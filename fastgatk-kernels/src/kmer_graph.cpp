#include "fastgatk/kernels/kmer_graph.hpp"

#include "fastgatk/core/plan.hpp"
#include "fastgatk/kernels/smith_waterman.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;
using KmerKey = std::string;

std::uint8_t encode(std::uint8_t base) {
    if (base < 4) return base;
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

KmerKey kmer_key(const std::vector<std::uint8_t>& bases,
                 std::size_t offset, std::uint32_t k) {
    static constexpr char symbols[] = {'A', 'C', 'G', 'T'};
    KmerKey key(k, 'N');
    for (std::uint32_t i = 0; i < k; ++i)
        key[i] = symbols[encode(bases[offset + i])];
    return key;
}

std::vector<std::uint8_t> orient_read(const std::vector<std::uint8_t>& bases,
                                      std::size_t begin, std::size_t end,
                                      bool reverse) {
    if (begin > end || end > bases.size()) return {};
    // ReadThreadingGraph.addRead() threads GATKRead.getBases() directly.
    // HTSlib's decoded BAM bases have the same stored record order, including
    // for FLAG 0x10 reads; the strand flag is alignment metadata rather than
    // an instruction to reverse-complement this graph payload.
    (void)reverse;
    std::vector<std::uint8_t> oriented(end - begin);
    for (std::size_t offset = 0; offset < oriented.size(); ++offset) {
        const auto value = encode(bases[begin + offset]);
        oriented[offset] = value < 4U
            ? value
            : static_cast<std::uint8_t>('N');
    }
    return oriented;
}

struct CigarElement {
    std::size_t length = 0;
    char operation = 0;
};

// ReadThreadingGraph's dangling-end merge is defined over a small CIGAR
// language rather than over a score threshold.  Keep that interpretation on
// Host with the variable-length traceback, while the corresponding affine DP
// score is checked against smith_waterman_score_kokkos below.
bool parse_cigar(const std::string& cigar, std::vector<CigarElement>& elements) {
    elements.clear();
    std::size_t index = 0;
    while (index < cigar.size()) {
        if (!std::isdigit(static_cast<unsigned char>(cigar[index]))) return false;
        std::size_t length = 0;
        while (index < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[index]))) {
            const auto digit = static_cast<std::size_t>(cigar[index] - '0');
            if (length > (std::numeric_limits<std::size_t>::max() - digit) / 10U)
                return false;
            length = length * 10U + digit;
            ++index;
        }
        if (length == 0 || index >= cigar.size()) return false;
        const auto operation = cigar[index++];
        if (operation != 'M' && operation != 'I' && operation != 'D') return false;
        elements.push_back(CigarElement{length, operation});
    }
    return !elements.empty();
}

std::size_t cigar_read_length(const std::vector<CigarElement>& elements) {
    std::size_t length = 0;
    for (const auto& element : elements)
        if (element.operation == 'M' || element.operation == 'I') length += element.length;
    return length;
}

std::size_t cigar_reference_length(const std::vector<CigarElement>& elements) {
    std::size_t length = 0;
    for (const auto& element : elements)
        if (element.operation == 'M' || element.operation == 'D') length += element.length;
    return length;
}

std::size_t longest_suffix_match(const std::string& reference, const std::string& alternate,
                                 const std::size_t reference_last) {
    if (reference.empty() || alternate.empty() || reference_last >= reference.size()) return 0;
    std::size_t length = 0;
    while (length < alternate.size() && length <= reference_last) {
        const auto reference_index = reference_last - length;
        const auto alternate_index = alternate.size() - 1U - length;
        if (encode(static_cast<std::uint8_t>(reference[reference_index])) !=
            encode(static_cast<std::uint8_t>(alternate[alternate_index])))
            break;
        ++length;
    }
    return length;
}

struct DanglingTailRecovery {
    std::string sequence;
    std::size_t alternate_merge_index = 0;
    std::size_t reference_merge_index = 0;
};

struct DanglingHeadMerge {
    std::size_t reference_index = 0;
    std::size_t alternate_index = 0;
};

// Mirrors AbstractReadThreadingGraph.mergeDanglingTail(): a valid
// LeadingIndel CIGAR must end in M, have at most three elements, and expose
// the configured exact suffix before reconnecting the dangling graph path to
// the reference suffix.  The returned sequence is the language accepted by
// the new graph edge; no numerical likelihood is calculated on Host here.
std::optional<DanglingTailRecovery> recover_dangling_tail_sequence(
    const std::string& alternate, const std::string& reference,
    const SmithWatermanAlignment& alignment, const std::int32_t minimum_matching_suffix) {
    std::vector<CigarElement> elements;
    if (!parse_cigar(alignment.cigar, elements)) return std::nullopt;
    // GATK AlignmentUtils.removeTrailingDeletions() happens before the
    // CIGAR complexity and terminal-M checks.
    if (!elements.empty() && elements.back().operation == 'D') elements.pop_back();
    if (elements.empty() || elements.size() > 3U || elements.back().operation != 'M')
        return std::nullopt;
    const auto read_length = cigar_read_length(elements);
    const auto reference_length = cigar_reference_length(elements);
    if (read_length == 0 || reference_length == 0 || read_length > alternate.size() ||
        reference_length > reference.size())
        return std::nullopt;
    const auto matching_suffix = std::min(
        longest_suffix_match(reference, alternate, reference_length - 1U),
        elements.back().length);
    // GATK's -1 legacy setting is shared by both dangling-end operations:
    // tails require a non-empty matching suffix, whereas non-negative values
    // require that exact threshold. Casting -1 to size_t here would reject
    // every recoverable tail and disconnect otherwise valid branches.
    if (minimum_matching_suffix < 0
            ? matching_suffix == 0U
            : matching_suffix < static_cast<std::size_t>(minimum_matching_suffix))
        return std::nullopt;
    // GATK uses Math.max(readLength - matchingSuffix - 1, 0): a complete
    // suffix match reconnects from the LCA instead of underflowing an
    // unsigned index and silently discarding the recovery.
    const auto alternate_merge_index = read_length <= matching_suffix
        ? std::size_t{0}
        : read_length - matching_suffix - 1U;
    auto reference_merge_index = reference_length - matching_suffix;
    // Preserve GATK's left-aligned leading-deletion correction. Without it a
    // deletion whose first reference base is the LCA loses one base when the
    // recovered edge is materialized as a Host path.
    if (elements.front().operation == 'D' &&
        elements.front().length + matching_suffix == reference_length)
        ++reference_merge_index;
    if (reference_merge_index == 0 || reference_merge_index >= reference.size() ||
        alternate_merge_index >= alternate.size())
        return std::nullopt;
    std::string recovered = alternate.substr(0, alternate_merge_index + 1U);
    recovered.append(reference, reference_merge_index, std::string::npos);
    // A full-M alignment can leave the emitted sequence text unchanged while
    // still requiring a new edge from the alternate path back to the
    // reference path.  `mergeDanglingTail()` always performs that topology
    // mutation after it has selected valid (non-zero) path indexes; rejecting
    // it merely because the two path languages spell the same bases drops a
    // legitimate rejoin and loses the corresponding K-best combinations.
    if (recovered.empty()) return std::nullopt;
    return DanglingTailRecovery{
        std::move(recovered), alternate_merge_index, reference_merge_index};
}

// Mirrors AbstractReadThreadingGraph.mergeDanglingHead().  GATK constructs
// both paths backwards from their highest common descendant, runs
// LEADING_INDEL SW, and then looks for an exact prefix in that backwards
// representation.  The returned indexes address the corresponding graph
// paths; translating those indexes into an edge remains Host topology work.
std::optional<DanglingHeadMerge> recover_dangling_head_merge(
    const std::string& alternate_upwards, const std::string& reference_upwards,
    const SmithWatermanAlignment& alignment, const std::size_t minimum_matching_bases) {
    std::vector<CigarElement> elements;
    if (!parse_cigar(alignment.cigar, elements)) return std::nullopt;
    // GATK removes terminal deletions before applying its bounded CIGAR
    // language.  For a head only the first element is required to be M.
    if (!elements.empty() && elements.back().operation == 'D') elements.pop_back();
    if (elements.empty() || elements.size() > 3U || elements.front().operation != 'M')
        return std::nullopt;
    const auto reference_length = cigar_reference_length(elements);
    if (reference_length == 0 || reference_length > reference_upwards.size() ||
        alternate_upwards.empty())
        return std::nullopt;

    std::size_t reference_index = reference_length - 1U;
    std::size_t alternate_index = alternate_upwards.size() - 1U;
    std::size_t matches = 0;
    bool mismatch = false;
    for (auto element = elements.rbegin(); element != elements.rend() && !mismatch; ++element) {
        if (element->operation != 'M') break;
        for (std::size_t base = 0; base < element->length; ++base) {
            if (reference_index >= reference_upwards.size() ||
                alternate_index >= alternate_upwards.size())
                return std::nullopt;
            if (encode(static_cast<std::uint8_t>(reference_upwards[reference_index])) !=
                encode(static_cast<std::uint8_t>(alternate_upwards[alternate_index]))) {
                mismatch = true;
                break;
            }
            ++matches;
            if (reference_index == 0 || alternate_index == 0) break;
            --reference_index;
            --alternate_index;
        }
    }
    if (matches < minimum_matching_bases || reference_index == 0 || alternate_index == 0)
        return std::nullopt;
    return DanglingHeadMerge{reference_index, alternate_index};
}

// Mirrors AbstractReadThreadingGraph.mergeDanglingHeadLegacy(), selected by
// GATK 4's default minMatchingBasesToDanglingEndRecovery == -1. It compares
// the leading M portion of the upward-path CIGAR and tolerates at most
// max(1, M/k) mismatches, returning the last mismatch position used to place
// the reconnecting edge. This is deliberately distinct from the modern
// exact-prefix algorithm above: changing the default silently admits noisy
// dangling heads that GATK rejects.
std::optional<DanglingHeadMerge> recover_dangling_head_legacy_merge(
    const std::string& alternate_upwards, const std::string& reference_upwards,
    const SmithWatermanAlignment& alignment, const std::uint32_t kmer_size) {
    std::vector<CigarElement> elements;
    if (!parse_cigar(alignment.cigar, elements)) return std::nullopt;
    if (!elements.empty() && elements.back().operation == 'D') elements.pop_back();
    if (elements.empty() || elements.size() > 3U || elements.front().operation != 'M')
        return std::nullopt;
    const auto maximum_index = static_cast<std::size_t>(elements.front().length);
    if (maximum_index == 0U || maximum_index > reference_upwards.size() ||
        maximum_index > alternate_upwards.size())
        return std::nullopt;
    const auto maximum_mismatches = std::max<std::size_t>(
        1U, maximum_index / std::max<std::uint32_t>(1U, kmer_size));
    std::size_t mismatches = 0;
    std::optional<std::size_t> last_mismatch;
    for (std::size_t index = 0; index < maximum_index; ++index) {
        if (reference_upwards[index] == alternate_upwards[index]) continue;
        if (++mismatches > maximum_mismatches) return std::nullopt;
        last_mismatch = index;
    }
    if (!last_mismatch.has_value() || *last_mismatch == 0U) return std::nullopt;
    return DanglingHeadMerge{*last_mismatch, *last_mismatch};
}

struct RewriteNode {
    std::string sequence;
    bool synthetic = false;
    bool alive = true;
    std::set<std::size_t> incoming;
    std::set<std::size_t> outgoing;
};

struct RewriteGraph {
    std::vector<RewriteNode> nodes;
    std::set<std::pair<std::size_t, std::size_t>> edges;

    std::size_t add_node(std::string sequence, const bool synthetic) {
        nodes.push_back(RewriteNode{std::move(sequence), synthetic, true, {}, {}});
        return nodes.size() - 1;
    }

    void add_edge(const std::size_t from, const std::size_t to) {
        if (from >= nodes.size() || to >= nodes.size() ||
            !nodes[from].alive || !nodes[to].alive) return;
        if (!edges.emplace(from, to).second) return;
        nodes[from].outgoing.insert(to);
        nodes[to].incoming.insert(from);
    }

    void remove_edge(const std::size_t from, const std::size_t to) {
        if (!edges.erase({from, to})) return;
        nodes[from].outgoing.erase(to);
        nodes[to].incoming.erase(from);
    }

    void remove_node(const std::size_t node) {
        if (node >= nodes.size() || !nodes[node].alive) return;
        const auto incoming = nodes[node].incoming;
        const auto outgoing = nodes[node].outgoing;
        for (const auto from : incoming) remove_edge(from, node);
        for (const auto to : outgoing) remove_edge(node, to);
        nodes[node].alive = false;
        nodes[node].incoming.clear();
        nodes[node].outgoing.clear();
    }
};

std::size_t common_prefix_length(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return 0;
    auto length = sequences.front().size();
    for (std::size_t index = 1; index < sequences.size(); ++index)
        length = std::min(length, sequences[index].size());
    std::size_t common = 0;
    while (common < length) {
        const auto base = sequences.front()[common];
        bool same = true;
        for (const auto& sequence : sequences)
            if (sequence[common] != base) { same = false; break; }
        if (!same) break;
        ++common;
    }
    return common;
}

std::size_t common_suffix_length(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return 0;
    auto length = sequences.front().size();
    for (std::size_t index = 1; index < sequences.size(); ++index)
        length = std::min(length, sequences[index].size());
    std::size_t common = 0;
    while (common < length) {
        const auto base = sequences.front()[sequences.front().size() - common - 1];
        bool same = true;
        for (const auto& sequence : sequences)
            if (sequence[sequence.size() - common - 1] != base) { same = false; break; }
        if (!same) break;
        ++common;
    }
    return common;
}

std::size_t active_sequence_node_count(const RewriteGraph& graph) {
    std::size_t count = 0;
    for (const auto& node : graph.nodes)
        if (node.alive && !node.synthetic) ++count;
    return count;
}

std::size_t active_sequence_edge_count(const RewriteGraph& graph) {
    std::size_t count = 0;
    for (const auto& edge : graph.edges)
        if (graph.nodes[edge.first].alive && graph.nodes[edge.second].alive &&
            !graph.nodes[edge.first].synthetic && !graph.nodes[edge.second].synthetic)
            ++count;
    return count;
}

// `ReadThreadingGraph.toSequenceGraph()` retains the multiplicity and
// reference state of every BaseEdge.  The earlier path-only rewriter below is
// useful telemetry, but cannot be used to choose K-best paths: after the 128
// path bound has been applied it no longer contains alternatives that a
// SeqGraph merge can promote.  Keep a small Host representation of the real
// SeqGraph here.  Kokkos still owns counting, degree calculation and pruning;
// these routines implement the variable-topology rewrites performed by the
// GATK Java host code.
struct WeightedSeqEdge {
    std::uint32_t multiplicity = 0;
    bool reference = false;
};

struct WeightedSeqNode {
    std::string sequence;
    std::uint32_t support = 0;
    bool alive = true;
    std::set<std::size_t> incoming;
    std::set<std::size_t> outgoing;
};

struct WeightedSeqGraph {
    std::vector<WeightedSeqNode> nodes;
    std::map<std::pair<std::size_t, std::size_t>, WeightedSeqEdge> edges;

    std::size_t add_node(std::string sequence, const std::uint32_t support = 0) {
        nodes.push_back(WeightedSeqNode{std::move(sequence), support, true, {}, {}});
        return nodes.size() - 1U;
    }

    const WeightedSeqEdge* edge(const std::size_t from, const std::size_t to) const {
        const auto found = edges.find({from, to});
        return found == edges.end() ? nullptr : &found->second;
    }

    bool add_edge(const std::size_t from, const std::size_t to, const WeightedSeqEdge value) {
        if (from >= nodes.size() || to >= nodes.size() || !nodes[from].alive || !nodes[to].alive)
            return false;
        if (!edges.emplace(std::make_pair(from, to), value).second) return false;
        nodes[from].outgoing.insert(to);
        nodes[to].incoming.insert(from);
        return true;
    }

    void add_or_update_edge(const std::size_t from, const std::size_t to,
                            const WeightedSeqEdge value) {
        const auto found = edges.find({from, to});
        if (found == edges.end()) {
            add_edge(from, to, value);
        } else {
            found->second.multiplicity += value.multiplicity;
            found->second.reference = found->second.reference || value.reference;
        }
    }

    void remove_edge(const std::size_t from, const std::size_t to) {
        if (!edges.erase({from, to})) return;
        nodes[from].outgoing.erase(to);
        nodes[to].incoming.erase(from);
    }

    void remove_node(const std::size_t node) {
        if (node >= nodes.size() || !nodes[node].alive) return;
        const auto incoming = nodes[node].incoming;
        const auto outgoing = nodes[node].outgoing;
        for (const auto from : incoming) remove_edge(from, node);
        for (const auto to : outgoing) remove_edge(node, to);
        nodes[node].alive = false;
        nodes[node].incoming.clear();
        nodes[node].outgoing.clear();
    }

    bool reference_node(const std::size_t node) const {
        if (node >= nodes.size() || !nodes[node].alive) return false;
        for (const auto from : nodes[node].incoming) {
            const auto* value = edge(from, node);
            if (value != nullptr && value->reference) return true;
        }
        for (const auto to : nodes[node].outgoing) {
            const auto* value = edge(node, to);
            if (value != nullptr && value->reference) return true;
        }
        std::size_t active = 0;
        for (const auto& value : nodes) if (value.alive) ++active;
        return active == 1U;
    }

    std::optional<std::size_t> reference_source() const {
        // SeqGraph.zipLinearChains() replaces a reference chain with one
        // vertex and leaves the old vertices in this backing store marked
        // dead.  Source/sink identity is a property of the live graph, not
        // the historical allocation count: a fully zipped reference path is
        // simultaneously its source and sink.
        const auto active_nodes = static_cast<std::size_t>(std::count_if(
            nodes.begin(), nodes.end(), [](const auto& node) { return node.alive; }));
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            if (!nodes[node].alive) continue;
            bool incoming_reference = false;
            bool outgoing_reference = false;
            for (const auto from : nodes[node].incoming) {
                const auto* value = edge(from, node);
                incoming_reference = incoming_reference || (value != nullptr && value->reference);
            }
            for (const auto to : nodes[node].outgoing) {
                const auto* value = edge(node, to);
                outgoing_reference = outgoing_reference || (value != nullptr && value->reference);
            }
            if (!incoming_reference && (outgoing_reference || active_nodes == 1U)) return node;
        }
        return std::nullopt;
    }

    std::optional<std::size_t> reference_sink() const {
        const auto active_nodes = static_cast<std::size_t>(std::count_if(
            nodes.begin(), nodes.end(), [](const auto& node) { return node.alive; }));
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            if (!nodes[node].alive) continue;
            bool incoming_reference = false;
            bool outgoing_reference = false;
            for (const auto from : nodes[node].incoming) {
                const auto* value = edge(from, node);
                incoming_reference = incoming_reference || (value != nullptr && value->reference);
            }
            for (const auto to : nodes[node].outgoing) {
                const auto* value = edge(node, to);
                outgoing_reference = outgoing_reference || (value != nullptr && value->reference);
            }
            if (!outgoing_reference && (incoming_reference || active_nodes == 1U)) return node;
        }
        return std::nullopt;
    }
};

std::size_t weighted_seqgraph_node_count(const WeightedSeqGraph& graph) {
    return static_cast<std::size_t>(std::count_if(graph.nodes.begin(), graph.nodes.end(),
        [](const auto& node) { return node.alive; }));
}

std::size_t weighted_seqgraph_edge_count(const WeightedSeqGraph& graph) {
    return graph.edges.size();
}

WeightedSeqEdge or_edge(const std::vector<WeightedSeqEdge>& edges,
                        const std::uint32_t multiplicity = 1U) {
    WeightedSeqEdge result{multiplicity, false};
    for (const auto& edge : edges) result.reference = result.reference || edge.reference;
    return result;
}

std::size_t common_suffix_length_limited(const std::vector<std::string>& sequences,
                                         const std::size_t limit) {
    if (sequences.empty()) return 0;
    std::size_t length = std::min(limit, sequences.front().size());
    for (std::size_t index = 1; index < sequences.size(); ++index)
        length = std::min(length, sequences[index].size());
    std::size_t common = 0;
    while (common < length) {
        const auto base = sequences.front()[sequences.front().size() - common - 1U];
        bool same = true;
        for (const auto& sequence : sequences)
            if (sequence[sequence.size() - common - 1U] != base) { same = false; break; }
        if (!same) break;
        ++common;
    }
    return common;
}

bool weighted_seqgraph_zip_linear_chains(WeightedSeqGraph& graph) {
    std::vector<std::size_t> starts;
    for (std::size_t node = 0; node < graph.nodes.size(); ++node) {
        if (!graph.nodes[node].alive || graph.nodes[node].outgoing.size() != 1U) continue;
        if (graph.nodes[node].incoming.size() != 1U ||
            graph.nodes[*graph.nodes[node].incoming.begin()].outgoing.size() > 1U)
            starts.push_back(node);
    }
    bool changed = false;
    for (const auto start : starts) {
        if (!graph.nodes[start].alive) continue;
        std::vector<std::size_t> chain{start};
        auto last = start;
        auto last_is_reference = graph.reference_node(last);
        while (graph.nodes[last].outgoing.size() == 1U) {
            const auto target = *graph.nodes[last].outgoing.begin();
            if (target == last || graph.nodes[target].incoming.size() != 1U ||
                graph.reference_node(target) != last_is_reference)
                break;
            chain.push_back(target);
            last = target;
            last_is_reference = graph.reference_node(last);
        }
        if (chain.size() <= 1U) continue;
        const auto first = chain.front();
        const auto final = chain.back();
        std::string sequence;
        std::uint32_t support = 0;
        for (const auto node : chain) {
            sequence += graph.nodes[node].sequence;
            support = std::max(support, graph.nodes[node].support);
        }
        const auto merged = graph.add_node(std::move(sequence), support);
        const auto outgoing = graph.nodes[final].outgoing;
        const auto incoming = graph.nodes[first].incoming;
        for (const auto to : outgoing) {
            const auto* edge = graph.edge(final, to);
            if (edge != nullptr) graph.add_edge(merged, to, *edge);
        }
        for (const auto from : incoming) {
            const auto* edge = graph.edge(from, first);
            if (edge != nullptr) graph.add_edge(from, merged, *edge);
        }
        for (const auto node : chain) graph.remove_node(node);
        changed = true;
    }
    return changed;
}

// Port of SharedVertexSequenceSplitter.splitAndUpdate().  `top`/`bottom`
// deliberately remain optional because MergeTails has no bottom vertex.
bool weighted_seqgraph_split_shared_vertices(WeightedSeqGraph& graph,
                                             const std::vector<std::size_t>& middles,
                                             const std::optional<std::size_t> top,
                                             const std::optional<std::size_t> bottom) {
    if (middles.size() < 2U || (!top.has_value() && !bottom.has_value())) return false;
    std::vector<std::string> sequences;
    sequences.reserve(middles.size());
    for (const auto middle : middles) {
        if (middle >= graph.nodes.size() || !graph.nodes[middle].alive) return false;
        sequences.push_back(graph.nodes[middle].sequence);
    }
    const auto prefix_length = common_prefix_length(sequences);
    const auto minimum_length = std::min_element(sequences.begin(), sequences.end(),
        [](const auto& left, const auto& right) { return left.size() < right.size(); })->size();
    const auto suffix_length = common_suffix_length_limited(sequences,
        minimum_length >= prefix_length ? minimum_length - prefix_length : 0U);
    const auto prefix = sequences.front().substr(0, prefix_length);
    const auto suffix = suffix_length == 0U ? std::string{}
        : sequences.front().substr(sequences.front().size() - suffix_length);

    struct Branch { std::string residual; WeightedSeqEdge to; WeightedSeqEdge from; std::uint32_t support; };
    std::vector<Branch> branches;
    branches.reserve(middles.size());
    for (std::size_t index = 0; index < middles.size(); ++index) {
        const auto middle = middles[index];
        WeightedSeqEdge to{0U, graph.reference_node(middle)};
        WeightedSeqEdge from{0U, graph.reference_node(middle)};
        if (!graph.nodes[middle].incoming.empty()) {
            const auto predecessor = *graph.nodes[middle].incoming.begin();
            if (const auto* edge = graph.edge(predecessor, middle); edge != nullptr) to = *edge;
        }
        if (!graph.nodes[middle].outgoing.empty()) {
            const auto successor = *graph.nodes[middle].outgoing.begin();
            if (const auto* edge = graph.edge(middle, successor); edge != nullptr) from = *edge;
        }
        const auto residual_length = sequences[index].size() - prefix_length - suffix_length;
        branches.push_back(Branch{sequences[index].substr(prefix_length, residual_length), to, from,
                                  graph.nodes[middle].support});
    }
    for (const auto middle : middles) graph.remove_node(middle);

    // Build the temporary split graph's prefix/suffix edge lists.  We only
    // need its edges, not a separate graph object, because updateGraph()
    // immediately reconnects them to the outer graph.
    struct ResidualBranch { std::optional<std::size_t> node; WeightedSeqEdge to; WeightedSeqEdge from; };
    std::vector<ResidualBranch> residuals;
    std::vector<WeightedSeqEdge> prefix_outgoing;
    std::vector<WeightedSeqEdge> suffix_incoming;
    std::optional<WeightedSeqEdge> direct;
    for (const auto& branch : branches) {
        if (branch.residual.empty()) {
            WeightedSeqEdge combined = branch.to;
            combined.multiplicity += branch.from.multiplicity;
            combined.reference = combined.reference || branch.from.reference;
            if (direct.has_value()) {
                direct->multiplicity += combined.multiplicity;
                direct->reference = direct->reference || combined.reference;
            } else direct = combined;
            prefix_outgoing.push_back(combined);
        } else {
            const auto node = graph.add_node(branch.residual, branch.support);
            residuals.push_back(ResidualBranch{node, branch.to, branch.from});
            prefix_outgoing.push_back(branch.to);
            suffix_incoming.push_back(branch.from);
        }
    }
    const bool has_direct = direct.has_value();
    if (direct.has_value()) suffix_incoming.push_back(*direct);
    const bool only_direct = has_direct && residuals.empty();
    const bool need_prefix = !prefix.empty() || (!top.has_value() && !only_direct);
    const bool need_suffix = !suffix.empty() || (!bottom.has_value() && !only_direct);
    const auto prefix_node = need_prefix ? std::optional<std::size_t>(graph.add_node(prefix)) : top;
    const auto suffix_node = need_suffix ? std::optional<std::size_t>(graph.add_node(suffix)) : bottom;
    if (!prefix_node.has_value() || !suffix_node.has_value()) return false;
    if (need_prefix && top.has_value()) graph.add_edge(*top, *prefix_node, or_edge(prefix_outgoing));
    if (need_suffix && bottom.has_value()) graph.add_edge(*suffix_node, *bottom, or_edge(suffix_incoming));
    for (const auto& branch : residuals) {
        graph.add_edge(*prefix_node, *branch.node, branch.to);
        graph.add_edge(*branch.node, *suffix_node, branch.from);
    }
    // In updateGraph() a prefix->suffix edge is connected to `botForConnect`.
    // With a dangling tail that endpoint *is* the newly created suffix, so
    // this edge must be retained even when the original bottom was null.
    if (direct.has_value())
        graph.add_edge(*prefix_node, *suffix_node, *direct);
    return true;
}

bool weighted_seqgraph_merge_diamonds_once(WeightedSeqGraph& graph) {
    for (std::size_t top = 0; top < graph.nodes.size(); ++top) {
        if (!graph.nodes[top].alive || graph.nodes[top].outgoing.size() <= 1U) continue;
        std::vector<std::size_t> middles(graph.nodes[top].outgoing.begin(), graph.nodes[top].outgoing.end());
        std::optional<std::size_t> bottom;
        bool valid = true;
        for (const auto middle : middles) {
            if (graph.nodes[middle].outgoing.empty() || graph.nodes[middle].incoming.size() != 1U) {
                valid = false; break;
            }
            for (const auto target : graph.nodes[middle].outgoing) {
                if (!bottom.has_value()) bottom = target;
                else if (*bottom != target) { valid = false; break; }
            }
            if (!valid) break;
        }
        if (!valid || !bottom.has_value() || graph.nodes[*bottom].incoming.size() != middles.size()) continue;
        std::vector<std::string> sequences;
        for (const auto middle : middles) sequences.push_back(graph.nodes[middle].sequence);
        if (common_prefix_length(sequences) == 0U && common_suffix_length(sequences) == 0U) continue;
        return weighted_seqgraph_split_shared_vertices(graph, middles, top, bottom);
    }
    return false;
}

bool weighted_seqgraph_merge_tails_once(WeightedSeqGraph& graph) {
    constexpr std::size_t min_common = 10U;
    for (std::size_t top = 0; top < graph.nodes.size(); ++top) {
        if (!graph.nodes[top].alive || graph.nodes[top].outgoing.size() <= 1U) continue;
        std::vector<std::size_t> tails(graph.nodes[top].outgoing.begin(), graph.nodes[top].outgoing.end());
        bool valid = true;
        for (const auto tail : tails) {
            if (!graph.nodes[tail].outgoing.empty() || graph.nodes[tail].incoming.size() > 1U) {
                valid = false; break;
            }
        }
        if (!valid) continue;
        std::vector<std::string> sequences;
        for (const auto tail : tails) sequences.push_back(graph.nodes[tail].sequence);
        if (common_suffix_length(sequences) < min_common) continue;
        return weighted_seqgraph_split_shared_vertices(graph, tails, top, std::nullopt);
    }
    return false;
}

bool weighted_seqgraph_split_common_suffix_once(WeightedSeqGraph& graph) {
    for (std::size_t bottom = 0; bottom < graph.nodes.size(); ++bottom) {
        if (!graph.nodes[bottom].alive || graph.nodes[bottom].incoming.size() < 2U) continue;
        std::vector<std::size_t> middles(graph.nodes[bottom].incoming.begin(), graph.nodes[bottom].incoming.end());
        bool safe = true;
        for (const auto middle : middles) {
            if (middle == bottom || graph.nodes[middle].outgoing.size() != 1U ||
                *graph.nodes[middle].outgoing.begin() != bottom ||
                graph.nodes[bottom].outgoing.count(middle) != 0U) {
                safe = false; break;
            }
        }
        if (!safe) continue;
        std::vector<std::string> sequences;
        for (const auto middle : middles) sequences.push_back(graph.nodes[middle].sequence);
        const auto suffix_length = common_suffix_length(sequences);
        if (suffix_length == 0U) continue;
        bool eliminates_ref_source = false;
        bool all_suffix = true;
        for (const auto middle : middles) {
            const bool incoming_reference = std::any_of(graph.nodes[middle].incoming.begin(), graph.nodes[middle].incoming.end(),
                [&](const auto from) { const auto* edge = graph.edge(from, middle); return edge != nullptr && edge->reference; });
            const bool outgoing_reference = std::any_of(graph.nodes[middle].outgoing.begin(), graph.nodes[middle].outgoing.end(),
                [&](const auto to) { const auto* edge = graph.edge(middle, to); return edge != nullptr && edge->reference; });
            eliminates_ref_source = eliminates_ref_source || (!incoming_reference && outgoing_reference &&
                graph.nodes[middle].sequence.size() == suffix_length);
            all_suffix = all_suffix && graph.nodes[middle].sequence.size() == suffix_length;
        }
        if (eliminates_ref_source || all_suffix) continue;
        const auto suffix = sequences.front().substr(sequences.front().size() - suffix_length);
        struct Branch { std::string prefix; WeightedSeqEdge out; std::vector<std::pair<std::size_t, WeightedSeqEdge>> incoming; std::uint32_t support; };
        std::vector<Branch> branches;
        for (const auto middle : middles) {
            const auto* out = graph.edge(middle, bottom);
            if (out == nullptr) { safe = false; break; }
            Branch branch{graph.nodes[middle].sequence.substr(0, graph.nodes[middle].sequence.size() - suffix_length),
                          *out, {}, graph.nodes[middle].support};
            for (const auto from : graph.nodes[middle].incoming) {
                const auto* in = graph.edge(from, middle);
                if (in != nullptr) branch.incoming.emplace_back(from, *in);
            }
            branches.push_back(std::move(branch));
        }
        if (!safe) continue;
        for (const auto middle : middles) graph.remove_node(middle);
        for (const auto& branch : branches) {
            const auto suffix_node = graph.add_node(suffix);
            std::size_t incoming_target = suffix_node;
            if (!branch.prefix.empty()) {
                incoming_target = graph.add_node(branch.prefix, branch.support);
                graph.add_edge(incoming_target, suffix_node, WeightedSeqEdge{1U, branch.out.reference});
            }
            graph.add_edge(suffix_node, bottom, branch.out);
            for (const auto& [from, edge] : branch.incoming) graph.add_edge(from, incoming_target, edge);
        }
        return true;
    }
    return false;
}

bool weighted_seqgraph_merge_common_suffix_once(WeightedSeqGraph& graph) {
    for (std::size_t bottom = 0; bottom < graph.nodes.size(); ++bottom) {
        if (!graph.nodes[bottom].alive || graph.nodes[bottom].incoming.empty()) continue;
        std::vector<std::size_t> previous(graph.nodes[bottom].incoming.begin(), graph.nodes[bottom].incoming.end());
        const auto& sequence = graph.nodes[previous.front()].sequence;
        bool valid = true;
        for (const auto prior : previous) {
            if (graph.nodes[prior].sequence != sequence || graph.nodes[prior].outgoing.size() != 1U ||
                *graph.nodes[prior].outgoing.begin() != bottom || graph.nodes[prior].incoming.empty()) {
                valid = false; break;
            }
        }
        if (!valid) continue;
        const auto merged = graph.add_node(sequence + graph.nodes[bottom].sequence);
        for (const auto prior : previous) {
            for (const auto from : graph.nodes[prior].incoming) {
                if (const auto* edge = graph.edge(from, prior); edge != nullptr)
                    graph.add_edge(from, merged, *edge);
            }
        }
        for (const auto to : graph.nodes[bottom].outgoing) {
            if (const auto* edge = graph.edge(bottom, to); edge != nullptr)
                graph.add_edge(merged, to, *edge);
        }
        for (const auto prior : previous) graph.remove_node(prior);
        graph.remove_node(bottom);
        return true;
    }
    return false;
}

void weighted_seqgraph_simplify(WeightedSeqGraph& graph) {
    weighted_seqgraph_zip_linear_chains(graph);
    for (std::size_t cycle = 0; cycle < 100U; ++cycle) {
        bool changed = false;
        while (weighted_seqgraph_merge_diamonds_once(graph)) changed = true;
        while (weighted_seqgraph_merge_tails_once(graph)) changed = true;
        while (weighted_seqgraph_split_common_suffix_once(graph)) changed = true;
        while (weighted_seqgraph_merge_common_suffix_once(graph)) changed = true;
        changed = weighted_seqgraph_zip_linear_chains(graph) || changed;
        if (!changed) break;
    }
}

void weighted_seqgraph_remove_undirected_nonreference(WeightedSeqGraph& graph) {
    const auto source = graph.reference_source();
    if (!source.has_value()) return;
    std::vector<bool> visited(graph.nodes.size(), false);
    std::deque<std::size_t> pending{*source};
    visited[*source] = true;
    while (!pending.empty()) {
        const auto node = pending.front();
        pending.pop_front();
        for (const auto neighbor : graph.nodes[node].incoming) {
            if (!visited[neighbor]) { visited[neighbor] = true; pending.push_back(neighbor); }
        }
        for (const auto neighbor : graph.nodes[node].outgoing) {
            if (!visited[neighbor]) { visited[neighbor] = true; pending.push_back(neighbor); }
        }
    }
    for (std::size_t node = 0; node < graph.nodes.size(); ++node)
        if (graph.nodes[node].alive && !visited[node]) graph.remove_node(node);
}

bool weighted_seqgraph_remove_paths_not_connected_to_ref(WeightedSeqGraph& graph) {
    const auto source = graph.reference_source();
    const auto sink = graph.reference_sink();
    if (!source.has_value() || !sink.has_value()) return false;
    std::vector<bool> forward(graph.nodes.size(), false), reverse(graph.nodes.size(), false);
    std::deque<std::size_t> pending{*source};
    forward[*source] = true;
    while (!pending.empty()) {
        const auto node = pending.front(); pending.pop_front();
        for (const auto to : graph.nodes[node].outgoing)
            if (!forward[to]) { forward[to] = true; pending.push_back(to); }
    }
    pending.push_back(*sink);
    reverse[*sink] = true;
    while (!pending.empty()) {
        const auto node = pending.front(); pending.pop_front();
        for (const auto from : graph.nodes[node].incoming)
            if (!reverse[from]) { reverse[from] = true; pending.push_back(from); }
    }
    for (std::size_t node = 0; node < graph.nodes.size(); ++node)
        if (graph.nodes[node].alive && (!forward[node] || !reverse[node])) graph.remove_node(node);
    return true;
}

struct ExactSeqGraphKBestPath {
    std::string sequence;
    double score = 0.0;
    std::uint32_t support = 0;
    bool alternate = false;
    std::size_t alt_begin = std::numeric_limits<std::size_t>::max();
    std::size_t alt_end = 0;
};

struct ExactSeqGraphKBestResult {
    std::vector<ExactSeqGraphKBestPath> paths;
    std::size_t nodes = 0;
    std::size_t edges = 0;
};

ExactSeqGraphKBestResult exact_seqgraph_kbest(
    const std::vector<std::string>& kmers,
    const std::vector<std::uint32_t>& node_support,
    const std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint8_t>& edge_kind,
    const std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t>& edge_support,
    const std::vector<bool>& keep,
    const std::uint32_t max_paths) {
    ExactSeqGraphKBestResult result;
    WeightedSeqGraph graph;
    std::vector<std::size_t> node_map(kmers.size(), std::numeric_limits<std::size_t>::max());
    std::vector<std::uint32_t> indegree(kmers.size(), 0U);
    for (const auto& [edge, kind] : edge_kind) {
        (void)kind;
        if (edge.second < indegree.size()) ++indegree[edge.second];
    }
    for (std::size_t node = 0; node < kmers.size(); ++node) {
        if (node >= keep.size() || !keep[node]) continue;
        // BaseGraph.toSequenceGraph() emits a full k-mer only for source
        // vertices; every other vertex contributes its additional suffix.
        const auto sequence = indegree[node] == 0U ? kmers[node]
            : (kmers[node].empty() ? std::string{} : kmers[node].substr(kmers[node].size() - 1U));
        node_map[node] = graph.add_node(sequence, node < node_support.size() ? node_support[node] : 0U);
    }
    for (const auto& [pair, kind] : edge_kind) {
        const auto from = pair.first;
        const auto to = pair.second;
        if (from >= node_map.size() || to >= node_map.size() ||
            node_map[from] == std::numeric_limits<std::size_t>::max() ||
            node_map[to] == std::numeric_limits<std::size_t>::max()) continue;
        const auto support = edge_support.find(pair);
        const auto multiplicity = static_cast<std::uint32_t>(
            (support == edge_support.end() ? 0U : support->second) + ((kind & 1U) != 0U ? 1U : 0U));
        graph.add_edge(node_map[from], node_map[to], WeightedSeqEdge{multiplicity, (kind & 1U) != 0U});
    }
    const auto debug_seqgraph = std::getenv("FASTGATK_DEBUG_EVENTMAP_SEQGRAPH") != nullptr;
    if (debug_seqgraph)
        std::cerr << "[FASTGATK_EVENTMAP_SEQGRAPH] stage=initial nodes="
                  << weighted_seqgraph_node_count(graph) << " edges="
                  << weighted_seqgraph_edge_count(graph) << '\n';
    // cleanupSeqGraph(): zip, remove weakly disconnected nodes, simplify,
    // then a directed ref-source/ref-sink trim and one final simplification.
    weighted_seqgraph_zip_linear_chains(graph);
    if (debug_seqgraph)
        std::cerr << "[FASTGATK_EVENTMAP_SEQGRAPH] stage=zipped nodes="
                  << weighted_seqgraph_node_count(graph) << " edges="
                  << weighted_seqgraph_edge_count(graph) << '\n';
    weighted_seqgraph_remove_undirected_nonreference(graph);
    weighted_seqgraph_simplify(graph);
    if (debug_seqgraph)
        std::cerr << "[FASTGATK_EVENTMAP_SEQGRAPH] stage=merged nodes="
                  << weighted_seqgraph_node_count(graph) << " edges="
                  << weighted_seqgraph_edge_count(graph) << '\n';
    if (weighted_seqgraph_remove_paths_not_connected_to_ref(graph)) weighted_seqgraph_simplify(graph);
    result.nodes = weighted_seqgraph_node_count(graph);
    result.edges = weighted_seqgraph_edge_count(graph);
    if (debug_seqgraph)
        std::cerr << "[FASTGATK_EVENTMAP_SEQGRAPH] stage=final nodes=" << result.nodes
                  << " edges=" << result.edges << '\n';
    const auto source = graph.reference_source();
    const auto sink = graph.reference_sink();
    if (!source.has_value() || !sink.has_value()) return result;

    struct PendingPath {
        std::size_t node = 0;
        std::string sequence;
        double score = 0.0;
        std::uint32_t support = 0;
        bool alternate = false;
        std::size_t alt_begin = std::numeric_limits<std::size_t>::max();
        std::size_t alt_end = 0;
        std::vector<std::size_t> vertices;
    };
    struct KBestOrder {
        bool operator()(const PendingPath& left, const PendingPath& right) const {
            if (left.score != right.score) return left.score < right.score;
            // GraphBasedKBestHaplotypeFinder reverses BaseUtils' lexical
            // comparator, so lexically larger bases break an exact score tie.
            if (left.sequence != right.sequence) return left.sequence < right.sequence;
            return left.node > right.node;
        }
    };
    std::priority_queue<PendingPath, std::vector<PendingPath>, KBestOrder> queue;
    queue.push(PendingPath{*source, graph.nodes[*source].sequence, 0.0,
                           graph.nodes[*source].support, false,
                           std::numeric_limits<std::size_t>::max(), 0U, {*source}});
    std::vector<std::size_t> vertex_visits(graph.nodes.size(), 0U);
    while (!queue.empty() && result.paths.size() < max_paths) {
        auto path = queue.top();
        queue.pop();
        if (path.node == *sink) {
            result.paths.push_back(ExactSeqGraphKBestPath{std::move(path.sequence), path.score,
                path.support, path.alternate, path.alt_begin, path.alt_end});
            continue;
        }
        if (vertex_visits[path.node]++ >= max_paths) continue;
        std::uint64_t total = 0;
        for (const auto to : graph.nodes[path.node].outgoing) {
            const auto* edge = graph.edge(path.node, to);
            if (edge != nullptr) total += edge->multiplicity;
        }
        if (total == 0U) continue;
        for (const auto to : graph.nodes[path.node].outgoing) {
            const auto* edge = graph.edge(path.node, to);
            if (edge == nullptr || edge->multiplicity == 0U) continue;
            // The ReadThreadingAssembler has already rejected cyclic input;
            // retain this guard for malformed direct API callers.
            if (std::find(path.vertices.begin(), path.vertices.end(), to) != path.vertices.end()) continue;
            auto next = path;
            next.node = to;
            const auto previous_length = next.sequence.size();
            next.sequence += graph.nodes[to].sequence;
            next.score += std::log10(static_cast<double>(edge->multiplicity)) -
                std::log10(static_cast<double>(total));
            next.support = std::max({next.support, graph.nodes[to].support, edge->multiplicity});
            if (!edge->reference) {
                next.alternate = true;
                if (next.alt_begin == std::numeric_limits<std::size_t>::max())
                    next.alt_begin = previous_length;
                // A dangling-tail recovery edge is deliberately
                // non-reference in GATK, but its target is the first
                // reference vertex at which the branch rejoins.  SeqGraph
                // may have zipped that target together with the complete
                // downstream reference suffix.  Its whole sequence is not
                // part of the read-derived branch: retain only the joining
                // base in the compact branch span.  This provenance is used
                // by the Host EventMap projection to distinguish a genuine
                // local indel from the reference halo.
                next.alt_end = graph.reference_node(to)
                    ? previous_length + std::min<std::size_t>(
                        1U, graph.nodes[to].sequence.size())
                    : next.sequence.size();
            }
            next.vertices.push_back(to);
            queue.push(std::move(next));
        }
    }
    return result;
}

}  // namespace

SeqGraphSimplificationResult simplify_seqgraph_paths(
    const std::vector<SeqGraphPath>& paths,
    const std::size_t min_common_tail_bases) {
    SeqGraphSimplificationResult result;
    RewriteGraph graph;
    // GATK's SeqGraph uses shared source/sink vertices to expose diamonds;
    // group by the coordinate frame carried by materialized paths.  A path
    // without coordinates gets a private source and is therefore treated as
    // an unanchored tail rather than accidentally merged with another locus.
    std::map<std::tuple<std::int32_t, std::int32_t, std::size_t>, std::size_t> sources;
    std::map<std::pair<std::int32_t, std::int32_t>, std::size_t> sinks;
    std::vector<std::size_t> source_nodes;
    for (std::size_t path_index = 0; path_index < paths.size(); ++path_index) {
        const auto& path = paths[path_index];
        const auto source_key = path.tid >= 0 && path.start >= 0
            ? std::make_tuple(path.tid, path.start, std::size_t{0})
            : std::make_tuple(std::int32_t{-1}, std::int32_t{-1}, path_index + 1);
        auto source_iter = sources.find(source_key);
        if (source_iter == sources.end()) {
            const auto source = graph.add_node({}, true);
            source_iter = sources.emplace(source_key, source).first;
            source_nodes.push_back(source);
        }
        const auto middle = graph.add_node(path.sequence, false);
        graph.add_edge(source_iter->second, middle);
        if (path.tid >= 0 && path.end >= 0) {
            const auto sink_key = std::make_pair(path.tid, path.end);
            auto sink_iter = sinks.find(sink_key);
            if (sink_iter == sinks.end())
                sink_iter = sinks.emplace(sink_key, graph.add_node({}, true)).first;
            graph.add_edge(middle, sink_iter->second);
        }
    }
    result.initial_nodes = active_sequence_node_count(graph);
    result.initial_edges = active_sequence_edge_count(graph);

    auto add_sequence_node = [&](const std::string& sequence) {
        return graph.add_node(sequence, false);
    };
    auto connect_branch = [&](const std::size_t top, const std::size_t bottom,
                              const std::string& prefix,
                              const std::string& residual,
                              const std::string& suffix) {
        std::size_t current = top;
        if (!prefix.empty()) {
            const auto node = add_sequence_node(prefix);
            graph.add_edge(current, node);
            current = node;
        }
        if (!residual.empty()) {
            const auto node = add_sequence_node(residual);
            graph.add_edge(current, node);
            current = node;
        }
        if (!suffix.empty()) {
            const auto node = add_sequence_node(suffix);
            graph.add_edge(current, node);
            current = node;
        }
        graph.add_edge(current, bottom);
    };

    // The operation order follows SeqGraph.simplifyGraph: common linear
    // chains are zipped, diamonds are split/merged, and tails are then
    // merged.  A bounded fixpoint handles a diamond that becomes visible only
    // after a neighbouring tail/suffix rewrite.
    for (std::size_t pass = 0; pass < 8; ++pass) {
        bool changed = false;
        for (std::size_t top = 0; top < graph.nodes.size(); ++top) {
            if (!graph.nodes[top].alive || graph.nodes[top].outgoing.size() < 2) continue;
            std::map<std::size_t, std::vector<std::size_t>> by_bottom;
            for (const auto middle : graph.nodes[top].outgoing) {
                if (!graph.nodes[middle].alive || graph.nodes[middle].synthetic ||
                    graph.nodes[middle].incoming.size() != 1 ||
                    graph.nodes[middle].outgoing.size() != 1) continue;
                const auto bottom = *graph.nodes[middle].outgoing.begin();
                by_bottom[bottom].push_back(middle);
            }
            for (const auto& [bottom, middles] : by_bottom) {
                if (middles.size() < 2 || !graph.nodes[bottom].alive ||
                    graph.nodes[bottom].incoming.size() != middles.size()) continue;
                std::vector<std::string> sequences;
                sequences.reserve(middles.size());
                for (const auto middle : middles) sequences.push_back(graph.nodes[middle].sequence);
                const auto prefix_length = common_prefix_length(sequences);
                const auto suffix_length = common_suffix_length(sequences);
                if (prefix_length == 0 && suffix_length == 0) continue;
                // The shared prefix/suffix must not overlap.  When a branch is
                // entirely shared, retain one base in the middle so the
                // branch vertex remains a valid path separator.
                auto safe_suffix = suffix_length;
                if (prefix_length + safe_suffix > sequences.front().size())
                    safe_suffix = sequences.front().size() - prefix_length;
                const auto prefix = sequences.front().substr(0, prefix_length);
                const auto suffix = safe_suffix == 0 ? std::string{}
                    : sequences.front().substr(sequences.front().size() - safe_suffix);
                for (const auto middle : middles) {
                    graph.remove_edge(top, middle);
                    graph.remove_edge(middle, bottom);
                    graph.remove_node(middle);
                }
                std::size_t prefix_node = std::numeric_limits<std::size_t>::max();
                std::size_t suffix_node = std::numeric_limits<std::size_t>::max();
                if (!prefix.empty()) {
                    prefix_node = add_sequence_node(prefix);
                    graph.add_edge(top, prefix_node);
                }
                if (!suffix.empty()) suffix_node = add_sequence_node(suffix);
                for (std::size_t index = 0; index < sequences.size(); ++index) {
                    const auto& sequence = sequences[index];
                    const auto residual_length = sequence.size() - prefix_length - safe_suffix;
                    const auto residual = sequence.substr(prefix_length, residual_length);
                    // Each branch gets its own residual; the common prefix and
                    // suffix are represented once by the shared topology.
                    std::size_t current = prefix_node == std::numeric_limits<std::size_t>::max()
                        ? top : prefix_node;
                    if (!residual.empty()) {
                        const auto node = add_sequence_node(residual);
                        graph.add_edge(current, node);
                        current = node;
                    }
                    if (!suffix.empty()) {
                        graph.add_edge(current, suffix_node);
                        current = suffix_node;
                    }
                    graph.add_edge(current, bottom);
                }
                ++result.diamond_merges;
                if (safe_suffix > 0) {
                    ++result.suffix_splits;
                    ++result.suffix_merges;
                }
                changed = true;
                break;
            }
            if (changed) break;
        }
        if (changed) continue;

        // Merge tails with the same source and a sufficiently long common
        // suffix, matching MergeTails' sink-only precondition.
        for (std::size_t top = 0; top < graph.nodes.size() && !changed; ++top) {
            if (!graph.nodes[top].alive || graph.nodes[top].outgoing.size() < 2) continue;
            std::vector<std::size_t> tails;
            for (const auto child : graph.nodes[top].outgoing)
                if (graph.nodes[child].alive && !graph.nodes[child].synthetic &&
                    graph.nodes[child].incoming.size() == 1 && graph.nodes[child].outgoing.empty())
                    tails.push_back(child);
            if (tails.size() < 2) continue;
            std::vector<std::string> sequences;
            for (const auto tail : tails) sequences.push_back(graph.nodes[tail].sequence);
            const auto suffix_length = common_suffix_length(sequences);
            if (suffix_length < min_common_tail_bases || suffix_length == 0) continue;
            const auto suffix = sequences.front().substr(sequences.front().size() - suffix_length);
            std::size_t suffix_node = std::numeric_limits<std::size_t>::max();
            for (const auto tail : tails) {
                graph.remove_edge(top, tail);
                graph.remove_node(tail);
            }
            suffix_node = add_sequence_node(suffix);
            for (const auto& sequence : sequences) {
                const auto residual = sequence.substr(0, sequence.size() - suffix_length);
                connect_branch(top, suffix_node, {}, residual, {});
            }
            ++result.tail_merges;
            ++result.suffix_splits;
            ++result.suffix_merges;
            changed = true;
        }
        if (changed) continue;

        // Zip only genuine sequence vertices.  Synthetic source/sink vertices
        // remain graph boundaries, so a complete linear path is represented
        // by one sequence vertex and still reports a single path.
        for (std::size_t node = 0; node < graph.nodes.size(); ++node) {
            if (!graph.nodes[node].alive || graph.nodes[node].synthetic ||
                graph.nodes[node].incoming.size() != 1 || graph.nodes[node].outgoing.size() != 1) continue;
            const auto predecessor = *graph.nodes[node].incoming.begin();
            const auto successor = *graph.nodes[node].outgoing.begin();
            if (!graph.nodes[predecessor].alive || !graph.nodes[successor].alive ||
                graph.nodes[predecessor].synthetic || graph.nodes[successor].synthetic ||
                graph.nodes[predecessor].outgoing.size() != 1 ||
                graph.nodes[successor].incoming.size() != 1) continue;
            graph.remove_edge(predecessor, node);
            graph.remove_edge(node, successor);
            graph.nodes[predecessor].sequence += graph.nodes[node].sequence;
            graph.remove_node(node);
            graph.add_edge(predecessor, successor);
            ++result.linear_chain_merges;
            changed = true;
            break;
        }
        if (!changed) break;
    }

    result.final_nodes = active_sequence_node_count(graph);
    result.final_edges = active_sequence_edge_count(graph);

    // Reconstruct the language of the transformed graph.  This makes the
    // rewrite an actual path/topology operation rather than a counter-only
    // probe; callers can feed these paths directly to PairHMM or SW.
    std::set<std::string> transformed_paths;
    std::function<void(std::size_t, const std::string&, std::set<std::size_t>&)> walk =
        [&](const std::size_t node, const std::string& sequence, std::set<std::size_t>& visiting) {
            if (!graph.nodes[node].alive || !visiting.insert(node).second) return;
            auto next_sequence = sequence;
            if (!graph.nodes[node].synthetic) next_sequence += graph.nodes[node].sequence;
            if (graph.nodes[node].outgoing.empty()) {
                if (!next_sequence.empty()) transformed_paths.insert(next_sequence);
            } else {
                for (const auto next : graph.nodes[node].outgoing)
                    walk(next, next_sequence, visiting);
            }
            visiting.erase(node);
        };
    for (const auto source : source_nodes) {
        std::set<std::size_t> visiting;
        walk(source, {}, visiting);
    }
    // If a rewrite removed every boundary path (which is only possible for an
    // empty input), preserve the API's deterministic empty result.
    result.path_sequences.assign(transformed_paths.begin(), transformed_paths.end());
    return result;
}

KmerGraphResult build_kmer_graph_impl(const KmerGraphInput& input,
                                      KmerGraphOptions options) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (options.k < 1 || options.min_count == 0 || options.max_paths == 0)
        throw std::invalid_argument(
            "k-mer graph requires positive k/min_count/max_paths");
    if (input.offsets.empty() || input.offsets.front() != 0 ||
        input.offsets.back() != input.bases.size())
        throw std::invalid_argument("k-mer graph offsets are malformed");
    for (std::size_t i = 1; i < input.offsets.size(); ++i)
        if (input.offsets[i] < input.offsets[i - 1])
            throw std::invalid_argument("k-mer graph offsets are not monotonic");
    const auto read_count = input.offsets.size() - 1;
    if (!input.tids.empty() && input.tids.size() != read_count)
        throw std::invalid_argument("k-mer graph read tids do not match offsets");
    if (!input.flags.empty() && input.flags.size() != read_count)
        throw std::invalid_argument("k-mer graph read flags do not match offsets");
    if (!input.sample_ids.empty() && input.sample_ids.size() != read_count)
        throw std::invalid_argument("k-mer graph read sample ids do not match offsets");
    if (options.num_pruning_samples == 0 ||
        options.max_unpruned_variants == 0 ||
        !std::isfinite(options.initial_error_rate_for_pruning) ||
        options.initial_error_rate_for_pruning <= 0.0 ||
        options.initial_error_rate_for_pruning >= 1.0 ||
        !std::isfinite(options.pruning_log_odds_threshold) ||
        !std::isfinite(options.pruning_seeding_log_odds_threshold))
        throw std::invalid_argument("k-mer graph pruning options are invalid (thresholds must be finite)");
    if (!input.reference_offsets.empty()) {
        if (input.reference_offsets.front() != 0 ||
            input.reference_offsets.back() != input.reference_bases.size())
            throw std::invalid_argument("k-mer graph reference offsets are malformed");
        for (std::size_t i = 1; i < input.reference_offsets.size(); ++i)
            if (input.reference_offsets[i] < input.reference_offsets[i - 1])
                throw std::invalid_argument("k-mer graph reference offsets are not monotonic");
        const auto reference_segments = input.reference_offsets.size() - 1;
        if (!input.reference_tids.empty() && input.reference_tids.size() != reference_segments)
            throw std::invalid_argument("k-mer graph reference tids do not match offsets");
        if (!input.reference_starts.empty() && input.reference_starts.size() != reference_segments)
            throw std::invalid_argument("k-mer graph reference starts do not match offsets");
        if (!input.reference_ends.empty() && input.reference_ends.size() != reference_segments)
            throw std::invalid_argument("k-mer graph reference ends do not match offsets");
    } else if (!input.reference_bases.empty()) {
        throw std::invalid_argument("k-mer graph reference bases require offsets");
    }

    // ReadThreadingGraph deliberately distinguishes node identity from k-mer
    // text.  A k-mer repeated inside one pending sequence is excluded from
    // kmerToVertexMap, so it can reuse an existing outgoing suffix edge but
    // cannot globally merge into an equal-text vertex.  Keep that identity
    // explicit instead of using KmerKey as the compact node ID.
    struct EdgeKey { std::uint32_t from, to; };
    std::vector<EdgeKey> encoded_edges;
    std::vector<bool> edge_is_reference;
    std::vector<std::int32_t> edge_sample_ids;
    std::vector<KmerKey> encoded_nodes;
    std::vector<std::uint32_t> node_read_occurrences;
    std::vector<bool> node_is_read;
    std::vector<bool> node_is_reference;
    std::set<KmerKey> non_unique_keys;
    std::set<KmerKey> reference_non_unique_keys;
    struct ReferenceStart {
        std::uint32_t node = 0;
        std::uint32_t sink = 0;
        std::int32_t tid = -1;
        std::int32_t start = -1;
        std::int32_t end = -1;
    };
    std::vector<ReferenceStart> reference_starts;
    std::vector<std::vector<std::uint32_t>> reference_node_paths;
    std::size_t observed_kmers = 0;

    // ReadThreadingGraph tracks non-unique kmers per source sequence rather
    // than treating the same kmer in two independent reads as a graph
    // repeat. Recompute the compact key set here from each read/reference
    // segment; this also supplies a deterministic dynamic-k selection hint
    // to callers without coupling the traversal to Java's object graph.
    const auto collect_non_unique = [&](const std::vector<std::uint8_t>& bases,
                                        const std::uint32_t begin,
                                        const std::uint32_t end,
                                        const bool reference_sequence) {
        if (end - begin < options.k) return;
        std::map<KmerKey, std::uint32_t> local_counts;
        for (std::size_t offset = begin; offset + options.k <= end; ++offset) {
            bool valid = true;
            for (std::uint32_t i = 0; i < options.k; ++i) {
                if (encode(bases[offset + i]) >= 4) { valid = false; break; }
            }
            if (valid) ++local_counts[kmer_key(bases, offset, options.k)];
        }
        for (const auto& count : local_counts) {
            if (count.second <= 1) continue;
            non_unique_keys.insert(count.first);
            if (reference_sequence) reference_non_unique_keys.insert(count.first);
        }
    };
    for (std::size_t read = 0; read < read_count; ++read) {
        const auto begin = input.offsets[read];
        const auto end = input.offsets[read + 1];
        const bool reverse = !input.flags.empty() && (input.flags[read] & 0x10U) != 0;
        const auto sequence = orient_read(input.bases, begin, end, reverse);
        collect_non_unique(sequence, 0, static_cast<std::uint32_t>(sequence.size()), false);
        for (std::size_t offset = 0; offset + options.k <= sequence.size(); ++offset) {
            bool valid = true;
            for (std::uint32_t i = 0; i < options.k; ++i)
                if (encode(sequence[offset + i]) >= 4) { valid = false; break; }
            if (valid) ++observed_kmers;
        }
    }
    for (std::size_t contig = 0; contig + 1 < input.reference_offsets.size(); ++contig)
        collect_non_unique(input.reference_bases, input.reference_offsets[contig],
                           input.reference_offsets[contig + 1], true);

    KmerGraphResult result;
    result.input_reads = read_count;
    result.input_kmers = observed_kmers;
    result.kmer_size = options.k;
    result.execution_space = ExecSpace::name();
    result.used = true;
    result.reference_non_unique_kmers = reference_non_unique_keys.size();
    result.reference_kmer_rejected = !options.allow_non_unique_kmers_in_ref &&
        !reference_non_unique_keys.empty();
    result.non_unique_kmers.reserve(non_unique_keys.size());
    for (const auto& key : non_unique_keys)
        result.non_unique_kmers.push_back(key);
    // ReadThreadingAssembler does not construct a graph for this k when the
    // reference contains repeated kmers. Keep deterministic telemetry so Host
    // orchestration can retry with a larger k.
    if (result.reference_kmer_rejected) return result;

    // Build in the same order as AbstractReadThreadingGraph: reference first,
    // then each first-seen sample's reads.  The only global merge table holds
    // unique kmers; each vertex still has a Host-side immutable ID for the
    // Kokkos count/degree kernel.
    std::map<KmerKey, std::uint32_t> unique_kmer_to_node;
    std::vector<std::vector<std::uint32_t>> threaded_outgoing, threaded_incoming;
    std::optional<KmerKey> reference_source_key;
    const auto create_node = [&](const KmerKey& key) -> std::uint32_t {
        const auto node = static_cast<std::uint32_t>(encoded_nodes.size());
        encoded_nodes.push_back(key);
        node_read_occurrences.push_back(0U);
        node_is_read.push_back(false);
        node_is_reference.push_back(false);
        threaded_outgoing.emplace_back();
        threaded_incoming.emplace_back();
        if (non_unique_keys.find(key) == non_unique_keys.end())
            unique_kmer_to_node.emplace(key, node);
        return node;
    };
    const auto add_topology_edge = [&](const std::uint32_t from, const std::uint32_t to) {
        auto& outgoing = threaded_outgoing[from];
        if (std::find(outgoing.begin(), outgoing.end(), to) == outgoing.end()) {
            outgoing.push_back(to);
            threaded_incoming[to].push_back(from);
        }
    };
    const auto add_edge_observation = [&](const std::uint32_t from, const std::uint32_t to,
                                          const bool is_reference, const std::int32_t sample_id) {
        encoded_edges.push_back(EdgeKey{from, to});
        edge_is_reference.push_back(is_reference);
        edge_sample_ids.push_back(sample_id);
        add_topology_edge(from, to);
    };
    const auto mark_node = [&](const std::uint32_t node, const bool is_reference) {
        if (is_reference) {
            node_is_reference[node] = true;
        } else {
            node_is_read[node] = true;
            ++node_read_occurrences[node];
        }
    };
    const auto thread_sequence = [&](const std::vector<std::uint8_t>& sequence,
                                     const std::size_t begin, const std::size_t end,
                                     const bool is_reference, const std::int32_t sample_id)
        -> std::vector<std::uint32_t> {
        std::vector<std::uint32_t> path;
        if (end < begin || end - begin < options.k) return path;
        std::size_t start = begin;
        if (!is_reference) {
            bool found_start = false;
            // findStart() excludes the final k-mer: a read has to supply an
            // extension edge before it is threaded into the graph.
            for (std::size_t offset = begin; offset < end - options.k; ++offset) {
                const auto key = kmer_key(sequence, offset, options.k);
                if (non_unique_keys.find(key) == non_unique_keys.end()) {
                    start = offset;
                    found_start = true;
                    break;
                }
            }
            if (!found_start) return path;
        }
        const auto start_key = kmer_key(sequence, start, options.k);
        const auto start_iter = unique_kmer_to_node.find(start_key);
        const auto current_start = start_iter == unique_kmer_to_node.end()
            ? create_node(start_key) : start_iter->second;
        if (is_reference && !reference_source_key.has_value())
            reference_source_key = start_key;
        mark_node(current_start, is_reference);
        path.push_back(current_start);

        // AbstractReadThreadingGraph.threadSequence() does more than extend
        // forward from findStart(): with INCREASE_COUNTS_BACKWARDS enabled it
        // walks at most k-2 already-matched incoming edges and increments
        // their MultiSampleEdge multiplicities.  This is especially relevant
        // when initial repeated kmers force findStart() into the middle of a
        // read: GATK retains evidence immediately before the first usable
        // kmer, while a forward-only port undercounts the competing branch.
        // `increaseCountsThroughBranches` is false in the production
        // assembler, so source increments only through a vertex with exactly
        // one incoming edge.
        if (!is_reference && options.k >= 2U) {
            auto matched_vertex = current_start;
            for (std::int64_t offset = static_cast<std::int64_t>(options.k) - 2;
                 offset >= 0 && threaded_incoming[matched_vertex].size() == 1U;
                 --offset) {
                const auto previous = threaded_incoming[matched_vertex].front();
                if (encoded_nodes[previous].back() !=
                    start_key[static_cast<std::size_t>(offset)])
                    break;
                // The existing edge may be a reference edge.  Recording the
                // read observation separately preserves its reference flag
                // while contributing the read's total/sample multiplicity,
                // exactly as MultiSampleEdge.incMultiplicity().
                add_edge_observation(previous, matched_vertex, false, sample_id);
                matched_vertex = previous;
            }
        }
        auto current = current_start;
        for (std::size_t offset = start + 1U; offset + options.k <= end; ++offset) {
            const auto next_key = kmer_key(sequence, offset, options.k);
            const auto suffix = next_key.back();
            std::optional<std::uint32_t> next;
            // extendChainByOne() first reuses an existing outgoing edge by
            // suffix. This is the only legal merge for non-unique kmers.
            for (const auto candidate : threaded_outgoing[current]) {
                if (encoded_nodes[candidate].back() == suffix) {
                    next = candidate;
                    break;
                }
            }
            if (!next.has_value() && !is_reference &&
                (!reference_source_key.has_value() || next_key != *reference_source_key)) {
                const auto merge = unique_kmer_to_node.find(next_key);
                if (merge != unique_kmer_to_node.end()) next = merge->second;
            }
            if (!next.has_value()) next = create_node(next_key);
            mark_node(*next, is_reference);
            add_edge_observation(current, *next, is_reference, sample_id);
            current = *next;
            path.push_back(current);
        }
        return path;
    };
    for (std::size_t contig = 0; contig + 1U < input.reference_offsets.size(); ++contig) {
        const auto begin = input.reference_offsets[contig];
        const auto end = input.reference_offsets[contig + 1U];
        const auto path = thread_sequence(input.reference_bases, begin, end, true, -1);
        reference_node_paths.push_back(path);
        if (path.empty()) continue;
        reference_starts.push_back(ReferenceStart{
            path.front(), path.back(),
            input.reference_tids.empty() ? static_cast<std::int32_t>(contig)
                                          : input.reference_tids[contig],
            input.reference_starts.empty() ? 0 : input.reference_starts[contig],
            input.reference_ends.empty()
                ? (input.reference_starts.empty() ? static_cast<std::int32_t>(end - begin)
                                                   : input.reference_starts[contig] + static_cast<std::int32_t>(end - begin))
                : input.reference_ends[contig]});
    }
    std::vector<std::int32_t> sample_order;
    std::map<std::int32_t, std::vector<std::size_t>> reads_by_sample;
    for (std::size_t record = 0; record < read_count; ++record) {
        const auto sample = input.sample_ids.empty() ? 0 : input.sample_ids[record];
        if (reads_by_sample.emplace(sample, std::vector<std::size_t>{}).second)
            sample_order.push_back(sample);
        reads_by_sample[sample].push_back(record);
    }
    for (const auto sample : sample_order) {
        for (const auto record : reads_by_sample[sample]) {
            const auto begin = input.offsets[record];
            const auto end = input.offsets[record + 1U];
            const bool reverse = !input.flags.empty() && (input.flags[record] & 0x10U) != 0;
            const auto sequence = orient_read(input.bases, begin, end, reverse);
            std::size_t run_begin = 0;
            for (std::size_t offset = 0; offset <= sequence.size(); ++offset) {
                const bool valid = offset < sequence.size() && encode(sequence[offset]) < 4U;
                if (valid) continue;
                if (offset - run_begin >= options.k)
                    thread_sequence(sequence, run_begin, offset, false, sample);
                run_begin = offset + 1U;
            }
        }
    }
    if (encoded_nodes.empty()) return result;
    result.reference_nodes = static_cast<std::size_t>(std::count(
        node_is_reference.begin(), node_is_reference.end(), true));
    std::vector<std::uint32_t> from_ids, to_ids;
    from_ids.reserve(encoded_edges.size());
    to_ids.reserve(encoded_edges.size());
    for (const auto& edge : encoded_edges) {
        from_ids.push_back(edge.from);
        to_ids.push_back(edge.to);
    }
    // Collapse duplicate read/reference observations into weighted graph
    // edges for ChainPruner.  The Kokkos degree kernel intentionally keeps
    // the occurrence arrays unchanged; this Host representation is only the
    // variable-topology pruning view and is deterministic by sorted node IDs.
    std::vector<SeqGraphEdge> pruning_edges;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> pruning_edge_index;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::map<std::int32_t, std::uint32_t>> pruning_sample_counts;
    for (std::size_t edge = 0; edge < encoded_edges.size(); ++edge) {
        const auto key = std::make_pair(from_ids[edge], to_ids[edge]);
        auto iter = pruning_edge_index.find(key);
        if (iter == pruning_edge_index.end()) {
            iter = pruning_edge_index.emplace(key, pruning_edges.size()).first;
            pruning_edges.push_back(SeqGraphEdge{key.first, key.second, 0U, false, {}});
        }
        auto& compact = pruning_edges[iter->second];
        if (edge_is_reference[edge]) {
            compact.reference = true;
            // ReadThreadingAssembler adds the reference first through
            // `addSequence("ref", ..., 1, true)`.  Its MultiSampleEdge is
            // therefore born with multiplicity one, and any read that later
            // follows the same edge increments that total.  Pruning's
            // AdaptiveChainPruner scores use getMultiplicity(), not merely
            // read support.  Leaving this seed out made a one-read sibling
            // branch look spuriously competitive with the reference path.
            ++compact.multiplicity;
        } else {
            ++compact.multiplicity;
            ++pruning_sample_counts[key][edge_sample_ids[edge]];
        }
    }
    for (auto& compact : pruning_edges) {
        const auto key = std::make_pair(compact.from, compact.to);
        const auto counts_iter = pruning_sample_counts.find(key);
        if (counts_iter != pruning_sample_counts.end()) {
            for (const auto& [sample, count] : counts_iter->second) {
                (void)sample;
                compact.sample_multiplicities.push_back(count);
            }
        }
    }
    fastgatk::core::HostBatch host("kmer-graph-v1");
    host.records = encoded_edges.size();
    host.bytes = encoded_nodes.size() * sizeof(std::uint32_t) +
                 (from_ids.size() + to_ids.size()) * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("kmer-graph");
    plan.begin_prepare(host);
    Kokkos::View<std::uint32_t*, MemorySpace> from("kmer_from", from_ids.size());
    Kokkos::View<std::uint32_t*, MemorySpace> to("kmer_to", to_ids.size());
    Kokkos::View<std::uint32_t*, MemorySpace> counts("kmer_counts", encoded_nodes.size());
    Kokkos::View<std::uint32_t*, MemorySpace> indegree("kmer_indegree", encoded_nodes.size());
    Kokkos::View<std::uint32_t*, MemorySpace> outdegree("kmer_outdegree", encoded_nodes.size());
    auto host_from = Kokkos::create_mirror_view(from);
    auto host_to = Kokkos::create_mirror_view(to);
    for (std::size_t i = 0; i < from_ids.size(); ++i) { host_from(i) = from_ids[i]; host_to(i) = to_ids[i]; }
    Kokkos::deep_copy(from, host_from);
    Kokkos::deep_copy(to, host_to);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(encoded_edges.size());
    device.bind("from", from); device.bind("to", to); device.bind("counts", counts);
    device.bind("indegree", indegree); device.bind("outdegree", outdegree);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    Kokkos::deep_copy(counts, static_cast<std::uint32_t>(0));
    Kokkos::deep_copy(indegree, static_cast<std::uint32_t>(0));
    Kokkos::deep_copy(outdegree, static_cast<std::uint32_t>(0));
    Kokkos::parallel_for("kmer_graph_count", Kokkos::RangePolicy<ExecSpace>(0, encoded_edges.size()),
        KOKKOS_LAMBDA(const std::size_t edge) {
            Kokkos::atomic_add(&counts(from(edge)), 1U);
            Kokkos::atomic_add(&counts(to(edge)), 1U);
            Kokkos::atomic_add(&outdegree(from(edge)), 1U);
            Kokkos::atomic_add(&indegree(to(edge)), 1U);
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_counts = Kokkos::create_mirror_view(counts);
    auto host_indegree = Kokkos::create_mirror_view(indegree);
    auto host_outdegree = Kokkos::create_mirror_view(outdegree);
    Kokkos::deep_copy(host_counts, counts);
    Kokkos::deep_copy(host_indegree, indegree);
    Kokkos::deep_copy(host_outdegree, outdegree);
    // `counts` includes the reference scaffold.  Recover a conservative
    // read-only support value on Host by subtracting each reference edge's
    // endpoint incidence.  This is deterministic and keeps variable-length
    // path traversal out of the device kernel while still exposing enough
    // evidence for caller-side path pruning.
    std::vector<std::uint32_t> reference_incidence(encoded_nodes.size(), 0);
    for (std::size_t edge = 0; edge < encoded_edges.size(); ++edge) {
        if (!edge_is_reference[edge]) continue;
        ++reference_incidence[from_ids[edge]];
        ++reference_incidence[to_ids[edge]];
    }
    std::vector<std::uint32_t> read_support(encoded_nodes.size(), 0);
    for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
        if (node_read_occurrences[node] != 0U)
            read_support[node] = node_read_occurrences[node];
        else
            read_support[node] = host_counts(node) > reference_incidence[node]
                ? host_counts(node) - reference_incidence[node] : 0;
    }
    const auto& is_read = node_is_read;
    const auto& is_reference = node_is_reference;
    const auto pruning = prune_seqgraph_chains(
        encoded_nodes.size(), pruning_edges, options, &encoded_nodes);
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint8_t> pruning_keep_by_pair;
    for (std::size_t edge = 0; edge < pruning_edges.size(); ++edge)
        pruning_keep_by_pair[std::make_pair(pruning_edges[edge].from, pruning_edges[edge].to)] =
            edge < pruning.keep_edges.size() ? pruning.keep_edges[edge] : 1U;
    if (std::getenv("FASTGATK_DEBUG_EVENTMAP_PRUNING") != nullptr) {
        // Keep this before the compact graph is filtered.  EventMap tracing
        // normally shows only surviving edges, which cannot distinguish a
        // read that never threaded from one that ChainPruner discarded.
        for (std::size_t edge = 0; edge < pruning_edges.size(); ++edge) {
            const auto& compact = pruning_edges[edge];
            if (compact.multiplicity == 0U) continue;
            std::cerr << "[FASTGATK_EVENTMAP_PRUNING_EDGE] from="
                      << encoded_nodes[compact.from] << " to="
                      << encoded_nodes[compact.to] << " support="
                      << compact.multiplicity << " reference=" << compact.reference
                      << " keep=" << static_cast<unsigned>(
                          edge < pruning.keep_edges.size() ? pruning.keep_edges[edge] : 1U)
                      << '\n';
        }
    }
    std::vector<std::uint8_t> node_has_kept_pruning_edge(encoded_nodes.size(), 0U);
    for (std::size_t edge = 0; edge < pruning_edges.size(); ++edge) {
        if (edge >= pruning.keep_edges.size() || pruning.keep_edges[edge] == 0U) continue;
        node_has_kept_pruning_edge[pruning_edges[edge].from] = 1U;
        node_has_kept_pruning_edge[pruning_edges[edge].to] = 1U;
    }
    std::vector<bool> keep(encoded_nodes.size(), false);
    for (std::size_t i = 0; i < encoded_nodes.size(); ++i) {
        // Reference nodes form the scaffold even when no read has enough
        // support. Read-only nodes survive only when at least one weighted
        // edge remains after ChainPruner (and still obey min_count).
        keep[i] = is_reference[i] ||
            (host_counts(i) >= options.min_count && node_has_kept_pruning_edge[i] != 0U);
        if (is_read[i] && !is_reference[i] && !keep[i] && options.use_adaptive_pruning)
            ++result.adaptive_pruned_nodes;
        if (is_read[i] && !keep[i]) ++result.pruned_nodes;
        if (keep[i] && (host_indegree(i) > 1 || host_outdegree(i) > 1)) ++result.branching_nodes;
        if (keep[i]) ++result.nodes;
    }
    std::vector<std::vector<std::uint32_t>> adjacency;
    std::vector<std::vector<std::uint32_t>> directed_adjacency;
    if (!reference_starts.empty()) adjacency.resize(encoded_nodes.size());
    if (!reference_starts.empty()) directed_adjacency.resize(encoded_nodes.size());
    for (std::size_t i = 0; i < encoded_edges.size(); ++i) {
        if (!keep[from_ids[i]] || !keep[to_ids[i]]) continue;
        const auto pruning_iter = pruning_keep_by_pair.find({from_ids[i], to_ids[i]});
        if (pruning_iter != pruning_keep_by_pair.end() && pruning_iter->second == 0U) continue;
        ++result.edges;
        if (edge_is_reference[i]) ++result.reference_edges;
        if (!adjacency.empty()) {
            adjacency[from_ids[i]].push_back(to_ids[i]);
            adjacency[to_ids[i]].push_back(from_ids[i]);
            directed_adjacency[from_ids[i]].push_back(to_ids[i]);
        }
    }
    if (!reference_starts.empty()) {
        // ReadThreadingAssembler passes SmithWatermanAlignmentConstants.STANDARD_NGS
        // to both dangling-end recovery operations.  These are deliberately
        // distinct from the 10/-15/-30/-5 read-to-haplotype parameters used
        // by the public SW API: the dangling CIGAR controls a Host topology
        // rewrite, so a different scoring regime can fabricate a graph edge.
        const SmithWatermanParameters dangling_end_sw_parameters{25, -50, -110, -6};
        // GATK keeps every ChainPruner-retained component while it runs
        // recoverDanglingTails() and recoverDanglingHeads().  A branch that
        // currently has no reference connection can acquire one through one
        // of those topology rewrites.  The source-to-sink intersection used
        // by removePathsNotConnectedToRef() is deliberately deferred until
        // after both recovery passes below.  In particular, an undirected
        // reachability pre-filter here would erase a valid dangling head
        // before its Smith-Waterman CIGAR has a chance to reconnect it.
        std::vector<bool> reachable(encoded_nodes.size(), false);

        // Deterministic source-to-sink traversal over the reference-connected
        // graph. GATK's K-best finder materializes only paths that reach a
        // reference sink. max_paths bounds the candidate population; an
        // optional positive max_depth is solely a guard on an unrejoined
        // alternate branch and never promotes a truncated prefix to a
        // haplotype.  In particular, the shared reference prefix/suffix can
        // be much longer than this native guard.
        for (auto& neighbors : adjacency) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
        for (auto& neighbors : directed_adjacency) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
        // Bit 0 = reference scaffold edge, bit 1 = read-only edge.  A pair
        // present in both is reference-connected and must not be labelled an
        // alternate merely because the read graph also observed it.
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint8_t> edge_kind_by_pair;
        // Keep edge multiplicity separate from node support. GATK's
        // MultiSampleEdge uses this evidence when choosing among equally
        // valid outgoing paths; a stable support-descending tie-break avoids
        // input-order dependent haplotype ordering while preserving the
        // alternate-before-reference policy below.
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> read_edge_support;
        // The dangling-end walk uses MultiSampleEdge's pruning multiplicity,
        // which is distinct from aggregate support used for K-best ranking.
        std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t>
            pruning_edge_support;
        for (const auto& [edge, sample_counts] : pruning_sample_counts) {
            std::vector<std::uint32_t> counts;
            counts.reserve(sample_counts.size());
            for (const auto& [sample, count] : sample_counts) {
                (void)sample;
                counts.push_back(count);
            }
            std::sort(counts.begin(), counts.end(), std::greater<>());
            std::uint32_t support = 0U;
            if (options.num_pruning_samples == 0U) {
                for (const auto count : counts) support += count;
            } else if (counts.size() >= options.num_pruning_samples) {
                support = counts[options.num_pruning_samples - 1U];
            }
            if (support != 0U) pruning_edge_support.emplace(edge, support);
        }
        for (std::size_t edge = 0; edge < encoded_edges.size(); ++edge) {
            if (!keep[from_ids[edge]] || !keep[to_ids[edge]]) continue;
            const auto pruning_iter = pruning_keep_by_pair.find({from_ids[edge], to_ids[edge]});
            if (pruning_iter != pruning_keep_by_pair.end() && pruning_iter->second == 0U) continue;
            const auto key = std::make_pair(from_ids[edge], to_ids[edge]);
            auto& kind = edge_kind_by_pair[key];
            kind |= edge_is_reference[edge] ? 1U : 2U;
            if (!edge_is_reference[edge]) ++read_edge_support[key];
        }
        // BaseGraph.isReferenceNode() is an edge property: a vertex belongs
        // to the reference graph if any incident edge is a reference edge.
        // `node_is_reference` records how the vertex was first threaded and
        // is deliberately not equivalent after a read merges into a unique
        // reference k-mer.  Dangling recovery must use the source predicate.
        const auto is_reference_node = [&](const std::uint32_t node) {
            for (const auto& [edge, kind] : edge_kind_by_pair) {
                if ((kind & 1U) != 0U && (edge.first == node || edge.second == node))
                    return true;
            }
            return false;
        };
        if (std::getenv("FASTGATK_DEBUG_EVENTMAP") != nullptr) {
            for (const auto& [edge, kind] : edge_kind_by_pair) {
                if ((kind & 2U) == 0U) continue;
                const auto support = read_edge_support.find(edge);
                std::cerr << "[FASTGATK_EVENTMAP_GRAPH_EDGE] from="
                          << encoded_nodes[edge.first] << " to="
                          << encoded_nodes[edge.second] << " kind="
                          << static_cast<unsigned>(kind) << " support="
                          << (support == read_edge_support.end() ? 0U : support->second)
                          << '\n';
            }
        }
        // Keep the graph topology canonical before running GATK's dangling
        // head walk: graph input keeps one edge record per read observation,
        // but a source with several observations of the *same* next k-mer is
        // still a single-outdegree chain in ReadThreadingGraph.
        for (auto& neighbors : directed_adjacency) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
        // ReadThreadingAssembler recovers tails before heads, while the graph
        // still represents the original dangling topology.  The recovery is
        // a Host graph rewrite, not a post-enumeration sequence patch: find
        // the source-equivalent LCA/reference walks, obtain every affine
        // score from the Kokkos batch, validate the Host CIGAR traceback, and
        // then add the one alternate-to-reference edge selected by
        // AbstractReadThreadingGraph.mergeDanglingTail().  Subsequent
        // SeqGraph conversion and K-best traversal therefore see the same
        // reconnected topology as GATK.
        {
            // `allowNonUniqueKmersInRef` is consulted only when
            // ReadThreadingAssembler decides whether this k-mer attempt may
            // be built.  Once built, GATK always runs recoverDanglingTails()
            // and recoverDanglingHeads() over its actual graph topology.
            // Skipping tail recovery here under the opt-in flag made a
            // repeated reference k-mer silently change the cleanup state
            // machine, rather than merely admitting the graph attempt.
            std::vector<std::vector<std::uint32_t>> incoming_adjacency(encoded_nodes.size());
            for (std::size_t from = 0; from < directed_adjacency.size(); ++from)
                for (const auto to : directed_adjacency[from])
                    incoming_adjacency[to].push_back(static_cast<std::uint32_t>(from));
            for (auto& neighbors : incoming_adjacency) {
                std::sort(neighbors.begin(), neighbors.end());
                neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
            }
            const auto has_incoming_reference_edge = [&](const std::uint32_t node) {
                if (node >= incoming_adjacency.size()) return false;
                for (const auto previous : incoming_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({previous, node});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return true;
                }
                return false;
            };
            // BaseGraph.isRefSource()/isRefSink() are edge predicates, not
            // construction-history flags.  A k-mer visited by a read can
            // still be a reference endpoint, and dangling recovery must not
            // use that endpoint as an alternate LCA or terminal.
            const auto is_reference_source = [&](const std::uint32_t node) {
                if (node >= incoming_adjacency.size() || node >= directed_adjacency.size())
                    return false;
                for (const auto previous : incoming_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({previous, node});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return false;
                }
                for (const auto next : directed_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({node, next});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return true;
                }
                return encoded_nodes.size() == 1U;
            };
            const auto is_reference_sink = [&](const std::uint32_t node) {
                if (node >= incoming_adjacency.size() || node >= directed_adjacency.size())
                    return false;
                for (const auto next : directed_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({node, next});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return false;
                }
                for (const auto previous : incoming_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({previous, node});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return true;
                }
                return encoded_nodes.size() == 1U;
            };
            // BaseGraph's dangling-tail reference walk is seeded with the
            // *edge* returned by getHeaviestIncomingEdge(altPath[1]).
            // MultiSampleEdge.getMultiplicity() includes the one reference
            // observation used to thread the scaffold, whereas
            // read_edge_support deliberately contains reads only.  Keep the
            // two notions separate: pruning below uses the per-sample value,
            // while this tie-break reproduces the graph edge multiplicity.
            const auto graph_edge_multiplicity = [&](const std::uint32_t from,
                                                     const std::uint32_t to) {
                const auto key = std::make_pair(from, to);
                const auto support = read_edge_support.find(key);
                std::uint32_t multiplicity =
                    support == read_edge_support.end() ? 0U : support->second;
                const auto kind = edge_kind_by_pair.find(key);
                if (kind != edge_kind_by_pair.end() && (kind->second & 1U) != 0U)
                    ++multiplicity;
                return std::max(1U, multiplicity);
            };
            const auto choose_heaviest_incoming = [&](const std::uint32_t node)
                -> std::optional<std::uint32_t> {
                if (node >= incoming_adjacency.size() || incoming_adjacency[node].empty())
                    return std::nullopt;
                auto best = incoming_adjacency[node].front();
                std::uint32_t best_support = 0U;
                for (const auto candidate : incoming_adjacency[node]) {
                    const auto count = graph_edge_multiplicity(candidate, node);
                    if (count > best_support || (count == best_support && candidate < best)) {
                        best = candidate;
                        best_support = count;
                    }
                }
                return best;
            };
            const auto pruning_multiplicity = [&](const std::uint32_t from,
                                                   const std::uint32_t to) {
                const auto support = pruning_edge_support.find({from, to});
                // Reference edges are constructed with multiplicity one in
                // GATK. The normal tail walk stops at the LCA before using
                // one, but recover-all traversal can inspect shared scaffold
                // topology as well.
                return support == pruning_edge_support.end() ? 1U : support->second;
            };
            // MutectReadThreadingAssemblerArgumentCollection delegates
            // adaptive pruning to AdaptiveChainPruner and constructs its
            // ReadThreadingAssembler with pruneFactor == 0.  That factor is
            // subsequently passed verbatim to recoverDanglingTails/Heads;
            // it must not be replaced with the CLI min-pruning floor after
            // the adaptive pruner has retained a heterogeneous chain.
            const auto dangling_prune_factor = options.use_adaptive_pruning
                ? 0U : options.min_pruning;
            const auto find_tail_path = [&](const std::uint32_t terminal)
                -> std::optional<std::vector<std::uint32_t>> {
                if (terminal >= encoded_nodes.size() || is_reference_sink(terminal) ||
                    !keep[terminal] ||
                    terminal >= directed_adjacency.size() ||
                    !directed_adjacency[terminal].empty())
                    return std::nullopt;
                std::vector<std::uint32_t> reversed;
                std::set<std::uint32_t> discarded;
                auto current = terminal;
                while (true) {
                    const auto indegree = incoming_adjacency[current].size();
                    const auto outdegree = directed_adjacency[current].size();
                    const bool done = options.recover_all_dangling_branches
                        ? (has_incoming_reference_edge(current) || indegree == 0U)
                        : (indegree != 1U || outdegree >= 2U);
                    if (done) {
                        const bool returns_path = options.recover_all_dangling_branches
                            ? (outdegree > 1U && has_incoming_reference_edge(current))
                            : outdegree > 1U;
                        if (!returns_path) return std::nullopt;
                        reversed.push_back(current);
                        break;
                    }
                    const auto previous = options.recover_all_dangling_branches
                        ? choose_heaviest_incoming(current)
                        : std::optional<std::uint32_t>(incoming_adjacency[current].front());
                    if (!previous.has_value()) return std::nullopt;
                    // GATK findPath() clears only the downstream suffix at
                    // a sub-threshold per-sample edge, then continues toward
                    // the LCA instead of discarding the full chain.
                    if (dangling_prune_factor != 0U &&
                        pruning_multiplicity(*previous, current) < dangling_prune_factor) {
                        discarded.insert(reversed.begin(), reversed.end());
                        reversed.clear();
                    } else {
                        reversed.push_back(current);
                    }
                    current = *previous;
                    if (discarded.find(current) != discarded.end() ||
                        std::find(reversed.begin(), reversed.end(), current) != reversed.end())
                        return std::nullopt;
                }
                std::reverse(reversed.begin(), reversed.end());
                const auto minimum_length = std::max<std::uint32_t>(
                    1U, options.min_dangling_branch_length);
                if (reversed.size() < static_cast<std::size_t>(minimum_length) + 1U)
                    return std::nullopt;
                return reversed;
            };
            const auto reference_path_from_lca = [&](const std::vector<std::uint32_t>& alternate)
                -> std::vector<std::uint32_t> {
                std::vector<std::uint32_t> reference;
                if (alternate.size() < 2U) return reference;
                // AbstractReadThreadingGraph does not blacklist the LCA's
                // first alternate edge by position.  It passes
                // getHeaviestIncomingEdge(altPath[1]) to
                // getReferencePath(), which may be a different incoming edge
                // when the first alternate vertex is revisited by another
                // read branch.  Treat the blacklist as an exact directed
                // edge, preserving BaseGraph's subsequent singleton walk.
                const auto blacklisted_to = alternate[1];
                const auto heaviest = choose_heaviest_incoming(blacklisted_to);
                const auto blacklisted_from = heaviest.value_or(alternate[0]);
                std::set<std::uint32_t> visited;
                auto current = alternate.front();
                while (current < directed_adjacency.size() && visited.insert(current).second) {
                    reference.push_back(current);
                    std::optional<std::uint32_t> next;
                    for (const auto candidate : directed_adjacency[current]) {
                        const auto edge = edge_kind_by_pair.find({current, candidate});
                        if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U) {
                            next = candidate;
                            break;
                        }
                    }
                    if (!next.has_value()) {
                        std::optional<std::uint32_t> singleton;
                        for (const auto candidate : directed_adjacency[current]) {
                            if (current == blacklisted_from && candidate == blacklisted_to) continue;
                            if (singleton.has_value()) { singleton.reset(); break; }
                            singleton = candidate;
                        }
                        next = singleton;
                    }
                    if (!next.has_value()) break;
                    current = *next;
                }
                return reference;
            };
            struct DanglingTailTopologyRequest {
                std::vector<std::uint32_t> alternate_path;
                std::vector<std::uint32_t> reference_path;
                std::string alternate_bases;
                std::string reference_bases;
            };
            std::vector<DanglingTailTopologyRequest> tail_requests;
            std::vector<SmithWatermanRequest> tail_score_requests;
            for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
                const auto alternate = find_tail_path(static_cast<std::uint32_t>(node));
                if (!alternate.has_value()) continue;
                if (std::getenv("FASTGATK_DEBUG_DANGLING_RECOVERY") != nullptr) {
                    std::cerr << "[FASTGATK_DANGLING_TAIL_PATH] terminal="
                              << encoded_nodes[node] << " path=";
                    for (const auto path_node : *alternate)
                        std::cerr << encoded_nodes[path_node] << ',';
                    std::cerr << " ref_source="
                              << (is_reference_source(alternate->front()) ? 1 : 0) << '\n';
                }
                if (is_reference_source(alternate->front())) continue;
                const auto reference = reference_path_from_lca(*alternate);
                if (reference.empty()) continue;
                std::string alternate_bases;
                std::string reference_bases;
                alternate_bases.reserve(alternate->size());
                reference_bases.reserve(reference.size());
                for (const auto path_node : *alternate)
                    alternate_bases.push_back(encoded_nodes[path_node].back());
                for (const auto path_node : reference)
                    reference_bases.push_back(encoded_nodes[path_node].back());
                if (alternate_bases.empty() || reference_bases.empty()) continue;
                tail_requests.push_back(DanglingTailTopologyRequest{
                    *alternate, reference, std::move(alternate_bases), std::move(reference_bases)});
                const auto& request = tail_requests.back();
                tail_score_requests.push_back(SmithWatermanRequest{
                    std::vector<std::uint8_t>(request.alternate_bases.begin(),
                                              request.alternate_bases.end()),
                    std::vector<std::uint8_t>(request.reference_bases.begin(),
                                              request.reference_bases.end())});
            }
            const auto tail_scores = smith_waterman_score_kokkos(
                tail_score_requests, dangling_end_sw_parameters,
                SmithWatermanOverhangStrategy::LeadingIndel);
            result.prepare_seconds += tail_scores.prepare_seconds;
            result.seconds += tail_scores.seconds;
            for (std::size_t request_index = 0; request_index < tail_requests.size(); ++request_index) {
                const auto& request = tail_requests[request_index];
                const auto alignment = smith_waterman_align_reference(
                    reinterpret_cast<const std::uint8_t*>(request.alternate_bases.data()),
                    request.alternate_bases.size(),
                    reinterpret_cast<const std::uint8_t*>(request.reference_bases.data()),
                    request.reference_bases.size(), dangling_end_sw_parameters,
                    SmithWatermanOverhangStrategy::LeadingIndel);
                if (request_index >= tail_scores.scores.size() ||
                    alignment.score != tail_scores.scores[request_index])
                    throw std::runtime_error(
                        "NUMERICAL_CONTRACT_FAILURE: dangling-tail SW score/traceback mismatch");
                const auto recovered = recover_dangling_tail_sequence(
                    request.alternate_bases, request.reference_bases, alignment,
                    options.min_dangling_matching_bases);
                if (std::getenv("FASTGATK_DEBUG_DANGLING_RECOVERY") != nullptr) {
                    std::cerr << "[FASTGATK_DANGLING_TAIL] alt="
                              << request.alternate_bases << " ref="
                              << request.reference_bases << " cigar=" << alignment.cigar
                              << " recovered=" << (recovered.has_value() ? 1 : 0);
                    if (recovered.has_value())
                        std::cerr << " alt_index=" << recovered->alternate_merge_index
                                  << " ref_index=" << recovered->reference_merge_index;
                    std::cerr << '\n';
                }
                if (!recovered.has_value() ||
                    recovered->alternate_merge_index >= request.alternate_path.size() ||
                    recovered->reference_merge_index >= request.reference_path.size())
                    continue;
                const auto from = request.alternate_path[recovered->alternate_merge_index];
                const auto to = request.reference_path[recovered->reference_merge_index];
                if (from == to) continue;
                const auto edge_key = std::make_pair(from, to);
                if (edge_kind_by_pair.find(edge_key) != edge_kind_by_pair.end()) continue;
                if (std::getenv("FASTGATK_DEBUG_DANGLING_RECOVERY") != nullptr) {
                    std::cerr << "[FASTGATK_DANGLING_TAIL_EDGE] from="
                              << encoded_nodes[from] << " to=" << encoded_nodes[to] << '\n';
                }
                directed_adjacency[from].push_back(to);
                edge_kind_by_pair[edge_key] = 2U;
                // AbstractReadThreadingGraph.mergeDanglingTail() inserts a
                // new MultiSampleEdge(false, 1), not an edge weighted by the
                // source vertex's aggregate coverage. Its multiplicity is
                // observable to GraphBasedKBestHaplotypeFinder.
                read_edge_support[edge_key] = 1U;
                pruning_edge_support[edge_key] = 1U;
                ++result.dangling_recovered_paths;
                result.dangling_recovered_bases += request.alternate_path.size() - 1U;
            }
            for (auto& neighbors : directed_adjacency) {
                std::sort(neighbors.begin(), neighbors.end());
                neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
            }
        }
        // AbstractReadThreadingGraph.recoverDanglingHeads() runs before the
        // SeqGraph conversion.  A recoverable head is an indegree-zero
        // non-reference chain that reaches a reference node; GATK scores its
        // upward (reversed) path with LEADING_INDEL then adds one graph edge
        // from the selected upstream reference vertex into the non-reference
        // chain.  When that target lies inside GATK's expanded source k-mer,
        // the source implementation first materializes synthetic k-mers from
        // the reference prefix; retain that path language for later Host path
        // materialization.  Every affine score remains in the Kokkos SW batch
        // and the Host traceback is checked against the device score.
        {
            // As for dangling tails above, GATK consults
            // allowNonUniqueKmersInRef only while deciding whether this
            // k-mer graph attempt is admissible.  Once built,
            // ReadThreadingAssembler always executes both
            // recoverDanglingTails() and recoverDanglingHeads().  A repeated
            // reference k-mer must not turn off the head topology rewrite.
            std::vector<std::uint32_t> directed_indegree_before_head(encoded_nodes.size(), 0U);
            for (const auto& neighbors : directed_adjacency)
                for (const auto next : neighbors) ++directed_indegree_before_head[next];
            // BaseGraph.isRefSource() is edge-based, not a label assigned
            // when the reference was first threaded. A non-reference vertex
            // with no incoming edge but an outgoing reference edge is the
            // reference source and must never be recovered as a dangling
            // head. This distinction matters for duplicated non-unique
            // k-mers that meet the reference at their first suffix edge.
            const auto is_reference_source = [&](const std::uint32_t node) {
                if (node >= directed_adjacency.size()) return false;
                for (const auto next : directed_adjacency[node]) {
                    const auto edge = edge_kind_by_pair.find({node, next});
                    if (edge != edge_kind_by_pair.end() && (edge->second & 1U) != 0U)
                        return true;
                }
                return false;
            };
            // Match BaseGraph.isRefSink(): a vertex is the terminal
            // reference vertex when a reference edge enters it but no
            // reference edge leaves it.  (The degenerate one-vertex graph
            // is also both source and sink.)  This must be tested before
            // constructing the dangling-head SW request.
            const auto is_reference_sink = [&](const std::uint32_t node) {
                if (node >= directed_adjacency.size()) return false;
                bool has_incoming_reference = false;
                for (const auto& [edge, kind] : edge_kind_by_pair) {
                    if ((kind & 1U) == 0U) continue;
                    if (edge.first == node) return false;
                    if (edge.second == node) has_incoming_reference = true;
                }
                return has_incoming_reference || encoded_nodes.size() == 1U;
            };
            struct DanglingHeadRequest {
                std::vector<std::uint32_t> path_nodes;
                std::size_t reference_segment = 0;
                std::size_t reference_hcd_start = 0;
                std::string alternate_upwards;
                std::string reference_upwards;
            };
            const auto minimum_matching_bases = options.min_dangling_matching_bases < 0
                ? std::size_t{0}
                : static_cast<std::size_t>(options.min_dangling_matching_bases);
            // The same ReadThreadingAssembler pruneFactor is passed to both
            // recoverDanglingTails() and recoverDanglingHeads().  In
            // particular, adaptive pruning has already decided which chains
            // survive, so its dangling-path walk uses a factor of zero.
            const auto dangling_prune_factor = options.use_adaptive_pruning
                ? 0U : options.min_pruning;
            const auto graph_edge_multiplicity = [&](const std::uint32_t from,
                                                     const std::uint32_t to) {
                const auto key = std::make_pair(from, to);
                const auto support = read_edge_support.find(key);
                std::uint32_t multiplicity =
                    support == read_edge_support.end() ? 0U : support->second;
                const auto kind = edge_kind_by_pair.find(key);
                if (kind != edge_kind_by_pair.end() && (kind->second & 1U) != 0U)
                    ++multiplicity;
                return std::max(1U, multiplicity);
            };
            const auto pruning_multiplicity = [&](const std::uint32_t from,
                                                   const std::uint32_t to) {
                const auto support = pruning_edge_support.find({from, to});
                return support == pruning_edge_support.end() ? 1U : support->second;
            };
            std::vector<DanglingHeadRequest> dangling_head_requests;
            std::vector<SmithWatermanRequest> dangling_head_score_requests;
            for (std::size_t source = 0; source < encoded_nodes.size(); ++source) {
                if (!keep[source] ||
                    directed_indegree_before_head[source] != 0U ||
                    is_reference_source(static_cast<std::uint32_t>(source)))
                    continue;
                std::vector<std::uint32_t> head_nodes{static_cast<std::uint32_t>(source)};
                std::set<std::uint32_t> visited{static_cast<std::uint32_t>(source)};
                auto current = static_cast<std::uint32_t>(source);
                bool valid_path = true;
                while (!is_reference_node(current)) {
                    const auto& next_nodes = directed_adjacency[current];
                    if (next_nodes.empty()) { valid_path = false; break; }
                    if (!options.recover_all_dangling_branches && next_nodes.size() != 1U) {
                        valid_path = false;
                        break;
                    }
                    auto next = next_nodes.front();
                    if (options.recover_all_dangling_branches && next_nodes.size() > 1U) {
                        std::uint32_t best_support = 0U;
                        for (const auto candidate : next_nodes) {
                            const auto support = graph_edge_multiplicity(current, candidate);
                            if (support > best_support ||
                                (support == best_support && candidate < next)) {
                                next = candidate;
                                best_support = support;
                            }
                        }
                    }
                    if (!visited.insert(next).second) { valid_path = false; break; }
                    // AbstractReadThreadingGraph.findPath() removes the
                    // source-side prefix after a sub-threshold edge, then
                    // continues toward the highest common descendant.  Do
                    // not let a retained heterogeneous chain turn that
                    // discarded prefix into a dangling-head CIGAR input.
                    if (dangling_prune_factor != 0U &&
                        pruning_multiplicity(current, next) < dangling_prune_factor)
                        head_nodes.clear();
                    head_nodes.push_back(next);
                    current = next;
                }
                // AbstractReadThreadingGraph rejects a dangling head whose
                // highest common descendant is itself the reference sink:
                // there is no downstream reference vertex available for the
                // reconnection edge.  In particular, do not synthesize a
                // new alternate extension past that sink; GATK returns
                // without attempting the legacy merge in this case.
                if (!valid_path || !is_reference_node(current) ||
                    is_reference_sink(current) ||
                    head_nodes.size() < static_cast<std::size_t>(options.min_dangling_branch_length) + 1U)
                    continue;

                std::size_t reference_segment = input.reference_offsets.size();
                std::size_t reference_hcd_start = 0;
                for (std::size_t segment = 0; segment + 1 < input.reference_offsets.size(); ++segment) {
                    const auto begin = input.reference_offsets[segment];
                    const auto end = input.reference_offsets[segment + 1];
                    if (end < begin || end - begin < options.k) continue;
                    if (segment >= reference_node_paths.size()) continue;
                    const auto& reference_nodes = reference_node_paths[segment];
                    const auto position = std::find(reference_nodes.begin(), reference_nodes.end(), current);
                    if (position != reference_nodes.end()) {
                        reference_segment = segment;
                        reference_hcd_start = static_cast<std::size_t>(
                            std::distance(reference_nodes.begin(), position));
                        break;
                    }
                }
                if (reference_segment >= input.reference_offsets.size()) continue;

                std::string head_sequence = encoded_nodes[head_nodes.front()];
                for (std::size_t node = 1; node < head_nodes.size(); ++node)
                    head_sequence.push_back(encoded_nodes[head_nodes[node]].back());
                const auto reference_begin = input.reference_offsets[reference_segment];
                const auto reference_end = input.reference_offsets[reference_segment + 1U];
                if (reference_end < reference_begin ||
                    reference_hcd_start + options.k > reference_end - reference_begin)
                    continue;
                std::string reference_upwards(
                    reinterpret_cast<const char*>(input.reference_bases.data() + reference_begin),
                    reference_hcd_start + options.k);
                std::reverse(head_sequence.begin(), head_sequence.end());
                std::reverse(reference_upwards.begin(), reference_upwards.end());
                dangling_head_requests.push_back(DanglingHeadRequest{
                    std::move(head_nodes), reference_segment, reference_hcd_start,
                    std::move(head_sequence), std::move(reference_upwards)});
                const auto& request = dangling_head_requests.back();
                dangling_head_score_requests.push_back(SmithWatermanRequest{
                    std::vector<std::uint8_t>(request.alternate_upwards.begin(),
                                              request.alternate_upwards.end()),
                    std::vector<std::uint8_t>(request.reference_upwards.begin(),
                                              request.reference_upwards.end())});
            }
            const auto dangling_head_scores = smith_waterman_score_kokkos(
                dangling_head_score_requests, dangling_end_sw_parameters,
                SmithWatermanOverhangStrategy::LeadingIndel);
            result.prepare_seconds += dangling_head_scores.prepare_seconds;
            result.seconds += dangling_head_scores.seconds;
            for (std::size_t request_index = 0; request_index < dangling_head_requests.size();
                 ++request_index) {
                const auto& request = dangling_head_requests[request_index];
                const auto alignment = smith_waterman_align_reference(
                    reinterpret_cast<const std::uint8_t*>(request.alternate_upwards.data()),
                    request.alternate_upwards.size(),
                    reinterpret_cast<const std::uint8_t*>(request.reference_upwards.data()),
                    request.reference_upwards.size(), dangling_end_sw_parameters,
                    SmithWatermanOverhangStrategy::LeadingIndel);
                if (request_index >= dangling_head_scores.scores.size() ||
                    alignment.score != dangling_head_scores.scores[request_index])
                    throw std::runtime_error(
                        "NUMERICAL_CONTRACT_FAILURE: dangling-head SW score/traceback mismatch");
                const auto merge = options.min_dangling_matching_bases < 0
                    ? recover_dangling_head_legacy_merge(
                        request.alternate_upwards, request.reference_upwards, alignment,
                        options.k)
                    : recover_dangling_head_merge(
                        request.alternate_upwards, request.reference_upwards, alignment,
                        minimum_matching_bases);
                if (std::getenv("FASTGATK_DEBUG_DANGLING_RECOVERY") != nullptr &&
                    !request.path_nodes.empty()) {
                    std::cerr << "[FASTGATK_DANGLING_HEAD] source="
                              << encoded_nodes[request.path_nodes.front()]
                              << " hcd=" << encoded_nodes[request.path_nodes.back()]
                              << " cigar=" << alignment.cigar
                              << " merge=" << (merge.has_value() ? 1 : 0);
                    if (merge.has_value())
                        std::cerr << " ref_index=" << merge->reference_index
                                  << " alt_index=" << merge->alternate_index;
                    std::cerr << '\n';
                }
                if (!merge.has_value() ||
                    merge->reference_index >= request.reference_hcd_start)
                    continue;
                if (merge->alternate_index >= request.path_nodes.size()) {
                    // AbstractReadThreadingGraph.extendDanglingPathAgainstReference:
                    // the CIGAR can point into the source vertex's expanded
                    // k-mer rather than at a materialized graph node.  GATK
                    // rewrites the *graph* in this case: it detaches the
                    // dangling source edge, creates a chain of synthetic
                    // k-mer vertices, and connects the selected reference
                    // node to that chain.  Do that Host topology mutation
                    // here before K-best traversal.  The degree/count stage
                    // has already run on immutable Kokkos IDs; these are
                    // recovery-only Host vertices and do not replace any
                    // numerical kernel.
                    std::vector<CigarElement> elements;
                    if (!parse_cigar(alignment.cigar, elements)) continue;
                    if (!elements.empty() && elements.back().operation == 'D')
                        elements.pop_back();
                    const auto reference_length = cigar_reference_length(elements);
                    const auto read_length = cigar_read_length(elements);
                    if (reference_length == 0 || read_length == 0) continue;
                    const auto nodes_to_extend = merge->alternate_index -
                        request.path_nodes.size() + 2U;
                    if (nodes_to_extend > options.k) continue;
                    const auto reference_offset = static_cast<std::int64_t>(reference_length) -
                        static_cast<std::int64_t>(read_length);
                    const auto reference_node_index = static_cast<std::int64_t>(
                        request.path_nodes.size() - 1U) + reference_offset +
                        static_cast<std::int64_t>(nodes_to_extend);
                    if (reference_node_index < 0 ||
                        static_cast<std::size_t>(reference_node_index) > request.reference_hcd_start)
                        continue;
                    const auto reference_source_start = request.reference_hcd_start -
                        static_cast<std::size_t>(reference_node_index);
                    const auto reference_anchor_start = request.reference_hcd_start -
                        merge->reference_index - 1U;
                    if (std::getenv("FASTGATK_DEBUG_EVENTMAP") != nullptr) {
                        std::cerr << "[FASTGATK_EVENTMAP_HEAD_EXTENSION]"
                                  << " hcd=" << request.reference_hcd_start
                                  << " nodes=" << request.path_nodes.size()
                                  << " ref_index=" << merge->reference_index
                                  << " alt_index=" << merge->alternate_index
                                  << " cigar=" << alignment.cigar
                                  << " extend=" << nodes_to_extend
                                  << " ref_node=" << reference_node_index
                                  << " ref_source=" << reference_source_start
                                  << " anchor=" << reference_anchor_start << '\n';
                    }
                    const auto reference_begin = input.reference_offsets[request.reference_segment];
                    const auto reference_end = input.reference_offsets[request.reference_segment + 1U];
                    if (reference_end < reference_begin ||
                        reference_source_start + options.k > reference_end - reference_begin ||
                        reference_anchor_start + options.k > reference_end - reference_begin)
                        continue;
                    if (request.reference_segment >= reference_node_paths.size() ||
                        reference_anchor_start >= reference_node_paths[request.reference_segment].size())
                        continue;
                    const auto reference_node =
                        reference_node_paths[request.reference_segment][reference_anchor_start];
                    const auto dangling_source = request.path_nodes.front();
                    if (request.path_nodes.size() < 2U || dangling_source >= directed_adjacency.size())
                        continue;
                    const auto original_target = request.path_nodes[1U];
                    const auto old_edge = std::make_pair(dangling_source, original_target);
                    const auto old_support = read_edge_support.find(old_edge);
                    const auto source_support = old_support == read_edge_support.end()
                        ? 1U : std::max(1U, old_support->second);
                    const auto old_kind = edge_kind_by_pair.find(old_edge);
                    if (old_kind == edge_kind_by_pair.end()) continue;

                    std::string sequence_to_extend(
                        reinterpret_cast<const char*>(input.reference_bases.data() +
                                                       reference_begin + reference_source_start),
                        nodes_to_extend);
                    sequence_to_extend += encoded_nodes[dangling_source];
                    if (sequence_to_extend.size() < options.k + nodes_to_extend) continue;

                    // removeEdge(danglingSource, prevV)
                    auto& source_outgoing = directed_adjacency[dangling_source];
                    source_outgoing.erase(std::remove(source_outgoing.begin(), source_outgoing.end(),
                                                      original_target), source_outgoing.end());
                    edge_kind_by_pair.erase(old_edge);
                    read_edge_support.erase(old_edge);

                    // extendDanglingPathAgainstReference() iterates from N
                    // down to one.  Preserve that order so the CIGAR-derived
                    // danglingPath index below selects the same new vertex.
                    std::vector<std::uint32_t> synthetic_nodes;
                    synthetic_nodes.reserve(nodes_to_extend);
                    auto previous = original_target;
                    for (std::size_t offset = nodes_to_extend; offset > 0U; --offset) {
                        const auto node = static_cast<std::uint32_t>(encoded_nodes.size());
                        encoded_nodes.push_back(sequence_to_extend.substr(offset, options.k));
                        node_read_occurrences.push_back(0U);
                        node_is_read.push_back(false);
                        node_is_reference.push_back(false);
                        read_support.push_back(source_support);
                        keep.push_back(true);
                        reachable.push_back(true);
                        directed_adjacency.emplace_back();
                        directed_adjacency[node].push_back(previous);
                        const auto edge = std::make_pair(node, previous);
                        edge_kind_by_pair[edge] = 2U;
                        read_edge_support[edge] = source_support;
                        synthetic_nodes.push_back(node);
                        previous = node;
                    }
                    // GATK removes the original dangling source before it
                    // appends these nodes. Its first synthetic vertex thus
                    // occupies `path_nodes.size() - 1` in danglingPath.
                    // `merge->alternate_index` is still expressed in the
                    // pre-removal path coordinate, so account for that one
                    // position here. Without it the bridge skips a k-mer and
                    // creates a false one-base deletion at the join.
                    const auto synthetic_index = merge->alternate_index -
                        request.path_nodes.size() + 1U;
                    if (synthetic_index >= synthetic_nodes.size())
                        continue;
                    const auto alternate_node = synthetic_nodes[synthetic_index];
                    const auto recovery_edge = std::make_pair(reference_node, alternate_node);
                    if (edge_kind_by_pair.find(recovery_edge) != edge_kind_by_pair.end())
                        continue;
                    directed_adjacency[reference_node].push_back(alternate_node);
                    edge_kind_by_pair[recovery_edge] = 2U;
                    read_edge_support[recovery_edge] = 1U;
                    ++result.dangling_recovered_paths;
                    result.dangling_recovered_bases += request.path_nodes.size() - 1U + nodes_to_extend;
                    continue;
                }
                const auto reference_kmer_start = request.reference_hcd_start -
                    merge->reference_index - 1U;
                if (request.reference_segment >= reference_node_paths.size() ||
                    reference_kmer_start >= reference_node_paths[request.reference_segment].size())
                    continue;
                const auto reference_node =
                    reference_node_paths[request.reference_segment][reference_kmer_start];
                const auto alternate_node = request.path_nodes[
                    request.path_nodes.size() - 1U - merge->alternate_index];
                if (reference_node == alternate_node || !is_reference_node(reference_node)) continue;
                const auto edge_key = std::make_pair(reference_node, alternate_node);
                if (edge_kind_by_pair.find(edge_key) != edge_kind_by_pair.end()) continue;
                if (std::getenv("FASTGATK_DEBUG_DANGLING_RECOVERY") != nullptr)
                    std::cerr << "[FASTGATK_DANGLING_HEAD_EDGE] from="
                              << encoded_nodes[reference_node] << " to="
                              << encoded_nodes[alternate_node] << '\n';
                directed_adjacency[reference_node].push_back(alternate_node);
                edge_kind_by_pair[edge_key] = 2U;
                // mergeDanglingHead(), like tail recovery, adds a fresh
                // non-reference edge with multiplicity one.
                read_edge_support[edge_key] = 1U;
                ++result.dangling_recovered_paths;
                result.dangling_recovered_bases += request.path_nodes.size() - 1U;
            }
        }
        // ReadThreadingAssembler.removePathsNotConnectedToRef() is a
        // *directed* source-to-sink intersection, not an undirected
        // connected-component test.  It runs after dangling-tail/head graph
        // rewrites, so a recovered tail gains reachability to the reference
        // sink and a recovered head gains reachability from the reference
        // source.  Conversely, a long read branch that only leaves the
        // scaffold is removed even though it is weakly connected to it.
        if (std::getenv("FASTGATK_DEBUG_EVENTMAP_TOPOLOGY_PRE_CLEAN") != nullptr) {
            for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
                if (!keep[node]) continue;
                std::cerr << "[FASTGATK_EVENTMAP_PRE_CLEAN_NODE] seq="
                          << encoded_nodes[node] << " reference="
                          << (is_reference[node] ? 1 : 0) << '\n';
            }
            for (const auto& [edge, kind] : edge_kind_by_pair) {
                std::cerr << "[FASTGATK_EVENTMAP_PRE_CLEAN_EDGE] from="
                          << encoded_nodes[edge.first] << " to="
                          << encoded_nodes[edge.second] << " kind="
                          << static_cast<unsigned>(kind) << " support="
                          << read_edge_support[edge] << '\n';
            }
        }
        {
            std::vector<std::vector<std::uint32_t>> reverse_adjacency(encoded_nodes.size());
            for (std::size_t from = 0; from < directed_adjacency.size(); ++from) {
                if (!keep[from]) continue;
                for (const auto to : directed_adjacency[from])
                    if (to < keep.size() && keep[to]) reverse_adjacency[to].push_back(
                        static_cast<std::uint32_t>(from));
            }
            std::vector<std::uint32_t> reference_sources;
            std::vector<std::uint32_t> reference_sinks;
            for (const auto& reference_start : reference_starts) {
                reference_sources.push_back(reference_start.node);
                reference_sinks.push_back(reference_start.sink);
            }
            std::sort(reference_sources.begin(), reference_sources.end());
            reference_sources.erase(std::unique(reference_sources.begin(), reference_sources.end()),
                                    reference_sources.end());
            std::sort(reference_sinks.begin(), reference_sinks.end());
            reference_sinks.erase(std::unique(reference_sinks.begin(), reference_sinks.end()),
                                  reference_sinks.end());
            std::vector<bool> from_reference_source(encoded_nodes.size(), false);
            std::vector<bool> to_reference_sink(encoded_nodes.size(), false);
            std::deque<std::uint32_t> pending;
            for (const auto source : reference_sources) {
                if (source >= keep.size() || !keep[source] || from_reference_source[source]) continue;
                from_reference_source[source] = true;
                pending.push_back(source);
            }
            while (!pending.empty()) {
                const auto node = pending.front();
                pending.pop_front();
                for (const auto next : directed_adjacency[node]) {
                    if (next >= keep.size() || !keep[next] || from_reference_source[next]) continue;
                    from_reference_source[next] = true;
                    pending.push_back(next);
                }
            }
            for (const auto sink : reference_sinks) {
                if (sink >= keep.size() || !keep[sink] || to_reference_sink[sink]) continue;
                to_reference_sink[sink] = true;
                pending.push_back(sink);
            }
            while (!pending.empty()) {
                const auto node = pending.front();
                pending.pop_front();
                for (const auto previous : reverse_adjacency[node]) {
                    if (previous >= keep.size() || !keep[previous] || to_reference_sink[previous]) continue;
                    to_reference_sink[previous] = true;
                    pending.push_back(previous);
                }
            }
            reachable.assign(encoded_nodes.size(), false);
            result.reference_connected_nodes = 0;
            result.dangling_nodes = 0;
            for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
                const bool on_reference_path = keep[node] && from_reference_source[node] &&
                    to_reference_sink[node];
                if (keep[node] && !on_reference_path) {
                    keep[node] = false;
                    if (is_read[node]) ++result.pruned_nodes;
                }
                reachable[node] = on_reference_path;
                if (is_read[node]) {
                    if (on_reference_path) ++result.reference_connected_nodes;
                    else ++result.dangling_nodes;
                }
            }
            for (auto& neighbors : directed_adjacency) {
                neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
                    [&](const auto node) { return node >= keep.size() || !keep[node]; }),
                    neighbors.end());
            }
            for (auto edge = edge_kind_by_pair.begin(); edge != edge_kind_by_pair.end();) {
                if (edge->first.first >= keep.size() || edge->first.second >= keep.size() ||
                    !keep[edge->first.first] || !keep[edge->first.second]) {
                    read_edge_support.erase(edge->first);
                    edge = edge_kind_by_pair.erase(edge);
                } else ++edge;
            }
            result.nodes = 0;
            result.edges = 0;
            result.reference_edges = 0;
            for (const auto value : keep) if (value) ++result.nodes;
            for (const auto& edge : edge_kind_by_pair) {
                ++result.edges;
                if ((edge.second & 1U) != 0U) ++result.reference_edges;
            }
        }
        if (std::getenv("FASTGATK_DEBUG_EVENTMAP_TOPOLOGY") != nullptr) {
            for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
                if (!keep[node]) continue;
                std::cerr << "[FASTGATK_EVENTMAP_TOPOLOGY_NODE] seq="
                          << encoded_nodes[node] << " reference="
                          << (is_reference[node] ? 1 : 0) << '\n';
            }
            for (const auto& [edge, kind] : edge_kind_by_pair) {
                std::cerr << "[FASTGATK_EVENTMAP_TOPOLOGY_EDGE] from="
                          << encoded_nodes[edge.first] << " to="
                          << encoded_nodes[edge.second] << " kind="
                          << static_cast<unsigned>(kind) << " support="
                          << read_edge_support[edge] << '\n';
            }
        }
        // Mirror SeqGraph.zipLinearChains on the compact directed topology.
        // A chain may be merged only when each adjacent vertex has a single
        // incoming/outgoing edge and both vertices have the same
        // reference/read classification.  The deterministic start rule is
        // the one used by GATK: a chain starts at an indegree-zero/junction
        // vertex, or immediately after a branching predecessor.  Cycles
        // deliberately have no start and are left intact for the cycle
        // retry policy above.
        result.seqgraph_nodes = 0;
        for (const auto value : keep) if (value) ++result.seqgraph_nodes;
        result.seqgraph_edges = edge_kind_by_pair.size();
        std::vector<std::uint32_t> compact_indegree(encoded_nodes.size(), 0);
        std::vector<std::uint32_t> compact_outdegree(encoded_nodes.size(), 0);
        std::vector<std::uint32_t> compact_predecessor(encoded_nodes.size(),
                                                       std::numeric_limits<std::uint32_t>::max());
        for (const auto& edge_entry : edge_kind_by_pair) {
            const auto from = edge_entry.first.first;
            const auto to = edge_entry.first.second;
            ++compact_outdegree[from];
            ++compact_indegree[to];
            compact_predecessor[to] = from;
        }
        std::vector<bool> compressed_node(encoded_nodes.size(), false);
        for (std::size_t node = 0; node < encoded_nodes.size(); ++node) {
            if (!keep[node] || compact_outdegree[node] != 1 ||
                (compact_indegree[node] == 1 &&
                 compact_predecessor[node] < encoded_nodes.size() &&
                 compact_outdegree[compact_predecessor[node]] <= 1))
                continue;
            std::vector<std::uint32_t> chain;
            auto current = static_cast<std::uint32_t>(node);
            chain.push_back(current);
            while (compact_outdegree[current] == 1) {
                const auto iter = edge_kind_by_pair.lower_bound(
                    std::make_pair(current, std::uint32_t{0}));
                if (iter == edge_kind_by_pair.end() || iter->first.first != current) break;
                const auto next = iter->first.second;
                if (next == current || compact_indegree[next] != 1 ||
                    is_reference[next] != is_reference[current]) break;
                if (compressed_node[next]) break;
                chain.push_back(next);
                current = next;
            }
            if (chain.size() <= 1) continue;
            for (const auto member : chain) compressed_node[member] = true;
            const auto removed = chain.size() - 1;
            result.seqgraph_nodes = result.seqgraph_nodes >= removed
                ? result.seqgraph_nodes - removed : 0;
            result.seqgraph_edges = result.seqgraph_edges >= removed
                ? result.seqgraph_edges - removed : 0;
        }
        // Detect alternate-containing directed cycles independently of path
        // materialization. A conventional three-colour DFS can miss a cycle
        // when the alternate edge enters a node that was already explored by
        // a different branch, so check reachability for every read-supported
        // edge instead: u->v participates in a cycle iff u is reachable from
        // v. Only a read-only (non-reference) edge can establish such a
        // cycle: an edge shared by the reference and reads is reference
        // topology, not an alternate branch. This is Host-only and bounded
        // by the local AssemblyRegion.
        std::vector<std::uint32_t> cycle_seen_stamp(encoded_nodes.size(), 0);
        std::uint32_t cycle_stamp = 0;
        for (const auto& edge_entry : edge_kind_by_pair) {
            if (edge_entry.second != 2U) continue;
            const auto from_node = edge_entry.first.first;
            const auto to_node = edge_entry.first.second;
            if (!keep[from_node] || !keep[to_node] ||
                !reachable[from_node] || !reachable[to_node]) continue;
            if (++cycle_stamp == 0U) {
                std::fill(cycle_seen_stamp.begin(), cycle_seen_stamp.end(), 0U);
                cycle_stamp = 1U;
            }
            std::deque<std::uint32_t> pending;
            cycle_seen_stamp[to_node] = cycle_stamp;
            pending.push_back(to_node);
            while (!pending.empty() && cycle_seen_stamp[from_node] != cycle_stamp) {
                const auto node = pending.front();
                pending.pop_front();
                for (const auto next : directed_adjacency[node]) {
                    if (!keep[next] || !reachable[next] ||
                        cycle_seen_stamp[next] == cycle_stamp) continue;
                    cycle_seen_stamp[next] = cycle_stamp;
                    pending.push_back(next);
                }
            }
            if (cycle_seen_stamp[from_node] == cycle_stamp) {
                if (std::getenv("FASTGATK_DEBUG_EVENTMAP_TOPOLOGY") != nullptr) {
                    std::cerr << "[FASTGATK_EVENTMAP_NON_REFERENCE_CYCLE] from="
                              << encoded_nodes[from_node] << " to="
                              << encoded_nodes[to_node] << '\n';
                }
                result.has_non_reference_cycles = true;
                break;
            }
        }
        // Explore alternate edges before the reference self-cycle.  In a
        // homopolymer scaffold the latter is a valid de Bruijn cycle and can
        // otherwise consume the bounded path budget before a supported branch
        // is ever visited.
        for (std::size_t node = 0; node < directed_adjacency.size(); ++node) {
            auto& neighbors = directed_adjacency[node];
            std::sort(neighbors.begin(), neighbors.end(), [&](const auto left, const auto right) {
                const auto left_iter = edge_kind_by_pair.find(std::make_pair(
                    static_cast<std::uint32_t>(node), left));
                const auto right_iter = edge_kind_by_pair.find(std::make_pair(
                    static_cast<std::uint32_t>(node), right));
                const bool left_alt = left_iter != edge_kind_by_pair.end() && left_iter->second == 2U;
                const bool right_alt = right_iter != edge_kind_by_pair.end() && right_iter->second == 2U;
                if (left_alt != right_alt) return left_alt > right_alt;
                const auto left_support = read_edge_support.find(std::make_pair(
                    static_cast<std::uint32_t>(node), left));
                const auto right_support = read_edge_support.find(std::make_pair(
                    static_cast<std::uint32_t>(node), right));
                const auto left_count = left_support == read_edge_support.end() ? 0U : left_support->second;
                const auto right_count = right_support == read_edge_support.end() ? 0U : right_support->second;
                if (left_count != right_count) return left_count > right_count;
                return left < right;
            });
        }
        std::vector<std::uint32_t> directed_indegree(encoded_nodes.size(), 0);
        for (const auto& neighbors : directed_adjacency)
            for (const auto next : neighbors) ++directed_indegree[next];
        std::vector<std::uint32_t> starts;
        std::vector<std::uint32_t> sinks;
        std::vector<std::int32_t> start_tids;
        std::vector<std::int32_t> start_positions;
        std::vector<std::int32_t> start_ends;
        std::set<std::uint32_t> seen_starts;
        // Prefer the first k-mer of each reference window.  Unlike a global
        // graph indegree scan this preserves the coordinate frame of local
        // AssemblyRegions even when overlapping reads add an incoming edge.
        for (const auto& reference_start : reference_starts) {
            const auto node = reference_start.node;
            const auto sink = reference_start.sink;
            if (!keep[node] || !reachable[node] || !keep[sink] || !reachable[sink] ||
                (node != sink && directed_adjacency[node].empty()) ||
                !seen_starts.insert(node).second) continue;
            starts.push_back(node);
            sinks.push_back(sink);
            start_tids.push_back(reference_start.tid);
            start_positions.push_back(reference_start.start);
            start_ends.push_back(reference_start.end);
        }
        // A reference-free graph (or a malformed local scaffold) retains the
        // historical deterministic indegree-zero fallback, but marks its
        // paths as unanchored so callers never mistake them for coordinates.
        if (starts.empty()) {
            for (std::size_t i = 0; i < encoded_nodes.size(); ++i) {
                if (keep[i] && reachable[i] && !directed_adjacency[i].empty() &&
                    directed_indegree[i] == 0) {
                    starts.push_back(static_cast<std::uint32_t>(i));
                    sinks.push_back(std::numeric_limits<std::uint32_t>::max());
                    start_tids.push_back(-1);
                    start_positions.push_back(-1);
                    start_ends.push_back(-1);
                }
            }
        }
        if (starts.empty()) {
            for (std::size_t i = 0; i < encoded_nodes.size(); ++i)
                if (keep[i] && is_reference[i]) {
                    starts.push_back(static_cast<std::uint32_t>(i));
                    sinks.push_back(std::numeric_limits<std::uint32_t>::max());
                    start_tids.push_back(-1);
                    start_positions.push_back(-1);
                    start_ends.push_back(-1);
                    break;
                }
        }
        const auto max_paths = options.max_paths;
        const auto max_depth = options.max_depth;
        // GraphBasedKBestHaplotypeFinder uses a priority queue ordered by
        // the accumulated log10(edgeMultiplicity / totalOutgoingMultiplicity).
        // A DFS is not equivalent once the 128-haplotype cap is reached: it
        // retains whichever lexical branch happens to be visited first and
        // can discard a better-supported EventMap haplotype.
        struct TraversalPath {
            std::uint32_t node = 0;
            std::uint32_t sink = std::numeric_limits<std::uint32_t>::max();
            std::size_t depth = 0;
            std::string sequence;
            std::int32_t tid = -1;
            std::int32_t start = -1;
            std::int32_t end = -1;
            std::uint32_t support = 0;
            bool has_alt = false;
            bool rejoined_reference = false;
            std::size_t alt_begin = std::numeric_limits<std::size_t>::max();
            std::size_t alt_end = 0;
            std::uint32_t alt_edge_from = std::numeric_limits<std::uint32_t>::max();
            std::uint32_t alt_edge_to = std::numeric_limits<std::uint32_t>::max();
            double score = 0.0;
            std::vector<std::uint32_t> nodes;
        };
        struct LowerPriority {
            bool operator()(const TraversalPath& left, const TraversalPath& right) const {
                if (left.score != right.score) return left.score < right.score;
                if (left.sequence != right.sequence) return left.sequence < right.sequence;
                return left.node > right.node;
            }
        };
        const auto edge_multiplicity = [&](const std::uint32_t from, const std::uint32_t to) {
            const auto key = std::make_pair(from, to);
            const auto support = read_edge_support.find(key);
            std::uint32_t multiplicity = support == read_edge_support.end() ? 0U : support->second;
            const auto kind = edge_kind_by_pair.find(key);
            if (kind != edge_kind_by_pair.end() && (kind->second & 1U) != 0U)
                ++multiplicity;  // ReadThreadingGraph seeds reference edges at one observation.
            return std::max(1U, multiplicity);
        };
        const auto record_haplotype = [&](const TraversalPath& path) {
            if (std::getenv("FASTGATK_DEBUG_EVENTMAP_KBEST") != nullptr) {
                std::cerr << "[FASTGATK_EVENTMAP_KBEST_PATH] index="
                          << result.haplotype_path_count
                          << " score=" << std::setprecision(17) << path.score
                          << " support=" << path.support
                          << " alternate=" << (path.has_alt ? 1 : 0)
                          << " sequence=" << path.sequence << '\n';
            }
            ++result.haplotype_path_count;
            result.haplotype_path_sequences.push_back(path.sequence);
            result.haplotype_path_tids.push_back(path.tid);
            result.haplotype_path_starts.push_back(path.start);
            result.haplotype_path_ends.push_back(path.end);
            result.haplotype_path_support.push_back(path.support);
            result.haplotype_path_scores.push_back(path.score);
            result.haplotype_path_has_non_reference_edge.push_back(path.has_alt ? 1U : 0U);
            result.haplotype_path_alt_read_starts.push_back(
                path.alt_begin == std::numeric_limits<std::size_t>::max()
                    ? 0U : static_cast<std::uint32_t>(std::min<std::size_t>(path.alt_begin, UINT32_MAX)));
            result.haplotype_path_alt_read_ends.push_back(
                path.alt_end == 0 ? 0U : static_cast<std::uint32_t>(std::min<std::size_t>(path.alt_end, UINT32_MAX)));
        };
        std::priority_queue<TraversalPath, std::vector<TraversalPath>, LowerPriority> pending_paths;
        std::vector<std::size_t> vertex_visits(encoded_nodes.size(), 0);
        for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
            const auto start = starts[start_index];
            pending_paths.push(TraversalPath{start,
                start_index < sinks.size() ? sinks[start_index] : std::numeric_limits<std::uint32_t>::max(),
                0, encoded_nodes[start], start_tids[start_index], start_positions[start_index],
                start_ends[start_index], read_support[start], false, false,
                std::numeric_limits<std::size_t>::max(), 0,
                std::numeric_limits<std::uint32_t>::max(), std::numeric_limits<std::uint32_t>::max(),
                0.0, {start}});
        }
        while (result.haplotype_path_count < max_paths && !pending_paths.empty()) {
            auto path = pending_paths.top();
            pending_paths.pop();
            const auto branch_length = path.alt_end > path.alt_begin ? path.alt_end - path.alt_begin : 0U;
            const bool short_dangling = path.has_alt && !path.rejoined_reference &&
                options.min_dangling_branch_length > 0 && branch_length < options.min_dangling_branch_length;
            if (path.sink != std::numeric_limits<std::uint32_t>::max() && path.node == path.sink) {
                if (short_dangling) {
                    ++result.dangling_branch_paths;
                    result.dangling_branch_bases += branch_length;
                } else {
                    record_haplotype(path);
                }
                continue;
            }
            const bool alternate_branch_exhausted =
                max_depth != 0 && path.has_alt && !path.rejoined_reference &&
                path.alt_begin != std::numeric_limits<std::size_t>::max() &&
                path.sequence.size() >= path.alt_begin &&
                path.sequence.size() - path.alt_begin >= max_depth;
            if (alternate_branch_exhausted || directed_adjacency[path.node].empty()) {
                if (short_dangling) {
                    ++result.dangling_branch_paths;
                    result.dangling_branch_bases += branch_length;
                }
                continue;
            }
            // This is the same per-vertex expansion cap used by GATK's
            // GraphBasedKBestHaplotypeFinder.  It bounds queue growth without
            // substituting traversal order for likelihood.
            if (vertex_visits[path.node]++ >= max_paths) continue;
            std::uint64_t total_outgoing = 0;
            for (const auto next : directed_adjacency[path.node])
                if (keep[next] && reachable[next])
                    total_outgoing += edge_multiplicity(path.node, next);
            if (total_outgoing == 0) continue;
            for (const auto next : directed_adjacency[path.node]) {
                if (!keep[next] || !reachable[next]) continue;
                const auto edge = edge_kind_by_pair.find(std::make_pair(path.node, next));
                const bool edge_is_alt = edge != edge_kind_by_pair.end() && edge->second == 2U;
                if (edge_is_alt && path.has_alt &&
                    path.node == path.alt_edge_from && next == path.alt_edge_to) continue;
                if (std::find(path.nodes.begin(), path.nodes.end(), next) != path.nodes.end()) {
                    if (path.has_alt || edge_is_alt) {
                        result.has_non_reference_cycles = true;
                        if (std::getenv("FASTGATK_DEBUG_KMER_CYCLES") != nullptr) {
                            std::cerr << "[FASTGATK_KMER_CYCLE]"
                                      << " from=" << encoded_nodes[path.node]
                                      << " to=" << encoded_nodes[next]
                                      << " path_has_alt=" << (path.has_alt ? 1 : 0)
                                      << " edge_is_alt=" << (edge_is_alt ? 1 : 0)
                                      << " depth=" << path.depth << '\n';
                        }
                    }
                    continue;
                }
                auto extended = path;
                extended.node = next;
                ++extended.depth;
                extended.sequence.push_back(encoded_nodes[next].back());
                extended.nodes.push_back(next);
                const auto multiplicity = edge_multiplicity(path.node, next);
                extended.score += std::log10(static_cast<double>(multiplicity)) -
                    std::log10(static_cast<double>(total_outgoing));
                extended.support = std::max({path.support, read_support[next], multiplicity});
                extended.has_alt = path.has_alt || edge_is_alt;
                extended.rejoined_reference = path.rejoined_reference ||
                    ((path.has_alt || edge_is_alt) && is_reference[next]);
                if (edge_is_alt) {
                    if (extended.alt_begin == std::numeric_limits<std::size_t>::max())
                        extended.alt_begin = path.sequence.size();
                    extended.alt_end = extended.sequence.size();
                    if (!path.has_alt) {
                        extended.alt_edge_from = path.node;
                        extended.alt_edge_to = next;
                    }
                }
                pending_paths.push(std::move(extended));
            }
        }
        if (std::getenv("FASTGATK_DEBUG_EVENTMAP") != nullptr &&
            (options.linked_de_bruijn_graph || options.disable_seqgraph_simplification)) {
            std::cerr << "[FASTGATK_EVENTMAP_KBEST_RAW] paths="
                      << result.haplotype_path_count << '\n';
        }
        // The normal ReadThreadingAssembler mode calls
        // GraphBasedKBestHaplotypeFinder on cleanupSeqGraph()'s final
        // SeqGraph, not on the raw k-mer topology.  Preserve the latter for
        // linked-de-Bruijn mode only; replace the bounded normal-mode list
        // with the source-equivalent, weighted SeqGraph traversal.
        if (!options.linked_de_bruijn_graph && !options.disable_seqgraph_simplification) {
            const auto seqgraph_kbest = exact_seqgraph_kbest(
                encoded_nodes, read_support, edge_kind_by_pair, read_edge_support, keep, max_paths);
            result.haplotype_path_count = 0;
            result.haplotype_path_sequences.clear();
            result.haplotype_path_tids.clear();
            result.haplotype_path_starts.clear();
            result.haplotype_path_ends.clear();
            result.haplotype_path_support.clear();
            result.haplotype_path_scores.clear();
            result.haplotype_path_has_non_reference_edge.clear();
            result.haplotype_path_alt_read_starts.clear();
            result.haplotype_path_alt_read_ends.clear();
            result.seqgraph_nodes = seqgraph_kbest.nodes;
            result.seqgraph_edges = seqgraph_kbest.edges;
            const auto tid = start_tids.empty() ? -1 : start_tids.front();
            const auto start_position = start_positions.empty() ? -1 : start_positions.front();
            const auto end_position = start_ends.empty() ? -1 : start_ends.front();
            for (const auto& path : seqgraph_kbest.paths) {
                if (std::getenv("FASTGATK_DEBUG_EVENTMAP_KBEST") != nullptr) {
                    std::cerr << "[FASTGATK_EVENTMAP_SEQGRAPH_KBEST_PATH] index="
                              << result.haplotype_path_count
                              << " score=" << std::setprecision(17) << path.score
                              << " support=" << path.support
                              << " alternate=" << (path.alternate ? 1 : 0)
                              << " sequence=" << path.sequence << '\n';
                }
                ++result.haplotype_path_count;
                result.haplotype_path_sequences.push_back(path.sequence);
                result.haplotype_path_tids.push_back(tid);
                result.haplotype_path_starts.push_back(start_position);
                result.haplotype_path_ends.push_back(end_position);
                result.haplotype_path_support.push_back(path.support);
                result.haplotype_path_scores.push_back(path.score);
                result.haplotype_path_has_non_reference_edge.push_back(path.alternate ? 1U : 0U);
                result.haplotype_path_alt_read_starts.push_back(
                    path.alt_begin == std::numeric_limits<std::size_t>::max()
                        ? 0U : static_cast<std::uint32_t>(std::min<std::size_t>(path.alt_begin, UINT32_MAX)));
                result.haplotype_path_alt_read_ends.push_back(
                    path.alt_end == 0U ? 0U : static_cast<std::uint32_t>(
                        std::min<std::size_t>(path.alt_end, UINT32_MAX)));
            }
        }
        // JunctionTreeLinkedDeBruijnGraph deliberately does not zip the
        // k-mer graph into a SeqGraph.  Its K-best finder can therefore leave
        // a pivotal edge uncovered when the bounded queue reaches
        // max_paths.  GATK's default linked mode recovers such an edge by
        // stapling it onto the best completed path that reaches the edge's
        // source.  Recreate that operation here from the raw topology so the
        // linked and regular modes share the same Kokkos-built graph while
        // retaining the linked-mode recovery toggle.
        if (options.linked_de_bruijn_graph &&
            !options.disable_artificial_haplotype_recovery &&
            result.haplotype_path_count < max_paths) {
            const auto graph_node_for_kmer = [&](const std::string& sequence,
                                                  const std::size_t offset) {
                if (offset + options.k > sequence.size()) return std::numeric_limits<std::uint32_t>::max();
                const auto key = sequence.substr(offset, options.k);
                const auto unique = unique_kmer_to_node.find(key);
                if (unique != unique_kmer_to_node.end()) return unique->second;
                const auto node = std::find(encoded_nodes.begin(), encoded_nodes.end(), key);
                return node == encoded_nodes.end()
                    ? std::numeric_limits<std::uint32_t>::max()
                    : static_cast<std::uint32_t>(node - encoded_nodes.begin());
            };
            std::vector<std::vector<std::uint32_t>> path_nodes(
                result.haplotype_path_sequences.size());
            std::set<std::pair<std::uint32_t, std::uint32_t>> covered_edges;
            std::vector<std::size_t> best_path_for_node(encoded_nodes.size(),
                                                        std::numeric_limits<std::size_t>::max());
            for (std::size_t path = 0; path < result.haplotype_path_sequences.size(); ++path) {
                const auto& sequence = result.haplotype_path_sequences[path];
                auto& nodes_for_path = path_nodes[path];
                if (sequence.size() < options.k) continue;
                nodes_for_path.reserve(sequence.size() - options.k + 1U);
                for (std::size_t offset = 0; offset + options.k <= sequence.size(); ++offset) {
                    const auto node = graph_node_for_kmer(sequence, offset);
                    if (node == std::numeric_limits<std::uint32_t>::max()) continue;
                    nodes_for_path.push_back(node);
                }
                for (std::size_t index = 1; index < nodes_for_path.size(); ++index)
                    covered_edges.emplace(nodes_for_path[index - 1], nodes_for_path[index]);
                const auto support = path < result.haplotype_path_support.size()
                    ? result.haplotype_path_support[path] : 0U;
                for (const auto node : nodes_for_path) {
                    const auto previous = best_path_for_node[node];
                    const auto previous_support = previous == std::numeric_limits<std::size_t>::max()
                        ? 0U : result.haplotype_path_support[previous];
                    if (previous == std::numeric_limits<std::size_t>::max() || support > previous_support ||
                        (support == previous_support && path < previous))
                        best_path_for_node[node] = path;
                }
            }
            std::vector<std::pair<std::uint32_t, std::uint32_t>> uncovered_edges;
            for (std::size_t from = 0; from < directed_adjacency.size(); ++from) {
                for (const auto to : directed_adjacency[from]) {
                    const auto edge_kind = edge_kind_by_pair.find({static_cast<std::uint32_t>(from), to});
                    if (edge_kind == edge_kind_by_pair.end() || edge_kind->second != 2U) continue;
                    if (!covered_edges.count({static_cast<std::uint32_t>(from), to}))
                        uncovered_edges.emplace_back(static_cast<std::uint32_t>(from), to);
                }
            }
            std::sort(uncovered_edges.begin(), uncovered_edges.end(), [&](const auto& left, const auto& right) {
                const auto left_support = read_edge_support.find(left);
                const auto right_support = read_edge_support.find(right);
                const auto left_count = left_support == read_edge_support.end() ? 0U : left_support->second;
                const auto right_count = right_support == read_edge_support.end() ? 0U : right_support->second;
                if (left_count != right_count) return left_count > right_count;
                return left < right;
            });
            std::set<std::string> existing_sequences(result.haplotype_path_sequences.begin(),
                                                      result.haplotype_path_sequences.end());
            for (const auto& [from, to] : uncovered_edges) {
                if (result.haplotype_path_count >= max_paths) break;
                const auto source_path = best_path_for_node[from];
                if (source_path == std::numeric_limits<std::size_t>::max() ||
                    source_path >= path_nodes.size()) continue;
                const auto& source_nodes = path_nodes[source_path];
                const auto source_iter = std::find(source_nodes.rbegin(), source_nodes.rend(), from);
                if (source_iter == source_nodes.rend()) continue;
                const auto source_index = static_cast<std::size_t>(source_iter - source_nodes.rbegin());
                const auto source_offset = source_nodes.size() - 1U - source_index;
                const auto& source_sequence = result.haplotype_path_sequences[source_path];
                const auto& target_kmer = encoded_nodes[to];
                auto artificial_sequence = source_sequence.substr(0, source_offset + options.k);
                artificial_sequence.push_back(target_kmer.back());
                std::set<std::uint32_t> visited_nodes(source_nodes.begin(),
                                                       source_nodes.begin() + source_offset + 1U);
                visited_nodes.insert(to);
                auto current = to;
                // `max_depth` protects the un-rejoined alternate branch; it
                // is not a cap on the completed reference prefix.  Preserve
                // the already materialized source path and allow the bounded
                // recovery branch beyond it, rather than truncating a long
                // AssemblyRegion merely because the user selected a small
                // graph resource guard.
                const auto artificial_depth_limit = options.max_depth == 0 ||
                    source_sequence.size() > std::numeric_limits<std::size_t>::max() - options.max_depth
                    ? std::numeric_limits<std::size_t>::max()
                    : source_sequence.size() + options.max_depth;
                while (artificial_sequence.size() < source_sequence.size() + options.k &&
                       artificial_sequence.size() < artificial_depth_limit &&
                       !directed_adjacency[current].empty()) {
                    const auto& candidates = directed_adjacency[current];
                    auto best_next = std::numeric_limits<std::uint32_t>::max();
                    std::uint32_t best_support = 0;
                    bool best_reference = false;
                    for (const auto next : candidates) {
                        if (visited_nodes.count(next)) continue;
                        const auto kind_iter = edge_kind_by_pair.find({current, next});
                        const auto is_reference_edge = kind_iter != edge_kind_by_pair.end() &&
                            kind_iter->second == 1U;
                        const auto support_iter = read_edge_support.find({current, next});
                        const auto support = support_iter == read_edge_support.end() ? 0U : support_iter->second;
                        if (best_next == std::numeric_limits<std::uint32_t>::max() ||
                            (is_reference_edge && !best_reference) ||
                            (is_reference_edge == best_reference && support > best_support) ||
                            (is_reference_edge == best_reference && support == best_support && next < best_next)) {
                            best_next = next;
                            best_support = support;
                            best_reference = is_reference_edge;
                        }
                    }
                    if (best_next == std::numeric_limits<std::uint32_t>::max()) break;
                    artificial_sequence.push_back(encoded_nodes[best_next].back());
                    visited_nodes.insert(best_next);
                    current = best_next;
                    if (is_reference[current] && artificial_sequence.size() >= options.k + 1U) break;
                }
                if (artificial_sequence.size() <= options.k || !existing_sequences.insert(artificial_sequence).second)
                    continue;
                const auto support_iter = read_edge_support.find({from, to});
                const auto edge_support = support_iter == read_edge_support.end() ? 0U : support_iter->second;
                const auto source_support = source_path < result.haplotype_path_support.size()
                    ? result.haplotype_path_support[source_path] : 0U;
                const auto source_score = source_path < result.haplotype_path_scores.size()
                    ? result.haplotype_path_scores[source_path]
                    : -std::numeric_limits<double>::infinity();
                result.haplotype_path_sequences.push_back(std::move(artificial_sequence));
                result.haplotype_path_tids.push_back(result.haplotype_path_tids[source_path]);
                result.haplotype_path_starts.push_back(result.haplotype_path_starts[source_path]);
                result.haplotype_path_ends.push_back(result.haplotype_path_ends[source_path]);
                result.haplotype_path_support.push_back(std::max(source_support, edge_support));
                result.haplotype_path_scores.push_back(source_score);
                result.haplotype_path_has_non_reference_edge.push_back(1U);
                result.haplotype_path_alt_read_starts.push_back(static_cast<std::uint32_t>(source_offset));
                result.haplotype_path_alt_read_ends.push_back(
                    static_cast<std::uint32_t>(result.haplotype_path_sequences.back().size()));
                ++result.haplotype_path_count;
                ++result.artificial_haplotype_recovery_paths;
                result.artificial_haplotype_recovery_bases += result.haplotype_path_sequences.back().size();
                covered_edges.emplace(from, to);
            }
        }
        result.reference_path_count = static_cast<std::size_t>(std::count_if(
            start_tids.begin(), start_tids.end(),
            [](const auto tid) { return tid >= 0; }));
    }
    // Run the actual SeqGraph rewrite over the materialized path language.
    // The raw de-Bruijn topology above remains available for support/cycle
    // accounting; this second representation is the sequence-level graph
    // consumed by HC/Mutect2 when they need shared diamond/tail/suffix nodes.
    std::vector<SeqGraphPath> seqgraph_paths;
    seqgraph_paths.reserve(result.haplotype_path_sequences.size());
    for (std::size_t path = 0; path < result.haplotype_path_sequences.size(); ++path) {
        seqgraph_paths.push_back(SeqGraphPath{
            result.haplotype_path_sequences[path],
            path < result.haplotype_path_tids.size() ? result.haplotype_path_tids[path] : -1,
            path < result.haplotype_path_starts.size() ? result.haplotype_path_starts[path] : -1,
            path < result.haplotype_path_ends.size() ? result.haplotype_path_ends[path] : -1,
            path < result.haplotype_path_has_non_reference_edge.size()
                ? result.haplotype_path_has_non_reference_edge[path] == 0 : false});
    }
    if (!options.disable_seqgraph_simplification && options.linked_de_bruijn_graph) {
        const auto seqgraph = simplify_seqgraph_paths(seqgraph_paths);
        result.seqgraph_nodes = seqgraph.final_nodes;
        result.seqgraph_edges = seqgraph.final_edges;
        result.seqgraph_linear_chain_merges = seqgraph.linear_chain_merges;
        result.seqgraph_diamond_merges = seqgraph.diamond_merges;
        result.seqgraph_tail_merges = seqgraph.tail_merges;
        result.seqgraph_suffix_splits = seqgraph.suffix_splits;
        result.seqgraph_suffix_merges = seqgraph.suffix_merges;
    }
    return result;
}

KmerGraphResult build_kmer_graph_kokkos(const KmerGraphInput& input,
                                        KmerGraphOptions options) {
    // GATK's ReadThreadingAssembler accepts a sorted, repeatable k-mer list
    // (the production default is [10,25]).  `assemble()` builds every
    // explicit graph, and `findBestPaths()` adds the paths from every
    // successful graph to one LinkedHashSet<Haplotype>.  A later k can thus
    // contribute an alternate the earlier graph did not resolve.  The k+10
    // expansion loop only runs when all explicit graphs fail.
    std::vector<std::uint32_t> requested_kmers = options.requested_kmer_sizes;
    if (requested_kmers.empty()) requested_kmers.push_back(options.k);
    for (const auto k : requested_kmers)
        if (k < 1)
            throw std::invalid_argument("requested k-mer sizes must be positive");
    std::sort(requested_kmers.begin(), requested_kmers.end());
    requested_kmers.erase(std::unique(requested_kmers.begin(), requested_kmers.end()),
                          requested_kmers.end());
    std::size_t maximum_sequence_length = 0;
    for (std::size_t read = 0; read + 1 < input.offsets.size(); ++read)
        maximum_sequence_length = std::max<std::size_t>(
            maximum_sequence_length, input.offsets[read + 1] - input.offsets[read]);
    for (std::size_t segment = 0; segment + 1 < input.reference_offsets.size(); ++segment)
        maximum_sequence_length = std::max<std::size_t>(
            maximum_sequence_length,
            input.reference_offsets[segment + 1] - input.reference_offsets[segment]);
    const auto maximum_k = static_cast<std::uint32_t>(std::min<std::size_t>(
        maximum_sequence_length, std::numeric_limits<std::uint32_t>::max()));
    // ReadThreadingAssembler adds exactly six k+10 attempts after the
    // explicit --kmer-size list. Its sixth expansion is special: GATK allows
    // a non-unique reference k-mer (and its low-complexity guard) only there.
    // Kokkos continues to process integer node IDs; Host keeps full strings
    // until that compact ID assignment, so this path is no longer limited to
    // the former 31-mer two-bit representation.
    constexpr std::size_t kMaxExpandedKmerAttempts = 6;
    KmerGraphResult selected;
    bool have_selected = false;
    bool saw_non_reference_cycle = false;
    std::size_t attempted_graphs = 0;

    const auto merge_haplotype_paths = [](KmerGraphResult& destination,
                                          const KmerGraphResult& source) {
        // Haplotype equality in the GATK result set is sequence based.  Walk
        // graphs in sorted k order and retain the first occurrence so the
        // Host coordinates and support match LinkedHashSet insertion order.
        std::set<std::string> seen(destination.haplotype_path_sequences.begin(),
                                   destination.haplotype_path_sequences.end());
        for (std::size_t path = 0; path < source.haplotype_path_sequences.size(); ++path) {
            const auto& sequence = source.haplotype_path_sequences[path];
            if (!seen.insert(sequence).second) continue;
            destination.haplotype_path_sequences.push_back(sequence);
            destination.haplotype_path_tids.push_back(
                path < source.haplotype_path_tids.size() ? source.haplotype_path_tids[path] : -1);
            destination.haplotype_path_starts.push_back(
                path < source.haplotype_path_starts.size() ? source.haplotype_path_starts[path] : -1);
            destination.haplotype_path_ends.push_back(
                path < source.haplotype_path_ends.size() ? source.haplotype_path_ends[path] : -1);
            destination.haplotype_path_support.push_back(
                path < source.haplotype_path_support.size() ? source.haplotype_path_support[path] : 0U);
            destination.haplotype_path_scores.push_back(
                path < source.haplotype_path_scores.size()
                    ? source.haplotype_path_scores[path]
                    : -std::numeric_limits<double>::infinity());
            destination.haplotype_path_has_non_reference_edge.push_back(
                path < source.haplotype_path_has_non_reference_edge.size()
                    ? source.haplotype_path_has_non_reference_edge[path] : 0U);
            destination.haplotype_path_alt_read_starts.push_back(
                path < source.haplotype_path_alt_read_starts.size()
                    ? source.haplotype_path_alt_read_starts[path] : 0U);
            destination.haplotype_path_alt_read_ends.push_back(
                path < source.haplotype_path_alt_read_ends.size()
                    ? source.haplotype_path_alt_read_ends[path] : 0U);
        }
        destination.haplotype_path_count = destination.haplotype_path_sequences.size();
    };

    const auto try_graph = [&](const std::uint32_t kmer_size,
                               const bool allow_non_unique_kmers) {
        KmerGraphOptions attempt_options = options;
        attempt_options.k = kmer_size;
        attempt_options.dont_increase_kmer_sizes_for_cycles = true;
        attempt_options.allow_non_unique_kmers_in_ref = allow_non_unique_kmers;
        auto attempt = build_kmer_graph_impl(input, attempt_options);
        ++attempted_graphs;
        attempt.kmer_iterations = attempted_graphs;
        const auto current_attempt_has_cycle = attempt.has_non_reference_cycles;
        if (std::getenv("FASTGATK_DEBUG_KMER_ATTEMPTS") != nullptr) {
            std::cerr << "[FASTGATK_KMER_ATTEMPT]"
                      << " k=" << kmer_size
                      << " reference_rejected=" << (attempt.reference_kmer_rejected ? 1 : 0)
                      << " non_reference_cycle=" << (current_attempt_has_cycle ? 1 : 0)
                      << " nodes=" << attempt.nodes
                      << " paths=" << attempt.haplotype_path_count
                      << '\n';
        }
        saw_non_reference_cycle = saw_non_reference_cycle || current_attempt_has_cycle;

        // createGraph() returns null for a rejected reference k-mer or a
        // cyclic ReadThreadingGraph.  A valid reference-only graph is still
        // an AssemblyResult and therefore prevents automatic k+10 retries.
        const bool accepted = !attempt.reference_kmer_rejected && !current_attempt_has_cycle;
        if (accepted) {
            if (!have_selected) {
                selected = std::move(attempt);
                have_selected = true;
            } else {
                merge_haplotype_paths(selected, attempt);
            }
        } else if (!have_selected) {
            // Keep diagnostics for an all-failed call without exposing its
            // paths as assembly candidates.
            selected = std::move(attempt);
        }
        return accepted;
    };

    // Explicit --kmer-size values are independent GATK graph attempts, not
    // a retry chain.  Do not stop after the first graph with an alternate.
    for (const auto kmer_size : requested_kmers)
        try_graph(kmer_size, options.allow_non_unique_kmers_in_ref);

    // Java enters the expansion loop only when its explicit result list is
    // empty, then stops at the first non-null result.  The sixth attempt is
    // the special low-complexity/non-unique-reference admission point; this
    // Host implementation carries the latter policy into the compact graph.
    if (!have_selected && !options.dont_increase_kmer_sizes_for_cycles) {
        auto current_k = requested_kmers.back();
        for (std::size_t expanded_attempts = 1;
             expanded_attempts <= kMaxExpandedKmerAttempts;
             ++expanded_attempts) {
            if (current_k >= maximum_k || current_k > maximum_k - 10U) break;
            current_k += 10U;
            const bool last_attempt = expanded_attempts == kMaxExpandedKmerAttempts;
            if (try_graph(current_k,
                          options.allow_non_unique_kmers_in_ref || last_attempt))
                break;
        }
    }

    selected.kmer_iterations = attempted_graphs;
    selected.has_non_reference_cycles = saw_non_reference_cycle;
    return selected;
}

}  // namespace fastgatk::kernels
