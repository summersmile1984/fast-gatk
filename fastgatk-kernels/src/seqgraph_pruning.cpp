#include "fastgatk/kernels/kmer_graph.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fastgatk::kernels {
namespace {

struct Chain {
    std::vector<std::size_t> edge_ids;
    std::uint32_t first_vertex = 0;
    std::uint32_t last_vertex = 0;
};

double fast_bernoulli_entropy(const double p) {
    const auto product = p * (1.0 - p);
    return product * (11.0 + 33.0 * product) / (2.0 + 20.0 * product);
}

double gatk_digamma_positive(const double input) {
    // This is Apache Commons Math 3.x Gamma.digamma (Bernardo AS 103), which
    // is the implementation used by MathUtils.digamma in GATK. The cutoff
    // and recursive subtraction order matter for the cached integer values;
    // a generic recurrence-to-eight approximation differs by several ulps at
    // the read depths where adaptive pruning makes a branch decision.
    constexpr double kEulerGamma = 0.577215664901532860606512090082;
    constexpr double kSmall = 1.0e-5;
    constexpr double kLarge = 49.0;
    if (std::isnan(input) || std::isinf(input)) return input;
    if (!(input > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (input <= kSmall) return -kEulerGamma - 1.0 / input;

    double x = input;
    // Java's implementation is recursive: digamma(x) = digamma(x+1)-1/x.
    // Accumulate the same terms in the unwind order (largest x first).
    while (x < kLarge) x += 1.0;
    double result;
    {
        const double inverse_squared = 1.0 / (x * x);
        result = std::log(x) - 0.5 / x - inverse_squared *
            (1.0 / 12.0 + inverse_squared * (1.0 / 120.0 - inverse_squared / 252.0));
    }
    double lower = x - 1.0;
    while (lower >= input) {
        result -= 1.0 / lower;
        lower -= 1.0;
    }
    return result;
}

std::uint64_t gcd_u64(std::uint64_t left, std::uint64_t right) {
    while (right != 0U) {
        const auto remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

double gatk_binomial_coefficient_log(const std::uint64_t n,
                                     std::uint64_t k) {
    if (k == 0U || k == n) return 0.0;
    if (k == 1U || k == n - 1U) return std::log(static_cast<double>(n));
    if (k > n / 2U) k = n - k;

    // Commons Math computes small coefficients exactly, then takes the log.
    // This avoids the lgamma rounding drift seen in the Java oracle.
    if (n < 67U) {
        std::uint64_t coefficient = 1U;
        std::uint64_t factor = n - k + 1U;
        for (std::uint64_t divisor = 1U; divisor <= k; ++divisor, ++factor) {
            const auto divisor_gcd = gcd_u64(factor, divisor);
            coefficient = (coefficient / (divisor / divisor_gcd)) *
                (factor / divisor_gcd);
        }
        return std::log(static_cast<double>(coefficient));
    }
    if (n < 1030U) {
        double coefficient = 1.0;
        for (std::uint64_t i = 1U; i <= k; ++i)
            coefficient *= static_cast<double>(n - k + i) / static_cast<double>(i);
        coefficient = std::floor(coefficient + 0.5);
        return std::log(coefficient);
    }

    // For values that overflow Commons Math's double coefficient, mirror its
    // explicit log-sum path rather than calling lgamma.
    double log_sum = 0.0;
    for (std::uint64_t i = n - k + 1U; i <= n; ++i)
        log_sum += std::log(static_cast<double>(i));
    for (std::uint64_t i = 2U; i <= k; ++i)
        log_sum -= std::log(static_cast<double>(i));
    return log_sum;
}

double constant_error_log_likelihood_ratio(const std::uint32_t n_ref,
                                           const std::uint32_t n_alt,
                                           const double error_probability) {
    if (n_alt == 0U) return 0.0;
    const auto bounded_error = std::clamp(error_probability,
                                          std::numeric_limits<double>::min(),
                                          1.0 - std::numeric_limits<double>::epsilon());
    // QualityUtils.errorProbToQual uses Math.round(-10*log10(p)) and caps the
    // resulting SAM quality to [1,93].  Reusing that quantization is required
    // for the same branch decisions as Mutect2Engine.logLikelihoodRatio.
    const auto quality = std::clamp(static_cast<int>(std::llround(
        -10.0 * std::log10(bounded_error))), 1, 93);
    // QualityUtils uses Math.pow for epsilon, while NaturalLogUtils computes
    // log(epsilon) from the fixed Phred factor and evaluates log(1-epsilon)
    // through its two-branch log1mexp implementation. Keep both paths rather
    // than deriving one value with log/pow, which changes the low-Q branch.
    const auto epsilon = std::pow(10.0, -static_cast<double>(quality) / 10.0);
    const auto f_tilde_ratio = std::exp(
        gatk_digamma_positive(static_cast<double>(n_ref) + 1.0) -
        gatk_digamma_positive(static_cast<double>(n_alt) + 1.0));
    const auto z_bar_alt = (1.0 - epsilon) /
        (1.0 - epsilon + epsilon * f_tilde_ratio);
    const auto log_epsilon = static_cast<double>(quality) * (-std::log(10.0) / 10.0);
    const auto log_one_minus_epsilon = log_epsilon < std::log(0.5) ?
        std::log1p(-std::exp(log_epsilon)) : std::log(-std::expm1(log_epsilon));
    const auto read_term = z_bar_alt *
        (log_one_minus_epsilon - log_epsilon) + fast_bernoulli_entropy(z_bar_alt);
    const auto n = static_cast<std::uint64_t>(n_ref) + static_cast<std::uint64_t>(n_alt);
    const auto beta_entropy = -std::log(static_cast<double>(n) + 1.0) -
        gatk_binomial_coefficient_log(n, n_alt);
    return beta_entropy + static_cast<double>(n_alt) * read_term;
}

std::uint32_t pruning_multiplicity(const SeqGraphEdge& edge,
                                   const std::uint32_t required_samples) {
    if (required_samples == 0U) return edge.multiplicity;
    if (edge.sample_multiplicities.size() < required_samples) return 0U;
    auto values = edge.sample_multiplicities;
    std::sort(values.begin(), values.end(), std::greater<>());
    return values[required_samples - 1U];
}

struct ChainScore {
    double left = 0.0;
    double right = 0.0;
    bool first_edge_reference = false;
};

std::vector<ChainScore> score_chains(
    const std::vector<Chain>& chains,
    const std::vector<SeqGraphEdge>& edges,
    const std::vector<std::vector<std::size_t>>& outgoing,
    const std::vector<std::vector<std::size_t>>& incoming,
    const double error_rate) {
    std::vector<ChainScore> scores(chains.size());
    for (std::size_t chain_index = 0; chain_index < chains.size(); ++chain_index) {
        const auto& chain = chains[chain_index];
        if (chain.edge_ids.empty()) continue;
        const auto first_edge = chain.edge_ids.front();
        const auto last_edge = chain.edge_ids.back();
        scores[chain_index].first_edge_reference = edges[first_edge].reference;
        std::uint32_t left_total = 0;
        for (const auto edge_id : outgoing[chain.first_vertex])
            left_total += edges[edge_id].multiplicity;
        std::uint32_t right_total = 0;
        for (const auto edge_id : incoming[chain.last_vertex])
            right_total += edges[edge_id].multiplicity;
        const auto left_multiplicity = edges[first_edge].multiplicity;
        const auto right_multiplicity = edges[last_edge].multiplicity;
        // A source/sink has no competing chain on that side, exactly as
        // AdaptiveChainPruner's graph.isSource/isSink special case.
        scores[chain_index].left = incoming[chain.first_vertex].empty() ? 0.0 :
            constant_error_log_likelihood_ratio(
                left_total >= left_multiplicity ? left_total - left_multiplicity : 0U,
                left_multiplicity, error_rate);
        scores[chain_index].right = outgoing[chain.last_vertex].empty() ? 0.0 :
            constant_error_log_likelihood_ratio(
                right_total >= right_multiplicity ? right_total - right_multiplicity : 0U,
                right_multiplicity, error_rate);
    }
    return scores;
}

std::vector<std::uint8_t> likely_error_chains(
    const std::vector<Chain>& chains,
    const std::vector<SeqGraphEdge>& edges,
    const std::vector<ChainScore>& scores,
    const std::vector<std::vector<std::size_t>>& chains_by_start,
    const std::vector<std::vector<std::size_t>>& chains_by_end,
    const double log_odds_threshold,
    const double seeding_log_odds_threshold,
    const std::uint32_t max_unpruned_variants,
    const std::vector<std::string>* chain_first_sequences,
    const std::vector<std::string>* chain_bases) {
    std::vector<std::uint8_t> good(chains.size(), 0U);
    if (chains.empty()) return good;
    auto is_good_incoming = [&](const std::size_t chain) {
        return scores[chain].right >= log_odds_threshold ||
            scores[chain].first_edge_reference;
    };
    auto is_good_outgoing = [&](const std::size_t chain) {
        return scores[chain].left >= log_odds_threshold ||
            scores[chain].first_edge_reference;
    };
    auto is_seedable = [&](const std::size_t chain) {
        return scores[chain].right >= seeding_log_odds_threshold &&
            scores[chain].left >= seeding_log_odds_threshold;
    };
    std::vector<std::vector<std::size_t>> seedable_by_vertex(chains_by_start.size());
    for (std::size_t chain = 0; chain < chains.size(); ++chain) {
        if (!is_seedable(chain)) continue;
        seedable_by_vertex[chains[chain].first_vertex].push_back(chain);
        seedable_by_vertex[chains[chain].last_vertex].push_back(chain);
    }

    struct QueueEntry {
        double score = 0.0;
        std::size_t chain = 0;
    };
    struct QueueOrder {
        const std::vector<std::string>* first_sequences = nullptr;
        const std::vector<std::string>* bases = nullptr;

        bool operator()(const QueueEntry& left, const QueueEntry& right) const {
            // Java's PriorityQueue removes the greatest score first (its
            // comparator stores -score), then uses BaseUtils.BASES_COMPARATOR
            // on the first vertex and full chain bases.  `priority_queue`
            // places the element for which this predicate is false at top.
            if (left.score != right.score) return left.score < right.score;
            if (first_sequences != nullptr && bases != nullptr &&
                left.chain < first_sequences->size() && right.chain < first_sequences->size() &&
                left.chain < bases->size() && right.chain < bases->size()) {
                const auto& left_first = (*first_sequences)[left.chain];
                const auto& right_first = (*first_sequences)[right.chain];
                if (left_first != right_first) return left_first > right_first;
                const auto& left_bases = (*bases)[left.chain];
                const auto& right_bases = (*bases)[right.chain];
                if (left_bases != right_bases) return left_bases > right_bases;
            }
            return left.chain > right.chain;
        }
    };
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, QueueOrder> queue(
        QueueOrder{chain_first_sequences, chain_bases});
    auto enqueue = [&](const std::size_t chain, const double score) {
        queue.push(QueueEntry{score, chain});
    };
    const auto max_chain = std::max_element(chains.begin(), chains.end(), [&](const auto& left, const auto& right) {
        // AdaptiveChainPruner.getMaxWeightChain orders by the greatest
        // multiplicity of any edge in the chain (then by path length), not
        // by the sum of all edge multiplicities.  Summing a long, moderate
        // chain can incorrectly seed it over a short, high-confidence chain
        // and changes which component survives the bounded traversal.
        const auto left_weight = std::accumulate(left.edge_ids.begin(), left.edge_ids.end(), std::uint32_t{0},
            [&](const auto value, const auto edge_id) {
                return std::max(value, edges[edge_id].multiplicity);
            });
        const auto right_weight = std::accumulate(right.edge_ids.begin(), right.edge_ids.end(), std::uint32_t{0},
            [&](const auto value, const auto edge_id) {
                return std::max(value, edges[edge_id].multiplicity);
            });
        if (left_weight != right_weight) return left_weight < right_weight;
        if (left.edge_ids.size() != right.edge_ids.size()) return left.edge_ids.size() < right.edge_ids.size();
        const auto left_index = static_cast<std::size_t>(&left - chains.data());
        const auto right_index = static_cast<std::size_t>(&right - chains.data());
        if (chain_first_sequences != nullptr &&
            left_index < chain_first_sequences->size() && right_index < chain_first_sequences->size())
            return (*chain_first_sequences)[left_index] < (*chain_first_sequences)[right_index];
        return left.edge_ids.front() < right.edge_ids.front();
    });
    enqueue(static_cast<std::size_t>(max_chain - chains.begin()), std::numeric_limits<double>::infinity());
    for (std::size_t vertex = 0; vertex < seedable_by_vertex.size(); ++vertex) {
        if (seedable_by_vertex[vertex].size() <= 2U) continue;
        for (const auto chain : chains_by_start[vertex])
            if (is_good_outgoing(chain)) enqueue(chain, scores[chain].left);
        for (const auto chain : chains_by_end[vertex])
            if (is_good_incoming(chain)) enqueue(chain, scores[chain].right);
    }

    std::vector<std::uint8_t> processed_vertex(chains_by_start.size(), 0U);
    std::vector<std::uint8_t> has_good_outgoing(chains_by_start.size(), 0U);
    std::uint32_t variant_count = 0U;
    // Preserve AdaptiveChainPruner's loop guard.  It intentionally stops
    // consuming the queue as soon as the bounded good subgraph has exceeded
    // maxUnprunedVariants; continuing to mark queued chains as good retains
    // hundreds of low-support branches in high-depth regions.
    while (!queue.empty() && variant_count <= max_unpruned_variants) {
        const auto chain = queue.top().chain;
        queue.pop();
        if (good[chain] != 0U) continue;
        good[chain] = 1U;
        const auto first = chains[chain].first_vertex;
        const auto new_variant = has_good_outgoing[first] != 0U;
        has_good_outgoing[first] = 1U;
        if (new_variant) {
            ++variant_count;
            if (variant_count > max_unpruned_variants) continue;
        }
        for (const auto vertex : {chains[chain].first_vertex, chains[chain].last_vertex}) {
            if (processed_vertex[vertex] != 0U) continue;
            processed_vertex[vertex] = 1U;
            for (const auto candidate : chains_by_start[vertex])
                if (is_good_outgoing(candidate)) enqueue(candidate, scores[candidate].left);
            for (const auto candidate : chains_by_end[vertex])
                if (is_good_incoming(candidate)) enqueue(candidate, scores[candidate].right);
        }
    }
    // Reference chains are never discarded, even if their endpoints were not
    // reached by the bounded priority traversal.
    for (std::size_t chain = 0; chain < chains.size(); ++chain)
        if (std::any_of(chains[chain].edge_ids.begin(), chains[chain].edge_ids.end(),
                        [&](const auto edge_id) { return edges[edge_id].reference; }))
            good[chain] = 1U;
    return good;
}

}  // namespace

double seqgraph_constant_error_log_likelihood_ratio(
    const std::uint32_t n_ref,
    const std::uint32_t n_alt,
    const double error_probability) {
    return constant_error_log_likelihood_ratio(n_ref, n_alt, error_probability);
}

SeqGraphPruningResult prune_seqgraph_chains(
    const std::size_t node_count,
    const std::vector<SeqGraphEdge>& edges,
    const KmerGraphOptions& options,
    const std::vector<std::string>* node_sequences) {
    SeqGraphPruningResult result;
    result.keep_edges.assign(edges.size(), 1U);
    if (node_count == 0U || edges.empty()) return result;
    std::vector<std::vector<std::size_t>> outgoing(node_count), incoming(node_count);
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        if (edges[edge].from >= node_count || edges[edge].to >= node_count)
            throw std::invalid_argument("SeqGraph pruning edge endpoint out of range");
        outgoing[edges[edge].from].push_back(edge);
        incoming[edges[edge].to].push_back(edge);
    }
    std::vector<Chain> chains;
    auto build_chain = [&](const std::size_t first_edge) {
        Chain chain;
        chain.first_vertex = edges[first_edge].from;
        auto edge_id = first_edge;
        while (true) {
            chain.edge_ids.push_back(edge_id);
            chain.last_vertex = edges[edge_id].to;
            const auto vertex = chain.last_vertex;
            if (vertex == chain.first_vertex || outgoing[vertex].size() != 1U ||
                incoming[vertex].size() > 1U) break;
            const auto next = outgoing[vertex].front();
            edge_id = next;
        }
        chains.push_back(std::move(chain));
    };
    // Match ChainPruner.findAllChains(): seed only source vertices, then
    // enqueue each previously unseen chain endpoint.  In particular, do not
    // synthesize fallback chains for cyclic components; the Java pruner does
    // not visit them, and ReadThreadingAssembler handles cycles separately.
    std::deque<std::uint32_t> chain_starts;
    std::vector<std::uint8_t> already_seen(node_count, 0U);
    for (std::size_t vertex = 0; vertex < node_count; ++vertex) {
        if (!incoming[vertex].empty()) continue;
        chain_starts.push_back(static_cast<std::uint32_t>(vertex));
        already_seen[vertex] = 1U;
    }
    while (!chain_starts.empty()) {
        const auto start = chain_starts.front();
        chain_starts.pop_front();
        for (const auto edge : outgoing[start]) {
            build_chain(edge);
            const auto end = chains.back().last_vertex;
            if (already_seen[end] == 0U) {
                chain_starts.push_back(end);
                already_seen[end] = 1U;
            }
        }
    }
    result.chain_count = chains.size();
    if (chains.empty()) return result;

    std::vector<std::string> chain_first_sequences;
    std::vector<std::string> chain_base_sequences;
    if (node_sequences != nullptr && node_sequences->size() >= node_count) {
        chain_first_sequences.reserve(chains.size());
        chain_base_sequences.reserve(chains.size());
        for (const auto& chain : chains) {
            auto bases = (*node_sequences)[chain.first_vertex];
            for (const auto edge_id : chain.edge_ids) {
                const auto target = edges[edge_id].to;
                if (target >= node_sequences->size() || (*node_sequences)[target].empty()) continue;
                bases.push_back((*node_sequences)[target].back());
            }
            chain_first_sequences.push_back((*node_sequences)[chain.first_vertex]);
            chain_base_sequences.push_back(std::move(bases));
        }
    }
    const auto* first_sequences = chain_first_sequences.empty() ? nullptr : &chain_first_sequences;
    const auto* bases_sequences = chain_base_sequences.empty() ? nullptr : &chain_base_sequences;

    std::vector<std::vector<std::size_t>> chains_by_start(node_count), chains_by_end(node_count);
    for (std::size_t chain = 0; chain < chains.size(); ++chain) {
        chains_by_start[chains[chain].first_vertex].push_back(chain);
        chains_by_end[chains[chain].last_vertex].push_back(chain);
    }
    if (!options.use_adaptive_pruning) {
        for (const auto& chain : chains) {
            const auto has_reference = std::any_of(chain.edge_ids.begin(), chain.edge_ids.end(),
                [&](const auto edge_id) { return edges[edge_id].reference; });
            const auto needs_pruning = !has_reference && options.min_pruning > 0U &&
                std::all_of(chain.edge_ids.begin(), chain.edge_ids.end(), [&](const auto edge_id) {
                    return pruning_multiplicity(edges[edge_id], options.num_pruning_samples) < options.min_pruning;
                });
            if (!needs_pruning) continue;
            ++result.pruned_chain_count;
            for (const auto edge_id : chain.edge_ids) result.keep_edges[edge_id] = 0U;
        }
        return result;
    }

    const auto initial_scores = score_chains(chains, edges, outgoing, incoming,
                                             options.initial_error_rate_for_pruning);
    const auto initial_bad = likely_error_chains(
        chains, edges, initial_scores, chains_by_start, chains_by_end,
        options.pruning_log_odds_threshold,
        options.pruning_seeding_log_odds_threshold,
        options.max_unpruned_variants, first_sequences, bases_sequences);
    std::uint64_t error_count = 0;
    std::uint64_t total_bases = 0;
    for (const auto& chain : chains) {
        for (const auto edge_id : chain.edge_ids) total_bases += edges[edge_id].multiplicity;
        // `likely_error_chains` returns the chains retained by the bounded
        // traversal (the same polarity as AdaptiveChainPruner's `good`
        // bitmap).  Only chains not retained contribute their terminal edge
        // multiplicity to the empirical error-rate estimate.
        if (initial_bad[&chain - chains.data()] == 0U)
            error_count += edges[chain.edge_ids.back()].multiplicity;
    }
    const auto initial_rate = options.initial_error_rate_for_pruning;
    result.estimated_error_rate = total_bases == 0U ? initial_rate :
        std::clamp(static_cast<double>(error_count) / static_cast<double>(total_bases),
                   std::numeric_limits<double>::min(), 1.0 - std::numeric_limits<double>::epsilon());
    const auto final_scores = score_chains(chains, edges, outgoing, incoming,
                                           result.estimated_error_rate);
    const auto final_bad = likely_error_chains(
        chains, edges, final_scores, chains_by_start, chains_by_end,
        options.pruning_log_odds_threshold,
        options.pruning_seeding_log_odds_threshold,
        options.max_unpruned_variants, first_sequences, bases_sequences);
    for (std::size_t chain = 0; chain < chains.size(); ++chain) {
        const auto has_reference = std::any_of(chains[chain].edge_ids.begin(), chains[chain].edge_ids.end(),
            [&](const auto edge_id) { return edges[edge_id].reference; });
        // `final_bad` is the retained/good bitmap returned by the traversal;
        // only a chain that was not retained is an adaptive-pruning target.
        if (final_bad[chain] != 0U || has_reference) continue;
        ++result.adaptive_pruned_chain_count;
        for (const auto edge_id : chains[chain].edge_ids) result.keep_edges[edge_id] = 0U;
    }
    return result;
}

}  // namespace fastgatk::kernels
