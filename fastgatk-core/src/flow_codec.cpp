#include "fastgatk/io/flow_codec.hpp"
#include "fastgatk/io/hts_reader.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fastgatk::io {
namespace {

constexpr std::size_t kProbabilityStride = 256;
constexpr double kMinimalCallProbability = 0.1;

std::uint8_t normalize_base(std::uint8_t value) {
    if (value >= static_cast<std::uint8_t>('a') && value <= static_cast<std::uint8_t>('z'))
        value = static_cast<std::uint8_t>(value - static_cast<std::uint8_t>('a') + static_cast<std::uint8_t>('A'));
    return value;
}

bool is_n(std::uint8_t value) {
    return normalize_base(value) == static_cast<std::uint8_t>('N');
}

double phred_probability(std::uint8_t quality) {
    return std::pow(10.0, -static_cast<double>(quality) / 10.0);
}

std::size_t clamp_hmer(std::int32_t value, std::size_t max_hmer) {
    if (value <= 0) return 0;
    return std::min<std::size_t>(static_cast<std::size_t>(value), max_hmer);
}

void spread_boundary(std::vector<double>& matrix, const std::vector<std::int32_t>& key,
                     std::size_t flow, std::size_t max_hmer, double filling_value) {
    if (flow >= key.size() || key[flow] <= 0) return;
    const auto call = clamp_hmer(key[flow], max_hmer);
    const auto count = max_hmer - call + 1;
    double total = 0.0;
    for (std::size_t hmer = call; hmer <= max_hmer; ++hmer)
        total += matrix[flow * kProbabilityStride + hmer];
    const double fill = std::max(total / static_cast<double>(count), filling_value);
    for (std::size_t hmer = call; hmer <= max_hmer; ++hmer)
        matrix[flow * kProbabilityStride + hmer] = fill;
}

}  // namespace

DecodedFlowRead decode_flow_read(const FlowReadTags& input) {
    if (input.bases.empty()) throw std::invalid_argument("flow codec requires non-empty bases");
    if (input.bases.size() != input.qualities.size())
        throw std::invalid_argument("flow codec bases/qualities length mismatch");
    if (input.tp.empty())
        throw std::invalid_argument("flow codec requires tp tag");
    if (input.tp.size() != input.bases.size())
        throw std::invalid_argument("flow codec tp/bases length mismatch");
    if (input.use_t0_tag && !input.t0_phred.empty() &&
        input.t0_phred.size() != input.bases.size())
        throw std::invalid_argument("flow codec t0/bases length mismatch");
    if (input.flow_order.empty()) throw std::invalid_argument("flow codec requires flow order");
    if (input.max_hmer == 0 || input.max_hmer >= kProbabilityStride)
        throw std::invalid_argument("flow codec max_hmer must be in [1,255]");
    if (!std::isfinite(input.filling_value) || input.filling_value < 0.0 || input.filling_value > 1.0)
        throw std::invalid_argument("flow codec filling value must be finite in [0,1]");
    if (input.use_t0_tag && input.t0_phred.empty())
        throw std::invalid_argument("flow codec requested t0 but tag is absent");

    const auto encoded = encode_flow_key(input.bases, input.flow_order);
    DecodedFlowRead decoded;
    decoded.key = encoded.key;
    decoded.flow_order = encoded.flow_order;
    decoded.max_hmer = input.max_hmer;
    decoded.filling_value = input.filling_value;
    if (decoded.key.empty()) throw std::invalid_argument("flow codec produced an empty key");
    decoded.probabilities.assign(decoded.key.size() * kProbabilityStride, input.filling_value);
    // GATK's default base-format decoder initializes all hmer rows to the
    // filling value, then accumulates tp observations into the selected row.
    std::size_t base_offset = 0;
    for (std::size_t flow = 0; flow < decoded.key.size(); ++flow) {
        const auto run = decoded.key[flow];
        if (run > 0) {
            for (std::size_t base = base_offset; base < base_offset + static_cast<std::size_t>(run); ++base) {
                const auto offset = input.tp[base];
                if (offset == 0) continue;
                const auto hmer = std::min<std::size_t>(
                    std::max<std::int32_t>(static_cast<std::int32_t>(run) + offset, 0),
                    static_cast<std::int32_t>(input.max_hmer));
                auto& cell = decoded.probabilities[flow * kProbabilityStride + hmer];
                const auto probability = phred_probability(input.qualities[base]);
                cell = (cell == input.filling_value) ? probability : cell + probability;
            }
        }
        if (run == 0 && input.use_t0_tag) {
            // This is FlowBasedRead.parseZeroQuals verbatim in data-flow form.
            if (base_offset != 0 && base_offset != input.bases.size()) {
                auto probability = std::min(phred_probability(input.t0_phred[base_offset - 1]),
                                            phred_probability(input.t0_phred[base_offset]));
                if (probability <= input.filling_value * 3.0) probability = 0.0;
                auto& cell = decoded.probabilities[flow * kProbabilityStride + 1];
                cell = std::max(cell, probability);
            }
        }
        double total_error = 0.0;
        for (std::size_t hmer = 0; hmer < input.max_hmer; ++hmer)
            total_error += decoded.probabilities[flow * kProbabilityStride + hmer];
        const auto call_probability = std::max(kMinimalCallProbability, 1.0 - total_error);
        decoded.probabilities[flow * kProbabilityStride + clamp_hmer(run, input.max_hmer)] = call_probability;
        // FlowBasedRead.applyFilteringFlowMatrix() always applies clipProbs
        // with the production defaults. Very small non-called probabilities
        // are quantization noise and are restored to the filling value; this
        // is observable in GATK's matrix oracle (for example a Q=26 error in
        // a non-called hmer row is clipped from 0.00251189 to 0.001).
        for (std::size_t hmer = 0; hmer < input.max_hmer; ++hmer) {
            if (hmer != clamp_hmer(run, input.max_hmer) &&
                decoded.probabilities[flow * kProbabilityStride + hmer] <= input.filling_value * 3.0)
                decoded.probabilities[flow * kProbabilityStride + hmer] = input.filling_value;
        }
        // FlowBasedRead.getProb(hmer) clamps hmer to maxHmer.  Materialize
        // that rule in the fixed 256-entry table so a haplotype containing a
        // longer homopolymer cannot accidentally read an uninitialized/zero
        // probability in the Kokkos kernel.
        const auto max_probability = decoded.probabilities[flow * kProbabilityStride + input.max_hmer];
        for (std::size_t hmer = input.max_hmer + 1; hmer < kProbabilityStride; ++hmer)
            decoded.probabilities[flow * kProbabilityStride + hmer] = max_probability;
        base_offset += static_cast<std::size_t>(std::max<std::int32_t>(run, 0));
    }

    if (!input.keep_boundary_flows) {
        auto first = std::find_if(decoded.key.begin(), decoded.key.end(), [](auto value) { return value != 0; });
        if (first != decoded.key.end()) spread_boundary(decoded.probabilities, decoded.key,
                                                         static_cast<std::size_t>(first - decoded.key.begin()),
                                                         input.max_hmer, input.filling_value);
        auto last = std::find_if(decoded.key.rbegin(), decoded.key.rend(), [](auto value) { return value != 0; });
        if (last != decoded.key.rend()) spread_boundary(decoded.probabilities, decoded.key,
                                                         decoded.key.size() - 1 - static_cast<std::size_t>(last - decoded.key.rbegin()),
                                                         input.max_hmer, input.filling_value);
    }
    for (std::size_t flow = 0; flow < decoded.key.size(); ++flow) {
        const auto max_probability = decoded.probabilities[flow * kProbabilityStride + input.max_hmer];
        for (std::size_t hmer = input.max_hmer + 1; hmer < kProbabilityStride; ++hmer)
            decoded.probabilities[flow * kProbabilityStride + hmer] = max_probability;
    }
    return decoded;
}

void clip_decoded_flow_read(DecodedFlowRead& decoded,
                            const std::size_t left_bases,
                            const std::size_t right_bases,
                            const bool spread_boundary_probabilities) {
    constexpr std::size_t probability_stride = kProbabilityStride;
    if (decoded.key.empty() || decoded.flow_order.size() != decoded.key.size() ||
        decoded.probabilities.size() != decoded.key.size() * probability_stride)
        throw std::invalid_argument("decoded flow read has inconsistent key/order/probability lengths");
    if (decoded.max_hmer == 0 || decoded.max_hmer >= probability_stride)
        throw std::invalid_argument("decoded flow max-hmer is outside probability matrix");
    if (!std::isfinite(decoded.filling_value) || decoded.filling_value < 0.0 ||
        decoded.filling_value > 1.0)
        throw std::invalid_argument("decoded flow filling value is invalid");
    std::size_t total_bases = 0;
    for (const auto value : decoded.key) {
        if (value < 0) throw std::invalid_argument("decoded flow key contains a negative hmer");
        total_bases += static_cast<std::size_t>(value);
    }
    if (left_bases > total_bases || right_bases > total_bases - left_bases ||
        left_bases + right_bases >= total_bases)
        throw std::invalid_argument("flow base clipping removes the complete read");

    // Reconstruct FlowBasedRead.flow2base/reverse-flow2base for the same
    // sentinel convention used by findLeft/RightClipping in GATK.
    EncodedFlowKey encoded;
    encoded.key = decoded.key;
    encoded.flow_order = decoded.flow_order;
    encoded.flow_to_base.resize(encoded.key.size());
    if (!encoded.flow_to_base.empty()) {
        encoded.flow_to_base[0] = -1;
        for (std::size_t flow = 1; flow < encoded.flow_to_base.size(); ++flow)
            encoded.flow_to_base[flow] = encoded.flow_to_base[flow - 1] + encoded.key[flow - 1];
    }
    encoded.reverse_key = encoded.key;
    std::reverse(encoded.reverse_key.begin(), encoded.reverse_key.end());
    encoded.reverse_flow_to_base.resize(encoded.reverse_key.size());
    if (!encoded.reverse_flow_to_base.empty()) {
        encoded.reverse_flow_to_base[0] = -1;
        for (std::size_t flow = 1; flow < encoded.reverse_flow_to_base.size(); ++flow)
            encoded.reverse_flow_to_base[flow] =
                encoded.reverse_flow_to_base[flow - 1] + encoded.reverse_key[flow - 1];
    }

    const auto left_pair = find_left_flow_clipping(left_bases, encoded);
    const auto right_pair = find_right_flow_clipping(right_bases, encoded);
    std::size_t clip_left = left_pair.first;
    std::size_t clip_right = right_pair.first;
    const auto left_hmer_clip = left_pair.second;
    const auto right_hmer_clip = right_pair.second;
    const auto original_length = decoded.key.size();
    if (clip_left >= original_length || clip_right >= original_length ||
        clip_left + clip_right >= original_length)
        throw std::invalid_argument("flow clipping resolved to an empty key");

    bool shift_left = true;
    if (left_hmer_clip > static_cast<std::size_t>(decoded.key[clip_left]))
        throw std::invalid_argument("left flow clipping exceeds boundary hmer");
    decoded.key[clip_left] -= static_cast<std::int32_t>(left_hmer_clip);
    while (decoded.key[clip_left] == 0) {
        ++clip_left;
        shift_left = false;
        if (clip_left >= original_length || clip_left + clip_right >= original_length)
            throw std::invalid_argument("left flow clipping removed the key");
    }

    bool shift_right = true;
    const auto right_index = original_length - 1 - clip_right;
    if (right_hmer_clip > static_cast<std::size_t>(decoded.key[right_index]))
        throw std::invalid_argument("right flow clipping exceeds boundary hmer");
    decoded.key[right_index] -= static_cast<std::int32_t>(right_hmer_clip);
    while (decoded.key[original_length - 1 - clip_right] == 0) {
        ++clip_right;
        shift_right = false;
        if (clip_left + clip_right >= original_length)
            throw std::invalid_argument("right flow clipping removed the key");
    }

    const auto retained_end = original_length - clip_right;
    if (clip_left >= retained_end)
        throw std::invalid_argument("flow clipping produced an empty retained key");
    std::vector<std::int32_t> new_key(
        decoded.key.begin() + static_cast<std::ptrdiff_t>(clip_left),
        decoded.key.begin() + static_cast<std::ptrdiff_t>(retained_end));
    std::vector<std::uint8_t> new_order(
        decoded.flow_order.begin() + static_cast<std::ptrdiff_t>(clip_left),
        decoded.flow_order.begin() + static_cast<std::ptrdiff_t>(retained_end));
    std::vector<double> new_probabilities(new_key.size() * probability_stride, 0.0);
    for (std::size_t flow = 0; flow < new_key.size(); ++flow) {
        const auto source = clip_left + flow;
        std::copy_n(decoded.probabilities.begin() +
                        static_cast<std::ptrdiff_t>(source * probability_stride),
                    probability_stride,
                    new_probabilities.begin() +
                        static_cast<std::ptrdiff_t>(flow * probability_stride));
    }
    const auto shift_column = [&](const std::size_t flow, const std::size_t shift) {
        if (shift == 0) return;
        if (shift >= probability_stride || flow >= new_key.size())
            throw std::invalid_argument("flow hmer shift exceeds probability matrix");
        auto* column = new_probabilities.data() + flow * probability_stride;
        for (std::size_t hmer = 0; hmer + shift < probability_stride; ++hmer)
            column[hmer] = column[hmer + shift];
        for (std::size_t hmer = probability_stride - shift; hmer < probability_stride; ++hmer)
            column[hmer] = 0.0;
    };
    if (shift_left) shift_column(0, left_hmer_clip);
    if (shift_right) shift_column(new_key.size() - 1, right_hmer_clip);

    decoded.key = std::move(new_key);
    decoded.flow_order = std::move(new_order);
    decoded.probabilities = std::move(new_probabilities);
    if (spread_boundary_probabilities) {
        const auto first = std::find_if(decoded.key.begin(), decoded.key.end(),
                                        [](const auto value) { return value != 0; });
        if (first != decoded.key.end())
            spread_boundary(decoded.probabilities, decoded.key,
                            static_cast<std::size_t>(first - decoded.key.begin()),
                            decoded.max_hmer, decoded.filling_value);
        const auto last = std::find_if(decoded.key.rbegin(), decoded.key.rend(),
                                       [](const auto value) { return value != 0; });
        if (last != decoded.key.rend())
            spread_boundary(decoded.probabilities, decoded.key,
                            decoded.key.size() - 1 - static_cast<std::size_t>(last - decoded.key.rbegin()),
                            decoded.max_hmer, decoded.filling_value);
    }
}

EncodedFlowKey encode_flow_key(const std::vector<std::uint8_t>& bases,
                               const std::string& flow_order) {
    if (bases.empty()) throw std::invalid_argument("flow key encoder requires non-empty bases");
    if (flow_order.empty()) throw std::invalid_argument("flow key encoder requires flow order");
    EncodedFlowKey encoded;
    std::size_t location = 0;
    std::size_t flow_number = 0;
    std::size_t period_guard = 0;
    while (location < bases.size()) {
        const auto flow_base = normalize_base(static_cast<std::uint8_t>(
            flow_order[flow_number % flow_order.size()]));
        const auto current = normalize_base(bases[location]);
        if (current != flow_base && !is_n(current)) {
            encoded.key.push_back(0);
            if (++period_guard > flow_order.size())
                throw std::invalid_argument("flow key base is absent from flow order");
        } else {
            std::size_t count = 0;
            while (location < bases.size()) {
                const auto base = normalize_base(bases[location]);
                if (base != flow_base && !is_n(base)) break;
                ++location;
                ++count;
            }
            encoded.key.push_back(static_cast<std::int32_t>(count));
            period_guard = 0;
        }
        ++flow_number;
    }
    encoded.flow_order.resize(encoded.key.size());
    for (std::size_t i = 0; i < encoded.flow_order.size(); ++i)
        encoded.flow_order[i] = static_cast<std::uint8_t>(flow_order[i % flow_order.size()]);
    encoded.flow_to_base.resize(encoded.key.size());
    if (!encoded.flow_to_base.empty()) {
        encoded.flow_to_base[0] = -1;
        for (std::size_t i = 1; i < encoded.flow_to_base.size(); ++i)
            encoded.flow_to_base[i] = encoded.flow_to_base[i - 1] + encoded.key[i - 1];
    }
    encoded.reverse_key = encoded.key;
    std::reverse(encoded.reverse_key.begin(), encoded.reverse_key.end());
    encoded.reverse_flow_to_base.resize(encoded.reverse_key.size());
    if (!encoded.reverse_flow_to_base.empty()) {
        encoded.reverse_flow_to_base[0] = -1;
        for (std::size_t i = 1; i < encoded.reverse_flow_to_base.size(); ++i)
            encoded.reverse_flow_to_base[i] = encoded.reverse_flow_to_base[i - 1] +
                encoded.reverse_key[i - 1];
    }
    return encoded;
}

std::vector<std::uint8_t> collapse_flow_homopolymers(
    const std::vector<std::uint8_t>& bases, const std::size_t hmer_size_threshold) {
    if (hmer_size_threshold == 0 || bases.empty()) return bases;
    if (hmer_size_threshold > std::numeric_limits<std::size_t>::max() / 2)
        throw std::invalid_argument("flow hmer collapse threshold is too large");

    // This is intentionally the same state machine as GATK's
    // LongHomopolymerHaplotypeCollapsingEngine.collapseBases(): the first
    // homopolymer is left untouched, while later runs retain at most the
    // configured threshold bases.  Keeping the first run is important because
    // the Java implementation uses it as the leading-boundary anchor.
    std::vector<std::uint8_t> collapsed;
    collapsed.reserve(bases.size());
    std::uint8_t last_base = 0;
    std::size_t base_same_count = 0;
    bool first_homopolymer = true;
    for (const auto base : bases) {
        if (base == last_base) {
            ++base_same_count;
            if (!first_homopolymer && base_same_count >= hmer_size_threshold)
                continue;
        } else {
            if (last_base != 0) first_homopolymer = false;
            last_base = base;
            base_same_count = 0;
        }
        collapsed.push_back(base);
    }
    return collapsed;
}

bool flow_homopolymer_exceeds_threshold(
    const std::vector<std::uint8_t>& bases, const std::size_t hmer_size_threshold) {
    if (hmer_size_threshold == 0 || bases.empty()) return false;
    std::uint8_t last_base = 0;
    std::size_t base_same_count = 0;
    for (const auto base : bases) {
        if (base == last_base) {
            ++base_same_count;
            // Java starts the counter at zero for the first base in a run, so
            // >= threshold means a run of threshold+1 bases.
            if (base_same_count >= hmer_size_threshold) return true;
        } else {
            last_base = base;
            base_same_count = 0;
        }
    }
    return false;
}

std::pair<std::size_t, std::size_t> find_left_flow_clipping(
    std::size_t base_clipping, const EncodedFlowKey& encoded) {
    if (base_clipping == 0) return {0, 0};
    if (encoded.key.empty() || encoded.flow_to_base.size() != encoded.key.size())
        throw std::invalid_argument("left flow clipping requires a complete key mapping");
    for (std::size_t flow = 0; flow < encoded.key.size(); ++flow) {
        const auto first = encoded.flow_to_base[flow];
        // GATK's getKeyToBase starts at -1.  That sentinel is intentional:
        // clipping inside the first homopolymer returns hmerClip=baseClip-1,
        // rather than skipping to the next spacer flow.
        const auto span_end = static_cast<std::int64_t>(first) +
            static_cast<std::int64_t>(std::max<std::int32_t>(encoded.key[flow], 0));
        if (span_end >= static_cast<std::int64_t>(base_clipping)) {
            const auto hmer_clip = static_cast<std::int64_t>(base_clipping) -
                static_cast<std::int64_t>(first) - 1;
            if (hmer_clip < 0) throw std::invalid_argument("left flow clipping produced a negative hmer clip");
            return {flow, static_cast<std::size_t>(hmer_clip)};
        }
    }
    throw std::invalid_argument("left flow clipping exceeds encoded read");
}

std::pair<std::size_t, std::size_t> find_right_flow_clipping(
    std::size_t base_clipping, const EncodedFlowKey& encoded) {
    if (base_clipping == 0) return {0, 0};
    if (encoded.reverse_key.empty() ||
        encoded.reverse_flow_to_base.size() != encoded.reverse_key.size())
        throw std::invalid_argument("right flow clipping requires a complete key mapping");
    for (std::size_t flow = 0; flow < encoded.reverse_key.size(); ++flow) {
        const auto first = encoded.reverse_flow_to_base[flow];
        const auto span_end = static_cast<std::int64_t>(first) +
            static_cast<std::int64_t>(std::max<std::int32_t>(encoded.reverse_key[flow], 0));
        if (span_end >= static_cast<std::int64_t>(base_clipping)) {
            const auto hmer_clip = static_cast<std::int64_t>(base_clipping) -
                static_cast<std::int64_t>(first) - 1;
            if (hmer_clip < 0) throw std::invalid_argument("right flow clipping produced a negative hmer clip");
            return {flow, static_cast<std::size_t>(hmer_clip)};
        }
    }
    throw std::invalid_argument("right flow clipping exceeds encoded read");
}

FlowReadTags flow_tags_from_batch(const ReadBatch& batch, std::size_t record,
                                  bool use_t0_tag, double filling_value,
                                  bool keep_boundary_flows) {
    const auto record_count = batch.positions.size();
    if (record >= record_count || batch.offsets.size() != record_count + 1)
        throw std::invalid_argument("flow batch record is out of range");
    if (batch.flow_tp_offsets.size() != record_count + 1 ||
        batch.flow_t0_offsets.size() != record_count + 1 ||
        batch.flow_order_offsets.size() != record_count + 1 ||
        batch.flow_max_hmer.size() != record_count)
        throw std::invalid_argument("flow metadata is not present for this ReadBatch");
    const auto begin = batch.offsets[record];
    const auto end = batch.offsets[record + 1];
    const auto tp_begin = batch.flow_tp_offsets[record];
    const auto tp_end = batch.flow_tp_offsets[record + 1];
    const auto t0_begin = batch.flow_t0_offsets[record];
    const auto t0_end = batch.flow_t0_offsets[record + 1];
    const auto order_begin = batch.flow_order_offsets[record];
    const auto order_end = batch.flow_order_offsets[record + 1];
    if (begin > end || end > batch.bases.size() || end > batch.qualities.size() ||
        tp_begin > tp_end || tp_end > batch.flow_tp.size() ||
        t0_begin > t0_end || t0_end > batch.flow_t0_phred.size() ||
        order_begin > order_end || order_end > batch.flow_orders.size())
        throw std::invalid_argument("flow metadata offsets are malformed");
    FlowReadTags result;
    result.bases.assign(batch.bases.begin() + begin, batch.bases.begin() + end);
    result.qualities.assign(batch.qualities.begin() + begin, batch.qualities.begin() + end);
    result.tp.assign(batch.flow_tp.begin() + tp_begin, batch.flow_tp.begin() + tp_end);
    result.t0_phred.assign(batch.flow_t0_phred.begin() + t0_begin, batch.flow_t0_phred.begin() + t0_end);
    result.flow_order.assign(batch.flow_orders.begin() + order_begin, batch.flow_orders.begin() + order_end);
    result.max_hmer = batch.flow_max_hmer[record] == 0 ? 12 : batch.flow_max_hmer[record];
    result.filling_value = filling_value;
    result.use_t0_tag = use_t0_tag;
    result.keep_boundary_flows = keep_boundary_flows;
    return result;
}

}  // namespace fastgatk::io
