#include "fastgatk/calling/pipeline.hpp"

#include "fastgatk/io/flow_codec.hpp"
#include "fastgatk/kernels/smith_waterman.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fastgatk::calling {
namespace {

std::uint8_t complement(const std::uint8_t base) {
    switch (base) {
        case 'A': case 'a': return static_cast<std::uint8_t>('T');
        case 'C': case 'c': return static_cast<std::uint8_t>('G');
        case 'G': case 'g': return static_cast<std::uint8_t>('C');
        case 'T': case 't': return static_cast<std::uint8_t>('A');
        default: return base;
    }
}

std::vector<std::uint8_t> reverse_complement(
    const std::vector<std::uint8_t>& input) {
    std::vector<std::uint8_t> output(input.size());
    for (std::size_t index = 0; index < input.size(); ++index)
        output[index] = complement(input[input.size() - 1 - index]);
    return output;
}

bool same_base(const std::vector<std::uint8_t>& sequence,
               const std::int64_t offset,
               const std::uint8_t base,
               const std::size_t length) {
    if (offset < 0 || static_cast<std::uint64_t>(offset) > sequence.size() ||
        length > sequence.size() - static_cast<std::size_t>(offset))
        return false;
    for (std::size_t index = 0; index < length; ++index)
        if (sequence[static_cast<std::size_t>(offset) + index] != base) return false;
    return true;
}

bool needs_collapsing(const std::vector<std::uint8_t>& sequence,
                      const std::size_t threshold) {
    if (sequence.empty()) return false;
    std::uint8_t last = 0;
    std::size_t same = 0;
    for (const auto base : sequence) {
        if (base == last) {
            if (++same >= threshold) return true;
        } else {
            last = base;
            same = 0;
        }
    }
    return false;
}

struct OneDirection {
    std::vector<std::uint8_t> bases;
    std::int32_t offset = 0;
    bool expanded = false;
};

OneDirection uncollapse_one_direction(
    const std::vector<std::uint8_t>& input_bases,
    const std::vector<std::uint8_t>& input_reference,
    const std::size_t threshold,
    const bool partial_mode) {
    OneDirection result;
    result.bases = input_bases;
    if (input_bases.empty() || input_reference.empty()) return result;
    const fastgatk::kernels::SmithWatermanParameters parameters{
        200, -150, -260, -11};
    const auto score = fastgatk::kernels::smith_waterman_score_kokkos(
        {fastgatk::kernels::SmithWatermanRequest{input_bases, input_reference}},
        parameters, fastgatk::kernels::SmithWatermanOverhangStrategy::Indel);
    if (score.scores.size() != 1)
        return result;
    const auto alignment = fastgatk::kernels::smith_waterman_align_reference(
        input_bases.data(), input_bases.size(), input_reference.data(),
        input_reference.size(), parameters,
        fastgatk::kernels::SmithWatermanOverhangStrategy::Indel);
    if (alignment.score != score.scores.front())
        throw std::runtime_error(
            "NUMERICAL_CONTRACT_FAILURE: flow HMER SW score/traceback mismatch");
    result.offset = alignment.alignment_offset;
    if (alignment.cigar.empty()) return result;

    std::size_t output_length = input_bases.size();
    std::size_t number = 0;
    for (const auto character : alignment.cigar) {
        if (std::isdigit(static_cast<unsigned char>(character))) {
            const auto digit = static_cast<std::size_t>(character - '0');
            if (number > (std::numeric_limits<std::size_t>::max() - digit) / 10)
                return OneDirection{input_bases, alignment.alignment_offset, false};
            number = number * 10 + digit;
            continue;
        }
        if (number == 0) return OneDirection{input_bases, alignment.alignment_offset, false};
        if (character == 'D') {
            if (output_length > std::numeric_limits<std::size_t>::max() - number)
                return OneDirection{input_bases, alignment.alignment_offset, false};
            output_length += number;
        }
        number = 0;
    }
    if (number != 0) return OneDirection{input_bases, alignment.alignment_offset, false};

    std::vector<std::uint8_t> output(output_length);
    std::size_t bases_offset = alignment.alignment_offset < 0
        ? 0U : static_cast<std::size_t>(alignment.alignment_offset);
    if (bases_offset > input_bases.size())
        return OneDirection{input_bases, alignment.alignment_offset, false};
    std::size_t reference_offset = 0;
    std::size_t output_offset = 0;
    number = 0;
    const auto fail = [&]() {
        return OneDirection{input_bases, alignment.alignment_offset, false};
    };
    for (const auto character : alignment.cigar) {
        if (std::isdigit(static_cast<unsigned char>(character))) {
            number = number * 10 + static_cast<std::size_t>(character - '0');
            continue;
        }
        if (number == 0) return fail();
        const auto consumes_read = character == 'M' || character == '=' ||
            character == 'X' || character == 'I' || character == 'S';
        const auto consumes_reference = character == 'M' || character == '=' ||
            character == 'X' || character == 'D';
        if (character != 'D') {
            if (consumes_read &&
                (number > input_bases.size() - bases_offset ||
                 number > output.size() - output_offset))
                return fail();
            if (consumes_read) {
                std::copy_n(input_bases.begin() + static_cast<std::ptrdiff_t>(bases_offset),
                            number, output.begin() + static_cast<std::ptrdiff_t>(output_offset));
                bases_offset += number;
                output_offset += number;
            }
        } else {
            bool incoming_hmer = false;
            const auto forward_end = std::min(input_bases.size(),
                bases_offset > input_bases.size() - std::min(threshold, input_bases.size())
                    ? input_bases.size() : bases_offset + threshold);
            std::vector<std::uint8_t> forward(
                input_bases.begin() + static_cast<std::ptrdiff_t>(bases_offset),
                input_bases.begin() + static_cast<std::ptrdiff_t>(forward_end));
            incoming_hmer = needs_collapsing(forward, threshold - 1);
            const auto backward_begin = bases_offset > threshold
                ? bases_offset - threshold : 0U;
            std::vector<std::uint8_t> backward(
                input_bases.begin() + static_cast<std::ptrdiff_t>(backward_begin),
                input_bases.begin() + static_cast<std::ptrdiff_t>(bases_offset));
            incoming_hmer = incoming_hmer || needs_collapsing(backward, threshold - 1);
            if (incoming_hmer && reference_offset < input_reference.size()) {
                const auto fill_from_left = reference_offset >= threshold &&
                    same_base(input_reference,
                              static_cast<std::int64_t>(reference_offset) -
                                  static_cast<std::int64_t>(threshold),
                              input_reference[reference_offset], threshold);
                const auto right_start = number > input_reference.size() - reference_offset
                    ? input_reference.size() : reference_offset + number;
                const auto fill_from_right = right_start > 0 &&
                    right_start <= input_reference.size() &&
                    same_base(input_reference,
                              static_cast<std::int64_t>(right_start) - 1,
                              input_reference[right_start - 1], threshold);
                if (fill_from_left || fill_from_right) {
                    const auto base = fill_from_left ? input_reference[reference_offset]
                                                     : input_reference[right_start - 1];
                    for (std::size_t fill = 0; fill < number; ++fill) {
                        const auto ref_index = fill_from_left
                            ? reference_offset + fill : right_start - 1 - fill;
                        if (partial_mode &&
                            (ref_index >= input_reference.size() ||
                             input_reference[ref_index] != base))
                            break;
                        if (output_offset >= output.size()) return fail();
                        output[output_offset++] = base;
                        result.expanded = true;
                    }
                }
            }
        }
        if (consumes_reference) {
            if (number > input_reference.size() - reference_offset) return fail();
            reference_offset += number;
        }
        number = 0;
    }
    output.resize(output_offset);
    result.bases = std::move(output);
    return result;
}

}  // namespace

FlowHmerUncollapseResult uncollapse_flow_hmers(
    const std::vector<std::uint8_t>& bases,
    const std::vector<std::uint8_t>& reference,
    const std::size_t hmer_size_threshold,
    const bool partial_mode,
    const bool limit_to_hmer_size_threshold) {
    FlowHmerUncollapseResult identity;
    identity.bases = bases;
    if (bases.empty() || reference.empty() || hmer_size_threshold == 0)
        return identity;
    if (hmer_size_threshold > std::numeric_limits<std::size_t>::max() / 2)
        throw std::invalid_argument("flow hmer uncollapse threshold is too large");
    auto forward = uncollapse_one_direction(bases, reference,
                                            hmer_size_threshold, partial_mode);
    auto reverse = uncollapse_one_direction(reverse_complement(bases),
                                             reverse_complement(reference),
                                             hmer_size_threshold, partial_mode);
    const bool use_reverse = reverse.bases.size() > forward.bases.size();
    auto selected = use_reverse ? std::move(reverse) : std::move(forward);
    if (use_reverse) selected.bases = reverse_complement(selected.bases);
    if (limit_to_hmer_size_threshold)
        selected.bases = fastgatk::io::collapse_flow_homopolymers(
            selected.bases, hmer_size_threshold);
    return FlowHmerUncollapseResult{
        std::move(selected.bases), selected.offset, selected.expanded};
}

}  // namespace fastgatk::calling
