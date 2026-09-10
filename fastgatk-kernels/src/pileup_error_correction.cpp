#include "fastgatk/kernels/pileup_error_correction.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

constexpr std::uint8_t kMatch = 0;
constexpr std::uint8_t kInsertion = 1;
constexpr std::uint8_t kDeletion = 2;
constexpr std::uint8_t kReferenceSkip = 3;
constexpr std::uint8_t kSoftClip = 4;
constexpr std::uint8_t kSequenceMatch = 7;
constexpr std::uint8_t kSequenceMismatch = 8;

KOKKOS_INLINE_FUNCTION
std::uint8_t encode_base(const std::uint8_t base) {
    switch (base) {
        case 0: case 'A': case 'a': return 0;
        case 1: case 'C': case 'c': return 1;
        case 2: case 'G': case 'g': return 2;
        case 3: case 'T': case 't': return 3;
        default: return 4;
    }
}

KOKKOS_INLINE_FUNCTION
std::uint8_t decode_base(const std::uint8_t base) {
    return base == 0 ? static_cast<std::uint8_t>('A') :
           base == 1 ? static_cast<std::uint8_t>('C') :
           base == 2 ? static_cast<std::uint8_t>('G') :
           base == 3 ? static_cast<std::uint8_t>('T') :
           static_cast<std::uint8_t>('N');
}

struct LocusKey {
    std::int32_t tid = -1;
    std::int32_t position = -1;
    bool operator<(const LocusKey& other) const noexcept {
        return tid < other.tid || (tid == other.tid && position < other.position);
    }
};

struct Observation {
    std::size_t absolute = 0;
    std::size_t record = 0;
    std::uint32_t offset = 0;
    std::uint8_t base = 4;
    std::uint8_t quality = 0;
    bool eligible = false;
};

struct LocusData {
    std::array<std::uint32_t, 4> counts{};
    std::vector<Observation> observations;
};

struct ReadEdit {
    std::uint32_t offset = 0;
    std::uint8_t target = 4;
};

double integer_digamma_difference(const std::size_t n_ref, const std::size_t n_alt) {
    // psi(n+1) = H_n - gamma; the Euler constant cancels in the difference.
    const auto harmonic = [](const std::size_t n) {
        double value = 0.0;
        for (std::size_t i = 1; i <= n; ++i) value += 1.0 / static_cast<double>(i);
        return value;
    };
    return harmonic(n_ref) - harmonic(n_alt);
}

double bernoulli_entropy(const double probability) {
    if (probability <= 0.0 || probability >= 1.0) return 0.0;
    return -probability * std::log(probability) -
           (1.0 - probability) * std::log1p(-probability);
}

double pileup_log_likelihood_ratio(const std::size_t n_ref,
                                   const std::vector<std::uint8_t>& alt_qualities) {
    const auto n_alt = alt_qualities.size();
    const auto n = n_ref + n_alt;
    if (n == 0 || n_alt == 0) return -std::numeric_limits<double>::infinity();
    const auto f_tilde_ratio = std::exp(integer_digamma_difference(n_ref, n_alt));
    double read_sum = 0.0;
    for (const auto quality : alt_qualities) {
        const auto epsilon = std::pow(10.0, -static_cast<double>(quality) / 10.0);
        const auto denominator = 1.0 - epsilon + epsilon * f_tilde_ratio;
        const auto z_bar_alt = denominator == 0.0 ? 0.0 :
            (1.0 - epsilon) / denominator;
        const auto log_epsilon = std::log(epsilon);
        const auto log_one_minus_epsilon = std::log1p(-epsilon);
        read_sum += z_bar_alt * (log_one_minus_epsilon - log_epsilon) +
                    bernoulli_entropy(z_bar_alt);
    }
    // Flat Beta(1,1) prior: -log(n+1) - log(C(n,nAlt)).  lgamma is the
    // same real-valued formulation used by Apache Commons Math's
    // CombinatoricsUtils.binomialCoefficientLog.
    const auto beta_entropy = -std::log(static_cast<double>(n + 1)) -
        (std::lgamma(static_cast<double>(n + 1)) -
         std::lgamma(static_cast<double>(n_alt + 1)) -
         std::lgamma(static_cast<double>(n_ref + 1)));
    return beta_entropy + read_sum;
}

bool is_match_like(const std::uint8_t code) {
    return code == kMatch || code == kSequenceMatch || code == kSequenceMismatch;
}

bool projects_adjacent_to_indel_or_softclip(const std::vector<std::uint8_t>& codes,
                                            const std::size_t index) {
    const auto structural = [](const std::uint8_t code) {
        return code == kInsertion || code == kDeletion || code == kSoftClip;
    };
    return (index > 0 && structural(codes[index - 1])) ||
           (index + 1 < codes.size() && structural(codes[index + 1]));
}

}  // namespace

PileupReadErrorCorrectionResult correct_reads_by_pileup_kokkos(
    const PileupReadErrorCorrectionInput& input,
    PileupReadErrorCorrectionOptions options) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (!std::isfinite(options.log_odds_threshold) &&
        options.log_odds_threshold != -std::numeric_limits<double>::infinity())
        throw std::invalid_argument("pileup error-correction threshold must be finite or -infinity");
    if (options.indel_span == 0 || options.indel_mismatches == 0 ||
        options.corrected_base_quality < 2)
        throw std::invalid_argument("pileup error-correction options are invalid");
    if (input.offsets.empty() || input.offsets.front() != 0 ||
        input.offsets.back() != input.bases.size() ||
        input.qualities.size() != input.bases.size() ||
        input.tids.size() != input.offsets.size() - 1 ||
        input.positions.size() != input.offsets.size() - 1)
        throw std::invalid_argument("pileup error-correction input layout is malformed");
    for (std::size_t i = 1; i < input.offsets.size(); ++i)
        if (input.offsets[i] < input.offsets[i - 1])
            throw std::invalid_argument("pileup error-correction offsets are not monotonic");
    const auto records = input.offsets.size() - 1;
    if (!input.cigar_offsets.empty() &&
        (input.cigar_offsets.size() != records + 1 ||
         input.cigar_offsets.back() != input.cigar_ops.size()))
        throw std::invalid_argument("pileup error-correction CIGAR layout is malformed");

    PileupReadErrorCorrectionResult result;
    result.input_reads = records;
    result.execution_space = ExecSpace::name();
    result.bases = input.bases;
    result.qualities = input.qualities;
    result.used = options.log_odds_threshold != -std::numeric_limits<double>::infinity();
    if (!result.used) return result;

    std::map<LocusKey, LocusData> loci;
    for (std::size_t record = 0; record < records; ++record) {
        const auto read_begin = input.offsets[record];
        const auto read_end = input.offsets[record + 1];
        std::vector<std::uint8_t> codes;
        std::vector<std::uint32_t> lengths;
        if (!input.cigar_offsets.empty()) {
            const auto cigar_begin = input.cigar_offsets[record];
            const auto cigar_end = input.cigar_offsets[record + 1];
            codes.reserve(cigar_end - cigar_begin);
            lengths.reserve(cigar_end - cigar_begin);
            for (std::size_t index = cigar_begin; index < cigar_end; ++index) {
                const auto packed = input.cigar_ops[index];
                const auto code = static_cast<std::uint8_t>(packed & 0x0fU);
                const auto length = packed >> 4U;
                if (length == 0 || code > kSequenceMismatch)
                    throw std::invalid_argument("pileup error-correction CIGAR op is invalid");
                codes.push_back(code);
                lengths.push_back(length);
            }
        }
        std::int64_t reference_cursor = input.positions[record];
        std::size_t read_cursor = read_begin;
        if (codes.empty()) {
            for (std::size_t absolute = read_begin; absolute < read_end; ++absolute) {
                const auto base = encode_base(input.bases[absolute]);
                if (base >= 4 || reference_cursor < 0 ||
                    reference_cursor > std::numeric_limits<std::int32_t>::max()) {
                    ++reference_cursor;
                    continue;
                }
                auto& locus = loci[LocusKey{input.tids[record],
                    static_cast<std::int32_t>(reference_cursor)}];
                ++locus.counts[base];
                locus.observations.push_back(Observation{
                    absolute, record, static_cast<std::uint32_t>(absolute - read_begin), base,
                    input.qualities[absolute], true});
                ++reference_cursor;
            }
            continue;
        }
        for (std::size_t segment = 0; segment < codes.size(); ++segment) {
            const auto code = codes[segment];
            const auto length = static_cast<std::size_t>(lengths[segment]);
            if (is_match_like(code)) {
                for (std::size_t offset = 0; offset < length; ++offset) {
                    if (read_cursor >= read_end || reference_cursor < 0 ||
                        reference_cursor > std::numeric_limits<std::int32_t>::max()) break;
                    const auto absolute = read_cursor++;
                    const auto base = encode_base(input.bases[absolute]);
                    if (base < 4) {
                        const bool eligible = !projects_adjacent_to_indel_or_softclip(codes, segment);
                        auto& locus = loci[LocusKey{input.tids[record],
                            static_cast<std::int32_t>(reference_cursor)}];
                        ++locus.counts[base];
                        locus.observations.push_back(Observation{
                            absolute, record, static_cast<std::uint32_t>(absolute - read_begin), base,
                            input.qualities[absolute], eligible});
                        if (!eligible) ++result.skipped_indel_adjacent_bases;
                    }
                    ++reference_cursor;
                }
            } else {
                if (code == kInsertion || code == kSoftClip) read_cursor += length;
                if (code == kDeletion || code == kReferenceSkip) reference_cursor += length;
            }
            if (read_cursor > read_end)
                throw std::invalid_argument("pileup error-correction CIGAR consumes too many bases");
        }
    }
    result.loci = loci.size();

    std::vector<std::uint8_t> target(input.bases.size(), 4);
    std::vector<std::vector<ReadEdit>> edits(records);
    for (const auto& entry : loci) {
        const auto& locus = entry.second;
        std::uint8_t plurality = 0;
        for (std::uint8_t base = 1; base < 4; ++base)
            if (locus.counts[base] > locus.counts[plurality]) plurality = base;
        std::vector<std::uint8_t> alt_qualities;
        alt_qualities.reserve(locus.observations.size());
        std::size_t ref_count = 0;
        for (const auto& observation : locus.observations) {
            if (observation.base == plurality) ++ref_count;
            else alt_qualities.push_back(observation.quality);
        }
        if (alt_qualities.empty() ||
            pileup_log_likelihood_ratio(ref_count, alt_qualities) >= options.log_odds_threshold)
            continue;
        bool locus_corrected = false;
        for (const auto& observation : locus.observations) {
            if (observation.base == plurality || !observation.eligible) continue;
            edits[observation.record].push_back(ReadEdit{observation.offset, plurality});
            locus_corrected = true;
        }
        if (locus_corrected) ++result.corrected_loci;
    }

    for (auto& read_edits : edits) {
        std::sort(read_edits.begin(), read_edits.end(),
                  [](const auto& left, const auto& right) { return left.offset < right.offset; });
        std::ptrdiff_t first_edit = 0;
        const auto mismatch_count = static_cast<std::size_t>(options.indel_mismatches);
        const auto span = static_cast<std::size_t>(options.indel_span);
        for (std::size_t n = 0; n + mismatch_count < read_edits.size(); ++n) {
            if (read_edits[n + mismatch_count - 1].offset - read_edits[n].offset < span)
                first_edit = static_cast<std::ptrdiff_t>(n + mismatch_count);
        }
        std::ptrdiff_t last_edit = read_edits.empty()
            ? -1 : static_cast<std::ptrdiff_t>(read_edits.size() - 1);
        if (read_edits.size() >= mismatch_count) {
            for (std::ptrdiff_t n = static_cast<std::ptrdiff_t>(read_edits.size()) - 1;
                 n >= static_cast<std::ptrdiff_t>(mismatch_count - 1); --n) {
                const auto current = static_cast<std::size_t>(n);
                const auto earlier = static_cast<std::size_t>(n - mismatch_count + 1);
                if (read_edits[current].offset - read_edits[earlier].offset < span)
                    last_edit = n - static_cast<std::ptrdiff_t>(mismatch_count);
            }
        }
        if (read_edits.empty() || first_edit > last_edit) continue;
        for (auto index = first_edit; index <= last_edit; ++index) {
            const auto absolute = input.offsets[&read_edits - edits.data()] + read_edits[index].offset;
            target[absolute] = read_edits[index].target;
        }
    }

    fastgatk::core::HostBatch host("pileup-error-correction-v1");
    host.records = input.bases.size();
    host.bytes = input.bases.size() * 3U + input.offsets.size() * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("pileup-error-correction");
    plan.begin_prepare(host);
    Kokkos::View<std::uint8_t*, MemorySpace> bases("pileup_correction_bases", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> corrected("pileup_correction_corrected", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> changed("pileup_correction_changed", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> target_view("pileup_correction_target", input.bases.size());
    auto host_bases = Kokkos::create_mirror_view(bases);
    auto host_target = Kokkos::create_mirror_view(target_view);
    for (std::size_t index = 0; index < input.bases.size(); ++index) {
        host_bases(index) = input.bases[index];
        host_target(index) = target[index];
    }
    Kokkos::deep_copy(bases, host_bases);
    Kokkos::deep_copy(target_view, host_target);
    fastgatk::core::DeviceBatch<ExecSpace> device(input.bases.size());
    device.bind("bases", bases);
    device.bind("corrected", corrected);
    device.bind("changed", changed);
    device.bind("target", target_view);
    ExecSpace().fence();
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    Kokkos::parallel_for("pileup_error_correction", Kokkos::RangePolicy<ExecSpace>(0, input.bases.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            corrected(index) = bases(index);
            changed(index) = 0;
            const auto desired = target_view(index);
            if (desired < 4 && encode_base(bases(index)) != desired) {
                corrected(index) = decode_base(desired);
                changed(index) = 1;
            }
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_corrected = Kokkos::create_mirror_view(corrected);
    auto host_changed = Kokkos::create_mirror_view(changed);
    Kokkos::deep_copy(host_corrected, corrected);
    Kokkos::deep_copy(host_changed, changed);
    std::vector<std::uint8_t> read_changed(records, 0);
    for (std::size_t index = 0; index < input.bases.size(); ++index) {
        result.bases[index] = host_corrected(index);
        if (host_changed(index) == 0) continue;
        result.qualities[index] = options.corrected_base_quality;
        ++result.corrected_bases;
        const auto record = static_cast<std::size_t>(
            std::upper_bound(input.offsets.begin(), input.offsets.end(),
                             static_cast<std::uint32_t>(index)) - input.offsets.begin() - 1);
        if (record < read_changed.size()) read_changed[record] = 1;
    }
    result.corrected_reads = std::count(read_changed.begin(), read_changed.end(), std::uint8_t{1});
    return result;
}

}  // namespace fastgatk::kernels
