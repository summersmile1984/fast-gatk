#include "fastgatk/kernels/read_error_correction.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

KOKKOS_INLINE_FUNCTION
std::uint8_t encode_base_device(const std::uint8_t base) {
    switch (base) {
        case 0: case 'A': case 'a': return 0;
        case 1: case 'C': case 'c': return 1;
        case 2: case 'G': case 'g': return 2;
        case 3: case 'T': case 't': return 3;
        default: return 4;
    }
}

KOKKOS_INLINE_FUNCTION
std::uint8_t decode_base_device(const std::uint8_t base) {
    switch (base) {
        case 0: return static_cast<std::uint8_t>('A');
        case 1: return static_cast<std::uint8_t>('C');
        case 2: return static_cast<std::uint8_t>('G');
        case 3: return static_cast<std::uint8_t>('T');
        default: return static_cast<std::uint8_t>('N');
    }
}

std::uint8_t encode_base_host(const std::uint8_t base) {
    switch (base) {
        case 0: return 0;
        case 1: return 1;
        case 2: return 2;
        case 3: return 3;
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

bool make_kmer_key(const std::vector<std::uint8_t>& bases,
                  const std::size_t begin, const std::size_t length,
                  std::uint64_t& key) {
    key = 0;
    for (std::size_t offset = 0; offset < length; ++offset) {
        const auto code = encode_base_host(bases[begin + offset]);
        if (code >= 4) return false;
        key = (key << 2U) | code;
    }
    return true;
}

}  // namespace

ReadErrorCorrectionResult correct_read_errors_kokkos(
    const ReadErrorCorrectionInput& input,
    ReadErrorCorrectionOptions options) {
    if (!Kokkos::is_initialized())
        throw std::runtime_error("Kokkos is not initialized");
    if (input.offsets.empty() || input.offsets.front() != 0 ||
        input.offsets.back() != input.bases.size())
        throw std::invalid_argument("read error correction offsets are malformed");
    for (std::size_t index = 1; index < input.offsets.size(); ++index)
        if (input.offsets[index] < input.offsets[index - 1])
            throw std::invalid_argument("read error correction offsets are not monotonic");
    if (!input.qualities.empty() && input.qualities.size() != input.bases.size())
        throw std::invalid_argument("read error correction qualities do not match bases");
    if (options.kmer_length == 0 || options.kmer_length > 31 ||
        options.min_observations_for_kmer_to_be_solid == 0 ||
        options.max_mismatches == 0 || options.max_mismatches > 2 ||
        options.max_observations_for_kmer_to_be_correctable == 0 ||
        options.corrected_base_quality < 2)
        throw std::invalid_argument(
            "read error correction requires kmer length [1,31], positive solid threshold, "
            "max mismatches [1,2], positive sparse-kmer floor and corrected quality >=2");

    ReadErrorCorrectionResult result;
    result.used = true;
    result.input_reads = input.offsets.size() - 1;
    result.execution_space = ExecSpace::name();
    result.bases = input.bases;
    result.qualities = input.qualities;
    if (result.qualities.empty()) result.qualities.assign(result.bases.size(), 0);

    // The exact key table is immutable during correction.  Host sorting keeps
    // key IDs stable across execution spaces; the device only performs binary
    // lookup and per-base candidate scoring.
    std::vector<std::uint64_t> observed_keys;
    for (std::size_t read = 0; read + 1 < input.offsets.size(); ++read) {
        const auto begin = input.offsets[read];
        const auto end = input.offsets[read + 1];
        if (end - begin < options.kmer_length) continue;
        for (std::size_t start = begin; start + options.kmer_length <= end; ++start) {
            std::uint64_t key = 0;
            if (make_kmer_key(input.bases, start, options.kmer_length, key))
                observed_keys.push_back(key);
        }
    }
    std::sort(observed_keys.begin(), observed_keys.end());
    observed_keys.erase(std::unique(observed_keys.begin(), observed_keys.end()), observed_keys.end());
    if (observed_keys.empty()) return result;

    std::vector<std::uint32_t> observed_counts(observed_keys.size(), 0);
    for (std::size_t read = 0; read + 1 < input.offsets.size(); ++read) {
        const auto begin = input.offsets[read];
        const auto end = input.offsets[read + 1];
        if (end - begin < options.kmer_length) continue;
        for (std::size_t start = begin; start + options.kmer_length <= end; ++start) {
            std::uint64_t key = 0;
            if (!make_kmer_key(input.bases, start, options.kmer_length, key)) continue;
            const auto iter = std::lower_bound(observed_keys.begin(), observed_keys.end(), key);
            if (iter != observed_keys.end() && *iter == key)
                ++observed_counts[static_cast<std::size_t>(iter - observed_keys.begin())];
        }
    }
    std::unordered_map<std::uint64_t, std::uint32_t> count_by_key;
    count_by_key.reserve(observed_keys.size() * 2U);
    for (std::size_t index = 0; index < observed_keys.size(); ++index)
        count_by_key.emplace(observed_keys[index], observed_counts[index]);
    std::unordered_map<std::uint64_t, std::uint64_t> correction_by_key;
    correction_by_key.reserve(observed_keys.size());
    const auto solid_threshold = options.min_observations_for_kmer_to_be_solid;
    const auto max_mismatches = options.max_mismatches;
    const auto max_sparse_observations = options.max_observations_for_kmer_to_be_correctable;
    for (std::size_t index = 0; index < observed_keys.size(); ++index) {
        const auto key = observed_keys[index];
        const auto count = observed_counts[index];
        if (count >= solid_threshold) {
            ++result.solid_kmers;
            continue;
        }
        if (count > max_sparse_observations) {
            ++result.uncorrectable_kmers;
            continue;
        }

        std::uint32_t best_distance = max_mismatches + 1U;
        std::uint64_t best_key = 0;
        bool found = false;
        const auto consider = [&](const std::uint64_t candidate,
                                  const std::uint32_t distance) {
            const auto iter = count_by_key.find(candidate);
            if (iter == count_by_key.end() || iter->second < solid_threshold) return;
            if (!found || distance < best_distance ||
                (distance == best_distance && candidate < best_key)) {
                found = true;
                best_distance = distance;
                best_key = candidate;
            }
        };
        // Enumerate Hamming-distance-one neighbors directly.  The key is
        // 2-bit encoded from left to right, so each position has exactly
        // three substitutions and no temporary strings are required.
        for (std::uint32_t position = 0; position < options.kmer_length; ++position) {
            const auto shift = 2U * (options.kmer_length - position - 1U);
            const auto original = static_cast<std::uint64_t>((key >> shift) & 3U);
            for (std::uint64_t base = 0; base < 4; ++base) {
                if (base == original) continue;
                consider(key ^ ((original ^ base) << shift), 1U);
            }
        }
        if (max_mismatches >= 2U) {
            for (std::uint32_t first = 0; first < options.kmer_length; ++first) {
                const auto first_shift = 2U * (options.kmer_length - first - 1U);
                const auto first_original = static_cast<std::uint64_t>((key >> first_shift) & 3U);
                for (std::uint32_t second = first + 1U; second < options.kmer_length; ++second) {
                    const auto second_shift = 2U * (options.kmer_length - second - 1U);
                    const auto second_original = static_cast<std::uint64_t>((key >> second_shift) & 3U);
                    for (std::uint64_t first_base = 0; first_base < 4; ++first_base) {
                        if (first_base == first_original) continue;
                        for (std::uint64_t second_base = 0; second_base < 4; ++second_base) {
                            if (second_base == second_original) continue;
                            const auto candidate = key ^
                                ((first_original ^ first_base) << first_shift) ^
                                ((second_original ^ second_base) << second_shift);
                            consider(candidate, 2U);
                        }
                    }
                }
            }
        }
        if (found) {
            correction_by_key.emplace(key, best_key);
            ++result.corrected_kmers;
        } else {
            ++result.uncorrectable_kmers;
        }
    }

    // Build the Java CorrectionSet equivalent on Host. Each base receives
    // all edits proposed by overlapping corrected kmers; conflicting edits
    // are rejected by strict consensus rather than resolved by support.
    std::vector<std::uint8_t> candidate_base(input.bases.size(), 0);
    std::vector<std::uint8_t> candidate_count(input.bases.size(), 0);
    std::vector<std::uint8_t> candidate_conflict(input.bases.size(), 0);
    for (std::size_t read = 0; read + 1 < input.offsets.size(); ++read) {
        const auto begin = input.offsets[read];
        const auto end = input.offsets[read + 1];
        if (end - begin < options.kmer_length) continue;
        for (std::size_t start = begin; start + options.kmer_length <= end; ++start) {
            std::uint64_t key = 0;
            if (!make_kmer_key(input.bases, start, options.kmer_length, key)) continue;
            const auto correction = correction_by_key.find(key);
            if (correction == correction_by_key.end() || correction->second == key) continue;
            for (std::uint32_t position = 0; position < options.kmer_length; ++position) {
                const auto shift = 2U * (options.kmer_length - position - 1U);
                const auto from = static_cast<std::uint8_t>((key >> shift) & 3U);
                const auto to = static_cast<std::uint8_t>((correction->second >> shift) & 3U);
                if (from == to) continue;
                const auto absolute = start + position;
                if (candidate_count[absolute] == 0) {
                    candidate_base[absolute] = to;
                } else if (candidate_base[absolute] != to) {
                    candidate_conflict[absolute] = 1;
                }
                if (candidate_count[absolute] != std::numeric_limits<std::uint8_t>::max())
                    ++candidate_count[absolute];
            }
        }
    }

    fastgatk::core::HostBatch host("read-error-correction-v1");
    host.records = input.bases.size();
    host.bytes = input.bases.size() + input.offsets.size() * sizeof(std::uint32_t) +
                 observed_keys.size() * (sizeof(std::uint64_t) + sizeof(std::uint32_t)) +
                 input.bases.size() * 3U;
    fastgatk::core::KernelPlan<ExecSpace> plan("read-error-correction");
    plan.begin_prepare(host);
    Kokkos::View<std::uint8_t*, MemorySpace> bases("error_correction_bases", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> corrected_bases("error_correction_corrected_bases", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> changed("error_correction_changed", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> candidate_base_view(
        "error_correction_candidate_base", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> candidate_count_view(
        "error_correction_candidate_count", input.bases.size());
    Kokkos::View<std::uint8_t*, MemorySpace> candidate_conflict_view(
        "error_correction_candidate_conflict", input.bases.size());
    auto host_bases = Kokkos::create_mirror_view(bases);
    auto host_candidate_base = Kokkos::create_mirror_view(candidate_base_view);
    auto host_candidate_count = Kokkos::create_mirror_view(candidate_count_view);
    auto host_candidate_conflict = Kokkos::create_mirror_view(candidate_conflict_view);
    for (std::size_t index = 0; index < input.bases.size(); ++index) host_bases(index) = input.bases[index];
    for (std::size_t index = 0; index < input.bases.size(); ++index) {
        host_candidate_base(index) = candidate_base[index];
        host_candidate_count(index) = candidate_count[index];
        host_candidate_conflict(index) = candidate_conflict[index];
    }
    Kokkos::deep_copy(bases, host_bases);
    Kokkos::deep_copy(candidate_base_view, host_candidate_base);
    Kokkos::deep_copy(candidate_count_view, host_candidate_count);
    Kokkos::deep_copy(candidate_conflict_view, host_candidate_conflict);
    fastgatk::core::DeviceBatch<ExecSpace> device(input.bases.size());
    device.bind("bases", bases);
    device.bind("corrected_bases", corrected_bases);
    device.bind("changed", changed);
    device.bind("candidate_base", candidate_base_view);
    device.bind("candidate_count", candidate_count_view);
    device.bind("candidate_conflict", candidate_conflict_view);
    ExecSpace().fence();
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    Kokkos::parallel_for("read_error_correction", Kokkos::RangePolicy<ExecSpace>(0, input.bases.size()),
        KOKKOS_LAMBDA(const std::size_t position) {
            const auto current = encode_base_device(bases(position));
            corrected_bases(position) = bases(position);
            changed(position) = 0;
            if (current >= 4) return;
            if (candidate_count_view(position) == 0 || candidate_conflict_view(position) != 0)
                return;
            const auto corrected = candidate_base_view(position);
            if (corrected == current) return;
            corrected_bases(position) = decode_base_device(corrected);
            changed(position) = 1;
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_corrected = Kokkos::create_mirror_view(corrected_bases);
    auto host_changed = Kokkos::create_mirror_view(changed);
    Kokkos::deep_copy(host_corrected, corrected_bases);
    Kokkos::deep_copy(host_changed, changed);
    std::vector<std::uint8_t> per_read_changed(result.input_reads, 0);
    for (std::size_t index = 0; index < input.bases.size(); ++index) {
        result.bases[index] = host_corrected(index);
        if (host_changed(index) == 0) continue;
        result.qualities[index] = options.corrected_base_quality;
        ++result.corrected_bases;
        std::size_t read = 0;
        while (read + 1 < input.offsets.size() &&
               !(index >= input.offsets[read] && index < input.offsets[read + 1])) ++read;
        if (read < per_read_changed.size()) per_read_changed[read] = 1;
    }
    result.corrected_reads = std::count(per_read_changed.begin(), per_read_changed.end(), std::uint8_t{1});
    return result;
}

}  // namespace fastgatk::kernels
