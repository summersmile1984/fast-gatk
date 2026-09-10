#include "fastgatk/kernels/bqsr.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace fastgatk::kernels {
namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

}  // namespace

BqsrCountResult count_bqsr_quality_kokkos(const BqsrObservationBatch& input) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (input.qualities.size() != input.mismatches.size())
        throw std::invalid_argument("BQSR observation arrays differ in length");

    BqsrCountResult result;
    result.bins.resize(94);
    result.observations = input.qualities.size();
    result.execution_space = ExecSpace::name();
    result.execution_policy = "TeamPolicy";
    if (input.qualities.empty()) return result;

    fastgatk::core::HostBatch host("bqsr-quality-v1");
    host.records = input.qualities.size();
    host.bytes = input.qualities.size() * 2 * sizeof(std::uint8_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("bqsr-quality");
    plan.begin_prepare(host);

    Kokkos::View<std::uint8_t*, MemorySpace> qualities("bqsr_qualities", input.qualities.size());
    Kokkos::View<std::uint8_t*, MemorySpace> mismatches("bqsr_mismatches", input.mismatches.size());
    auto host_qualities = Kokkos::create_mirror_view(qualities);
    auto host_mismatches = Kokkos::create_mirror_view(mismatches);
    for (std::size_t index = 0; index < input.qualities.size(); ++index) {
        host_qualities(index) = std::min<std::uint8_t>(input.qualities[index], 93);
        host_mismatches(index) = input.mismatches[index] == 0 ? 0 : 1;
    }
    Kokkos::deep_copy(qualities, host_qualities);
    Kokkos::deep_copy(mismatches, host_mismatches);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(input.qualities.size());
    device.bind("qualities", qualities);
    device.bind("mismatches", mismatches);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;

    plan.begin_execute();
    // Use a deterministic team-local histogram rather than a global bin
    // atomic.  Each team owns a contiguous observation tile; the Host sums
    // complete team rows in league order after the kernel, so integer counts
    // are backend-independent while high-core CPU/GPU contention is bounded.
    constexpr std::size_t observations_per_team = 256;
    constexpr std::size_t histogram_width = 188; // 94 counts + 94 errors
    const auto observation_count = input.qualities.size();
    const auto league_size = (observation_count + observations_per_team - 1) /
                             observations_per_team;
    Kokkos::View<std::uint64_t**, Kokkos::LayoutRight, MemorySpace> team_histogram(
        "bqsr_team_histogram", league_size, histogram_width);
    Kokkos::deep_copy(team_histogram, std::uint64_t{0});
    using TeamMember = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
    const Kokkos::TeamPolicy<ExecSpace> policy(league_size, Kokkos::AUTO());
    Kokkos::parallel_for("bqsr_quality_team_histogram", policy,
        KOKKOS_LAMBDA(const TeamMember& team) {
            const auto league = static_cast<std::size_t>(team.league_rank());
            const auto begin = league * observations_per_team;
            const auto end = (begin + observations_per_team < observation_count)
                ? begin + observations_per_team : observation_count;
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, begin, end),
                [&](const std::size_t index) {
                    const auto quality = static_cast<std::size_t>(qualities(index));
                    Kokkos::atomic_add(&team_histogram(league, quality), std::uint64_t{1});
                    if (mismatches(index) != 0)
                        Kokkos::atomic_add(&team_histogram(league, 94 + quality), std::uint64_t{1});
                });
        });
    ExecSpace().fence();
    plan.end_execute();
    result.execute_seconds = plan.telemetry().execute_seconds;
    auto host_team_histogram = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), team_histogram);
    std::array<std::uint64_t, histogram_width> histogram{};
    for (std::size_t league = 0; league < league_size; ++league)
        for (std::size_t bin = 0; bin < histogram_width; ++bin)
            histogram[bin] += host_team_histogram(league, bin);
    for (std::size_t quality = 0; quality < result.bins.size(); ++quality)
        result.bins[quality] = BqsrQualityCounts{histogram[quality], histogram[94 + quality]};
    return result;
}

BqsrCovariateCountResult count_bqsr_covariates_kokkos(
    const BqsrCovariateObservationBatch& input) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (input.covariate_ids.size() != input.mismatches.size())
        throw std::invalid_argument("BQSR covariate id and mismatch arrays differ in length");

    BqsrCovariateCountResult result;
    result.observations = input.covariate_ids.size();
    result.execution_space = ExecSpace::name();
    result.execution_policy = "TeamPolicy";
    if (input.covariate_ids.empty()) return result;

    const auto max_id = *std::max_element(input.covariate_ids.begin(),
                                          input.covariate_ids.end());
    // The Host contract assigns dense [0,n) ids for the current batch.  Fail
    // closed on a sparse/untrusted id rather than allocating an attacker-sized
    // device histogram for a single observation.
    if (static_cast<std::size_t>(max_id) >= input.covariate_ids.size())
        throw std::invalid_argument("BQSR covariate ids must be dense per batch");
    const auto covariate_count = static_cast<std::size_t>(max_id) + 1;
    result.covariates.resize(covariate_count);

    fastgatk::core::HostBatch host("bqsr-covariate-v1");
    host.records = input.covariate_ids.size();
    host.bytes = input.covariate_ids.size() *
                 (sizeof(std::uint32_t) + sizeof(std::uint8_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("bqsr-covariate");
    plan.begin_prepare(host);

    Kokkos::View<std::uint32_t*, MemorySpace> ids(
        "bqsr_covariate_ids", input.covariate_ids.size());
    Kokkos::View<std::uint8_t*, MemorySpace> mismatches(
        "bqsr_covariate_mismatches", input.mismatches.size());
    auto host_ids = Kokkos::create_mirror_view(ids);
    auto host_mismatches = Kokkos::create_mirror_view(mismatches);
    for (std::size_t index = 0; index < input.covariate_ids.size(); ++index) {
        host_ids(index) = input.covariate_ids[index];
        host_mismatches(index) = input.mismatches[index] == 0 ? 0 : 1;
    }
    Kokkos::deep_copy(ids, host_ids);
    Kokkos::deep_copy(mismatches, host_mismatches);
    ExecSpace().fence();

    // A dense per-team table avoids global atomic contention for normal BQSR
    // batches.  Keep the allocation bounded because a batch containing a
    // nearly-unique covariate key per observation would otherwise multiply
    // the key width by the number of teams.
    constexpr std::size_t observations_per_team = 256;
    constexpr std::uint64_t max_team_histogram_entries = 4ULL * 1024ULL * 1024ULL;
    const auto observation_count = input.covariate_ids.size();
    const auto league_size = (observation_count + observations_per_team - 1) /
                             observations_per_team;
    const auto histogram_width = covariate_count > std::numeric_limits<std::size_t>::max() / 2
        ? std::numeric_limits<std::size_t>::max() : covariate_count * 2;
    const bool team_histogram_fits = histogram_width != std::numeric_limits<std::size_t>::max() &&
        league_size != 0 && histogram_width <= std::numeric_limits<std::size_t>::max() / league_size &&
        static_cast<std::uint64_t>(histogram_width * league_size) <= max_team_histogram_entries;
    result.team_local_histogram = team_histogram_fits;
    result.execution_policy = team_histogram_fits ? "TeamPolicy" : "RangePolicy";
    result.workspace_bytes = static_cast<std::uint64_t>(
        (team_histogram_fits ? histogram_width * league_size : histogram_width) *
        sizeof(std::uint64_t));

    if (team_histogram_fits) {
        Kokkos::View<std::uint64_t**, Kokkos::LayoutRight, MemorySpace> team_counts(
            "bqsr_covariate_team_counts", league_size, histogram_width);
        Kokkos::deep_copy(team_counts, std::uint64_t{0});
        fastgatk::core::DeviceBatch<ExecSpace> device(input.covariate_ids.size());
        device.bind("covariate_ids", ids);
        device.bind("mismatches", mismatches);
        device.bind("team_counts", team_counts);
        plan.end_prepare(device);
        result.prepare_seconds = plan.telemetry().prepare_seconds;
        plan.begin_execute();
        using TeamMember = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
        Kokkos::parallel_for("bqsr_covariate_team_count",
            Kokkos::TeamPolicy<ExecSpace>(league_size, Kokkos::AUTO()),
            KOKKOS_LAMBDA(const TeamMember& team) {
                const auto league = static_cast<std::size_t>(team.league_rank());
                const auto begin = league * observations_per_team;
                const auto end = (begin + observations_per_team < observation_count)
                    ? begin + observations_per_team : observation_count;
                Kokkos::parallel_for(Kokkos::TeamThreadRange(team, begin, end),
                    [&](const std::size_t index) {
                        const auto id = static_cast<std::size_t>(ids(index));
                        Kokkos::atomic_add(&team_counts(league, id * 2), std::uint64_t{1});
                        if (mismatches(index) != 0)
                            Kokkos::atomic_add(&team_counts(league, id * 2 + 1), std::uint64_t{1});
                    });
            });
        ExecSpace().fence();
        plan.end_execute();
        result.execute_seconds = plan.telemetry().execute_seconds;
        const auto host_counts = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), team_counts);
        for (std::size_t league = 0; league < league_size; ++league)
            for (std::size_t id = 0; id < covariate_count; ++id) {
                result.covariates[id].count += host_counts(league, id * 2);
                result.covariates[id].mismatches += host_counts(league, id * 2 + 1);
            }
    } else {
        Kokkos::View<std::uint64_t**, Kokkos::LayoutRight, MemorySpace> counts(
            "bqsr_covariate_counts", covariate_count, 2);
        Kokkos::deep_copy(counts, std::uint64_t{0});
        fastgatk::core::DeviceBatch<ExecSpace> device(input.covariate_ids.size());
        device.bind("covariate_ids", ids);
        device.bind("mismatches", mismatches);
        device.bind("counts", counts);
        plan.end_prepare(device);
        result.prepare_seconds = plan.telemetry().prepare_seconds;
        plan.begin_execute();
        Kokkos::parallel_for("bqsr_covariate_count",
            Kokkos::RangePolicy<ExecSpace>(0, input.covariate_ids.size()),
            KOKKOS_LAMBDA(const std::size_t index) {
                const auto id = static_cast<std::size_t>(ids(index));
                Kokkos::atomic_add(&counts(id, 0), std::uint64_t{1});
                if (mismatches(index) != 0)
                    Kokkos::atomic_add(&counts(id, 1), std::uint64_t{1});
            });
        ExecSpace().fence();
        plan.end_execute();
        result.execute_seconds = plan.telemetry().execute_seconds;
        const auto host_counts = Kokkos::create_mirror_view_and_copy(
            Kokkos::HostSpace(), counts);
        for (std::size_t id = 0; id < covariate_count; ++id)
            result.covariates[id] = BqsrQualityCounts{host_counts(id, 0), host_counts(id, 1)};
    }
    return result;
}

BqsrQualityTransformResult apply_bqsr_quality_kokkos(
    const std::vector<std::uint8_t>& input,
    const std::vector<std::int16_t>& deltas) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (input.size() != deltas.size())
        throw std::invalid_argument("BQSR quality and delta arrays differ in length");

    BqsrQualityTransformResult result;
    result.adjusted.resize(input.size());
    result.execution_space = ExecSpace::name();
    result.execution_policy = "RangePolicy";
    if (input.empty()) return result;

    fastgatk::core::HostBatch host("bqsr-quality-transform-v1");
    host.records = input.size();
    host.bytes = input.size() * (sizeof(std::uint8_t) + sizeof(std::int16_t));
    fastgatk::core::KernelPlan<ExecSpace> plan("bqsr-quality-transform");
    plan.begin_prepare(host);
    Kokkos::View<std::uint8_t*, MemorySpace> qualities("bqsr_transform_qualities", input.size());
    Kokkos::View<std::int16_t*, MemorySpace> delta_view("bqsr_transform_deltas", input.size());
    Kokkos::View<std::uint8_t*, MemorySpace> adjusted("bqsr_transform_adjusted", input.size());
    auto host_qualities = Kokkos::create_mirror_view(qualities);
    auto host_deltas = Kokkos::create_mirror_view(delta_view);
    for (std::size_t index = 0; index < input.size(); ++index) {
        host_qualities(index) = std::min<std::uint8_t>(input[index], 93);
        host_deltas(index) = deltas[index];
    }
    Kokkos::deep_copy(qualities, host_qualities);
    Kokkos::deep_copy(delta_view, host_deltas);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(input.size());
    device.bind("qualities", qualities);
    device.bind("deltas", delta_view);
    device.bind("adjusted", adjusted);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    Kokkos::parallel_for("bqsr_quality_transform",
        Kokkos::RangePolicy<ExecSpace>(0, input.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            const int value = static_cast<int>(qualities(index)) + delta_view(index);
            // QualityUtils.boundQual(): recalibrated qualities are encoded in
            // [1,93]. A raw Q0/Q1 base can still be preserved by passing a
            // zero delta from Host before this transform.
            adjusted(index) = static_cast<std::uint8_t>(value < 1 ? 1 : value > 93 ? 93 : value);
        });
    ExecSpace().fence();
    plan.end_execute();
    result.execute_seconds = plan.telemetry().execute_seconds;
    auto host_adjusted = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), adjusted);
    for (std::size_t index = 0; index < input.size(); ++index) result.adjusted[index] = host_adjusted(index);
    return result;
}

}  // namespace fastgatk::kernels
