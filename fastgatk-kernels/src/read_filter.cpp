#include "fastgatk/kernels/read_filter.hpp"

#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace fastgatk::kernels {
namespace {
using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecSpace::memory_space;

// AssemblyRegionWalker delegates --max-reads-per-alignment-start to
// PositionalDownsampler, which uses java.util.Random's Algorithm R reservoir.
// Keep the stream-ordered reservoir on Host after Kokkos has produced the
// stable filter mask; this irregular selection must not be approximated by a
// per-record hash and does not belong in the PairHMM device kernel.
class JavaRandom {
public:
    explicit JavaRandom(const std::uint64_t seed)
        : seed_((seed ^ kMultiplier) & kMask) {}

    std::uint32_t next_int(const std::uint32_t bound) {
        if (bound == 0) throw std::invalid_argument("JavaRandom bound must be positive");
        if ((bound & (bound - 1U)) == 0U)
            return static_cast<std::uint32_t>((
                static_cast<std::uint64_t>(bound) * next_bits(31U)) >> 31U);
        while (true) {
            const auto bits = static_cast<std::uint32_t>(next_bits(31U));
            const auto value = bits % bound;
            // java.util.Random retries when signed 32-bit
            // `bits - value + (bound - 1)` overflows. Avoid signed overflow
            // while retaining that exact acceptance condition.
            if (bits - value <= std::numeric_limits<std::int32_t>::max() -
                                    static_cast<std::int64_t>(bound - 1U))
                return value;
        }
    }

private:
    static constexpr std::uint64_t kMultiplier = 0x5DEECE66DULL;
    static constexpr std::uint64_t kAddend = 0xBULL;
    static constexpr std::uint64_t kMask = (1ULL << 48U) - 1ULL;

    std::uint64_t next_bits(const unsigned int bits) {
        seed_ = (seed_ * kMultiplier + kAddend) & kMask;
        return seed_ >> (48U - bits);
    }

    std::uint64_t seed_;
};
}  // namespace

ReadFilterResult filter_reads_kokkos(const ReadFilterInput& input,
                                     ReadFilterOptions options) {
    if (!Kokkos::is_initialized()) throw std::runtime_error("Kokkos is not initialized");
    if (input.offsets.empty() || input.offsets.front() != 0 ||
        input.offsets.back() < input.offsets.front())
        throw std::invalid_argument("read filter offsets are malformed");
    for (std::size_t i = 1; i < input.offsets.size(); ++i)
        if (input.offsets[i] < input.offsets[i - 1])
            throw std::invalid_argument("read filter offsets are not monotonic");
    const auto reads = input.offsets.size() - 1;
    if (input.tids.size() != reads || input.positions.size() != reads || input.mapq.size() != reads)
        throw std::invalid_argument("read filter metadata does not match offsets");
    if (!input.flags.empty() && input.flags.size() != reads)
        throw std::invalid_argument("read filter flags do not match offsets");
    const bool has_cigar = !input.cigar_offsets.empty() || !input.cigar_ops.empty();
    if (has_cigar && (input.cigar_offsets.size() != reads + 1 || input.cigar_offsets.front() != 0 ||
                      input.cigar_offsets.back() != input.cigar_ops.size()))
        throw std::invalid_argument("read filter CIGAR metadata does not match offsets");
    if (has_cigar)
        for (std::size_t i = 1; i < input.cigar_offsets.size(); ++i)
            if (input.cigar_offsets[i] < input.cigar_offsets[i - 1])
                throw std::invalid_argument("read filter CIGAR offsets are not monotonic");
    const bool has_original_alignment = !input.original_alignment_offsets.empty() ||
        !input.original_alignments.empty();
    const bool has_mate_contig = !input.mate_contig_offsets.empty() ||
        !input.mate_contigs.empty();
    const bool has_original_alignment_presence = !input.original_alignment_present.empty();
    const bool has_mate_contig_presence = !input.mate_contig_present.empty();
    const bool has_read_group = !input.read_group_offsets.empty() ||
        !input.read_groups.empty();
    if (has_read_group &&
        (input.read_group_offsets.size() != reads + 1 ||
         input.read_group_offsets.front() != 0 ||
         input.read_group_offsets.back() != input.read_groups.size()))
        throw std::invalid_argument("read filter RG metadata does not match offsets");
    if (has_read_group)
        for (std::size_t i = 1; i < input.read_group_offsets.size(); ++i)
            if (input.read_group_offsets[i] < input.read_group_offsets[i - 1])
                throw std::invalid_argument("read filter RG offsets are not monotonic");
    if (has_original_alignment &&
        (input.original_alignment_offsets.size() != reads + 1 ||
         input.original_alignment_offsets.front() != 0 ||
         input.original_alignment_offsets.back() != input.original_alignments.size()))
        throw std::invalid_argument("read filter OA metadata does not match offsets");
    if (has_mate_contig &&
        (input.mate_contig_offsets.size() != reads + 1 ||
         input.mate_contig_offsets.front() != 0 ||
         input.mate_contig_offsets.back() != input.mate_contigs.size()))
        throw std::invalid_argument("read filter XM metadata does not match offsets");
    if (has_original_alignment_presence && input.original_alignment_present.size() != reads)
        throw std::invalid_argument("read filter OA presence metadata does not match offsets");
    if (has_mate_contig_presence && input.mate_contig_present.size() != reads)
        throw std::invalid_argument("read filter XM presence metadata does not match offsets");
    if (has_original_alignment)
        for (std::size_t i = 1; i < input.original_alignment_offsets.size(); ++i)
            if (input.original_alignment_offsets[i] < input.original_alignment_offsets[i - 1])
                throw std::invalid_argument("read filter OA offsets are not monotonic");
    if (has_mate_contig)
        for (std::size_t i = 1; i < input.mate_contig_offsets.size(); ++i)
            if (input.mate_contig_offsets[i] < input.mate_contig_offsets[i - 1])
                throw std::invalid_argument("read filter XM offsets are not monotonic");
    // An empty tile has no records to evaluate, so a selected CIGAR filter is
    // vacuously satisfied.  This matters for indexed region streaming, where
    // a bounded tile can legitimately contain no alignments; requiring a
    // metadata payload in that case would turn an otherwise valid empty
    // batch into a hard failure.
    if (reads != 0 && (options.require_good_cigar || options.require_nonzero_reference_span ||
                       options.require_no_n_cigar) && !has_cigar)
        throw std::invalid_argument("read filter CIGAR metadata is required by the selected filter");
    if (options.require_read_length && options.min_read_length > options.max_read_length)
        throw std::invalid_argument("read length filter minimum exceeds maximum");

    ReadFilterResult result;
    result.used = true;
    result.input_reads = reads;
    result.keep.assign(reads, 0);
    result.execution_space = ExecSpace::name();
    if (reads == 0) return result;

    fastgatk::core::HostBatch host("read-filter-v1");
    host.records = reads;
    host.bytes = input.offsets.size() * sizeof(std::uint32_t) +
                 input.tids.size() * sizeof(std::int32_t) * 2 + input.mapq.size();
    fastgatk::core::KernelPlan<ExecSpace> plan("read-filter");
    plan.begin_prepare(host);
    Kokkos::View<std::uint32_t*, MemorySpace> lengths("read_filter_lengths", reads);
    Kokkos::View<std::int32_t*, MemorySpace> tids("read_filter_tids", reads);
    Kokkos::View<std::int32_t*, MemorySpace> positions("read_filter_positions", reads);
    Kokkos::View<std::uint8_t*, MemorySpace> mapq("read_filter_mapq", reads);
    Kokkos::View<std::uint16_t*, MemorySpace> flags("read_filter_flags", reads);
    Kokkos::View<std::uint32_t*, MemorySpace> cigar_offsets("read_filter_cigar_offsets", reads + 1);
    Kokkos::View<std::uint32_t*, MemorySpace> cigar_ops("read_filter_cigar_ops", input.cigar_ops.size());
    Kokkos::View<std::uint32_t*, MemorySpace> oa_offsets("read_filter_oa_offsets", reads + 1);
    Kokkos::View<std::uint8_t*, MemorySpace> oa_payload("read_filter_oa_payload", input.original_alignments.size());
    Kokkos::View<std::uint8_t*, MemorySpace> oa_present("read_filter_oa_present", reads);
    Kokkos::View<std::uint32_t*, MemorySpace> xm_offsets("read_filter_xm_offsets", reads + 1);
    Kokkos::View<std::uint8_t*, MemorySpace> xm_payload("read_filter_xm_payload", input.mate_contigs.size());
    Kokkos::View<std::uint8_t*, MemorySpace> xm_present("read_filter_xm_present", reads);
    Kokkos::View<std::uint32_t*, MemorySpace> rg_offsets("read_filter_rg_offsets", reads + 1);
    Kokkos::View<std::uint8_t*, MemorySpace> rg_payload("read_filter_rg_payload", input.read_groups.size());
    Kokkos::View<std::uint8_t*, MemorySpace> keep("read_filter_keep", reads);
    auto host_lengths = Kokkos::create_mirror_view(lengths);
    auto host_tids = Kokkos::create_mirror_view(tids);
    auto host_positions = Kokkos::create_mirror_view(positions);
    auto host_mapq = Kokkos::create_mirror_view(mapq);
    auto host_flags = Kokkos::create_mirror_view(flags);
    auto host_cigar_offsets = Kokkos::create_mirror_view(cigar_offsets);
    auto host_cigar_ops = Kokkos::create_mirror_view(cigar_ops);
    auto host_oa_offsets = Kokkos::create_mirror_view(oa_offsets);
    auto host_oa_payload = Kokkos::create_mirror_view(oa_payload);
    auto host_oa_present = Kokkos::create_mirror_view(oa_present);
    auto host_xm_offsets = Kokkos::create_mirror_view(xm_offsets);
    auto host_xm_payload = Kokkos::create_mirror_view(xm_payload);
    auto host_xm_present = Kokkos::create_mirror_view(xm_present);
    auto host_rg_offsets = Kokkos::create_mirror_view(rg_offsets);
    auto host_rg_payload = Kokkos::create_mirror_view(rg_payload);
    for (std::size_t i = 0; i < reads; ++i) {
        host_lengths(i) = input.offsets[i + 1] - input.offsets[i];
        host_tids(i) = input.tids[i];
        host_positions(i) = input.positions[i];
        host_mapq(i) = input.mapq[i];
        host_flags(i) = input.flags.empty() ? 0 : input.flags[i];
    }
    for (std::size_t i = 0; i < input.cigar_offsets.size(); ++i) host_cigar_offsets(i) = input.cigar_offsets[i];
    for (std::size_t i = 0; i < input.cigar_ops.size(); ++i) host_cigar_ops(i) = input.cigar_ops[i];
    for (std::size_t i = 0; i <= reads; ++i) host_oa_offsets(i) = 0;
    for (std::size_t i = 0; i < input.original_alignment_offsets.size(); ++i)
        host_oa_offsets(i) = input.original_alignment_offsets[i];
    for (std::size_t i = 0; i < input.original_alignments.size(); ++i)
        host_oa_payload(i) = input.original_alignments[i];
    for (std::size_t i = 0; i < reads; ++i) {
        const auto begin = input.original_alignment_offsets.size() == reads + 1
            ? input.original_alignment_offsets[i] : 0U;
        const auto end = input.original_alignment_offsets.size() == reads + 1
            ? input.original_alignment_offsets[i + 1] : begin;
        host_oa_present(i) = has_original_alignment_presence
            ? input.original_alignment_present[i]
            : (end > begin ? 1U : 0U);
    }
    for (std::size_t i = 0; i <= reads; ++i) host_xm_offsets(i) = 0;
    for (std::size_t i = 0; i < input.mate_contig_offsets.size(); ++i)
        host_xm_offsets(i) = input.mate_contig_offsets[i];
    for (std::size_t i = 0; i < input.mate_contigs.size(); ++i)
        host_xm_payload(i) = input.mate_contigs[i];
    for (std::size_t i = 0; i < reads; ++i) {
        const auto begin = input.mate_contig_offsets.size() == reads + 1
            ? input.mate_contig_offsets[i] : 0U;
        const auto end = input.mate_contig_offsets.size() == reads + 1
            ? input.mate_contig_offsets[i + 1] : begin;
        host_xm_present(i) = has_mate_contig_presence
            ? input.mate_contig_present[i]
            : (end > begin ? 1U : 0U);
    }
    for (std::size_t i = 0; i <= reads; ++i) host_rg_offsets(i) = 0;
    for (std::size_t i = 0; i < input.read_group_offsets.size(); ++i)
        host_rg_offsets(i) = input.read_group_offsets[i];
    for (std::size_t i = 0; i < input.read_groups.size(); ++i)
        host_rg_payload(i) = input.read_groups[i];
    Kokkos::deep_copy(lengths, host_lengths);
    Kokkos::deep_copy(tids, host_tids);
    Kokkos::deep_copy(positions, host_positions);
    Kokkos::deep_copy(mapq, host_mapq);
    Kokkos::deep_copy(flags, host_flags);
    Kokkos::deep_copy(cigar_offsets, host_cigar_offsets);
    Kokkos::deep_copy(cigar_ops, host_cigar_ops);
    Kokkos::deep_copy(oa_offsets, host_oa_offsets);
    Kokkos::deep_copy(oa_payload, host_oa_payload);
    Kokkos::deep_copy(oa_present, host_oa_present);
    Kokkos::deep_copy(xm_offsets, host_xm_offsets);
    Kokkos::deep_copy(xm_payload, host_xm_payload);
    Kokkos::deep_copy(xm_present, host_xm_present);
    Kokkos::deep_copy(rg_offsets, host_rg_offsets);
    Kokkos::deep_copy(rg_payload, host_rg_payload);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device(reads);
    device.bind("lengths", lengths); device.bind("tids", tids);
    device.bind("positions", positions); device.bind("mapq", mapq);
    device.bind("flags", flags); device.bind("keep", keep);
    device.bind("cigar_offsets", cigar_offsets); device.bind("cigar_ops", cigar_ops);
    device.bind("oa_offsets", oa_offsets); device.bind("oa_payload", oa_payload);
    device.bind("oa_present", oa_present);
    device.bind("xm_offsets", xm_offsets); device.bind("xm_payload", xm_payload);
    device.bind("xm_present", xm_present);
    device.bind("rg_offsets", rg_offsets); device.bind("rg_payload", rg_payload);
    plan.end_prepare(device);
    result.prepare_seconds = plan.telemetry().prepare_seconds;
    plan.begin_execute();
    const auto min_mapq = options.min_mapq;
    const auto exclude_mapping_quality_unavailable = options.exclude_mapping_quality_unavailable;
    const auto exclude_mapping_quality_zero = options.exclude_mapping_quality_zero;
    const auto exclude_unmapped = options.exclude_unmapped;
    const auto exclude_secondary = options.exclude_secondary;
    const auto exclude_supplementary = options.exclude_supplementary;
    const auto exclude_duplicates = options.exclude_duplicates;
    const auto exclude_qcfail = options.exclude_qcfail;
    const auto require_good_cigar = options.require_good_cigar;
    const auto require_nonzero_reference_span = options.require_nonzero_reference_span;
    const auto require_non_chimeric_original_alignment =
        options.require_non_chimeric_original_alignment;
    const auto require_no_n_cigar = options.require_no_n_cigar;
    const auto require_read_group = options.require_read_group;
    const auto require_read_length = options.require_read_length;
    const auto min_read_length = options.min_read_length;
    const auto max_read_length = options.max_read_length;
    Kokkos::parallel_for("read_filter", Kokkos::RangePolicy<ExecSpace>(0, reads),
        KOKKOS_LAMBDA(const std::size_t i) {
            const auto flag = flags(i);
            bool accepted = lengths(i) > 0 && mapq(i) >= min_mapq && positions(i) >= 0;
            accepted = accepted && (!exclude_mapping_quality_unavailable || mapq(i) != 255U);
            accepted = accepted && (!exclude_mapping_quality_zero || mapq(i) != 0U);
            accepted = accepted && (!require_read_length ||
                                    (lengths(i) >= min_read_length && lengths(i) <= max_read_length));
            accepted = accepted && (!exclude_unmapped || (flag & 0x4U) == 0);
            accepted = accepted && (!exclude_secondary || (flag & 0x100U) == 0);
            accepted = accepted && (!exclude_supplementary || (flag & 0x800U) == 0);
            accepted = accepted && (!exclude_duplicates || (flag & 0x400U) == 0);
            accepted = accepted && (!exclude_qcfail || (flag & 0x200U) == 0);
            accepted = accepted && (!require_read_group ||
                                    (has_read_group && rg_offsets(i + 1) > rg_offsets(i)));
            if (require_non_chimeric_original_alignment) {
                const auto oa_begin = oa_offsets(i);
                const auto oa_end = oa_offsets(i + 1);
                const auto xm_begin = xm_offsets(i);
                const auto xm_end = xm_offsets(i + 1);
                // GATK accepts reads missing either tag. When both are
                // present, OA's contig prefix must equal XM exactly.
                if (oa_present(i) != 0U && xm_present(i) != 0U) {
                    std::uint32_t oa_contig_end = oa_begin;
                    while (oa_contig_end < oa_end && oa_payload(oa_contig_end) != ',')
                        ++oa_contig_end;
                    const auto oa_length = oa_contig_end - oa_begin;
                    const auto xm_length = xm_end - xm_begin;
                    bool same = oa_length == xm_length;
                    for (std::uint32_t cursor = 0; same && cursor < oa_length; ++cursor)
                        same = oa_payload(oa_begin + cursor) == xm_payload(xm_begin + cursor);
                    accepted = accepted && same;
                }
            }
            if (require_good_cigar || require_nonzero_reference_span || require_no_n_cigar) {
                bool good = true;
                bool consumes_reference = false;
                std::uint64_t read_consumed = 0;
                const auto cigar_begin = cigar_offsets(i);
                const auto cigar_end = cigar_offsets(i + 1);
                if (cigar_end <= cigar_begin) good = false;
                for (std::uint32_t cursor = cigar_begin; good && cursor < cigar_end; ++cursor) {
                    const auto packed = cigar_ops(cursor);
                    const auto length = packed >> 4;
                    const auto code = packed & 0x0fU;
                    if (length == 0 || code > 8U) { good = false; break; }
                    if (require_no_n_cigar && code == 3U) { good = false; break; }
                    const bool consumes_read = code == 0U || code == 1U || code == 4U || code == 7U || code == 8U;
                    const bool consumes_ref = code == 0U || code == 2U || code == 3U || code == 7U || code == 8U;
                    if (consumes_read) read_consumed += length;
                    if (consumes_ref) consumes_reference = true;
                }
                const bool unmapped = (flag & 0x4U) != 0;
                if (good && !(unmapped && read_consumed == 0) && read_consumed != lengths(i)) good = false;
                accepted = accepted && good && (!require_nonzero_reference_span || consumes_reference);
            }
            keep(i) = accepted ? 1 : 0;
        });
    ExecSpace().fence();
    plan.end_execute();
    result.seconds = plan.telemetry().execute_seconds;
    auto host_keep = Kokkos::create_mirror_view(keep);
    Kokkos::deep_copy(host_keep, keep);
    for (std::size_t i = 0; i < reads; ++i) result.keep[i] = host_keep(i);

    for (const auto value : result.keep) if (value) ++result.passed_reads;
    result.filtered_reads = result.input_reads - result.passed_reads;

    struct Key { std::int32_t tid, position; bool operator<(const Key& other) const {
        return tid < other.tid || (tid == other.tid && position < other.position);
    }};
    std::map<Key, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < reads; ++i)
        if (result.keep[i]) groups[Key{input.tids[i], input.positions[i]}].push_back(i);

    if (options.mutect2_downsampling) {
        // Mutect2 creates MutectDownsampler instead of the ordinary
        // AssemblyRegionWalker PositionalDownsampler.  Its stride begins at
        // the first *observed* alignment start (rather than a genomic window
        // boundary), and it maintains the shared GATK Java RNG across all
        // finalized strides.
        const auto stride = static_cast<std::size_t>(options.mutect2_downsampling_stride);
        if (stride == 0)
            throw std::invalid_argument("MutectDownsampler stride must be positive");
        const auto suspicious_per_start =
            options.mutect2_max_suspicious_reads_per_alignment_start;
        if (options.max_reads_per_locus == 0 && suspicious_per_start <= 0)
            return result;
        constexpr auto java_int_max = static_cast<std::size_t>(
            std::numeric_limits<std::int32_t>::max());
        std::size_t max_coverage = java_int_max;
        if (options.max_reads_per_locus != 0) {
            if (options.max_reads_per_locus > java_int_max ||
                stride > java_int_max / options.max_reads_per_locus)
                throw std::invalid_argument("MutectDownsampler coverage exceeds GATK int range");
            max_coverage = options.max_reads_per_locus * stride;
        }
        std::size_t suspicious_per_stride = java_int_max;
        if (suspicious_per_start > 0) {
            const auto positive_limit = static_cast<std::size_t>(suspicious_per_start);
            if (positive_limit > java_int_max || stride > java_int_max / positive_limit)
                throw std::invalid_argument("MutectDownsampler suspicious-read limit exceeds GATK int range");
            suspicious_per_stride = positive_limit * stride;
        }

        struct StrideGroup {
            std::int32_t tid = -1;
            std::int32_t first_position = -1;
            std::vector<std::size_t> members;
        };
        std::vector<StrideGroup> stride_groups;
        for (const auto& [key, members] : groups) {
            const bool starts_new_stride = stride_groups.empty() ||
                key.tid != stride_groups.back().tid ||
                static_cast<std::size_t>(key.position - stride_groups.back().first_position) >= stride;
            if (starts_new_stride)
                stride_groups.push_back(StrideGroup{key.tid, key.position, {}});
            auto& target = stride_groups.back().members;
            target.insert(target.end(), members.begin(), members.end());
        }

        JavaRandom random(options.seed);
        for (const auto& group : stride_groups) {
            const auto& members = group.members;
            std::size_t suspicious_reads = 0;
            bool reject_all = false;
            for (const auto member : members) {
                if (input.mapq[member] <= 50U) ++suspicious_reads;
                if (suspicious_reads >= suspicious_per_stride) {
                    reject_all = true;
                    break;
                }
            }
            if (reject_all) {
                for (const auto member : members) result.keep[member] = 0;
                result.downsampled_reads += members.size();
                continue;
            }
            if (members.size() <= max_coverage) continue;
            if (members.size() > java_int_max)
                throw std::invalid_argument("MutectDownsampler reservoir exceeds GATK int range");
            std::vector<std::size_t> reservoir;
            reservoir.reserve(max_coverage);
            std::size_t submitted = 0;
            for (const auto member : members) {
                // MutectDownsampler's overloaded-pool path discards
                // suspicious reads and feeds only MAPQ > 50 records to its
                // Java ReservoirDownsampler.
                if (input.mapq[member] <= 50U) continue;
                if (submitted < max_coverage) {
                    reservoir.push_back(member);
                } else {
                    const auto slot = random.next_int(
                        static_cast<std::uint32_t>(submitted + 1U));
                    if (slot < max_coverage) reservoir[slot] = member;
                }
                ++submitted;
            }
            for (const auto member : members) result.keep[member] = 0;
            for (const auto member : reservoir) result.keep[member] = 1;
            result.downsampled_reads += members.size() - reservoir.size();
        }
        result.passed_reads -= result.downsampled_reads;
        return result;
    }

    if (options.max_reads_per_locus == 0) return result;
    JavaRandom random(options.seed);
    for (const auto& [key, members] : groups) {
        static_cast<void>(key);
        if (members.size() <= options.max_reads_per_locus) continue;
        if (options.max_reads_per_locus >
                static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) ||
            members.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            throw std::invalid_argument("read downsampling reservoir exceeds GATK int range");
        std::vector<std::size_t> reservoir;
        reservoir.reserve(options.max_reads_per_locus);
        std::size_t submitted = 0;
        for (const auto member : members) {
            if (submitted < options.max_reads_per_locus) {
                reservoir.push_back(member);
                ++submitted;
                continue;
            }
            const auto slot = random.next_int(static_cast<std::uint32_t>(submitted + 1U));
            if (slot < options.max_reads_per_locus)
                reservoir[slot] = member;
            ++submitted;
        }
        for (const auto member : members) result.keep[member] = 0;
        for (const auto member : reservoir) result.keep[member] = 1;
        result.downsampled_reads += members.size() - reservoir.size();
    }
    result.passed_reads -= result.downsampled_reads;
    return result;
}

}  // namespace fastgatk::kernels
