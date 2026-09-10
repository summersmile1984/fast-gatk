#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/faidx.h>
#endif

namespace {

constexpr std::array<char, 4> BASES{'A', 'C', 'G', 'T'};

struct Options {
    std::string input;
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    std::string output;
    std::string manifest;
    std::string sample;
    std::size_t batch_records = 4096;
    int minimum_mapping_quality = 30;
    int minimum_base_quality = 20;
    int threads = 1;
    bool include_duplicates = false;
    // CollectAllelicCounts inherits LocusWalker defaults (Wellformed and
    // Mapped) and adds Mapped, NonZeroReferenceLength, NotDuplicate and
    // MappingQuality(>=30).  Keep this state explicit so the common GATK
    // read-filter controls are real controls rather than ignored flags.
    bool disable_tool_default_read_filters = false;
    bool wellformed_filter = true;
    bool mapped_filter = true;
    bool nonzero_reference_filter = true;
    bool duplicate_filter = true;
    bool mapping_quality_filter = true;
    bool vendor_quality_filter = false;
    bool secondary_filter = false;
    bool supplementary_filter = false;
    std::vector<std::string> disabled_read_filters;
    std::vector<std::string> enabled_read_filters;
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    fastgatk::io::HtsIntervalSetRule interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
};

struct Locus {
    std::int32_t tid = -1;
    std::int64_t position = 0;
    char reference = 'N';
};

// One result per decoded ReadBatch keeps all mutable read counters on the
// compute stage. The Kokkos count view is shared by that single stage, while
// the result itself crosses the normal encode/sink boundary deterministically.
struct CountBatchResult {
    std::uint64_t reads_seen = 0;
    std::uint64_t reads_used = 0;
    std::uint64_t reads_filtered = 0;
    std::uint64_t projected_bases = 0;
    std::uint64_t kernel_batches = 0;
    std::uint64_t kernel_observations = 0;
    double kernel_prepare_seconds = 0.0;
    double kernel_execute_seconds = 0.0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

enum class ReadFilterKind : std::uint8_t {
    Wellformed,
    Mapped,
    NonZeroReferenceLength,
    NotDuplicate,
    MappingQuality,
    PassesVendorQuality,
    NotSecondary,
    NotSupplementary,
};

ReadFilterKind read_filter_kind(const std::string& name) {
    if (name == "WellformedReadFilter") return ReadFilterKind::Wellformed;
    if (name == "MappedReadFilter") return ReadFilterKind::Mapped;
    if (name == "NonZeroReferenceLengthAlignmentReadFilter")
        return ReadFilterKind::NonZeroReferenceLength;
    if (name == "NotDuplicateReadFilter") return ReadFilterKind::NotDuplicate;
    if (name == "MappingQualityReadFilter") return ReadFilterKind::MappingQuality;
    if (name == "PassesVendorQualityCheckReadFilter")
        return ReadFilterKind::PassesVendorQuality;
    if (name == "NotSecondaryAlignmentReadFilter") return ReadFilterKind::NotSecondary;
    if (name == "NotSupplementaryAlignmentReadFilter") return ReadFilterKind::NotSupplementary;
    throw std::invalid_argument("UNSUPPORTED_PARAMETER: read filter " + name);
}

void set_read_filter(Options& options, ReadFilterKind kind, bool enabled) {
    switch (kind) {
        case ReadFilterKind::Wellformed: options.wellformed_filter = enabled; break;
        case ReadFilterKind::Mapped: options.mapped_filter = enabled; break;
        case ReadFilterKind::NonZeroReferenceLength: options.nonzero_reference_filter = enabled; break;
        case ReadFilterKind::NotDuplicate:
            options.duplicate_filter = enabled;
            options.include_duplicates = !enabled;
            break;
        case ReadFilterKind::MappingQuality: options.mapping_quality_filter = enabled; break;
        case ReadFilterKind::PassesVendorQuality: options.vendor_quality_filter = enabled; break;
        case ReadFilterKind::NotSecondary: options.secondary_filter = enabled; break;
        case ReadFilterKind::NotSupplementary: options.supplementary_filter = enabled; break;
    }
}

bool is_tool_default_filter(ReadFilterKind kind) {
    switch (kind) {
        case ReadFilterKind::Wellformed:
        case ReadFilterKind::Mapped:
        case ReadFilterKind::NonZeroReferenceLength:
        case ReadFilterKind::NotDuplicate:
        case ReadFilterKind::MappingQuality:
            return true;
        case ReadFilterKind::PassesVendorQuality:
        case ReadFilterKind::NotSecondary:
        case ReadFilterKind::NotSupplementary:
            return false;
    }
    return false;
}

void disable_tool_defaults(Options& options) {
    // This is deliberately the complete GATK default list for this walker:
    // Wellformed + Mapped from LocusWalker, then the four tool additions.
    options.wellformed_filter = false;
    options.mapped_filter = false;
    options.nonzero_reference_filter = false;
    options.duplicate_filter = false;
    options.mapping_quality_filter = false;
    options.include_duplicates = true;
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-collect-allelic-counts (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                    BAM/CRAM/SAM input\n"
                         "  -R, --reference FILE                reference FASTA\n"
                         "  -L, --intervals REGION               interval/file selector (repeatable)\n"
                         "  -XL, --exclude-intervals REGION     subtract interval/file selector (repeatable)\n"
                         "      --interval-padding N            pad included intervals\n"
                         "      --interval-exclusion-padding N  pad exclusions before subtraction\n"
                         "      --interval-set-rule RULE        UNION (default) or INTERSECTION\n"
                         "  -O, --output FILE                   allelic-count TSV\n"
                         "      --minimum-mapping-quality N     default 30\n"
                         "      --minimum-base-quality N        default 20\n"
                         "      --include-duplicates            retain duplicate reads\n"
                         "  -DF, --disable-read-filter NAME    disable a default read filter\n"
                         "      --read-filter NAME             enable a read filter\n"
                         "      --disable-tool-default-read-filters [BOOL]\n"
                         "                                      disable all tool default filters\n"
                         "      --sample NAME                   sample metadata\n"
                         "      --output-manifest FILE          OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                has_option(argument, "--region") ? "--region" :
                has_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals")) {
            options.exclusions.push_back(require_value(index, argc, argv, argument,
                                                        "--exclude-intervals", "-XL"));
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--minimum-mapping-quality")) {
            options.minimum_mapping_quality = std::stoi(require_value(
                index, argc, argv, argument, "--minimum-mapping-quality"));
        } else if (has_option(argument, "--minimum-base-quality")) {
            options.minimum_base_quality = std::stoi(require_value(
                index, argc, argv, argument, "--minimum-base-quality"));
        } else if (has_option(argument, "--batch-records")) {
            options.batch_records = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--batch-records")));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--sample")) {
            options.sample = require_value(index, argc, argv, argument, "--sample");
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "-imr" || has_option(argument, "--interval-merging-rule")) {
            const auto rule = require_value(index, argc, argv, argument, "--interval-merging-rule", "-imr");
            if (rule != "OVERLAPPING_ONLY")
                throw std::invalid_argument("native CollectAllelicCounts requires --interval-merging-rule OVERLAPPING_ONLY");
        } else if (argument == "-ip" || has_option(argument, "--interval-padding") ||
                   argument == "-ixp" || has_option(argument, "--interval-exclusion-padding")) {
            const bool exclusion = argument == "-ixp" ||
                argument.rfind("--interval-exclusion-padding", 0) == 0;
            const char* name = exclusion ? "--interval-exclusion-padding" : "--interval-padding";
            const auto value = std::stoll(require_value(index, argc, argv, argument, name,
                                                         exclusion ? "-ixp" : "-ip"));
            if (value < 0) throw std::invalid_argument(std::string(name) + " must be non-negative");
            if (exclusion) options.exclusion_padding = value;
            else options.interval_padding = value;
        } else if (argument == "-isr" || has_option(argument, "--interval-set-rule")) {
            auto rule = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char c) {
                return static_cast<char>(std::toupper(c));
            });
            if (rule == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (rule == "INTERSECTION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + rule);
        } else if (argument == "--include-duplicates") {
            options.include_duplicates = true;
            options.duplicate_filter = false;
        } else if (argument == "--disable-tool-default-read-filters" ||
                   has_option(argument, "--disable-tool-default-read-filters")) {
            options.disable_tool_default_read_filters =
                fastgatk::native::parse_optional_boolean(
                    index, argc, argv, argument, "--disable-tool-default-read-filters");
            if (options.disable_tool_default_read_filters) disable_tool_defaults(options);
        } else if (argument == "-DF" || has_option(argument, "--disable-read-filter")) {
            const auto filter = require_value(index, argc, argv, argument,
                                               "--disable-read-filter", "-DF");
            // GATK only permits disabling filters in the tool's resolved
            // default list.  Rejecting a non-default name is safer than
            // silently running with a different pileup.
            const auto kind = read_filter_kind(filter);
            if (!is_tool_default_filter(kind))
                throw std::invalid_argument("UNSUPPORTED_PARAMETER: --disable-read-filter " + filter);
            options.disabled_read_filters.push_back(filter);
        } else if (argument == "-RF" || has_option(argument, "--read-filter")) {
            const auto filter = require_value(index, argc, argv, argument,
                                               "--read-filter", "-RF");
            (void)read_filter_kind(filter);
            options.enabled_read_filters.push_back(filter);
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility switches; the native filter policy is recorded in the manifest.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    // Barclay resolves filter directives after parsing, independently of
    // argument order: disabling tool defaults removes the complete default
    // list, then explicit --read-filter directives are added back, and
    // --disable-read-filter removes defaults by class name.  Reproduce that
    // ordering instead of making command-line ordering observable.
    if (options.disable_tool_default_read_filters) disable_tool_defaults(options);
    for (const auto& filter : options.disabled_read_filters)
        set_read_filter(options, read_filter_kind(filter), false);
    for (const auto& filter : options.enabled_read_filters)
        set_read_filter(options, read_filter_kind(filter), true);
    if (options.include_duplicates) options.duplicate_filter = false;
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.regions.empty()) throw std::invalid_argument("-L/--intervals is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.minimum_mapping_quality < 0 || options.minimum_mapping_quality > 255 ||
        options.minimum_base_quality < 0 || options.minimum_base_quality > 93 ||
        options.batch_records == 0 || options.threads < 1)
        throw std::invalid_argument("invalid allelic-count resource/filter parameters");
    return options;
}

int base_index(char base) {
    base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    for (std::size_t i = 0; i < BASES.size(); ++i)
        if (BASES[i] == base) return static_cast<int>(i);
    return -1;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto c : text) {
        if (c == '"' || c == '\\') escaped << '\\';
        if (c == '\n') escaped << "\\n";
        else if (c == '\r') escaped << "\\r";
        else if (c == '\t') escaped << "\\t";
        else escaped << c;
    }
    return escaped.str();
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

#if FASTGATK_HAS_HTSLIB
std::vector<char> fetch_reference(faidx_t* reference, const std::string& contig,
                                  std::int64_t start, std::int64_t end) {
    if (end <= start) return {};
    int fetched_length = 0;
    char* sequence = faidx_fetch_seq(reference, contig.c_str(), static_cast<int>(start),
                                     static_cast<int>(end - 1), &fetched_length);
    if (sequence == nullptr || fetched_length != end - start) {
        free(sequence);
        throw std::runtime_error("BAD_INPUT: reference FASTA does not cover interval " +
                                 contig + ":" + std::to_string(start + 1) + "-" + std::to_string(end));
    }
    std::vector<char> result(sequence, sequence + fetched_length);
    free(sequence);
    for (auto& base : result) base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    return result;
}
#endif

std::size_t find_locus(const std::vector<fastgatk::io::HtsInterval>& intervals,
                       const std::vector<std::size_t>& offsets,
                       std::int32_t tid, std::int64_t position) {
    std::size_t low = 0;
    std::size_t high = intervals.size();
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        const auto& interval = intervals[middle];
        if (interval.tid < tid || (interval.tid == tid && interval.end <= position)) low = middle + 1;
        else high = middle;
    }
    if (low >= intervals.size() || intervals[low].tid != tid ||
        position < intervals[low].start || position >= intervals[low].end)
        return std::numeric_limits<std::size_t>::max();
    return offsets[low] + static_cast<std::size_t>(position - intervals[low].start);
}

std::vector<fastgatk::io::HtsInterval> subtract_intervals(
    const std::vector<fastgatk::io::HtsInterval>& input,
    const std::vector<fastgatk::io::HtsInterval>& exclusions) {
    if (exclusions.empty()) return input;
    std::vector<fastgatk::io::HtsInterval> result;
    for (const auto& source : input) {
        std::int64_t cursor = source.start;
        for (const auto& exclusion : exclusions) {
            if (exclusion.tid != source.tid) continue;
            if (exclusion.end <= cursor) continue;
            if (exclusion.start >= source.end) break;
            if (exclusion.start > cursor)
                result.push_back({source.tid, cursor, std::min(source.end, exclusion.start)});
            if (exclusion.end >= source.end) { cursor = source.end; break; }
            cursor = std::max(cursor, exclusion.end);
        }
        if (cursor < source.end) result.push_back({source.tid, cursor, source.end});
    }
    return result;
}

bool wellformed_record(const fastgatk::io::ReadBatch& batch, std::size_t record) {
    if (record >= batch.positions.size() || record >= batch.tids.size() ||
        batch.positions[record] < 0 || batch.tids[record] < 0 || !batch.has_cigar() ||
        !batch.cigar_record_layout_valid(record))
        return false;
    // GATK's WellformedReadFilter includes CigarContainsNoNOperator.  A
    // skipped reference block is therefore not part of the default
    // CollectAllelicCounts pileup, even though HTSlib can project bases past
    // the N operation safely.
    for (std::size_t i = batch.cigar_offsets[record]; i < batch.cigar_offsets[record + 1]; ++i) {
        if (fastgatk::io::CigarOp::unpack(batch.cigar_ops[i]).code ==
            fastgatk::io::CigarOpCode::ReferenceSkip)
            return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
#if !FASTGATK_HAS_HTSLIB
    (void)argc;
    (void)argv;
    std::cerr << "BACKEND_UNAVAILABLE: build with HTSlib to run CollectAllelicCounts\n";
    return 2;
#else
    bool initialized = false;
    faidx_t* reference = nullptr;
    try {
        const auto started = std::chrono::steady_clock::now();
        auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;

        fastgatk::io::HtsReader reader(options.input, options.reference, options.regions,
                                       options.batch_records, {}, options.exclusions,
                                       options.interval_padding, options.exclusion_padding,
                                       options.interval_set_rule);
        if (!reader.compiled()) throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib is not available");
        const auto& header_samples = reader.header().samples;
        std::string sample_source = "cli";
        if (options.sample.empty()) {
            if (header_samples.size() != 1)
                throw std::runtime_error(
                    "BAD_INPUT: CollectAllelicCounts requires --sample when BAM header has zero or multiple SM samples");
            options.sample = header_samples.front();
            sample_source = "read-group-header";
        }
        reference = fai_load(options.reference.c_str());
        if (!reference) throw std::runtime_error("BAD_INPUT: cannot load reference FASTA/index: " + options.reference);
        const auto excluded_interval_count = reader.exclusion_intervals().size();
        const auto intervals = subtract_intervals(reader.intervals(), reader.exclusion_intervals());
        if (intervals.empty()) throw std::runtime_error("BAD_INPUT: no intervals were selected");

        std::vector<std::size_t> offsets;
        offsets.reserve(intervals.size() + 1);
        offsets.push_back(0);
        std::vector<Locus> loci;
        for (const auto& interval : intervals) {
            if (interval.end - interval.start > 50'000'000 ||
                loci.size() > 50'000'000 - static_cast<std::size_t>(interval.end - interval.start))
                throw std::runtime_error("RESOURCE_EXHAUSTED: selected loci exceed native memory safety limit");
            const auto bases = fetch_reference(reference,
                reader.header().contigs.at(static_cast<std::size_t>(interval.tid)), interval.start, interval.end);
            for (std::size_t i = 0; i < bases.size(); ++i)
                loci.push_back(Locus{interval.tid, interval.start + static_cast<std::int64_t>(i), bases[i]});
            offsets.push_back(loci.size());
        }
        std::uint64_t reads_seen = 0;
        std::uint64_t reads_used = 0;
        std::uint64_t reads_filtered = 0;
        std::uint64_t projected_bases = 0;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        std::uint64_t count_kernel_batches = 0;
        std::uint64_t count_kernel_observations = 0;
        double count_kernel_prepare_seconds = 0.0;
        double count_kernel_execute_seconds = 0.0;
        fastgatk::runtime::ThreeStagePipeline<fastgatk::io::ReadBatch, CountBatchResult, CountBatchResult>::Metrics pipeline_metrics;
        std::vector<std::uint64_t> final_counts;
        {
            Kokkos::View<std::uint64_t*> device_counts("collect_allelic_counts", loci.size() * BASES.size());
            Kokkos::deep_copy(device_counts, std::uint64_t{0});
            using Pipeline = fastgatk::runtime::ThreeStagePipeline<fastgatk::io::ReadBatch, CountBatchResult, CountBatchResult>;
            const auto safe_budget = resources.safe_memory_budget_bytes();
            const auto decoded_capacity = safe_budget > 0
                ? std::max<std::uint64_t>(safe_budget / 4, 1)
                : std::max<std::uint64_t>(64ULL * 1024ULL * 1024ULL,
                    static_cast<std::uint64_t>(options.batch_records) * 256ULL);
            const auto result_capacity = std::max<std::uint64_t>(sizeof(CountBatchResult) * 8ULL, 4096ULL);
            Pipeline pipeline(
                Pipeline::Limits{decoded_capacity, result_capacity, result_capacity},
                [&reader]() -> std::optional<fastgatk::io::ReadBatch> {
                    fastgatk::io::ReadBatch batch;
                    if (!reader.next(batch)) return std::nullopt;
                    return batch;
                },
                [&](fastgatk::io::ReadBatch batch) -> std::optional<CountBatchResult> {
                    CountBatchResult result;
                    result.reads_seen = batch.records();
                    std::vector<std::size_t> observation_loci;
                    std::vector<std::uint8_t> observation_bases;
                    for (std::size_t record = 0; record < batch.records(); ++record) {
                        const auto flags = record < batch.flags.size() ? batch.flags[record] : 0xffffU;
                        const auto mapq = record < batch.mapq.size() ? batch.mapq[record] : 0;
                        const auto tid = record < batch.tids.size() ? batch.tids[record] : -1;
                        const auto position = record < batch.positions.size() ? batch.positions[record] : -1;
                        const auto reference_end = fastgatk::io::reference_end(batch, record);
                        const bool filtered =
                            (options.wellformed_filter && !wellformed_record(batch, record)) ||
                            (options.mapped_filter && ((flags & 0x4U) != 0 || tid < 0 || position < 0)) ||
                            (options.nonzero_reference_filter && reference_end <= position) ||
                            (options.duplicate_filter && (flags & 0x400U) != 0) ||
                            (options.mapping_quality_filter && mapq < options.minimum_mapping_quality) ||
                            (options.vendor_quality_filter && (flags & 0x200U) != 0) ||
                            (options.secondary_filter && (flags & 0x100U) != 0) ||
                            (options.supplementary_filter && (flags & 0x800U) != 0);
                        if (filtered) {
                            ++result.reads_filtered;
                            continue;
                        }
                        ++result.reads_used;
                        const auto begin = batch.offsets[record];
                        const auto end = batch.offsets[record + 1];
                        for (std::size_t read_offset = 0; read_offset < end - begin; ++read_offset) {
                            if (batch.qualities[begin + read_offset] < options.minimum_base_quality) continue;
                            const int base = base_index(static_cast<char>(batch.bases[begin + read_offset]));
                            if (base < 0) continue;
                            fastgatk::io::ReadProjection projection;
                            if (!fastgatk::io::project_read_offset(batch, record, read_offset, projection) ||
                                projection.operation == fastgatk::io::CigarOpCode::Insertion ||
                                projection.operation == fastgatk::io::CigarOpCode::SoftClip)
                                continue;
                            const auto locus = find_locus(intervals, offsets, batch.tids[record], projection.reference_position);
                            if (locus == std::numeric_limits<std::size_t>::max()) continue;
                            observation_loci.push_back(locus);
                            observation_bases.push_back(static_cast<std::uint8_t>(base));
                            ++result.projected_bases;
                        }
                    }
                    if (observation_loci.empty()) return result;
                    Kokkos::View<std::size_t*> loci_view("allelic_observation_loci", observation_loci.size());
                    Kokkos::View<std::uint8_t*> bases_view("allelic_observation_bases", observation_bases.size());
                    auto host_loci = Kokkos::create_mirror_view(loci_view);
                    auto host_bases = Kokkos::create_mirror_view(bases_view);
                    for (std::size_t i = 0; i < observation_loci.size(); ++i) {
                        host_loci(i) = observation_loci[i];
                        host_bases(i) = observation_bases[i];
                    }
                    Kokkos::deep_copy(loci_view, host_loci);
                    Kokkos::deep_copy(bases_view, host_bases);
                    fastgatk::core::HostBatch count_host("collect-allelic-counts-v1");
                    count_host.records = observation_loci.size();
                    count_host.bytes = observation_loci.size() * (sizeof(std::size_t) + sizeof(std::uint8_t));
                    fastgatk::core::KernelPlan<ExecSpace> count_plan("collect-allelic-counts");
                    count_plan.begin_prepare(count_host);
                    fastgatk::core::DeviceBatch<ExecSpace> count_device(observation_loci.size());
                    count_device.bind("loci", loci_view);
                    count_device.bind("bases", bases_view);
                    count_device.bind("counts", device_counts);
                    ExecSpace().fence();
                    count_plan.end_prepare(count_device);
                    count_plan.begin_execute();
                    Kokkos::parallel_for("collect_allelic_counts_kernel", Kokkos::RangePolicy<ExecSpace>(0, observation_loci.size()),
                        KOKKOS_LAMBDA(const std::size_t i) {
                            Kokkos::atomic_add(&device_counts(loci_view(i) * BASES.size() + bases_view(i)),
                                               std::uint64_t{1});
                        });
                    ExecSpace().fence();
                    count_plan.end_execute();
                    result.kernel_batches = 1;
                    result.kernel_observations = observation_loci.size();
                    result.kernel_prepare_seconds = count_plan.telemetry().prepare_seconds;
                    result.kernel_execute_seconds = count_plan.telemetry().execute_seconds;
                    return result;
                },
                [](CountBatchResult result) -> std::optional<CountBatchResult> { return result; },
                [&](CountBatchResult result) {
                    reads_seen += result.reads_seen;
                    reads_used += result.reads_used;
                    reads_filtered += result.reads_filtered;
                    projected_bases += result.projected_bases;
                    count_kernel_batches += result.kernel_batches;
                    count_kernel_observations += result.kernel_observations;
                    count_kernel_prepare_seconds += result.kernel_prepare_seconds;
                    count_kernel_execute_seconds += result.kernel_execute_seconds;
                },
                [](const fastgatk::io::ReadBatch& batch) {
                    return std::max<std::uint64_t>(1, batch.bytes());
                },
                [](const CountBatchResult&) { return static_cast<std::uint64_t>(sizeof(CountBatchResult)); },
                [](const CountBatchResult&) { return static_cast<std::uint64_t>(sizeof(CountBatchResult)); });
            pipeline_metrics = pipeline.run();
            auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device_counts);
            final_counts.resize(loci.size() * BASES.size());
            for (std::size_t i = 0; i < final_counts.size(); ++i) final_counts[i] = host_counts(i);
        }
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write allelic counts: " + options.output);
        output << "@HD\tVN:1.6\n";
        for (std::size_t tid = 0; tid < reader.header().contigs.size(); ++tid)
            output << "@SQ\tSN:" << reader.header().contigs[tid] << "\tLN:" << reader.header().contig_lengths[tid] << '\n';
        output << "@RG\tID:GATKCopyNumber\tSM:" << options.sample << '\n';
        output << "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n";
        std::size_t rows = 0;
        for (std::size_t locus = 0; locus < loci.size(); ++locus) {
            const int reference_base = base_index(loci[locus].reference);
            if (reference_base < 0) continue;
            std::uint64_t total = 0;
            std::array<std::uint64_t, 4> counts{};
            for (std::size_t base = 0; base < BASES.size(); ++base) {
                counts[base] = final_counts[locus * BASES.size() + base];
                total += counts[base];
            }
            const auto reference_count = counts[static_cast<std::size_t>(reference_base)];
            const auto alternate_count = total - reference_count;
            char alternate_base = 'N';
            if (alternate_count != 0) {
                std::uint64_t best = 0;
                for (std::size_t base = 0; base < BASES.size(); ++base) {
                    if (static_cast<int>(base) == reference_base) continue;
                    if (counts[base] > best) {
                        best = counts[base];
                        alternate_base = BASES[base];
                    }
                }
            }
            output << reader.header().contigs.at(static_cast<std::size_t>(loci[locus].tid)) << '\t'
                   << (loci[locus].position + 1) << '\t' << reference_count << '\t'
                   << alternate_count << '\t' << loci[locus].reference << '\t' << alternate_base << '\n';
            ++rows;
        }
        output.close();
        if (!output || !file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete allelic-count output: " + options.output);
        if (!options.manifest.empty()) {
            std::error_code output_error;
            const auto output_bytes = std::filesystem::file_size(options.output, output_error);
            const auto wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"CollectAllelicCounts\","
                     << "\"implementation\":\"fastgatk-collect-allelic-counts\",\"status\":\"prototype\","
                     << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                     << "\",\"determinism\":\"strict\","
                     << "\"input\":\"" << json_escape(options.input) << "\",\"reference\":\""
                     << json_escape(options.reference) << "\",\"output\":\"" << json_escape(options.output)
                     << "\",\"sample\":\"" << json_escape(options.sample) << "\",\"intervals\":"
                     << intervals.size() << ",\"excluded_intervals\":" << excluded_interval_count
                     << ",\"loci\":" << loci.size() << ",\"rows\":" << rows
                     << ",\"reads_seen\":" << reads_seen << ",\"reads_used\":" << reads_used
                     << ",\"reads_filtered\":" << reads_filtered << ",\"projected_bases\":" << projected_bases
                     << ",\"minimum_mapping_quality\":" << options.minimum_mapping_quality
                     << ",\"minimum_base_quality\":" << options.minimum_base_quality
                     << ",\"include_duplicates\":" << (options.include_duplicates ? "true" : "false")
                     << ",\"read_filter_policy\":{\"wellformed\":"
                     << (options.wellformed_filter ? "true" : "false")
                     << ",\"mapped\":" << (options.mapped_filter ? "true" : "false")
                     << ",\"nonzero_reference_length\":"
                     << (options.nonzero_reference_filter ? "true" : "false")
                     << ",\"not_duplicate\":" << (options.duplicate_filter ? "true" : "false")
                     << ",\"mapping_quality\":"
                     << (options.mapping_quality_filter ? "true" : "false")
                     << ",\"passes_vendor_quality\":"
                     << (options.vendor_quality_filter ? "true" : "false")
                     << ",\"not_secondary\":" << (options.secondary_filter ? "true" : "false")
                     << ",\"not_supplementary\":"
                     << (options.supplementary_filter ? "true" : "false") << "}"
                     << ",\"disable_tool_default_read_filters\":"
                     << (options.disable_tool_default_read_filters ? "true" : "false")
                     << ",\"sample_source\":\"" << sample_source << "\"," 
                     << "\"htslib_backend\":\"" << json_escape(reader.backend_description()) << "\","
                     << "\"reference_n_only_rows_omitted\":true,"
                     << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                     << "\",\"kind\":\"allelic-counts\",\"complete\":"
                     << (file_complete(options.output) ? "true" : "false") << "}],"
                     << "\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\"" << ExecSpace::name()
                     << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_batches\":" << count_kernel_batches
                     << ",\"kernel_observations\":" << count_kernel_observations
                     << ",\"kernel_prepare_seconds\":" << count_kernel_prepare_seconds
                     << ",\"kernel_execute_seconds\":" << count_kernel_execute_seconds
                     << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\""
                     << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                     << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                     << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                     << ",\"pipeline_decoded_bytes\":" << pipeline_metrics.decoded_bytes
                     << ",\"pipeline_computed_bytes\":" << pipeline_metrics.computed_bytes
                     << ",\"pipeline_encoded_bytes\":" << pipeline_metrics.encoded_bytes
                     << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes
                     << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes
                     << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes
                     << ",\"output_bytes\":"
                     << (output_error ? 0 : output_bytes)
                     << ",\"wall_seconds\":" << wall_seconds << "}}\n";
            manifest.close();
            if (!manifest || !file_complete(options.manifest))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest: " + options.manifest);
        }
        fai_destroy(reference);
        reference = nullptr;
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"CollectAllelicCounts\",\"status\":\"prototype\",\"sample\":\""
                  << json_escape(options.sample) << "\",\"sample_source\":\"" << sample_source
                  << "\",\"rows\":"
                  << rows << ",\"loci\":" << loci.size() << ",\"excluded_intervals\":" << excluded_interval_count
                  << ",\"reads_seen\":" << reads_seen
                  << ",\"reads_used\":" << reads_used << ",\"projected_bases\":" << projected_bases
                  << ",\"kernel_batches\":" << count_kernel_batches
                  << ",\"kernel_observations\":" << count_kernel_observations
                  << ",\"kernel_execution_space\":\"" << ExecSpace::name() << "\""
                  << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\""
                  << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                  << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                  << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items << "}\n";
        return 0;
    } catch (const std::exception& error) {
        if (reference) fai_destroy(reference);
        if (initialized) Kokkos::finalize();
        std::cerr << error.what() << '\n';
        return 2;
    }
#endif
}
