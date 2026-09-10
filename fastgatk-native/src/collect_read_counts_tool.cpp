#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/hdf5_count_collection.hpp"
#include "optional_boolean.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    std::string output;
    std::string manifest;
    std::string sample = "UNKNOWN";
    // GATK CollectReadCounts defaults to the HDF5 SimpleCountCollection.  TSV
    // remains an explicit portable interchange mode, matching --format TSV.
    std::string format = "HDF5";
    int minimum_mapping_quality = 30;
    std::size_t batch_records = 4096;
    int threads = 1;
    bool include_duplicates = false;
    std::string interval_merging_rule = "OVERLAPPING_ONLY";
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    std::string read_validation_stringency = "SILENT";
    std::string tmp_dir;
    std::string verbosity = "INFO";
    bool quiet = false;
};

struct Interval {
    int tid = -1;
    std::int64_t begin = 0;
    std::int64_t end = 0;
    std::uint64_t count = 0;
};

// Result owned by the compute stage.  Counts stay on Host after the Kokkos
// reduction so the single sink can merge them in input order without sharing
// mutable interval state with the worker.
struct CountBatchResult {
    std::vector<std::uint64_t> counts;
    std::uint64_t records = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
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
            std::cout << "fastgatk-collect-read-counts (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 BAM/CRAM/SAM input\n"
                         "  -L, --intervals REGION           interval/file selector (repeatable)\n"
                         "  -XL, --exclude-intervals REGION  subtract interval/file selector (repeatable)\n"
                         "  -R, --reference FILE              CRAM reference\n"
                         "  -O, --output FILE                TSV/HDF5 read-count output\n"
                         "      --format TSV|HDF5           GATK SimpleCountCollection format (default HDF5)\n"
                         "      --sample NAME                sample metadata\n"
                         "      --minimum-mapping-quality N  default 30\n"
                         "      --interval-merging-rule R    OVERLAPPING_ONLY (default) or ALL\n"
                         "      --interval-padding N          must be 0 for CollectReadCounts\n"
                         "      --interval-exclusion-padding N must be 0 for CollectReadCounts\n"
                         "      --include-duplicates         retain duplicate reads\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || is_option(argument, "--input"))
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                is_option(argument, "--region") ? "--region" :
                is_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (argument == "-XL" || is_option(argument, "--exclude-intervals"))
            options.exclusions.push_back(require_value(index, argc, argv, argument,
                                                        "--exclude-intervals", "-XL"));
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--format"))
            options.format = require_value(index, argc, argv, argument, "--format");
        else if (is_option(argument, "--sample"))
            options.sample = require_value(index, argc, argv, argument, "--sample");
        else if (is_option(argument, "--minimum-mapping-quality"))
            options.minimum_mapping_quality = std::stoi(require_value(
                index, argc, argv, argument, "--minimum-mapping-quality"));
        else if (is_option(argument, "--batch-records"))
            options.batch_records = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--batch-records")));
        else if (is_option(argument, "--threads"))
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--include-duplicates")
            options.include_duplicates = true;
        else if (argument == "-imr" || is_option(argument, "--interval-merging-rule")) {
            options.interval_merging_rule = require_value(
                index, argc, argv, argument, "--interval-merging-rule", "-imr");
            for (auto& character : options.interval_merging_rule)
                character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
            if (options.interval_merging_rule != "OVERLAPPING_ONLY" &&
                options.interval_merging_rule != "ALL")
                throw std::invalid_argument(
                    "--interval-merging-rule must be OVERLAPPING_ONLY or ALL");
        } else if (argument == "-ip" || is_option(argument, "--interval-padding") ||
                   argument == "-ixp" || is_option(argument, "--interval-exclusion-padding")) {
            const bool exclusion = argument == "-ixp" ||
                argument.rfind("--interval-exclusion-padding", 0) == 0;
            const char* name = exclusion ? "--interval-exclusion-padding" : "--interval-padding";
            const auto value = std::stoll(require_value(index, argc, argv, argument, name,
                                                         exclusion ? "-ixp" : "-ip"));
            if (value < 0) throw std::invalid_argument(std::string(name) + " must be non-negative");
            if (value != 0)
                throw std::invalid_argument(exclusion
                    ? "Interval exclusion padding must be set to 0."
                    : "Interval padding must be set to 0.");
            if (exclusion) options.exclusion_padding = value;
            else options.interval_padding = value;
        } else if (argument == "--QUIET" || argument == "--quiet" ||
                   argument.rfind("--QUIET=", 0) == 0 || argument.rfind("--quiet=", 0) == 0) {
            options.quiet = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, argument.rfind("--QUIET", 0) == 0 ? "--QUIET" : "--quiet");
        } else if (argument == "--read-validation-stringency" || argument == "-VS" ||
                   is_option(argument, "--read-validation-stringency")) {
            options.read_validation_stringency = require_value(index, argc, argv, argument,
                "--read-validation-stringency", "-VS");
            std::transform(options.read_validation_stringency.begin(), options.read_validation_stringency.end(),
                           options.read_validation_stringency.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
            if (options.read_validation_stringency != "STRICT" &&
                options.read_validation_stringency != "LENIENT" &&
                options.read_validation_stringency != "SILENT")
                throw std::invalid_argument("--read-validation-stringency must be STRICT, LENIENT, or SILENT");
        } else if (argument == "--disable-bam-index-caching" || argument == "-DBIC" ||
                   argument.rfind("--disable-bam-index-caching=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(index, argc, argv, argument,
                argument == "-DBIC" ? "-DBIC" : "--disable-bam-index-caching");
        } else if (argument == "--tmp-dir" || is_option(argument, "--tmp-dir")) {
            options.tmp_dir = require_value(index, argc, argv, argument, "--tmp-dir");
        } else if (argument == "--use-jdk-deflater" || argument == "-jdk-deflater" ||
                   argument == "--use-jdk-inflater" || argument == "-jdk-inflater" ||
                   argument.rfind("--use-jdk-deflater=", 0) == 0 || argument.rfind("--use-jdk-inflater=", 0) == 0) {
            const bool deflater = argument.rfind("--use-jdk-deflater", 0) == 0 || argument == "-jdk-deflater";
            (void)fastgatk::native::parse_optional_boolean(index, argc, argv, argument,
                deflater ? (argument == "-jdk-deflater" ? "-jdk-deflater" : "--use-jdk-deflater")
                         : (argument == "-jdk-inflater" ? "-jdk-inflater" : "--use-jdk-inflater"));
        } else if (argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--disable-tool-default-read-filters") {
            // Compatibility flags.  The native default read-filter mask remains
            // explicit in the manifest; disabling it is not silently supported.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--VERBOSITY") ||
                   is_option(argument, "--seconds-between-progress-updates")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--VERBOSITY", 0) == 0 ? "--VERBOSITY" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
            if (argument == "--verbosity" || argument.rfind("--verbosity=", 0) == 0 ||
                argument == "--VERBOSITY" || argument.rfind("--VERBOSITY=", 0) == 0)
                options.verbosity = argument.find('=') == std::string::npos ? argv[index] : option_value(argument, argument.rfind("--VERBOSITY", 0) == 0 ? "--VERBOSITY" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.regions.empty()) throw std::invalid_argument("-L/--intervals is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    for (auto& character : options.format)
        character = static_cast<char>(std::toupper(static_cast<unsigned char>(character)));
    if (options.format != "TSV" && options.format != "HDF5")
        throw std::invalid_argument("--format must be TSV or HDF5");
    if (options.format == "HDF5" && !fastgatk::hdf5::available())
        throw std::invalid_argument("BACKEND_UNAVAILABLE: HDF5 support was not built; use --format TSV or --fallback");
    if (options.minimum_mapping_quality < 0 || options.minimum_mapping_quality > 255 ||
        options.batch_records == 0 || options.threads < 1)
        throw std::invalid_argument("invalid read-count resource/filter parameters");
    if (options.interval_merging_rule != "OVERLAPPING_ONLY")
        throw std::invalid_argument("Interval merging rule must be set to OVERLAPPING_ONLY.");
    return options;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto character : text) {
        if (character == '"' || character == '\\') escaped << '\\';
        if (character == '\n') escaped << "\\n";
        else if (character == '\r') escaped << "\\r";
        else if (character == '\t') escaped << "\\t";
        else escaped << character;
    }
    return escaped.str();
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

std::vector<Interval> merge_intervals(std::vector<Interval> intervals,
                                      const std::string& rule,
                                      std::size_t& merge_count) {
    std::stable_sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        if (left.tid != right.tid) return left.tid < right.tid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Interval> merged;
    merged.reserve(intervals.size());
    for (auto interval : intervals) {
        if (interval.end <= interval.begin) continue;
        if (!merged.empty() && merged.back().tid == interval.tid &&
            (interval.begin < merged.back().end ||
             (rule == "ALL" && interval.begin == merged.back().end))) {
            merged.back().end = std::max(merged.back().end, interval.end);
            ++merge_count;
        } else {
            interval.count = 0;
            merged.push_back(std::move(interval));
        }
    }
    return merged;
}

std::vector<Interval> subtract_intervals(const std::vector<Interval>& input,
                                         const std::vector<fastgatk::io::HtsInterval>& exclusions) {
    if (exclusions.empty()) return input;
    std::vector<Interval> result;
    for (const auto& source : input) {
        std::int64_t cursor = source.begin;
        for (const auto& exclusion : exclusions) {
            if (exclusion.tid != source.tid) continue;
            if (exclusion.end <= cursor) continue;
            if (exclusion.start >= source.end) break;
            if (exclusion.start > cursor)
                result.push_back(Interval{source.tid, cursor,
                                          std::min(source.end, exclusion.start), 0});
            if (exclusion.end >= source.end) { cursor = source.end; break; }
            cursor = std::max(cursor, exclusion.end);
        }
        if (cursor < source.end) result.push_back(Interval{source.tid, cursor, source.end, 0});
    }
    return result;
}

void write_tsv(const std::string& path, const Options& options,
               const fastgatk::io::HeaderSummary& header,
               const std::vector<Interval>& intervals) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write counts: " + path);
    output << "@HD\tVN:1.6\tSO:coordinate\n";
    for (std::size_t tid = 0; tid < header.contigs.size(); ++tid)
        output << "@SQ\tSN:" << header.contigs[tid] << "\tLN:" << header.contig_lengths[tid] << '\n';
    output << "@RG\tID:fastgatk\tSM:" << options.sample << '\n';
    output << "CONTIG\tSTART\tEND\tCOUNT\n";
    for (const auto& interval : intervals)
        output << header.contigs.at(static_cast<std::size_t>(interval.tid)) << '\t'
               << (interval.begin + 1) << '\t' << interval.end << '\t'
               << interval.count << '\n';
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize counts: " + path);
}

void write_hdf5(const std::string& path, const Options& options,
                const fastgatk::io::HeaderSummary& header,
                const std::vector<Interval>& intervals) {
    fastgatk::hdf5::CountCollection collection;
    collection.sample = options.sample;
    // HDF5SimpleCountCollection stores the SAMTextHeaderCodec rendering of a
    // sequence dictionary, not the complete input BAM header.  HtsReader
    // preserves SQ tags (AS/M5/UR/SP) and normalizes the HD line to the
    // HTSJDK VN:1.6 form, so Java GATK can round-trip the same dictionary.
    collection.sequence_dictionary = header.sequence_dictionary;
    if (collection.sequence_dictionary.empty()) {
        collection.sequence_dictionary = "@HD\tVN:1.6\n";
        for (std::size_t tid = 0; tid < header.contigs.size(); ++tid)
            collection.sequence_dictionary += "@SQ\tSN:" + header.contigs[tid] + "\tLN:" +
                std::to_string(header.contig_lengths[tid]) + "\n";
    }
    collection.intervals.reserve(intervals.size());
    collection.counts.reserve(intervals.size());
    for (const auto& interval : intervals) {
        collection.intervals.push_back({header.contigs.at(static_cast<std::size_t>(interval.tid)),
                                        interval.begin + 1, interval.end});
        collection.counts.push_back(static_cast<double>(interval.count));
    }
    fastgatk::hdf5::write_simple_count_collection(path, collection);
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const auto wall_begin = std::chrono::steady_clock::now();
        const auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;

        const fastgatk::runtime::AdaptiveController controller;
        const auto initial_limits = controller.initial(
            resources.budget(), fastgatk::runtime::WorkEstimate{
                static_cast<std::uint32_t>(std::min<std::size_t>(
                    options.batch_records, std::numeric_limits<std::uint32_t>::max())),
                1, 2048, 0, 2048});
        const auto initial_batch_records = std::min<std::size_t>(
            options.batch_records, initial_limits.max_reads);
        fastgatk::io::HtsReader reader(options.input, options.reference, options.regions,
                                       initial_batch_records, {}, options.exclusions,
                                       options.interval_padding, options.exclusion_padding);
        if (!reader.compiled()) throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib is not available");
        std::vector<Interval> requested_intervals;
        for (const auto& interval : reader.requested_intervals())
            requested_intervals.push_back(Interval{interval.tid, interval.start, interval.end, 0});
        if (requested_intervals.empty()) throw std::runtime_error("BAD_INPUT: no intervals were selected");
        const auto requested_interval_count = requested_intervals.size();
        const auto excluded_interval_count = reader.exclusion_intervals().size();
        requested_intervals = subtract_intervals(requested_intervals, reader.exclusion_intervals());
        if (requested_intervals.empty())
            throw std::runtime_error("BAD_INPUT: no intervals remain after applying --exclude-intervals");
        std::size_t interval_merges = 0;
        auto intervals = merge_intervals(
            std::move(requested_intervals), options.interval_merging_rule, interval_merges);
        if (intervals.empty()) throw std::runtime_error("BAD_INPUT: no non-empty intervals were selected");

        // Interval metadata is invariant across read batches.  Keep one
        // device-resident copy so the hot loop only stages the changing read
        // columns and the per-batch histogram; this matters for small batches
        // and avoids allocator/PCIe churn on accelerator backends.
        Kokkos::View<std::int32_t*> interval_tids("read_count_interval_tids", intervals.size());
        Kokkos::View<std::int64_t*> interval_begins("read_count_interval_begins", intervals.size());
        Kokkos::View<std::int64_t*> interval_ends("read_count_interval_ends", intervals.size());
        auto host_interval_tids = Kokkos::create_mirror_view(interval_tids);
        auto host_interval_begins = Kokkos::create_mirror_view(interval_begins);
        auto host_interval_ends = Kokkos::create_mirror_view(interval_ends);
        for (std::size_t i = 0; i < intervals.size(); ++i) {
            host_interval_tids(i) = intervals[i].tid;
            host_interval_begins(i) = intervals[i].begin;
            host_interval_ends(i) = intervals[i].end;
        }
        Kokkos::deep_copy(interval_tids, host_interval_tids);
        Kokkos::deep_copy(interval_begins, host_interval_begins);
        Kokkos::deep_copy(interval_ends, host_interval_ends);

        std::uint64_t reads_seen = 0;
        std::uint64_t reads_used = 0;
        std::uint64_t reads_filtered = 0;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        using ReadBatch = fastgatk::io::ReadBatch;
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<ReadBatch, CountBatchResult, CountBatchResult>;
        std::uint64_t count_kernel_batches = 0;
        std::uint64_t count_kernel_records = 0;
        double count_kernel_prepare_seconds = 0.0;
        double count_kernel_execute_seconds = 0.0;
        auto batch_limits = initial_limits;
        batch_limits.max_reads = static_cast<std::uint32_t>(std::min<std::size_t>(
            reader.batch_records(), std::numeric_limits<std::uint32_t>::max()));
        std::uint64_t adaptive_reductions = 0;
        const auto saturating_multiply = [](std::uint64_t left, std::uint64_t right) {
            return right != 0 && left > std::numeric_limits<std::uint64_t>::max() / right
                ? std::numeric_limits<std::uint64_t>::max() : left * right;
        };
        const auto saturating_add = [](std::uint64_t left, std::uint64_t right) {
            return right > std::numeric_limits<std::uint64_t>::max() - left
                ? std::numeric_limits<std::uint64_t>::max() : left + right;
        };
        const auto safe_budget = resources.safe_memory_budget_bytes();
        const auto decoded_capacity = safe_budget != 0
            ? safe_budget
            : std::max<std::uint64_t>(64ULL << 20,
                saturating_multiply(static_cast<std::uint64_t>(initial_batch_records), 16384));
        const auto result_bytes = saturating_add(
            static_cast<std::uint64_t>(sizeof(CountBatchResult)),
            saturating_multiply(static_cast<std::uint64_t>(intervals.size()), sizeof(std::uint64_t)));
        const auto result_capacity = std::max<std::uint64_t>(sizeof(CountBatchResult), result_bytes);
        const auto pipeline_limits = Pipeline::Limits{
            decoded_capacity, result_capacity, result_capacity};
        const auto pipeline_metrics = Pipeline(
            pipeline_limits,
            [&]() -> std::optional<ReadBatch> {
                ReadBatch batch;
                if (!reader.next(batch)) return std::nullopt;
                reads_seen += batch.records();
                fastgatk::runtime::RuntimeTelemetry runtime;
                runtime.host_bytes = batch.bytes();
                runtime.inflight_bytes = runtime.host_bytes;
                runtime.compute_queue_empty = true;
                const auto pressure = safe_budget != 0 && runtime.host_bytes > safe_budget
                    ? fastgatk::runtime::Pressure::HostMemory
                    : fastgatk::runtime::Pressure::Normal;
                const auto next_limits = controller.next(runtime, batch_limits, pressure);
                if (next_limits.max_reads < reader.batch_records()) {
                    reader.set_batch_records(next_limits.max_reads);
                    ++adaptive_reductions;
                }
                batch_limits = next_limits;
                return batch;
            },
            [&](ReadBatch batch) -> std::optional<CountBatchResult> {
                const auto records = batch.records();
                fastgatk::core::HostBatch count_host("collect-read-counts-v1");
                count_host.records = records;
                count_host.bytes = records * (sizeof(std::int32_t) + sizeof(std::int64_t) +
                                              sizeof(std::uint16_t) + sizeof(std::uint8_t));
                fastgatk::core::KernelPlan<ExecSpace> count_plan("collect-read-counts");
                count_plan.begin_prepare(count_host);
                Kokkos::View<std::int32_t*> tids("read_count_tids", records);
                Kokkos::View<std::int64_t*> positions("read_count_positions", records);
                Kokkos::View<std::uint16_t*> flags("read_count_flags", records);
                Kokkos::View<std::uint8_t*> mapq("read_count_mapq", records);
                Kokkos::View<std::uint64_t*> counts("read_count_intervals", intervals.size());
                auto host_tids = Kokkos::create_mirror_view(tids);
                auto host_positions = Kokkos::create_mirror_view(positions);
                auto host_flags = Kokkos::create_mirror_view(flags);
                auto host_mapq = Kokkos::create_mirror_view(mapq);
                for (std::size_t record = 0; record < records; ++record) {
                    host_tids(record) = record < batch.tids.size() ? batch.tids[record] : -1;
                    host_positions(record) = record < batch.positions.size() ? batch.positions[record] : -1;
                    host_flags(record) = record < batch.flags.size() ? batch.flags[record] : 0xffffU;
                    host_mapq(record) = record < batch.mapq.size() ? batch.mapq[record] : 0;
                }
                Kokkos::deep_copy(tids, host_tids);
                Kokkos::deep_copy(positions, host_positions);
                Kokkos::deep_copy(flags, host_flags);
                Kokkos::deep_copy(mapq, host_mapq);
                Kokkos::deep_copy(counts, std::uint64_t{0});
                fastgatk::core::DeviceBatch<ExecSpace> count_device(records);
                count_device.bind("tids", tids);
                count_device.bind("positions", positions);
                count_device.bind("flags", flags);
                count_device.bind("mapq", mapq);
                count_device.bind("interval_tids", interval_tids);
                count_device.bind("interval_begins", interval_begins);
                count_device.bind("interval_ends", interval_ends);
                count_device.bind("counts", counts);
                ExecSpace().fence();
                count_plan.end_prepare(count_device);
                const auto interval_count = intervals.size();
                const auto minimum_mapping_quality = static_cast<std::uint8_t>(options.minimum_mapping_quality);
                const auto include_duplicates = options.include_duplicates;
                count_plan.begin_execute();
                Kokkos::parallel_for("collect_read_counts", Kokkos::RangePolicy<ExecSpace>(0, records), KOKKOS_LAMBDA(const std::size_t record) {
                    const auto flag = flags(record);
                    if (tids(record) < 0 || positions(record) < 0 ||
                        (flag & (0x4U | 0x100U | 0x200U)) != 0 ||
                        (!include_duplicates && (flag & 0x400U) != 0) ||
                        mapq(record) < minimum_mapping_quality || mapq(record) == 255) return;
                    for (std::size_t interval = 0; interval < interval_count; ++interval) {
                        if (tids(record) == interval_tids(interval) &&
                            positions(record) >= interval_begins(interval) &&
                            positions(record) < interval_ends(interval)) {
                            Kokkos::atomic_add(&counts(interval), std::uint64_t{1});
                            break;
                        }
                    }
                });
                ExecSpace().fence();
                count_plan.end_execute();
                auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, counts);
                CountBatchResult result;
                result.records = records;
                result.counts.resize(intervals.size());
                for (std::size_t interval = 0; interval < intervals.size(); ++interval)
                    result.counts[interval] = host_counts(interval);
                result.prepare_seconds = count_plan.telemetry().prepare_seconds;
                result.execute_seconds = count_plan.telemetry().execute_seconds;
                return result;
            },
            [](CountBatchResult result) -> std::optional<CountBatchResult> {
                return result;
            },
            [&](CountBatchResult result) {
                ++count_kernel_batches;
                count_kernel_records += result.records;
                count_kernel_prepare_seconds += result.prepare_seconds;
                count_kernel_execute_seconds += result.execute_seconds;
                for (std::size_t interval = 0; interval < intervals.size(); ++interval) {
                    intervals[interval].count += result.counts[interval];
                    reads_used += result.counts[interval];
                }
            },
            [](const ReadBatch& batch) { return batch.bytes(); },
            [&](const CountBatchResult& result) {
                return saturating_add(static_cast<std::uint64_t>(sizeof(CountBatchResult)),
                                      saturating_multiply(static_cast<std::uint64_t>(result.counts.size()), sizeof(std::uint64_t)));
            },
            [&](const CountBatchResult& result) {
                return saturating_add(static_cast<std::uint64_t>(sizeof(CountBatchResult)),
                                      saturating_multiply(static_cast<std::uint64_t>(result.counts.size()), sizeof(std::uint64_t)));
            }).run();
        reads_filtered = reads_seen > reads_used ? reads_seen - reads_used : 0;
        if (options.format == "HDF5") write_hdf5(options.output, options, reader.header(), intervals);
        else write_tsv(options.output, options, reader.header(), intervals);
        if (!file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: counts output is missing or empty");
        std::error_code output_error;
        const auto output_bytes = std::filesystem::file_size(options.output, output_error);
        if (output_error || output_bytes == 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: counts output size is unavailable");
        bool hdf5_roundtrip = false;
        if (options.format == "HDF5") {
            const auto roundtrip = fastgatk::hdf5::read_simple_count_collection(options.output);
            if (roundtrip.sample != options.sample || roundtrip.intervals.size() != intervals.size() ||
                roundtrip.counts.size() != intervals.size())
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: HDF5 count metadata round-trip mismatch");
            for (std::size_t index = 0; index < intervals.size(); ++index) {
                const auto& interval = intervals[index];
                const auto& stored = roundtrip.intervals[index];
                const auto contig = reader.header().contigs.at(static_cast<std::size_t>(interval.tid));
                if (stored.contig != contig || stored.start != interval.begin + 1 ||
                    stored.end != interval.end || roundtrip.counts[index] != static_cast<double>(interval.count))
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: HDF5 interval/count round-trip mismatch");
            }
            hdf5_roundtrip = true;
        }
        const auto wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_begin).count();

        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"CollectReadCounts\","
                 << "\"implementation\":\"fastgatk-collect-read-counts\",\"status\":\"prototype\","
                 << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\","
                 << "\"determinism\":\"strict\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\","
                 << "\"primary_output_kind\":\"" << (options.format == "HDF5" ? "hdf5-simple-count-collection" : "tsv-simple-count-collection") << "\","
                 << "\"compatibility\":{\"tsv\":true,\"hdf5\":" << (options.format == "HDF5" ? "true" : "false") << ",\"read_start_semantics\":true,"
                 << "\"interval_merging_rule\":\"" << json_escape(options.interval_merging_rule) << "\","
                 << "\"exclude_intervals\":true,\"zero_interval_padding_required\":true,"
                 << "\"quiet\":" << (options.quiet ? "true" : "false")
                 << ",\"read_validation_stringency\":\""
                 << json_escape(options.read_validation_stringency) << "\",\"tmp_dir\":\""
                 << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                 << json_escape(options.verbosity) << "\","
                 << "\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"" << (options.format == "HDF5" ? "hdf5" : "tsv") << "\",\"complete\":true}],"
                 << "\"telemetry\":{\"resources\":{\"backend\":\""
                 << json_escape(reader.backend_description()) << "\"},\"intervals\":" << intervals.size()
                 << ",\"requested_intervals\":" << requested_interval_count
                 << ",\"excluded_intervals\":" << excluded_interval_count
                 << ",\"interval_merges\":" << interval_merges
                 << ",\"interval_file_inputs\":" << reader.interval_file_inputs()
                 << ",\"interval_file_records\":" << reader.interval_file_records()
                 << ",\"reads_seen\":" << reads_seen << ",\"reads_used\":" << reads_used
                 << ",\"reads_filtered\":" << reads_filtered
                 << ",\"quiet\":" << (options.quiet ? "true" : "false")
                 << ",\"read_validation_stringency\":\""
                 << json_escape(options.read_validation_stringency) << "\",\"tmp_dir\":\""
                 << json_escape(options.tmp_dir) << "\",\"verbosity\":\""
                 << json_escape(options.verbosity) << "\""
                 << ",\"count_kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"count_kernel_execution_space\":\"" << ExecSpace::name()
                 << "\",\"count_kernel_execution_policy\":\"RangePolicy\""
                 << ",\"count_kernel_batches\":" << count_kernel_batches
                 << ",\"count_kernel_records\":" << count_kernel_records
                 << ",\"count_kernel_prepare_seconds\":" << count_kernel_prepare_seconds
                 << ",\"count_kernel_execute_seconds\":" << count_kernel_execute_seconds
                 << ",\"pipeline_lifecycle\":\"decode->compute->encode->sink\""
                 << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                 << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                 << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                 << ",\"pipeline_decoded_bytes\":" << pipeline_metrics.decoded_bytes
                 << ",\"pipeline_computed_bytes\":" << pipeline_metrics.computed_bytes
                 << ",\"pipeline_encoded_bytes\":" << pipeline_metrics.encoded_bytes
                 << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes
                 << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes
                 << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes
                 << ",\"initial_batch_records\":" << initial_batch_records
                 << ",\"effective_batch_records\":" << reader.batch_records()
                 << ",\"adaptive_batch_reductions\":" << adaptive_reductions
                 << ",\"minimum_mapping_quality\":" << options.minimum_mapping_quality
                 << ",\"exclude_duplicates\":" << (options.include_duplicates ? "false" : "true")
                 << ",\"output_bytes\":" << output_bytes
                 << ",\"hdf5_roundtrip\":" << (hdf5_roundtrip ? "true" : "false")
                 << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}}\n";
        manifest.close();
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize manifest");
        std::cout << "{\"tool\":\"CollectReadCounts\",\"status\":\"prototype\","
                  << "\"format\":\"" << options.format << "\",\"intervals\":" << intervals.size()
                  << ",\"interval_merging_rule\":\"" << json_escape(options.interval_merging_rule) << "\""
                  << ",\"requested_intervals\":" << requested_interval_count
                  << ",\"excluded_intervals\":" << excluded_interval_count
                  << ",\"interval_merges\":" << interval_merges
                  << ",\"reads_seen\":" << reads_seen << ",\"reads_used\":" << reads_used
                  << ",\"reads_filtered\":" << reads_filtered
                  << ",\"count_kernel_batches\":" << count_kernel_batches
                  << ",\"count_kernel_records\":" << count_kernel_records
                  << ",\"count_kernel_execution_space\":\"" << ExecSpace::name() << "\""
                  << ",\"pipeline_lifecycle\":\"decode->compute->encode->sink\""
                  << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                  << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                  << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                  << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes
                  << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes
                  << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes
                  << ",\"output_bytes\":" << output_bytes
                  << ",\"hdf5_roundtrip\":" << (hdf5_roundtrip ? "true" : "false")
                  << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}\n";
        // Kokkos requires every View allocation to be released before
        // finalize.  These interval views deliberately outlive the batch loop
        // for reuse, so drop both device and mirror handles explicitly here.
        interval_tids = {};
        interval_begins = {};
        interval_ends = {};
        host_interval_tids = {};
        host_interval_begins = {};
        host_interval_ends = {};
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
