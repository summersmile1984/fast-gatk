#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#if FASTGATK_HAS_HTSLIB
#include <htslib/faidx.h>
#endif

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    std::string output;
    std::string manifest;
    std::string sample;
    std::size_t batch_records = 4096;
    std::size_t max_loci = 50000000;
    int min_base_quality = 0;
    int max_base_quality = 127;
    int start = 1;
    int stop = 500;
    int n_bins = 499;
    int summary_coverage_threshold = 15;
    int threads = 1;
    bool print_base_counts = false;
    bool omit_locus_table = false;
    bool omit_interval_statistics = false;
    bool omit_per_sample_statistics = false;
    bool omit_depth_output_at_each_base = false;
    bool include_deletions = false;
    bool ignore_deletion_sites = false;
    bool include_ref_n_sites = false;
    bool count_fragments = false;
};

struct Target {
    std::int32_t tid = -1;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

// Per-ReadBatch counters cross the encode/sink boundary so all aggregate
// telemetry remains deterministic. The compute stage owns the Host CIGAR /
// fragment state and submits only fixed observation arrays to Kokkos.
struct CoverageBatchResult {
    std::uint64_t reads_seen = 0;
    std::uint64_t reads_used = 0;
    std::uint64_t reads_filtered = 0;
    std::uint64_t observations = 0;
    std::uint64_t kernel_batches = 0;
    std::uint64_t kernel_observations = 0;
    double kernel_prepare_seconds = 0.0;
    double kernel_execute_seconds = 0.0;
};

// COUNT_FRAGMENTS needs a stable Host-side identity across read batches.  The
// locus is part of the key because one fragment contributes once at every
// covered locus; the read name is the same for both mates in a pair.
struct FragmentObservationKey {
    std::uint32_t sample = 0;
    std::uint64_t locus = 0;
    std::string name;

    bool operator==(const FragmentObservationKey& other) const {
        return sample == other.sample && locus == other.locus && name == other.name;
    }
};

struct FragmentObservationKeyHash {
    std::size_t operator()(const FragmentObservationKey& key) const noexcept {
        std::size_t hash = std::hash<std::string>{}(key.name);
        hash ^= std::hash<std::uint64_t>{}(key.locus) + static_cast<std::size_t>(0x9e3779b9) +
                (hash << 6) + (hash >> 2);
        hash ^= std::hash<std::uint32_t>{}(key.sample) + static_cast<std::size_t>(0x9e3779b9) +
                (hash << 6) + (hash >> 2);
        return hash;
    }
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
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
            std::cout << "fastgatk-depth-of-coverage (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                    BAM/CRAM/SAM input\n"
                         "  -R, --reference FILE                reference FASTA\n"
                         "  -L, --intervals REGION               interval selector (repeatable)\n"
                         "  -XL, --exclude-intervals REGION      exclude selector (repeatable)\n"
                         "      --interval-set-rule RULE         UNION (default) or INTERSECTION\n"
                         "      --interval-padding N             padding for included intervals\n"
                         "      --interval-exclusion-padding N   padding for excluded intervals\n"
                         "  -O, --output FILE                   locus output prefix/path\n"
                         "      --sample NAME                   sample name override\n"
                         "      --min-base-quality N            default 0\n"
                         "      --max-base-quality N            default 127\n"
                         "      --include-deletions             count deletion pileups\n"
                         "      --ignore-deletion-sites         do not count deletion pileups\n"
                         "      --include-ref-n-sites          include loci whose reference base is N\n"
                         "      --count-type TYPE              COUNT_READS or COUNT_FRAGMENTS\n"
                         "      --print-base-counts             append A/C/G/T/N counts\n"
                         "      --omit-locus-table              omit the main locus table\n"
                         "      --omit-depth-output-at-each-base alias for --omit-locus-table\n"
                         "      --omit-interval-statistics       omit interval histogram\n"
                         "      --omit-per-sample-statistics     omit sample histogram\n"
                         "      --omit-intervals                legacy alias for --omit-interval-statistics\n"
                         "      --omit-sample-summary           legacy alias for --omit-per-sample-statistics\n"
                         "      --omit-genes                    accepted legacy no-op (no gene sidecar)\n"
                         "      --summary-coverage-threshold N   default 15\n"
                         "      --start N --stop N --nBins N     cumulative histogram range\n"
                         "      --max-loci N                    memory guard (default 50000000)\n"
                         "      --output-manifest FILE          OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || has_option(argument, "--input")) {
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                has_option(argument, "--region") ? "--region" : "--interval";
            options.regions.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals") ||
                   has_option(argument, "--exclude-interval") || has_option(argument, "--exclude-region")) {
            const char* name = argument == "-XL" ? "--exclude-intervals" :
                has_option(argument, "--exclude-region") ? "--exclude-region" : "--exclude-interval";
            options.exclusions.push_back(require_value(index, argc, argv, argument, name, "-XL"));
        } else if (argument == "-isr" || has_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        } else if (argument == "-ip" || has_option(argument, "--interval-padding")) {
            const auto value = require_value(index, argc, argv, argument, "--interval-padding", "-ip");
            try { options.interval_padding = std::stoll(value); }
            catch (...) { throw std::invalid_argument("invalid --interval-padding: " + value); }
            if (options.interval_padding < 0) throw std::invalid_argument("--interval-padding must be non-negative");
        } else if (argument == "-ixp" || has_option(argument, "--interval-exclusion-padding")) {
            const auto value = require_value(index, argc, argv, argument,
                "--interval-exclusion-padding", "-ixp");
            try { options.exclusion_padding = std::stoll(value); }
            catch (...) { throw std::invalid_argument("invalid --interval-exclusion-padding: " + value); }
            if (options.exclusion_padding < 0)
                throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (has_option(argument, "--sample")) {
            options.sample = require_value(index, argc, argv, argument, "--sample");
        } else if (has_option(argument, "--min-base-quality") || has_option(argument, "--minimum-base-quality")) {
            options.min_base_quality = std::stoi(require_value(index, argc, argv, argument,
                argument.rfind("--minimum", 0) == 0 ? "--minimum-base-quality" : "--min-base-quality"));
        } else if (has_option(argument, "--max-base-quality")) {
            options.max_base_quality = std::stoi(require_value(index, argc, argv, argument, "--max-base-quality"));
        } else if (has_option(argument, "--batch-records")) {
            options.batch_records = static_cast<std::size_t>(std::stoull(require_value(index, argc, argv, argument, "--batch-records")));
        } else if (has_option(argument, "--max-loci")) {
            options.max_loci = static_cast<std::size_t>(std::stoull(require_value(index, argc, argv, argument, "--max-loci")));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--start")) {
            options.start = std::stoi(require_value(index, argc, argv, argument, "--start"));
        } else if (has_option(argument, "--stop")) {
            options.stop = std::stoi(require_value(index, argc, argv, argument, "--stop"));
        } else if (has_option(argument, "--nBins") || has_option(argument, "--nbins")) {
            options.n_bins = std::stoi(require_value(index, argc, argv, argument, "--nBins"));
        } else if (has_option(argument, "--summary-coverage-threshold")) {
            options.summary_coverage_threshold = std::stoi(require_value(index, argc, argv, argument, "--summary-coverage-threshold"));
        } else if (argument == "--print-base-counts") {
            options.print_base_counts = true;
        } else if (argument == "--include-deletions") {
            options.include_deletions = true;
        } else if (argument == "--ignore-deletion-sites" ||
                   argument.rfind("--ignore-deletion-sites=", 0) == 0) {
            options.ignore_deletion_sites = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--ignore-deletion-sites");
        } else if (argument == "--omit-locus-table" ||
                   argument.rfind("--omit-locus-table=", 0) == 0) {
            options.omit_locus_table = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-locus-table");
        } else if (argument == "--omit-depth-output-at-each-base" ||
                   argument.rfind("--omit-depth-output-at-each-base=", 0) == 0) {
            options.omit_depth_output_at_each_base = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-depth-output-at-each-base");
        } else if (argument == "--omit-interval-statistics" ||
                   argument.rfind("--omit-interval-statistics=", 0) == 0) {
            options.omit_interval_statistics = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-interval-statistics");
        } else if (argument == "--omit-per-sample-statistics" ||
                   argument.rfind("--omit-per-sample-statistics=", 0) == 0) {
            options.omit_per_sample_statistics = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-per-sample-statistics");
        } else if (argument == "--include-ref-n-sites") {
            options.include_ref_n_sites = true;
        } else if (argument == "--omit-genes") {
            // Gene aggregation is not part of the bounded native output
            // contract; accepting this legacy switch is therefore a safe
            // no-op (there is no gene sidecar to suppress).
        } else if (argument == "--omit-intervals" ||
                   argument.rfind("--omit-intervals=", 0) == 0) {
            // Older GATK wrappers used --omit-intervals for the same output
            // family now named --omit-interval-statistics.
            options.omit_interval_statistics = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-intervals");
        } else if (argument == "--omit-sample-summary" ||
                   argument.rfind("--omit-sample-summary=", 0) == 0) {
            // Compatibility alias used by pre-4.x DepthOfCoverage wrappers.
            options.omit_per_sample_statistics = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--omit-sample-summary");
        } else if (argument == "--disable-sequence-dictionary-validation" || argument == "--quiet") {
            // The bounded native path has no additional dictionary validation
            // or progress stream to disable.
        } else if (has_option(argument, "--partition-type")) {
            const auto partition = require_value(index, argc, argv, argument, "--partition-type");
            if (partition != "sample")
                throw std::invalid_argument("UNSUPPORTED: native DepthOfCoverage currently supports --partition-type sample only");
        } else if (has_option(argument, "--count-type")) {
            const auto type = require_value(index, argc, argv, argument, "--count-type");
            if (type == "COUNT_READS") options.count_fragments = false;
            else if (type == "COUNT_FRAGMENTS") options.count_fragments = true;
            else throw std::invalid_argument("UNSUPPORTED: native DepthOfCoverage supports COUNT_READS or COUNT_FRAGMENTS");
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
    if (options.input.empty() || options.reference.empty() || options.output.empty())
        throw std::invalid_argument("-I/--input, -R/--reference, and -O/--output are required");
    if (options.min_base_quality < 0 || options.max_base_quality > 127 ||
        options.min_base_quality > options.max_base_quality || options.batch_records == 0 ||
        options.max_loci == 0 || options.threads < 1 || options.start < 0 ||
        options.stop < options.start || options.n_bins < 1 || options.summary_coverage_threshold < 0)
        throw std::invalid_argument("invalid DepthOfCoverage resource/filter parameters");
    options.omit_locus_table = options.omit_locus_table || options.omit_depth_output_at_each_base;
    return options;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const char c : text) {
        if (c == '"' || c == '\\') escaped << '\\';
        if (c == '\n') escaped << "\\n";
        else if (c == '\r') escaped << "\\r";
        else if (c == '\t') escaped << "\\t";
        else escaped << c;
    }
    return escaped.str();
}

std::uintmax_t file_bytes(const std::string& path) {
    std::error_code error;
    const auto result = std::filesystem::file_size(path, error);
    return error ? 0 : result;
}

std::size_t find_locus(const std::vector<Target>& targets,
                       const std::vector<std::size_t>& offsets,
                       std::int32_t tid, std::int64_t position) {
    std::size_t low = 0, high = targets.size();
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        const auto& target = targets[middle];
        if (target.tid < tid || (target.tid == tid && target.end <= position)) low = middle + 1;
        else high = middle;
    }
    if (low >= targets.size() || targets[low].tid != tid ||
        position < targets[low].start || position >= targets[low].end) return std::numeric_limits<std::size_t>::max();
    return offsets[low] + static_cast<std::size_t>(position - targets[low].start);
}

// Apply GATK's include-then-exclude traversal boundary before allocating the
// dense locus matrices.  Keeping excluded bases out of the materialized target
// list (rather than merely zeroing their counters) is important: locus output,
// interval summaries, quantiles and cumulative denominators must all observe
// the same set of coordinates.  The input vectors are sorted/normalized by
// HtsReader, so this remains linear in the number of include/exclude spans.
std::vector<Target> subtract_exclusions(const std::vector<Target>& includes,
                                        const std::vector<fastgatk::io::HtsInterval>& exclusions) {
    std::vector<Target> result;
    for (const auto& include : includes) {
        std::int64_t cursor = include.start;
        for (const auto& exclusion : exclusions) {
            if (exclusion.tid < include.tid) continue;
            if (exclusion.tid > include.tid) break;
            if (exclusion.end <= cursor) continue;
            if (exclusion.start >= include.end) break;
            if (exclusion.start > cursor)
                result.push_back(Target{include.tid, cursor,
                                        std::min(include.end, exclusion.start)});
            cursor = std::max(cursor, exclusion.end);
            if (cursor >= include.end) break;
        }
        if (cursor < include.end)
            result.push_back(Target{include.tid, cursor, include.end});
    }
    return result;
}

std::string base_counts(const std::uint64_t* counts, bool include_deletions) {
    std::ostringstream out;
    const char* names = "ACGTN";
    for (int i = 0; i < 5; ++i) {
        if (i) out << ' ';
        out << names[i] << ':' << counts[i];
    }
    if (include_deletions) out << " D:" << counts[5];
    return out.str();
}

std::string format_two(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string format_one(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << value;
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
#if !FASTGATK_HAS_HTSLIB
    (void)argc; (void)argv;
    std::cerr << "BACKEND_UNAVAILABLE: build with HTSlib to run DepthOfCoverage\n";
    return 2;
#else
    bool initialized = false;
    faidx_t* reference_fai = nullptr;
    try {
        const auto wall_begin = std::chrono::steady_clock::now();
        auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;

        fastgatk::io::HtsReader reader(
            options.input, options.reference, options.regions, options.batch_records,
            options.sample, options.exclusions, options.interval_padding,
            options.exclusion_padding, options.interval_set_rule);
        if (!reader.compiled()) throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib is not available");
        std::vector<std::string> samples;
        if (!options.sample.empty()) {
            if (!reader.header().samples.empty() &&
                std::find(reader.header().samples.begin(), reader.header().samples.end(), options.sample) ==
                    reader.header().samples.end())
                throw std::runtime_error("BAD_INPUT: --sample is not present in the BAM @RG SM header: " + options.sample);
            samples.push_back(options.sample);
        } else if (!reader.header().samples.empty()) {
            samples = reader.header().samples;
        } else {
            samples.push_back("UNKNOWN");
        }
        const std::size_t sample_count = samples.size();
        std::map<std::string, std::size_t> sample_index;
        for (std::size_t index = 0; index < samples.size(); ++index)
            sample_index.emplace(samples[index], index);

        std::vector<Target> included_targets;
        if (!options.regions.empty()) {
            for (const auto& interval : reader.intervals())
                included_targets.push_back(Target{interval.tid, interval.start, interval.end});
        } else {
            for (std::size_t tid = 0; tid < reader.header().contigs.size(); ++tid)
                included_targets.push_back(Target{static_cast<std::int32_t>(tid), 0, reader.header().contig_lengths[tid]});
        }
        const auto targets = subtract_exclusions(included_targets, reader.exclusion_intervals());
        if (targets.empty()) throw std::runtime_error("BAD_INPUT: no loci were selected");
        std::vector<std::size_t> offsets(targets.size() + 1, 0);
        for (std::size_t i = 0; i < targets.size(); ++i) {
            if (targets[i].end < targets[i].start ||
                static_cast<std::uint64_t>(targets[i].end - targets[i].start) > options.max_loci - offsets[i])
                throw std::runtime_error("RESOURCE_EXHAUSTED: selected loci exceed --max-loci memory guard");
            offsets[i + 1] = offsets[i] + static_cast<std::size_t>(targets[i].end - targets[i].start);
        }
        const std::size_t loci = offsets.back();
        if (loci == 0) throw std::runtime_error("BAD_INPUT: selected intervals have zero length");

        // GATK excludes reference-N loci by default.  Build one flat mask in
        // Host memory so both the observation projection and every summary
        // denominator use exactly the same eligibility rule.  The mask is
        // deliberately a byte vector: it is cheap to retain for a batch and
        // avoids copying reference strings into the Kokkos counting kernel.
        reference_fai = fai_load(options.reference.c_str());
        if (!reference_fai)
            throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
        std::vector<std::uint8_t> canonical_locus(loci, 1);
        if (!options.include_ref_n_sites) {
            for (std::size_t target_index = 0; target_index < targets.size(); ++target_index) {
                const auto& target = targets[target_index];
                hts_pos_t fetched_length = 0;
                char* sequence = faidx_fetch_seq64(
                    reference_fai, reader.header().contigs[static_cast<std::size_t>(target.tid)].c_str(),
                    target.start, target.end - 1, &fetched_length);
                if (!sequence || fetched_length != target.end - target.start) {
                    free(sequence);
                    throw std::runtime_error("BAD_INPUT: failed to fetch reference bases for DepthOfCoverage");
                }
                for (std::int64_t offset = 0; offset < target.end - target.start; ++offset) {
                    const auto base = static_cast<unsigned char>(sequence[offset]);
                    const auto upper = static_cast<unsigned char>(std::toupper(base));
                    canonical_locus[offsets[target_index] + static_cast<std::size_t>(offset)] =
                        (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') ? 1 : 0;
                }
                free(sequence);
            }
        }
        const std::size_t counted_loci = options.include_ref_n_sites
            ? loci
            : static_cast<std::size_t>(std::count(canonical_locus.begin(), canonical_locus.end(), std::uint8_t{1}));
        if (counted_loci == 0)
            throw std::runtime_error("BAD_INPUT: no canonical reference loci remain; use --include-ref-n-sites to include N sites");

        // LayoutRight keeps the locus dimension contiguous for each sample;
        // this is the access pattern used by the Host output pass and avoids
        // a second transpose when a multi-sample BAM is partitioned.
        Kokkos::View<std::uint64_t**, Kokkos::LayoutRight> depth("doc_depth", sample_count, loci);
        Kokkos::View<std::uint64_t***, Kokkos::LayoutRight> bases("doc_bases_count", sample_count, loci, 6);
        Kokkos::deep_copy(depth, static_cast<std::uint64_t>(0));
        Kokkos::deep_copy(bases, static_cast<std::uint64_t>(0));
        std::unordered_set<FragmentObservationKey, FragmentObservationKeyHash> seen_fragments;
        if (options.count_fragments) seen_fragments.reserve(1024);
        std::uint64_t reads_seen = 0, reads_used = 0, reads_filtered = 0, observations = 0;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        std::uint64_t count_kernel_batches = 0;
        std::uint64_t count_kernel_observations = 0;
        double count_kernel_prepare_seconds = 0.0;
        double count_kernel_execute_seconds = 0.0;
        const auto read_group_samples = reader.header().read_group_samples;
        const bool count_deletions = options.include_deletions && !options.ignore_deletion_sites;
        const bool header_has_no_samples = reader.header().samples.empty();
        std::uint64_t next_record_ordinal = 0;
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<fastgatk::io::ReadBatch, CoverageBatchResult, CoverageBatchResult>;
        const auto safe_budget = resources.safe_memory_budget_bytes();
        const auto decoded_capacity = safe_budget > 0
            ? std::max<std::uint64_t>(safe_budget / 4, 1)
            : std::max<std::uint64_t>(64ULL * 1024ULL * 1024ULL,
                static_cast<std::uint64_t>(options.batch_records) * 512ULL);
        const auto result_capacity = std::max<std::uint64_t>(sizeof(CoverageBatchResult) * 8ULL, 4096ULL);
        Pipeline pipeline(
            Pipeline::Limits{decoded_capacity, result_capacity, result_capacity},
            [&reader]() -> std::optional<fastgatk::io::ReadBatch> {
                fastgatk::io::ReadBatch batch;
                if (!reader.next(batch)) return std::nullopt;
                return batch;
            },
            [&](fastgatk::io::ReadBatch batch) -> std::optional<CoverageBatchResult> {
                CoverageBatchResult result;
                result.reads_seen = batch.records();
                const auto batch_start_ordinal = next_record_ordinal;
                next_record_ordinal += batch.records();
                std::vector<std::uint64_t> observation_loci;
                std::vector<std::uint8_t> observation_bases;
                std::vector<std::uint32_t> observation_samples;
                for (std::size_t record = 0; record < batch.records(); ++record) {
                    const std::uint16_t flags = batch.flags.size() == batch.records() ? batch.flags[record] : 0;
                    // DepthOfCoverage's GATK 4.6.2.0 default filter set is
                    // WellformedReadFilter + NotDuplicateReadFilter +
                    // NotSecondaryAlignmentReadFilter + MappedReadFilter.
                    // QC-fail (0x200) and supplementary (0x800) reads are
                    // intentionally *not* part of that set and must remain
                    // countable.  The former used to be dropped here,
                    // producing a silent depth mismatch on QC-fail BAMs.
                    if ((flags & 0x4U) || (flags & 0x100U) || (flags & 0x400U)) {
                        ++result.reads_filtered;
                        continue;
                    }
                    std::string read_sample;
                    if (batch.read_group_offsets.size() == batch.records() + 1) {
                        const auto rg_begin = batch.read_group_offsets[record];
                        const auto rg_end = batch.read_group_offsets[record + 1];
                        if (rg_begin <= rg_end && rg_end <= batch.read_groups.size())
                            read_sample.assign(batch.read_groups.begin() + rg_begin,
                                               batch.read_groups.begin() + rg_end);
                    }
                    if (!read_sample.empty()) {
                        const auto rg = read_group_samples.find(read_sample);
                        read_sample = rg == read_group_samples.end() ? std::string{} : rg->second;
                    }
                    std::size_t sample_slot = std::numeric_limits<std::size_t>::max();
                    if (!read_sample.empty()) {
                        const auto selected = sample_index.find(read_sample);
                        if (selected != sample_index.end()) sample_slot = selected->second;
                    } else if (sample_count == 1 && header_has_no_samples) {
                        sample_slot = 0;
                    }
                    if (sample_slot == std::numeric_limits<std::size_t>::max()) continue;
                    if (batch.has_cigar() && !batch.cigar_layout_valid())
                        throw std::runtime_error("BAD_INPUT: malformed CIGAR metadata");
                    // WellformedReadFilter rejects any CIGAR containing the
                    // reference-skip (N) operator.  It is not equivalent to
                    // silently advancing the reference cursor: doing so
                    // would count the two M islands of a spliced read even
                    // though GATK drops the complete read.  Keep this check
                    // at the record boundary so other reads in the batch are
                    // still processed.
                    if (batch.has_cigar()) {
                        bool contains_reference_skip = false;
                        for (std::size_t cigar = batch.cigar_offsets[record];
                             cigar < batch.cigar_offsets[record + 1]; ++cigar) {
                            if (fastgatk::io::CigarOp::unpack(batch.cigar_ops[cigar]).code ==
                                fastgatk::io::CigarOpCode::ReferenceSkip) {
                                contains_reference_skip = true;
                                break;
                            }
                        }
                        if (contains_reference_skip) {
                            ++result.reads_filtered;
                            continue;
                        }
                    }
                    const std::size_t begin = batch.offsets[record], end = batch.offsets[record + 1];
                    bool used = false;
                    std::int64_t reference_position = batch.positions[record];
                    std::size_t read_cursor = 0;
                    std::string fragment_name;
                    if (options.count_fragments && batch.name_offsets.size() == batch.records() + 1) {
                        const auto name_begin = batch.name_offsets[record];
                        const auto name_end = batch.name_offsets[record + 1];
                        if (name_begin <= name_end && name_end <= batch.names.size())
                            fragment_name.assign(batch.names.begin() + name_begin, batch.names.begin() + name_end);
                    }
                    if (fragment_name.empty())
                        fragment_name = "__record_" + std::to_string(batch_start_ordinal + record);
                    const auto emit = [&](std::int64_t position, std::size_t read_offset, std::uint8_t base) {
                        const auto locus = find_locus(targets, offsets, batch.tids[record], position);
                        if (locus == std::numeric_limits<std::size_t>::max()) return;
                        if (!options.include_ref_n_sites && canonical_locus[locus] == 0) return;
                        if (base != 5) {
                            if (read_offset >= end - begin) return;
                            const auto quality = batch.qualities[begin + read_offset];
                            if (quality < options.min_base_quality || quality > options.max_base_quality) return;
                        }
                        if (options.count_fragments &&
                            !seen_fragments.emplace(FragmentObservationKey{
                                static_cast<std::uint32_t>(sample_slot), static_cast<std::uint64_t>(locus), fragment_name
                            }).second)
                            return;
                        observation_loci.push_back(static_cast<std::uint64_t>(locus));
                        observation_bases.push_back(base);
                        observation_samples.push_back(static_cast<std::uint32_t>(sample_slot));
                        used = true;
                    };
                    if (!batch.has_cigar()) {
                        for (; read_cursor < end - begin; ++read_cursor) {
                            const auto base = static_cast<char>(batch.bases[begin + read_cursor]);
                            const char* names = "ACGTN";
                            int index = 4;
                            for (int i = 0; i < 4; ++i)
                                if (std::toupper(static_cast<unsigned char>(base)) == names[i]) index = i;
                            emit(reference_position + static_cast<std::int64_t>(read_cursor), read_cursor,
                                 static_cast<std::uint8_t>(index));
                        }
                    } else {
                        for (std::size_t cigar = batch.cigar_offsets[record];
                             cigar < batch.cigar_offsets[record + 1]; ++cigar) {
                            const auto op = fastgatk::io::CigarOp::unpack(batch.cigar_ops[cigar]);
                            if (op.projects_base()) {
                                for (std::uint32_t j = 0; j < op.length; ++j) {
                                    const auto base = static_cast<char>(batch.bases[begin + read_cursor + j]);
                                    const char* names = "ACGTN";
                                    int index = 4;
                                    for (int i = 0; i < 4; ++i)
                                        if (std::toupper(static_cast<unsigned char>(base)) == names[i]) index = i;
                                    emit(reference_position + static_cast<std::int64_t>(j), read_cursor + j,
                                         static_cast<std::uint8_t>(index));
                                }
                            } else if (op.code == fastgatk::io::CigarOpCode::Deletion && count_deletions) {
                                for (std::uint32_t j = 0; j < op.length; ++j)
                                    emit(reference_position + static_cast<std::int64_t>(j), 0, 5);
                            }
                            if (op.consumes_read()) read_cursor += op.length;
                            if (op.consumes_reference()) reference_position += op.length;
                        }
                    }
                    if (used) ++result.reads_used;
                }
                result.observations = observation_loci.size();
                if (observation_loci.empty()) return result;
                fastgatk::core::HostBatch count_host("depth-of-coverage-count-v1");
                count_host.records = observation_loci.size();
                count_host.bytes = observation_loci.size() * sizeof(std::uint64_t) +
                                   observation_bases.size() * sizeof(std::uint8_t) +
                                   observation_samples.size() * sizeof(std::uint32_t);
                fastgatk::core::KernelPlan<ExecSpace> count_plan("depth-of-coverage-count");
                count_plan.begin_prepare(count_host);
                Kokkos::View<std::uint64_t*> loci_view("doc_loci", observation_loci.size());
                Kokkos::View<std::uint8_t*> base_view("doc_bases", observation_bases.size());
                Kokkos::View<std::uint32_t*> sample_view("doc_samples", observation_samples.size());
                auto host_loci = Kokkos::create_mirror_view(loci_view);
                auto host_base = Kokkos::create_mirror_view(base_view);
                auto host_sample = Kokkos::create_mirror_view(sample_view);
                for (std::size_t i = 0; i < observation_loci.size(); ++i) {
                    host_loci(i) = observation_loci[i];
                    host_base(i) = observation_bases[i];
                    host_sample(i) = observation_samples[i];
                }
                Kokkos::deep_copy(loci_view, host_loci);
                Kokkos::deep_copy(base_view, host_base);
                Kokkos::deep_copy(sample_view, host_sample);
                fastgatk::core::DeviceBatch<ExecSpace> count_device(observation_loci.size());
                count_device.bind("loci", loci_view);
                count_device.bind("bases", base_view);
                count_device.bind("samples", sample_view);
                count_device.bind("depth", depth);
                count_device.bind("base_counts", bases);
                ExecSpace().fence();
                count_plan.end_prepare(count_device);
                count_plan.begin_execute();
                Kokkos::parallel_for("depth_of_coverage_count",
                    Kokkos::RangePolicy<ExecSpace>(0, observation_loci.size()),
                    KOKKOS_LAMBDA(const std::size_t i) {
                        const auto locus = loci_view(i);
                        const auto base = base_view(i);
                        const auto sample = sample_view(i);
                        Kokkos::atomic_add(&depth(sample, locus), static_cast<std::uint64_t>(1));
                        Kokkos::atomic_add(&bases(sample, locus, base), static_cast<std::uint64_t>(1));
                    });
                ExecSpace().fence();
                count_plan.end_execute();
                result.kernel_batches = 1;
                result.kernel_observations = observation_loci.size();
                result.kernel_prepare_seconds = count_plan.telemetry().prepare_seconds;
                result.kernel_execute_seconds = count_plan.telemetry().execute_seconds;
                return result;
            },
            [](CoverageBatchResult result) -> std::optional<CoverageBatchResult> { return result; },
            [&](CoverageBatchResult result) {
                reads_seen += result.reads_seen;
                reads_used += result.reads_used;
                reads_filtered += result.reads_filtered;
                observations += result.observations;
                count_kernel_batches += result.kernel_batches;
                count_kernel_observations += result.kernel_observations;
                count_kernel_prepare_seconds += result.kernel_prepare_seconds;
                count_kernel_execute_seconds += result.kernel_execute_seconds;
            },
            [](const fastgatk::io::ReadBatch& batch) {
                return std::max<std::uint64_t>(1, batch.bytes());
            },
            [](const CoverageBatchResult&) { return static_cast<std::uint64_t>(sizeof(CoverageBatchResult)); },
            [](const CoverageBatchResult&) { return static_cast<std::uint64_t>(sizeof(CoverageBatchResult)); });
        const auto pipeline_metrics = pipeline.run();
        std::vector<std::uint64_t> host_depth(sample_count * loci, 0);
        std::vector<std::uint64_t> host_bases(sample_count * loci * 6, 0);
        auto depth_host = Kokkos::create_mirror_view(depth);
        auto bases_host = Kokkos::create_mirror_view(bases);
        Kokkos::deep_copy(depth_host, depth);
        Kokkos::deep_copy(bases_host, bases);
        for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
            for (std::size_t i = 0; i < loci; ++i) {
                host_depth[sample_slot * loci + i] = depth_host(sample_slot, i);
                for (int j = 0; j < 6; ++j)
                    host_bases[(sample_slot * loci + i) * 6 + j] = bases_host(sample_slot, i, j);
            }
        }

        std::vector<std::string> outputs;
        const auto add_output = [&](const std::string& path) { outputs.push_back(path); };
        if (!options.omit_locus_table) {
            std::ofstream out(options.output);
            if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write " + options.output);
            out << "Locus,Total_Depth";
            if (sample_count == 1)
                out << ",Average_Depth_sample,Depth_for_" << samples.front();
            else
                for (const auto& sample : samples) out << ",Average_Depth_" << sample << ",Depth_for_" << sample;
            if (options.print_base_counts)
                for (const auto& sample : samples) out << ',' << sample << "_base_counts";
            out << '\n';
            for (std::size_t target = 0; target < targets.size(); ++target) {
                for (std::int64_t p = targets[target].start; p < targets[target].end; ++p) {
                    const auto locus = offsets[target] + static_cast<std::size_t>(p - targets[target].start);
                    if (!options.include_ref_n_sites && canonical_locus[locus] == 0) continue;
                    std::uint64_t total_locus_depth = 0;
                    for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot)
                        total_locus_depth += host_depth[sample_slot * loci + locus];
                    out << reader.header().contigs[static_cast<std::size_t>(targets[target].tid)] << ':' << (p + 1) << ','
                        << total_locus_depth;
                    for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
                        const auto value = host_depth[sample_slot * loci + locus];
                        out << ',' << format_two(static_cast<double>(value)) << ',' << value;
                    }
                    if (options.print_base_counts)
                        for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot)
                            out << ',' << base_counts(host_bases.data() + (sample_slot * loci + locus) * 6,
                                                       count_deletions);
                    out << '\n';
                }
            }
            add_output(options.output);
        }
        std::vector<std::uint64_t> sample_total(sample_count, 0), sample_above(sample_count, 0);
        std::vector<double> sample_mean(sample_count, 0.0);
        std::vector<int> sample_q1(sample_count, 0), sample_median(sample_count, 0), sample_q3(sample_count, 0);
        const auto quantile_bin = [&](std::size_t sample_slot, double fraction) {
            if (counted_loci == 0) return 0;
            std::vector<std::uint64_t> sorted;
            sorted.reserve(counted_loci);
            for (std::size_t locus = 0; locus < loci; ++locus)
                if (options.include_ref_n_sites || canonical_locus[locus] != 0)
                    sorted.push_back(host_depth[sample_slot * loci + locus]);
            std::sort(sorted.begin(), sorted.end());
            const auto rank = static_cast<std::size_t>(fraction * static_cast<double>(sorted.size() - 1));
            const auto value = sorted[rank];
            return value >= static_cast<std::uint64_t>(options.stop) ? options.stop + 1 : static_cast<int>(value + 1);
        };
        for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
            for (std::size_t locus = 0; locus < loci; ++locus) {
                if (!options.include_ref_n_sites && canonical_locus[locus] == 0) continue;
                const auto value = host_depth[sample_slot * loci + locus];
                sample_total[sample_slot] += value;
                if (value >= static_cast<std::uint64_t>(options.summary_coverage_threshold)) ++sample_above[sample_slot];
            }
            sample_mean[sample_slot] = counted_loci == 0 ? 0.0 : static_cast<double>(sample_total[sample_slot]) / counted_loci;
            sample_q1[sample_slot] = quantile_bin(sample_slot, 0.25);
            sample_median[sample_slot] = quantile_bin(sample_slot, 0.5);
            sample_q3[sample_slot] = quantile_bin(sample_slot, 0.75);
        }
        std::uint64_t total_depth = std::accumulate(sample_total.begin(), sample_total.end(), std::uint64_t{0});
        const double mean = sample_count == 0 || counted_loci == 0 ? 0.0 : static_cast<double>(total_depth) / (sample_count * counted_loci);
        const auto interval_summary = options.output + ".sample_interval_summary";
        if (!options.omit_interval_statistics) {
          std::ofstream out(interval_summary); if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write interval summary");
          out << "Target,total_coverage,average_coverage";
          for (const auto& sample : samples)
              out << ',' << sample << "_total_cvg," << sample << "_mean_cvg," << sample << "_granular_Q1,"
                  << sample << "_granular_median," << sample << "_granular_Q3," << sample << "_%_above_"
                  << options.summary_coverage_threshold;
          out << '\n';
          for (std::size_t i = 0; i < targets.size(); ++i) {
              out << reader.header().contigs[static_cast<std::size_t>(targets[i].tid)] << ':' << targets[i].start + 1 << '-' << targets[i].end
                  << ',' << total_depth << ',' << format_two(mean);
              for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot)
                    out << ',' << sample_total[sample_slot] << ',' << format_two(sample_mean[sample_slot])
                      << ',' << sample_q1[sample_slot] << ',' << sample_median[sample_slot] << ',' << sample_q3[sample_slot]
                      << ',' << format_one(counted_loci == 0 ? 0.0 : 100.0 * static_cast<double>(sample_above[sample_slot]) / counted_loci);
              out << '\n';
          }
          add_output(interval_summary);
        }
        const auto interval_statistics = options.output + ".sample_interval_statistics";
        if (!options.omit_interval_statistics) {
            std::ofstream out(interval_statistics);
            if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write interval statistics");
            out << "Number_of_sources";
            for (int i = 0; i <= options.stop; ++i) out << ",depth>=" << i;
            out << '\n' << "At_least_" << sample_count << "_samples";
            for (int i = 0; i <= options.stop; ++i) {
                std::size_t count = 0;
                for (const auto value : sample_mean) if (value >= static_cast<double>(i)) ++count;
                out << ',' << count;
            }
            out << '\n';
            add_output(interval_statistics);
        }
        const auto sample_summary = options.output + ".sample_summary";
        if (!options.omit_per_sample_statistics) {
          std::ofstream out(sample_summary); if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sample summary");
          out << "sample_id,total,mean,granular_third_quartile,granular_median,granular_first_quartile,%_bases_above_" << options.summary_coverage_threshold << '\n';
          for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot)
              out << samples[sample_slot] << ',' << sample_total[sample_slot] << ',' << format_two(sample_mean[sample_slot])
                  << ',' << sample_q3[sample_slot] << ',' << sample_median[sample_slot] << ',' << sample_q1[sample_slot]
                  << ',' << format_one(counted_loci == 0 ? 0.0 : 100.0 * static_cast<double>(sample_above[sample_slot]) / counted_loci) << '\n';
          out << "Total," << total_depth << ',' << format_two(mean) << ",N/A,N/A,N/A,N/A\n";
          add_output(sample_summary);
        }
        const auto sample_statistics = options.output + ".sample_statistics";
        if (!options.omit_per_sample_statistics) {
            std::ofstream out(sample_statistics);
            if (!out) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sample statistics");
            out << "Source_of_reads";
            for (int i = 0; i < options.stop; ++i) out << ",from_" << i << "_to_" << (i + 1) << ')';
            out << ",from_" << options.stop << "_to_inf\n";
            for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
                out << "sample_" << samples[sample_slot];
                for (int i = 0; i < options.stop; ++i) {
                    std::uint64_t count = 0;
                    for (std::size_t locus = 0; locus < loci; ++locus)
                        if ((options.include_ref_n_sites || canonical_locus[locus] != 0) &&
                            host_depth[sample_slot * loci + locus] == static_cast<std::uint64_t>(i)) ++count;
                    out << ',' << count;
                }
                std::uint64_t tail = 0;
                for (std::size_t locus = 0; locus < loci; ++locus)
                    if ((options.include_ref_n_sites || canonical_locus[locus] != 0) &&
                        host_depth[sample_slot * loci + locus] >= static_cast<std::uint64_t>(options.stop)) ++tail;
                out << ',' << tail << '\n';
            }
            add_output(sample_statistics);
        }
        const auto cumulative_counts = options.output + ".sample_cumulative_coverage_counts";
        const auto cumulative_proportions = options.output + ".sample_cumulative_coverage_proportions";
        {
            std::ofstream counts(cumulative_counts), proportions(cumulative_proportions);
            if (!counts || !proportions) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write cumulative coverage outputs");
            counts << "sample"; proportions << "sample";
            for (int i = 0; i <= options.stop; ++i) { counts << ",gte_" << i; proportions << ",gte_" << i; }
            counts << '\n'; proportions << '\n';
            for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
                counts << (sample_count == 1 ? "NSamples_1" : samples[sample_slot]);
                proportions << samples[sample_slot];
                for (int i = 0; i <= options.stop; ++i) {
                    std::uint64_t at_least = 0;
                    for (std::size_t locus = 0; locus < loci; ++locus)
                        if ((options.include_ref_n_sites || canonical_locus[locus] != 0) &&
                            host_depth[sample_slot * loci + locus] >= static_cast<std::uint64_t>(i)) ++at_least;
                    counts << ',' << at_least;
                    proportions << ',' << format_two(counted_loci == 0 ? 0.0 : static_cast<double>(at_least) / counted_loci);
                }
                counts << '\n'; proportions << '\n';
            }
        }
        add_output(cumulative_counts); add_output(cumulative_proportions);
        for (const auto& path : outputs) if (file_bytes(path) == 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: output is missing or empty: " + path);

        const auto wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_begin).count();
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path); if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest");
        manifest << "{\"schema_version\":1,\"tool\":\"DepthOfCoverage\",\"implementation\":\"fastgatk-depth-of-coverage\",\"status\":\"prototype\",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\",\"determinism\":\"strict\",\"primary_output\":\"" << json_escape(options.output) << "\",\"samples\":[";
        for (std::size_t sample_slot = 0; sample_slot < sample_count; ++sample_slot) {
            if (sample_slot) manifest << ',';
            manifest << '\"' << json_escape(samples[sample_slot]) << '\"';
        }
        manifest << "],\"compatibility\":{\"partition_type\":\"sample\",\"count_type\":\"" << (options.count_fragments ? "COUNT_FRAGMENTS" : "COUNT_READS") << "\",\"bit_identical_to_gatk\":false,\"include_deletions\":" << (options.include_deletions ? "true" : "false") << ",\"ignore_deletion_sites\":" << (options.ignore_deletion_sites ? "true" : "false") << ",\"effective_include_deletions\":" << (count_deletions ? "true" : "false") << ",\"multi_sample_read_group_partition\":"
                 << (sample_count > 1 ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection ? "INTERSECTION" : "UNION")
                 << "\",\"interval_padding\":" << options.interval_padding
                 << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                 << ",\"excluded_intervals\":" << reader.exclusion_intervals().size()
                 << "},\"outputs\":[";
        for (std::size_t i = 0; i < outputs.size(); ++i) { if (i) manifest << ','; manifest << "{\"path\":\"" << json_escape(outputs[i]) << "\",\"complete\":true,\"bytes\":" << file_bytes(outputs[i]) << "}"; }
        manifest << "],\"telemetry\":{\"resources\":{\"backend\":\"" << json_escape(reader.backend_description()) << "\"},\"loci\":" << loci << ",\"counted_loci\":" << counted_loci << ",\"include_ref_n_sites\":" << (options.include_ref_n_sites ? "true" : "false") << ",\"include_deletions\":" << (options.include_deletions ? "true" : "false") << ",\"ignore_deletion_sites\":" << (options.ignore_deletion_sites ? "true" : "false") << ",\"effective_include_deletions\":" << (count_deletions ? "true" : "false") << ",\"count_type\":\"" << (options.count_fragments ? "COUNT_FRAGMENTS" : "COUNT_READS") << "\",\"fragment_keys\":" << seen_fragments.size() << ",\"targets\":" << targets.size() << ",\"requested_interval_files\":" << reader.interval_file_inputs() << ",\"requested_interval_records\":" << reader.interval_file_records() << ",\"excluded_intervals\":" << reader.exclusion_intervals().size() << ",\"interval_set_rule\":\"" << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection ? "INTERSECTION" : "UNION") << "\",\"interval_padding\":" << options.interval_padding << ",\"interval_exclusion_padding\":" << options.exclusion_padding << ",\"reads_seen\":" << reads_seen << ",\"reads_used\":" << reads_used << ",\"reads_filtered\":" << reads_filtered << ",\"observations\":" << observations << ",\"initial_batch_records\":" << options.batch_records << ",\"effective_batch_records\":" << reader.batch_records() << ",\"count_kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"count_kernel_execution_space\":\"" << ExecSpace::name() << "\",\"count_kernel_execution_policy\":\"RangePolicy\",\"count_kernel_batches\":" << count_kernel_batches << ",\"count_kernel_observations\":" << count_kernel_observations << ",\"count_kernel_prepare_seconds\":" << count_kernel_prepare_seconds << ",\"count_kernel_execute_seconds\":" << count_kernel_execute_seconds << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items << ",\"pipeline_decoded_bytes\":" << pipeline_metrics.decoded_bytes << ",\"pipeline_computed_bytes\":" << pipeline_metrics.computed_bytes << ",\"pipeline_encoded_bytes\":" << pipeline_metrics.encoded_bytes << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes << ",\"output_bytes\":" << file_bytes(options.output) << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}}\n";
        std::cout << "{\"tool\":\"DepthOfCoverage\",\"status\":\"prototype\",\"sample_count\":" << sample_count
                  << ",\"loci\":" << loci << ",\"counted_loci\":" << counted_loci
                  << ",\"include_ref_n_sites\":" << (options.include_ref_n_sites ? "true" : "false")
                  << ",\"count_type\":\"" << (options.count_fragments ? "COUNT_FRAGMENTS" : "COUNT_READS") << "\""
                  << ",\"reads_seen\":" << reads_seen << ",\"reads_used\":" << reads_used
                  << ",\"observations\":" << observations << ",\"count_kernel_batches\":" << count_kernel_batches
                  << ",\"count_kernel_observations\":" << count_kernel_observations
                  << ",\"count_kernel_execution_space\":\"" << ExecSpace::name() << "\",\"output_bytes\":" << file_bytes(options.output)
                  << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}\n";
        // Kokkos allocations must be released before finalize; otherwise the
        // runtime diagnoses a post-finalize View destructor as a hard error.
        depth_host = decltype(depth_host)();
        bases_host = decltype(bases_host)();
        depth = decltype(depth)();
        bases = decltype(bases)();
        fai_destroy(reference_fai);
        reference_fai = nullptr;
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        if (reference_fai != nullptr) fai_destroy(reference_fai);
        if (initialized) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
#endif
}
