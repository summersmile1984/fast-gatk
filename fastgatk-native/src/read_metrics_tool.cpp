#include "fastgatk/core/plan.hpp"
#include "fastgatk/runtime/pipeline.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include "fastgatk/io/hts_read_guard.hpp"

#include <array>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cmath>
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
#include <string_view>
#include <vector>

namespace {

#ifdef FASTGATK_FLAGSTAT
constexpr bool kFlagStat = true;
constexpr const char* kToolName = "FlagStat";
constexpr const char* kImplementation = "fastgatk-flag-stat";
#else
constexpr bool kFlagStat = false;
constexpr const char* kToolName = "CountReads";
constexpr const char* kImplementation = "fastgatk-count-reads";
#endif

struct Options {
    std::vector<std::string> inputs;
    std::string reference;
    std::string output;
    std::string manifest;
    bool quiet = false;
    std::string tmp_dir;
    std::string verbosity;
    std::vector<std::string> intervals;
    std::vector<std::string> exclude_intervals;
    enum class IntervalSetRule : std::uint8_t { Union, Intersection };
    IntervalSetRule interval_set_rule = IntervalSetRule::Union;
    enum class IntervalMergingRule : std::uint8_t { All, OverlappingOnly };
    IntervalMergingRule interval_merging_rule = IntervalMergingRule::All;
    std::int64_t interval_padding = 0;
    std::int64_t interval_exclusion_padding = 0;
    std::vector<std::string> read_filters;
    std::vector<std::string> inverted_read_filters;
    std::vector<std::string> disabled_read_filters;
    std::vector<std::string> read_names;
    std::vector<std::string> keep_read_groups;
    std::vector<std::pair<std::string, std::string>> read_group_blacklist;
    enum class ReadTagOperator : std::uint8_t {
        Less, LessOrEqual, Greater, GreaterOrEqual, Equal, NotEqual
    };
    std::string read_filter_tag;
    double read_filter_tag_comp = 0.0;
    ReadTagOperator read_filter_tag_op = ReadTagOperator::Equal;
    bool read_filter_tag_comp_seen = false;
    bool read_filter_tag_op_seen = false;
    bool disable_tool_default_read_filters = false;
    bool require_read_length = false;
    std::uint32_t min_read_length = 1;
    std::uint32_t max_read_length = std::numeric_limits<std::uint32_t>::max();
    bool mapping_quality_argument_seen = false;
    int min_mapping_quality = 10;
    std::optional<int> max_mapping_quality;
    bool fragment_length_argument_seen = false;
    int min_fragment_length = 0;
    int max_fragment_length = 1000000;
    std::size_t batch_records = 65536;
    int threads = 1;
};

struct Interval {
    int tid = -1;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

struct ReadProjection {
    std::uint16_t flags = 0;
    std::int32_t tid = -1;
    std::int32_t mtid = -1;
    std::int64_t pos = -1;
    std::int64_t mpos = -1;
    std::int32_t fragment_length = 0;
    std::uint8_t mapq = 0;
};

struct KernelTotals {
    std::uint64_t records = 0;
    std::uint64_t batches = 0;
    std::size_t prepare_bytes = 0;
    std::size_t device_bytes = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::array<std::uint64_t, 12> counters{};
    std::size_t buffer_capacity_records = 0;
    std::size_t buffer_allocations = 0;
    std::size_t buffer_reuses = 0;
    std::uint64_t pipeline_decoded_items = 0;
    std::uint64_t pipeline_computed_items = 0;
    std::uint64_t pipeline_encoded_items = 0;
    std::uint64_t pipeline_decoded_bytes = 0;
    std::uint64_t pipeline_computed_bytes = 0;
    std::uint64_t pipeline_encoded_bytes = 0;
    std::uint64_t pipeline_peak_decoded_bytes = 0;
    std::uint64_t pipeline_peak_computed_bytes = 0;
    std::uint64_t pipeline_peak_encoded_bytes = 0;
};

// Read batches arrive in a bounded stream.  Keep the Kokkos staging views
// alive across batches and grow them only when a larger batch is observed;
// this turns allocator churn into a one-time capacity change while preserving
// the exact logical batch length in DeviceBatch telemetry.
template<class ExecSpace>
class ReadMetricsBuffers {
public:
    using FlagsView = Kokkos::View<std::uint16_t*>;
    FlagsView flags;
#ifdef FASTGATK_FLAGSTAT
    Kokkos::View<std::int32_t*> tids;
    Kokkos::View<std::int32_t*> mtids;
    Kokkos::View<std::int64_t*> positions;
    Kokkos::View<std::int64_t*> mate_positions;
    Kokkos::View<std::uint8_t*> mapqs;
    Kokkos::View<std::uint64_t*> counters;
#endif
    std::size_t capacity = 0;
    std::size_t allocations = 0;
    std::size_t reuses = 0;

    void ensure(std::size_t records) {
        if (records <= capacity) {
            ++reuses;
            return;
        }
        const auto next = std::max(records, capacity == 0 ? std::size_t{1} : capacity * 2);
        flags = FlagsView("read_flags_persistent", next);
#ifdef FASTGATK_FLAGSTAT
        tids = Kokkos::View<std::int32_t*>("read_tids_persistent", next);
        mtids = Kokkos::View<std::int32_t*>("read_mtids_persistent", next);
        positions = Kokkos::View<std::int64_t*>("read_positions_persistent", next);
        mate_positions = Kokkos::View<std::int64_t*>("read_mate_positions_persistent", next);
        mapqs = Kokkos::View<std::uint8_t*>("read_mapqs_persistent", next);
        counters = Kokkos::View<std::uint64_t*>("flagstat_counters_persistent", 12);
#endif
        capacity = next;
        ++allocations;
    }

    void clear() {
        flags = FlagsView{};
#ifdef FASTGATK_FLAGSTAT
        tids = Kokkos::View<std::int32_t*>{};
        mtids = Kokkos::View<std::int32_t*>{};
        positions = Kokkos::View<std::int64_t*>{};
        mate_positions = Kokkos::View<std::int64_t*>{};
        mapqs = Kokkos::View<std::uint8_t*>{};
        counters = Kokkos::View<std::uint64_t*>{};
#endif
        capacity = 0;
    }
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string value_for(int& index, int argc, char** argv, const std::string& argument,
                     std::string_view name, const char* short_name = nullptr) {
    if (const auto value = inline_value(argument, name); !value.empty()) return value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument("missing value for " + std::string(name));
}

bool is_supported_read_filter(const std::string& name);

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-" << (kFlagStat ? "flag-stat" : "count-reads")
                      << " (GATK-compatible native; contract-compatible)\n"
                         "  -I, --input FILE             input SAM/BAM/CRAM (repeatable)\n"
                         "  -R, --reference FILE        optional reference FASTA (required for some CRAMs)\n"
                         "  -O, --output FILE            optional text output\n"
                         "  -L, --intervals REGION       repeatable contig[:start-end] or interval file\n"
                         "  -isr, --interval-set-rule R UNION (default) or INTERSECTION across -L selectors\n"
                         "  -imr, --interval-merging-rule R ALL (default) or OVERLAPPING_ONLY\n"
                         "  -ip, --interval-padding N    padding (bp) for included intervals\n"
                         "  -ixp, --interval-exclusion-padding N  padding (bp) for excluded intervals\n"
                         "  -XL, --exclude-intervals REGION  repeatable excluded region/file\n"
                         "  -RF, --read-filter NAME      repeatable common read filter\n"
                         "  -XRF, --inverted-read-filter NAME  retain reads rejected by this filter\n"
                         "      supported: MappingQuality, fragment-length, mate/strand and CIGAR/RG filters\n"
                         "  -DF, --disable-read-filter NAME  disable a default read filter\n"
                         "      --disable-tool-default-read-filters  remove implicit tool filters\n"
                         "      --min-read-length N      ReadLengthReadFilter lower bound\n"
                         "      --max-read-length N      ReadLengthReadFilter upper bound\n"
                         "      --minimum-mapping-quality N  MappingQualityReadFilter lower bound (inclusive)\n"
                         "      --maximum-mapping-quality N  MappingQualityReadFilter upper bound (inclusive)\n"
                         "      --min-fragment-length N    FragmentLengthReadFilter lower bound\n"
                         "      --max-fragment-length N    FragmentLengthReadFilter upper bound\n"
                         "      --read-name NAME         ReadNameReadFilter exact name (repeatable)\n"
                         "      --keep-read-group NAME  ReadGroupReadFilter exact RG (repeatable)\n"
                         "      --read-group-black-list ATTR:VALUE  ReadGroupBlackListReadFilter (repeatable)\n"
                         "      --read-filter-tag TAG   ReadTagValueFilter two-character tag\n"
                         "      --read-filter-tag-comp N  numeric comparison value\n"
                         "      --read-filter-tag-op OP  LESS/LESS_OR_EQUAL/GREATER/GREATER_OR_EQUAL/EQUAL/NOT_EQUAL\n"
                         "      --batch-records N        bounded HTSlib/Kokkos batch size (default 65536)\n"
                         "      --threads N              HTSlib/Kokkos threads\n"
                         "      --output-manifest FILE   OutputManifest JSON\n"
                         "      --QUIET[=true|false]     Picard/GATK quiet switch\n"
                         "      --tmp-dir DIR            temporary directory\n"
                         "      --verbosity LEVEL        ERROR/WARNING/INFO/DEBUG\n";
            std::exit(0);
        } else if (argument == "-I" || argument == "--input" ||
                   !inline_value(argument, "--input").empty()) {
            options.inputs.push_back(value_for(index, argc, argv, argument, "--input", "-I"));
        } else if (argument == "-R" || argument == "--reference" ||
                   !inline_value(argument, "--reference").empty()) {
            options.reference = value_for(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || argument == "--intervals" || argument == "--interval" ||
                   argument == "--region" || !inline_value(argument, "--intervals").empty() ||
                   !inline_value(argument, "--interval").empty() ||
                   !inline_value(argument, "--region").empty()) {
            const auto name = argument.rfind("--interval", 0) == 0 ? "--intervals" :
                (argument.rfind("--region", 0) == 0 ? "--region" : "--intervals");
            options.intervals.push_back(value_for(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-isr" || argument == "--interval-set-rule" ||
                   !inline_value(argument, "--interval-set-rule").empty()) {
            auto value = value_for(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
                return static_cast<char>(std::toupper(character));
            });
            if (value == "UNION") options.interval_set_rule = Options::IntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = Options::IntervalSetRule::Intersection;
            else throw std::invalid_argument("BAD_INPUT: --interval-set-rule must be UNION or INTERSECTION");
        } else if (argument == "-imr" || argument == "--interval-merging-rule" ||
                   !inline_value(argument, "--interval-merging-rule").empty()) {
            auto value = value_for(index, argc, argv, argument, "--interval-merging-rule", "-imr");
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
                return static_cast<char>(std::toupper(character));
            });
            if (value == "ALL") options.interval_merging_rule = Options::IntervalMergingRule::All;
            else if (value == "OVERLAPPING_ONLY")
                options.interval_merging_rule = Options::IntervalMergingRule::OverlappingOnly;
            else throw std::invalid_argument(
                "BAD_INPUT: --interval-merging-rule must be ALL or OVERLAPPING_ONLY");
        } else if (argument == "-ip" || argument == "--interval-padding" ||
                   !inline_value(argument, "--interval-padding").empty()) {
            const auto value = value_for(index, argc, argv, argument, "--interval-padding", "-ip");
            std::size_t consumed = 0;
            try { options.interval_padding = std::stoll(value, &consumed); }
            catch (...) { throw std::invalid_argument("BAD_INPUT: invalid --interval-padding: " + value); }
            if (consumed != value.size() || options.interval_padding < 0)
                throw std::invalid_argument("BAD_INPUT: --interval-padding must be non-negative");
        } else if (argument == "-ixp" || argument == "--interval-exclusion-padding" ||
                   !inline_value(argument, "--interval-exclusion-padding").empty()) {
            const auto value = value_for(index, argc, argv, argument,
                                         "--interval-exclusion-padding", "-ixp");
            std::size_t consumed = 0;
            try { options.interval_exclusion_padding = std::stoll(value, &consumed); }
            catch (...) {
                throw std::invalid_argument(
                    "BAD_INPUT: invalid --interval-exclusion-padding: " + value);
            }
            if (consumed != value.size() || options.interval_exclusion_padding < 0)
                throw std::invalid_argument(
                    "BAD_INPUT: --interval-exclusion-padding must be non-negative");
        } else if (argument == "-XL" || argument == "--exclude-intervals" ||
                   !inline_value(argument, "--exclude-intervals").empty()) {
            options.exclude_intervals.push_back(value_for(index, argc, argv, argument,
                                                           "--exclude-intervals", "-XL"));
        } else if (argument == "-RF" || argument == "--read-filter" ||
                   !inline_value(argument, "--read-filter").empty()) {
            options.read_filters.push_back(value_for(index, argc, argv, argument,
                                                      "--read-filter", "-RF"));
            if (options.read_filters.back() == "ReadLengthReadFilter")
                options.require_read_length = true;
        } else if (argument == "-XRF" || argument == "--inverted-read-filter" ||
                   !inline_value(argument, "--inverted-read-filter").empty()) {
            options.inverted_read_filters.push_back(value_for(
                index, argc, argv, argument, "--inverted-read-filter", "-XRF"));
            if (options.inverted_read_filters.back() == "ReadLengthReadFilter")
                options.require_read_length = true;
        } else if (argument == "-DF" || argument == "--disable-read-filter" ||
                   !inline_value(argument, "--disable-read-filter").empty()) {
            options.disabled_read_filters.push_back(value_for(index, argc, argv, argument,
                                                               "--disable-read-filter", "-DF"));
            if (options.disabled_read_filters.back() == "ReadLengthReadFilter")
                options.require_read_length = false;
        } else if (argument == "--read-name" || !inline_value(argument, "--read-name").empty()) {
            options.read_names.push_back(value_for(index, argc, argv, argument, "--read-name"));
        } else if (argument == "--keep-read-group" ||
                   !inline_value(argument, "--keep-read-group").empty()) {
            options.keep_read_groups.push_back(
                value_for(index, argc, argv, argument, "--keep-read-group"));
        } else if (argument == "--read-group-black-list" ||
                   !inline_value(argument, "--read-group-black-list").empty()) {
            const auto expression = value_for(index, argc, argv, argument,
                                              "--read-group-black-list");
            const auto separator = expression.find(':');
            if (separator != 2 || separator + 1 >= expression.size())
                throw std::invalid_argument(
                    "BAD_INPUT: --read-group-black-list must be ATTR:VALUE");
            options.read_group_blacklist.emplace_back(
                expression.substr(0, separator), expression.substr(separator + 1));
        } else if (argument == "--read-filter-tag" ||
                   !inline_value(argument, "--read-filter-tag").empty()) {
            options.read_filter_tag = value_for(index, argc, argv, argument, "--read-filter-tag");
        } else if (argument == "--read-filter-tag-comp" ||
                   !inline_value(argument, "--read-filter-tag-comp").empty()) {
            const auto value = value_for(index, argc, argv, argument, "--read-filter-tag-comp");
            std::size_t consumed = 0;
            try { options.read_filter_tag_comp = std::stod(value, &consumed); }
            catch (...) { throw std::invalid_argument("BAD_INPUT: invalid --read-filter-tag-comp: " + value); }
            if (consumed != value.size() || !std::isfinite(options.read_filter_tag_comp))
                throw std::invalid_argument("BAD_INPUT: --read-filter-tag-comp must be finite");
            options.read_filter_tag_comp_seen = true;
        } else if (argument == "--read-filter-tag-op" ||
                   !inline_value(argument, "--read-filter-tag-op").empty()) {
            auto value = value_for(index, argc, argv, argument, "--read-filter-tag-op");
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
                return static_cast<char>(std::toupper(character));
            });
            if (value == "LESS") options.read_filter_tag_op = Options::ReadTagOperator::Less;
            else if (value == "LESS_OR_EQUAL")
                options.read_filter_tag_op = Options::ReadTagOperator::LessOrEqual;
            else if (value == "GREATER") options.read_filter_tag_op = Options::ReadTagOperator::Greater;
            else if (value == "GREATER_OR_EQUAL")
                options.read_filter_tag_op = Options::ReadTagOperator::GreaterOrEqual;
            else if (value == "EQUAL") options.read_filter_tag_op = Options::ReadTagOperator::Equal;
            else if (value == "NOT_EQUAL") options.read_filter_tag_op = Options::ReadTagOperator::NotEqual;
            else throw std::invalid_argument(
                "BAD_INPUT: --read-filter-tag-op must be LESS, LESS_OR_EQUAL, GREATER, "
                "GREATER_OR_EQUAL, EQUAL or NOT_EQUAL");
            options.read_filter_tag_op_seen = true;
        } else if (argument == "--disable-tool-default-read-filters") {
            // There is no extra implicit predicate in this read-QC adapter;
            // accept the standard GATK switch for command-line compatibility
            // and preserve all explicitly requested -RF predicates.
            options.disable_tool_default_read_filters = true;
        } else if (argument == "--min-read-length" ||
                   !inline_value(argument, "--min-read-length").empty()) {
            const auto value = std::stoull(value_for(index, argc, argv, argument, "--min-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-read-length is out of range");
            options.min_read_length = static_cast<std::uint32_t>(value);
            options.require_read_length = true;
        } else if (argument == "--max-read-length" ||
                   !inline_value(argument, "--max-read-length").empty()) {
            const auto value = std::stoull(value_for(index, argc, argv, argument, "--max-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-read-length is out of range");
            options.max_read_length = static_cast<std::uint32_t>(value);
            options.require_read_length = true;
        } else if (argument == "--minimum-mapping-quality" ||
                   !inline_value(argument, "--minimum-mapping-quality").empty()) {
            const auto value = std::stoll(value_for(index, argc, argv, argument,
                                                    "--minimum-mapping-quality"));
            if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
                throw std::invalid_argument("--minimum-mapping-quality is out of range");
            options.min_mapping_quality = static_cast<int>(value);
            options.mapping_quality_argument_seen = true;
        } else if (argument == "--maximum-mapping-quality" ||
                   !inline_value(argument, "--maximum-mapping-quality").empty()) {
            const auto value = std::stoll(value_for(index, argc, argv, argument,
                                                    "--maximum-mapping-quality"));
            if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
                throw std::invalid_argument("--maximum-mapping-quality is out of range");
            options.max_mapping_quality = static_cast<int>(value);
            options.mapping_quality_argument_seen = true;
        } else if (argument == "--min-fragment-length" ||
                   !inline_value(argument, "--min-fragment-length").empty()) {
            const auto value = std::stoll(value_for(index, argc, argv, argument,
                                                    "--min-fragment-length"));
            if (value < 0 || value > std::numeric_limits<int>::max())
                throw std::invalid_argument("--min-fragment-length is out of range");
            options.min_fragment_length = static_cast<int>(value);
            options.fragment_length_argument_seen = true;
        } else if (argument == "--max-fragment-length" ||
                   !inline_value(argument, "--max-fragment-length").empty()) {
            const auto value = std::stoll(value_for(index, argc, argv, argument,
                                                    "--max-fragment-length"));
            if (value < 0 || value > std::numeric_limits<int>::max())
                throw std::invalid_argument("--max-fragment-length is out of range");
            options.max_fragment_length = static_cast<int>(value);
            options.fragment_length_argument_seen = true;
        } else if (argument == "-O" || argument == "--output" ||
                   !inline_value(argument, "--output").empty()) {
            options.output = value_for(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--batch-records" || !inline_value(argument, "--batch-records").empty()) {
            const auto value = value_for(index, argc, argv, argument, "--batch-records");
            const auto parsed = std::stoull(value);
            if (parsed == 0) throw std::invalid_argument("--batch-records must be positive");
            options.batch_records = static_cast<std::size_t>(parsed);
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(value_for(index, argc, argv, argument, "--threads"));
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() ||
                   !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = value_for(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--seconds-between-progress-updates", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                     "--seconds-between-progress-updates");
                (void)value_for(index, argc, argv, argument, name);
            }
        } else if (argument == "--QUIET" || argument == "--quiet" ||
                   argument.rfind("--QUIET=", 0) == 0 ||
                   argument.rfind("--quiet=", 0) == 0) {
            // Picard/GATK launcher quiet switch (P1 CLI parity).
            options.quiet = true;
        } else if (argument == "--tmp-dir" || !inline_value(argument, "--tmp-dir").empty()) {
            options.tmp_dir = value_for(index, argc, argv, argument, "--tmp-dir");
        } else if (argument == "--verbosity" ||
                   !inline_value(argument, "--verbosity").empty()) {
            options.verbosity = value_for(index, argc, argv, argument, "--verbosity");
        } else if (argument == "--java-options" ||
                   !inline_value(argument, "--java-options").empty() ||
                   argument == "-java-options" ||
                   !inline_value(argument, "-java-options").empty()) {
            // Launcher-level JVM passthrough: accepted for CLI parity; the
            // native binary has no JVM.
            (void)value_for(index, argc, argv, argument,
                            argument.rfind("-java-", 0) == 0 ? "-java-options"
                                                             : "--java-options");
        } else if (argument == "--gatk-config-file" ||
                   !inline_value(argument, "--gatk-config-file").empty()) {
            (void)value_for(index, argc, argv, argument, "--gatk-config-file");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    if (options.require_read_length && options.min_read_length > options.max_read_length)
        throw std::invalid_argument("read length filter minimum exceeds maximum");
    for (const auto& filter : options.read_filters)
        if (!is_supported_read_filter(filter))
            throw std::invalid_argument("UNSUPPORTED_PARAMETER: unsupported read filter: " + filter);
    for (const auto& filter : options.inverted_read_filters)
        if (!is_supported_read_filter(filter))
            throw std::invalid_argument("UNSUPPORTED_PARAMETER: unsupported inverted read filter: " + filter);
    for (const auto& filter : options.disabled_read_filters)
        if (!is_supported_read_filter(filter))
            throw std::invalid_argument("UNSUPPORTED_PARAMETER: unsupported disabled read filter: " + filter);
    const auto has_filter = [&](const char* name) {
        return std::find(options.read_filters.begin(), options.read_filters.end(), name) !=
                   options.read_filters.end() ||
               std::find(options.inverted_read_filters.begin(), options.inverted_read_filters.end(), name) !=
                   options.inverted_read_filters.end();
    };
    const bool has_mapping_quality_filter = has_filter("MappingQualityReadFilter");
    if (options.mapping_quality_argument_seen && !has_mapping_quality_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --minimum/maximum-mapping-quality requires MappingQualityReadFilter");
    const bool has_fragment_length_filter = has_filter("FragmentLengthReadFilter");
    if (options.fragment_length_argument_seen && !has_fragment_length_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --min/max-fragment-length requires FragmentLengthReadFilter");
    const bool has_read_name_filter = has_filter("ReadNameReadFilter");
    if (!options.read_names.empty() && !has_read_name_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --read-name requires ReadNameReadFilter");
    const bool has_read_group_filter = has_filter("ReadGroupReadFilter");
    if (!options.keep_read_groups.empty() && !has_read_group_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --keep-read-group requires ReadGroupReadFilter");
    const bool has_read_group_blacklist_filter = has_filter("ReadGroupBlackListReadFilter");
    if (!options.read_group_blacklist.empty() && !has_read_group_blacklist_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --read-group-black-list requires ReadGroupBlackListReadFilter");
    const bool has_read_tag_filter = has_filter("ReadTagValueFilter");
    if ((options.read_filter_tag.empty() || options.read_filter_tag.size() != 2) && has_read_tag_filter)
        throw std::invalid_argument(
            "BAD_INPUT: ReadTagValueFilter requires a two-character --read-filter-tag");
    if (!options.read_filter_tag.empty() && !has_read_tag_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --read-filter-tag requires ReadTagValueFilter");
    if ((options.read_filter_tag_comp_seen || options.read_filter_tag_op_seen) && !has_read_tag_filter)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: --read-filter-tag-comp/op requires ReadTagValueFilter");
    for (const auto& input : options.inputs)
        if (!std::filesystem::is_regular_file(input))
            throw std::runtime_error("BAD_INPUT: input is not a readable regular file: " + input);
    return options;
}

bool is_supported_read_filter(const std::string& name) {
    static const std::array<std::string_view, 30> supported{
        "AllowAllReadsReadFilter", "WellformedReadFilter", "MappedReadFilter",
        "NotDuplicateReadFilter", "NotSecondaryAlignmentReadFilter",
        "NotSupplementaryAlignmentReadFilter", "PrimaryLineReadFilter",
        "PassesVendorQualityCheckReadFilter", "ProperlyPairedReadFilter",
        "PairedReadFilter", "FirstOfPairReadFilter", "SecondOfPairReadFilter",
        "MappingQualityNotZeroReadFilter", "MappingQualityAvailableReadFilter", "MappingQualityReadFilter",
        "NonZeroReferenceLengthAlignmentReadFilter", "ReadLengthReadFilter", "NotProperlyPairedReadFilter",
        "NotOpticalDuplicateReadFilter", "NonZeroFragmentLengthReadFilter",
        "MateOnSameContigOrNoMappedMateReadFilter", "CigarContainsNoNOperator",
        "GoodCigarReadFilter", "HasReadGroupReadFilter", "MateDifferentStrandReadFilter",
        "FragmentLengthReadFilter", "ReadNameReadFilter", "ReadGroupReadFilter",
        "ReadTagValueFilter", "ReadGroupBlackListReadFilter"};
    return std::find(supported.begin(), supported.end(), name) != supported.end();
}

bool good_cigar(const bam1_t* record) {
    const auto* cigar = bam_get_cigar(record);
    bool previous_indel = false;
    for (std::uint32_t index = 0; index < record->core.n_cigar; ++index) {
        const auto element = cigar[index];
        const auto operation = bam_cigar_op(element);
        if (bam_cigar_oplen(element) == 0) return false;
        if (operation == BAM_CREF_SKIP) return false;
        const bool indel = operation == BAM_CINS || operation == BAM_CDEL;
        if (indel && previous_indel) return false;
        previous_indel = indel;
    }
    // GATK's CigarUtils.isGood rejects a deletion as the first/last
    // non-clipping operation, including when clips precede/follow it.
    for (int direction = 0; direction < 2; ++direction) {
        const auto count = record->core.n_cigar;
        for (std::uint32_t offset = 0; offset < count; ++offset) {
            const auto index = direction == 0 ? offset : count - offset - 1;
            const auto operation = bam_cigar_op(cigar[index]);
            const bool clipping = operation == BAM_CSOFT_CLIP ||
                                  operation == BAM_CHARD_CLIP || operation == BAM_CPAD;
            if (operation == BAM_CDEL) return false;
            if (!clipping) break;
        }
    }
    return true;
}

bool wellformed_read(const bam1_t* record, const sam_hdr_t* header) {
    // This mirrors the engine-level WellformedReadFilter composite used by
    // CountReads and FlagStat: mapped coordinates agree with the header,
    // sequence/RG data are present, read length matches CIGAR, and skipped
    // reference regions are rejected.  HTSlib has already validated the SAM
    // lexical/CIGAR encoding before this boundary is reached.
    if (record == nullptr || header == nullptr || record->core.l_qseq <= 0 ||
        bam_aux_get(record, "RG") == nullptr)
        return false;
    const bool unmapped = (record->core.flag & BAM_FUNMAP) != 0 ||
                          record->core.tid < 0 || record->core.pos < 0;
    if (!unmapped) {
        if (record->core.tid >= header->n_targets || record->core.pos < 0 ||
            bam_endpos(record) <= record->core.pos || record->core.n_cigar == 0)
            return false;
    }
    if (!unmapped && static_cast<int>(bam_cigar2qlen(record->core.n_cigar,
                                                      bam_get_cigar(record))) != record->core.l_qseq)
        return false;
    return good_cigar(record);
}

bool passes_read_filter(const Options& options, const std::string& name,
                        const bam1_t* record, const sam_hdr_t* header) {
    const auto flag = record->core.flag;
    if (name == "AllowAllReadsReadFilter") return true;
    if (name == "WellformedReadFilter") return wellformed_read(record, header);
    if (name == "MappedReadFilter") return (flag & BAM_FUNMAP) == 0 && record->core.tid >= 0;
    if (name == "NotDuplicateReadFilter") return (flag & BAM_FDUP) == 0;
    if (name == "NotSecondaryAlignmentReadFilter") return (flag & BAM_FSECONDARY) == 0;
    if (name == "NotSupplementaryAlignmentReadFilter") return (flag & BAM_FSUPPLEMENTARY) == 0;
    if (name == "PrimaryLineReadFilter") return (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
    if (name == "PassesVendorQualityCheckReadFilter") return (flag & BAM_FQCFAIL) == 0;
    if (name == "ProperlyPairedReadFilter") return (flag & (BAM_FPAIRED | BAM_FPROPER_PAIR)) ==
                                                             (BAM_FPAIRED | BAM_FPROPER_PAIR);
    if (name == "PairedReadFilter") return (flag & BAM_FPAIRED) != 0;
    if (name == "FirstOfPairReadFilter") return (flag & (BAM_FPAIRED | BAM_FREAD1)) ==
                                                        (BAM_FPAIRED | BAM_FREAD1);
    if (name == "SecondOfPairReadFilter") return (flag & (BAM_FPAIRED | BAM_FREAD2)) ==
                                                         (BAM_FPAIRED | BAM_FREAD2);
    if (name == "MappingQualityNotZeroReadFilter") return record->core.qual != 0;
    if (name == "MappingQualityAvailableReadFilter") return record->core.qual != 255;
    if (name == "MappingQualityReadFilter")
        return static_cast<int>(record->core.qual) >= options.min_mapping_quality &&
               (!options.max_mapping_quality.has_value() ||
                static_cast<int>(record->core.qual) <= *options.max_mapping_quality);
    if (name == "NonZeroFragmentLengthReadFilter") return record->core.isize != 0;
    if (name == "FragmentLengthReadFilter") {
        if ((flag & BAM_FPAIRED) == 0) return true;
        const auto fragment = std::llabs(static_cast<long long>(record->core.isize));
        return fragment >= options.min_fragment_length && fragment <= options.max_fragment_length;
    }
    if (name == "ReadNameReadFilter") {
        const auto* qname = bam_get_qname(record);
        return qname != nullptr && std::find(options.read_names.begin(), options.read_names.end(), qname) !=
            options.read_names.end();
    }
    if (name == "ReadGroupReadFilter") {
        const auto* tag = bam_aux_get(record, "RG");
        const auto* value = tag == nullptr ? nullptr : bam_aux2Z(tag);
        return value != nullptr && std::find(options.keep_read_groups.begin(),
                                             options.keep_read_groups.end(), value) !=
            options.keep_read_groups.end();
    }
    if (name == "ReadGroupBlackListReadFilter") {
        const auto* rg_tag = bam_aux_get(record, "RG");
        if (rg_tag == nullptr) return true;
        const auto* rg_id = bam_aux2Z(rg_tag);
        if (rg_id == nullptr) return true;
        for (const auto& [attribute, blocked] : options.read_group_blacklist) {
            if (attribute == "RG" && blocked == rg_id) return false;
            kstring_t value{0, 0, nullptr};
            const auto status = sam_hdr_find_tag_id(
                const_cast<sam_hdr_t*>(header), "RG", "ID", rg_id,
                attribute.c_str(), &value);
            const bool match = status == 0 && value.s != nullptr && blocked == value.s;
            std::free(value.s);
            if (match) return false;
        }
        return true;
    }
    if (name == "ReadTagValueFilter") {
        const auto* tag = bam_aux_get(record, options.read_filter_tag.c_str());
        if (tag == nullptr) return false;
        const auto value = bam_aux2f(tag);
        switch (options.read_filter_tag_op) {
        case Options::ReadTagOperator::Less: return value < options.read_filter_tag_comp;
        case Options::ReadTagOperator::LessOrEqual: return value <= options.read_filter_tag_comp;
        case Options::ReadTagOperator::Greater: return value > options.read_filter_tag_comp;
        case Options::ReadTagOperator::GreaterOrEqual: return value >= options.read_filter_tag_comp;
        case Options::ReadTagOperator::Equal: return value == options.read_filter_tag_comp;
        case Options::ReadTagOperator::NotEqual: return value != options.read_filter_tag_comp;
        }
        return false;
    }
    if (name == "MateOnSameContigOrNoMappedMateReadFilter") {
        const bool paired = (flag & BAM_FPAIRED) != 0;
        const bool unmapped = (flag & BAM_FUNMAP) != 0 || record->core.tid < 0 || record->core.pos < 0;
        const bool mate_unmapped = (flag & BAM_FMUNMAP) != 0 || record->core.mtid < 0 ||
                                   record->core.mpos < 0;
        return !paired || mate_unmapped ||
               (!unmapped && record->core.tid == record->core.mtid);
    }
    if (name == "CigarContainsNoNOperator") {
        const auto* cigar = bam_get_cigar(record);
        for (std::uint32_t index = 0; index < record->core.n_cigar; ++index)
            if (bam_cigar_op(cigar[index]) == BAM_CREF_SKIP) return false;
        return true;
    }
    if (name == "GoodCigarReadFilter") return good_cigar(record);
    if (name == "HasReadGroupReadFilter") return bam_aux_get(record, "RG") != nullptr;
    if (name == "MateDifferentStrandReadFilter") {
        const bool paired = (flag & BAM_FPAIRED) != 0;
        const bool unmapped = (flag & BAM_FUNMAP) != 0 || record->core.tid < 0 || record->core.pos < 0;
        const bool mate_unmapped = (flag & BAM_FMUNMAP) != 0 || record->core.mtid < 0 ||
                                   record->core.mpos < 0;
        const bool reverse = (flag & BAM_FREVERSE) != 0;
        const bool mate_reverse = (flag & BAM_FMREVERSE) != 0;
        return paired && !unmapped && !mate_unmapped && reverse != mate_reverse;
    }
    if (name == "NonZeroReferenceLengthAlignmentReadFilter")
        return (flag & BAM_FUNMAP) == 0 && record->core.tid >= 0 && record->core.pos >= 0 &&
               bam_endpos(record) > record->core.pos;
    if (name == "ReadLengthReadFilter")
        return !options.require_read_length ||
               (record->core.l_qseq >= static_cast<int>(options.min_read_length) &&
                record->core.l_qseq <= static_cast<int>(options.max_read_length));
    if (name == "NotProperlyPairedReadFilter")
        return (flag & BAM_FPAIRED) != 0 && (flag & BAM_FPROPER_PAIR) == 0;
    if (name == "NotOpticalDuplicateReadFilter") {
        // HTSlib exposes the duplicate bit but not Picard's optical-duplicate
        // classification without parsing the OQ/read-name metadata.  Keep
        // the common non-optical boundary conservative: all records pass.
        return true;
    }
    return false;
}

bool passes_read_filters(const Options& options, const bam1_t* record,
                         const sam_hdr_t* header) {
    const bool default_disabled = std::find(options.disabled_read_filters.begin(),
                                            options.disabled_read_filters.end(),
                                            "WellformedReadFilter") != options.disabled_read_filters.end();
    if (!options.disable_tool_default_read_filters && !default_disabled &&
        !passes_read_filter(options, "WellformedReadFilter", record, header))
        return false;
    for (const auto& filter : options.read_filters) {
        if (std::find(options.disabled_read_filters.begin(), options.disabled_read_filters.end(), filter) !=
            options.disabled_read_filters.end()) continue;
        if (!passes_read_filter(options, filter, record, header)) return false;
    }
    // GATK wraps each -XRF predicate in InvertedReadFilter after ordinary
    // default/-RF filters have been applied.
    for (const auto& filter : options.inverted_read_filters) {
        if (std::find(options.disabled_read_filters.begin(), options.disabled_read_filters.end(), filter) !=
            options.disabled_read_filters.end()) continue;
        if (passes_read_filter(options, filter, record, header)) return false;
    }
    return true;
}

std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    const auto offset = value.size() - ending.size();
    for (std::size_t index = 0; index < ending.size(); ++index) {
        if (std::tolower(static_cast<unsigned char>(value[offset + index])) !=
            std::tolower(static_cast<unsigned char>(ending[index]))) return false;
    }
    return true;
}

int header_tid(const sam_hdr_t* header, const std::string& contig) {
    for (int tid = 0; tid < header->n_targets; ++tid)
        if (contig == header->target_name[tid]) return tid;
    return -1;
}

std::int64_t parse_coordinate(const std::string& text, const std::string& original) {
    if (text.empty()) throw std::runtime_error("BAD_INPUT: empty interval coordinate: " + original);
    std::size_t consumed = 0;
    const auto value = std::stoll(text, &consumed);
    if (consumed != text.size() || value < 1)
        throw std::runtime_error("BAD_INPUT: invalid interval coordinate: " + original);
    return value;
}

Interval parse_literal_interval(const sam_hdr_t* header, const std::string& raw) {
    const auto value = trim_copy(raw);
    const auto colon = value.find(':');
    const auto contig = value.substr(0, colon);
    const auto tid = header_tid(header, contig);
    if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in header: " + contig);
    if (colon == std::string::npos)
        return Interval{tid, 0, static_cast<std::int64_t>(header->target_len[tid])};
    const auto coordinates = value.substr(colon + 1);
    const auto dash = coordinates.find('-');
    const auto start = parse_coordinate(coordinates.substr(0, dash), value);
    const auto end = dash == std::string::npos
        ? start : parse_coordinate(coordinates.substr(dash + 1), value);
    if (end < start || end > static_cast<std::int64_t>(header->target_len[tid]))
        throw std::runtime_error("BAD_INPUT: interval outside contig: " + value);
    return Interval{tid, start - 1, end};
}

void append_interval_file(const sam_hdr_t* header, const std::string& path,
                          std::vector<Interval>& output) {
    const bool bed = suffix_ci(path, ".bed") || suffix_ci(path, ".bed.gz");
    std::size_t parsed = 0;
    const auto process_line = [&](std::string line) {
        line = trim_copy(std::move(line));
        if (line.empty() || line.front() == '#' || line.front() == '@' ||
            line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) return;
        std::istringstream fields(line);
        std::string first;
        fields >> first;
        if (first.find(':') != std::string::npos) {
            output.push_back(parse_literal_interval(header, first));
            ++parsed;
            return;
        }
        std::string second, third;
        if (!(fields >> second >> third))
            throw std::runtime_error("BAD_INPUT: malformed interval file line: " + line);
        const auto tid = header_tid(header, first);
        if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in header: " + first);
        std::size_t consumed_start = 0, consumed_end = 0;
        const auto start = std::stoll(second, &consumed_start);
        const auto end = std::stoll(third, &consumed_end);
        if (consumed_start != second.size() || consumed_end != third.size())
            throw std::runtime_error("BAD_INPUT: malformed interval coordinates: " + line);
        const auto begin = bed ? start : start - 1;
        const auto finish = end;
        if (begin < 0 || finish <= begin || finish > static_cast<std::int64_t>(header->target_len[tid]))
            throw std::runtime_error("BAD_INPUT: interval outside contig: " + line);
        output.push_back(Interval{tid, begin, finish});
        ++parsed;
    };
    if (suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file)
            throw std::runtime_error("BAD_INPUT: cannot open compressed interval file: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (fastgatk::io::read_text_line(file, &line, path) >= 0)
                process_line(line.s == nullptr ? std::string{} : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(file);
            throw;
        }
        const auto close_status = hts_close(file);
        free(line.s);
        if (close_status != 0)
            throw std::runtime_error("BAD_INPUT: failed reading compressed interval file: " + path);
    } else {
        std::ifstream stream(path);
        if (!stream) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        for (std::string line; std::getline(stream, line);)
            process_line(std::move(line));
        if (stream.bad()) throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
    }
    if (parsed == 0) throw std::runtime_error("BAD_INPUT: interval file contains no intervals: " + path);
}

std::vector<Interval> materialize_intervals(const sam_hdr_t* header,
                                            const std::vector<std::string>& selectors,
                                            Options::IntervalSetRule set_rule =
                                                Options::IntervalSetRule::Union,
                                            Options::IntervalMergingRule merging_rule =
                                                Options::IntervalMergingRule::All,
                                            std::int64_t padding = 0) {
    if (selectors.empty()) return {};
    // A literal or interval-list file is one GATK interval selector.  Build
    // each selector independently before applying INTERSECTION, otherwise
    // intersecting the flattened records would incorrectly require a read to
    // overlap every line in one interval file.
    std::vector<std::vector<Interval>> selector_sets;
    selector_sets.reserve(selectors.size());
    for (const auto& selector : selectors) {
        std::vector<Interval> parsed;
        std::error_code error;
        if (std::filesystem::is_regular_file(selector, error) && !error)
            append_interval_file(header, selector, parsed);
        else
            parsed.push_back(parse_literal_interval(header, selector));
        if (padding != 0) {
            for (auto& interval : parsed) {
                interval.start = std::max<std::int64_t>(0, interval.start - padding);
                const auto length = static_cast<std::int64_t>(header->target_len[interval.tid]);
                interval.end = interval.end > std::numeric_limits<std::int64_t>::max() - padding
                    ? length : std::min(length, interval.end + padding);
            }
        }
        std::sort(parsed.begin(), parsed.end(), [](const Interval& left, const Interval& right) {
            if (left.tid != right.tid) return left.tid < right.tid;
            if (left.start != right.start) return left.start < right.start;
            return left.end < right.end;
        });
        std::vector<Interval> normalized;
        for (const auto& interval : parsed) {
            const bool merge = merging_rule == Options::IntervalMergingRule::All
                ? interval.start <= (normalized.empty() ? interval.start : normalized.back().end)
                : interval.start < (normalized.empty() ? interval.start : normalized.back().end);
            if (!normalized.empty() && normalized.back().tid == interval.tid &&
                merge)
                normalized.back().end = std::max(normalized.back().end, interval.end);
            else normalized.push_back(interval);
        }
        selector_sets.push_back(std::move(normalized));
    }
    if (set_rule == Options::IntervalSetRule::Intersection) {
        std::vector<Interval> result = selector_sets.front();
        for (std::size_t set = 1; set < selector_sets.size() && !result.empty(); ++set) {
            std::vector<Interval> intersection;
            for (const auto& left : result) {
                for (const auto& right : selector_sets[set]) {
                    if (left.tid != right.tid) continue;
                    const auto start = std::max(left.start, right.start);
                    const auto end = std::min(left.end, right.end);
                    if (start < end) intersection.push_back(Interval{left.tid, start, end});
                }
            }
            result.clear();
            for (const auto& interval : intersection) {
                const bool merge = merging_rule == Options::IntervalMergingRule::All
                    ? interval.start <= (result.empty() ? interval.start : result.back().end)
                    : interval.start < (result.empty() ? interval.start : result.back().end);
                if (!result.empty() && result.back().tid == interval.tid && merge)
                    result.back().end = std::max(result.back().end, interval.end);
                else result.push_back(interval);
            }
        }
        return result;
    }
    std::vector<Interval> intervals;
    for (auto& selector : selector_sets)
        intervals.insert(intervals.end(), selector.begin(), selector.end());
    std::sort(intervals.begin(), intervals.end(), [](const Interval& left, const Interval& right) {
        if (left.tid != right.tid) return left.tid < right.tid;
        if (left.start != right.start) return left.start < right.start;
        return left.end < right.end;
    });
    std::vector<Interval> merged;
    for (const auto& interval : intervals) {
        const bool merge = merging_rule == Options::IntervalMergingRule::All
            ? interval.start <= (merged.empty() ? interval.start : merged.back().end)
            : interval.start < (merged.empty() ? interval.start : merged.back().end);
        if (!merged.empty() && merged.back().tid == interval.tid && merge) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }
    return merged;
}

// HTSlib's multi-region iterator performs the same read de-duplication that
// GATK's interval traversal guarantees when an alignment overlaps multiple
// interval fragments.  Build its owned region list from our normalized
// zero-based half-open intervals so OVERLAPPING_ONLY does not double count
// records at interval boundaries.
hts_reglist_t* make_indexed_region_list(const sam_hdr_t* header,
                                        const std::vector<Interval>& intervals,
                                        unsigned int& count) {
    count = 0;
    if (intervals.empty()) return nullptr;
    std::vector<int> tids;
    std::vector<std::size_t> sizes;
    for (const auto& interval : intervals) {
        if (tids.empty() || tids.back() != interval.tid) {
            tids.push_back(interval.tid);
            sizes.push_back(0);
        }
        ++sizes.back();
    }
    auto* regions = static_cast<hts_reglist_t*>(
        std::calloc(tids.size(), sizeof(hts_reglist_t)));
    if (!regions) throw std::bad_alloc();
    std::size_t offset = 0;
    for (std::size_t index = 0; index < tids.size(); ++index) {
        regions[index].reg = sam_hdr_tid2name(header, tids[index]);
        regions[index].tid = tids[index];
        regions[index].count = static_cast<std::uint32_t>(sizes[index]);
        regions[index].intervals = static_cast<hts_pair_pos_t*>(
            std::malloc(sizes[index] * sizeof(hts_pair_pos_t)));
        if (!regions[index].intervals) {
            hts_reglist_free(regions, static_cast<int>(tids.size()));
            throw std::bad_alloc();
        }
        for (std::size_t item = 0; item < sizes[index]; ++item) {
            regions[index].intervals[item].beg = static_cast<hts_pos_t>(
                intervals[offset + item].start);
            regions[index].intervals[item].end = static_cast<hts_pos_t>(
                intervals[offset + item].end);
        }
        offset += sizes[index];
    }
    count = static_cast<unsigned int>(tids.size());
    return regions;
}

bool record_overlaps(const bam1_t* record, const std::vector<Interval>& intervals) {
    if (intervals.empty() || record->core.tid < 0 || record->core.pos < 0) return false;
    const auto start = static_cast<std::int64_t>(record->core.pos);
    auto end = static_cast<std::int64_t>(bam_endpos(record));
    if (end <= start) end = start + 1;
    return std::any_of(intervals.begin(), intervals.end(), [&](const Interval& interval) {
        return interval.tid == record->core.tid && start < interval.end && end > interval.start;
    });
}

// When no -L is provided GATK traverses the whole reference and subtracts
// -XL territory.  A read remains eligible if any part of its CIGAR span lies
// outside the excluded union; a one-base exclusion therefore must not drop
// every read spanning that base.
bool record_has_nonexcluded_span(const bam1_t* record,
                                 const std::vector<Interval>& intervals) {
    if (intervals.empty() || record->core.tid < 0 || record->core.pos < 0) return false;
    const auto start = static_cast<std::int64_t>(record->core.pos);
    auto end = static_cast<std::int64_t>(bam_endpos(record));
    if (end <= start) end = start + 1;
    std::int64_t cursor = start;
    for (const auto& interval : intervals) {
        if (interval.tid != record->core.tid) continue;
        if (interval.end <= cursor) continue;
        if (interval.start >= end) break;
        if (interval.start > cursor) return true;
        cursor = std::max(cursor, interval.end);
        if (cursor >= end) break;
    }
    return cursor < end;
}

bool has_remaining_interval_territory(const std::vector<Interval>& included,
                                      const std::vector<Interval>& excluded) {
    if (included.empty()) return true;
    for (const auto& interval : included) {
        std::int64_t cursor = interval.start;
        for (const auto& removal : excluded) {
            if (removal.tid != interval.tid) continue;
            if (removal.end <= cursor) continue;
            if (removal.start >= interval.end) break;
            if (removal.start > cursor) return true;
            cursor = std::max(cursor, removal.end);
            if (cursor >= interval.end) break;
        }
        if (cursor < interval.end) return true;
    }
    return false;
}

std::vector<Interval> subtract_intervals(const std::vector<Interval>& included,
                                         const std::vector<Interval>& excluded) {
    if (included.empty() || excluded.empty()) return included;
    std::vector<Interval> remaining;
    for (const auto& interval : included) {
        std::int64_t cursor = interval.start;
        for (const auto& removal : excluded) {
            if (removal.tid != interval.tid) continue;
            if (removal.end <= cursor) continue;
            if (removal.start >= interval.end) break;
            if (removal.start > cursor)
                remaining.push_back(Interval{interval.tid, cursor,
                                             std::min(removal.start, interval.end)});
            cursor = std::max(cursor, removal.end);
            if (cursor >= interval.end) break;
        }
        if (cursor < interval.end)
            remaining.push_back(Interval{interval.tid, cursor, interval.end});
    }
    return remaining;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char character : value) {
        if (character == '"' || character == '\\') out << '\\';
        if (character == '\n') out << "\\n";
        else if (character == '\r') out << "\\r";
        else if (character == '\t') out << "\\t";
        else out << character;
    }
    return out.str();
}

template<class ExecSpace>
KernelTotals process_batch(const std::vector<ReadProjection>& batch,
                           ReadMetricsBuffers<ExecSpace>& buffers) {
    KernelTotals totals;
    fastgatk::core::HostBatch host(std::string(kImplementation) + "-v1");
    host.records = batch.size();
    host.bytes = batch.size() * sizeof(ReadProjection);
    fastgatk::core::KernelPlan<ExecSpace> plan(kImplementation);
    plan.begin_prepare(host);

    buffers.ensure(batch.size());
    auto flags = buffers.flags;
    auto host_flags = Kokkos::create_mirror_view(Kokkos::subview(flags,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    for (std::size_t index = 0; index < batch.size(); ++index)
        host_flags(index) = batch[index].flags;
    auto flags_batch = Kokkos::subview(flags,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    Kokkos::deep_copy(flags_batch, host_flags);

    fastgatk::core::DeviceBatch<ExecSpace> device(batch.size());
    device.bind("flags", flags, batch.size());
#ifdef FASTGATK_FLAGSTAT
    auto tids = buffers.tids;
    auto mtids = buffers.mtids;
    auto positions = buffers.positions;
    auto mate_positions = buffers.mate_positions;
    auto mapqs = buffers.mapqs;
    auto counters = buffers.counters;
    auto tids_batch = Kokkos::subview(tids,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    auto mtids_batch = Kokkos::subview(mtids,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    auto positions_batch = Kokkos::subview(positions,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    auto mate_positions_batch = Kokkos::subview(mate_positions,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    auto mapqs_batch = Kokkos::subview(mapqs,
        std::make_pair<std::size_t, std::size_t>(0, batch.size()));
    auto host_tids = Kokkos::create_mirror_view(tids);
    auto host_mtids = Kokkos::create_mirror_view(mtids);
    auto host_positions = Kokkos::create_mirror_view(positions);
    auto host_mate_positions = Kokkos::create_mirror_view(mate_positions);
    auto host_mapqs = Kokkos::create_mirror_view(mapqs);
    for (std::size_t index = 0; index < batch.size(); ++index) {
        host_tids(index) = batch[index].tid;
        host_mtids(index) = batch[index].mtid;
        host_positions(index) = batch[index].pos;
        host_mate_positions(index) = batch[index].mpos;
        host_mapqs(index) = batch[index].mapq;
    }
    Kokkos::deep_copy(tids_batch, Kokkos::subview(host_tids,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    Kokkos::deep_copy(mtids_batch, Kokkos::subview(host_mtids,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    Kokkos::deep_copy(positions_batch, Kokkos::subview(host_positions,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    Kokkos::deep_copy(mate_positions_batch, Kokkos::subview(host_mate_positions,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    Kokkos::deep_copy(mapqs_batch, Kokkos::subview(host_mapqs,
        std::make_pair<std::size_t, std::size_t>(0, batch.size())));
    Kokkos::deep_copy(counters, std::uint64_t{0});
    device.bind("tids", tids, batch.size());
    device.bind("mtids", mtids, batch.size());
    device.bind("positions", positions, batch.size());
    device.bind("mate_positions", mate_positions, batch.size());
    device.bind("mapqs", mapqs, batch.size());
    device.bind("counters", counters);
#endif
    ExecSpace().fence();
    plan.end_prepare(device);
    plan.begin_execute();
#ifdef FASTGATK_FLAGSTAT
    Kokkos::parallel_for("flagstat_accumulate", Kokkos::RangePolicy<ExecSpace>(0, batch.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            const auto flag = flags(index);
            const bool unmapped = (flag & BAM_FUNMAP) != 0 || tids(index) < 0 || positions(index) < 0;
            const bool mate_unmapped = (flag & BAM_FMUNMAP) != 0 || mtids(index) < 0 || mate_positions(index) < 0;
            const bool paired = (flag & BAM_FPAIRED) != 0;
            Kokkos::atomic_fetch_add(&counters(0), std::uint64_t{1});
            if ((flag & BAM_FQCFAIL) != 0) Kokkos::atomic_fetch_add(&counters(1), std::uint64_t{1});
            if ((flag & BAM_FDUP) != 0) Kokkos::atomic_fetch_add(&counters(2), std::uint64_t{1});
            if (!unmapped) Kokkos::atomic_fetch_add(&counters(3), std::uint64_t{1});
            if (!paired) return;
            Kokkos::atomic_fetch_add(&counters(4), std::uint64_t{1});
            if ((flag & BAM_FREAD2) != 0) Kokkos::atomic_fetch_add(&counters(5), std::uint64_t{1});
            else if ((flag & BAM_FREAD1) != 0) Kokkos::atomic_fetch_add(&counters(6), std::uint64_t{1});
            if ((flag & BAM_FPROPER_PAIR) != 0) Kokkos::atomic_fetch_add(&counters(7), std::uint64_t{1});
            if (!unmapped && !mate_unmapped) {
                Kokkos::atomic_fetch_add(&counters(8), std::uint64_t{1});
                if (tids(index) != mtids(index)) {
                    Kokkos::atomic_fetch_add(&counters(10), std::uint64_t{1});
                    if (mapqs(index) >= 5) Kokkos::atomic_fetch_add(&counters(11), std::uint64_t{1});
                }
            }
            if (!unmapped && mate_unmapped) Kokkos::atomic_fetch_add(&counters(9), std::uint64_t{1});
        });
    ExecSpace().fence();
    plan.end_execute();
    auto host_counters = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counters);
    for (std::size_t index = 0; index < totals.counters.size(); ++index)
        totals.counters[index] += host_counters(index);
#else
    std::uint64_t count = 0;
    Kokkos::parallel_reduce("count_reads", Kokkos::RangePolicy<ExecSpace>(0, batch.size()),
        KOKKOS_LAMBDA(const std::size_t, std::uint64_t& value) { value += 1; }, count);
    ExecSpace().fence();
    plan.end_execute();
    totals.counters[0] += count;
#endif
    const auto& telemetry = plan.telemetry();
    totals.records += batch.size();
    ++totals.batches;
    totals.prepare_bytes += telemetry.prepare_bytes;
    totals.device_bytes += telemetry.device_bytes;
    totals.prepare_seconds += telemetry.prepare_seconds;
    totals.execute_seconds += telemetry.execute_seconds;
    totals.buffer_capacity_records = buffers.capacity;
    totals.buffer_allocations = buffers.allocations;
    totals.buffer_reuses = buffers.reuses;
    return totals;
}

#ifdef FASTGATK_FLAGSTAT
std::string percent(std::uint64_t numerator, std::uint64_t denominator) {
    if (denominator == 0) return "NaN";
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << (static_cast<double>(numerator) / static_cast<double>(denominator)) * 100.0;
    return out.str();
}
#endif

std::string render_output(const KernelTotals& totals) {
#ifdef FASTGATK_FLAGSTAT
    const auto& count = totals.counters;
    std::ostringstream out;
    out << count[0] << " in total\n"
        << count[1] << " QC failure\n"
        << count[2] << " duplicates\n"
        << count[3] << " mapped (" << percent(count[3], count[0]) << "%)\n"
        << count[4] << " paired in sequencing\n"
        << count[6] << " read1\n"
        << count[5] << " read2\n"
        << count[7] << " properly paired (" << percent(count[7], count[0]) << "%)\n"
        << count[8] << " with itself and mate mapped\n"
        << count[9] << " singletons (" << percent(count[9], count[0]) << "%)\n"
        << count[10] << " with mate mapped to a different chr\n"
        << count[11] << " with mate mapped to a different chr (mapQ>=5)\n";
    return out.str();
#else
    return std::to_string(totals.counters[0]) + "\n";
#endif
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
    const auto effective_threads = resources.effective_threads(
        static_cast<std::size_t>(options.threads));
    auto batch_records = options.batch_records;
    const auto initial_limits = fastgatk::runtime::AdaptiveController{}.initial(
        resources.budget(), fastgatk::runtime::WorkEstimate{
            static_cast<std::uint32_t>(std::min<std::size_t>(
                options.batch_records, std::numeric_limits<std::uint32_t>::max())),
            1, sizeof(ReadProjection), sizeof(ReadProjection), sizeof(ReadProjection)});
    if (initial_limits.max_reads != 0)
        batch_records = std::min<std::size_t>(batch_records, initial_limits.max_reads);
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(effective_threads);
    Kokkos::initialize(settings);
    bool initialized = true;
    KernelTotals totals;
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    ReadMetricsBuffers<ExecSpace> buffers;
    samFile* input = nullptr;
    sam_hdr_t* header = nullptr;
    bam1_t* record = nullptr;
    hts_idx_t* interval_index = nullptr;
    hts_itr_t* interval_iterator = nullptr;
    hts_reglist_t* indexed_region_list = nullptr;
    unsigned int indexed_region_count = 0;
    std::size_t interval_count = 0;
    std::size_t exclude_interval_count = 0;
    bool indexed_intervals = true;
    try {
        for (std::size_t input_index = 0; input_index < options.inputs.size(); ++input_index) {
            const auto& input_path = options.inputs[input_index];
            input = sam_open(input_path.c_str(), "r");
            if (!input)
                throw std::runtime_error("BAD_INPUT: HTSlib cannot open input: " + input_path);
            if (!options.reference.empty() && hts_set_fai_filename(input, options.reference.c_str()) != 0)
                throw std::runtime_error("BAD_INPUT: HTSlib cannot configure reference FASTA: " + options.reference);
            if (effective_threads > 1 && hts_set_threads(input, static_cast<int>(effective_threads)) != 0)
                throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib thread decoder setup failed");
            header = sam_hdr_read(input);
            record = bam_init1();
            if (!header || !record)
                throw std::runtime_error("BAD_INPUT: HTSlib cannot read input header");
        const auto include_intervals = materialize_intervals(
            header, options.intervals, options.interval_set_rule,
            options.interval_merging_rule, options.interval_padding);
        const auto exclude_intervals = materialize_intervals(
            header, options.exclude_intervals, Options::IntervalSetRule::Union,
            options.interval_merging_rule, options.interval_exclusion_padding);
        if (!include_intervals.empty() && !has_remaining_interval_territory(
                include_intervals, exclude_intervals))
            throw std::invalid_argument(
                "BAD_INPUT: intervals specified for exclusion with -XL removed all "
                "territory specified by -L");
        const auto effective_intervals = subtract_intervals(include_intervals, exclude_intervals);
        interval_count += include_intervals.size();
        exclude_interval_count += exclude_intervals.size();
        if (!effective_intervals.empty()) {
            interval_index = sam_index_load(input, input_path.c_str());
            if (interval_index != nullptr) {
                indexed_region_list = make_indexed_region_list(
                    header, effective_intervals, indexed_region_count);
                interval_iterator = sam_itr_regions(
                    interval_index, header, indexed_region_list, indexed_region_count);
                if (interval_iterator != nullptr) {
                    // Ownership of indexed_region_list transfers to the iterator.
                    indexed_region_list = nullptr;
                } else {
                    hts_reglist_free(indexed_region_list,
                                     static_cast<int>(indexed_region_count));
                    indexed_region_list = nullptr;
                }
            }
        }
        indexed_intervals = indexed_intervals &&
            (!include_intervals.empty() && interval_iterator != nullptr);

        // Keep decode, Kokkos compute and output aggregation in separate
        // stages.  The final sink is the only writer of totals, so output
        // ordering remains deterministic while HTSlib decode overlaps the
        // previous Kokkos batch.  Queue capacity is expressed in bytes and
        // is bounded to one logical input batch (plus small result queues).
        using ReadBatch = std::vector<ReadProjection>;
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<ReadBatch, KernelTotals, KernelTotals>;
        const auto saturating_multiply = [](std::uint64_t left, std::uint64_t right) {
            return right != 0 && left > std::numeric_limits<std::uint64_t>::max() / right
                ? std::numeric_limits<std::uint64_t>::max() : left * right;
        };
        const auto decoded_capacity = std::max<std::uint64_t>(
            sizeof(ReadProjection), saturating_multiply(
                static_cast<std::uint64_t>(batch_records), sizeof(ReadProjection)));
        const auto result_capacity = std::max<std::uint64_t>(sizeof(KernelTotals) * 2U, 1U);
        const auto pipeline_limits = Pipeline::Limits{
            decoded_capacity, result_capacity, result_capacity};
        bool input_exhausted = false;
        int status = 0;
        const auto next_record = [&]() -> bool {
            if (input_exhausted) return false;
            if (!effective_intervals.empty() && interval_index != nullptr) {
                if (interval_iterator != nullptr) {
                    status = sam_itr_next(input, interval_iterator, record);
                    if (status >= 0) return true;
                    hts_itr_destroy(interval_iterator);
                    interval_iterator = nullptr;
                    if (status < -1)
                        throw std::runtime_error(
                            "BAD_INPUT: HTSlib failed while reading indexed interval");
                    input_exhausted = true;
                    return false;
                }
            }
            status = sam_read1(input, header, record);
            if (status >= 0) return true;
            if (status < -1)
                throw std::runtime_error("BAD_INPUT: HTSlib failed while reading records");
            input_exhausted = true;
            return false;
        };
        Pipeline pipeline(
            pipeline_limits,
            [&]() -> std::optional<ReadBatch> {
                ReadBatch selected;
                selected.reserve(batch_records);
                while (selected.size() < batch_records && next_record()) {
                    const bool included = include_intervals.empty() ||
                        record_overlaps(record, effective_intervals);
                    const bool excluded = include_intervals.empty() &&
                        !exclude_intervals.empty() &&
                        !record_has_nonexcluded_span(record, exclude_intervals);
                    if (!included || excluded || !passes_read_filters(options, record, header))
                        continue;
                    selected.push_back(ReadProjection{
                        record->core.flag, record->core.tid, record->core.mtid,
                        static_cast<std::int64_t>(record->core.pos),
                        static_cast<std::int64_t>(record->core.mpos),
                        static_cast<std::int32_t>(record->core.isize),
                        record->core.qual});
                }
                if (selected.empty()) return std::nullopt;
                return selected;
            },
            [&](ReadBatch batch) -> std::optional<KernelTotals> {
                return process_batch<ExecSpace>(batch, buffers);
            },
            [](KernelTotals result) -> std::optional<KernelTotals> {
                return result;
            },
            [&](KernelTotals result) {
                totals.records += result.records;
                totals.batches += result.batches;
                totals.prepare_bytes += result.prepare_bytes;
                totals.device_bytes += result.device_bytes;
                totals.prepare_seconds += result.prepare_seconds;
                totals.execute_seconds += result.execute_seconds;
                for (std::size_t index = 0; index < totals.counters.size(); ++index)
                    totals.counters[index] += result.counters[index];
                totals.buffer_capacity_records = std::max(
                    totals.buffer_capacity_records, result.buffer_capacity_records);
                totals.buffer_allocations = result.buffer_allocations;
                totals.buffer_reuses = result.buffer_reuses;
            },
            [&](const ReadBatch& value) {
                return saturating_multiply(static_cast<std::uint64_t>(value.size()),
                                           sizeof(ReadProjection));
            },
            [](const KernelTotals&) { return static_cast<std::uint64_t>(sizeof(KernelTotals)); },
            [](const KernelTotals&) { return static_cast<std::uint64_t>(sizeof(KernelTotals)); });
        const auto pipeline_metrics = pipeline.run();
        totals.pipeline_decoded_items += pipeline_metrics.decoded_items;
        totals.pipeline_computed_items += pipeline_metrics.computed_items;
        totals.pipeline_encoded_items += pipeline_metrics.encoded_items;
        totals.pipeline_decoded_bytes += pipeline_metrics.decoded_bytes;
        totals.pipeline_computed_bytes += pipeline_metrics.computed_bytes;
        totals.pipeline_encoded_bytes += pipeline_metrics.encoded_bytes;
        totals.pipeline_peak_decoded_bytes = std::max(
            totals.pipeline_peak_decoded_bytes, pipeline_metrics.peak_decoded_bytes);
        totals.pipeline_peak_computed_bytes = std::max(
            totals.pipeline_peak_computed_bytes, pipeline_metrics.peak_computed_bytes);
        totals.pipeline_peak_encoded_bytes = std::max(
            totals.pipeline_peak_encoded_bytes, pipeline_metrics.peak_encoded_bytes);
        if (input_index + 1 == options.inputs.size()) {
        const auto text = render_output(totals);
        std::cout << text;
        if (!options.output.empty()) {
            std::ofstream output(options.output);
            if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write output: " + options.output);
            output << text;
            if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize output: " + options.output);
        }
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"" << kToolName
                     << "\",\"implementation\":\"" << kImplementation
                     << "\",\"status\":\"contract-compatible\",\"input\":\"" << json_escape(options.inputs.front())
                     << "\",\"inputs\":[";
            for (std::size_t index = 0; index < options.inputs.size(); ++index) {
                if (index != 0) manifest << ',';
                manifest << '"' << json_escape(options.inputs[index]) << '"';
            }
            manifest << "],\"input_count\":" << options.inputs.size()
                     << ",\"records\":" << totals.records << ",\"batches\":" << totals.batches
                     << ",\"batch_records\":" << batch_records
                     << ",\"requested_threads\":" << options.threads
                     << ",\"effective_threads\":" << effective_threads
                     << ",\"reference\":\"" << json_escape(options.reference) << "\""
                     << ",\"interval_selectors\":" << options.intervals.size()
                     << ",\"intervals\":" << interval_count
                     << ",\"exclude_interval_selectors\":" << options.exclude_intervals.size()
                     << ",\"exclude_intervals\":" << exclude_interval_count
                     << ",\"interval_set_rule\":\""
                     << (options.interval_set_rule == Options::IntervalSetRule::Intersection
                             ? "INTERSECTION" : "UNION") << "\""
                     << ",\"interval_merging_rule\":\""
                     << (options.interval_merging_rule == Options::IntervalMergingRule::OverlappingOnly
                             ? "OVERLAPPING_ONLY" : "ALL") << "\""
                     << ",\"interval_padding\":" << options.interval_padding
                     << ",\"interval_exclusion_padding\":"
                     << options.interval_exclusion_padding
                     << ",\"indexed_intervals\":"
                     << ((!options.intervals.empty() && indexed_intervals) ? "true" : "false")
                     << ",\"read_filter_count\":" << options.read_filters.size()
                     << ",\"inverted_read_filter_count\":" << options.inverted_read_filters.size()
                     << ",\"disabled_read_filter_count\":" << options.disabled_read_filters.size()
                     << ",\"read_name_count\":" << options.read_names.size()
                     << ",\"keep_read_group_count\":" << options.keep_read_groups.size()
                     << ",\"read_group_blacklist_count\":" << options.read_group_blacklist.size()
                     << ",\"read_filter_tag\":\"" << json_escape(options.read_filter_tag) << "\""
                     << ",\"read_filter_tag_comp\":" << options.read_filter_tag_comp
                     << ",\"read_filter_tag_op\":\""
                     << (options.read_filter_tag_op == Options::ReadTagOperator::Less ? "LESS" :
                         options.read_filter_tag_op == Options::ReadTagOperator::LessOrEqual ? "LESS_OR_EQUAL" :
                         options.read_filter_tag_op == Options::ReadTagOperator::Greater ? "GREATER" :
                         options.read_filter_tag_op == Options::ReadTagOperator::GreaterOrEqual ? "GREATER_OR_EQUAL" :
                         options.read_filter_tag_op == Options::ReadTagOperator::NotEqual ? "NOT_EQUAL" : "EQUAL")
                     << "\""
                     << ",\"default_read_filter_count\":"
                     << ((!options.disable_tool_default_read_filters &&
                          std::find(options.disabled_read_filters.begin(),
                                    options.disabled_read_filters.end(),
                                    "WellformedReadFilter") == options.disabled_read_filters.end()) ? 1 : 0)
                     << ",\"default_read_filters\":["
                     << ((!options.disable_tool_default_read_filters &&
                          std::find(options.disabled_read_filters.begin(),
                                    options.disabled_read_filters.end(),
                                    "WellformedReadFilter") == options.disabled_read_filters.end())
                             ? "\"WellformedReadFilter\"" : "") << "]"
                     << ",\"disable_tool_default_read_filters\":"
                     << (options.disable_tool_default_read_filters ? "true" : "false")
                     << ",\"read_filter_require_read_length\":"
                     << (options.require_read_length ? "true" : "false")
                     << ",\"read_filter_min_read_length\":" << options.min_read_length
                     << ",\"read_filter_max_read_length\":" << options.max_read_length
                     << ",\"read_filter_min_mapping_quality\":" << options.min_mapping_quality
                     << ",\"read_filter_max_mapping_quality\":"
                     << (options.max_mapping_quality.has_value()
                             ? std::to_string(*options.max_mapping_quality) : "null")
                     << ",\"read_filter_min_fragment_length\":" << options.min_fragment_length
                     << ",\"read_filter_max_fragment_length\":" << options.max_fragment_length
                     << ",\"output\":\"" << json_escape(options.output)
                     << "\",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\""
                     << ExecSpace::name() << "\",\"kernel_execution_policy\":\""
#ifdef FASTGATK_FLAGSTAT
                     << "RangePolicy+atomic\""
#else
                     << "RangePolicy+reduction\""
#endif
                     << ",\"kernel_batches\":" << totals.batches
                     << ",\"kernel_prepare_bytes\":" << totals.prepare_bytes
                     << ",\"kernel_device_bytes\":" << totals.device_bytes
                     << ",\"kernel_prepare_seconds\":" << totals.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << totals.execute_seconds
                     << ",\"persistent_buffer_capacity_records\":" << totals.buffer_capacity_records
                     << ",\"persistent_buffer_allocations\":" << totals.buffer_allocations
                     << ",\"persistent_buffer_reuses\":" << totals.buffer_reuses
                     << ",\"pipeline_lifecycle\":\"decode->compute->encode->sink\""
                     << ",\"pipeline_decoded_items\":" << totals.pipeline_decoded_items
                     << ",\"pipeline_computed_items\":" << totals.pipeline_computed_items
                     << ",\"pipeline_encoded_items\":" << totals.pipeline_encoded_items
                     << ",\"pipeline_decoded_bytes\":" << totals.pipeline_decoded_bytes
                     << ",\"pipeline_computed_bytes\":" << totals.pipeline_computed_bytes
                     << ",\"pipeline_encoded_bytes\":" << totals.pipeline_encoded_bytes
                     << ",\"pipeline_peak_decoded_bytes\":" << totals.pipeline_peak_decoded_bytes
                     << ",\"pipeline_peak_computed_bytes\":" << totals.pipeline_peak_computed_bytes
                     << ",\"pipeline_peak_encoded_bytes\":" << totals.pipeline_peak_encoded_bytes << "}}\n";
        }
        std::cerr << "{\"tool\":\"" << kToolName << "\",\"status\":\"contract-compatible\",\"records\":"
                  << totals.records << ",\"batches\":" << totals.batches
                  << ",\"requested_threads\":" << options.threads
                  << ",\"effective_threads\":" << effective_threads
                  << ",\"batch_records\":" << batch_records
                  << ",\"kernel_execution_space\":\"" << ExecSpace::name()
                  << "\",\"kernel_execute_seconds\":" << totals.execute_seconds
                  << ",\"persistent_buffer_allocations\":" << totals.buffer_allocations
                  << ",\"persistent_buffer_reuses\":" << totals.buffer_reuses
                  << ",\"pipeline_decoded_items\":" << totals.pipeline_decoded_items
                  << ",\"pipeline_computed_items\":" << totals.pipeline_computed_items
                  << ",\"pipeline_encoded_items\":" << totals.pipeline_encoded_items << "}\n";
        }  // final-input output and telemetry
        if (interval_iterator) {
            hts_itr_destroy(interval_iterator);
            interval_iterator = nullptr;
        }
        if (indexed_region_list) {
            hts_reglist_free(indexed_region_list,
                             static_cast<int>(indexed_region_count));
            indexed_region_list = nullptr;
            indexed_region_count = 0;
        }
        if (interval_index) {
            hts_idx_destroy(interval_index);
            interval_index = nullptr;
        }
        if (record) {
            bam_destroy1(record);
            record = nullptr;
        }
        if (header) {
            sam_hdr_destroy(header);
            header = nullptr;
        }
        if (input) {
            const int close_status = sam_close(input);
            input = nullptr;
            if (close_status != 0)
                throw std::runtime_error("BAD_INPUT: HTSlib failed closing input");
        }
        }  // input loop
        buffers.clear();
        Kokkos::finalize();
        initialized = false;
    } catch (...) {
        if (interval_iterator) hts_itr_destroy(interval_iterator);
        if (indexed_region_list)
            hts_reglist_free(indexed_region_list,
                             static_cast<int>(indexed_region_count));
        if (interval_index) hts_idx_destroy(interval_index);
        buffers.clear();
        if (initialized) Kokkos::finalize();
        if (record) bam_destroy1(record);
        if (header) sam_hdr_destroy(header);
        if (input) sam_close(input);
        throw;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return main_impl(argc, argv);
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
