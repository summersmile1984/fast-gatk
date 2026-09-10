#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/kernels/bqsr.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/sam.h>
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string reference;
    std::vector<std::string> known_sites;
    std::string recalibration;
    std::string output;
    std::string manifest;
    std::string index_path;
    // GATK's output-index switch is an optional boolean.  Keep indexing
    // enabled by default, but honor an explicit false/0 value so a native
    // output has the same sidecar contract as the Java tool.
    bool create_output_bam_index = true;
    std::vector<std::string> regions;
    std::size_t interval_file_inputs = 0;
    std::size_t interval_file_records = 0;
    std::size_t batch_records = 4096;
    int preserve_qualities_less_than = 6;
    int quantization_levels = 0;
    std::vector<int> static_quantized_quals;
    bool round_down_quantized = false;
    bool allow_missing_read_group = false;
    double global_qscore_prior = -1.0;
    bool use_original_qualities = false;
    bool emit_original_qualities = false;
    bool compute_indel_bqsr_tables = false;
    // BaseRecalibrator inherits GATK's read-filter plugin descriptor.  Keep
    // the resolved default mask explicit so -RF/-DF and
    // --disable-tool-default-read-filters have the same meaning as Java,
    // rather than accepting compatibility switches and silently ignoring
    // them.
    bool disable_tool_default_read_filters = false;
    bool mapping_quality_not_zero_filter = true;
    bool mapping_quality_available_filter = true;
    bool mapped_filter = true;
    bool not_secondary_filter = true;
    bool not_duplicate_filter = true;
    bool passes_vendor_quality_filter = true;
    bool wellformed_filter = true;
    bool not_supplementary_filter = false;
    std::vector<std::string> enabled_read_filters;
    std::vector<std::string> disabled_read_filters;
    // The GATK plugin descriptor applies inverted filters after the regular
    // filter conjunction. Store names for manifest provenance and a compact
    // resolved mask for the per-record Host predicate.
    std::vector<std::string> inverted_read_filters;
    std::uint16_t inverted_read_filter_mask = 0;
    int indels_context_size = 3;
    int mismatches_context_size = 2;
    // CycleCovariate has a bounded signed cycle key. GATK defaults to 500
    // and rejects a long read above that bound; expose the same model
    // dimension so long reads can opt into a larger deterministic domain.
    int maximum_cycle_value = 500;
    int insertions_default_quality = 45;
    int deletions_default_quality = 45;
    std::string checkpoint;
    bool resume_checkpoint = false;
    std::size_t checkpoint_every_batches = 1;
    // Resource-aware decode telemetry.  Checkpointed runs keep a stable
    // batch boundary; non-checkpointed runs may reduce the next decode batch
    // after observing HostBatch bytes under SLURM/cgroup pressure.
    std::uint64_t adaptive_batch_reductions = 0;
    // The first decode batch after applying the allocation cap.  Keep this
    // separate from batch_records (the user's requested schedule) so the
    // manifest is auditable when a constrained worker starts below that cap.
    std::size_t initial_batch_records = 0;
    std::size_t effective_batch_records = 0;
};

// A small, tool-neutral projection of ThreeStagePipeline metrics.  Keeping
// this outside the HTSlib-only code lets the manifest schema describe the
// lifecycle for both BaseRecalibrator and ApplyBQSR without leaking BAM
// ownership into the runtime library.
struct PipelineSummary {
    bool used = false;
    std::uint64_t decoded_items = 0;
    std::uint64_t computed_items = 0;
    std::uint64_t encoded_items = 0;
    std::uint64_t decoded_bytes = 0;
    std::uint64_t computed_bytes = 0;
    std::uint64_t encoded_bytes = 0;
    std::uint64_t peak_decoded_bytes = 0;
    std::uint64_t peak_computed_bytes = 0;
    std::uint64_t peak_encoded_bytes = 0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* long_name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, long_name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == long_name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

enum class BqsrReadFilterKind : std::uint8_t {
    MappingQualityNotZero,
    MappingQualityAvailable,
    Mapped,
    NotSecondary,
    NotDuplicate,
    PassesVendorQuality,
    Wellformed,
    NotSupplementary,
};

constexpr std::uint16_t bqsr_read_filter_bit(BqsrReadFilterKind kind) noexcept {
    return static_cast<std::uint16_t>(1U << static_cast<unsigned>(kind));
}

[[maybe_unused]] BqsrReadFilterKind bqsr_read_filter_kind(const std::string& name) {
    if (name == "MappingQualityNotZeroReadFilter")
        return BqsrReadFilterKind::MappingQualityNotZero;
    if (name == "MappingQualityAvailableReadFilter")
        return BqsrReadFilterKind::MappingQualityAvailable;
    if (name == "MappedReadFilter") return BqsrReadFilterKind::Mapped;
    if (name == "NotSecondaryAlignmentReadFilter")
        return BqsrReadFilterKind::NotSecondary;
    if (name == "NotDuplicateReadFilter") return BqsrReadFilterKind::NotDuplicate;
    if (name == "PassesVendorQualityCheckReadFilter")
        return BqsrReadFilterKind::PassesVendorQuality;
    if (name == "WellformedReadFilter") return BqsrReadFilterKind::Wellformed;
    if (name == "NotSupplementaryAlignmentReadFilter")
        return BqsrReadFilterKind::NotSupplementary;
    throw std::invalid_argument("UNSUPPORTED_PARAMETER: read filter " + name);
}

[[maybe_unused]] bool bqsr_tool_default_filter(BqsrReadFilterKind kind) noexcept {
    switch (kind) {
        case BqsrReadFilterKind::MappingQualityNotZero:
        case BqsrReadFilterKind::MappingQualityAvailable:
        case BqsrReadFilterKind::Mapped:
        case BqsrReadFilterKind::NotSecondary:
        case BqsrReadFilterKind::NotDuplicate:
        case BqsrReadFilterKind::PassesVendorQuality:
        case BqsrReadFilterKind::Wellformed:
            return true;
        case BqsrReadFilterKind::NotSupplementary:
            return false;
    }
    return false;
}

[[maybe_unused]] void set_bqsr_read_filter(Options& options, BqsrReadFilterKind kind, bool enabled) {
    switch (kind) {
        case BqsrReadFilterKind::MappingQualityNotZero:
            options.mapping_quality_not_zero_filter = enabled; break;
        case BqsrReadFilterKind::MappingQualityAvailable:
            options.mapping_quality_available_filter = enabled; break;
        case BqsrReadFilterKind::Mapped:
            options.mapped_filter = enabled; break;
        case BqsrReadFilterKind::NotSecondary:
            options.not_secondary_filter = enabled; break;
        case BqsrReadFilterKind::NotDuplicate:
            options.not_duplicate_filter = enabled; break;
        case BqsrReadFilterKind::PassesVendorQuality:
            options.passes_vendor_quality_filter = enabled; break;
        case BqsrReadFilterKind::Wellformed:
            options.wellformed_filter = enabled; break;
        case BqsrReadFilterKind::NotSupplementary:
            options.not_supplementary_filter = enabled; break;
    }
}

[[maybe_unused]] void disable_bqsr_tool_defaults(Options& options) {
    options.mapping_quality_not_zero_filter = false;
    options.mapping_quality_available_filter = false;
    options.mapped_filter = false;
    options.not_secondary_filter = false;
    options.not_duplicate_filter = false;
    options.passes_vendor_quality_filter = false;
    options.wellformed_filter = false;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto ch : text) {
        if (ch == '"' || ch == '\\') escaped << '\\';
        if (ch == '\n') escaped << "\\n";
        else if (ch == '\r') escaped << "\\r";
        else if (ch == '\t') escaped << "\\t";
        else escaped << ch;
    }
    return escaped.str();
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
#if defined(FASTGATK_BQSR_APPLY)
            std::cout << "fastgatk-apply-bqsr (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 input BAM/CRAM\n"
                         "  --bqsr-recal-file, -bqsr FILE      FASTGATK-BQSR report\n"
                         "  --preserve-qscores-less-than N   leave qualities below N unchanged\n"
                         "  --use-original-qualities, -OQ  reset from OQ tag before recalibration\n"
                         "  --emit-original-quals           write original QUAL under OQ when absent\n"
                         "  --quantize-quals N               0=no quantization, -1=report map\n"
                         "  --static-quantized-quals N       repeatable static quality bins\n"
                         "  --round-down-quantized          round static bins down\n"
                         "  --allow-missing-read-group      use quantized qualities for unknown RG\n"
                         "  --create-output-bam-index [BOOL]  write BAM/CRAM index (default true)\n"
                         "  --global-qscore-prior Q           constant Bayesian prior (>0)\n"
                         "  --checkpoint FILE                atomic ApplyBQSR restart state\n"
                         "  --resume-checkpoint FILE        validate and restart from checkpoint\n"
                         "  --checkpoint-every-batches N    checkpoint cadence (default 1)\n"
                         "  -R, --reference FILE              CRAM reference FASTA\n"
                         "  -O, --output FILE                 output BAM/CRAM\n";
#else
            std::cout << "fastgatk-bqsr (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 input BAM/CRAM\n"
                         "  -R, --reference FILE              reference FASTA\n"
                         "  --known-sites FILE               VCF known sites (optional)\n"
                         "  --compute-indel-bqsr-tables     emit GATK I/D recalibration tables\n"
                         "  -RF, --read-filter NAME         add a supported GATK read filter\n"
                         "  -XRF, --inverted-read-filter NAME\n"
                         "                                  add the inverse of a supported read filter\n"
                         "  -DF, --disable-read-filter NAME disable one default read filter\n"
                         "      --disable-tool-default-read-filters [BOOL]\n"
                         "                                  disable all default read filters\n"
                         "  --indels-context-size N         indel context size (default 3)\n"
                         "  --mismatches-context-size N    mismatch context size (default 2)\n"
                         "  --maximum-cycle-value N        maximum absolute cycle (default 500)\n"
                         "  --insertions-default-quality N  indel insertion quality (default 45)\n"
                         "  --deletions-default-quality N   indel deletion quality (default 45)\n"
                         "  --checkpoint FILE                atomic BaseRecalibrator checkpoint\n"
                         "  --resume-checkpoint FILE        validate and resume checkpoint\n"
                         "  --checkpoint-every-batches N    checkpoint cadence (default 1)\n"
                         "  -L, --intervals REGION           contig:start-end\n"
                         "  -O, --output FILE                recalibration report\n";
#endif
            std::exit(0);
        } else if (is_option(argument, "--input") || argument == "-I")
            options.input = require_value(i, argc, argv, argument, "--input", "-I");
        else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(i, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--known-sites"))
            options.known_sites.push_back(require_value(i, argc, argv, argument, "--known-sites"));
        else if (
#if defined(FASTGATK_BQSR_APPLY)
                 argument == "-bqsr" ||
#endif
                 is_option(argument, "--bqsr-recal-file") ||
                 is_option(argument, "--recal-file"))
#if defined(FASTGATK_BQSR_APPLY)
            options.recalibration = require_value(i, argc, argv, argument,
                argument == "-bqsr" ? "--bqsr-recal-file" :
                argument.rfind("--recal-file", 0) == 0 ? "--recal-file" : "--bqsr-recal-file",
                "-bqsr");
#else
            options.recalibration = require_value(i, argc, argv, argument,
                argument.rfind("--recal-file", 0) == 0 ? "--recal-file" : "--bqsr-recal-file");
#endif
        else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(i, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(i, argc, argv, argument,
                argument == "--manifest" ? "--manifest" : "--output-manifest");
        else if (is_option(argument, "--intervals") || is_option(argument, "--region") || argument == "-L")
            options.regions.push_back(require_value(i, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                argument.rfind("--region", 0) == 0 ? "--region" : "--intervals", "-L"));
        else if (is_option(argument, "--batch-records"))
            options.batch_records = std::stoull(require_value(i, argc, argv, argument, "--batch-records"));
        else if (is_option(argument, "--preserve-qscores-less-than"))
            options.preserve_qualities_less_than = std::stoi(require_value(
                i, argc, argv, argument, "--preserve-qscores-less-than"));
        else if (is_option(argument, "--quantize-quals"))
            options.quantization_levels = std::stoi(require_value(
                i, argc, argv, argument, "--quantize-quals"));
        else if (is_option(argument, "--static-quantized-quals"))
            options.static_quantized_quals.push_back(std::stoi(require_value(
                i, argc, argv, argument, "--static-quantized-quals")));
        else if (argument == "--round-down-quantized" ||
                 argument.rfind("--round-down-quantized=", 0) == 0)
            options.round_down_quantized = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--round-down-quantized");
        else if (argument == "--allow-missing-read-group" ||
                 argument.rfind("--allow-missing-read-group=", 0) == 0)
            options.allow_missing_read_group = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--allow-missing-read-group");
        else if (is_option(argument, "--global-qscore-prior"))
            options.global_qscore_prior = std::stod(require_value(
                i, argc, argv, argument, "--global-qscore-prior"));
        else if (argument == "--use-original-qualities" || argument == "-OQ" ||
                 argument.rfind("--use-original-qualities=", 0) == 0) {
            const char* name = argument == "-OQ" ? "-OQ" : "--use-original-qualities";
            options.use_original_qualities = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, name);
        }
        else if (argument == "--emit-original-quals" ||
                 argument.rfind("--emit-original-quals=", 0) == 0)
            options.emit_original_qualities = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--emit-original-quals");
        else if (argument == "--compute-indel-bqsr-tables" || argument == "--indels")
            options.compute_indel_bqsr_tables = true;
#if !defined(FASTGATK_BQSR_APPLY)
        else if (is_option(argument, "--mismatches-context-size") || is_option(argument, "-mcs"))
            options.mismatches_context_size = std::stoi(require_value(
                i, argc, argv, argument,
                argument.rfind("-mcs", 0) == 0 ? "-mcs" : "--mismatches-context-size"));
        else if (argument == "--disable-tool-default-read-filters" ||
                 argument.rfind("--disable-tool-default-read-filters=", 0) == 0) {
            options.disable_tool_default_read_filters = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--disable-tool-default-read-filters");
            if (options.disable_tool_default_read_filters) disable_bqsr_tool_defaults(options);
        }
        else if (argument == "-RF" || is_option(argument, "--read-filter")) {
            const auto filter = require_value(i, argc, argv, argument, "--read-filter", "-RF");
            (void)bqsr_read_filter_kind(filter);
            options.enabled_read_filters.push_back(filter);
        }
        else if (argument == "-DF" || is_option(argument, "--disable-read-filter")) {
            const auto filter = require_value(i, argc, argv, argument,
                                               "--disable-read-filter", "-DF");
            (void)bqsr_read_filter_kind(filter);
            options.disabled_read_filters.push_back(filter);
        }
        else if (argument == "-XRF" || is_option(argument, "--inverted-read-filter")) {
            const auto filter = require_value(i, argc, argv, argument,
                                               "--inverted-read-filter", "-XRF");
            (void)bqsr_read_filter_kind(filter);
            options.inverted_read_filters.push_back(filter);
        }
        else if (is_option(argument, "--maximum-cycle-value") || is_option(argument, "--max-cycle"))
            options.maximum_cycle_value = std::stoi(require_value(
                i, argc, argv, argument,
                argument.rfind("--max-cycle", 0) == 0 ? "--max-cycle" : "--maximum-cycle-value"));
#endif
        else if (is_option(argument, "--indels-context-size") || is_option(argument, "-ics"))
            options.indels_context_size = std::stoi(require_value(
                i, argc, argv, argument,
                argument.rfind("-ics", 0) == 0 ? "-ics" : "--indels-context-size"));
        else if (is_option(argument, "--insertions-default-quality"))
            options.insertions_default_quality = std::stoi(require_value(
                i, argc, argv, argument, "--insertions-default-quality"));
        else if (is_option(argument, "--deletions-default-quality"))
            options.deletions_default_quality = std::stoi(require_value(
                i, argc, argv, argument, "--deletions-default-quality"));
        else if (is_option(argument, "--checkpoint"))
            options.checkpoint = require_value(i, argc, argv, argument, "--checkpoint");
        else if (is_option(argument, "--resume-checkpoint")) {
            options.checkpoint = require_value(i, argc, argv, argument, "--resume-checkpoint");
            options.resume_checkpoint = true;
        } else if (is_option(argument, "--checkpoint-every-batches"))
            options.checkpoint_every_batches = std::stoull(require_value(
                i, argc, argv, argument, "--checkpoint-every-batches"));
        else if (argument == "--create-output-bam-index" ||
                 argument.rfind("--create-output-bam-index=", 0) == 0)
            options.create_output_bam_index = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--create-output-bam-index");
        else if (argument == "--quiet") {
            // Standalone flag; do not consume the next tool argument.
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("--input is required");
    if (options.output.empty()) throw std::invalid_argument("--output is required");
    if (options.batch_records == 0) throw std::invalid_argument("--batch-records must be positive");
    if (options.checkpoint_every_batches == 0)
        throw std::invalid_argument("--checkpoint-every-batches must be positive");
    if (options.preserve_qualities_less_than < 0 || options.preserve_qualities_less_than > 93)
        throw std::invalid_argument("--preserve-qscores-less-than must be in [0,93]");
    if (options.quantization_levels > 93)
        throw std::invalid_argument("--quantize-quals must be <= 93 (0 disables quantization)");
    if (!options.static_quantized_quals.empty()) {
        if (options.quantization_levels != 0)
            throw std::invalid_argument("--static-quantized-quals cannot be combined with --quantize-quals");
        for (const auto quality : options.static_quantized_quals)
            if (quality < 0 || quality > 254)
                throw std::invalid_argument("--static-quantized-quals values must be in [0,254]");
    } else if (options.round_down_quantized) {
        throw std::invalid_argument("--round-down-quantized requires --static-quantized-quals");
    }
#if !defined(FASTGATK_BQSR_APPLY)
    if (!options.static_quantized_quals.empty() || options.round_down_quantized ||
        options.allow_missing_read_group)
        throw std::invalid_argument("static quantization and missing-read-group options require ApplyBQSR");
    // Barclay resolves plugin directives after argument parsing.  The native
    // implementation therefore applies the same order: start with either
    // the complete default list or an empty list, remove -DF defaults, then
    // append explicit -RF filters and finally apply -XRF inversions. Mirror
    // the descriptor's mutually-exclusive directive checks rather than
    // silently choosing one side of a contradictory command line.
    const auto duplicate_name = [](const std::vector<std::string>& names) {
        std::set<std::string> seen;
        for (const auto& name : names)
            if (!seen.insert(name).second) return name;
        return std::string{};
    };
    if (const auto duplicate = duplicate_name(options.enabled_read_filters); !duplicate.empty())
        throw std::invalid_argument("BAD_INPUT: read filter is enabled more than once: " + duplicate);
    if (const auto duplicate = duplicate_name(options.disabled_read_filters); !duplicate.empty())
        throw std::invalid_argument("BAD_INPUT: read filter is disabled more than once: " + duplicate);
    const auto contains = [](const std::vector<std::string>& names, const std::string& name) {
        return std::find(names.begin(), names.end(), name) != names.end();
    };
    for (const auto& filter : options.enabled_read_filters) {
        if (contains(options.disabled_read_filters, filter))
            throw std::invalid_argument("BAD_INPUT: read filter is both enabled and disabled: " + filter);
        if (contains(options.inverted_read_filters, filter))
            throw std::invalid_argument("BAD_INPUT: read filter is both enabled and inverted: " + filter);
    }
    for (const auto& filter : options.inverted_read_filters) {
        if (contains(options.disabled_read_filters, filter))
            throw std::invalid_argument("BAD_INPUT: read filter is both inverted and disabled: " + filter);
        const auto kind = bqsr_read_filter_kind(filter);
        if (!options.disable_tool_default_read_filters && bqsr_tool_default_filter(kind))
            throw std::invalid_argument(
                "BAD_INPUT: inverted default read filter requires --disable-tool-default-read-filters: " + filter);
        options.inverted_read_filter_mask = static_cast<std::uint16_t>(
            options.inverted_read_filter_mask | bqsr_read_filter_bit(kind));
    }
    if (options.disable_tool_default_read_filters) disable_bqsr_tool_defaults(options);
    for (const auto& filter : options.disabled_read_filters)
        set_bqsr_read_filter(options, bqsr_read_filter_kind(filter), false);
    for (const auto& filter : options.enabled_read_filters)
        set_bqsr_read_filter(options, bqsr_read_filter_kind(filter), true);
#endif
    if (options.global_qscore_prior == 0.0 || options.global_qscore_prior < -1.0 ||
        options.global_qscore_prior > 93.0)
        throw std::invalid_argument("--global-qscore-prior must be -1 or in (0,93]");
    if (options.indels_context_size <= 0 || options.indels_context_size > 13)
        throw std::invalid_argument("--indels-context-size must be in [1,13]");
    if (options.mismatches_context_size <= 0 || options.mismatches_context_size > 13)
        throw std::invalid_argument("--mismatches-context-size must be in [1,13]");
    if (options.maximum_cycle_value <= 0)
        throw std::invalid_argument("--maximum-cycle-value must be positive");
    if (options.insertions_default_quality < 0 || options.insertions_default_quality > 93 ||
        options.deletions_default_quality < 0 || options.deletions_default_quality > 93)
        throw std::invalid_argument("indel default qualities must be in [0,93]");
#if defined(FASTGATK_BQSR_APPLY)
    if (options.recalibration.empty()) throw std::invalid_argument("--bqsr-recal-file is required");
#endif
    return options;
}

std::uint8_t encode_base(std::uint8_t value) {
    switch (value) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

std::vector<std::string> load_reference(const std::string& path,
                                        const fastgatk::io::HeaderSummary& header) {
    std::vector<std::string> sequences(header.contigs.size());
    if (path.empty()) return sequences;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open reference FASTA: " + path);
    std::map<std::string, std::size_t> ids;
    for (std::size_t i = 0; i < header.contigs.size(); ++i) ids.emplace(header.contigs[i], i);
    std::string line, current;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            std::istringstream name(line.substr(1));
            name >> current;
            continue;
        }
        const auto found = ids.find(current);
        if (found != ids.end()) sequences[found->second] += line;
    }
    return sequences;
}

// Java's KnownSites mask is interval based: a VariantContext covers the
// complete REF span, not just the first base.  Keep intervals compact and
// merge them instead of expanding every long deletion into a potentially
// enormous point set.
struct KnownSiteInterval {
    std::int32_t tid = -1;
    std::int32_t start = -1;  // zero based, inclusive
    std::int32_t end = -1;    // zero based, inclusive
    bool operator<(const KnownSiteInterval& other) const noexcept {
        if (tid != other.tid) return tid < other.tid;
        if (start != other.start) return start < other.start;
        return end < other.end;
    }
};

std::string known_sites_signature(const Options& options) {
    std::ostringstream value;
    for (std::size_t index = 0; index < options.known_sites.size(); ++index) {
        if (index != 0) value << '\x1f';
        value << options.known_sites[index];
    }
    return value.str();
}

std::string static_quantized_signature(const Options& options) {
    std::ostringstream value;
    for (std::size_t index = 0; index < options.static_quantized_quals.size(); ++index) {
        if (index != 0) value << ',';
        value << options.static_quantized_quals[index];
    }
    return value.str();
}

std::string read_filter_signature(const std::vector<std::string>& filters) {
    std::ostringstream value;
    for (std::size_t index = 0; index < filters.size(); ++index) {
        if (index != 0) value << ',';
        value << filters[index];
    }
    return value.str();
}

std::vector<KnownSiteInterval> load_known_sites(
    const std::vector<std::string>& paths, const fastgatk::io::HeaderSummary& header) {
    std::vector<KnownSiteInterval> sites;
    if (paths.empty()) return sites;
#if FASTGATK_HAS_HTSLIB
    for (const auto& path : paths) {
        htsFile* input = bcf_open(path.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open known-sites VCF/BCF: " + path);
        bcf_hdr_t* header_vcf = bcf_hdr_read(input);
        bcf1_t* record = bcf_init();
        if (!header_vcf || !record) {
            if (record) bcf_destroy(record);
            if (header_vcf) bcf_hdr_destroy(header_vcf);
            bcf_close(input);
            throw std::runtime_error("BAD_INPUT: cannot read known-sites VCF/BCF header: " + path);
        }
        while (true) {
            const auto status = bcf_read(input, header_vcf, record);
            if (status < 0) {
                if (status != -1) {
                    bcf_destroy(record); bcf_hdr_destroy(header_vcf); bcf_close(input);
                    throw std::runtime_error("BAD_INPUT: failed reading known-sites VCF/BCF: " + path);
                }
                break;
            }
            if (record->rid < 0 || record->pos < 0) continue;
            const char* contig = bcf_hdr_id2name(header_vcf, record->rid);
            if (!contig) continue;
            const auto found = std::find(header.contigs.begin(), header.contigs.end(), contig);
            if (found == header.contigs.end()) continue;
            const auto tid = static_cast<std::int32_t>(std::distance(header.contigs.begin(), found));
            if (record->pos > std::numeric_limits<std::int32_t>::max()) continue;
            std::int64_t end = static_cast<std::int64_t>(record->pos) +
                std::max<std::int32_t>(record->rlen, 1) - 1;
            // Symbolic records can carry an explicit 1-based inclusive END;
            // VariantContext uses that span for KnownSites masking.
            int32_t* end_values = nullptr;
            int end_count = 0;
            const int end_status = bcf_get_info_int32(
                header_vcf, record, "END", &end_values, &end_count);
            if (end_status > 0 && end_values[0] > 0)
                end = static_cast<std::int64_t>(end_values[0]) - 1;
            free(end_values);
            if (end < record->pos) end = record->pos;
            if (end > std::numeric_limits<std::int32_t>::max())
                end = std::numeric_limits<std::int32_t>::max();
            sites.push_back(KnownSiteInterval{
                tid, static_cast<std::int32_t>(record->pos), static_cast<std::int32_t>(end)});
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header_vcf);
        if (bcf_close(input) != 0)
            throw std::runtime_error("BAD_INPUT: failed closing known-sites VCF/BCF: " + path);
    }
#else
    (void)header;
    throw std::runtime_error("BACKEND_UNAVAILABLE: known-sites requires HTSlib");
#endif
    std::sort(sites.begin(), sites.end());
    std::vector<KnownSiteInterval> merged;
    merged.reserve(sites.size());
    for (const auto& interval : sites) {
        if (!merged.empty() && merged.back().tid == interval.tid &&
            interval.start <= merged.back().end + 1) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }
    return merged;
}

bool known_site_contains(const std::vector<KnownSiteInterval>& sites,
                         std::int32_t tid, std::int64_t position) noexcept {
    if (tid < 0 || position < 0 || position > std::numeric_limits<std::int32_t>::max())
        return false;
    const auto needle = KnownSiteInterval{tid, static_cast<std::int32_t>(position),
                                          std::numeric_limits<std::int32_t>::max()};
    const auto found = std::upper_bound(sites.begin(), sites.end(), needle,
        [](const KnownSiteInterval& value, const KnownSiteInterval& item) {
            if (value.tid != item.tid) return value.tid < item.tid;
            return value.start < item.start;
        });
    if (found == sites.begin()) return false;
    const auto& candidate = *(found - 1);
    return candidate.tid == tid && candidate.start <= position && candidate.end >= position;
}

struct QualityBin {
    std::uint64_t count = 0;
    std::uint64_t mismatches = 0;
};

struct CovariateKey {
    std::string read_group;
    std::int32_t cycle = 0;
    std::string context;
    std::int32_t quality = 0;
    char event = 'M';
    bool operator<(const CovariateKey& other) const {
        if (read_group != other.read_group) return read_group < other.read_group;
        if (event != other.event) return event < other.event;
        if (cycle != other.cycle) return cycle < other.cycle;
        if (context != other.context) return context < other.context;
        return quality < other.quality;
    }
};

using CovariateTable = std::map<CovariateKey, QualityBin>;

struct BqsrCheckpointState {
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
    std::array<QualityBin, 94> bins{};
    CovariateTable covariates;
    std::size_t effective_batch_records = 0;
};

struct ApplyBqsrCheckpointState {
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
};

std::uintmax_t file_size_or_zero(const std::string& path) {
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    return error ? 0 : size;
}

std::intmax_t file_mtime_or_zero(const std::string& path) {
    std::error_code error;
    const auto stamp = std::filesystem::last_write_time(path, error);
    if (error) return 0;
    return static_cast<std::intmax_t>(stamp.time_since_epoch().count());
}

std::string regions_signature(const Options& options) {
    std::ostringstream value;
    for (std::size_t index = 0; index < options.regions.size(); ++index) {
        if (index != 0) value << '\x1f';
        value << options.regions[index];
    }
    return value.str();
}

void write_bqsr_checkpoint(const std::string& path, const Options& options,
                           const BqsrCheckpointState& state) {
    if (path.empty()) return;
    const std::string temporary = path + ".tmp";
    {
        std::ofstream output(temporary, std::ios::trunc);
        if (!output) throw std::runtime_error("cannot write BQSR checkpoint: " + temporary);
        output << "FASTGATK_BQSR_CHECKPOINT 1\n"
               << "tool BaseRecalibrator\n"
               << "input " << std::quoted(options.input) << "\n"
               << "reference " << std::quoted(options.reference) << "\n"
               << "known_sites " << std::quoted(known_sites_signature(options)) << "\n"
               << "regions " << std::quoted(regions_signature(options)) << "\n"
               << "input_size " << file_size_or_zero(options.input) << "\n"
               << "input_mtime " << file_mtime_or_zero(options.input) << "\n"
               << "batch_records " << options.batch_records << "\n"
               << "compute_indel " << (options.compute_indel_bqsr_tables ? 1 : 0) << "\n"
               << "indels_context_size " << options.indels_context_size << "\n"
               << "mismatches_context_size " << options.mismatches_context_size << "\n"
               << "maximum_cycle_value " << options.maximum_cycle_value << "\n"
               << "insertions_default_quality " << options.insertions_default_quality << "\n"
               << "deletions_default_quality " << options.deletions_default_quality << "\n"
               << "effective_batch_records " << options.effective_batch_records << "\n"
               << "records " << state.records << "\n"
               << "bases " << state.bases << "\n";
        for (std::size_t quality = 0; quality < state.bins.size(); ++quality)
            output << "bin " << quality << ' ' << state.bins[quality].count << ' '
                   << state.bins[quality].mismatches << "\n";
        for (const auto& [key, bin] : state.covariates)
            output << "cov " << std::quoted(key.read_group) << ' ' << key.cycle << ' '
                   << std::quoted(key.context) << ' ' << key.quality << ' ' << key.event << ' '
                   << bin.count << ' ' << bin.mismatches << "\n";
        output.flush();
        if (!output) throw std::runtime_error("cannot finalize BQSR checkpoint: " + temporary);
    }
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) throw std::runtime_error("cannot atomically publish BQSR checkpoint: " + path);
}

BqsrCheckpointState read_bqsr_checkpoint(const std::string& path, const Options& options) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open BQSR checkpoint: " + path);
    BqsrCheckpointState state;
    std::string line;
    bool header_seen = false;
    bool tool_seen = false;
    bool input_seen = false;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream fields(line);
        std::string tag;
        fields >> tag;
        if (tag == "FASTGATK_BQSR_CHECKPOINT") {
            int version = 0;
            fields >> version;
            if (version != 1) throw std::runtime_error("BAD_INPUT: unsupported BQSR checkpoint version");
            header_seen = true;
        } else if (tag == "tool") {
            std::string value;
            fields >> value;
            if (value != "BaseRecalibrator") throw std::runtime_error("BAD_INPUT: checkpoint tool mismatch");
            tool_seen = true;
        } else if (tag == "input" || tag == "reference" || tag == "known_sites" ||
                   tag == "region" || tag == "regions") {
            std::string value;
            fields >> std::quoted(value);
            const auto expected = tag == "input" ? options.input :
                tag == "reference" ? options.reference : known_sites_signature(options);
            const auto expected_regions = regions_signature(options);
            if (tag == "region") {
                const auto legacy_region = options.regions.size() == 1 ? options.regions.front() : std::string{};
                if (value != legacy_region)
                    throw std::runtime_error("BAD_INPUT: BQSR checkpoint option mismatch: region");
                continue;
            }
            if (tag == "regions") {
                if (value != expected_regions)
                    throw std::runtime_error("BAD_INPUT: BQSR checkpoint option mismatch: regions");
                continue;
            }
            if (value != expected) throw std::runtime_error("BAD_INPUT: BQSR checkpoint option mismatch: " + tag);
            if (tag == "input") input_seen = true;
        } else if (tag == "input_size" || tag == "input_mtime") {
            std::intmax_t value = 0;
            fields >> value;
            const auto expected = tag == "input_size"
                ? static_cast<std::intmax_t>(file_size_or_zero(options.input))
                : file_mtime_or_zero(options.input);
            if (value != expected) throw std::runtime_error("BAD_INPUT: BQSR input changed since checkpoint");
        } else if (tag == "batch_records") {
            std::size_t value = 0;
            fields >> value;
            if (value != options.batch_records) throw std::runtime_error("BAD_INPUT: BQSR batch size mismatch");
        } else if (tag == "compute_indel" || tag == "indels_context_size" ||
                   tag == "mismatches_context_size" ||
                   tag == "maximum_cycle_value" ||
                   tag == "insertions_default_quality" || tag == "deletions_default_quality") {
            int value = 0;
            fields >> value;
            const int expected = tag == "compute_indel" ? (options.compute_indel_bqsr_tables ? 1 : 0) :
                tag == "indels_context_size" ? options.indels_context_size :
                tag == "mismatches_context_size" ? options.mismatches_context_size :
                tag == "maximum_cycle_value" ? options.maximum_cycle_value :
                tag == "insertions_default_quality" ? options.insertions_default_quality :
                options.deletions_default_quality;
            if (value != expected) throw std::runtime_error("BAD_INPUT: BQSR checkpoint option mismatch: " + tag);
        } else if (tag == "records") {
            fields >> state.records;
        } else if (tag == "bases") {
            fields >> state.bases;
        } else if (tag == "effective_batch_records") {
            fields >> state.effective_batch_records;
        } else if (tag == "bin") {
            std::size_t quality = 0;
            QualityBin bin;
            fields >> quality >> bin.count >> bin.mismatches;
            if (quality >= state.bins.size()) throw std::runtime_error("BAD_INPUT: invalid BQSR checkpoint quality");
            state.bins[quality] = bin;
        } else if (tag == "cov") {
            CovariateKey key;
            QualityBin bin;
            fields >> std::quoted(key.read_group) >> key.cycle >> std::quoted(key.context) >>
                key.quality >> key.event >> bin.count >> bin.mismatches;
            if (!fields || (key.event != 'M' && key.event != 'I' && key.event != 'D'))
                throw std::runtime_error("BAD_INPUT: malformed BQSR checkpoint covariate");
            state.covariates[key] = bin;
        } else {
            throw std::runtime_error("BAD_INPUT: unknown BQSR checkpoint field: " + tag);
        }
    }
    if (!header_seen || !tool_seen || !input_seen || state.records == 0)
        throw std::runtime_error("BAD_INPUT: incomplete BQSR checkpoint");
    if (state.effective_batch_records != 0 &&
        state.effective_batch_records != options.effective_batch_records)
        throw std::runtime_error("BAD_INPUT: BQSR effective batch size mismatch");
    return state;
}

void write_apply_bqsr_checkpoint(const std::string& path, const Options& options,
                                 const ApplyBqsrCheckpointState& state) {
    if (path.empty()) return;
    const std::string temporary = path + ".tmp";
    {
        std::ofstream output(temporary, std::ios::trunc);
        if (!output) throw std::runtime_error("cannot write ApplyBQSR checkpoint: " + temporary);
        output << "FASTGATK_BQSR_CHECKPOINT 1\n"
               << "tool ApplyBQSR\n"
               << "input " << std::quoted(options.input) << "\n"
               << "reference " << std::quoted(options.reference) << "\n"
               << "recalibration " << std::quoted(options.recalibration) << "\n"
               << "input_size " << file_size_or_zero(options.input) << "\n"
               << "input_mtime " << file_mtime_or_zero(options.input) << "\n"
               << "recalibration_size " << file_size_or_zero(options.recalibration) << "\n"
               << "recalibration_mtime " << file_mtime_or_zero(options.recalibration) << "\n"
               << "batch_records " << options.batch_records << "\n"
               << "preserve_qualities_less_than " << options.preserve_qualities_less_than << "\n"
               << "quantization_levels " << options.quantization_levels << "\n"
               << "static_quantized_quals " << std::quoted(static_quantized_signature(options)) << "\n"
               << "round_down_quantized " << (options.round_down_quantized ? 1 : 0) << "\n"
               << "allow_missing_read_group " << (options.allow_missing_read_group ? 1 : 0) << "\n"
               << "use_original_qualities " << (options.use_original_qualities ? 1 : 0) << "\n"
               << "emit_original_qualities " << (options.emit_original_qualities ? 1 : 0) << "\n"
               << "create_output_bam_index " << (options.create_output_bam_index ? 1 : 0) << "\n"
               << "global_qscore_prior " << std::setprecision(17) << options.global_qscore_prior << "\n"
               << "records " << state.records << "\n"
               << "bases " << state.bases << "\n";
        output.flush();
        if (!output) throw std::runtime_error("cannot finalize ApplyBQSR checkpoint: " + temporary);
    }
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) throw std::runtime_error("cannot atomically publish ApplyBQSR checkpoint: " + path);
}

ApplyBqsrCheckpointState read_apply_bqsr_checkpoint(const std::string& path,
                                                    const Options& options) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open ApplyBQSR checkpoint: " + path);
    ApplyBqsrCheckpointState state;
    std::string line;
    bool header_seen = false;
    bool tool_seen = false;
    bool input_seen = false;
    bool recalibration_seen = false;
    bool use_original_seen = false;
    bool emit_original_seen = false;
    bool static_quantized_seen = false;
    bool round_down_seen = false;
    bool allow_missing_seen = false;
    bool create_index_seen = false;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream fields(line);
        std::string tag;
        fields >> tag;
        if (tag == "FASTGATK_BQSR_CHECKPOINT") {
            int version = 0;
            fields >> version;
            if (version != 1) throw std::runtime_error("BAD_INPUT: unsupported ApplyBQSR checkpoint version");
            header_seen = true;
        } else if (tag == "tool") {
            std::string value;
            fields >> value;
            if (value != "ApplyBQSR") throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint tool mismatch");
            tool_seen = true;
        } else if (tag == "input" || tag == "reference" || tag == "recalibration") {
            std::string value;
            fields >> std::quoted(value);
            const auto* expected = tag == "input" ? &options.input :
                tag == "reference" ? &options.reference : &options.recalibration;
            if (value != *expected)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: " + tag);
            if (tag == "input") input_seen = true;
            if (tag == "recalibration") recalibration_seen = true;
        } else if (tag == "input_size" || tag == "input_mtime" ||
                   tag == "recalibration_size" || tag == "recalibration_mtime") {
            std::intmax_t value = 0;
            fields >> value;
            const auto expected = tag == "input_size"
                ? static_cast<std::intmax_t>(file_size_or_zero(options.input))
                : tag == "input_mtime"
                ? file_mtime_or_zero(options.input)
                : tag == "recalibration_size"
                ? static_cast<std::intmax_t>(file_size_or_zero(options.recalibration))
                : file_mtime_or_zero(options.recalibration);
            if (value != expected)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR input/recalibration changed since checkpoint");
        } else if (tag == "batch_records") {
            std::size_t value = 0;
            fields >> value;
            if (value != options.batch_records)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR batch size mismatch");
        } else if (tag == "preserve_qualities_less_than" || tag == "quantization_levels") {
            int value = 0;
            fields >> value;
            const int expected = tag == "preserve_qualities_less_than"
                ? options.preserve_qualities_less_than : options.quantization_levels;
            if (value != expected)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: " + tag);
        } else if (tag == "static_quantized_quals") {
            std::string value;
            fields >> std::quoted(value);
            if (value != static_quantized_signature(options))
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: static_quantized_quals");
            static_quantized_seen = true;
        } else if (tag == "round_down_quantized" || tag == "allow_missing_read_group") {
            int value = 0;
            fields >> value;
            const int expected = tag == "round_down_quantized"
                ? (options.round_down_quantized ? 1 : 0)
                : (options.allow_missing_read_group ? 1 : 0);
            if (value != expected)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: " + tag);
            if (tag == "round_down_quantized") round_down_seen = true;
            else allow_missing_seen = true;
        } else if (tag == "use_original_qualities") {
            int value = 0;
            fields >> value;
            if (value != (options.use_original_qualities ? 1 : 0))
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: use_original_qualities");
            use_original_seen = true;
        } else if (tag == "emit_original_qualities") {
            int value = 0;
            fields >> value;
            if (value != (options.emit_original_qualities ? 1 : 0))
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: emit_original_qualities");
            emit_original_seen = true;
        } else if (tag == "create_output_bam_index") {
            int value = 0;
            fields >> value;
            if (value != (options.create_output_bam_index ? 1 : 0))
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: create_output_bam_index");
            create_index_seen = true;
        } else if (tag == "global_qscore_prior") {
            double value = 0.0;
            fields >> value;
            if (value != options.global_qscore_prior)
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint option mismatch: global_qscore_prior");
        } else if (tag == "records") {
            fields >> state.records;
        } else if (tag == "bases") {
            fields >> state.bases;
        } else {
            throw std::runtime_error("BAD_INPUT: unknown ApplyBQSR checkpoint field: " + tag);
        }
    }
    if (!header_seen || !tool_seen || !input_seen || !recalibration_seen || state.records == 0)
        throw std::runtime_error("BAD_INPUT: incomplete ApplyBQSR checkpoint");
    if (options.use_original_qualities && !use_original_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks use_original_qualities signature");
    if (options.emit_original_qualities && !emit_original_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks emit_original_qualities signature");
    if (!options.static_quantized_quals.empty() && !static_quantized_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks static_quantized_quals signature");
    if (options.round_down_quantized && !round_down_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks round_down_quantized signature");
    if (options.allow_missing_read_group && !allow_missing_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks allow_missing_read_group signature");
    // A pre-index-switch checkpoint historically implied the default (index
    // enabled).  Refuse to resume it with an explicit disable request because
    // that changes the durable output artifact contract.
    if (!options.create_output_bam_index && !create_index_seen)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint lacks create_output_bam_index signature");
    return state;
}

struct CovariateModel {
    std::map<std::string, QualityBin> read_group;
    std::map<std::string, double> read_group_reported_quality;
    std::map<std::pair<std::string, int>, QualityBin> read_group_quality;
    std::map<std::tuple<std::string, std::string, int>, QualityBin> context_quality;
    std::map<std::tuple<std::string, int, int>, QualityBin> cycle_quality;
};

void add_quality_bin(QualityBin& destination, const QualityBin& source) {
    destination.count += source.count;
    destination.mismatches += source.mismatches;
}

std::string read_group_at(const fastgatk::io::ReadBatch& batch, std::size_t record) {
    if (batch.read_group_offsets.size() != batch.records() + 1) return {};
    const auto begin = batch.read_group_offsets[record];
    const auto end = batch.read_group_offsets[record + 1];
    if (begin > end || end > batch.read_groups.size()) return {};
    return std::string(batch.read_groups.begin() + begin, batch.read_groups.begin() + end);
}

bool bqsr_wellformed_read(const fastgatk::io::ReadBatch& batch, std::size_t record,
                          const fastgatk::io::HeaderSummary& header) {
    if (record >= batch.records() || batch.offsets.size() != batch.records() + 1 ||
        batch.flags.size() != batch.records() || batch.tids.size() != batch.records() ||
        batch.positions.size() != batch.records() ||
        batch.read_group_offsets.size() != batch.records() + 1 ||
        batch.offsets[record] > batch.offsets[record + 1] ||
        batch.offsets[record + 1] > batch.bases.size() ||
        batch.offsets[record + 1] > batch.qualities.size() ||
        batch.offsets[record] == batch.offsets[record + 1] ||
        read_group_at(batch, record).empty())
        return false;
    const auto flags = batch.flags[record];
    const bool unmapped = (flags & 0x4U) != 0;
    if (!unmapped) {
        const auto tid = batch.tids[record];
        const auto position = batch.positions[record];
        if (tid < 0 || static_cast<std::size_t>(tid) >= header.contigs.size() || position < 0)
            return false;
        if (static_cast<std::size_t>(tid) < header.contig_lengths.size() &&
            position >= header.contig_lengths[static_cast<std::size_t>(tid)])
            return false;
    }
    // The HTSlib reader keeps a CIGAR offset span for every record.  Mirror
    // WellformedReadFilter's no-N and read-length checks without carrying a
    // bam1_t object into the Kokkos stage.
    if (!batch.cigar_layout_valid()) return false;
    std::uint64_t read_consumed = 0;
    for (std::size_t index = batch.cigar_offsets[record];
         index < batch.cigar_offsets[record + 1]; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid() || operation.code == fastgatk::io::CigarOpCode::ReferenceSkip)
            return false;
        if (operation.consumes_read()) read_consumed += operation.length;
    }
    if (!unmapped && read_consumed != batch.offsets[record + 1] - batch.offsets[record])
        return false;
    return true;
}

bool bqsr_read_filter_passes(const fastgatk::io::ReadBatch& batch, std::size_t record,
                             const fastgatk::io::HeaderSummary& header,
                             BqsrReadFilterKind kind) {
    if (record >= batch.records() || batch.flags.size() != batch.records() ||
        batch.mapq.size() != batch.records())
        return false;
    const auto flags = batch.flags[record];
    switch (kind) {
        case BqsrReadFilterKind::MappingQualityNotZero:
            return batch.mapq[record] != 0;
        case BqsrReadFilterKind::MappingQualityAvailable:
            return batch.mapq[record] != 255;
        case BqsrReadFilterKind::Mapped:
            return (flags & 0x4U) == 0;
        case BqsrReadFilterKind::NotSecondary:
            return (flags & 0x100U) == 0;
        case BqsrReadFilterKind::NotDuplicate:
            return (flags & 0x400U) == 0;
        case BqsrReadFilterKind::PassesVendorQuality:
            return (flags & 0x200U) == 0;
        case BqsrReadFilterKind::Wellformed:
            return bqsr_wellformed_read(batch, record, header);
        case BqsrReadFilterKind::NotSupplementary:
            return (flags & 0x800U) == 0;
    }
    return false;
}

bool bqsr_default_read_filter(const fastgatk::io::ReadBatch& batch, std::size_t record,
                              const fastgatk::io::HeaderSummary& header,
                              const Options& options) {
    if (record >= batch.records() || batch.flags.size() != batch.records() ||
        batch.mapq.size() != batch.records()) return false;
    const auto flags = batch.flags[record];
    // BaseRecalibrator's standard filter list is MappingQualityNotZero,
    // MappingQualityAvailable, Mapped, NotSecondary, NotDuplicate,
    // PassesVendorQualityCheck and WellformedReadFilter.  Keep each
    // predicate independently switchable so the generic GATK -RF/-DF
    // plugin controls retain their Java semantics.
    if (options.mapping_quality_not_zero_filter && batch.mapq[record] == 0)
        return false;
    if (options.mapping_quality_available_filter && batch.mapq[record] == 255)
        return false;
    if (options.mapped_filter && (flags & 0x4U) != 0) return false;
    if (options.not_secondary_filter && (flags & 0x100U) != 0) return false;
    if (options.not_duplicate_filter && (flags & 0x400U) != 0) return false;
    if (options.passes_vendor_quality_filter && (flags & 0x200U) != 0) return false;
    if (options.not_supplementary_filter && (flags & 0x800U) != 0) return false;
    if (options.wellformed_filter && !bqsr_wellformed_read(batch, record, header)) return false;
    constexpr std::array inverted_kinds{
        BqsrReadFilterKind::MappingQualityNotZero,
        BqsrReadFilterKind::MappingQualityAvailable,
        BqsrReadFilterKind::Mapped,
        BqsrReadFilterKind::NotSecondary,
        BqsrReadFilterKind::NotDuplicate,
        BqsrReadFilterKind::PassesVendorQuality,
        BqsrReadFilterKind::Wellformed,
        BqsrReadFilterKind::NotSupplementary,
    };
    for (const auto kind : inverted_kinds) {
        if ((options.inverted_read_filter_mask & bqsr_read_filter_bit(kind)) != 0 &&
            bqsr_read_filter_passes(batch, record, header, kind))
            return false;
    }
    return true;
}

std::pair<std::size_t, std::size_t> bqsr_soft_clip_bounds(
    const fastgatk::io::ReadBatch& batch, std::size_t record) {
    if (record >= batch.records() || batch.cigar_offsets.size() != batch.records() + 1)
        return {0, batch.offsets[record + 1] - batch.offsets[record]};
    const auto begin = batch.cigar_offsets[record];
    const auto end = batch.cigar_offsets[record + 1];
    std::size_t left = 0, right = 0;
    bool seen_alignment = false;
    for (std::size_t index = begin; index < end; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) break;
        if (!seen_alignment && operation.code == fastgatk::io::CigarOpCode::SoftClip)
            left += operation.length;
        else if (operation.consumes_read() || operation.consumes_reference())
            seen_alignment = true;
    }
    seen_alignment = false;
    for (std::size_t index = end; index > begin; ) {
        --index;
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) break;
        if (!seen_alignment && operation.code == fastgatk::io::CigarOpCode::SoftClip)
            right += operation.length;
        else if (operation.consumes_read() || operation.consumes_reference())
            seen_alignment = true;
    }
    const auto length = batch.offsets[record + 1] - batch.offsets[record];
    if (left + right > length) return {0, length};
    return {left, length - right};
}

char complement_base(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

std::string read_context_from(std::span<const std::uint8_t> bases,
                              std::span<const std::uint8_t> qualities,
                              std::uint16_t flags, std::size_t offset,
                              std::size_t context_size = 2) {
    if (offset >= bases.size() || qualities.size() != bases.size()) return {};
    const bool reverse = (flags & 0x10U) != 0;
    const auto stranded_offset = reverse ? bases.size() - offset - 1 : offset;
    if (stranded_offset + 1 < context_size) return {};
    std::size_t left_clip = 0, right_clip = bases.size();
    while (left_clip < right_clip && qualities[left_clip] <= 2) ++left_clip;
    while (right_clip > left_clip && qualities[right_clip - 1] <= 2) --right_clip;
    std::string context;
    context.reserve(context_size);
    for (std::size_t position = stranded_offset + 1 - context_size;
         position <= stranded_offset; ++position) {
        const auto original = reverse ? bases.size() - position - 1 : position;
        char base = 'N';
        // GATK's ContextCovariate only rewrites low-quality *tails* to N.
        // Low-quality bases in the interior remain their observed base and
        // are still valid context symbols; checking qualities here would
        // incorrectly erase those contexts.
        if (original >= left_clip && original < right_clip) {
            base = static_cast<char>(bases[original]);
            if (reverse) base = complement_base(base);
            else if (base >= 'a' && base <= 'z') base = static_cast<char>(base - 'a' + 'A');
            if (encode_base(static_cast<std::uint8_t>(base)) > 3) base = 'N';
        }
        if (base == 'N') return {};
        context.push_back(base);
    }
    return context;
}

std::int32_t read_cycle(std::uint16_t flags, std::size_t offset, std::size_t length) {
    const bool reverse = (flags & 0x10U) != 0;
    const bool second = (flags & 0x80U) != 0;
    const int order_factor = second ? -1 : 1;
    if (reverse) return static_cast<std::int32_t>(length * order_factor - offset * order_factor);
    return static_cast<std::int32_t>(order_factor + offset * order_factor);
}

std::int32_t checked_cycle(std::uint16_t flags, std::size_t offset, std::size_t length,
                           int maximum_cycle_value) {
    const auto cycle = read_cycle(flags, offset, length);
    if (std::abs(static_cast<std::int64_t>(cycle)) > maximum_cycle_value)
        throw std::runtime_error(
            "BAD_INPUT: The maximum allowed value for the cycle is " +
            std::to_string(maximum_cycle_value) + ", but a larger cycle (" +
            std::to_string(std::abs(static_cast<std::int64_t>(cycle))) +
            ") was detected. Please use the --maximum-cycle-value argument "
            "(when creating the recalibration table in BaseRecalibrator) to increase this value");
    return cycle;
}

// BaseRecalibrationEngine computes the substitution error array for every
// read base, including bases consumed by an insertion.  Those insertion
// bases have an isSNP value of zero (there is no reference base to compare),
// but they still contribute to the quality/read-group/cycle covariates.  The
// generic ReadProjection API intentionally treats insertions as non-projected
// for callers that build a reference pileup, so BQSR needs this local
// projection that preserves the insertion operation and its left anchor.
bool bqsr_project_read_offset(const fastgatk::io::ReadBatch& batch,
                              std::size_t record, std::size_t read_offset,
                              fastgatk::io::ReadProjection& projection) noexcept {
    if (fastgatk::io::project_read_offset(batch, record, read_offset, projection))
        return true;
    if (record >= batch.records() || batch.positions.size() <= record ||
        batch.positions[record] < 0 || !batch.cigar_record_layout_valid(record) ||
        record + 1 >= batch.offsets.size() ||
        read_offset >= batch.offsets[record + 1] - batch.offsets[record]) return false;
    std::size_t read_cursor = 0;
    std::int64_t reference_position = batch.positions[record];
    for (std::size_t index = batch.cigar_offsets[record];
         index < batch.cigar_offsets[record + 1]; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) return false;
        if (operation.code == fastgatk::io::CigarOpCode::Insertion &&
            read_offset >= read_cursor && read_offset - read_cursor < operation.length) {
            projection.read_offset = static_cast<std::uint32_t>(read_offset);
            projection.operation = operation.code;
            // The insertion is anchored immediately before the next
            // reference base in the forward CIGAR coordinate system.  The
            // anchor is used only for known-site masking; substitution error
            // remains zero for every inserted base.
            projection.reference_position = reference_position > 0
                ? reference_position - 1 : reference_position;
            return true;
        }
        if (operation.consumes_read()) read_cursor += operation.length;
        if (operation.consumes_reference()) reference_position += operation.length;
        if (read_cursor > batch.offsets[record + 1] - batch.offsets[record]) return false;
    }
    return false;
}

struct IndelMasks {
    std::vector<std::uint8_t> insertion;
    std::vector<std::uint8_t> deletion;
    std::vector<std::uint8_t> read_bases;
    std::size_t clipped_left = 0;
    std::size_t clipped_right = 0;
};

// BaseRecalibrator's optional I/D error models use the per-base BI/BD
// qualities when those SAM tags are present.  HTSlib keeps the tags as flat
// spans in ReadBatch, with an empty span representing a missing tag; in the
// latter case GATK's ReadUtils supplies the configured Q45 (or explicit
// --insertions/deletions-default-quality) fallback.  Keep the validation at
// the Host boundary so a truncated tag cannot silently shift event rows.
int indel_quality_at(const fastgatk::io::ReadBatch& batch, std::size_t record,
                     std::size_t read_offset, char event, int default_quality) {
    const auto& offsets = event == 'I' ? batch.insertion_quality_offsets
                                       : batch.deletion_quality_offsets;
    const auto& qualities = event == 'I' ? batch.insertion_qualities
                                         : batch.deletion_qualities;
    if (offsets.empty()) return default_quality;
    if (offsets.size() != batch.records() + 1 || record >= batch.records() ||
        offsets[record] > offsets[record + 1] || offsets[record + 1] > qualities.size())
        throw std::runtime_error("BAD_INPUT: malformed BI/BD indel quality layout");
    const auto begin = static_cast<std::size_t>(batch.offsets[record]);
    const auto end = static_cast<std::size_t>(batch.offsets[record + 1]);
    const auto tag_begin = static_cast<std::size_t>(offsets[record]);
    const auto tag_end = static_cast<std::size_t>(offsets[record + 1]);
    // HTSlib preserves one empty offset span for a missing optional tag so
    // BI/BD records stay aligned with the read batch.  An empty span is the
    // GATK ReadUtils fallback case, not a malformed tag; only a non-empty
    // span must match the complete read length.
    if (tag_end == tag_begin) return default_quality;
    if (tag_end - tag_begin != end - begin)
        throw std::runtime_error(std::string("BAD_INPUT: ") + (event == 'I' ? "BI" : "BD") +
                                 " indel quality length does not match read length");
    if (read_offset >= end - begin)
        throw std::runtime_error("BAD_INPUT: indel quality offset exceeds read length");
    return std::min<int>(qualities[tag_begin + read_offset], 93);
}

IndelMasks indel_masks_at(const fastgatk::io::ReadBatch& batch, std::size_t record) {
    const auto begin = batch.offsets[record];
    const auto end = batch.offsets[record + 1];
    IndelMasks masks{std::vector<std::uint8_t>(end - begin, 0),
                     std::vector<std::uint8_t>(end - begin, 0),
                     std::vector<std::uint8_t>(end - begin, 0), 0, 0};
    if (batch.cigar_offsets.size() != batch.records() + 1) return masks;
    bool seen_alignment = false;
    for (std::size_t index = batch.cigar_offsets[record];
         index < batch.cigar_offsets[record + 1]; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) break;
        if (!seen_alignment && operation.code == fastgatk::io::CigarOpCode::SoftClip) {
            masks.clipped_left += operation.length;
        } else {
            seen_alignment = true;
        }
    }
    seen_alignment = false;
    for (std::size_t index = batch.cigar_offsets[record + 1];
         index > batch.cigar_offsets[record]; ) {
        --index;
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) break;
        if (!seen_alignment && operation.code == fastgatk::io::CigarOpCode::SoftClip) {
            masks.clipped_right += operation.length;
        } else {
            seen_alignment = true;
        }
    }
    const auto flags = record < batch.flags.size() ? batch.flags[record] : 0U;
    const bool reverse = (flags & 0x10U) != 0;
    std::size_t read_position = 0;
    const auto mark = [&](std::vector<std::uint8_t>& mask, std::int64_t index) {
        if (index >= 0 && static_cast<std::size_t>(index) < mask.size())
            mask[static_cast<std::size_t>(index)] = 1;
    };
    for (std::size_t index = batch.cigar_offsets[record];
         index < batch.cigar_offsets[record + 1]; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
        if (!operation.valid()) break;
        switch (operation.code) {
            case fastgatk::io::CigarOpCode::Match:
            case fastgatk::io::CigarOpCode::SequenceMatch:
            case fastgatk::io::CigarOpCode::SequenceMismatch:
                for (std::size_t base = 0; base < operation.length &&
                     read_position + base < masks.read_bases.size(); ++base)
                    masks.read_bases[read_position + base] = 1;
                read_position += operation.length;
                break;
            case fastgatk::io::CigarOpCode::SoftClip:
                read_position += operation.length;
                break;
            case fastgatk::io::CigarOpCode::Insertion:
                for (std::size_t base = 0; base < operation.length &&
                     read_position + base < masks.read_bases.size(); ++base)
                    masks.read_bases[read_position + base] = 1;
                if (!reverse) mark(masks.insertion, static_cast<std::int64_t>(read_position) - 1);
                read_position += operation.length;
                if (reverse) mark(masks.insertion, static_cast<std::int64_t>(read_position));
                break;
            case fastgatk::io::CigarOpCode::Deletion:
                mark(masks.deletion, reverse ? static_cast<std::int64_t>(read_position)
                                             : static_cast<std::int64_t>(read_position) - 1);
                break;
            case fastgatk::io::CigarOpCode::ReferenceSkip:
            case fastgatk::io::CigarOpCode::HardClip:
            case fastgatk::io::CigarOpCode::Padding:
                break;
        }
    }
    return masks;
}

CovariateTable read_covariates(const std::string& path) {
    CovariateTable table;
    std::ifstream input(path);
    if (!input) return table;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        CovariateKey key;
        QualityBin bin;
        double empirical = 0.0;
        int delta = 0;
        if (!(fields >> key.read_group >> key.cycle >> key.context >> key.quality >>
              bin.count >> bin.mismatches >> empirical >> delta))
            throw std::runtime_error("BAD_INPUT: malformed BQSR covariate row: " + path);
        std::string event;
        if (fields >> event && !event.empty()) key.event = event.front();
        if (key.read_group == "-") key.read_group.clear();
        if (key.context == "-") key.context.clear();
        // GatherBQSRReports may legitimately produce duplicate descriptor rows
        // when reports were created from independent shards.  Merge integer
        // counters instead of silently discarding all but the first row.
        add_quality_bin(table[key], bin);
    }
    return table;
}

CovariateModel build_covariate_model(const CovariateTable& table) {
    CovariateModel model;
    std::map<std::string, double> expected_errors;
    std::map<std::string, std::uint64_t> observations;
    for (const auto& [key, bin] : table) {
        // ApplyBQSR recalibrates base qualities from the substitution model.
        // I/D rows are retained in the report for GATK compatibility but must
        // not be folded into the M posterior.
        if (key.event != 'M') continue;
        add_quality_bin(model.read_group[key.read_group], bin);
        expected_errors[key.read_group] += static_cast<double>(bin.count) *
            std::pow(10.0, -static_cast<double>(key.quality) / 10.0);
        observations[key.read_group] += bin.count;
        add_quality_bin(model.read_group_quality[{key.read_group, key.quality}], bin);
        add_quality_bin(model.context_quality[{key.read_group, key.context, key.quality}], bin);
        add_quality_bin(model.cycle_quality[{key.read_group, key.cycle, key.quality}], bin);
    }
    for (const auto& [read_group, count] : observations) {
        const auto expected = expected_errors[read_group];
        model.read_group_reported_quality[read_group] =
            count == 0 || expected <= 0.0 ? 0.0 : -10.0 * std::log10(expected / count);
    }
    return model;
}

struct BqsrReportModel {
    // ContextCovariate is parameterized by the report that produced the
    // recalibration tables.  ApplyBQSR must reconstruct that same k-mer
    // width; its own CLI intentionally has no context-size override.
    int mismatches_context_size = 2;
    std::array<int, 94> quality_deltas{};
    std::array<int, 94> quantized_qualities{};
    std::array<std::uint64_t, 94> empirical_quality_counts{};
    std::map<std::string, QualityBin> read_group;
    std::map<std::string, double> read_group_reported_quality;
    std::map<std::pair<std::string, int>, QualityBin> read_group_quality;
    std::map<std::tuple<std::string, std::string, int>, QualityBin> context_quality;
    std::map<std::tuple<std::string, int, int>, QualityBin> cycle_quality;
};

int bayesian_empirical_quality(const QualityBin& bin, double prior_quality);
int covariate_delta(const QualityBin& bin, int quality);

int model_delta(const CovariateTable& exact,
                const CovariateModel& model,
                const BqsrReportModel& report,
                const CovariateKey& key,
                double global_qscore_prior = -1.0) {
    // Match BQSRReadTransformer's hierarchical Bayesian estimate:
    //   read-group -> (read-group, reported quality) -> special covariates.
    // A special covariate contributes a delta relative to the same posterior
    // quality-score datum; the deltas are not applied iteratively.
    const auto read_group_bin = [&]() -> const QualityBin* {
        if (const auto found = model.read_group.find(key.read_group);
            found != model.read_group.end()) return &found->second;
        if (const auto found = report.read_group.find(key.read_group);
            found != report.read_group.end()) return &found->second;
        return nullptr;
    }();
    const auto quality_bin = [&]() -> const QualityBin* {
        if (const auto found = model.read_group_quality.find({key.read_group, key.quality});
            found != model.read_group_quality.end()) return &found->second;
        if (const auto found = report.read_group_quality.find({key.read_group, key.quality});
            found != report.read_group_quality.end()) return &found->second;
        return nullptr;
    }();

    double read_group_prior = static_cast<double>(key.quality);
    if (global_qscore_prior > 0.0) {
        read_group_prior = global_qscore_prior;
    } else if (const auto found = model.read_group_reported_quality.find(key.read_group);
               found != model.read_group_reported_quality.end()) {
        read_group_prior = found->second;
    } else if (const auto found = report.read_group_reported_quality.find(key.read_group);
               found != report.read_group_reported_quality.end()) {
        read_group_prior = found->second;
    }
    double posterior = read_group_bin == nullptr
        ? read_group_prior
        : static_cast<double>(bayesian_empirical_quality(*read_group_bin, read_group_prior));
    if (quality_bin != nullptr) {
        posterior = bayesian_empirical_quality(*quality_bin, posterior);
    } else {
        posterior += report.quality_deltas[static_cast<std::size_t>(
            std::clamp(key.quality, 0, 93))];
    }
    const auto special_delta = [&](const auto& model_table, const auto& report_table,
                                   const auto& lookup) {
        if (const auto found = model_table.find(lookup); found != model_table.end())
            return static_cast<double>(bayesian_empirical_quality(found->second, posterior)) - posterior;
        if (const auto found = report_table.find(lookup); found != report_table.end())
            return static_cast<double>(bayesian_empirical_quality(found->second, posterior)) - posterior;
        return 0.0;
    };
    const auto corrected = posterior
        + special_delta(model.context_quality, report.context_quality,
                        std::make_tuple(key.read_group, key.context, key.quality))
        + special_delta(model.cycle_quality, report.cycle_quality,
                        std::make_tuple(key.read_group, key.cycle, key.quality));
    if (read_group_bin == nullptr && quality_bin == nullptr &&
        model.context_quality.find({key.read_group, key.context, key.quality}) == model.context_quality.end() &&
        model.cycle_quality.find({key.read_group, key.cycle, key.quality}) == model.cycle_quality.end() &&
        report.context_quality.find({key.read_group, key.context, key.quality}) == report.context_quality.end() &&
        report.cycle_quality.find({key.read_group, key.cycle, key.quality}) == report.cycle_quality.end() &&
        exact.empty()) {
        return report.quality_deltas[static_cast<std::size_t>(std::clamp(key.quality, 0, 93))];
    }
    return static_cast<int>(std::lround(std::clamp(
        corrected - static_cast<double>(key.quality), -20.0, 20.0)));
}

std::string normalize_read_group(std::string read_group) {
    if (read_group == "UNKNOWN" || read_group == "-") read_group.clear();
    return read_group;
}

BqsrReportModel read_report(const std::string& path) {
    BqsrReportModel model;
    for (std::size_t quality = 0; quality < model.quantized_qualities.size(); ++quality)
        model.quantized_qualities[quality] = static_cast<int>(quality);
    std::array<QualityBin, 94> global_bins{};
    std::map<std::string, double> expected_errors;
    std::map<std::string, std::uint64_t> reported_observations;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open recalibration report: " + path);
    std::string line;
    enum class Table { None, Arguments, Quantized, RecalTable0, RecalTable1, RecalTable2 } table = Table::None;
    bool saw_gatk_table = false;
    while (std::getline(input, line)) {
        if (line.rfind("#:GATKTable:Arguments:", 0) == 0) {
            table = Table::Arguments;
            continue;
        }
        if (line.rfind("#:GATKTable:RecalTable0:", 0) == 0) {
            table = Table::RecalTable0;
            saw_gatk_table = true;
            continue;
        }
        if (line.rfind("#:GATKTable:RecalTable1:", 0) == 0) {
            table = Table::RecalTable1;
            saw_gatk_table = true;
            continue;
        }
        if (line.rfind("#:GATKTable:Quantized:", 0) == 0) {
            table = Table::Quantized;
            continue;
        }
        if (line.rfind("#:GATKTable:RecalTable2:", 0) == 0) {
            table = Table::RecalTable2;
            saw_gatk_table = true;
            continue;
        }
        if (line.rfind("#:GATKTable:", 0) == 0 &&
            line.rfind("#:GATKTable:Arguments:", 0) != 0 &&
            line.rfind("#:GATKTable:RecalTable0:", 0) != 0 &&
            line.rfind("#:GATKTable:RecalTable1:", 0) != 0 &&
            line.rfind("#:GATKTable:RecalTable2:", 0) != 0) {
            table = Table::None;
            continue;
        }
        if (table == Table::Arguments) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            std::string argument, value;
            if (fields >> argument >> value && argument == "mismatches_context_size") {
                try {
                    std::size_t consumed = 0;
                    const auto parsed = std::stoi(value, &consumed);
                    if (consumed != value.size() || parsed < 1 || parsed > 13)
                        throw std::invalid_argument("outside ContextCovariate domain");
                    model.mismatches_context_size = parsed;
                } catch (const std::exception&) {
                    throw std::runtime_error(
                        "BAD_INPUT: invalid mismatches_context_size in recalibration report");
                }
            }
            continue;
        }
        if (table == Table::RecalTable0) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            std::string read_group, event;
            double empirical = 0.0, estimated = 0.0, errors = 0.0;
            std::uint64_t observations = 0;
            if (fields >> read_group >> event >> empirical >> estimated >> observations >> errors &&
                event == "M") {
                const auto normalized_group = normalize_read_group(read_group);
                model.read_group_reported_quality[normalized_group] = estimated;
            }
            continue;
        }
        if (table == Table::RecalTable1) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            std::string read_group, event;
            int quality = 0;
            double empirical = 0.0;
            std::uint64_t observations = 0;
            double errors = 0.0;
            if (fields >> read_group >> quality >> event >> empirical >> observations >> errors &&
                event == "M" && quality >= 0 && quality < static_cast<int>(global_bins.size())) {
                const QualityBin bin{observations, static_cast<std::uint64_t>(std::llround(std::max(errors, 0.0)))};
                add_quality_bin(global_bins[static_cast<std::size_t>(quality)], bin);
                const auto normalized_group = normalize_read_group(read_group);
                add_quality_bin(model.read_group[normalized_group], bin);
                add_quality_bin(model.read_group_quality[{normalized_group, quality}], bin);
                expected_errors[normalized_group] += static_cast<double>(observations) *
                    std::pow(10.0, -static_cast<double>(quality) / 10.0);
                reported_observations[normalized_group] += observations;
            }
            continue;
        }
        if (table == Table::Quantized) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            int quality = 0, quantized = 0;
            std::uint64_t count = 0;
            if (fields >> quality >> count >> quantized &&
                quality >= 0 && quality < static_cast<int>(model.quantized_qualities.size()) &&
                quantized >= 0 && quantized <= 93) {
                model.quantized_qualities[static_cast<std::size_t>(quality)] = quantized;
                model.empirical_quality_counts[static_cast<std::size_t>(quality)] = count;
            }
            continue;
        }
        if (table == Table::RecalTable2) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream fields(line);
            std::string read_group, covariate_value, covariate_name, event;
            int quality = 0;
            double empirical = 0.0;
            std::uint64_t observations = 0;
            double errors = 0.0;
            if (!(fields >> read_group >> quality >> covariate_value >> covariate_name >> event >>
                  empirical >> observations >> errors) || event != "M" || quality < 0 || quality >= 94)
                continue;
            const QualityBin bin{observations, static_cast<std::uint64_t>(std::llround(std::max(errors, 0.0)))};
            const auto normalized_group = normalize_read_group(read_group);
            if (covariate_value == "-") covariate_value.clear();
            if (covariate_name == "Context") {
                add_quality_bin(model.context_quality[{normalized_group, covariate_value, quality}], bin);
            } else if (covariate_name == "Cycle") {
                try {
                    const auto cycle = std::stoi(covariate_value);
                    add_quality_bin(model.cycle_quality[{normalized_group, cycle, quality}], bin);
                } catch (const std::exception&) {
                    // Ignore malformed optional covariate rows; RecalTable1 remains usable.
                }
            }
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        int quality = 0, delta = 0;
        std::uint64_t count = 0, mismatches = 0;
        double empirical = 0.0;
        if (fields >> quality >> count >> mismatches >> empirical >> delta) {
            if (quality >= 0 && quality < static_cast<int>(model.quality_deltas.size()))
                model.quality_deltas[static_cast<std::size_t>(quality)] = delta;
        }
    }
    for (const auto& [read_group, count] : reported_observations) {
        if (model.read_group_reported_quality.find(read_group) != model.read_group_reported_quality.end())
            continue;
        const auto expected = expected_errors[read_group];
        model.read_group_reported_quality[read_group] =
            count == 0 || expected <= 0.0 ? 0.0 : -10.0 * std::log10(expected / count);
    }
    if (saw_gatk_table) {
        for (std::size_t quality = 0; quality < global_bins.size(); ++quality)
            model.quality_deltas[quality] = covariate_delta(global_bins[quality], static_cast<int>(quality));
    }
    return model;
}

// Mirror BQSRReadTransformer.constructStaticQuantizedMapping from the pinned
// GATK release.  Static bins are deliberately separate from the report's
// dynamic Quantized table: Java first applies the report quantized value and
// then maps that value into the requested static bins.  Quality values below
// MIN_USABLE_Q_SCORE (Q6) are reserved/special and remain unchanged.
std::vector<int> construct_static_quantized_mapping(
    const std::vector<int>& requested, bool round_down) {
    constexpr int kMinUsableQuality = 6;
    constexpr int kMaxQuality = 254;
    std::vector<int> mapping(static_cast<std::size_t>(kMaxQuality + 1));
    for (int quality = 0; quality < kMaxQuality + 1; ++quality)
        mapping[static_cast<std::size_t>(quality)] = quality;
    if (requested.empty()) return mapping;

    auto bins = requested;
    std::sort(bins.begin(), bins.end());
    bins.erase(std::unique(bins.begin(), bins.end()), bins.end());
    if (bins.size() == 1) {
        for (int quality = kMinUsableQuality; quality < kMaxQuality + 1; ++quality)
            mapping[static_cast<std::size_t>(quality)] = bins.front();
        return mapping;
    }

    auto quality_probability = [](int quality) {
        return 1.0 - std::pow(10.0, -static_cast<double>(quality) / 10.0);
    };
    int previous_quality = kMinUsableQuality;
    double previous_probability = quality_probability(previous_quality);
    for (const int next_quality : bins) {
        const double next_probability = quality_probability(next_quality);
        for (int quality = previous_quality; quality < next_quality; ++quality) {
            if (round_down) {
                mapping[static_cast<std::size_t>(quality)] = previous_quality;
            } else {
                const double quality_probability_value = quality_probability(quality);
                mapping[static_cast<std::size_t>(quality)] =
                    quality_probability_value - previous_probability >
                            next_probability - quality_probability_value
                        ? next_quality : previous_quality;
            }
        }
        previous_quality = next_quality;
        previous_probability = next_probability;
    }
    for (int quality = previous_quality; quality < kMaxQuality + 1; ++quality)
        mapping[static_cast<std::size_t>(quality)] = previous_quality;
    return mapping;
}

// Reproduce the pinned GATK QualQuantizer greedy adjacent-interval merge.
// The report's Quantized table is an empirical-quality histogram, so positive
// --quantize-quals values can be recalculated without rescanning the BAM.
std::array<int, 94> construct_dynamic_quantized_mapping(
    const BqsrReportModel& report, int requested_levels) {
    constexpr int kMinInterestingQuality = 6;
    constexpr int kMaxQuality = 93;
    struct Interval {
        int begin = 0;
        int end = 0;
        std::uint64_t observations = 0;
        std::uint64_t errors = 0;
        bool fixed = true;
        int fixed_quality = 0;
    };
    std::array<int, 94> identity{};
    for (int quality = 0; quality <= kMaxQuality; ++quality)
        identity[static_cast<std::size_t>(quality)] = quality;
    if (requested_levels <= 0 || requested_levels >= 94) return identity;

    auto error_probability = [](int quality) {
        return std::pow(10.0, -static_cast<double>(quality) / 10.0);
    };
    auto interval_error_rate = [](const Interval& interval) {
        if (interval.fixed)
            return std::pow(10.0, -static_cast<double>(interval.fixed_quality) / 10.0);
        if (interval.observations == 0) return 0.0;
        return (static_cast<double>(interval.errors) + 1.0) /
               (static_cast<double>(interval.observations) + 1.0);
    };
    auto interval_quality = [&](const Interval& interval) {
        const auto rate = interval_error_rate(interval);
        if (rate <= 0.0) return kMaxQuality;
        const auto rounded = static_cast<int>(std::floor(-10.0 * std::log10(rate) + 0.5));
        return std::clamp(rounded, 1, kMaxQuality);
    };
    auto penalty = [&](const Interval& interval, double global_rate) {
        if (global_rate == 0.0) return 0.0;
        double total = 0.0;
        for (int quality = interval.begin; quality <= interval.end; ++quality) {
            if (quality <= kMinInterestingQuality) continue;
            const auto observations = report.empirical_quality_counts[static_cast<std::size_t>(quality)];
            if (observations == 0) continue;
            const auto leaf_rate = error_probability(quality);
            total += std::abs(std::log10(leaf_rate) - std::log10(global_rate)) *
                     static_cast<double>(observations);
        }
        return total;
    };

    std::vector<Interval> intervals;
    intervals.reserve(94);
    for (int quality = 0; quality <= kMaxQuality; ++quality) {
        const auto observations = report.empirical_quality_counts[static_cast<std::size_t>(quality)];
        const auto expected_errors = static_cast<double>(observations) *
                                     error_probability(quality);
        intervals.push_back(Interval{quality, quality, observations,
                                     static_cast<std::uint64_t>(std::floor(expected_errors)),
                                     true, quality});
    }
    while (static_cast<int>(intervals.size()) > requested_levels) {
        std::size_t best = 0;
        double best_penalty = std::numeric_limits<double>::infinity();
        for (std::size_t index = 0; index + 1 < intervals.size(); ++index) {
            const auto& left = intervals[index];
            const auto& right = intervals[index + 1];
            Interval merged{left.begin, right.end,
                             left.observations + right.observations,
                             left.errors + right.errors, false, 0};
            const auto candidate_penalty = penalty(merged, interval_error_rate(merged));
            // Java keeps the first candidate on an exact tie.
            if (candidate_penalty < best_penalty) {
                best_penalty = candidate_penalty;
                best = index;
            }
        }
        const auto& left = intervals[best];
        const auto& right = intervals[best + 1];
        intervals[best] = Interval{left.begin, right.end,
                                   left.observations + right.observations,
                                   left.errors + right.errors, false, 0};
        intervals.erase(intervals.begin() + static_cast<std::ptrdiff_t>(best + 1));
    }
    std::array<int, 94> mapping{};
    for (const auto& interval : intervals) {
        const auto quality = interval.fixed ? interval.fixed_quality : interval_quality(interval);
        for (int value = interval.begin; value <= interval.end; ++value)
            mapping[static_cast<std::size_t>(value)] = quality;
    }
    return mapping;
}

int bayesian_empirical_quality(const QualityBin& bin, double prior_quality) {
    // RecalDatum uses one pseudo-error and two pseudo-observations, then a
    // Gaussian prior N(prior, 0.5) and a binomial likelihood.  The search is
    // over GATK's reasonable quality range [0, 60], with the first maximum
    // retained to match MathUtils.maxElementIndex tie behavior.
    long double observations = static_cast<long double>(bin.count) + 2.0L;
    long double errors = static_cast<long double>(std::llround(
        static_cast<long double>(bin.mismatches))) + 1.0L;
    if (observations > static_cast<long double>(std::numeric_limits<std::int32_t>::max() - 1)) {
        const auto fraction = static_cast<long double>(std::numeric_limits<std::int32_t>::max() - 1) /
                              observations;
        errors = std::llround(errors * fraction);
        observations = static_cast<long double>(std::numeric_limits<std::int32_t>::max() - 1);
    }
    errors = std::clamp(errors, 0.0L, observations);
    int best_quality = 0;
    long double best_posterior = -std::numeric_limits<long double>::max();
    for (int quality = 0; quality <= 60; ++quality) {
        // GATK's RecalDatum.getLogPrior deliberately truncates the
        // quality/prior difference to an int before taking abs (rather than
        // using the continuous distance).  This creates equal prior mass for
        // adjacent integer bins around a fractional EstimatedQReported and
        // affects the deterministic MathUtils.maxElementIndex tie-break.
        const auto difference = std::min(
            std::abs(static_cast<int>(static_cast<long double>(quality) - prior_quality)), 40);
        const auto log_prior = -0.5L * (static_cast<long double>(difference) / 0.5L) *
                               (static_cast<long double>(difference) / 0.5L);
        const auto probability = std::pow(10.0L, -static_cast<long double>(quality) / 10.0L);
        long double log_likelihood = 0.0L;
        if (probability <= 0.0L || probability >= 1.0L) {
            if (probability >= 1.0L && errors < observations) log_likelihood =
                -std::numeric_limits<long double>::max();
            else if (probability <= 0.0L && errors > 0.0L) log_likelihood =
                -std::numeric_limits<long double>::max();
        } else {
            log_likelihood = std::lgammal(observations + 1.0L)
                - std::lgammal(errors + 1.0L)
                - std::lgammal(observations - errors + 1.0L)
                + errors * std::log(probability)
                + (observations - errors) * std::log1pl(-probability);
        }
        const auto posterior = log_prior + log_likelihood;
        if (posterior > best_posterior) {
            best_posterior = posterior;
            best_quality = quality;
        }
    }
    return best_quality;
}

double empirical_quality(const QualityBin& bin, int default_quality) {
    return static_cast<double>(bayesian_empirical_quality(bin, default_quality));
}

std::string display_read_group(const std::string& read_group) {
    return read_group.empty() ? "UNKNOWN" : read_group;
}

struct ReadGroupSummary {
    std::array<QualityBin, 94> bins{};
};

std::map<std::string, ReadGroupSummary> summarize_read_groups(
    const CovariateTable& covariates, const std::array<QualityBin, 94>& bins) {
    std::map<std::string, ReadGroupSummary> groups;
    for (const auto& [key, value] : covariates) {
        if (key.event != 'M') continue;
        auto& destination = groups[key.read_group];
        if (key.quality >= 0 && key.quality < static_cast<int>(destination.bins.size())) {
            destination.bins[static_cast<std::size_t>(key.quality)].count += value.count;
            destination.bins[static_cast<std::size_t>(key.quality)].mismatches += value.mismatches;
        }
    }
    if (groups.empty()) groups.emplace(std::string{}, ReadGroupSummary{});
    std::uint64_t grouped_observations = 0;
    for (const auto& [_, group] : groups)
        for (const auto& bin : group.bins) grouped_observations += bin.count;
    if (grouped_observations == 0) {
        auto& destination = groups.begin()->second;
        destination.bins = bins;
    }
    return groups;
}

struct QuantizedInterval {
    int begin = 0;
    int end = 0;
    std::uint64_t observations = 0;
    std::uint64_t errors = 0;
};

double quantizer_penalty(const QuantizedInterval& interval,
                         const std::array<std::uint64_t, 94>& histogram) {
    if (interval.observations == 0) return 0.0;
    const auto global_error = static_cast<long double>(interval.errors + 1.0L) /
                              static_cast<long double>(interval.observations + 1.0L);
    if (global_error <= 0.0L) return 0.0;
    double penalty = 0.0;
    for (int quality = interval.begin; quality <= interval.end; ++quality) {
        if (quality <= 6 || histogram[static_cast<std::size_t>(quality)] == 0) continue;
        const auto leaf_error = std::pow(10.0L, -static_cast<long double>(quality) / 10.0L);
        penalty += std::abs(std::log10(static_cast<double>(leaf_error)) -
                            std::log10(static_cast<double>(global_error))) *
                   static_cast<double>(histogram[static_cast<std::size_t>(quality)]);
    }
    return penalty;
}

std::array<int, 94> build_quantized_map(const std::array<std::uint64_t, 94>& histogram,
                                        int levels = 16) {
    std::vector<QuantizedInterval> intervals;
    intervals.reserve(histogram.size());
    for (int quality = 0; quality < static_cast<int>(histogram.size()); ++quality) {
        const auto observations = histogram[static_cast<std::size_t>(quality)];
        const auto error_probability = std::pow(10.0L, -static_cast<long double>(quality) / 10.0L);
        intervals.push_back(QuantizedInterval{
            quality, quality, observations,
            static_cast<std::uint64_t>(std::floor(static_cast<long double>(observations) *
                                                   error_probability))});
    }
    levels = std::max(1, std::min(levels, static_cast<int>(intervals.size())));
    while (static_cast<int>(intervals.size()) > levels) {
        std::size_t best = 0;
        double best_penalty = std::numeric_limits<double>::infinity();
        for (std::size_t index = 0; index + 1 < intervals.size(); ++index) {
            QuantizedInterval merged{
                intervals[index].begin, intervals[index + 1].end,
                intervals[index].observations + intervals[index + 1].observations,
                intervals[index].errors + intervals[index + 1].errors};
            const auto penalty = quantizer_penalty(merged, histogram);
            if (penalty < best_penalty) {
                best_penalty = penalty;
                best = index;
            }
        }
        intervals[best].end = intervals[best + 1].end;
        intervals[best].observations += intervals[best + 1].observations;
        intervals[best].errors += intervals[best + 1].errors;
        intervals.erase(intervals.begin() + static_cast<std::ptrdiff_t>(best + 1));
    }
    std::array<int, 94> mapping{};
    for (const auto& interval : intervals) {
        int quantized = interval.begin;
        if (interval.begin != interval.end) {
            // QualQuantizer treats an all-zero interval as error rate 0 and
            // therefore emits the SAM maximum Q93.  Applying the usual
            // pseudo-count (0+1)/(0+1)=1 here would incorrectly map an
            // unobserved high-quality tail to Q0/Q1.
            if (interval.observations == 0) {
                quantized = 93;
            } else {
                const auto error_rate = static_cast<double>(interval.errors + 1.0L) /
                                        static_cast<double>(interval.observations + 1.0L);
                quantized = static_cast<int>(std::lround(-10.0 * std::log10(
                    std::clamp(error_rate, 1.0e-300, 1.0))));
                quantized = std::clamp(quantized, 1, 93);
            }
        }
        for (int quality = interval.begin; quality <= interval.end; ++quality)
            mapping[static_cast<std::size_t>(quality)] = quantized;
    }
    // QualQuantizer includes the complete [0,93] histogram, including
    // zero-observation bins.  Do not special-case the high-quality tail:
    // those bins participate in the same greedy contiguous merge and may map
    // to a lower representative when all observations stop below Q93.
    return mapping;
}

std::vector<std::string> split_report_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos
            ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

std::vector<std::string> split_report_colons(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find(':', begin);
        fields.push_back(line.substr(begin, end == std::string::npos
            ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

bool report_numeric_column(const std::string& table, std::size_t column) {
    if (table == "Quantized") return true;
    if (table == "RecalTable0") return column >= 2;
    if (table == "RecalTable1") return column == 1 || column >= 3;
    if (table == "RecalTable2") return column == 1 || column >= 5;
    return false;
}

std::string pad_report_cell(const std::string& value, std::size_t width, bool right) {
    if (value.size() >= width) return value;
    const auto padding = std::string(width - value.size(), ' ');
    return right ? padding + value : value + padding;
}

// GATKReportTable's v1.1 reader is fixed-width: it determines column starts
// from the header line and slices every row at those offsets.  Tabs are valid
// whitespace but are not sufficient when a value is wider than its column
// name (eg. Arguments/binary_tag_name).  Native reports are initially built
// with simple tab delimiters; normalize the completed report once, preserving
// the values while emitting the same padded shape expected by HTSJDK.
void rewrite_gatk_report_fixed_width(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot reopen recalibration report: " + path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) lines.push_back(line);
    input.close();

    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot normalize recalibration report: " + path);
    for (std::size_t index = 0; index < lines.size();) {
        const auto& definition = lines[index];
        // A table definition starts with a numeric column count; table-name
        // and description lines also carry the GATKTable prefix but do not.
        if (definition.rfind("#:GATKTable:", 0) != 0 ||
            definition.size() <= 12 || definition[12] < '0' || definition[12] > '9' ||
            index + 2 >= lines.size()) {
            output << lines[index++] << '\n';
            continue;
        }
        const auto definition_fields = split_report_colons(definition);
        if (definition_fields.size() < 4) {
            output << lines[index++] << '\n';
            continue;
        }
        const auto columns = static_cast<std::size_t>(std::stoul(definition_fields[2]));
        const auto rows = static_cast<std::size_t>(std::stoul(definition_fields[3]));
        const auto name_fields = split_report_colons(lines[index + 1]);
        const auto table = name_fields.size() > 2 ? name_fields[2] : std::string{};
        const auto header = split_report_fields(lines[index + 2]);
        if (header.size() != columns || index + 3 + rows > lines.size()) {
            output << lines[index++] << '\n';
            continue;
        }
        std::vector<std::vector<std::string>> data;
        data.reserve(rows);
        std::vector<std::size_t> widths(columns, 0);
        for (std::size_t column = 0; column < columns; ++column)
            widths[column] = header[column].size();
        for (std::size_t row = 0; row < rows; ++row) {
            auto fields = split_report_fields(lines[index + 3 + row]);
            if (fields.size() != columns) {
                // Keep malformed input untouched; the native writer should
                // never generate this branch, but avoiding a silent rewrite
                // makes future format extensions fail visibly in Java.
                output << lines[index++] << '\n';
                data.clear();
                break;
            }
            for (std::size_t column = 0; column < columns; ++column)
                widths[column] = std::max(widths[column], fields[column].size());
            data.push_back(std::move(fields));
        }
        if (data.size() != rows) continue;
        output << definition << '\n' << lines[index + 1] << '\n';
        for (std::size_t column = 0; column < columns; ++column) {
            if (column != 0) output << "  ";
            output << pad_report_cell(header[column], widths[column], false);
        }
        output << '\n';
        for (const auto& fields : data) {
            for (std::size_t column = 0; column < columns; ++column) {
                if (column != 0) output << "  ";
                output << pad_report_cell(fields[column], widths[column],
                                          report_numeric_column(table, column));
            }
            output << '\n';
        }
        output << '\n';
        index += 4 + rows;
    }
    if (!output) throw std::runtime_error("cannot finalize normalized recalibration report: " + path);
}

void write_gatk_report(const std::string& path,
                       const std::array<QualityBin, 94>& bins,
                       const CovariateTable& covariates,
                       const Options& options) {
    const auto groups = summarize_read_groups(covariates, bins);
    std::array<std::uint64_t, 94> empirical_histogram{};
    // Quantization is built from the quality-score table (read-group,
    // reported-quality), not from the already-collapsed global histogram.
    // This distinction matters for low-depth groups because each datum gets
    // its own Bayesian prior before empirical qualities are pooled.
    std::map<std::pair<std::string, int>, QualityBin> quality_table;
    std::map<std::tuple<std::string, int, char>, QualityBin> indel_quality_table;
    for (const auto& [key, bin] : covariates) {
        if (key.event == 'M') {
            add_quality_bin(quality_table[{key.read_group, key.quality}], bin);
        } else if (options.compute_indel_bqsr_tables &&
                   (key.event == 'I' || key.event == 'D')) {
            add_quality_bin(indel_quality_table[{key.read_group, key.quality, key.event}], bin);
        }
    }
    std::map<std::tuple<std::string, std::string, int>, QualityBin> context_table;
    std::map<std::tuple<std::string, int, int>, QualityBin> cycle_table;
    std::map<std::tuple<std::string, std::string, int, char>, QualityBin> indel_context_table;
    std::map<std::tuple<std::string, int, int, char>, QualityBin> indel_cycle_table;
    for (const auto& [key, bin] : covariates) {
        if (key.event != 'M') {
            if (!options.compute_indel_bqsr_tables) continue;
            if (!key.context.empty())
                add_quality_bin(indel_context_table[{key.read_group, key.context, key.quality, key.event}], bin);
            if (key.cycle != std::numeric_limits<std::int32_t>::min())
                add_quality_bin(indel_cycle_table[{key.read_group, key.cycle, key.quality, key.event}], bin);
            continue;
        }
        if (!key.context.empty())
            add_quality_bin(context_table[{key.read_group, key.context, key.quality}], bin);
        add_quality_bin(cycle_table[{key.read_group, key.cycle, key.quality}], bin);
    }
    for (const auto& [key, bin] : quality_table) {
        const auto empirical = bayesian_empirical_quality(bin, static_cast<double>(key.second));
        empirical_histogram[static_cast<std::size_t>(std::clamp(empirical, 0, 93))] += bin.count;
    }
    // GATK's quantizer includes the optional I/D quality-score tables in the
    // histogram.  Their reported quality is the per-read BI/BD quality when
    // present (and the configured default otherwise); using the configured
    // fallback for every row would silently distort reports carrying indel
    // quality tags.
    if (options.compute_indel_bqsr_tables) {
        for (const auto& [key, bin] : indel_quality_table) {
            const auto empirical = bayesian_empirical_quality(
                bin, static_cast<double>(std::get<1>(key)));
            empirical_histogram[static_cast<std::size_t>(std::clamp(
                empirical, 0, 93))] += bin.count;
        }
    }
    if (quality_table.empty()) {
        for (std::size_t quality = 0; quality < bins.size(); ++quality)
            empirical_histogram[quality] = bins[quality].count;
    }
    auto quantized_map = build_quantized_map(empirical_histogram, 16);
    if (options.compute_indel_bqsr_tables) {
        int last_substitution_quality = 0;
        for (const auto& [key, bin] : quality_table)
            if (bin.count != 0) last_substitution_quality = std::max(last_substitution_quality, key.second);
        int first_indel_quality = 94;
        for (const auto& [key, bin] : indel_quality_table) {
            if (bin.count != 0)
                first_indel_quality = std::min(first_indel_quality, std::get<1>(key));
        }
        // QualQuantizer leaves an unobserved high-quality gap immediately
        // below the first observed I/D quality in the unquantized tail.  The
        // substitution-only path has no such internal gap.
        for (int quality = last_substitution_quality + 1;
             quality < first_indel_quality && quality < static_cast<int>(quantized_map.size()); ++quality)
            if (empirical_histogram[static_cast<std::size_t>(quality)] == 0)
                quantized_map[static_cast<std::size_t>(quality)] = 93;
    }
    std::size_t table1_rows = 0;
    for (const auto& [_, group] : groups) {
        for (const auto& bin : group.bins) if (bin.count != 0) ++table1_rows;
    }
    std::map<std::pair<std::string, char>, QualityBin> indel_groups;
    std::map<std::tuple<std::string, int, char>, QualityBin> indel_quality_groups;
    for (const auto& [key, bin] : covariates) {
        if (key.event == 'M' || !options.compute_indel_bqsr_tables) continue;
        add_quality_bin(indel_groups[{key.read_group, key.event}], bin);
        add_quality_bin(indel_quality_groups[{key.read_group, key.quality, key.event}], bin);
    }
    for (const auto& [_, bin] : indel_quality_groups) if (bin.count != 0) ++table1_rows;
    const std::size_t table0_rows = groups.size() + indel_groups.size();
    const std::size_t table2_rows = context_table.size() + cycle_table.size() +
                                    indel_context_table.size() + indel_cycle_table.size();
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open recalibration report: " + path);
    output << "#:GATKReport.v1.1:5\n"
           << "#:GATKTable:2:17:%s:%s:;\n"
           << "#:GATKTable:Arguments:Recalibration argument collection values used in this run\n"
           << "Argument\tValue\n"
           << "binary_tag_name\tnull\n"
           << "covariate\tReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate\n"
           << "default_platform\tnull\n"
           << "deletions_default_quality\t" << options.deletions_default_quality << "\n"
           << "force_platform\tnull\n"
           << "indels_context_size\t" << options.indels_context_size << "\n"
           << "insertions_default_quality\t" << options.insertions_default_quality << "\n"
           << "low_quality_tail\t2\n"
           << "maximum_cycle_value\t" << options.maximum_cycle_value << "\n"
           << "mismatches_context_size\t" << options.mismatches_context_size << "\n"
           << "mismatches_default_quality\t-1\n"
           << "no_standard_covs\tfalse\n"
           << "quantizing_levels\t16\n"
           << "recalibration_report\tnull\n"
           << "run_without_dbsnp\tfalse\n"
           << "solid_nocall_strategy\tTHROW_EXCEPTION\n"
           << "solid_recal_mode\tSET_Q_ZERO\n\n"
           << "#:GATKTable:3:94:%d:%d:%d:;\n"
           << "#:GATKTable:Quantized:Quality quantization map\n"
           << "QualityScore\tCount\tQuantizedScore\n";
    for (std::size_t quality = 0; quality < bins.size(); ++quality)
        output << quality << '\t' << empirical_histogram[quality] << '\t'
               << quantized_map[quality] << '\n';
    output << "\n#:GATKTable:6:" << table0_rows
           << ":%s:%s:%.4f:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable0:\n"
           << "ReadGroup\tEventType\tEmpiricalQuality\tEstimatedQReported\tObservations\tErrors\n";
    for (const auto& [read_group, group] : groups) {
        std::uint64_t observations = 0, mismatches = 0;
        double expected_errors = 0.0;
        for (std::size_t quality = 0; quality < group.bins.size(); ++quality) {
            observations += group.bins[quality].count;
            mismatches += group.bins[quality].mismatches;
            expected_errors += static_cast<double>(group.bins[quality].count) *
                std::pow(10.0, -static_cast<double>(quality) / 10.0);
        }
        const double estimated = observations == 0 || expected_errors <= 0.0 ? 0.0 :
            -10.0 * std::log10(expected_errors / static_cast<double>(observations));
        // RecalTable0's estimated Q reported is a continuous Bayesian prior;
        // GATK passes that value directly into the empirical-quality solver
        // rather than rounding it to an integer first.  Keeping the
        // fractional prior is observable when multiple known-sites masks
        // change the aggregate error rate.
        const double empirical = static_cast<double>(bayesian_empirical_quality(
            QualityBin{observations, mismatches}, estimated));
        const auto rg = display_read_group(read_group);
        output << rg << "\tM\t" << std::fixed << std::setprecision(4) << empirical << '\t'
               << estimated << '\t' << observations << '\t' << std::setprecision(2)
               << static_cast<double>(mismatches) << '\n';
    }
    for (const auto& [key, bin] : indel_groups) {
        const auto rg = display_read_group(key.first);
        double expected_errors = 0.0;
        for (const auto& [quality_key, quality_bin] : covariates) {
            if (quality_key.read_group == key.first && quality_key.event == key.second)
                expected_errors += static_cast<double>(quality_bin.count) *
                    std::pow(10.0, -static_cast<double>(quality_key.quality) / 10.0);
        }
        const auto estimated = bin.count == 0 || expected_errors <= 0.0 ? 0.0 :
            -10.0 * std::log10(expected_errors / static_cast<double>(bin.count));
        const auto empirical = static_cast<double>(bayesian_empirical_quality(bin, estimated));
        output << rg << '\t' << key.second << '\t' << std::fixed << std::setprecision(4)
               << empirical << '\t' << std::fixed << std::setprecision(4)
               << estimated << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    output << "\n#:GATKTable:6:" << table1_rows
           << ":%s:%s:%s:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable1:\n"
           << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations\tErrors\n";
    for (const auto& [read_group, group] : groups) {
        const auto rg = display_read_group(read_group);
        for (std::size_t quality = 0; quality < group.bins.size(); ++quality) {
            const auto& bin = group.bins[quality];
            if (bin.count == 0) continue;
            output << rg << '\t' << quality << "\tM\t" << std::fixed << std::setprecision(4)
                   << empirical_quality(bin, static_cast<int>(quality)) << '\t' << bin.count
                   << '\t' << std::setprecision(2) << static_cast<double>(bin.mismatches) << '\n';
        }
    }
    for (const auto& [key, bin] : indel_quality_groups) {
        if (bin.count == 0) continue;
        const auto quality = std::get<1>(key);
        const auto event = std::get<2>(key);
        const auto empirical = empirical_quality(bin, quality);
        output << display_read_group(std::get<0>(key)) << '\t' << quality << '\t' << event
               << '\t' << std::fixed << std::setprecision(4) << empirical << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    output << "\n#:GATKTable:8:" << table2_rows
           << ":%s:%s:%s:%s:%s:%.4f:%d:%.2f:;\n"
           << "#:GATKTable:RecalTable2:\n"
           << "ReadGroup\tQualityScore\tCovariateValue\tCovariateName\tEventType\tEmpiricalQuality\tObservations\tErrors\n";
    for (const auto& [key, bin] : context_table) {
        const auto& [read_group, context, quality] = key;
        const auto rg = display_read_group(read_group);
        output << rg << '\t' << quality << '\t' << context << "\tContext\tM\t"
               << std::fixed << std::setprecision(4) << empirical_quality(bin, quality) << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    for (const auto& [key, bin] : cycle_table) {
        const auto& [read_group, cycle, quality] = key;
        const auto rg = display_read_group(read_group);
        output << rg << '\t' << quality << '\t' << cycle << "\tCycle\tM\t"
               << std::fixed << std::setprecision(4) << empirical_quality(bin, quality) << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    for (const auto& [key, bin] : indel_context_table) {
        const auto& [read_group, context, quality, event] = key;
        const auto rg = display_read_group(read_group);
        const auto empirical = empirical_quality(bin, quality);
        output << rg << '\t' << quality << '\t' << context << "\tContext\t" << event << '\t'
               << std::fixed << std::setprecision(4) << empirical << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    for (const auto& [key, bin] : indel_cycle_table) {
        const auto& [read_group, cycle, quality, event] = key;
        const auto rg = display_read_group(read_group);
        const auto empirical = empirical_quality(bin, quality);
        output << rg << '\t' << quality << '\t' << cycle << "\tCycle\t" << event << '\t'
               << std::fixed << std::setprecision(4) << empirical << '\t'
               << bin.count << '\t' << std::setprecision(2)
               << static_cast<double>(bin.mismatches) << '\n';
    }
    if (!output) throw std::runtime_error("cannot write recalibration report: " + path);
    output.flush();
    output.close();
    rewrite_gatk_report_fixed_width(path);
}

std::string manifest_json(const Options& options, const char* tool,
                          const fastgatk::runtime::ResourceSnapshot& resources,
                          std::uint64_t records, std::uint64_t bases,
                          double kernel_prepare_seconds = 0.0,
                          double kernel_execute_seconds = 0.0,
                          const std::string& kernel_execution_space = {},
                          const std::string& kernel_execution_policy = {},
                          bool checkpoint_resumed = false,
                          std::uint64_t checkpoint_records = 0,
                          bool checkpoint_prefix_reused = false,
                          double wall_seconds = 0.0,
                          const PipelineSummary& pipeline = {},
                          std::uint64_t covariate_kernel_observations = 0,
                          double covariate_kernel_prepare_seconds = 0.0,
                          double covariate_kernel_execute_seconds = 0.0,
                          const std::string& covariate_kernel_execution_space = {},
                          const std::string& covariate_kernel_execution_policy = {},
                          std::uint64_t covariate_kernel_workspace_bytes = 0,
                          bool covariate_kernel_team_local_histogram = false) {
    std::ostringstream out;
    const bool report = std::string(tool) == "BaseRecalibrator";
    const auto covariate_path = options.output + ".covariates.tsv";
    const bool primary_complete = file_complete(options.output);
    const bool covariate_complete = !report || file_complete(covariate_path);
    const bool index_complete = options.index_path.empty() || file_complete(options.index_path);
    const auto file_bytes = [](const std::string& path) -> std::uintmax_t {
        std::error_code error;
        const auto bytes = std::filesystem::file_size(path, error);
        return error ? 0 : bytes;
    };
    const auto execution_space = kernel_execution_space.empty()
        ? std::string("Host") : kernel_execution_space;
    out << "{\"schema_version\":1,\"tool\":\"" << tool
        << "\",\"implementation\":\"fastgatk-bqsr\",\"status\":\"prototype\""
        << ",\"execution_space\":\"" << json_escape(execution_space)
        << "\",\"determinism\":\"strict\""
        << ",\"primary_output\":\"" << json_escape(options.output) << "\""
        << ",\"primary_output_kind\":\"" << (report ? "recalibration-report" : "bam-or-cram") << "\""
        << ",\"compatibility\":{\"gatk_parameter_aliases\":true,"
        << "\"integer_table\":true,\"gatk_report_v1_1\":" << (report ? "true" : "false")
        << ",\"covariate_table\":" << (report ? "true" : "false")
        << ",\"hierarchical_covariate_model\":true"
        << ",\"covariate_integer_reduction\":" << (report ? "true" : "false")
        << ",\"report_recal_table_fallback\":" << (!report ? "true" : "false")
        << ",\"read_filter_controls\":" << (report ? "true" : "false")
        << ",\"disable_tool_default_read_filters\":"
        << (options.disable_tool_default_read_filters ? "true" : "false")
        << ",\"known_sites_filter\":"
        << (!options.known_sites.empty() ? "true" : "false")
        << ",\"known_sites_count\":" << options.known_sites.size()
        << ",\"preserve_qscores_less_than\":true"
        << ",\"use_original_qualities\":true"
        << ",\"emit_original_quals\":true"
        << ",\"quantize_quals\":true"
        << ",\"static_quantized_quals\":true"
        << ",\"round_down_quantized\":true"
        << ",\"allow_missing_read_group\":true"
        << ",\"global_qscore_prior\":true"
        << ",\"compute_indel_bqsr_tables\":true"
        << ",\"mismatches_context_size\":true"
        << ",\"maximum_cycle_value\":true"
        << ",\"resource_adaptive_batch\":true"
        << ",\"checkpoint_resume\":true"
        << ",\"create_output_bam_index\":"
        << (options.create_output_bam_index ? "true" : "false")
        << ",\"bit_identical_to_gatk\":false,\"vcf_index\":false,\"hts_index\":"
        << (!options.index_path.empty() ? "true" : "false")
        << "},\"outputs\":[{\"path\":\"" << json_escape(options.output)
        << "\",\"kind\":\"" << (report ? "recalibration-report" : "bam-or-cram")
        << "\",\"complete\":" << (primary_complete ? "true" : "false") << "}"
        << (report ? ",{\"path\":\"" + json_escape(covariate_path) +
             "\",\"kind\":\"bqsr-covariate-table\",\"complete\":" +
             (covariate_complete ? "true}" : "false}") : "")
        << (options.index_path.empty() ? "" : ",{\"path\":\"" + json_escape(options.index_path) +
             "\",\"kind\":\"bam-or-cram-index\",\"complete\":" +
             (index_complete ? "true}" : "false}"))
        << "]"
        << ",\"telemetry\":{\"resources\":" << resources.to_json()
        << ",\"records\":" << records << ",\"bases\":" << bases
        << ",\"output_bytes\":" << file_bytes(options.output)
        << ",\"covariate_bytes\":" << (report ? file_bytes(covariate_path) : 0)
        << ",\"index_bytes\":" << (options.index_path.empty() ? 0 : file_bytes(options.index_path))
        << ",\"wall_seconds\":" << wall_seconds
        << ",\"preserve_qualities_less_than\":" << options.preserve_qualities_less_than
        << ",\"use_original_qualities\":" << (options.use_original_qualities ? "true" : "false")
        << ",\"emit_original_qualities\":" << (options.emit_original_qualities ? "true" : "false")
        << ",\"quantization_levels\":" << options.quantization_levels
        << ",\"static_quantized_quals\":\""
        << json_escape(static_quantized_signature(options)) << "\""
        << ",\"round_down_quantized\":" << (options.round_down_quantized ? "true" : "false")
        << ",\"allow_missing_read_group\":" << (options.allow_missing_read_group ? "true" : "false")
        << ",\"create_output_bam_index\":"
        << (options.create_output_bam_index ? "true" : "false")
        << ",\"global_qscore_prior\":" << options.global_qscore_prior
        << ",\"compute_indel_bqsr_tables\":" << (options.compute_indel_bqsr_tables ? "true" : "false")
        << ",\"mismatches_context_size\":" << options.mismatches_context_size
        << ",\"maximum_cycle_value\":" << options.maximum_cycle_value
        << ",\"read_filter_controls\":" << (report ? "true" : "false")
        << ",\"disable_tool_default_read_filters\":"
        << (options.disable_tool_default_read_filters ? "true" : "false")
        << ",\"enabled_read_filters\":\""
        << json_escape(read_filter_signature(options.enabled_read_filters)) << "\""
        << ",\"disabled_read_filters\":\""
        << json_escape(read_filter_signature(options.disabled_read_filters)) << "\""
        << ",\"inverted_read_filters\":\""
        << json_escape(read_filter_signature(options.inverted_read_filters)) << "\""
        << ",\"indels_context_size\":" << options.indels_context_size
        << ",\"insertions_default_quality\":" << options.insertions_default_quality
        << ",\"deletions_default_quality\":" << options.deletions_default_quality
        << ",\"kernel_execution_space\":\"" << json_escape(kernel_execution_space)
        << "\",\"kernel_execution_policy\":\""
        << json_escape(kernel_execution_policy)
        << "\",\"kernel_prepare_seconds\":" << kernel_prepare_seconds
        << ",\"kernel_execute_seconds\":" << kernel_execute_seconds
        << ",\"covariate_kernel_observations\":" << covariate_kernel_observations
        << ",\"covariate_kernel_execution_space\":\""
        << json_escape(covariate_kernel_execution_space)
        << "\",\"covariate_kernel_execution_policy\":\""
        << json_escape(covariate_kernel_execution_policy)
        << "\",\"covariate_kernel_prepare_seconds\":"
        << covariate_kernel_prepare_seconds
        << ",\"covariate_kernel_execute_seconds\":"
        << covariate_kernel_execute_seconds
        << ",\"covariate_kernel_workspace_bytes\":"
        << covariate_kernel_workspace_bytes
        << ",\"covariate_kernel_team_local_histogram\":"
        << (covariate_kernel_team_local_histogram ? "true" : "false")
        << ",\"checkpoint_enabled\":" << (!options.checkpoint.empty() ? "true" : "false")
        << ",\"checkpoint_path\":\"" << json_escape(options.checkpoint)
        << "\",\"checkpoint_resumed\":" << (checkpoint_resumed ? "true" : "false")
        << ",\"checkpoint_records\":" << checkpoint_records
        << ",\"checkpoint_prefix_reused\":"
        << (checkpoint_prefix_reused ? "true" : "false")
        << ",\"requested_batch_records\":" << options.batch_records
        << ",\"initial_batch_records\":" << options.initial_batch_records
        << ",\"effective_batch_records\":" << options.effective_batch_records
        << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
        << ",\"interval_file_inputs\":" << options.interval_file_inputs
        << ",\"interval_file_records\":" << options.interval_file_records
        << ",\"pipeline_lifecycle\":\""
        << (pipeline.used ? "Host decode->bounded queue->Kokkos compute->encode->sink" : "not-used")
        << "\",\"pipeline_decoded_items\":" << pipeline.decoded_items
        << ",\"pipeline_computed_items\":" << pipeline.computed_items
        << ",\"pipeline_encoded_items\":" << pipeline.encoded_items
        << ",\"pipeline_decoded_bytes\":" << pipeline.decoded_bytes
        << ",\"pipeline_computed_bytes\":" << pipeline.computed_bytes
        << ",\"pipeline_encoded_bytes\":" << pipeline.encoded_bytes
        << ",\"pipeline_peak_decoded_bytes\":" << pipeline.peak_decoded_bytes
        << ",\"pipeline_peak_computed_bytes\":" << pipeline.peak_computed_bytes
        << ",\"pipeline_peak_encoded_bytes\":" << pipeline.peak_encoded_bytes
        << "}}\n";
    return out.str();
}

#if !FASTGATK_HAS_HTSLIB
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with FAST_GATK_HTSLIB_ROOT for BQSR");
}
#else

// BAM records are variable-sized HTSlib objects.  Own them explicitly across
// pipeline queues so a decode/compute/encode exception cannot leak a staged
// record or double-free it during checkpoint cleanup.
struct OwnedBam {
    bam1_t* value = nullptr;

    OwnedBam() = default;
    explicit OwnedBam(bam1_t* record) : value(record) {}
    OwnedBam(const OwnedBam&) = delete;
    OwnedBam& operator=(const OwnedBam&) = delete;
    OwnedBam(OwnedBam&& other) noexcept : value(other.value) { other.value = nullptr; }
    OwnedBam& operator=(OwnedBam&& other) noexcept {
        if (this == &other) return *this;
        if (value != nullptr) bam_destroy1(value);
        value = other.value;
        other.value = nullptr;
        return *this;
    }
    ~OwnedBam() {
        if (value != nullptr) bam_destroy1(value);
    }
};

struct ApplyDecodedBatch {
    std::vector<OwnedBam> records;
    std::vector<std::uint8_t> raw_qualities;
    std::vector<std::int16_t> deltas;
    std::vector<std::uint8_t> bypass_preserve_threshold;
};

struct ApplyComputedBatch {
    std::vector<OwnedBam> records;
    std::vector<std::uint8_t> raw_qualities;
    std::vector<std::uint8_t> adjusted;
    std::vector<std::uint8_t> bypass_preserve_threshold;
};

struct ApplyEncodedBatch {
    std::vector<OwnedBam> records;
};

// GATK's ApplyBQSR checks the reported quality before applying either the
// hierarchical model or a quantizer. Keep this predicate shared by decode
// (which suppresses the model delta) and encode (which suppresses the final
// quantizer) so low-quality sentinel bases cannot be altered by the Kokkos
// clamp or by a later quantization branch.
bool preserve_quality_score(int raw_quality, int threshold) noexcept {
    return raw_quality < threshold;
}

// BaseRecalibrator keeps BAM/CRAM decoding and covariate projection on Host,
// while the compute stage owns only one decoded ReadBatch and its local
// integer reduction.  Returning local tables instead of mutating the global
// map in the worker preserves the sink's deterministic merge order.
struct BaseDecodedBatch {
    fastgatk::io::ReadBatch batch;
};

struct BaseComputedBatch {
    std::array<QualityBin, 94> bins{};
    CovariateTable covariates;
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
};

struct BaseEncodedBatch {
    std::array<QualityBin, 94> bins{};
    CovariateTable covariates;
    std::uint64_t records = 0;
    std::uint64_t bases = 0;
};

std::uint64_t owned_bam_bytes(const OwnedBam& record) {
    if (record.value == nullptr) return sizeof(OwnedBam);
    return static_cast<std::uint64_t>(sizeof(OwnedBam) + sizeof(bam1_t) +
                                      std::max<int>(0, record.value->l_data));
}

std::uint64_t apply_decoded_batch_bytes(const ApplyDecodedBatch& batch) {
    std::uint64_t bytes = sizeof(ApplyDecodedBatch) + batch.raw_qualities.size() +
                          batch.deltas.size() + batch.bypass_preserve_threshold.size();
    for (const auto& record : batch.records) bytes += owned_bam_bytes(record);
    return std::max<std::uint64_t>(1, bytes);
}

std::uint64_t apply_computed_batch_bytes(const ApplyComputedBatch& batch) {
    std::uint64_t bytes = sizeof(ApplyComputedBatch) + batch.raw_qualities.size() +
                          batch.adjusted.size() + batch.bypass_preserve_threshold.size();
    for (const auto& record : batch.records) bytes += owned_bam_bytes(record);
    return std::max<std::uint64_t>(1, bytes);
}

std::uint64_t apply_encoded_batch_bytes(const ApplyEncodedBatch& batch) {
    std::uint64_t bytes = sizeof(ApplyEncodedBatch);
    for (const auto& record : batch.records) bytes += owned_bam_bytes(record);
    return std::max<std::uint64_t>(1, bytes);
}

std::uint64_t covariate_table_bytes(const CovariateTable& table) {
    std::uint64_t bytes = sizeof(CovariateTable);
    for (const auto& [key, value] : table) {
        bytes += sizeof(key) + sizeof(value) + key.read_group.size() + key.context.size();
    }
    return bytes;
}

std::uint64_t base_decoded_batch_bytes(const BaseDecodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(BaseDecodedBatch) + batch.batch.bytes());
}

std::uint64_t base_computed_batch_bytes(const BaseComputedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(BaseComputedBatch) + covariate_table_bytes(batch.covariates));
}

std::uint64_t base_encoded_batch_bytes(const BaseEncodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(BaseEncodedBatch) + covariate_table_bytes(batch.covariates));
}

std::map<std::string, std::string> read_group_platform_units(const std::string& path,
                                                              const std::string& reference) {
    std::map<std::string, std::string> identifiers;
    htsFile* input = hts_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open reads for RG/PU header: " + path);
    if (!reference.empty() && hts_set_fai_filename(input, reference.c_str()) != 0) {
        hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot configure reference for RG/PU header");
    }
    bam_hdr_t* header = sam_hdr_read(input);
    if (!header) {
        hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read BAM header for RG/PU covariate");
    }
    std::istringstream lines(header->text == nullptr ? std::string{} : header->text);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.rfind("@RG\t", 0) != 0) continue;
        std::string id, platform_unit;
        std::istringstream fields(line);
        std::string field;
        while (std::getline(fields, field, '\t')) {
            if (field.rfind("ID:", 0) == 0) id = field.substr(3);
            else if (field.rfind("PU:", 0) == 0) platform_unit = field.substr(3);
        }
        if (!id.empty()) identifiers[id] = platform_unit.empty() ? id : platform_unit;
    }
    bam_hdr_destroy(header);
    hts_close(input);
    return identifiers;
}

std::string read_group_identifier(const std::string& read_group,
                                  const std::map<std::string, std::string>& identifiers) {
    if (read_group.empty()) return {};
    if (const auto found = identifiers.find(read_group); found != identifiers.end())
        return found->second;
    return read_group;
}

void require_compatible_sequence_dictionary(const bam_hdr_t* expected,
                                             const bam_hdr_t* observed,
                                             const std::string& path) {
    if (expected == nullptr || observed == nullptr ||
        sam_hdr_nref(expected) != sam_hdr_nref(observed))
        throw std::runtime_error("BAD_INPUT: ApplyBQSR resume output header is incompatible: " + path);
    for (int tid = 0; tid < sam_hdr_nref(expected); ++tid) {
        const auto* expected_name = sam_hdr_tid2name(expected, tid);
        const auto* observed_name = sam_hdr_tid2name(observed, tid);
        if (expected_name == nullptr || observed_name == nullptr ||
            std::string(expected_name) != observed_name ||
            sam_hdr_tid2len(expected, tid) != sam_hdr_tid2len(observed, tid))
            throw std::runtime_error("BAD_INPUT: ApplyBQSR resume output sequence dictionary mismatch: " + path);
    }
}

std::string bam_read_group(const bam1_t* record) {
    const auto* tag = bam_aux_get(record, "RG");
    if (tag == nullptr || *tag != 'Z') return {};
    const auto* value = bam_aux2Z(tag);
    return value == nullptr ? std::string{} : std::string(value);
}

// ApplyBQSR's --use-original-qualities is implemented by BQSRReadTransformer
// as a reset from the SAM OQ tag before looking up covariates.  Keep the
// decoded qualities separate from the mutable BAM QUAL field so the output
// still carries the original OQ tag unchanged, matching HTSJDK behavior.
std::vector<std::uint8_t> bam_apply_qualities(const bam1_t* record,
                                              bool use_original_qualities) {
    const auto length = static_cast<std::size_t>(record->core.l_qseq);
    std::vector<std::uint8_t> qualities(length, 0);
    const auto* current = bam_get_qual(const_cast<bam1_t*>(record));
    for (std::size_t index = 0; index < length; ++index)
        qualities[index] = current[index] == 0xff
            ? 0 : static_cast<std::uint8_t>(std::min<int>(current[index], 93));
    if (!use_original_qualities) return qualities;

    const auto* tag = bam_aux_get(record, "OQ");
    if (tag == nullptr) return qualities;
    if (tag[0] != 'Z')
        throw std::runtime_error("BAD_INPUT: ApplyBQSR OQ tag must use SAM Z type");
    const auto* encoded = bam_aux2Z(tag);
    if (encoded == nullptr || std::strlen(encoded) != length)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR OQ tag length does not match read length");
    for (std::size_t index = 0; index < length; ++index) {
        const int phred = static_cast<unsigned char>(encoded[index]) - 33;
        if (phred < 0 || phred > 93)
            throw std::runtime_error("BAD_INPUT: ApplyBQSR OQ tag contains invalid FASTQ quality");
        qualities[index] = static_cast<std::uint8_t>(phred);
    }
    return qualities;
}

void bam_emit_original_qualities_if_absent(bam1_t* record) {
    if (bam_aux_get(record, "OQ") != nullptr) return;
    const auto length = static_cast<std::size_t>(record->core.l_qseq);
    const auto* current = bam_get_qual(record);
    std::string encoded(length, '!');
    for (std::size_t index = 0; index < length; ++index) {
        const int phred = current[index] == 0xff
            ? 0 : std::min<int>(current[index], 93);
        encoded[index] = static_cast<char>(phred + 33);
    }
    if (bam_aux_update_str(record, "OQ", static_cast<int>(encoded.size()), encoded.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot emit ApplyBQSR OQ tag");
}

std::vector<std::uint8_t> bam_read_bases(const bam1_t* record) {
    const char* alphabet = "=ACMGRSVTWYHKDBN";
    std::vector<std::uint8_t> bases(static_cast<std::size_t>(record->core.l_qseq));
    const auto* sequence = bam_get_seq(record);
    for (int index = 0; index < record->core.l_qseq; ++index)
        bases[static_cast<std::size_t>(index)] =
            static_cast<std::uint8_t>(alphabet[bam_seqi(sequence, index) & 15]);
    return bases;
}

int covariate_delta(const QualityBin& bin, int quality) {
    const auto empirical = bayesian_empirical_quality(bin, quality);
    return std::clamp(empirical - quality, -20, 20);
}

[[maybe_unused]] int run_base_recalibrator(Options& options,
                          const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto started = std::chrono::steady_clock::now();
    const fastgatk::runtime::AdaptiveController controller;
    const auto initial_limits = controller.initial(
        resources.budget(), fastgatk::runtime::WorkEstimate{
            static_cast<std::uint32_t>(std::min<std::size_t>(
                options.batch_records, std::numeric_limits<std::uint32_t>::max())),
            1, 2048, 1024, 2048});
    // Apply the allocation cap before the first HTSlib decode.  This is
    // checkpoint-safe because the effective initial batch is deterministic
    // for a given allocation and requested batch size.
    const auto initial_batch_records = std::min<std::size_t>(
        options.batch_records, initial_limits.max_reads);
    options.initial_batch_records = initial_batch_records;
    options.effective_batch_records = initial_batch_records;
    fastgatk::io::HtsReader reader(options.input, options.reference, options.regions,
                                   initial_batch_records);
    auto batch_limits = initial_limits;
    batch_limits.max_reads = static_cast<std::uint32_t>(std::min<std::size_t>(
        reader.batch_records(), std::numeric_limits<std::uint32_t>::max()));
    options.interval_file_inputs = reader.interval_file_inputs();
    options.interval_file_records = reader.interval_file_records();
    // GATK's ReadGroupCovariate is keyed by platform unit (PU) when one is
    // present, falling back to the read-group ID.  HTSlib batches expose the
    // RG tag, so normalize it at the Host boundary before building tables.
    const auto platform_units = read_group_platform_units(options.input, options.reference);
    const auto reference = load_reference(options.reference, reader.header());
    const auto known_sites = load_known_sites(options.known_sites, reader.header());
    const auto header_summary = reader.header();
    std::array<QualityBin, 94> bins{};
    CovariateTable covariates;
    std::uint64_t records = 0, bases = 0;
    std::uint64_t checkpoint_records = 0;
    bool checkpoint_resumed = false;
    if (options.resume_checkpoint) {
        const auto checkpoint = read_bqsr_checkpoint(options.checkpoint, options);
        records = checkpoint.records;
        bases = checkpoint.bases;
        bins = checkpoint.bins;
        covariates = checkpoint.covariates;
        checkpoint_records = checkpoint.records;
        checkpoint_resumed = true;
    }
    double kernel_prepare_seconds = 0.0;
    double kernel_execute_seconds = 0.0;
    std::string kernel_execution_space;
    std::string kernel_execution_policy;
    double covariate_kernel_prepare_seconds = 0.0;
    double covariate_kernel_execute_seconds = 0.0;
    std::string covariate_kernel_execution_space;
    std::string covariate_kernel_execution_policy;
    std::uint64_t covariate_kernel_observations = 0;
    std::uint64_t covariate_kernel_workspace_bytes = 0;
    bool covariate_kernel_team_local_histogram = false;
    std::uint64_t input_records_seen = 0;
    std::size_t completed_batches = 0;
    std::size_t observed_effective_batch_records = reader.batch_records();
    using Pipeline = fastgatk::runtime::ThreeStagePipeline<
        BaseDecodedBatch, BaseComputedBatch, BaseEncodedBatch>;
    auto decode = [&]() -> std::optional<BaseDecodedBatch> {
        fastgatk::io::ReadBatch batch;
        while (reader.next(batch)) {
            if (input_records_seen < checkpoint_records) {
                input_records_seen += batch.records();
                if (input_records_seen > checkpoint_records)
                    throw std::runtime_error("BAD_INPUT: BQSR checkpoint does not end at a batch boundary");
                continue;
            }
            input_records_seen += batch.records();
            // Adapt only at a decoded batch boundary.  The reader is owned by
            // this stage, so changing its next batch size is race-free.
            if (options.checkpoint.empty() && !options.resume_checkpoint) {
                fastgatk::runtime::RuntimeTelemetry runtime;
                runtime.host_bytes = batch.bytes();
                runtime.inflight_bytes = runtime.host_bytes;
                runtime.compute_queue_empty = true;
                const auto safe_budget = resources.safe_memory_budget_bytes();
                const auto pressure = safe_budget != 0 && runtime.host_bytes > safe_budget
                    ? fastgatk::runtime::Pressure::HostMemory
                    : fastgatk::runtime::Pressure::Normal;
                const auto next_limits = controller.next(runtime, batch_limits, pressure);
                if (next_limits.max_reads < reader.batch_records()) {
                    reader.set_batch_records(next_limits.max_reads);
                    ++options.adaptive_batch_reductions;
                }
                batch_limits = next_limits;
                observed_effective_batch_records = reader.batch_records();
            }
            return BaseDecodedBatch{std::move(batch)};
        }
        return std::nullopt;
    };
    auto compute = [&](BaseDecodedBatch decoded) -> std::optional<BaseComputedBatch> {
        auto batch = std::move(decoded.batch);
        BaseComputedBatch result;
        result.records = 0;
        result.bases = 0;
        fastgatk::kernels::BqsrObservationBatch observations;
        observations.qualities.reserve(batch.bases.size());
        observations.mismatches.reserve(batch.bases.size());
        fastgatk::kernels::BqsrCovariateObservationBatch covariate_observations;
        covariate_observations.covariate_ids.reserve(batch.bases.size());
        covariate_observations.mismatches.reserve(batch.bases.size());
        std::map<CovariateKey, std::uint32_t> covariate_ids;
        std::vector<CovariateKey> covariate_keys;
        covariate_keys.reserve(batch.bases.size());
        for (std::size_t record = 0; record < batch.records(); ++record) {
            if (!bqsr_default_read_filter(batch, record, header_summary, options)) continue;
            const auto begin = batch.offsets[record];
            const auto end = batch.offsets[record + 1];
            const auto [clip_left, clip_end] = bqsr_soft_clip_bounds(batch, record);
            if (clip_end <= clip_left) continue;
            ++result.records;
            result.bases += end - begin;
            const auto flags = record < batch.flags.size() ? batch.flags[record] : 0U;
            const std::span<const std::uint8_t> read_bases(
                batch.bases.data() + begin + clip_left, clip_end - clip_left);
            const std::span<const std::uint8_t> read_qualities(
                batch.qualities.data() + begin + clip_left, clip_end - clip_left);
            const auto indel_masks = options.compute_indel_bqsr_tables
                ? indel_masks_at(batch, record)
                : IndelMasks{};
            for (std::size_t offset = begin + clip_left; offset < begin + clip_end; ++offset) {
                const auto quality = std::min<std::size_t>(batch.qualities[offset], 93);
                if (quality < 6 || encode_base(batch.bases[offset]) > 3) continue;
                fastgatk::io::ReadProjection projection;
                const bool projected = bqsr_project_read_offset(
                    batch, record, offset - begin, projection);
                const bool has_reference = projected &&
                    projection.operation != fastgatk::io::CigarOpCode::Insertion &&
                    batch.tids[record] >= 0 &&
                    projection.reference_position >= 0 &&
                    static_cast<std::size_t>(batch.tids[record]) < reference.size() &&
                    static_cast<std::size_t>(projection.reference_position) <
                        reference[static_cast<std::size_t>(batch.tids[record])].size();
                const bool known = projected && known_site_contains(
                    known_sites, batch.tids[record], projection.reference_position);
                const auto read_group = read_group_identifier(
                    read_group_at(batch, record), platform_units);
                const auto clipped_offset = offset - begin - clip_left;
                const auto cycle = checked_cycle(flags, clipped_offset, clip_end - clip_left,
                                                 options.maximum_cycle_value);
                const auto full_read_offset = offset - begin;
                const auto add_event = [&](char event, int event_quality,
                                           std::size_t context_size, bool error) {
                    std::size_t indel_offset = 0;
                    const auto full_length = end - begin;
                    const auto effective_length = full_length >=
                        indel_masks.clipped_left + indel_masks.clipped_right
                        ? full_length - indel_masks.clipped_left - indel_masks.clipped_right : 0;
                    const bool inside_clipped_read = full_read_offset >= indel_masks.clipped_left &&
                        full_read_offset < indel_masks.clipped_left + effective_length;
                    if (inside_clipped_read) indel_offset = full_read_offset - indel_masks.clipped_left;
                    const bool valid_indel_cycle = inside_clipped_read && indel_offset >= 4 &&
                        indel_offset + 4 < effective_length;
                    const auto indel_cycle = valid_indel_cycle
                        ? checked_cycle(flags, indel_offset, effective_length,
                                        options.maximum_cycle_value)
                        : std::numeric_limits<std::int32_t>::min();
                    const auto covariate_cycle = event == 'M' ? cycle : indel_cycle;
                    CovariateKey key{read_group, covariate_cycle,
                        read_context_from(read_bases, read_qualities, flags,
                        clipped_offset, context_size),
                        event_quality, event};
                    auto found = covariate_ids.find(key);
                    std::uint32_t id = 0;
                    if (found == covariate_ids.end()) {
                        if (covariate_keys.size() >= std::numeric_limits<std::uint32_t>::max())
                            throw std::runtime_error("RESOURCE_EXHAUSTED: BQSR covariate id space exceeded");
                        id = static_cast<std::uint32_t>(covariate_keys.size());
                        covariate_ids.emplace(key, id);
                        covariate_keys.push_back(std::move(key));
                    } else {
                        id = found->second;
                    }
                    covariate_observations.covariate_ids.push_back(id);
                    covariate_observations.mismatches.push_back(error ? 1U : 0U);
                };
                // Java's isSNP array is initialized to zero for insertion
                // read bases. They are therefore valid substitution
                // observations (unless masked by a known site), even though
                // no reference base is consumed by the insertion.
                const bool insertion_observation = projected &&
                    projection.operation == fastgatk::io::CigarOpCode::Insertion;
                if (!known && (insertion_observation || has_reference) &&
                    (!has_reference || encode_base(static_cast<std::uint8_t>(
                        reference[static_cast<std::size_t>(batch.tids[record])]
                        [static_cast<std::size_t>(projection.reference_position)])) <= 3)) {
                    const bool mismatch = has_reference &&
                        encode_base(batch.bases[offset]) != encode_base(
                            static_cast<std::uint8_t>(reference[static_cast<std::size_t>(batch.tids[record])]
                                                       [static_cast<std::size_t>(projection.reference_position)]));
                    observations.qualities.push_back(static_cast<std::uint8_t>(quality));
                    observations.mismatches.push_back(mismatch ? 1 : 0);
                    add_event('M', static_cast<int>(quality),
                              static_cast<std::size_t>(options.mismatches_context_size), mismatch);
                }
                if (options.compute_indel_bqsr_tables && !known &&
                    full_read_offset < indel_masks.read_bases.size() && indel_masks.read_bases[full_read_offset] != 0) {
                    const auto insertion = full_read_offset < indel_masks.insertion.size()
                        ? indel_masks.insertion[full_read_offset] != 0 : false;
                    const auto deletion = full_read_offset < indel_masks.deletion.size()
                        ? indel_masks.deletion[full_read_offset] != 0 : false;
                    const auto insertion_quality = indel_quality_at(
                        batch, record, full_read_offset, 'I', options.insertions_default_quality);
                    const auto deletion_quality = indel_quality_at(
                        batch, record, full_read_offset, 'D', options.deletions_default_quality);
                    add_event('I', insertion_quality,
                              static_cast<std::size_t>(options.indels_context_size), insertion);
                    add_event('D', deletion_quality,
                              static_cast<std::size_t>(options.indels_context_size), deletion);
                }
            }
        }
        const auto batch_counts = fastgatk::kernels::count_bqsr_quality_kokkos(observations);
        kernel_prepare_seconds += batch_counts.prepare_seconds;
        kernel_execute_seconds += batch_counts.execute_seconds;
        kernel_execution_space = batch_counts.execution_space;
        kernel_execution_policy = batch_counts.execution_policy;
        for (std::size_t quality = 0; quality < result.bins.size(); ++quality) {
            result.bins[quality].count += batch_counts.bins[quality].count;
            result.bins[quality].mismatches += batch_counts.bins[quality].mismatches;
        }
        const auto covariate_counts = fastgatk::kernels::count_bqsr_covariates_kokkos(
            covariate_observations);
        covariate_kernel_prepare_seconds += covariate_counts.prepare_seconds;
        covariate_kernel_execute_seconds += covariate_counts.execute_seconds;
        covariate_kernel_execution_space = covariate_counts.execution_space;
        covariate_kernel_execution_policy = covariate_counts.execution_policy;
        covariate_kernel_observations += covariate_counts.observations;
        covariate_kernel_workspace_bytes = std::max(
            covariate_kernel_workspace_bytes, covariate_counts.workspace_bytes);
        covariate_kernel_team_local_histogram =
            covariate_kernel_team_local_histogram || covariate_counts.team_local_histogram;
        for (std::size_t id = 0; id < covariate_keys.size(); ++id) {
            const auto& key = covariate_keys[id];
            const auto& count = covariate_counts.covariates[id];
            auto& covariate = result.covariates[key];
            covariate.count += count.count;
            covariate.mismatches += count.mismatches;
        }
        return result;
    };
    auto encode = [](BaseComputedBatch batch) -> std::optional<BaseEncodedBatch> {
        return BaseEncodedBatch{std::move(batch.bins), std::move(batch.covariates),
                                batch.records, batch.bases};
    };
    auto sink = [&](BaseEncodedBatch batch) {
        records += batch.records;
        bases += batch.bases;
        for (std::size_t quality = 0; quality < bins.size(); ++quality)
            add_quality_bin(bins[quality], batch.bins[quality]);
        for (const auto& [key, value] : batch.covariates)
            add_quality_bin(covariates[key], value);
        ++completed_batches;
        if (!options.checkpoint.empty() &&
            completed_batches % options.checkpoint_every_batches == 0) {
            BqsrCheckpointState checkpoint{records, bases, bins, covariates,
                                           options.effective_batch_records};
            write_bqsr_checkpoint(options.checkpoint, options, checkpoint);
        }
    };
    const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    const auto nominal_batch_bytes = std::max<std::uint64_t>(
        4096ULL, static_cast<std::uint64_t>(initial_batch_records) * 256ULL);
    const auto stage_capacity = safe_memory_budget > 0
        ? std::max<std::uint64_t>(safe_memory_budget / 3ULL, nominal_batch_bytes * 2ULL)
        : std::max<std::uint64_t>(64ULL * 1024ULL * 1024ULL, nominal_batch_bytes * 2ULL);
    Pipeline pipeline(
        Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
        std::move(decode), std::move(compute), std::move(encode), std::move(sink),
        base_decoded_batch_bytes, base_computed_batch_bytes, base_encoded_batch_bytes);
    const auto pipeline_metrics = pipeline.run();
    // Publish decode-stage scheduling telemetry only after all worker/sink
    // callbacks have joined.  Checkpoint writes in sink therefore never race
    // with the adaptive decoder's mutable bookkeeping.
    options.effective_batch_records = observed_effective_batch_records;
    const PipelineSummary pipeline_summary{
        true,
        pipeline_metrics.decoded_items,
        pipeline_metrics.computed_items,
        pipeline_metrics.encoded_items,
        pipeline_metrics.decoded_bytes,
        pipeline_metrics.computed_bytes,
        pipeline_metrics.encoded_bytes,
        pipeline_metrics.peak_decoded_bytes,
        pipeline_metrics.peak_computed_bytes,
        pipeline_metrics.peak_encoded_bytes};
    if (input_records_seen < checkpoint_records)
        throw std::runtime_error("BAD_INPUT: BQSR checkpoint is beyond end of input");
    if (!options.checkpoint.empty()) {
        BqsrCheckpointState checkpoint{records, bases, bins, covariates,
                                       options.effective_batch_records};
        write_bqsr_checkpoint(options.checkpoint, options, checkpoint);
    }
    write_gatk_report(options.output, bins, covariates, options);
    const auto covariate_path = options.output + ".covariates.tsv";
    std::ofstream covariate_output(covariate_path);
    if (!covariate_output) throw std::runtime_error("cannot write BQSR covariate table: " + covariate_path);
    covariate_output << "# FASTGATK-BQSR-COVARIATES v2\n"
                     << "# read_group\tcycle\tcontext\traw_quality\tcount\tmismatches\tempirical_quality\tdelta\tevent_type\n"
                     << std::setprecision(12);
    for (const auto& [key, bin] : covariates) {
        const double empirical = empirical_quality(bin, key.quality);
        const auto correction = static_cast<int>(std::lround(std::clamp(
            empirical - static_cast<double>(key.quality), -20.0, 20.0)));
        covariate_output << (key.read_group.empty() ? "-" : key.read_group) << '\t'
                         << key.cycle << '\t' << (key.context.empty() ? "-" : key.context)
                         << '\t' << key.quality << '\t'
                         << bin.count << '\t' << bin.mismatches << '\t' << empirical << '\t'
                         << correction << '\t' << key.event << '\n';
    }
    if (!covariate_output) throw std::runtime_error("cannot write BQSR covariate table: " + covariate_path);
    covariate_output.flush();
    if (!file_complete(options.output) || !file_complete(covariate_path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: BQSR report or covariate table is missing/empty");
    const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("cannot write BQSR manifest: " + manifest_path);
    manifest << manifest_json(options, "BaseRecalibrator", resources, records, bases,
                              kernel_prepare_seconds, kernel_execute_seconds,
                              kernel_execution_space, kernel_execution_policy,
                              checkpoint_resumed, checkpoint_records, false,
                              std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count(),
                              pipeline_summary, covariate_kernel_observations,
                              covariate_kernel_prepare_seconds,
                              covariate_kernel_execute_seconds,
                              covariate_kernel_execution_space,
                              covariate_kernel_execution_policy,
                              covariate_kernel_workspace_bytes,
                              covariate_kernel_team_local_histogram);
    std::cout << "{\"tool\":\"BaseRecalibrator\",\"status\":\"prototype\",\"records\":"
              << records << ",\"bases\":" << bases << ",\"covariate_rows\":"
              << covariates.size() << ",\"covariate_observations\":"
              << covariate_kernel_observations << "}\n";
    return 0;
}

[[maybe_unused]] int run_apply_bqsr(Options& options,
                   const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto started = std::chrono::steady_clock::now();
    if (options.input == options.output) throw std::invalid_argument("input and output must differ");
    const fastgatk::runtime::AdaptiveController controller;
    const auto initial_limits = controller.initial(
        resources.budget(), fastgatk::runtime::WorkEstimate{
            static_cast<std::uint32_t>(std::min<std::size_t>(
                options.batch_records, std::numeric_limits<std::uint32_t>::max())),
            1, 2048, 1024, 2048});
    // A checkpointed ApplyBQSR stream keeps the requested batch schedule so
    // the durable prefix can be resumed at the same record boundary.  Fresh
    // runs are capped before staging the first BAM records.
    const bool adaptive_decode = options.checkpoint.empty() && !options.resume_checkpoint;
    const auto decode_batch_records = adaptive_decode
        ? std::min<std::size_t>(options.batch_records, initial_limits.max_reads)
        : options.batch_records;
    options.initial_batch_records = decode_batch_records;
    options.effective_batch_records = decode_batch_records;
    std::uint64_t checkpoint_records = 0;
    bool checkpoint_resumed = false;
    if (options.resume_checkpoint) {
        const auto checkpoint = read_apply_bqsr_checkpoint(options.checkpoint, options);
        checkpoint_records = checkpoint.records;
        checkpoint_resumed = true;
        if (!std::filesystem::is_regular_file(options.output))
            throw std::runtime_error("BAD_INPUT: ApplyBQSR resume requires the checkpoint output prefix: " + options.output);
    }
    const auto report_model = read_report(options.recalibration);
    // Manifest telemetry must describe the effective report-owned model
    // dimension rather than ApplyBQSR's parser default.
    options.mismatches_context_size = report_model.mismatches_context_size;
    const auto static_quantized_mapping = construct_static_quantized_mapping(
        options.static_quantized_quals, options.round_down_quantized);
    const auto dynamic_quantized_mapping = construct_dynamic_quantized_mapping(
        report_model, options.quantization_levels);
    const auto covariates = read_covariates(options.recalibration + ".covariates.tsv");
    const auto covariate_model = build_covariate_model(covariates);
    const char* input_mode = suffix(options.input, ".cram") ? "rc" : "r";
    const char* output_mode = suffix(options.output, ".cram") ? "wc" :
                              suffix(options.output, ".sam") ? "w" : "wb";
    htsFile* input = hts_open(options.input.c_str(), input_mode);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open input: " + options.input);
    if (!options.reference.empty() && hts_set_fai_filename(input, options.reference.c_str()) != 0) {
        hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot configure input reference");
    }
    bam_hdr_t* header = sam_hdr_read(input);
    if (!header) { hts_close(input); throw std::runtime_error("BAD_INPUT: cannot read input header"); }
    const auto platform_units = read_group_platform_units(options.input, options.reference);
    const bool resume_prefix = checkpoint_resumed;
    const std::string temporary_output = options.output + ".resume.tmp";
    htsFile* prefix_input = nullptr;
    bam_hdr_t* prefix_header = nullptr;
    if (resume_prefix) {
        prefix_input = hts_open(options.output.c_str(), "r");
        if (!prefix_input) {
            bam_hdr_destroy(header); hts_close(input);
            throw std::runtime_error("BAD_INPUT: cannot open ApplyBQSR checkpoint output: " + options.output);
        }
        if (!options.reference.empty() && hts_set_fai_filename(prefix_input, options.reference.c_str()) != 0) {
            hts_close(prefix_input); bam_hdr_destroy(header); hts_close(input);
            throw std::runtime_error("BAD_INPUT: cannot configure checkpoint output reference");
        }
        prefix_header = sam_hdr_read(prefix_input);
        if (!prefix_header) {
            hts_close(prefix_input); bam_hdr_destroy(header); hts_close(input);
            throw std::runtime_error("BAD_INPUT: cannot read ApplyBQSR checkpoint output header");
        }
        require_compatible_sequence_dictionary(header, prefix_header, options.output);
    }
    htsFile* output = hts_open(resume_prefix ? temporary_output.c_str() : options.output.c_str(), output_mode);
    if (!output) { bam_hdr_destroy(header); hts_close(input); throw std::runtime_error("cannot open output: " + options.output); }
    // HTSlib currently defaults to CRAM 3.1, while the pinned GATK/HTSJDK
    // 4.6.2.0 reader accepts CRAM 3.0 but rejects 3.1.  Select the portable
    // 3.0 container at the writer boundary so ApplyBQSR CRAM output remains
    // directly reopenable by the Java tool.  BAM/SAM output is unaffected.
    if (suffix(options.output, ".cram") &&
        hts_set_opt(output, CRAM_OPT_VERSION, "3.0") != 0) {
        hts_close(output); if (prefix_header) bam_hdr_destroy(prefix_header);
        if (prefix_input) hts_close(prefix_input);
        bam_hdr_destroy(header); hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot configure CRAM 3.0 output");
    }
    if (suffix(options.output, ".cram") &&
        hts_set_opt(output, CRAM_OPT_IGNORE_MD5, 1) != 0) {
        hts_close(output); if (prefix_header) bam_hdr_destroy(prefix_header);
        if (prefix_input) hts_close(prefix_input);
        bam_hdr_destroy(header); hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot configure CRAM reference-MD5 compatibility");
    }
    if (!options.reference.empty() && hts_set_fai_filename(output, options.reference.c_str()) != 0) {
        hts_close(output); if (prefix_header) bam_hdr_destroy(prefix_header);
        if (prefix_input) hts_close(prefix_input);
        bam_hdr_destroy(header); hts_close(input);
        throw std::runtime_error("BAD_INPUT: cannot configure output reference");
    }
    if (sam_hdr_write(output, header) < 0) {
        hts_close(output); if (prefix_header) bam_hdr_destroy(prefix_header);
        if (prefix_input) hts_close(prefix_input);
        bam_hdr_destroy(header); hts_close(input);
        throw std::runtime_error("cannot write output header");
    }
    bool checkpoint_prefix_reused = false;
    if (resume_prefix) {
        bam1_t* prefix_record = bam_init1();
        if (!prefix_record) {
            hts_close(output); bam_hdr_destroy(prefix_header); hts_close(prefix_input);
            bam_hdr_destroy(header); hts_close(input);
            throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed for ApplyBQSR checkpoint prefix");
        }
        for (std::uint64_t copied = 0; copied < checkpoint_records; ++copied) {
            if (sam_read1(prefix_input, prefix_header, prefix_record) < 0 ||
                sam_write1(output, header, prefix_record) < 0) {
                bam_destroy1(prefix_record); hts_close(output);
                bam_hdr_destroy(prefix_header); hts_close(prefix_input);
                bam_hdr_destroy(header); hts_close(input);
                throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint output is shorter than checkpoint");
            }
        }
        bam_destroy1(prefix_record);
        bam_hdr_destroy(prefix_header);
        hts_close(prefix_input);
        prefix_header = nullptr;
        prefix_input = nullptr;
        checkpoint_prefix_reused = true;
    }
    bam1_t* record = bam_init1();
    if (!record) { hts_close(output); bam_hdr_destroy(header); hts_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed"); }
    std::uint64_t records = checkpoint_prefix_reused ? checkpoint_records : 0;
    std::uint64_t bases = 0;
    if (checkpoint_prefix_reused) {
        const auto checkpoint = read_apply_bqsr_checkpoint(options.checkpoint, options);
        bases = checkpoint.bases;
    }
    std::uint64_t input_records_seen = 0;
    std::uint64_t completed_batches = checkpoint_prefix_reused
        ? checkpoint_records / options.batch_records : 0;
    double kernel_prepare_seconds = 0.0;
    double kernel_execute_seconds = 0.0;
    std::string kernel_execution_space;
    std::string kernel_execution_policy;
    // ApplyBQSR uses the same bounded decode -> Kokkos transform -> encode
    // lifecycle as the read-count tools.  BAM ownership remains Host-side;
    // only flat quality/delta arrays cross into the Kokkos kernel.
    using Pipeline = fastgatk::runtime::ThreeStagePipeline<
        ApplyDecodedBatch, ApplyComputedBatch, ApplyEncodedBatch>;
    Pipeline::Metrics pipeline_metrics;
    const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    const auto nominal_batch_bytes = std::max<std::uint64_t>(
        4096ULL, static_cast<std::uint64_t>(decode_batch_records) * 1024ULL);
    const auto stage_capacity = safe_memory_budget > 0
        ? std::max<std::uint64_t>(safe_memory_budget / 3ULL, nominal_batch_bytes * 2ULL)
        : std::max<std::uint64_t>(64ULL * 1024ULL * 1024ULL, nominal_batch_bytes * 2ULL);
    auto decode = [&]() -> std::optional<ApplyDecodedBatch> {
        ApplyDecodedBatch batch;
        batch.records.reserve(decode_batch_records);
        batch.raw_qualities.reserve(decode_batch_records * 150ULL);
        batch.deltas.reserve(decode_batch_records * 150ULL);
        while (batch.records.size() < decode_batch_records && sam_read1(input, header, record) >= 0) {
            if (input_records_seen < checkpoint_records) {
                ++input_records_seen;
                continue;
            }
            ++input_records_seen;
            OwnedBam staged_owner{bam_dup1(record)};
            if (staged_owner.value == nullptr)
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot stage ApplyBQSR batch record");
            if (options.emit_original_qualities)
                bam_emit_original_qualities_if_absent(staged_owner.value);
            const auto read_group = read_group_identifier(bam_read_group(record), platform_units);
            const auto read_bases = bam_read_bases(record);
            const auto read_qualities = bam_apply_qualities(record, options.use_original_qualities);
            const bool read_group_known =
                report_model.read_group.find(read_group) != report_model.read_group.end() ||
                report_model.read_group_reported_quality.find(read_group) !=
                    report_model.read_group_reported_quality.end() ||
                covariate_model.read_group.find(read_group) != covariate_model.read_group.end() ||
                covariate_model.read_group_reported_quality.find(read_group) !=
                    covariate_model.read_group_reported_quality.end();
            if (!read_group_known && !options.allow_missing_read_group)
                throw std::runtime_error(
                    "BAD_INPUT: read group not found in recalibration table: " +
                    (read_group.empty() ? std::string("UNKNOWN") : read_group) +
                    "; use --allow-missing-read-group to continue");
            for (int index = 0; index < record->core.l_qseq; ++index) {
                const auto raw = read_qualities[static_cast<std::size_t>(index)];
                if (!read_group_known && options.allow_missing_read_group) {
                    // GATK's missing-RG escape hatch bypasses model lookup and
                    // maps the original (or OQ-reset) quality directly through
                    // the selected quantizer.  A zero delta keeps the shared
                    // Kokkos transform in the same bounded pipeline.
                    batch.raw_qualities.push_back(raw);
                    batch.deltas.push_back(0);
                    batch.bypass_preserve_threshold.push_back(1);
                    continue;
                }
                const auto context = read_context_from(
                    read_bases, read_qualities, record->core.flag,
                    static_cast<std::size_t>(index),
                    static_cast<std::size_t>(report_model.mismatches_context_size));
                const CovariateKey key{
                    read_group,
                    read_cycle(record->core.flag, static_cast<std::size_t>(index),
                               static_cast<std::size_t>(record->core.l_qseq)),
                    context, raw};
                // GATK preserves low-quality bases before consulting the
                // recalibration tables (default threshold is Q6).
                const int delta = preserve_quality_score(raw, options.preserve_qualities_less_than)
                    ? 0 : model_delta(covariates, covariate_model, report_model, key,
                                       options.global_qscore_prior);
                batch.raw_qualities.push_back(static_cast<std::uint8_t>(raw));
                batch.deltas.push_back(static_cast<std::int16_t>(delta));
                batch.bypass_preserve_threshold.push_back(0);
            }
            batch.records.emplace_back(std::move(staged_owner));
        }
        if (batch.records.empty()) return std::nullopt;
        return batch;
    };
    auto compute = [&](ApplyDecodedBatch batch) -> std::optional<ApplyComputedBatch> {
        ApplyComputedBatch result;
        result.records = std::move(batch.records);
        result.raw_qualities = std::move(batch.raw_qualities);
        result.bypass_preserve_threshold = std::move(batch.bypass_preserve_threshold);
        if (result.raw_qualities.size() != batch.deltas.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: ApplyBQSR quality/delta width mismatch");
        if (result.raw_qualities.size() != result.bypass_preserve_threshold.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: ApplyBQSR preservation mask width mismatch");
        if (!batch.deltas.empty()) {
            const auto transformed = fastgatk::kernels::apply_bqsr_quality_kokkos(
                result.raw_qualities, batch.deltas);
            kernel_prepare_seconds += transformed.prepare_seconds;
            kernel_execute_seconds += transformed.execute_seconds;
            kernel_execution_space = transformed.execution_space;
            kernel_execution_policy = transformed.execution_policy;
            result.adjusted = transformed.adjusted;
        }
        return result;
    };
    auto encode = [&](ApplyComputedBatch batch) -> std::optional<ApplyEncodedBatch> {
        ApplyEncodedBatch result;
        result.records = std::move(batch.records);
        std::size_t cursor = 0;
        for (auto& owned : result.records) {
            if (owned.value == nullptr)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: null ApplyBQSR staged record");
            auto* staged_qualities = bam_get_qual(owned.value);
            for (int index = 0; index < owned.value->core.l_qseq; ++index) {
                if (cursor >= batch.raw_qualities.size() || cursor >= batch.adjusted.size())
                    throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: ApplyBQSR transformed width mismatch");
                const auto raw = batch.raw_qualities[cursor];
                const auto adjusted = batch.adjusted[cursor];
                const bool bypass_preserve = batch.bypass_preserve_threshold[cursor] != 0;
                // GATK applies the dynamic/static quantization map only to
                // bases that pass preserve-qscores-less-than.  Missing/Q0
                // qualities therefore remain untouched.
                const auto bounded_quality = std::min<std::size_t>(
                    adjusted, report_model.quantized_qualities.size() - 1);
                const auto dynamic_quantized = options.quantization_levels == 0
                    ? (bypass_preserve ? static_cast<int>(raw) : static_cast<int>(adjusted))
                    : options.quantization_levels < 0
                        ? report_model.quantized_qualities[bounded_quality]
                        : dynamic_quantized_mapping[bounded_quality];
                const auto quantized = options.static_quantized_quals.empty()
                    ? dynamic_quantized
                    : static_quantized_mapping[
                        std::min<std::size_t>(static_cast<std::size_t>(dynamic_quantized),
                                              static_quantized_mapping.size() - 1)];
                staged_qualities[index] = !bypass_preserve &&
                    preserve_quality_score(raw, options.preserve_qualities_less_than)
                    ? raw : quantized;
                ++cursor;
            }
        }
        if (cursor != batch.raw_qualities.size() || cursor != batch.adjusted.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: ApplyBQSR encoded width mismatch");
        return result;
    };
    auto sink = [&](ApplyEncodedBatch batch) {
        for (auto& owned : batch.records) {
            if (owned.value == nullptr || sam_write1(output, header, owned.value) < 0)
                throw std::runtime_error("cannot write ApplyBQSR output record");
            ++records;
            bases += static_cast<std::uint64_t>(owned.value->core.l_qseq);
        }
        ++completed_batches;
        if (!options.checkpoint.empty() &&
            completed_batches % options.checkpoint_every_batches == 0)
            write_apply_bqsr_checkpoint(options.checkpoint, options,
                                        ApplyBqsrCheckpointState{records, bases});
    };
    Pipeline pipeline(
        Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
        std::move(decode), std::move(compute), std::move(encode), std::move(sink),
        apply_decoded_batch_bytes, apply_computed_batch_bytes, apply_encoded_batch_bytes);
    pipeline_metrics = pipeline.run();
    const PipelineSummary pipeline_summary{
        true,
        pipeline_metrics.decoded_items,
        pipeline_metrics.computed_items,
        pipeline_metrics.encoded_items,
        pipeline_metrics.decoded_bytes,
        pipeline_metrics.computed_bytes,
        pipeline_metrics.encoded_bytes,
        pipeline_metrics.peak_decoded_bytes,
        pipeline_metrics.peak_computed_bytes,
        pipeline_metrics.peak_encoded_bytes};
    bam_destroy1(record);
    if (hts_close(output) != 0) {
        if (prefix_header) bam_hdr_destroy(prefix_header);
        if (prefix_input) hts_close(prefix_input);
        bam_hdr_destroy(header);
        hts_close(input);
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize BAM/CRAM output");
    }
    bam_hdr_destroy(header);
    hts_close(input);
    if (input_records_seen < checkpoint_records)
        throw std::runtime_error("BAD_INPUT: ApplyBQSR checkpoint is beyond end of input");
    if (resume_prefix) {
        std::error_code error;
        std::filesystem::rename(temporary_output, options.output, error);
        if (error)
            throw std::runtime_error("cannot atomically publish resumed ApplyBQSR output: " +
                                     error.message());
    }
    if (options.create_output_bam_index &&
        (suffix(options.output, ".bam") || suffix(options.output, ".cram"))) {
        options.index_path = options.output + (suffix(options.output, ".cram") ? ".crai" : ".bai");
        if (sam_index_build3(options.output.c_str(), options.index_path.c_str(), 0, 0) != 0)
            throw std::runtime_error("cannot build BAM/CRAM index: " + options.output);
    }
    if (!file_complete(options.output) ||
        (!options.index_path.empty() && !file_complete(options.index_path)))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: ApplyBQSR output or index is missing/empty");
    const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("cannot write BQSR manifest: " + manifest_path);
    manifest << manifest_json(options, "ApplyBQSR", resources, records, bases,
                              kernel_prepare_seconds, kernel_execute_seconds,
                              kernel_execution_space, kernel_execution_policy,
                              checkpoint_resumed,
                              checkpoint_records, checkpoint_prefix_reused,
                              std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count(),
                              pipeline_summary);
    std::cout << "{\"tool\":\"ApplyBQSR\",\"status\":\"prototype\",\"records\":"
              << records << ",\"bases\":" << bases << ",\"covariate_rows\":"
              << covariates.size() << "}\n";
    return 0;
}

int run_tool(Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
#if defined(FASTGATK_BQSR_APPLY)
    return run_apply_bqsr(options, resources);
#else
    return run_base_recalibrator(options, resources);
#endif
}
#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        Kokkos::initialize();
        initialized = true;
        auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        const auto status = run_tool(options, resources);
        Kokkos::finalize();
        initialized = false;
        return status;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
