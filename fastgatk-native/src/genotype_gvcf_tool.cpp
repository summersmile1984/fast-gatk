#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <set>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string reference;
    std::vector<std::string> regions;
    // GATK -XL/--exclude-intervals is applied after include-set
    // construction. Keep exclusions separate so an exclusion-only request
    // still traverses the complete input while removing overlapping records.
    std::vector<std::string> excluded_regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::string output;
    std::string manifest;
    bool create_index = true;
    // Match GATK GenotypeCalculationArgumentCollection defaults.  Priors are
    // only used by the optional posterior-QUAL path; PL-based GT/GQ
    // assignment remains unchanged by default.
    double snp_heterozygosity = 1.0e-3;
    double indel_heterozygosity = 1.0 / 8000.0;
    double heterozygosity_stdev = 0.01;
    // GenotypeGVCFs emits an ALT only when its independent AF posterior
    // passes this phred threshold.  Keeping the threshold explicit makes the
    // native output-allele subset match GATK's standard-min-confidence-
    // threshold-for-calling instead of merely retaining every ALT present in
    // a source shard.
    double standard_confidence_for_calling = 30.0;
    // GATK's GenotypingEngine reduces the proper ALT set before AF
    // calculation.  The symbolic <NON_REF> sentinel is never counted toward
    // this bound; the default is the GATK 4.6.2.0 value of six.
    int max_alternate_alleles = 6;
    // Optional GATK GenotypingEngine annotation.  NDA records the number of
    // concrete ALT alleles discovered before max-ALT subsetting, so callers
    // can distinguish an intentionally reduced site from a site that never
    // contained the dropped alleles.
    bool annotate_with_num_discovered_alleles = false;
    bool use_genotype_priors = true;
    bool use_posteriors_to_calculate_qual = false;
    // GATK's optional --include-non-variant-sites emits REF-only records from
    // reference-confidence blocks instead of dropping groups with no concrete
    // ALT.  The default remains variant-only materialization.
    bool include_non_variant_sites = false;
    // Deprecated GATK spelling for writer-side STARTS_IN interval filtering.
    // Traversal must still retain records that merely overlap -L so a long
    // deletion can contribute joint likelihoods at an interior locus; only
    // the finalized VariantContext is filtered by its start coordinate.
    bool only_output_calls_starting_in_intervals = false;
    std::uint64_t output_interval_skipped = 0;
    // GATK's GenotypingEngine uses the new AF calculator by default.  The
    // explicit negation is retained for compatibility/debugging and leaves
    // the existing PL/GT/AC path untouched.
    bool use_new_qual_calculator = true;
    // GenotypeGVCFs itself forces PREFER_PLS in the minimal genotyping
    // engine.  Keep that as the native default and expose the same enum for
    // callers that need posterior/no-call assignment semantics.
    std::string genotype_assignment_method = "PREFER_PLS";
    // Keep native diagnostic annotations available by default, but expose an
    // explicit GATK-compatible profile that removes raw/intermediate INFO and
    // caller-specific FORMAT fields after the standard annotations have been
    // finalized.  This profile is semantic/field compatible; HTSlib text
    // formatting is still reported separately by the oracle.
    bool gatk_annotation_compatibility = false;
    // Explicit bounded joint-locus mode.  The default aggregate path remains
    // available for byte-for-byte regression work; this switch makes input
    // decode and cross-shard merge lazy, retaining at most one locus per
    // source plus the three-stage pipeline queues.
    bool stream_by_locus = false;
    std::uint64_t streamed_loci = 0;
    std::uint64_t streamed_peak_host_bytes = 0;
    std::uint64_t stream_max_inflight_records = 0;
    bool pipeline_used = false;
    std::uint64_t pipeline_decoded_items = 0;
    std::uint64_t pipeline_computed_items = 0;
    std::uint64_t pipeline_encoded_items = 0;
    std::uint64_t pipeline_decoded_bytes = 0;
    std::uint64_t pipeline_computed_bytes = 0;
    std::uint64_t pipeline_encoded_bytes = 0;
    std::uint64_t pipeline_peak_decoded_bytes = 0;
    std::uint64_t pipeline_peak_computed_bytes = 0;
    std::uint64_t pipeline_peak_encoded_bytes = 0;
    std::uint64_t pipeline_stage_capacity_bytes = 0;
    // GATK does not materialize an orphan spanning-deletion '*' ALT when no
    // concrete deletion record covers the locus.  Count the bounded cleanup
    // at the joint-locus union boundary for compatibility diagnostics.
    std::uint64_t orphan_spanning_deletion_loci = 0;
};

struct GenotypeKernelTelemetry {
    std::uint64_t calls = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::string execution_space;
    std::uint64_t pl_remap_calls = 0;
    double pl_remap_prepare_seconds = 0.0;
    double pl_remap_execute_seconds = 0.0;
    std::string pl_remap_execution_space;
    std::uint64_t allele_field_remap_calls = 0;
    double allele_field_remap_prepare_seconds = 0.0;
    double allele_field_remap_execute_seconds = 0.0;
    std::string allele_field_remap_execution_space;
    std::uint64_t allele_count_calls = 0;
    double allele_count_prepare_seconds = 0.0;
    double allele_count_execute_seconds = 0.0;
    std::string allele_count_execution_space;
    std::uint64_t posterior_calls = 0;
    double posterior_prepare_seconds = 0.0;
    double posterior_execute_seconds = 0.0;
    std::string posterior_execution_space;
    std::uint64_t posterior_samples = 0;
    std::uint64_t cross_sample_reference_calls = 0;
    double cross_sample_reference_prepare_seconds = 0.0;
    double cross_sample_reference_execute_seconds = 0.0;
    std::string cross_sample_reference_execution_space;
    std::uint64_t cross_sample_reference_samples = 0;
    std::uint64_t cohort_calls = 0;
    double cohort_prepare_seconds = 0.0;
    double cohort_execute_seconds = 0.0;
    std::string cohort_execution_space;
    std::uint64_t cohort_samples = 0;
    std::uint64_t cohort_approximate_gq_samples = 0;
    std::uint64_t cohort_iterations = 0;
    std::uint64_t cohort_converged = 0;
    std::uint64_t excess_het_calls = 0;
    std::uint64_t inbreeding_coeff_calls = 0;
    std::uint64_t posterior_assignment_calls = 0;
    double posterior_assignment_prepare_seconds = 0.0;
    double posterior_assignment_execute_seconds = 0.0;
    std::string posterior_assignment_execution_space;
    std::uint64_t posterior_assignment_samples = 0;
    std::uint64_t posterior_annotation_assignment_calls = 0;
    double posterior_annotation_assignment_prepare_seconds = 0.0;
    double posterior_annotation_assignment_execute_seconds = 0.0;
    std::string posterior_annotation_assignment_execution_space;
    std::uint64_t posterior_annotation_assignment_samples = 0;
    std::uint64_t genotype_prior_calls = 0;
    double genotype_prior_prepare_seconds = 0.0;
    double genotype_prior_execute_seconds = 0.0;
    std::string genotype_prior_execution_space;
    // Number of FORMAT annotations removed by
    // SET_TO_NO_CALL_NO_ANNOTATIONS. GT is retained as a no-call; every other
    // genotype-level field is cleared, including caller-specific FORMAT IDs.
    std::uint64_t no_annotation_format_fields = 0;
    std::uint64_t output_allele_pruning_calls = 0;
    std::uint64_t output_alleles_pruned = 0;
    std::uint64_t max_alt_pruning_calls = 0;
    std::uint64_t max_alt_alleles_pruned = 0;
    std::uint64_t max_alt_score_kernel_calls = 0;
    double max_alt_score_prepare_seconds = 0.0;
    double max_alt_score_execute_seconds = 0.0;
    std::string max_alt_score_execution_space;
    std::uint64_t num_discovered_alleles_calls = 0;
};

GenotypeKernelTelemetry genotype_kernel_telemetry;

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

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() &&
           path.compare(path.size() - value.size(), value.size(), ending) == 0;
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

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

// HTSJDK's VCF encoder emits site QUAL with two decimal places.  Keep the
// cohort/posterior calculation in double precision, then apply the output
// contract's positive half-up rounding before handing the value to HTSlib.
// Unlike PL, GenotypingEngine does not impose a 999 cap on site QUAL: a
// multi-sample cohort can legitimately exceed it.
// Without this boundary native would expose the internal 33.4823... value as
// QUAL=33.4823 while GATK writes QUAL=33.48, even when the underlying model is
// otherwise identical.
float gatk_qual_output(double value) {
    if (!std::isfinite(value)) return 0.0F;
    const auto rounded = std::round(std::max(0.0, value) * 100.0) / 100.0;
    return static_cast<float>(std::min(
        rounded, static_cast<double>(std::numeric_limits<float>::max())));
}

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

#if FASTGATK_HAS_HTSLIB
float htslib_float_missing() {
    float value = 0.0F;
    const std::uint32_t bits = bcf_float_missing;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

bool has_header_field(const bcf_hdr_t* header, int header_type, const char* id) {
    return bcf_hdr_get_hrec(header, header_type, "ID", id, nullptr) != nullptr;
}
#endif

// Native GenomicsDBImport publishes a small JSON sidecar next to the
// portable input index.  Keep the reader dependency-free, but validate the
// fields that define the workspace contract before expanding it.  Older
// adapter workspaces may only contain fastgatk-inputs.tsv; those remain
// accepted for backwards compatibility.
std::size_t metadata_unsigned_field(const std::string& json, const char* key) {
    const std::string marker = std::string("\"") + key + "\":";
    const auto begin = json.find(marker);
    if (begin == std::string::npos)
        throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata is missing ") + key);
    auto cursor = begin + marker.size();
    while (cursor < json.size() && std::isspace(static_cast<unsigned char>(json[cursor]))) ++cursor;
    if (cursor == json.size() || !std::isdigit(static_cast<unsigned char>(json[cursor])))
        throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata has invalid ") + key);
    std::size_t value = 0;
    while (cursor < json.size() && std::isdigit(static_cast<unsigned char>(json[cursor]))) {
        const auto digit = static_cast<std::size_t>(json[cursor] - '0');
        if (value > (std::numeric_limits<std::size_t>::max() - digit) / 10U)
            throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata overflows ") + key);
        value = value * 10U + digit;
        ++cursor;
    }
    return value;
}

std::string metadata_string_field(const std::string& json, const char* key) {
    const std::string marker = std::string("\"") + key + "\"";
    const auto begin = json.find(marker);
    if (begin == std::string::npos)
        throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata is missing ") + key);
    auto cursor = json.find(':', begin + marker.size());
    if (cursor == std::string::npos)
        throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata has invalid ") + key);
    ++cursor;
    while (cursor < json.size() && std::isspace(static_cast<unsigned char>(json[cursor]))) ++cursor;
    if (cursor == json.size() || json[cursor] != '"')
        throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata has invalid ") + key);
    ++cursor;
    std::string value;
    while (cursor < json.size()) {
        const char ch = json[cursor++];
        if (ch == '"') return value;
        if (ch == '\\' && cursor < json.size()) value.push_back(json[cursor++]);
        else value.push_back(ch);
    }
    throw std::runtime_error(std::string("BAD_INPUT: native GenomicsDB metadata has unterminated ") + key);
}

void validate_native_workspace_metadata(const std::filesystem::path& workspace,
                                        const std::vector<std::string>& indexed_inputs) {
    const auto metadata_path = workspace / "fastgatk-workspace.json";
    std::ifstream metadata_stream(metadata_path);
    if (!metadata_stream) return;  // pre-metadata adapter workspace
    const std::string metadata((std::istreambuf_iterator<char>(metadata_stream)),
                               std::istreambuf_iterator<char>());
    if (metadata_unsigned_field(metadata, "schema_version") != 1 ||
        metadata_string_field(metadata, "backend") != "fastgatk-sparse-index")
        throw std::runtime_error("BAD_INPUT: unsupported or malformed native GenomicsDB workspace metadata: " +
                                 metadata_path.string());
    const auto expected = metadata_unsigned_field(metadata, "input_count");
    if (expected != indexed_inputs.size())
        throw std::runtime_error("BAD_INPUT: native GenomicsDB workspace input_count does not match fastgatk-inputs.tsv: " +
                                 workspace.string());
    for (const auto& input : indexed_inputs) {
        if (input.rfind("http://", 0) == 0 || input.rfind("https://", 0) == 0 ||
            input.rfind("s3://", 0) == 0 || input.rfind("gs://", 0) == 0) {
            throw std::runtime_error("BACKEND_UNAVAILABLE: native GenomicsDB workspace contains a remote input: " + input);
        }
        std::error_code error;
        if (!std::filesystem::is_regular_file(input, error) || error)
            throw std::runtime_error("BAD_INPUT: native GenomicsDB workspace input is missing: " + input);
    }
}

// True GATK GenomicsDB workspaces contain TileDB arrays rather than the
// fastgatk-inputs.tsv sidecar used by the portable sparse workspace.  The
// optional exporter is kept outside this process's ABI: it runs as a child,
// writes one plain-VCF shard per array, and returns a manifest of those shards.
// This makes direct replacement possible without linking the legacy GATK
// libstdc++ ABI into the Kokkos/HTSlib executable.
std::vector<std::filesystem::path> genomicsdb_bridge_temporary_files;
std::uint64_t genomicsdb_bridge_sequence = 0;
bool genomicsdb_bridge_used = false;

struct GenomicsDBBridgeCleanup {
    ~GenomicsDBBridgeCleanup() {
        for (const auto& path : genomicsdb_bridge_temporary_files) {
            std::error_code error;
            std::filesystem::remove(path, error);
        }
    }
};

GenomicsDBBridgeCleanup genomicsdb_bridge_cleanup;

std::string genomicsdb_bridge_executable() {
    if (const char* configured = std::getenv("FASTGATK_GENOMICSDB_BRIDGE");
        configured != nullptr && *configured != '\0')
        return configured;
#ifdef FASTGATK_GENOMICSDB_BRIDGE_DEFAULT
    return FASTGATK_GENOMICSDB_BRIDGE_DEFAULT;
#else
    return {};
#endif
}

bool has_regular_file(const std::filesystem::path& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error;
}

bool looks_like_true_genomicsdb_workspace(const std::filesystem::path& workspace) {
    std::error_code error;
    if (!std::filesystem::is_directory(workspace, error) || error) return false;
    if (!(has_regular_file(workspace / "callset.json") ||
          has_regular_file(workspace / "callsets.json"))) return false;
    if (!(has_regular_file(workspace / "vidmap.json") ||
          has_regular_file(workspace / "vid_mapping.json"))) return false;
    for (const auto& entry : std::filesystem::directory_iterator(workspace, error)) {
        if (error) return false;
        std::error_code schema_error;
        if (entry.is_directory(schema_error) && !schema_error &&
            has_regular_file(entry.path() / "__array_schema.tdb")) return true;
    }
    return false;
}

#if defined(__unix__) || defined(__APPLE__)
int run_genomicsdb_bridge(const std::vector<std::string>& command) {
    std::vector<char*> argv;
    argv.reserve(command.size() + 1);
    for (const auto& value : command) argv.push_back(const_cast<char*>(value.c_str()));
    argv.push_back(nullptr);
    const auto pid = fork();
    if (pid < 0) throw std::runtime_error("BACKEND_UNAVAILABLE: cannot fork GenomicsDB bridge");
    if (pid == 0) {
        execv(argv[0], argv.data());
        execvp(argv[0], argv.data());
        _exit(127);
    }
    int status = 0;
    if (waitpid(pid, &status, 0) < 0)
        throw std::runtime_error("BACKEND_UNAVAILABLE: cannot wait for GenomicsDB bridge");
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    if (WIFSIGNALED(status)) return 128 + WTERMSIG(status);
    return 1;
}
#else
int run_genomicsdb_bridge(const std::vector<std::string>&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: true GenomicsDB bridge requires POSIX process support");
}
#endif

std::vector<std::string> export_true_genomicsdb_workspace(
    const std::filesystem::path& workspace, const std::vector<std::string>& regions) {
    const auto executable = genomicsdb_bridge_executable();
    if (executable.empty())
        throw std::runtime_error(
            "BACKEND_UNAVAILABLE: true GenomicsDB workspace requires fastgatk-genomicsdb-export; "
            "set FASTGATK_GENOMICSDB_BRIDGE or build the optional bridge");
    std::error_code executable_error;
    if (!std::filesystem::is_regular_file(executable, executable_error) || executable_error)
        throw std::runtime_error("BACKEND_UNAVAILABLE: GenomicsDB bridge executable is missing: " + executable);
    const auto sequence = ++genomicsdb_bridge_sequence;
#if defined(__unix__) || defined(__APPLE__)
    const auto pid = std::to_string(static_cast<unsigned long long>(getpid()));
#else
    const auto pid = std::string("0");
#endif
    const auto manifest = std::filesystem::temp_directory_path() /
        ("fastgatk-genomicsdb-export-" + pid + "-" + std::to_string(sequence) + ".tsv");
    genomicsdb_bridge_temporary_files.push_back(manifest);
    std::vector<std::string> command{executable, workspace.string(), manifest.string()};
    // The bridge translates canonical contig:start-end selectors into the
    // workspace's global TileDB column range. Interval-list/BED selectors are
    // intentionally omitted here; HTSlib still applies them after the bridge
    // export, so unsupported selector syntax never widens a query silently.
    for (const auto& region : regions) {
        if (region.find(':') == std::string::npos || region.find('-', region.find(':') + 1) == std::string::npos)
            continue;
        command.emplace_back("--interval");
        command.push_back(region);
    }
    const auto status = run_genomicsdb_bridge(command);
    if (status != 0)
        throw std::runtime_error("BACKEND_UNAVAILABLE: true GenomicsDB bridge exited with status " +
                                 std::to_string(status) + " for " + workspace.string());
    std::ifstream stream(manifest);
    if (!stream)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GenomicsDB bridge did not publish a manifest: " +
                                 manifest.string());
    std::string line;
    std::vector<std::string> outputs;
    while (std::getline(stream, line)) {
        if (line.empty() || line.front() == '#') continue;
        const auto delimiter = line.find('\t');
        const auto output = delimiter == std::string::npos ? line : line.substr(delimiter + 1);
        if (output.empty())
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: malformed GenomicsDB bridge manifest: " +
                                     manifest.string());
        std::error_code output_error;
        if (!std::filesystem::is_regular_file(output, output_error) || output_error ||
            std::filesystem::file_size(output, output_error) == 0 || output_error)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GenomicsDB bridge output is missing: " + output);
        outputs.push_back(output);
        genomicsdb_bridge_temporary_files.emplace_back(output);
#if FASTGATK_HAS_HTSLIB
        const auto index = output + ".tbi";
        if (tbx_index_build3(output.c_str(), index.c_str(), 0, 0, &tbx_conf_vcf) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot index GenomicsDB bridge output: " + output);
        genomicsdb_bridge_temporary_files.emplace_back(index);
#endif
    }
    if (outputs.empty())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GenomicsDB bridge manifest is empty: " +
                                 manifest.string());
    genomicsdb_bridge_used = true;
    return outputs;
}

std::vector<std::string> expand_genomicsdb_inputs(
    const std::vector<std::string>& inputs, const std::vector<std::string>& regions) {
    std::vector<std::string> expanded;
    for (const auto& input : inputs) {
        constexpr const char* prefix = "gendb://";
        if (input.rfind(prefix, 0) != 0) {
            expanded.push_back(input);
            continue;
        }
        const auto workspace = std::filesystem::path(input.substr(std::strlen(prefix)));
        if (workspace.empty())
            throw std::invalid_argument("BAD_INPUT: empty gendb:// workspace");
        auto index = workspace / "fastgatk-inputs.tsv";
        // Native GenomicsDBImport publishes a self-contained materialized
        // index.  Prefer it when present; retaining the legacy index fallback
        // keeps pre-materialization adapter workspaces readable.
        const auto native_index = workspace / "fastgatk-native-inputs.tsv";
        std::error_code native_index_error;
        if (std::filesystem::is_regular_file(native_index, native_index_error) &&
            !native_index_error)
            index = native_index;
        std::ifstream stream(index);
        if (!stream)
        {
            if (looks_like_true_genomicsdb_workspace(workspace)) {
                const auto exported = export_true_genomicsdb_workspace(workspace, regions);
                expanded.insert(expanded.end(), exported.begin(), exported.end());
                continue;
            }
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: opaque GenomicsDB workspace requires the GATK GenotypeGVCFs backend: " +
                workspace.string());
        }
        std::string line;
        std::size_t count = 0;
        std::vector<std::string> workspace_inputs;
        while (std::getline(stream, line)) {
            if (line.empty() || line.front() == '#') continue;
            expanded.push_back(line);
            workspace_inputs.push_back(line);
            ++count;
        }
        if (count == 0)
            throw std::runtime_error("BAD_INPUT: gendb:// workspace input index is empty: " +
                                     workspace.string());
        // Validate the native sidecar only when present.  This preserves the
        // legacy index-only bridge while preventing a partially published or
        // cross-backend workspace from being silently interpreted as VCF.
        validate_native_workspace_metadata(workspace, workspace_inputs);
    }
    if (expanded.empty()) throw std::invalid_argument("--variant/-V is required");
    return expanded;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-genotype-gvcf (GATK-compatible native contract subset)\n"
                         "  -V, --variant FILE              input VCF/GVCF (repeatable)\n"
                         "  -R, --reference FILE            reference FASTA (accepted)\n"
                         "  -L, --intervals REGION          contig:start-end\n"
                         "  -XL, --exclude-intervals REGION exclude genomic intervals\n"
                         "      --interval-set-rule RULE   UNION (default) or INTERSECTION\n"
                         "  -O, --output FILE               output VCF/VCF.GZ\n"
                         "      --output-manifest FILE      OutputManifest JSON\n"
                         "      --heterozygosity FLOAT      SNP heterozygosity (default 1e-3)\n"
                         "      --indel-heterozygosity FLOAT\n"
                         "                                     indel heterozygosity (default 1/8000)\n"
                         "      --heterozygosity-stdev FLOAT  AF prior stdev (default 0.01)\n"
                         "      --standard-min-confidence-threshold-for-calling FLOAT\n"
                         "                                     ALT posterior threshold (default 30)\n"
                         "      --max-alternate-alleles INTEGER\n"
                         "                                     maximum proper ALT alleles (default 6)\n"
                         "      --annotate-with-num-discovered-alleles [true|false]\n"
                         "                                     emit pre-subset ALT count (default false)\n"
                         "      --use-new-qual-calculator    GATK cohort AF calculator (default)\n"
                         "      --no-use-new-qual-calculator disable cohort AF calculator\n"
                         "      --use-posteriors-to-calculate-qual\n"
                         "                                     posterior QUAL (GATK --gp-qual alias)\n"
                         "      --genotype-assignment-method MODE\n"
                         "                                     PREFER_PLS (default), USE_PLS_TO_ASSIGN,\n"
                         "                                     USE_POSTERIOR_PROBABILITIES, SET_TO_NO_CALL,\n"
                         "                                     SET_TO_NO_CALL_NO_ANNOTATIONS, BEST_MATCH_TO_ORIGINAL,\n"
                         "                                     USE_POSTERIORS_ANNOTATION, DO_NOT_ASSIGN_GENOTYPES\n"
                         "      --include-non-variant-sites  emit REF-only reference-confidence sites\n"
                         "      --only-output-calls-starting-in-intervals [true|false]\n"
                         "                                     emit finalized calls whose POS is inside -L\n"
                         "      --gatk-compatible-annotations  emit GATK standard INFO/FORMAT set\n"
                         "      --stream-by-locus            bounded k-way joint-locus merge\n"
                         "      --no-genotype-priors         flat prior for posterior QUAL\n";
            std::exit(0);
        } else if (is_option(argument, "--variant") || argument == "-V")
            options.inputs.push_back(require_value(i, argc, argv, argument, "--variant", "-V"));
        else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(i, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--intervals") || is_option(argument, "--interval") ||
                 is_option(argument, "--region") || argument == "-L") {
            const char* interval_option = argument == "-L" ? "--intervals" :
                argument.rfind("--intervals", 0) == 0 ? "--intervals" :
                argument.rfind("--interval", 0) == 0 ? "--interval" : "--region";
            options.regions.push_back(require_value(i, argc, argv, argument, interval_option, "-L"));
        }
        else if (is_option(argument, "--exclude-intervals") || argument == "-XL")
            options.excluded_regions.push_back(require_value(
                i, argc, argv, argument, "--exclude-intervals", "-XL"));
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(i, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else
                throw std::invalid_argument("BAD_INPUT: --interval-set-rule must be UNION or INTERSECTION");
        }
        else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(i, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(i, argc, argv, argument,
                argument == "--manifest" ? "--manifest" : "--output-manifest");
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--create-output-variant-index");
        } else if (is_option(argument, "--heterozygosity")) {
            options.snp_heterozygosity = std::stod(require_value(
                i, argc, argv, argument, "--heterozygosity"));
        } else if (is_option(argument, "--indel-heterozygosity")) {
            options.indel_heterozygosity = std::stod(require_value(
                i, argc, argv, argument, "--indel-heterozygosity"));
        } else if (is_option(argument, "--heterozygosity-stdev")) {
            options.heterozygosity_stdev = std::stod(require_value(
                i, argc, argv, argument, "--heterozygosity-stdev"));
        } else if (is_option(argument, "--standard-min-confidence-threshold-for-calling") ||
                   is_option(argument, "--stand-call-conf")) {
            const char* name = argument.rfind("--stand-call-conf", 0) == 0
                ? "--stand-call-conf" : "--standard-min-confidence-threshold-for-calling";
            options.standard_confidence_for_calling = std::stod(require_value(
                i, argc, argv, argument, name));
        } else if (is_option(argument, "--max-alternate-alleles")) {
            options.max_alternate_alleles = std::stoi(require_value(
                i, argc, argv, argument, "--max-alternate-alleles"));
        } else if (argument == "--annotate-with-num-discovered-alleles" ||
                   argument.rfind("--annotate-with-num-discovered-alleles=", 0) == 0) {
            options.annotate_with_num_discovered_alleles = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--annotate-with-num-discovered-alleles");
        } else if (argument == "--no-genotype-priors") {
            options.use_genotype_priors = false;
        } else if (argument == "--include-non-variant-sites") {
            options.include_non_variant_sites = true;
        } else if (argument == "--only-output-calls-starting-in-intervals" ||
                   argument.rfind("--only-output-calls-starting-in-intervals=", 0) == 0) {
            options.only_output_calls_starting_in_intervals =
                fastgatk::native::parse_optional_boolean(
                    i, argc, argv, argument,
                    "--only-output-calls-starting-in-intervals");
        } else if (argument == "--gatk-compatible-annotations" ||
                   argument == "--strict-gatk-annotations") {
            options.gatk_annotation_compatibility = true;
        } else if (argument == "--stream-by-locus" || argument == "--stream-loci") {
            options.stream_by_locus = true;
        } else if (argument == "--no-use-new-qual-calculator") {
            options.use_new_qual_calculator = false;
        } else if (argument == "--use-new-qual-calculator" ||
                   argument.rfind("--use-new-qual-calculator=", 0) == 0 ||
                   argument == "--new-qual" || argument.rfind("--new-qual=", 0) == 0) {
            const char* name = argument.rfind("--new-qual", 0) == 0
                ? "--new-qual" : "--use-new-qual-calculator";
            options.use_new_qual_calculator = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, name);
        } else if (argument == "--use-posteriors-to-calculate-qual" ||
                   argument.rfind("--use-posteriors-to-calculate-qual=", 0) == 0 ||
                   argument == "--gp-qual" || argument.rfind("--gp-qual=", 0) == 0) {
            const char* name = argument.rfind("--gp-qual", 0) == 0
                ? "--gp-qual" : "--use-posteriors-to-calculate-qual";
            options.use_posteriors_to_calculate_qual =
                fastgatk::native::parse_optional_boolean(i, argc, argv, argument, name);
        } else if (is_option(argument, "--genotype-assignment-method") ||
                   argument == "--gam") {
            options.genotype_assignment_method = require_value(
                i, argc, argv, argument, "--genotype-assignment-method", "--gam");
            std::transform(options.genotype_assignment_method.begin(),
                           options.genotype_assignment_method.end(),
                           options.genotype_assignment_method.begin(),
                           [](unsigned char value) { return static_cast<char>(std::toupper(value)); });
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Accepted launcher-compatible flags; the native VCF reader still
            // validates that all records have a header contig.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(i, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("--variant/-V is required");
    if (options.output.empty()) throw std::invalid_argument("--output/-O is required");
    if (!(options.snp_heterozygosity > 0.0 && options.snp_heterozygosity <= 1.0) ||
        !(options.indel_heterozygosity > 0.0 && options.indel_heterozygosity <= 1.0) ||
        !(options.heterozygosity_stdev > 0.0 && std::isfinite(options.heterozygosity_stdev)))
        throw std::invalid_argument("heterozygosity values must be in (0,1]");
    if (!(options.standard_confidence_for_calling >= 0.0) ||
        !std::isfinite(options.standard_confidence_for_calling))
        throw std::invalid_argument("standard confidence threshold must be finite and non-negative");
    if (options.max_alternate_alleles <= 0)
        throw std::invalid_argument("max alternate alleles must be positive");
    if (options.only_output_calls_starting_in_intervals && options.regions.empty())
        throw std::invalid_argument(
            "BAD_INPUT: --only-output-calls-starting-in-intervals requires -L/--intervals");
    const std::vector<std::string> supported_assignment_methods{
        "PREFER_PLS", "USE_PLS_TO_ASSIGN", "USE_POSTERIOR_PROBABILITIES",
        "SET_TO_NO_CALL", "SET_TO_NO_CALL_NO_ANNOTATIONS", "BEST_MATCH_TO_ORIGINAL",
        "USE_POSTERIORS_ANNOTATION",
        "DO_NOT_ASSIGN_GENOTYPES"};
    if (std::find(supported_assignment_methods.begin(), supported_assignment_methods.end(),
                  options.genotype_assignment_method) == supported_assignment_methods.end())
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: genotype-assignment-method must be one of " +
            std::string("PREFER_PLS, USE_PLS_TO_ASSIGN, USE_POSTERIOR_PROBABILITIES, SET_TO_NO_CALL, SET_TO_NO_CALL_NO_ANNOTATIONS, BEST_MATCH_TO_ORIGINAL, USE_POSTERIORS_ANNOTATION, DO_NOT_ASSIGN_GENOTYPES"));
    if (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
        !options.use_new_qual_calculator)
        throw std::invalid_argument(
            "UNSUPPORTED_PARAMETER: USE_POSTERIOR_PROBABILITIES requires the cohort AF calculator");
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct Region {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
    // Numeric contig IDs are local to each VCF header.  Keep the canonical
    // name from the first selector header so the same interval remains
    // correct when a multi-shard workspace presents disjoint dictionaries.
    std::string contig;
};

struct NativeRecordRange {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::uint64_t records = 0;
};

using NativeRecordIndex = std::map<std::string, std::vector<NativeRecordRange>>;

NativeRecordIndex load_native_record_index(const std::filesystem::path& workspace) {
    NativeRecordIndex result;
    const auto path = workspace / "fastgatk-record-index.tsv";
    std::ifstream stream(path);
    if (!stream) return result;  // legacy native workspace without span index
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::istringstream fields(line);
        std::string input;
        NativeRecordRange range;
        if (!std::getline(fields, input, '\t') || input.empty() ||
            !std::getline(fields, range.contig, '\t') || range.contig.empty())
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB record index: " + path.string());
        std::string start, end, records;
        if (!std::getline(fields, start, '\t') || !std::getline(fields, end, '\t') ||
            !std::getline(fields, records))
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB record index row: " + line);
        try {
            std::size_t consumed = 0;
            range.start = std::stoll(start, &consumed);
            if (consumed != start.size()) throw std::invalid_argument("start");
            consumed = 0;
            range.end = std::stoll(end, &consumed);
            if (consumed != end.size()) throw std::invalid_argument("end");
            consumed = 0;
            range.records = std::stoull(records, &consumed);
            if (consumed != records.size() || range.start < 0 || range.end <= range.start ||
                range.records == 0) throw std::invalid_argument("range");
        } catch (const std::exception&) {
            throw std::runtime_error("BAD_INPUT: malformed native GenomicsDB record index row: " + line);
        }
        result[input].push_back(std::move(range));
    }
    if (stream.bad())
        throw std::runtime_error("BAD_INPUT: failed reading native GenomicsDB record index: " + path.string());
    return result;
}

bool native_record_index_overlaps(const std::vector<NativeRecordRange>& ranges,
                                  const bcf_hdr_t* header,
                                  const std::vector<Region>& regions) {
    for (const auto& range : ranges) {
        for (const auto& region : regions) {
            // Region numeric IDs are header-local.  A native workspace can
            // legitimately contain shards with disjoint dictionaries (for
            // example one shard only has chr1 and another only chr2), so
            // comparing RIDs across headers would make an unrelated shard
            // look overlapping.  Resolve the query RID back to its contig
            // name before intersecting the name-keyed sidecar span.
            const auto query_contig = region.contig.empty() && region.rid >= 0
                ? bcf_hdr_id2name(header, region.rid) : region.contig.c_str();
            if (query_contig != nullptr && range.contig == query_contig &&
                range.start < region.end && range.end > region.begin)
                return true;
        }
    }
    return false;
}

// A gVCF workspace may contain either BGZF VCF/BCF plus a CSI index or a
// tabix-indexed VCF.GZ.  Keep the index wrapper local to the Host traversal
// boundary so the downstream staging/merge code remains independent of the
// on-disk index namespace.  In particular, TBI sequence IDs are not required
// to match the VCF header's rid values and must be translated by contig name.
struct GenotypeTraversalIndex {
    hts_idx_t* index = nullptr;
    tbx_t* tabix = nullptr;

    GenotypeTraversalIndex() = default;
    GenotypeTraversalIndex(const GenotypeTraversalIndex&) = delete;
    GenotypeTraversalIndex& operator=(const GenotypeTraversalIndex&) = delete;
    ~GenotypeTraversalIndex() {
        if (tabix != nullptr) tbx_destroy(tabix);
        else if (index != nullptr) hts_idx_destroy(index);
    }

    int tid_for(const bcf_hdr_t* header, int rid, const std::string& contig = {}) const {
        // CSI stores numeric sequence IDs in the order of the current VCF
        // header.  A multi-shard workspace may present disjoint dictionaries,
        // so the query RID from the first shard cannot be reused blindly.
        // Resolve by canonical contig name whenever one is available.
        if (tabix == nullptr) {
            if (!contig.empty()) return bcf_hdr_name2id(header, contig.c_str());
            return rid;
        }
        const auto* name = !contig.empty() ? contig.c_str() :
            (rid >= 0 ? bcf_hdr_id2name(header, rid) : nullptr);
        return name == nullptr ? -1 : tbx_name2id(tabix, name);
    }
};

std::unique_ptr<GenotypeTraversalIndex> load_genotype_traversal_index(
    const std::string& filename) {
    const std::filesystem::path csi_path = filename + ".csi";
    const std::filesystem::path tbi_path = filename + ".tbi";
    if (std::filesystem::exists(csi_path)) {
        if (auto* index = hts_idx_load3(filename.c_str(), csi_path.c_str(),
                                         HTS_FMT_CSI, HTS_IDX_SILENT_FAIL)) {
            auto result = std::make_unique<GenotypeTraversalIndex>();
            result->index = index;
            return result;
        }
    }
    if (std::filesystem::exists(tbi_path)) {
        if (auto* tabix = tbx_index_load2(filename.c_str(), tbi_path.c_str())) {
            auto result = std::make_unique<GenotypeTraversalIndex>();
            result->index = tabix->idx;
            result->tabix = tabix;
            return result;
        }
    }
    return {};
}

void append_region_selector(const std::string& selector, const bcf_hdr_t* header,
                            std::vector<Region>& regions,
                            fastgatk::io::IntervalFileStats& stats) {
    std::vector<fastgatk::io::IndexedInterval> parsed;
    // GenotypeGVCFs may consume a native workspace whose shards have
    // disjoint contig dictionaries.  Preserve selectors absent from the
    // first shard and resolve them by name against each input header below.
    fastgatk::io::append_interval_selector(selector, header, parsed, stats, true);
    for (const auto& interval : parsed) {
        const auto* contig = interval.contig.empty() && interval.rid >= 0
            ? bcf_hdr_id2name(header, interval.rid) : interval.contig.c_str();
        regions.push_back(Region{interval.rid, interval.begin, interval.end,
                                 contig == nullptr ? std::string{} : std::string(contig)});
    }
}

void normalize_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        if (left.contig != right.contig) return left.contig < right.contig;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.contig.empty() || region.end <= region.begin) continue;
        if (!merged.empty() && merged.back().contig == region.contig &&
            region.begin <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, region.end);
        } else {
            merged.push_back(region);
        }
    }
    regions.swap(merged);
}

// Apply GATK's interval-set rule across repeatable -L selectors.  INTERSECTION
// is a set operation between selector groups (not merely a post-hoc filter on
// the union), so each selector is parsed into a temporary normalized set and
// intersected with the accumulated half-open intervals.  Keeping contig names
// in the Region object avoids header-local RID mismatches for multi-shard VCFs.
void append_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    intersection.reserve(std::min(regions.size(), incoming.size()));
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].contig < incoming[right].contig) {
            ++left;
            continue;
        }
        if (incoming[right].contig < regions[left].contig) {
            ++right;
            continue;
        }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end) {
            Region overlap = regions[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (regions[left].end < incoming[right].end) ++left;
        else ++right;
    }
    regions.swap(intersection);
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, static_cast<int>(values[0]));
    free(values);
    return end;
}

bool overlaps_regions(const bcf_hdr_t* header, bcf1_t* record,
                      const std::vector<Region>& regions) {
    if (regions.empty()) return false;
    const auto record_end = record_end_exclusive(header, record);
    const auto* record_contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        const auto same_contig = !region.contig.empty()
            ? (record_contig != nullptr && region.contig == record_contig)
            : record->rid == region.rid;
        return same_contig && record_end > region.begin && record->pos < region.end;
    });
}

bool starts_in_regions(const bcf_hdr_t* header, const bcf1_t* record,
                       const std::vector<Region>& regions) {
    if (regions.empty()) return false;
    const auto* record_contig = record->rid >= 0
        ? bcf_hdr_id2name(header, record->rid) : nullptr;
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        const auto same_contig = !region.contig.empty()
            ? (record_contig != nullptr && region.contig == record_contig)
            : record->rid == region.rid;
        return same_contig && record->pos >= region.begin && record->pos < region.end;
    });
}

// Include intervals are evaluated using the existing GVCF span-overlap
// predicate, then exclusions are applied as a second set operation.  This is
// the same ordering as GATK's IntervalArgumentCollection and matters for a
// reference-confidence block whose END crosses an excluded coordinate.
bool in_regions(const bcf_hdr_t* header, bcf1_t* record,
                const std::vector<Region>& regions,
                const std::vector<Region>& excluded_regions = {}) {
    if (!regions.empty() && !overlaps_regions(header, record, regions)) return false;
    return !overlaps_regions(header, record, excluded_regions);
}

void merge_contig_dictionary(bcf_hdr_t* output_header, const bcf_hdr_t* input_header) {
    if (output_header == nullptr || input_header == nullptr) return;
    const auto contig_length = [](const bcf_hdr_t* header, const char* name) -> std::int64_t {
        if (header == nullptr || name == nullptr) return 0;
        const auto* hrec = bcf_hdr_get_hrec(header, BCF_HL_CTG, "ID", name, nullptr);
        if (hrec == nullptr) return 0;
        for (int index = 0; index < hrec->nkeys; ++index) {
            if (hrec->keys[index] == nullptr || hrec->vals[index] == nullptr ||
                std::strcmp(hrec->keys[index], "length") != 0)
                continue;
            char* end = nullptr;
            const auto value = std::strtoll(hrec->vals[index], &end, 10);
            if (end != hrec->vals[index] && *end == '\0' && value > 0) return value;
        }
        return 0;
    };
    for (int rid = 0; rid < input_header->n[BCF_DT_CTG]; ++rid) {
        const char* name = bcf_hdr_int2id(input_header, BCF_DT_CTG, rid);
        if (name == nullptr || *name == '\0') continue;
        const int existing = bcf_hdr_name2id(output_header, name);
        const auto input_length = contig_length(input_header, name);
        if (existing >= 0) {
            const auto output_length = contig_length(output_header, name);
            if (input_length > 0 && output_length > 0 && input_length != output_length)
                throw std::runtime_error(
                    std::string("BAD_INPUT: contig length mismatch across GenotypeGVCFs shards: ") + name);
            continue;
        }
        std::ostringstream line;
        line << "##contig=<ID=" << name;
        if (input_length > 0) line << ",length=" << input_length;
        line << ">";
        if (bcf_hdr_append(output_header, line.str().c_str()) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge VCF contig header");
    }
    if (bcf_hdr_sync(output_header) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync merged VCF contig header");
}

struct Record {
    bcf1_t* value = nullptr;
    int rid = -1;
    int pos = -1;
    std::string key;
    // REF followed by the concrete ALT alleles represented by this staged
    // record.  Keeping the names (rather than only integer indices) lets the
    // joint merge form an ALT union across shards before remapping FORMAT
    // fields.
    std::vector<std::string> alleles;
    std::vector<int> output_samples;
    int ploidy = 0;
    std::vector<int32_t> gt;
    // The merged source genotype in the union allele space, captured before
    // GenotypingEngine's PL-derived assignment replaces `gt`.  GATK's
    // PREFER_PLS path falls back to GATKVariantContextUtils.bestMatchToOriginalGT()
    // on this call, so it must survive the union/sample merge.  `source_alleles`
    // is the union allele list that `source_gt` indexes into; every record of a
    // merged group shares it, so the merge keeps the first one.
    std::vector<int32_t> source_gt;
    std::vector<std::string> source_alleles;
    std::vector<int32_t> dp;
    std::vector<int32_t> gq;
    // HaplotypeCaller's physical-phasing annotations are ordinary
    // genotype-level attributes to GenotypeGVCFs.  Preserve them through the
    // Host-side ALT union/sample merge instead of leaving whichever source
    // bcf1_t happened to be first.  The Kokkos PL/AF/GT numerical kernels do
    // not consume these strings; they remain at the GATK GenotypeBuilder
    // compatibility boundary.
    std::vector<std::string> pgt;
    std::vector<std::string> pid;
    std::vector<int32_t> ps;
    // GATK moves the genotype quality to RGQ for monomorphic reference-only
    // records emitted by --include-non-variant-sites.  Keep this separate
    // from the ordinary variant GQ so the native diagnostic profile can
    // retain its historical fields while the compatibility profile matches
    // GenotypeGVCFs' writer contract.
    std::vector<int32_t> rgq;
    std::vector<int32_t> min_dp;
    std::vector<int32_t> ad;
    // FORMAT/SB is the four-entry ref-forward/ref-reverse/alt-forward/
    // alt-reverse table emitted by HaplotypeCaller.  Keep it in the staged
    // representation so joint-shard merges can reproduce GATK's site-level
    // FS/SOR annotations instead of accidentally using only the first shard.
    std::vector<int32_t> sb;
    std::vector<int32_t> pl;
    // GATK PHRED_SCALED_POSTERIORS_KEY (FORMAT/PP), stored in the same
    // sample-major Number=G layout as PL for posterior-based assignment.
    std::vector<int32_t> pp;
    // True only when this native pass materialized GP/PG from the explicit
    // USE_POSTERIOR_PROBABILITIES assignment path. GATK's --gp-qual option is
    // conditional: with default PL assignment and no posterior field it must
    // leave the PL-derived QUAL unchanged.
    bool posterior_quality_available = false;
    int allele_count = 0;
    bool reference_block = false;
    // Preserve the original REF/<NON_REF> likelihood rows while a
    // reference-confidence block is compacted to a REF-only output record.
    // The compacted PL=[0] is sufficient for the VCF writer, but cannot be
    // used to compute a true cross-sample reference posterior.  Keeping this
    // bounded source-width copy lets the post-merge Kokkos kernel combine all
    // contributing samples before the diagnostic RCQ/RCP fields are emitted.
    std::vector<int32_t> reference_confidence_pl;
    int reference_confidence_allele_count = 0;
    // If an input genotype used an ALT removed from the joint union (currently
    // the bounded orphan '*' cleanup), GATK preserves the likelihood fields
    // but publishes that sample as an explicit no-call.  Keep this state
    // through the PL-derived assignment boundary.
    std::vector<std::uint8_t> force_no_call_samples;
    bool orphan_spanning_deletion = false;
    // True only for a locus this tool SYNTHESIZED because a spanning record
    // covers it (dense mode only).  Such a locus has no input record of its
    // own, which is what makes its site confidence exactly zero rather than a
    // ~1e-16 residue, and that in turn makes GATK's QualByDepth numerator the
    // sign-preserving -0.0 (see update_gatk_standard_annotations).
    bool materialized_spanning_locus = false;
    // True when the source FORMAT carried GQ for this record, i.e. when GATK's
    // merged genotype would report hasGQ().  Distinct from the native PL-derived
    // GQ vector, which is synthesized even for inputs without FORMAT/GQ.
    bool source_has_gq = false;
    // Set when --include-non-variant-sites keeps a locus whose alternate
    // alleles were all pruned.  GATK emits the resulting REF-only no-call
    // straight out of GenotypeGVCFsEngine.regenotypeVC, without any further
    // site or standard annotation, so the shared annotation stages must be
    // skipped for that record.
    bool finalized_monomorphic_ref = false;
    // AFCalculator returns per-ALT posterior absence probabilities.  Keep
    // this bounded vector through the staging boundary so the final output
    // allele subset can apply GATK's standard confidence threshold after
    // the cohort merge, before INFO/FORMAT annotations are finalized.
    std::vector<double> cohort_log10_p_allele_absent;
    std::vector<int32_t> cohort_integer_allele_counts;
    // AlleleFrequencyResult.log10_p_no_variant.  GenotypingEngine reports the
    // *complementary* posterior as the site QUAL for a locus that turned
    // monomorphic (MathUtils.log10OneMinusPow10), which is only reachable in
    // --include-non-variant-sites mode.  Keep the raw posterior alongside the
    // already-derived `cohort.qual` so that boundary can be reproduced.
    double cohort_log10_p_no_variant = 0.0;
    bool cohort_quality_available = false;
    // GenotypingEngine.calculateGenotypes() tests the RECOMPUTED site
    // confidence, not the source FILTER and not the published QUAL token:
    //
    //   final double phredScaledConfidence = (-10.0 * log10Confidence) + 0.0;
    //   ...
    //   if ( ! passesCallThreshold(phredScaledConfidence) ) {
    //       builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
    //   }
    //
    // (GenotypingEngine.java:163 and :184-186, with passesCallThreshold()
    // defined as `conf >= standardConfidenceForCalling` at :430-432).  The same
    // double is handed to builder.log10PError() at :183, but htsjdk renders
    // QUAL with two decimals, so the token can read at or above the cutoff
    // while the filter is applied: measured on the weak-locus fixture, the
    // confidence lies inside [23.135, 23.1358) while the published token is
    // 23.14.  Keep the raw double so the writer boundary can reproduce the
    // decision exactly, and keep the PRE-posterior one, because :184-186 runs
    // before the optional posterior QUAL update at :192-199.
    // `call_confidence_available` separates "the genotyping engine computed a
    // confidence" from "the record never reached it" -- a passthrough
    // reference-confidence block has a missing QUAL and GATK applies no filter
    // to it.
    double call_confidence = 0.0;
    bool call_confidence_available = false;
};

// The joint-materialization pass is Host/HTSlib heavy, but its final
// per-locus annotation and serialization pass is also a natural portable
// pipeline boundary.  Keep the record itself move-only across the queues so
// the caller thread remains the single VCF writer while the worker invokes
// the Kokkos genotype kernels.
struct GenotypeDecoded {
    Record record;
};

struct GenotypeComputed {
    Record record;
};

struct GenotypeEncoded {
    Record record;
    std::string text;
};

std::uint64_t genotype_record_bytes(const Record& record) {
    std::uint64_t bytes = 4096;
    const auto add = [&](std::size_t count, std::size_t width) {
        if (count > (std::numeric_limits<std::uint64_t>::max() - bytes) / width)
            bytes = std::numeric_limits<std::uint64_t>::max();
        else
            bytes += static_cast<std::uint64_t>(count) * width;
    };
    add(record.alleles.size(), sizeof(std::string));
    add(record.output_samples.size(), sizeof(int));
    add(record.gt.size(), sizeof(int32_t));
    add(record.source_gt.size(), sizeof(int32_t));
    add(record.source_alleles.size(), sizeof(std::string));
    add(record.dp.size(), sizeof(int32_t));
    add(record.gq.size(), sizeof(int32_t));
    add(record.rgq.size(), sizeof(int32_t));
    add(record.min_dp.size(), sizeof(int32_t));
    add(record.ad.size(), sizeof(int32_t));
    add(record.sb.size(), sizeof(int32_t));
    add(record.pl.size(), sizeof(int32_t));
    add(record.pp.size(), sizeof(int32_t));
    add(record.reference_confidence_pl.size(), sizeof(int32_t));
    add(record.cohort_log10_p_allele_absent.size(), sizeof(double));
    add(record.cohort_integer_allele_counts.size(), sizeof(int32_t));
    add(record.key.size(), sizeof(char));
    for (const auto& allele : record.alleles) add(allele.size(), sizeof(char));
    return bytes;
}

void destroy_records(std::vector<Record>& records) {
    for (auto& record : records) bcf_destroy(record.value);
    records.clear();
}

std::string record_key(const bcf1_t* record) {
    // Records at one locus must be coalesced even when different input
    // shards carried different concrete ALT subsets.  The union/remapping is
    // performed after grouping, so ALT names do not belong in this key.
    const char* ref = record->n_allele > 0 ? record->d.allele[0] : "";
    return std::to_string(record->rid) + ":" + std::to_string(record->pos) + ":" + ref;
}

int record_span_end(const bcf_hdr_t* header, const Record& record) {
    int end = record.pos + 1;
    if (!record.alleles.empty() && record.alleles.front() != "<NON_REF>")
        end = std::max(end, record.pos + static_cast<int>(record.alleles.front().size()));
    int32_t* end_value = nullptr;
    int end_count = 0;
    const auto length = bcf_get_info_int32(header, record.value, "END",
                                           &end_value, &end_count);
    if (length > 0 && end_count > 0 && end_value[0] > end)
        end = end_value[0];  // END is 1-based inclusive, hence exclusive here.
    free(end_value);
    return end;
}

struct VariantSpan {
    int rid = -1;
    int begin = -1;
    int end = -1;
};

// GATK's per-traversal record of the deletions it has EMITTED, i.e.
// ``GenotypingEngine.upstreamDeletionsLoc`` (GenotypingEngine.java:52), a
// PriorityQueue<Locatable> ordered by interval end whose only production
// mutation sites are ``recordDeletions()`` (:343-357) and its lazy cull
// (:344-346).  ``clearUpstreamDeletionsLoc()`` (:330-332) is
// ``@VisibleForTesting`` and is never called from production code, so the state
// accumulates across the whole ordered walk.
//
// ``recordDeletions(vc, emittedAlleles)`` runs at :179, i.e. AFTER the output
// allele subset was decided at :155, and appends
// ``new SimpleInterval(contig, vc.getStart(), vc.getStart() + deletionSize)``
// with ``deletionSize = vc.getReference().length() - allele.length()`` for every
// EMITTED allele with a positive deletion size.  Three consequences the previous
// input-derived span list could not reproduce:
//   * a deletion allele of a locus GATK dropped (or pruned) is never recorded;
//   * the recorded interval ends at ``start + (refLen - altLen)``, not at
//     ``start + refLen``;
//   * the symbolic '*' is an ordinary allele here -- its length is 1 -- so a
//     locus whose emitted set contains a '*' records an interval of its own.
// ``deletionSize`` (and hence the recorded interval) is invariant under
// ``GATKVariantContextUtils.reverseTrimAlleles()``, which only clips trailing
// bases, so the untrimmed merged record gives the same interval GATK computes.
//
// ``isVcCoveredByDeletion()`` (:365-371) then requires
//     loc.getStart() < vc.getStart() && vc.getStart() <= loc.getEnd()
// which in native's 0-based, half-open coordinates is
//     span.begin < record.pos && record.pos < span.end
// with ``end = pos + (refLen - altLen) + 1``.  The start test is STRICT, so the
// very locus being genotyped can never be covered by its own deletion.
struct EmittedDeletions {
    // Min-heap by interval end: GATK's PriorityQueue iteration order aside, the
    // heap is what makes the lazy cull cheap while the coverage test scans it.
    std::vector<VariantSpan> heap;

    static bool later_end(const VariantSpan& left, const VariantSpan& right) {
        return left.end > right.end;
    }

    // GenotypingEngine.java:344-346: drop the head while it cannot cover the
    // current locus any more (or sits on another contig).
    void cull(int rid, int pos) {
        while (!heap.empty() &&
               (heap.front().rid != rid || heap.front().end <= pos)) {
            std::pop_heap(heap.begin(), heap.end(), later_end);
            heap.pop_back();
        }
    }

    bool covers(int rid, int pos) const {
        for (const auto& span : heap)
            if (span.rid == rid && span.begin < pos && pos < span.end) return true;
        return false;
    }

    // recordDeletions() over the alleles a locus actually emitted.
    void record(const Record& record, const std::vector<std::string>& emitted) {
        if (record.rid < 0 || record.alleles.empty()) return;
        cull(record.rid, record.pos);
        const int reference_length = static_cast<int>(record.alleles.front().size());
        for (const auto& allele : emitted) {
            const int deletion_size = reference_length - static_cast<int>(allele.size());
            if (deletion_size <= 0) continue;
            heap.push_back(VariantSpan{record.rid, record.pos,
                                       record.pos + deletion_size + 1});
            std::push_heap(heap.begin(), heap.end(), later_end);
        }
    }
};

std::string reference_base(faidx_t* reference_index, const bcf_hdr_t* header,
                           int rid, int position) {
    if (reference_index == nullptr) return "N";
    const char* contig = bcf_hdr_id2name(header, rid);
    if (contig == nullptr) return "N";
    int fetched_length = 0;
    char* fetched = faidx_fetch_seq(reference_index, contig, position, position,
                                    &fetched_length);
    std::string result = fetched != nullptr && fetched_length == 1
        ? std::string(1, fetched[0]) : std::string("N");
    free(fetched);
    return result;
}

[[maybe_unused]] int genotype_index(int first, int second) {
    if (first > second) std::swap(first, second);
    return second * (second + 1) / 2 + first;
}

std::size_t genotype_width(int allele_count, int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t result = 1;
    for (int step = 1; step <= ploidy; ++step) {
        const auto numerator = static_cast<std::size_t>(allele_count + step - 1);
        if (result > std::numeric_limits<std::size_t>::max() / numerator) return 0;
        result *= numerator;
        result /= static_cast<std::size_t>(step);
    }
    return result;
}

void unrank_genotype(std::size_t index, int allele_count, int ploidy,
                     std::vector<int>& alleles) {
    alleles.assign(static_cast<std::size_t>(ploidy), 0);
    int available = allele_count;
    for (int position = ploidy - 1; position >= 0; --position) {
        int selected = available - 1;
        for (int candidate = 0; candidate < available; ++candidate) {
            const auto group = genotype_width(candidate + 1, position);
            if (index < group) {
                selected = candidate;
                break;
            }
            index -= group;
        }
        alleles[static_cast<std::size_t>(position)] = selected;
        available = selected + 1;
    }
}

std::size_t rank_genotype(const std::vector<int>& alleles, int allele_count) {
    (void)allele_count;
    std::size_t rank = 0;
    for (int position = static_cast<int>(alleles.size()) - 1; position >= 0; --position) {
        const int selected = alleles[static_cast<std::size_t>(position)];
        for (int candidate = 0; candidate < selected; ++candidate)
            rank += genotype_width(candidate + 1, position);
    }
    return rank;
}

[[maybe_unused]] void compact_genotype_fields(const bcf_hdr_t* input_header,
                             const bcf_hdr_t* output_header,
                             bcf1_t* record,
                             int selected,
                             int original_alleles) {
    const int sample_count = output_header->n[BCF_DT_SAMPLE];
    if (sample_count <= 0) return;

    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "GT") >= 0) {
        int32_t* genotypes = nullptr;
        int genotype_count = 0;
        const auto length = bcf_get_genotypes(input_header, record, &genotypes, &genotype_count);
        if (length > 0 && genotype_count >= sample_count && genotype_count % sample_count == 0) {
            const int ploidy = genotype_count / sample_count;
            for (int i = 0; i < genotype_count; ++i) {
                if (bcf_gt_is_missing(genotypes[i]) || genotypes[i] == bcf_int32_vector_end) {
                    genotypes[i] = bcf_gt_missing;
                    continue;
                }
                const int allele = bcf_gt_allele(genotypes[i]);
                const bool phased = bcf_gt_is_phased(genotypes[i]);
                if (allele == 0) genotypes[i] = phased ? bcf_gt_phased(0) : bcf_gt_unphased(0);
                else if (allele == selected) genotypes[i] = phased ? bcf_gt_phased(1) : bcf_gt_unphased(1);
                else genotypes[i] = bcf_gt_missing;
            }
            if (bcf_update_genotypes(output_header, record, genotypes,
                                     genotype_count) != 0) {
                free(genotypes);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact GT FORMAT field");
            }
            (void)ploidy;
        }
        free(genotypes);
    }

    const auto compact_per_sample = [&](const char* tag, int old_width, int new_width,
                                        auto index_for) {
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0) return;
        int32_t* values = nullptr;
        int value_count = 0;
        const auto length = bcf_get_format_int32(input_header, record, tag, &values, &value_count);
        if (length <= 0 || old_width <= 0 || value_count < sample_count * old_width ||
            value_count % sample_count != 0) {
            free(values);
            return;
        }
        std::vector<int32_t> compact(static_cast<std::size_t>(sample_count) * new_width,
                                     bcf_int32_missing);
        for (int sample = 0; sample < sample_count; ++sample) {
            const auto source = values + sample * old_width;
            auto destination = compact.data() + sample * new_width;
            for (int i = 0; i < new_width; ++i) {
                const int source_index = index_for(i);
                if (source_index >= 0 && source_index < old_width)
                    destination[i] = source[source_index];
            }
        }
        if (bcf_update_format_int32(output_header, record, tag, compact.data(),
                                    static_cast<int>(compact.size())) != 0) {
            free(values);
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot compact ") + tag + " FORMAT field");
        }
        free(values);
    };
    compact_per_sample("AD", original_alleles, 2,
                       [&](int index) { return index == 0 ? 0 : selected; });
    const int old_pl_width = original_alleles * (original_alleles + 1) / 2;
    const int zero_alt = selected * (selected + 1) / 2;
    const int hom_alt = zero_alt + selected;
    compact_per_sample("PL", old_pl_width, 3,
                       [&](int index) { return index == 0 ? 0 : index == 1 ? zero_alt : hom_alt; });
}

void extract_materialized_fields(const bcf_hdr_t* input_header, const bcf1_t* record,
                                 const std::vector<int>& selected, int original_alleles,
                                 const std::vector<int>& output_samples,
                                 Record& destination) {
    const int source_samples = input_header->n[BCF_DT_SAMPLE];
    destination.allele_count = static_cast<int>(selected.size()) + 1;
    if (source_samples <= 0 || static_cast<int>(output_samples.size()) != source_samples) return;
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto genotype_length = bcf_get_genotypes(input_header, const_cast<bcf1_t*>(record),
                                                    &genotypes, &genotype_count);
    if (genotype_length > 0 && genotype_count % source_samples == 0) {
        destination.ploidy = genotype_count / source_samples;
        destination.gt.assign(static_cast<std::size_t>(source_samples) * destination.ploidy,
                              bcf_gt_missing);
        for (int sample = 0; sample < source_samples; ++sample) {
            for (int allele_index = 0; allele_index < destination.ploidy; ++allele_index) {
                const auto source = genotypes[sample * destination.ploidy + allele_index];
                if (source == bcf_int32_vector_end) continue;
                // CombineGVCFs represents a phased no-call as .|.  Its
                // second missing allele carries HTSlib's phase bit, and
                // GenotypeBuilder retains that bit when it later assigns the
                // PL-derived call.  Do not collapse it to ./.: the physical
                // PGT/PID/PS attributes alone are not enough for the final
                // GenotypeGVCFs GT writer to recover the separator.
                if (bcf_gt_is_missing(source)) {
                    destination.gt[static_cast<std::size_t>(sample) * destination.ploidy +
                                   allele_index] = bcf_gt_is_phased(source)
                        ? bcf_gt_phased(-1) : bcf_gt_missing;
                    continue;
                }
                const auto allele = bcf_gt_allele(source);
                const auto phased = bcf_gt_is_phased(source);
                int output_allele = -1;
                if (allele == 0) output_allele = 0;
                else {
                    const auto found = std::find(selected.begin(), selected.end(), allele);
                    if (found != selected.end())
                        output_allele = static_cast<int>(found - selected.begin()) + 1;
                }
                if (output_allele >= 0)
                    destination.gt[static_cast<std::size_t>(sample) * destination.ploidy + allele_index] =
                        phased ? bcf_gt_phased(output_allele) : bcf_gt_unphased(output_allele);
            }
        }
    }
    free(genotypes);

    const auto read_scalar_format = [&](const char* tag, std::vector<int32_t>& destination) {
        int32_t* scalar = nullptr;
        int scalar_count = 0;
        const auto length = bcf_get_format_int32(input_header, const_cast<bcf1_t*>(record),
                                                 tag, &scalar, &scalar_count);
        if (length > 0 && scalar_count % source_samples == 0) {
            const int width = scalar_count / source_samples;
            if (width > 0) {
                destination.assign(static_cast<std::size_t>(source_samples), bcf_int32_missing);
                for (int sample = 0; sample < source_samples; ++sample)
                    destination[static_cast<std::size_t>(sample)] = scalar[sample * width];
            }
        }
        free(scalar);
    };
    read_scalar_format("DP", destination.dp);
    read_scalar_format("GQ", destination.gq);
    // Whether the source FORMAT actually carried GQ is a distinct contract from
    // the later PL-derived GQ: GATK's cleanupGenotypeAnnotations() branches on
    // the merged genotype's hasGQ().  read_scalar_format leaves the vector
    // untouched when the tag is absent, so emptiness records that here.
    destination.source_has_gq = !destination.gq.empty();
    read_scalar_format("MIN_DP", destination.min_dp);
    read_scalar_format("PS", destination.ps);

    const auto read_string_format = [&](const char* tag,
                                        std::vector<std::string>& output) {
        char** values = nullptr;
        int value_count = 0;
        const auto length = bcf_get_format_string(input_header,
                                                  const_cast<bcf1_t*>(record),
                                                  tag, &values, &value_count);
        if (length > 0 && values != nullptr && value_count >= source_samples) {
            output.assign(static_cast<std::size_t>(source_samples), ".");
            for (int sample = 0; sample < source_samples; ++sample) {
                if (values[sample] != nullptr && values[sample][0] != '\0')
                    output[static_cast<std::size_t>(sample)] = values[sample];
            }
        }
        // HTSlib allocates FORMAT strings as one contiguous block followed
        // by the pointer table.  This is the ownership rule used throughout
        // the native VCF tools.
        if (values != nullptr) {
            if (value_count > 0 && values[0] != nullptr) free(values[0]);
            free(values);
        }
    };
    read_string_format("PGT", destination.pgt);
    read_string_format("PID", destination.pid);

    // FORMAT/SB is Number=4 in GATK gVCFs.  Some legacy inputs encode it as a
    // wider vector; retain the first four entries per sample and fail closed
    // for malformed/missing rows rather than inventing strand evidence.
    int32_t* values = nullptr;
    int value_count = 0;
    const auto sb_length = bcf_get_format_int32(input_header,
                                                const_cast<bcf1_t*>(record),
                                                "SB", &values, &value_count);
    if (sb_length > 0 && value_count % source_samples == 0) {
        const int old_width = value_count / source_samples;
        if (old_width >= 4) {
            destination.sb.assign(static_cast<std::size_t>(source_samples) * 4,
                                  bcf_int32_missing);
            for (int sample = 0; sample < source_samples; ++sample) {
                const auto source = values + sample * old_width;
                auto* target = destination.sb.data() +
                    static_cast<std::size_t>(sample) * 4;
                for (int index = 0; index < 4; ++index)
                    target[index] = source[index];
            }
        }
    }
    free(values);

    values = nullptr;
    value_count = 0;
    const auto ad_length = bcf_get_format_int32(input_header, const_cast<bcf1_t*>(record),
                                                "AD", &values, &value_count);
    if (ad_length > 0 && value_count % source_samples == 0) {
        const int old_width = value_count / source_samples;
        destination.ad.assign(static_cast<std::size_t>(source_samples) * destination.allele_count,
                              bcf_int32_missing);
        for (int sample = 0; sample < source_samples; ++sample) {
            const auto source = values + sample * old_width;
            auto* destination_sample = destination.ad.data() +
                static_cast<std::size_t>(sample) * destination.allele_count;
            if (old_width > 0) destination_sample[0] = source[0];
            for (std::size_t index = 0; index < selected.size(); ++index) {
                const int source_index = selected[index];
                if (source_index >= 0 && source_index < old_width)
                    destination_sample[index + 1] = source[source_index];
            }
        }
    }
    free(values);

    values = nullptr;
    value_count = 0;
    const auto pl_length = bcf_get_format_int32(input_header, const_cast<bcf1_t*>(record),
                                                "PL", &values, &value_count);
    if (pl_length > 0 && value_count % source_samples == 0) {
        const int old_width = value_count / source_samples;
        const auto new_width = genotype_width(destination.allele_count, destination.ploidy);
        if (new_width == 0) {
            free(values);
            throw std::runtime_error("BAD_INPUT: unsupported genotype PL dimensions");
        }
        destination.pl.assign(static_cast<std::size_t>(source_samples) * new_width,
                              bcf_int32_missing);
        for (int sample = 0; sample < source_samples; ++sample) {
            const auto source = values + sample * old_width;
            auto* destination_sample = destination.pl.data() +
                static_cast<std::size_t>(sample) * new_width;
            std::vector<int> output_genotype;
            std::vector<int> source_genotype;
            output_genotype.reserve(destination.ploidy);
            source_genotype.reserve(destination.ploidy);
            for (std::size_t output_index = 0; output_index < new_width; ++output_index) {
                unrank_genotype(output_index, destination.allele_count, destination.ploidy,
                                output_genotype);
                source_genotype = output_genotype;
                for (auto& allele : source_genotype)
                    allele = allele == 0 ? 0 : selected[static_cast<std::size_t>(allele - 1)];
                const auto source_index = rank_genotype(source_genotype, original_alleles);
                if (source_index < static_cast<std::size_t>(old_width))
                    destination_sample[output_index] = source[source_index];
            }
        }
    }
    free(values);

    // FORMAT/PP is GATK's phred-scaled posterior annotation.  It follows the
    // same VCF Number=G rank as PL, so project it through the identical
    // source-allele tuple mapping before the joint ALT union is formed.
    values = nullptr;
    value_count = 0;
    const auto pp_length = bcf_get_format_int32(input_header, const_cast<bcf1_t*>(record),
                                                "PP", &values, &value_count);
    if (pp_length > 0 && value_count % source_samples == 0) {
        const int old_width = value_count / source_samples;
        const auto new_width = genotype_width(destination.allele_count, destination.ploidy);
        if (new_width == 0) {
            free(values);
            throw std::runtime_error("BAD_INPUT: unsupported genotype PP dimensions");
        }
        destination.pp.assign(static_cast<std::size_t>(source_samples) * new_width,
                              bcf_int32_missing);
        for (int sample = 0; sample < source_samples; ++sample) {
            const auto source = values + sample * old_width;
            auto* destination_sample = destination.pp.data() +
                static_cast<std::size_t>(sample) * new_width;
            std::vector<int> output_genotype;
            std::vector<int> source_genotype;
            output_genotype.reserve(destination.ploidy);
            source_genotype.reserve(destination.ploidy);
            for (std::size_t output_index = 0; output_index < new_width; ++output_index) {
                unrank_genotype(output_index, destination.allele_count, destination.ploidy,
                                output_genotype);
                source_genotype = output_genotype;
                for (auto& allele : source_genotype)
                    allele = allele == 0 ? 0 : selected[static_cast<std::size_t>(allele - 1)];
                const auto source_index = rank_genotype(source_genotype, original_alleles);
                if (source_index < static_cast<std::size_t>(old_width))
                    destination_sample[output_index] = source[source_index];
            }
        }
    }
    free(values);
    destination.output_samples = output_samples;
    (void)original_alleles;
}

std::string allele_list_string(const std::vector<std::string>& alleles) {
    std::string value;
    for (std::size_t index = 0; index < alleles.size(); ++index) {
        if (index != 0) value.push_back(',');
        value += alleles[index].empty() ? "N" : alleles[index];
    }
    return value;
}

enum class PriorAlleleType { Ref, Snp, Indel, Other };

std::vector<double> build_log10_genotype_priors(const std::vector<std::string>& alleles,
                                                int ploidy,
                                                double snp_heterozygosity,
                                                double indel_heterozygosity,
                                                bool use_priors) {
    const auto width = genotype_width(static_cast<int>(alleles.size()), ploidy);
    if (width == 0) throw std::runtime_error("BAD_INPUT: unsupported genotype prior dimensions");
    std::vector<double> priors(width, 0.0);
    if (!use_priors) return priors;
    if (alleles.empty()) throw std::runtime_error("BAD_INPUT: empty allele list for genotype priors");
    const auto& reference = alleles.front();
    const auto log10_snp = std::log10(snp_heterozygosity);
    const auto log10_indel = std::log10(indel_heterozygosity);
    // This is GenotypePriorCalculator.assumingHW from GATK. SNP priors are
    // divided across the three possible non-reference bases; indels and
    // symbolic/other alleles use their event-type prior directly.
    constexpr double log10_snp_normalization = 0.47712125471966244; // log10(3)
    std::vector<PriorAlleleType> types;
    types.reserve(alleles.size());
    types.push_back(PriorAlleleType::Ref);
    for (std::size_t index = 1; index < alleles.size(); ++index) {
        const auto& allele = alleles[index];
        if (allele == "*" || (!allele.empty() && allele.front() == '<')) {
            types.push_back(PriorAlleleType::Other);
        } else if (allele.size() == reference.size()) {
            types.push_back(PriorAlleleType::Snp);
        } else {
            types.push_back(PriorAlleleType::Indel);
        }
    }
    std::vector<double> het(alleles.size(), 0.0);
    std::vector<double> hom(alleles.size(), 0.0);
    for (std::size_t index = 1; index < alleles.size(); ++index) {
        switch (types[index]) {
            case PriorAlleleType::Snp:
                het[index] = log10_snp - log10_snp_normalization;
                hom[index] = 2.0 * log10_snp - log10_snp_normalization;
                break;
            case PriorAlleleType::Indel:
                het[index] = log10_indel;
                hom[index] = 2.0 * log10_indel;
                break;
            case PriorAlleleType::Other: {
                const auto log10_other = std::max(log10_snp, log10_indel);
                het[index] = log10_other;
                hom[index] = 2.0 * log10_other;
                break;
            }
            case PriorAlleleType::Ref:
                break;
        }
    }
    const auto calculated = fastgatk::kernels::calculate_genotype_priors_kokkos(
        het, hom, ploidy);
    ++genotype_kernel_telemetry.genotype_prior_calls;
    genotype_kernel_telemetry.genotype_prior_prepare_seconds += calculated.prepare_seconds;
    genotype_kernel_telemetry.genotype_prior_execute_seconds += calculated.seconds;
    genotype_kernel_telemetry.genotype_prior_execution_space = calculated.execution_space;
    return calculated.log10_priors;
}

void derive_gt_gq_from_pl(Record& record) {
    if (record.pl.empty() || record.ploidy <= 0 || record.allele_count < 2 ||
        record.output_samples.empty()) return;
    const auto derived = fastgatk::kernels::derive_genotype_gt_gq_kokkos(
        record.pl, record.output_samples.size(), record.allele_count, record.ploidy);
    ++genotype_kernel_telemetry.calls;
    genotype_kernel_telemetry.prepare_seconds += derived.prepare_seconds;
    genotype_kernel_telemetry.execute_seconds += derived.seconds;
    genotype_kernel_telemetry.execution_space = derived.execution_space;
    // GenotypeBuilder starts from the input genotype, then substitutes the
    // PL-derived allele list.  Its existing phased flag therefore survives
    // the call.  Preserve the corresponding HTSlib phase separators before
    // replacing the staged GT array; PGT/PID/PS are carried independently
    // below.
    const auto input_gt = record.gt;
    const bool input_phase_shape = input_gt.size() ==
        record.output_samples.size() * static_cast<std::size_t>(record.ploidy);
    record.gt.assign(record.output_samples.size() * static_cast<std::size_t>(record.ploidy),
                     bcf_gt_missing);
    record.gq.assign(record.output_samples.size(), bcf_int32_missing);
    for (std::size_t sample = 0; sample < record.output_samples.size(); ++sample) {
        bool missing = false;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto allele = derived.alleles[sample * static_cast<std::size_t>(record.ploidy) +
                                                static_cast<std::size_t>(copy)];
            if (allele < 0) {
                missing = true;
                break;
            }
            const auto index = sample * static_cast<std::size_t>(record.ploidy) +
                               static_cast<std::size_t>(copy);
            const bool phased = input_phase_shape && bcf_gt_is_phased(input_gt[index]);
            record.gt[index] = phased ? bcf_gt_phased(allele) : bcf_gt_unphased(allele);
        }
        if (missing) continue;
        if (derived.gq[sample] >= 0) record.gq[sample] = derived.gq[sample];
        if (sample < record.force_no_call_samples.size() &&
            record.force_no_call_samples[sample] != 0) {
            for (int copy = 0; copy < record.ploidy; ++copy)
                record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                          static_cast<std::size_t>(copy)] = bcf_gt_missing;
            record.gq[sample] = 0;
        }
    }
}

void apply_genotype_assignment(const bcf_hdr_t* output_header, Record& record,
                               const Options& options,
                               const std::vector<double>& posterior_priors) {
    // GenotypeGVCFs' Java engine clones its arguments and unconditionally
    // forces PREFER_PLS in createMinimalArgs().  The native diagnostic profile
    // keeps the richer assignment modes for focused kernel tests, while the
    // GATK-compatible writer must preserve that effective Java behavior.
    const std::string effective_method = options.gatk_annotation_compatibility
        ? "PREFER_PLS" : options.genotype_assignment_method;
    const auto& method = effective_method;
    if (method == "PREFER_PLS" || method == "USE_PLS_TO_ASSIGN" ||
        method == "DO_NOT_ASSIGN_GENOTYPES") return;
    if (method == "SET_TO_NO_CALL_NO_ANNOTATIONS") {
        if (record.ploidy > 0 && record.gt.empty() && !record.pl.empty()) {
            const auto width = genotype_width(record.allele_count, record.ploidy);
            if (width != 0 && record.pl.size() % width == 0)
                record.gt.assign((record.pl.size() / width) *
                                     static_cast<std::size_t>(record.ploidy),
                                 bcf_gt_missing);
        } else if (record.ploidy > 0 && !record.gt.empty()) {
            std::fill(record.gt.begin(), record.gt.end(), bcf_gt_missing);
        }
        if (!record.gt.empty() && bcf_update_genotypes(
                output_header, record.value, record.gt.data(),
                static_cast<int>(record.gt.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write no-call GTs");
        // Java's GenotypeBuilder.noAnnotations() clears every genotype-level
        // attribute, not only the fields materialized by Record (GQ/DP/AD/PL/
        // PP). Enumerate the FORMAT dictionary and remove fields present on
        // this record, retaining GT so the no-call remains representable.
        bcf_unpack(record.value, BCF_UN_FMT);
        std::vector<std::string> format_ids;
        format_ids.reserve(static_cast<std::size_t>(output_header->nhrec));
        for (int index = 0; index < output_header->nhrec; ++index) {
            const auto* hrec = output_header->hrec[index];
            if (hrec == nullptr || hrec->type != BCF_HL_FMT) continue;
            for (int key = 0; key < hrec->nkeys; ++key) {
                if (std::strcmp(hrec->keys[key], "ID") != 0 || hrec->vals[key] == nullptr) continue;
                const std::string id(hrec->vals[key]);
                if (id != "GT" && std::find(format_ids.begin(), format_ids.end(), id) == format_ids.end())
                    format_ids.push_back(id);
            }
        }
        const auto has_format = [&](const std::string& id) {
            const auto key = bcf_hdr_id2int(output_header, BCF_DT_ID, id.c_str());
            if (key < 0) return false;
            for (int index = 0; index < record.value->n_fmt; ++index)
                if (record.value->d.fmt[index].id == key) return true;
            return false;
        };
        for (const auto& id : format_ids) {
            if (!has_format(id)) continue;
            if (bcf_update_format(output_header, record.value, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove " + id + " FORMAT field");
            ++genotype_kernel_telemetry.no_annotation_format_fields;
        }
        record.gq.clear();
        record.ad.clear();
        record.pl.clear();
        record.pp.clear();
        record.dp.clear();
        return;
    }
    if (method == "BEST_MATCH_TO_ORIGINAL") {
        // The staged GT is the original call after allele-name projection.
        // GATK only converts an original GQ=0 hom-ref-like call to no-call
        // when PL[0] is also zero; otherwise it preserves the best match.
        const auto width = record.ploidy > 0
            ? genotype_width(record.allele_count, record.ploidy) : 0;
        const auto sample_count = record.gt.empty() || record.ploidy <= 0
            ? 0 : record.gt.size() / static_cast<std::size_t>(record.ploidy);
        if (sample_count > 0 && record.gq.size() == sample_count) {
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                const bool pl_zero = record.pl.empty() ||
                    (width != 0 && record.pl.size() >= (sample + 1) * width &&
                     record.pl[sample * width] == 0);
                if (record.gq[sample] != 0 || !pl_zero) continue;
                for (int copy = 0; copy < record.ploidy; ++copy)
                    record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                              static_cast<std::size_t>(copy)] = bcf_gt_missing;
            }
        }
        if (!record.gt.empty() && bcf_update_genotypes(
                output_header, record.value, record.gt.data(),
                static_cast<int>(record.gt.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write best-match GTs");
        return;
    }
    if (method == "USE_POSTERIORS_ANNOTATION") {
        if (record.ploidy <= 0 || record.allele_count < 2 || record.pp.empty())
            throw std::runtime_error(
                "BAD_INPUT: USE_POSTERIORS_ANNOTATION requires FORMAT/PP for every merged record");
        const auto width = genotype_width(record.allele_count, record.ploidy);
        if (width == 0 || record.pp.size() % width != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid merged PP dimensions for genotype assignment");
        const auto sample_count = record.pp.size() / width;
        const auto derived = fastgatk::kernels::derive_genotype_gt_gq_kokkos(
            record.pp, sample_count, record.allele_count, record.ploidy);
        ++genotype_kernel_telemetry.posterior_annotation_assignment_calls;
        genotype_kernel_telemetry.posterior_annotation_assignment_prepare_seconds += derived.prepare_seconds;
        genotype_kernel_telemetry.posterior_annotation_assignment_execute_seconds += derived.seconds;
        genotype_kernel_telemetry.posterior_annotation_assignment_execution_space = derived.execution_space;
        genotype_kernel_telemetry.posterior_annotation_assignment_samples += sample_count;
        const auto input_gt = record.gt;
        const bool input_phase_shape = input_gt.size() ==
            sample_count * static_cast<std::size_t>(record.ploidy);
        record.gt.assign(sample_count * static_cast<std::size_t>(record.ploidy), bcf_gt_missing);
        record.gq.assign(sample_count, bcf_int32_missing);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            bool missing = false;
            for (int copy = 0; copy < record.ploidy; ++copy) {
                const auto allele = derived.alleles[
                    sample * static_cast<std::size_t>(record.ploidy) +
                    static_cast<std::size_t>(copy)];
                if (allele < 0) {
                    missing = true;
                    break;
                }
                const auto index = sample * static_cast<std::size_t>(record.ploidy) +
                                   static_cast<std::size_t>(copy);
                const bool phased = input_phase_shape && bcf_gt_is_phased(input_gt[index]);
                record.gt[index] = phased ? bcf_gt_phased(allele) : bcf_gt_unphased(allele);
            }
            if (!missing && derived.gq[sample] >= 0)
                record.gq[sample] = derived.gq[sample];
        }
        if (bcf_update_genotypes(output_header, record.value, record.gt.data(),
                                 static_cast<int>(record.gt.size())) != 0 ||
            bcf_update_format_int32(output_header, record.value, "GQ", record.gq.data(),
                                    static_cast<int>(record.gq.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write posterior-annotation genotype fields");
        return;
    }
    if (method == "SET_TO_NO_CALL" && record.ploidy > 0 && !record.gt.empty()) {
        record.gt.assign(record.gt.size(), bcf_gt_missing);
        if (bcf_update_genotypes(output_header, record.value, record.gt.data(),
                                 static_cast<int>(record.gt.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write no-call GTs");
        return;
    }
    if (record.ploidy <= 0 || record.allele_count < 2 || record.pl.empty()) return;
    const auto width = genotype_width(record.allele_count, record.ploidy);
    if (width == 0 || record.pl.size() % width != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid merged PL dimensions for genotype assignment");
    const auto sample_count = record.pl.size() / width;
    if (method == "SET_TO_NO_CALL") {
        record.gt.assign(sample_count * static_cast<std::size_t>(record.ploidy), bcf_gt_missing);
        if (bcf_update_genotypes(output_header, record.value, record.gt.data(),
                                 static_cast<int>(record.gt.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write no-call GTs");
        return;
    }
    if (method != "USE_POSTERIOR_PROBABILITIES")
        throw std::runtime_error("UNSUPPORTED_PARAMETER: genotype assignment method is not implemented");
    if (posterior_priors.size() != width)
        throw std::runtime_error(
            "BACKEND_UNAVAILABLE: cohort posterior genotype priors were not materialized");
    const auto derived = fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
        record.pl, sample_count, record.allele_count, record.ploidy, posterior_priors);
    ++genotype_kernel_telemetry.posterior_assignment_calls;
    genotype_kernel_telemetry.posterior_assignment_prepare_seconds += derived.prepare_seconds;
    genotype_kernel_telemetry.posterior_assignment_execute_seconds += derived.seconds;
    genotype_kernel_telemetry.posterior_assignment_execution_space = derived.execution_space;
    genotype_kernel_telemetry.posterior_assignment_samples += sample_count;
    const auto input_gt = record.gt;
    const bool input_phase_shape = input_gt.size() ==
        sample_count * static_cast<std::size_t>(record.ploidy);
    record.gt.assign(sample_count * static_cast<std::size_t>(record.ploidy), bcf_gt_missing);
    record.gq.assign(sample_count, bcf_int32_missing);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        bool missing = false;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto allele = derived.alleles[
                sample * static_cast<std::size_t>(record.ploidy) +
                static_cast<std::size_t>(copy)];
            if (allele < 0) {
                missing = true;
                break;
            }
            const auto index = sample * static_cast<std::size_t>(record.ploidy) +
                               static_cast<std::size_t>(copy);
            const bool phased = input_phase_shape && bcf_gt_is_phased(input_gt[index]);
            record.gt[index] = phased ? bcf_gt_phased(allele) : bcf_gt_unphased(allele);
        }
        if (!missing && derived.gq[sample] >= 0)
            record.gq[sample] = derived.gq[sample];
    }
    if (derived.posterior_phred.size() != sample_count * width ||
        derived.prior_phred.size() != width)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: posterior GP/PG dimensions are incomplete");
    const float missing_float = htslib_float_missing();
    std::vector<float> gp(sample_count * width, missing_float);
    std::vector<float> pg(sample_count * width, missing_float);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        for (std::size_t genotype = 0; genotype < width; ++genotype) {
            const auto index = sample * width + genotype;
            if (derived.posterior_phred[index] >= 0.0)
                gp[index] = static_cast<float>(derived.posterior_phred[index]);
            if (derived.prior_phred[genotype] >= 0.0)
                pg[index] = static_cast<float>(derived.prior_phred[genotype]);
        }
    }
    if (bcf_update_genotypes(output_header, record.value, record.gt.data(),
                             static_cast<int>(record.gt.size())) != 0 ||
        bcf_update_format_int32(output_header, record.value, "GQ", record.gq.data(),
                                static_cast<int>(record.gq.size())) != 0 ||
        bcf_update_format_float(output_header, record.value, "GP", gp.data(),
                                static_cast<int>(gp.size())) != 0 ||
        bcf_update_format_float(output_header, record.value, "PG", pg.data(),
                                static_cast<int>(pg.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write posterior genotype fields");
    record.posterior_quality_available = true;
}

void account_pl_remap(const fastgatk::kernels::RemapGenotypePlResult& result) {
    ++genotype_kernel_telemetry.pl_remap_calls;
    genotype_kernel_telemetry.pl_remap_prepare_seconds += result.prepare_seconds;
    genotype_kernel_telemetry.pl_remap_execute_seconds += result.seconds;
    genotype_kernel_telemetry.pl_remap_execution_space = result.execution_space;
}

void account_allele_field_remap(const fastgatk::kernels::RemapAlleleFieldResult& result) {
    ++genotype_kernel_telemetry.allele_field_remap_calls;
    genotype_kernel_telemetry.allele_field_remap_prepare_seconds += result.prepare_seconds;
    genotype_kernel_telemetry.allele_field_remap_execute_seconds += result.seconds;
    genotype_kernel_telemetry.allele_field_remap_execution_space = result.execution_space;
}

void normalize_haplotype_caller_phasing(const bcf_hdr_t* output_header,
                                        Record& record) {
    // Direct translation of GenotypeGVCFsEngine.cleanupGenotypeAnnotations():
    // a GenotypeBuilder retains HaplotypeCaller's physical phase attributes,
    // and a final homozygous-variant call changes PGT from its heterozygous
    // representation to the Java constant "1|1".  GT/PGT/PID/PS are Host VCF
    // attributes; their computation must not move the Kokkos PL/AF kernels
    // out of their numerical role.
    if (record.pgt.empty() || record.ploidy <= 0 || record.gt.empty() ||
        record.output_samples.empty() ||
        bcf_hdr_id2int(output_header, BCF_DT_ID, "PGT") < 0)
        return;
    const auto sample_count = std::min(record.pgt.size(), record.output_samples.size());
    if (record.gt.size() < sample_count * static_cast<std::size_t>(record.ploidy)) return;
    bool changed = false;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        if (record.pgt[sample].empty() || record.pgt[sample] == ".") continue;
        int hom_allele = -1;
        bool hom_var = true;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto encoded = record.gt[
                sample * static_cast<std::size_t>(record.ploidy) +
                static_cast<std::size_t>(copy)];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                hom_var = false;
                break;
            }
            const auto allele = bcf_gt_allele(encoded);
            if (allele <= 0 || (hom_allele >= 0 && allele != hom_allele)) {
                hom_var = false;
                break;
            }
            hom_allele = allele;
        }
        if (hom_var && record.pgt[sample] != "1|1") {
            record.pgt[sample] = "1|1";
            changed = true;
        }
    }
    if (!changed) return;
    std::vector<const char*> values;
    values.reserve(record.pgt.size());
    for (const auto& value : record.pgt) values.push_back(value.c_str());
    if (bcf_update_format_string(output_header, record.value, "PGT", values.data(),
                                 static_cast<int>(values.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot normalize PGT FORMAT field");
}

void update_site_quality_from_posteriors(Record& record, const Options& options) {
    if (!options.use_posteriors_to_calculate_qual || !record.posterior_quality_available ||
        record.pl.empty() ||
        record.ploidy <= 0 || record.allele_count < 2) return;
    const auto width = genotype_width(record.allele_count, record.ploidy);
    if (width == 0 || record.pl.size() % width != 0) return;
    const auto sample_count = record.pl.size() / width;
    const auto priors = build_log10_genotype_priors(
        record.alleles, record.ploidy, options.snp_heterozygosity,
        options.indel_heterozygosity, options.use_genotype_priors);
    const auto spanning_deletion = std::find(record.alleles.begin(), record.alleles.end(), "*");
    const int spanning_deletion_index = spanning_deletion == record.alleles.end()
        ? -1 : static_cast<int>(std::distance(record.alleles.begin(), spanning_deletion));
    const auto posterior = fastgatk::kernels::calculate_site_posterior_kokkos(
        record.pl, sample_count, record.allele_count, record.ploidy, priors,
        spanning_deletion_index);
    ++genotype_kernel_telemetry.posterior_calls;
    genotype_kernel_telemetry.posterior_prepare_seconds += posterior.prepare_seconds;
    genotype_kernel_telemetry.posterior_execute_seconds += posterior.seconds;
    genotype_kernel_telemetry.posterior_execution_space = posterior.execution_space;
    genotype_kernel_telemetry.posterior_samples += posterior.samples_with_likelihoods;
    if (posterior.samples_with_likelihoods == 0) return;
    // GATK writes a phred-scaled site confidence. Keep the VCF float finite
    // and deterministic for extremely confident cohorts.
    record.value->qual = gatk_qual_output(posterior.qual);
}

void update_cross_sample_reference_confidence(const bcf_hdr_t* output_header,
                                              Record& record,
                                              const Options& options) {
    // Pure REF-only materialization has no alternate genotype row and is
    // already represented by the input GQ/PL.  Concrete joint records are the
    // sites where every contributing sample's REF confidence must be combined.
    const std::vector<int32_t>* posterior_pl = &record.pl;
    int posterior_allele_count = record.allele_count;
    std::vector<std::string> posterior_alleles = record.alleles;
    // A pure reference block has been compacted to REF-only before reaching
    // this stage.  Use the preserved source-width REF/<NON_REF> rows to
    // calculate a genuine cross-sample posterior rather than treating the
    // compact PL=[0] as certainty.  This branch is also what enables
    // --include-non-variant-sites to retain meaningful RCQ/RCP diagnostics.
    if ((record.pl.empty() || record.allele_count < 2) &&
        record.reference_confidence_allele_count >= 2 &&
        !record.reference_confidence_pl.empty()) {
        posterior_pl = &record.reference_confidence_pl;
        posterior_allele_count = record.reference_confidence_allele_count;
        if (posterior_alleles.size() < 2)
            posterior_alleles.push_back("<NON_REF>");
    }
    if (posterior_pl->empty() || record.ploidy <= 0 || posterior_allele_count < 2)
        return;
    const auto width = genotype_width(posterior_allele_count, record.ploidy);
    if (width == 0 || posterior_pl->size() % width != 0) return;
    const auto sample_count = posterior_pl->size() / width;
    const auto priors = build_log10_genotype_priors(
        posterior_alleles, record.ploidy, options.snp_heterozygosity,
        options.indel_heterozygosity, options.use_genotype_priors);
    const auto posterior =
        fastgatk::kernels::calculate_cross_sample_reference_confidence_kokkos(
            *posterior_pl, sample_count, posterior_allele_count, record.ploidy, priors);
    ++genotype_kernel_telemetry.cross_sample_reference_calls;
    genotype_kernel_telemetry.cross_sample_reference_prepare_seconds +=
        posterior.prepare_seconds;
    genotype_kernel_telemetry.cross_sample_reference_execute_seconds +=
        posterior.seconds;
    genotype_kernel_telemetry.cross_sample_reference_execution_space =
        posterior.execution_space;
    genotype_kernel_telemetry.cross_sample_reference_samples +=
        posterior.samples_with_likelihoods;
    const float joint_qual = static_cast<float>(std::min(999.0, posterior.joint_qual));
    const float joint_probability = static_cast<float>(std::clamp(
        std::pow(10.0, posterior.joint_log10_p_reference), 0.0, 1.0));
    if (bcf_update_info_float(output_header, record.value, "RCQ", &joint_qual, 1) != 0 ||
        bcf_update_info_float(output_header, record.value, "RCP", &joint_probability, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write cross-sample reference confidence");
}

std::vector<double> build_cohort_prior_pseudocounts(const std::vector<std::string>& alleles,
                                                    double snp_heterozygosity,
                                                    double indel_heterozygosity,
                                                    double heterozygosity_stdev) {
    if (alleles.empty() || !(heterozygosity_stdev > 0.0))
        throw std::invalid_argument("invalid cohort allele prior configuration");
    const auto reference_length = alleles.front().size();
    const auto ref_pseudocount = snp_heterozygosity /
                                 (heterozygosity_stdev * heterozygosity_stdev);
    const auto snp_pseudocount = snp_heterozygosity * ref_pseudocount;
    const auto indel_pseudocount = indel_heterozygosity * ref_pseudocount;
    std::vector<double> result(alleles.size(), indel_pseudocount);
    result[0] = ref_pseudocount;
    for (std::size_t allele = 1; allele < alleles.size(); ++allele)
        result[allele] = alleles[allele].size() == reference_length
            ? snp_pseudocount : indel_pseudocount;
    return result;
}

void update_cohort_af_annotations(const bcf_hdr_t* output_header, Record& record,
                                  const Options& options) {
    if (!options.use_new_qual_calculator || record.ploidy <= 0 ||
        record.allele_count < 2) return;
    const auto width = genotype_width(record.allele_count, record.ploidy);
    if (width == 0) return;
    std::size_t sample_count = 0;
    if (!record.pl.empty()) {
        if (record.pl.size() % width != 0) return;
        sample_count = record.pl.size() / width;
    } else if (!record.gt.empty() &&
               record.gt.size() % static_cast<std::size_t>(record.ploidy) == 0) {
        // Reblocked/reference-confidence inputs may omit FORMAT/PL entirely.
        // Keep a missing matrix so the shared Kokkos AF API can apply GATK's
        // diploid hom-ref + GQ approximation where it is legal.
        sample_count = record.gt.size() / static_cast<std::size_t>(record.ploidy);
    } else {
        return;
    }
    if (sample_count == 0) return;
    std::vector<int32_t> cohort_pl = record.pl;
    if (cohort_pl.empty())
        cohort_pl.assign(sample_count * width, bcf_int32_missing);
    std::vector<int32_t> allele_indices(record.gt.size(), -1);
    for (std::size_t index = 0; index < record.gt.size(); ++index) {
        const auto encoded = record.gt[index];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
        const auto allele = bcf_gt_allele(encoded);
        if (allele >= 0 && allele < record.allele_count) allele_indices[index] = allele;
    }
    auto prior_pseudocounts = build_cohort_prior_pseudocounts(
        record.alleles, options.snp_heterozygosity, options.indel_heterozygosity,
        options.heterozygosity_stdev);
    if (!options.use_genotype_priors)
        std::fill(prior_pseudocounts.begin(), prior_pseudocounts.end(), 1.0);
    const auto cohort = fastgatk::kernels::calculate_allele_frequency_kokkos(
        cohort_pl, sample_count, record.allele_count, record.ploidy,
        prior_pseudocounts, record.gq, allele_indices,
        [&]() {
            for (int allele = 1; allele < record.allele_count; ++allele)
                if (record.alleles[static_cast<std::size_t>(allele)] == "*") return allele;
            return -1;
        }());
    ++genotype_kernel_telemetry.cohort_calls;
    genotype_kernel_telemetry.cohort_prepare_seconds += cohort.prepare_seconds;
    genotype_kernel_telemetry.cohort_execute_seconds += cohort.seconds;
    genotype_kernel_telemetry.cohort_execution_space = cohort.execution_space;
    genotype_kernel_telemetry.cohort_samples += cohort.samples_with_likelihoods;
    genotype_kernel_telemetry.cohort_approximate_gq_samples += cohort.approximate_gq_samples;
    genotype_kernel_telemetry.cohort_iterations += static_cast<std::uint64_t>(cohort.iterations);
    if (cohort.converged) ++genotype_kernel_telemetry.cohort_converged;
    // Preserve the AFCalculator output needed by the subsequent GATK
    // output-allele subset step.  This is deliberately kept on Host as a
    // small per-ALT vector; the expensive posterior/EM work remains in the
    // shared Kokkos kernel above.
    record.cohort_log10_p_allele_absent = cohort.log10_p_allele_absent;
    record.cohort_integer_allele_counts = cohort.integer_allele_counts;
    if (cohort.samples_with_likelihoods == 0) return;

    // GATK's GenotypingEngine uses AFResult for site confidence and emits the
    // MLE allele count/frequency annotations separately from the genotype-call
    // AC/AN/AF fields.  AN is the number of called genotype alleles, not the
    // rounded EM count, so derive it from the authoritative merged GT vector.
    record.value->qual = gatk_qual_output(cohort.qual);
    record.cohort_log10_p_no_variant = cohort.log10_p_no_variant;
    record.cohort_quality_available = true;
    // The site confidence GenotypingEngine tests at :184-186 (and publishes as
    // QUAL at :183) for a non-monomorphic site is
    // -10 * AFresult.log10ProbOnlyRefAlleleExists(), i.e. the raw
    // `cohort.qual` before gatk_qual_output()'s two-decimal rendering.
    record.call_confidence = -10.0 * cohort.log10_p_no_variant;
    record.call_confidence_available = true;
    std::vector<int32_t> mle_ac(static_cast<std::size_t>(record.allele_count - 1), 0);
    for (std::size_t allele = 1; allele < cohort.integer_allele_counts.size(); ++allele)
        mle_ac[allele - 1] = cohort.integer_allele_counts[allele];
    if (bcf_update_info_int32(output_header, record.value, "MLEAC", mle_ac.data(),
                              static_cast<int>(mle_ac.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write MLEAC annotation");
    int32_t an = 0;
    if (!record.gt.empty()) {
        for (const auto encoded : record.gt)
            if (encoded != bcf_int32_vector_end && !bcf_gt_is_missing(encoded)) ++an;
    }
    if (an <= 0) {
        for (const auto count : cohort.integer_allele_counts) an += count;
    }
    std::vector<float> mle_af(mle_ac.size(), 0.0F);
    for (std::size_t allele = 0; allele < mle_ac.size(); ++allele)
        mle_af[allele] = an == 0 ? 0.0F : static_cast<float>(mle_ac[allele]) / static_cast<float>(an);
    if (bcf_update_info_float(output_header, record.value, "MLEAF", mle_af.data(),
                              static_cast<int>(mle_af.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write MLEAF annotation");
}

// HTSJDK's QualByDepth intentionally jitters raw QD values >= 35 with
// Utils.getRandomGenerator().nextGaussian().  GATK's utility uses the fixed
// release seed 47382911, so this is deterministic across JVM invocations and
// must be reproduced at the Host writer boundary for multi-sample cohorts.
// Keep the Java 48-bit LCG here instead of std::rand()/std::normal_distribution
// so the result is independent of libc and C++ standard-library versions.
class GatkJavaRandom {
public:
    explicit GatkJavaRandom(const std::uint64_t seed = 47382911ULL)
        : state_((seed ^ kMultiplier) & kMask) {}

    double next_gaussian() {
        if (have_next_gaussian_) {
            have_next_gaussian_ = false;
            return next_gaussian_;
        }
        double first = 0.0;
        double second = 0.0;
        double radius = 0.0;
        do {
            first = 2.0 * next_double() - 1.0;
            second = 2.0 * next_double() - 1.0;
            radius = first * first + second * second;
        } while (radius >= 1.0 || radius == 0.0);
        const auto multiplier = std::sqrt(-2.0 * std::log(radius) / radius);
        next_gaussian_ = second * multiplier;
        have_next_gaussian_ = true;
        return first * multiplier;
    }

private:
    static constexpr std::uint64_t kMultiplier = 0x5DEECE66DULL;
    static constexpr std::uint64_t kAddend = 0xBULL;
    static constexpr std::uint64_t kMask = (1ULL << 48) - 1ULL;

    std::uint64_t state_ = 0;
    bool have_next_gaussian_ = false;
    double next_gaussian_ = 0.0;

    std::uint32_t next_bits(const int bits) {
        state_ = (state_ * kMultiplier + kAddend) & kMask;
        return static_cast<std::uint32_t>(state_ >> (48 - bits));
    }

    double next_double() {
        const auto high = static_cast<std::uint64_t>(next_bits(26));
        const auto low = static_cast<std::uint64_t>(next_bits(27));
        return static_cast<double>((high << 27) + low) /
               static_cast<double>(1ULL << 53);
    }
};

double gatk_fix_high_qd(const double raw_qd, GatkJavaRandom& random) {
    constexpr double max_qd_before_fixing = 35.0;
    constexpr double ideal_high_qd = 30.0;
    constexpr double jitter_sigma = 3.0;
    if (raw_qd < max_qd_before_fixing) return raw_qd;
    return ideal_high_qd + random.next_gaussian() * jitter_sigma;
}

// GVCFs produced by HaplotypeCaller already carry Number=A MLEAC/MLEAF
// values for the concrete ALT plus <NON_REF>.  GenotypeGVCFs must project
// those stale per-ALT values whenever it forms the joint ALT union or removes
// <NON_REF>; otherwise --no-use-new-qual-calculator leaks a second value into
// the final VCF (for example MLEAC=1,0 instead of MLEAC=1).  Number=A excludes
// REF, so its source/target map is the allele map with the REF entry removed.
std::vector<int32_t> number_a_target_to_source(
    const std::vector<int32_t>& target_to_old,
    int target_allele_count) {
    if (target_allele_count < 2 ||
        target_to_old.size() != static_cast<std::size_t>(target_allele_count))
        throw std::invalid_argument("invalid Number=A allele projection dimensions");
    std::vector<int32_t> result(static_cast<std::size_t>(target_allele_count - 1), -1);
    for (int target = 1; target < target_allele_count; ++target) {
        const auto source = target_to_old[static_cast<std::size_t>(target)];
        result[static_cast<std::size_t>(target - 1)] = source > 0 ? source - 1 : -1;
    }
    return result;
}

void remap_info_number_a_int(const bcf_hdr_t* output_header, bcf1_t* value,
                             const char* tag, int old_allele_count,
                             int new_allele_count,
                             const std::vector<int32_t>& target_to_old) {
    if (old_allele_count < 2 || new_allele_count < 2 ||
        bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0)
        return;
    int32_t* source_values = nullptr;
    int source_count = 0;
    const auto length = bcf_get_info_int32(output_header, value, tag,
                                           &source_values, &source_count);
    if (length <= 0 || source_count < old_allele_count - 1) {
        free(source_values);
        return;
    }
    const std::vector<int32_t> source(source_values,
                                      source_values + old_allele_count - 1);
    free(source_values);
    const auto target_map = number_a_target_to_source(target_to_old, new_allele_count);
    const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
        source, 1, old_allele_count - 1, new_allele_count - 1, target_map);
    account_allele_field_remap(remapped);
    if (bcf_update_info_int32(output_header, value, tag, remapped.values.data(),
                              new_allele_count - 1) != 0)
        throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot remap INFO/") + tag);
}

void remap_info_number_a_float(const bcf_hdr_t* output_header, bcf1_t* value,
                               const char* tag, int old_allele_count,
                               int new_allele_count,
                               const std::vector<int32_t>& target_to_old) {
    if (old_allele_count < 2 || new_allele_count < 2 ||
        bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0)
        return;
    float* source_values = nullptr;
    int source_count = 0;
    const auto length = bcf_get_info_float(output_header, value, tag,
                                           &source_values, &source_count);
    if (length <= 0 || source_count < old_allele_count - 1) {
        free(source_values);
        return;
    }
    std::vector<float> remapped(static_cast<std::size_t>(new_allele_count - 1),
                                htslib_float_missing());
    for (int target = 1; target < new_allele_count; ++target) {
        const auto source = target_to_old[static_cast<std::size_t>(target)];
        if (source > 0 && source - 1 < old_allele_count - 1)
            remapped[static_cast<std::size_t>(target - 1)] = source_values[source - 1];
    }
    free(source_values);
    if (bcf_update_info_float(output_header, value, tag, remapped.data(),
                              new_allele_count - 1) != 0)
        throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot remap INFO/") + tag);
}

void remap_record_to_allele_union(const bcf_hdr_t* output_header, Record& record,
                                  const std::vector<std::string>& union_alleles,
                                  bool force_no_call_on_dropped_gt = true) {
    if (record.alleles == union_alleles) return;
    if (record.alleles.empty() || union_alleles.empty() ||
        record.allele_count != static_cast<int>(record.alleles.size()))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid staged allele metadata");

    std::vector<int> old_to_new(record.alleles.size(), -1);
    for (std::size_t old = 0; old < record.alleles.size(); ++old) {
        const auto found = std::find(union_alleles.begin(), union_alleles.end(), record.alleles[old]);
        // The same remap kernel serves both ALT-union expansion and the
        // post-AF output subset.  An old concrete ALT absent from the target
        // list is intentionally dropped in the latter case.
        if (found != union_alleles.end())
            old_to_new[old] = static_cast<int>(found - union_alleles.begin());
    }
    const int old_alleles = static_cast<int>(record.alleles.size());
    const int new_alleles = static_cast<int>(union_alleles.size());
    std::vector<int32_t> target_to_old(static_cast<std::size_t>(new_alleles), -1);
    for (int target = 0; target < new_alleles; ++target) {
        const auto found = std::find(record.alleles.begin(), record.alleles.end(),
                                     union_alleles[static_cast<std::size_t>(target)]);
        if (found != record.alleles.end())
            target_to_old[static_cast<std::size_t>(target)] =
                static_cast<int32_t>(found - record.alleles.begin());
    }

    // Project existing Number=A fields before HTSlib sees the new ALT list.
    // MLEAC/MLEAF are especially important for the legacy AF-calculator path:
    // they are inherited from HaplotypeCaller and are not overwritten when
    // --no-use-new-qual-calculator is selected.
    remap_info_number_a_int(output_header, record.value, "MLEAC",
                            old_alleles, new_alleles, target_to_old);
    remap_info_number_a_float(output_header, record.value, "MLEAF",
                              old_alleles, new_alleles, target_to_old);
    remap_info_number_a_float(output_header, record.value, "AF",
                              old_alleles, new_alleles, target_to_old);

    if (!record.gt.empty() && record.ploidy > 0) {
        const auto sample_count = record.gt.size() / static_cast<std::size_t>(record.ploidy);
        if (record.force_no_call_samples.size() < sample_count)
            record.force_no_call_samples.resize(sample_count, 0);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            bool dropped_gt_allele = false;
            for (int copy = 0; copy < record.ploidy; ++copy) {
                auto& encoded = record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                                           static_cast<std::size_t>(copy)];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            const int old_index = bcf_gt_allele(encoded);
            if (old_index < 0 || old_index >= old_alleles || old_to_new[old_index] < 0) {
                encoded = bcf_gt_missing;
                dropped_gt_allele = true;
                continue;
            }
            const int new_index = old_to_new[old_index];
            encoded = bcf_gt_is_phased(encoded) ? bcf_gt_phased(new_index)
                                                 : bcf_gt_unphased(new_index);
            }
            // ALT-union construction and max-ALT clipping preserve GATK's
            // no-call contract for a source GT that mentions a removed
            // allele.  The final AF-based output subset is different:
            // AlleleSubsettingUtils.subsetAlleles() rebuilds a PREFER_PLS
            // call from the projected Number=G row.  In particular, an
            // orphan spanning-deletion '*' may disappear from the emitted
            // ALT list while the projected triploid PLs still call 0/1/1.
            // Let that later Kokkos GT/GQ derivation replace the temporary
            // missing source GT instead of pinning the entire sample to ./.
            if (dropped_gt_allele && force_no_call_on_dropped_gt) {
                record.force_no_call_samples[sample] = 1;
                for (int copy = 0; copy < record.ploidy; ++copy)
                    record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                              static_cast<std::size_t>(copy)] = bcf_gt_missing;
                if (sample < record.gq.size()) record.gq[sample] = 0;
            }
        }
    }

    if (!record.ad.empty()) {
        const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
            record.ad, record.output_samples.size(), old_alleles, new_alleles,
            target_to_old);
        account_allele_field_remap(remapped);
        record.ad = remapped.values;
    }

    if (!record.pl.empty()) {
        const auto old_width = genotype_width(old_alleles, record.ploidy);
        const auto new_width = genotype_width(new_alleles, record.ploidy);
        if (old_width == 0 || new_width == 0)
            throw std::runtime_error("BAD_INPUT: unsupported genotype PL dimensions during ALT union");
        const auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
            record.pl, record.output_samples.size(), old_alleles, new_alleles,
            record.ploidy, target_to_old);
        account_pl_remap(remapped);
        record.pl = remapped.pl;
    }
    if (!record.pp.empty()) {
        const auto old_width = genotype_width(old_alleles, record.ploidy);
        const auto new_width = genotype_width(new_alleles, record.ploidy);
        if (old_width == 0 || new_width == 0)
            throw std::runtime_error("BAD_INPUT: unsupported genotype PP dimensions during ALT union");
        const auto remapped = fastgatk::kernels::remap_genotype_pl_kokkos(
            record.pp, record.output_samples.size(), old_alleles, new_alleles,
            record.ploidy, target_to_old);
        account_pl_remap(remapped);
        record.pp = remapped.pl;
    }

    int32_t* info_ad = nullptr;
    int info_ad_count = 0;
    const auto info_length = bcf_get_info_int32(output_header, record.value, "AD",
                                                &info_ad, &info_ad_count);
    if (info_length > 0 && info_ad_count >= old_alleles) {
        const std::vector<int32_t> source(info_ad, info_ad + old_alleles);
        const auto remapped = fastgatk::kernels::remap_allele_field_kokkos(
            source, 1, old_alleles, new_alleles, target_to_old);
        account_allele_field_remap(remapped);
        if (bcf_update_info_int32(output_header, record.value, "AD", remapped.values.data(), new_alleles) != 0) {
            free(info_ad);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remap INFO/AD field");
        }
    }
    free(info_ad);

    if (bcf_update_alleles_str(output_header, record.value,
                               allele_list_string(union_alleles).c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update joint ALT union");
    record.alleles = union_alleles;
    record.allele_count = new_alleles;
}

// GenotypingEngine applies --max-alternate-alleles before the cohort AF
// calculator.  AlleleSubsettingUtils ranks proper ALTs by the sum, over
// samples, of the PL distance between hom-ref and each sample's most likely
// genotype when that genotype contains the ALT.  Keep REF and <NON_REF>,
// retain original allele order after the stable score sort, and project all
// Number=A/Number=G fields through the common Kokkos remap path.
bool apply_gatk_max_alternate_alleles(const bcf_hdr_t* output_header,
                                      Record& record,
                                      const Options& options) {
    if (!options.gatk_annotation_compatibility ||
        options.max_alternate_alleles <= 0 || record.allele_count < 2 ||
        record.alleles.size() < 2)
        return false;
    int proper_alt_count = 0;
    for (std::size_t index = 1; index < record.alleles.size(); ++index)
        if (record.alleles[index] != "<NON_REF>") ++proper_alt_count;
    if (proper_alt_count <= options.max_alternate_alleles) return false;

    const auto width = genotype_width(record.allele_count, record.ploidy);
    const auto sample_count = width > 0 && !record.pl.empty() &&
        record.pl.size() % width == 0 ? record.pl.size() / width : 0;
    std::vector<double> scores(record.alleles.size(), 0.0);
    if (sample_count != 0) {
        const auto calculated = fastgatk::kernels::calculate_allele_likelihood_scores_kokkos(
            record.pl, sample_count, record.allele_count, record.ploidy);
        ++genotype_kernel_telemetry.max_alt_score_kernel_calls;
        genotype_kernel_telemetry.max_alt_score_prepare_seconds += calculated.prepare_seconds;
        genotype_kernel_telemetry.max_alt_score_execute_seconds += calculated.seconds;
        genotype_kernel_telemetry.max_alt_score_execution_space = calculated.execution_space;
        for (std::size_t index = 0; index < scores.size() &&
             index < calculated.scores.size(); ++index)
            scores[index] = calculated.scores[index];
    }

    std::vector<std::size_t> proper_indices;
    proper_indices.reserve(static_cast<std::size_t>(proper_alt_count));
    for (std::size_t index = 1; index < record.alleles.size(); ++index)
        if (record.alleles[index] != "<NON_REF>") proper_indices.push_back(index);
    // Java's stream.sorted comparator compares only the score.  Java's TimSort
    // is stable, so equal scores retain the input allele index; spell that
    // tie-break out to make the result independent of std::sort's choice.
    std::stable_sort(proper_indices.begin(), proper_indices.end(),
                     [&](const std::size_t left, const std::size_t right) {
                         return scores[left] > scores[right];
                     });
    std::vector<bool> keep(record.alleles.size(), false);
    keep[0] = true;
    int retained = 0;
    for (const auto index : proper_indices) {
        if (retained++ >= options.max_alternate_alleles) break;
        keep[index] = true;
    }
    for (std::size_t index = 1; index < record.alleles.size(); ++index)
        if (record.alleles[index] == "<NON_REF>") keep[index] = true;

    std::vector<std::string> selected;
    selected.reserve(static_cast<std::size_t>(options.max_alternate_alleles + 2));
    for (std::size_t index = 0; index < record.alleles.size(); ++index)
        if (keep[index]) selected.push_back(record.alleles[index]);
    if (selected.size() == record.alleles.size()) return false;
    const auto pruned = record.alleles.size() - selected.size();
    ++genotype_kernel_telemetry.max_alt_pruning_calls;
    genotype_kernel_telemetry.max_alt_alleles_pruned += pruned;
    remap_record_to_allele_union(output_header, record, selected);
    if (!record.gt.empty() && bcf_update_genotypes(
            output_header, record.value, record.gt.data(),
            static_cast<int>(record.gt.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish max-ALT GT subset");
    // GenotypingEngine's allele subset preserves a called genotype only when
    // every allele remains in the projection.  A genotype containing a
    // dropped ALT becomes a no-call and its scalar GQ is cleared; do not
    // re-derive it from an all-99 projected PL row (that would incorrectly
    // turn the no-call back into 0/0 with GQ=99).
    if (record.ploidy > 0 && record.gq.size() ==
            record.gt.size() / static_cast<std::size_t>(record.ploidy)) {
        for (std::size_t sample = 0; sample < record.gq.size(); ++sample) {
            bool missing = false;
            for (int copy = 0; copy < record.ploidy; ++copy) {
                const auto encoded = record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                                               static_cast<std::size_t>(copy)];
                if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                    missing = true;
                    break;
                }
            }
            if (missing) record.gq[sample] = bcf_int32_missing;
        }
    }
    // AlleleSubsettingUtils.scaleLogSpaceArrayForNumericalStability shifts
    // every retained PL row so its minimum is zero.  This is observable for
    // a sample whose best genotype used a dropped ALT: the projected row is
    // all 99 in the source scale but must be all zero in the new scale.
    const auto projected_width = genotype_width(record.allele_count, record.ploidy);
    if (projected_width != 0 && !record.pl.empty() &&
        record.pl.size() % projected_width == 0) {
        const auto projected_samples = record.pl.size() / projected_width;
        for (std::size_t sample = 0; sample < projected_samples; ++sample) {
            std::int32_t minimum = std::numeric_limits<std::int32_t>::max();
            for (std::size_t genotype = 0; genotype < projected_width; ++genotype) {
                const auto value = record.pl[sample * projected_width + genotype];
                if (value >= 0 && value < minimum) minimum = value;
            }
            if (minimum == std::numeric_limits<std::int32_t>::max() || minimum == 0) continue;
            for (std::size_t genotype = 0; genotype < projected_width; ++genotype) {
                auto& value = record.pl[sample * projected_width + genotype];
                if (value >= 0) value -= minimum;
            }
        }
    }
    const auto publish_format = [&](const char* tag, const std::vector<int32_t>& values) {
        if (values.empty()) return;
        if (bcf_update_format_int32(output_header, record.value, tag, values.data(),
                                    static_cast<int>(values.size())) != 0)
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot publish max-ALT ") + tag);
    };
    publish_format("DP", record.dp);
    publish_format("GQ", record.gq);
    publish_format("RGQ", record.rgq);
    publish_format("MIN_DP", record.min_dp);
    publish_format("AD", record.ad);
    publish_format("SB", record.sb);
    publish_format("PL", record.pl);
    publish_format("PP", record.pp);
    return true;
}

// NDA is intentionally captured before max-ALT reduction, and from the MERGED
// INPUT allele list rather than the published one.  That is exactly GATK's
// rule: GenotypingEngine.composeCallAttributes() stores
// `vc.getAlternateAlleles().size()` of the VariantContext handed to
// calculateGenotypes() (GenotypingEngine.java:464-465), i.e. the merged record
// -- before `*` ownership pruning (:312-316), before the confidence test and
// before the max-ALT reduction into `reducedVC` (:137-144).  Two details of
// that value are load-bearing:
//   * <NON_REF> is NOT counted, because the merger that built `vc` already
//     removed it (GenotypeGVCFsEngine.java:136 passes
//     removeNonRefSymbolicAllele = true;
//     ReferenceConfidenceVariantContextMerger.java:339-345);
//   * a symbolic '*' IS counted, because the merger re-adds Allele.SPAN_DEL
//     for a spanning event (ReferenceConfidenceVariantContextMerger.java:340-342).
// So the value is 1 for a single-concrete-ALT locus, and the old
// `record.alleles.size() < 3` guard (>= two ALTs) silently omitted NDA on every
// one-ALT shape -- including the REF-only dense row materialized from a covered
// '*', whose merged input set is [ref, *].  A locus whose merged ALT set is
// genuinely empty carries no NDA on either side: that denormalized REF-only
// record comes from regenotypeVC's non-variant branch -- the method entered at
// GenotypeGVCFsEngine.java:148 whose `originalVC.isVariant()` test at :154 fails
// -- and never reaches composeCallAttributes().
void annotate_num_discovered_alleles(const bcf_hdr_t* output_header,
                                     Record& record,
                                     const Options& options) {
    if (!options.annotate_with_num_discovered_alleles ||
        record.alleles.size() < 2 ||
        bcf_hdr_id2int(output_header, BCF_DT_ID, "NDA") < 0)
        return;
    int discovered = 0;
    for (std::size_t index = 1; index < record.alleles.size(); ++index)
        if (record.alleles[index] != "<NON_REF>") ++discovered;
    if (discovered <= 0) return;
    const auto value = static_cast<std::int32_t>(discovered);
    if (bcf_update_info_int32(output_header, record.value, "NDA", &value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write INFO/NDA");
    ++genotype_kernel_telemetry.num_discovered_alleles_calls;
}

// Remove the symbolic <NON_REF> allele after all shards have contributed their
// concrete ALT union.  A reference-only block has likelihoods for REF versus
// any-unseen-ALT; when materializing it at a concrete site, every concrete
// ALT inherits the REF/<NON_REF> likelihood and every ALT/ALT genotype
// inherits <NON_REF>/<NON_REF>.  Variant records keep their concrete PL rows
// and simply discard symbolic rows.  This is the deterministic single-sample
// projection used before the cohort-level allele-count kernel.
void remove_non_ref_allele(const bcf_hdr_t* output_header, Record& record,
                           const std::vector<std::string>& concrete_alleles,
                           bool best_match_to_original) {
    const auto non_ref_iter = std::find(record.alleles.begin(), record.alleles.end(), "<NON_REF>");
    if (non_ref_iter == record.alleles.end()) return;
    const int non_ref = static_cast<int>(non_ref_iter - record.alleles.begin());
    const int old_alleles = static_cast<int>(record.alleles.size());
    const int new_alleles = static_cast<int>(concrete_alleles.size());
    if (new_alleles < 2 || record.allele_count != old_alleles)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid <NON_REF> projection");
    const auto old_pl_width = genotype_width(old_alleles, record.ploidy);
    const auto new_pl_width = genotype_width(new_alleles, record.ploidy);
    if (old_pl_width == 0 || new_pl_width == 0)
        throw std::runtime_error("BAD_INPUT: unsupported genotype dimensions while removing <NON_REF>");

    std::vector<int32_t> target_to_old(static_cast<std::size_t>(new_alleles), -1);
    for (int target = 0; target < new_alleles; ++target) {
        if (record.reference_block) {
            target_to_old[static_cast<std::size_t>(target)] = target == 0 ? 0 : non_ref;
            continue;
        }
        const auto found = std::find(record.alleles.begin(), record.alleles.end(),
                                     concrete_alleles[static_cast<std::size_t>(target)]);
        if (found != record.alleles.end()) {
            const auto source_index = static_cast<int>(found - record.alleles.begin());
            if (source_index < non_ref)
            target_to_old[static_cast<std::size_t>(target)] = source_index;
        }
    }

    // MLEAC/MLEAF and any carried AF values use Number=A (ALT-only) indexing.
    // Compact them along with PL/AD before dropping the symbolic allele so a
    // caller opting out of the new AF calculator still gets GATK's projected
    // legacy annotations.
    remap_info_number_a_int(output_header, record.value, "MLEAC",
                            old_alleles, new_alleles, target_to_old);
    remap_info_number_a_float(output_header, record.value, "MLEAF",
                              old_alleles, new_alleles, target_to_old);
    remap_info_number_a_float(output_header, record.value, "AF",
                              old_alleles, new_alleles, target_to_old);

    if (!record.pl.empty()) {
        const auto projected = fastgatk::kernels::remap_genotype_pl_kokkos(
            record.pl, record.output_samples.size(), old_alleles, new_alleles,
            record.ploidy, target_to_old);
        account_pl_remap(projected);
        record.pl = projected.pl;
    }
    if (!record.pp.empty()) {
        const auto projected = fastgatk::kernels::remap_genotype_pl_kokkos(
            record.pp, record.output_samples.size(), old_alleles, new_alleles,
            record.ploidy, target_to_old);
        account_pl_remap(projected);
        record.pp = projected.pl;
    }

    if (!record.ad.empty()) {
        const auto projected = fastgatk::kernels::remap_allele_field_kokkos(
            record.ad, record.output_samples.size(), old_alleles, new_alleles,
            target_to_old);
        account_allele_field_remap(projected);
        record.ad = projected.values;
    }

    if (!record.gt.empty() && !record.reference_block) {
        for (auto& encoded : record.gt) {
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            const int allele = bcf_gt_allele(encoded);
            if (allele == non_ref && best_match_to_original)
                encoded = bcf_gt_is_phased(encoded) ? bcf_gt_phased(0) : bcf_gt_unphased(0);
            else if (allele < 0 || allele >= new_alleles) encoded = bcf_gt_missing;
        }
    }
    if (bcf_update_alleles_str(output_header, record.value,
                               allele_list_string(concrete_alleles).c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove <NON_REF> allele");
    // FORMAT arrays are kept in the staged Host representation here.  A
    // staged bcf1_t may still carry only its source-shard samples while the
    // merged output header already contains samples from later shards; the
    // final merge_sample_fields call writes the projected global arrays with
    // the correct output-header sample count.
    record.alleles = concrete_alleles;
    record.allele_count = new_alleles;
    record.reference_block = false;
}

// GATK MathUtils.log10OneMinusPow10 (via NaturalLogUtils.log1mexp): the log10
// of `1 - 10^log10_value`, computed without cancelling away the tiny
// complement.  GenotypingEngine reports this complementary posterior as the
// site QUAL of a locus that turned monomorphic
// (GenotypingEngine.java:158-163).
double gatk_log10_one_minus_pow10(const double log10_value) {
    if (log10_value > 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (log10_value == 0.0) return -std::numeric_limits<double>::infinity();
    constexpr double kLog10 = 2.30258509299404568402;
    constexpr double kLog2 = 0.69314718055994530942;
    const double natural = log10_value * kLog10;
    const double log1mexp = natural > -kLog2
        ? std::log(-std::expm1(natural)) : std::log1p(-std::exp(natural));
    return log1mexp / kLog10;
}

// GATK keeps a locus in --include-non-variant-sites mode even when every
// alternate allele failed the standard confidence threshold.  The rebuilt
// VariantContext is REF-only: GenotypingEngine.calculateGenotypes() builds its
// genotypes with GATKVariantContextUtils.subsetToRefOnly()
// (GenotypingEngine.java:190), then GenotypeGVCFsEngine.regenotypeVC()
// runs cleanupGenotypeAnnotations(result, createRefGTs=true, keepSB=false)
// (GenotypeGVCFsEngine.java:191-194, defined at :440).  That helper re-installs
// a hom-ref call only when the projected genotype still carries a positive
// depth *and* a GQ -- in which case the GQ becomes RGQ and DP is the MIN_DP
// (else the source DP) -- and otherwise publishes './.' with DP/GQ/PL removed
// (GenotypeGVCFsEngine.java:479-491).  The locus is written because it is not
// "spanning-deletion only": it has no ALT at all
// (GenotypeGVCFs.java:324-330, GATKVariantContextUtils.java:2089-2091).
//
// This must be materialized before the shared Number=G remap, which requires at
// least two target alleles and therefore cannot represent a REF-only output.
void materialize_gatk_monomorphic_ref_call(const bcf_hdr_t* output_header,
                                           Record& record) {
    if (record.value == nullptr || record.alleles.empty())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid REF-only call");
    const auto ploidy = record.ploidy > 0 ? record.ploidy : 2;
    const auto sample_count = static_cast<std::size_t>(record.value->n_sample);
    if (sample_count == 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: REF-only call has no samples");
    const auto reference = record.alleles.front();

    const auto present = [](const std::vector<int32_t>& field, const std::size_t index) {
        return field.size() > index && field[index] != bcf_int32_missing &&
               field[index] != bcf_int32_vector_end;
    };

    std::vector<int32_t> gt(sample_count * static_cast<std::size_t>(ploidy),
                            bcf_gt_missing);
    std::vector<int32_t> rgq(sample_count, bcf_int32_missing);
    std::vector<int32_t> dp(sample_count, bcf_int32_missing);
    bool any_rgq = false;
    // cleanupGenotypeAnnotations() keys the hom-ref conversion on the merged
    // genotype's GQ (GenotypeGVCFsEngine.java:485-488), so use the recorded
    // source-FORMAT presence rather than the native PL-derived GQ vector.
    int32_t* source_gq = nullptr;
    int source_gq_count = 0;
    if (record.source_has_gq)
        (void)bcf_get_format_int32(output_header, record.value, "GQ", &source_gq,
                                   &source_gq_count);
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        int32_t depth = 0;
        if (present(record.min_dp, sample)) depth = record.min_dp[sample];
        else if (present(record.dp, sample)) depth = record.dp[sample];
        const auto gq = source_gq != nullptr && static_cast<int>(sample) < source_gq_count
            ? source_gq[sample] : bcf_int32_missing;
        if (depth <= 0 || gq == bcf_int32_missing || gq == bcf_int32_vector_end)
            continue;
        if (gq > 0)
            for (int copy = 0; copy < ploidy; ++copy)
                gt[sample * static_cast<std::size_t>(ploidy) +
                   static_cast<std::size_t>(copy)] = bcf_gt_unphased(0);
        rgq[sample] = gq;
        dp[sample] = depth;
        any_rgq = true;
    }
    free(source_gq);

    // cleanupGenotypeAnnotations() drops AD/PL/PP/MIN_DP/SB and, for the
    // no-call shape, DP/GQ as well; a G0 call keeps DP and RGQ instead.
    const char* const dropped_format[] = {"AD", "PL", "PP", "MIN_DP", "SB",
                                          "GP", "PG", "GQ"};
    for (const auto* tag : dropped_format) {
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0) continue;
        if (std::strcmp(tag, "GP") == 0 || std::strcmp(tag, "PG") == 0)
            (void)bcf_update_format_float(output_header, record.value, tag, nullptr, 0);
        else
            (void)bcf_update_format_int32(output_header, record.value, tag, nullptr, 0);
    }
    for (const auto* tag : {"DP", "RGQ"}) {
        if (any_rgq || bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0) continue;
        (void)bcf_update_format_int32(output_header, record.value, tag, nullptr, 0);
    }
    if (bcf_update_genotypes(output_header, record.value, gt.data(),
                             static_cast<int>(gt.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write REF-only GT");
    if (any_rgq) {
        if (bcf_update_format_int32(output_header, record.value, "RGQ", rgq.data(),
                                    static_cast<int>(rgq.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write REF-only RGQ");
        if (bcf_update_format_int32(output_header, record.value, "DP", dp.data(),
                                    static_cast<int>(dp.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write REF-only DP");
    }
    // The Number=A MLEAC/MLEAF vectors survive as missing values on the empty
    // ALT list, which is how GenotypingEngine.composeCallAttributes()
    // (:437-441) renders a site with no emitted allele counts.
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MLEAC") >= 0) {
        const int32_t missing = bcf_int32_missing;
        if (bcf_update_info_int32(output_header, record.value, "MLEAC", &missing, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write REF-only MLEAC");
    }
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MLEAF") >= 0) {
        float missing = 0.0F;
        bcf_float_set_missing(missing);
        if (bcf_update_info_float(output_header, record.value, "MLEAF", &missing, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write REF-only MLEAF");
    }
    if (bcf_update_alleles_str(output_header, record.value, reference.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize REF-only alleles");

    // GenotypingEngine treats the locus as monomorphic here (no ALT passed the
    // standard confidence threshold), so the site confidence is the
    // complementary posterior rather than the AFCalculator's own QUAL
    // (GenotypingEngine.java:158-163); both formulas were already distinguished
    // for HaplotypeCaller's gVCF blocks.
    if (record.cohort_quality_available) {
        const auto complement = gatk_log10_one_minus_pow10(record.cohort_log10_p_no_variant);
        if (std::isfinite(complement)) {
            record.value->qual = gatk_qual_output(-10.0 * complement);
            // The confidence GenotypingEngine tests at :184-186 for a
            // monomorphic site is -10 * log10ProbVariantPresent(), i.e. the
            // same double, before the two-decimal rendering above.  A NaN
            // complement (log10Confidence > 0, unreachable for a real
            // posterior) keeps the previous leave-unchanged behaviour and is
            // NOT recorded, so that degenerate shape stays exactly where the
            // measured contract had it.
            record.call_confidence = -10.0 * complement;
            record.call_confidence_available = true;
        } else if (std::isinf(complement) && complement < 0.0) {
            // GenotypingEngine.java:158-163 assigns that same
            // log10ProbVariantPresent() straight into the record
            // (builder.log10PError(log10Confidence) at :183), so when the
            // posterior puts ALL of its mass on "no variant present",
            // MathUtils.log10OneMinusPow10(0.0) is Double.NEGATIVE_INFINITY
            // and GATK's QUAL really is Double.POSITIVE_INFINITY.  It is not a
            // sentinel: htsjdk renders it with ordinary decimal formatting, so
            // the published token is Java's `%.2f` spelling of an infinite
            // double, "Infinity".  gatk_qual_output() would flatten that to 0,
            // so keep the infinite value here and let
            // format_gatk_qual_value() supply the token.  A NaN complement
            // (log10Confidence > 0) keeps the previous leave-unchanged
            // behaviour.
            record.value->qual = std::numeric_limits<float>::infinity();
            // +Infinity passes every finite --standard-min-confidence-
            // threshold-for-calling, so the filter decision is recorded as
            // well: GATK's passesCallThreshold(-10 * -Infinity) is true.
            record.call_confidence = std::numeric_limits<double>::infinity();
            record.call_confidence_available = true;
        }
    }

    record.alleles = {reference};
    record.allele_count = 1;
    record.ploidy = ploidy;
    record.gt = std::move(gt);
    record.rgq = std::move(rgq);
    record.dp = std::move(dp);
    record.gq.clear();
    record.min_dp.clear();
    record.ad.clear();
    record.pl.clear();
    record.pp.clear();
    record.sb.clear();
    record.force_no_call_samples.clear();
    record.orphan_spanning_deletion = false;
    record.reference_block = false;
    record.finalized_monomorphic_ref = true;
}

// GATKVariantContextUtils.makeGenotypeCall()'s PREFER_PLS branch is a two-way
// split (GATKVariantContextUtils.java:331-352).  Only its first half derives
// the call from the likelihood row that AlleleSubsettingUtils.subsetAlleles()
// projected onto the KEPT alleles; when that row is not informative the row is
// ignored and the call is projected from the SOURCE genotype instead
// (bestMatchToOriginalGT: calls at :333-338, definition at :397-403).
//
// isInformative() is sum(log10GL) < SUM_GL_THRESH_NOCALL = -0.1 (:54-58), and
// the row the builder stores has already been shifted by its minimum
// (scaleLogSpaceArrayForNumericalStability, AlleleSubsettingUtils.java:90-95,
// followed by GenotypeLikelihoods.GLsToPLs).  On integer PLs that makes the
// test exactly sum(pl_i - min_pl) <= 1, which is shift-invariant and therefore
// also correct for a row native has not normalized yet.
bool gatk_prefer_pls_row_is_uninformative(const std::vector<int32_t>& pl,
                                          const std::size_t begin,
                                          const std::size_t end) {
    int64_t sum = 0;
    std::int32_t minimum = std::numeric_limits<std::int32_t>::max();
    std::size_t count = 0;
    for (std::size_t index = begin; index < end; ++index) {
        const auto value = pl[index];
        if (value < 0 || value == bcf_int32_vector_end) continue;
        sum += value;
        if (value < minimum) minimum = value;
        ++count;
    }
    // An entirely absent row is the `genotypeLikelihoods == null` arm of the
    // same Java condition.  Native keeps its existing no-call for that input
    // rather than projecting the source genotype.
    if (count == 0) return false;
    return sum - static_cast<int64_t>(minimum) * static_cast<int64_t>(count) <= 1;
}

// GenotypingEngine does not publish every ALT that arrived from the merged
// gVCF.  It asks AFCalculator whether each allele independently passes the
// standard confidence threshold, then subsets the VariantContext and all
// Number=A/Number=G fields before final annotation.  Retaining every source
// ALT is observably wrong for multi-ALT sites: an unsupported sibling changes
// the ALT list, PL width, AC/AF/MLEAC/MLEAF vectors and downstream QUAL/QD.
// Apply the same output-allele boundary after the shared Kokkos AF result and
// before the ordinary site/standard-annotation passes.
// GenotypeGVCFs' merged-depth precondition (GATK's regenotypeVC gate).
//
// GenotypeGVCFsEngine.java:157-174 enters the whole regenotyping block -- the
// output-allele subset, the reverse trim at :167 and every site annotation --
// only for a record that is a variant with a positive merged depth:
//
//     if ( originalVC.isVariant() && originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 ) {
//         ... calculateGenotypes / finalizeAnnotations / reverseTrimAlleles ...
//     } else {
//         result = originalVC;                       // :174 passthrough
//     }
//
// and :181 applies the same depth test again before anything is annotated or
// emitted, so that a passthrough record falls through to `return null` at :198
// unless --include-non-variant-sites is set.
//
// The depth is that of the MERGED record, and
// ReferenceConfidenceVariantContextMerger.calculateVCDepth() (:352-360) returns
// INFO/DP whenever the key is present -- never falling back to the genotypes --
// and otherwise sums getBestDepthValue() over the samples, which is MIN_DP when
// the genotype carries it and DP otherwise.
std::int64_t gatk_merged_record_depth(const bcf_hdr_t* output_header,
                                      const Record& record) {
    int32_t* info_dp = nullptr;
    int info_dp_count = 0;
    const auto info_length = bcf_get_info_int32(output_header, record.value, "DP",
                                                &info_dp, &info_dp_count);
    const bool has_info_dp = info_length > 0 && info_dp_count >= 1 &&
        info_dp[0] != bcf_int32_missing && info_dp[0] != bcf_int32_vector_end;
    std::int64_t depth = has_info_dp ? info_dp[0] : 0;
    free(info_dp);
    if (has_info_dp) return depth;
    for (std::size_t sample = 0; sample < record.dp.size(); ++sample) {
        const auto usable = [](const int32_t value) {
            return value != bcf_int32_missing && value != bcf_int32_vector_end &&
                   value >= 0;
        };
        const auto min_dp = sample < record.min_dp.size() ? record.min_dp[sample]
                                                          : bcf_int32_missing;
        const auto dp = record.dp[sample];
        const auto value = usable(min_dp) ? min_dp : (usable(dp) ? dp : bcf_int32_missing);
        if (usable(value)) depth += value;
    }
    return depth;
}

// True when GATK would not regenotype this locus at all, i.e. when native must
// publish nothing for it in the default (non dense) traversal.  The dense arm
// of the same rule is GATK's passthrough record (FORMAT `GT:AD` hom-ref no-call
// with INFO/DP=0) plus ALT='.' rows for the covered positions; that shape is
// entangled with the tracked covered-locus materialization gap and is left
// untouched here on purpose.
bool gatk_skips_regenotyping(const bcf_hdr_t* output_header, const Record& record,
                             const Options& options) {
    if (!options.gatk_annotation_compatibility) return false;
    if (options.include_non_variant_sites) return false;
    if (record.alleles.size() < 2) return false;  // originalVC.isVariant()
    return gatk_merged_record_depth(output_header, record) <= 0;
}

// GenotypeGVCFs' reverse allele trim (GATK's reverseTrimAlleles).
//
// GenotypeGVCFsEngine.java:160-169 regenotypes a locus, then rewrites the
// record before any site annotation runs:
//
//     finalizeAnnotations(...);                                     // :165
//     regenotypedVC = GATKVariantContextUtils.reverseTrimAlleles(regenotypedVC); // :167
//     annotateContext(regenotypedVC, ...);                          // :189
//
// GATKVariantContextUtils.reverseTrimAlleles() (:1443-1445) is
// trimAlleles(vc, trimForward=false, trimReverse=true), so only the trailing
// bases move:
//
//   * :1458 -- the guard.  A record with one allele, or with ANY allele whose
//     length() is 1 that is not the symbolic '*' (Allele.SPAN_DEL), is returned
//     untouched.  Because GenotypingEngine.calculateOutputAlleleSubset() has
//     already run (:155, mirrored by the caller of this helper), the guard sees
//     the SURVIVING allele list -- a one-base ALT that the standard-confidence
//     threshold pruned no longer protects the locus.
//   * :1462-1467 -- the candidates are the non-symbolic, non-'*' alleles (the
//     REF is one of them), and AlignmentUtils.normalizeAlleles(..., maxShift=0,
//     trim=true) (AlignmentUtils.java:818-845) consumes the common trailing run
//     while the shortest candidate still has a base left, so the run is capped
//     by the shortest candidate and endShift is that capped value.
//   * :1469-1475 -- if the cap was reached a whole allele was consumed, and
//     because nothing was clipped from the front (trimForward is false) exactly
//     one base is restored, so an allele can never be emptied.
//   * :1489-1521 -- every non-symbolic, non-'*' allele loses those trailing
//     bases (:1504) while symbolic alleles and '*' are copied untouched
//     (:1497-1501); genotypes are remapped by allele identity, so GT indices,
//     phase and the likelihood fields are unchanged; and the record is rebuilt
//     at the SAME start with stop = start + REF.length() - 1 (:1515-1518).
//
// Measured against pinned GATK 4.6.2.0 (see
// scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py): AAAA/AACA -> AAA/AAC
// at POS 2, ACGTACGT/ACGT -> ACGTA/A (the common run is four, three are
// clipped), AAAA/AACC unchanged, and AA/* -> A/* with the '*' intact.
//
// Placement: this helper must run *after* the emitted-deletion bookkeeping of
// EmittedDeletions::record(), because GATK records the deletion size of the
// alleles the locus emits at GenotypingEngine.java:178-179 -- before :167 --
// and a 'AA'->'A' rewrite would otherwise change that size.
void apply_gatk_reverse_trim(const bcf_hdr_t* output_header, Record& record) {
    if (record.alleles.size() <= 1) return;  // :1458 (getNAlleles() <= 1)
    const auto symbolic = [](const std::string& allele) {
        return !allele.empty() && allele.front() == '<';
    };
    for (const auto& allele : record.alleles) {
        // :1458.  '*' is Allele.SPAN_DEL and is explicitly exempt; a symbolic
        // allele is longer than one base, so it never triggers the guard.
        if (!symbolic(allele) && allele != "*" && allele.size() == 1) return;
    }
    std::vector<std::size_t> candidates;
    std::size_t shortest = std::numeric_limits<std::size_t>::max();
    for (std::size_t index = 0; index < record.alleles.size(); ++index) {
        const auto& allele = record.alleles[index];
        if (symbolic(allele) || allele == "*") continue;
        candidates.push_back(index);
        shortest = std::min(shortest, allele.size());
    }
    // Unreachable in practice (the REF is always a concrete allele, and a
    // one-base REF is removed by the guard above), but AlignmentUtils
    // .normalizeAlleles() requires a non-empty sequence list.
    if (candidates.empty() || shortest == 0) return;
    std::size_t end_trim = 0;
    while (end_trim < shortest && [&] {
               const auto& first = record.alleles[candidates.front()];
               const char probe = first[first.size() - 1 - end_trim];
               for (const auto index : candidates) {
                   const auto& allele = record.alleles[index];
                   if (allele[allele.size() - 1 - end_trim] != probe) return false;
               }
               return true;
           }())
        ++end_trim;
    // AlignmentUtils.java:818-838 stops when a candidate would be emptied, and
    // GATKVariantContextUtils.java:1469-1475 restores one base for exactly that
    // case (startTrim is 0 there, since the front loop only runs while
    // minSize > 0).
    const std::size_t rev_trim = end_trim == shortest ? end_trim - 1 : end_trim;
    if (rev_trim == 0) return;  // :1492-1493, nothing to do
    std::vector<std::string> trimmed;
    trimmed.reserve(record.alleles.size());
    for (const auto& allele : record.alleles) {
        if (symbolic(allele) || allele == "*") {
            trimmed.push_back(allele);
            continue;
        }
        trimmed.push_back(allele.substr(0, allele.size() - rev_trim));
    }
    std::string joined = trimmed.front();
    for (std::size_t index = 1; index < trimmed.size(); ++index) {
        joined += ',';
        joined += trimmed[index];
    }
    // htslib recomputes rlen from the new REF, which is GATK's
    // builder.stop(start + REF.length() - 1).
    if (bcf_update_alleles_str(output_header, record.value, joined.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot apply the reverse allele trim");
    record.alleles = std::move(trimmed);
}

// htsjdk rebuilds the subsetted genotypes through GenotypeLikelihoods, which
// normalises each sample's phred-scaled likelihood row so that its minimum is 0.
// The projection itself preserves the source row's offset, so a multi-allelic
// source whose best genotype does not survive the output-allele subset would
// otherwise be published shifted by that offset.  Measured on GATK's own chr20
// corpus: `20:10002458` GATK `2090,167,0` vs native `2172,249,82` (a uniform +82,
// the projected row's own minimum) and `20:10024300` +45.
void normalize_sample_pl(Record& record) {
    if (record.pl.empty() || record.ploidy <= 0 || record.allele_count < 2) return;
    const auto width = genotype_width(record.allele_count, record.ploidy);
    if (width == 0 || record.pl.size() % width != 0) return;
    const auto sample_count = record.pl.size() / width;
    const auto usable = [](const int32_t value) {
        return value != bcf_int32_missing && value != bcf_int32_vector_end && value >= 0;
    };
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto begin = record.pl.begin() + static_cast<std::ptrdiff_t>(sample * width);
        const auto end = begin + static_cast<std::ptrdiff_t>(width);
        int32_t minimum = std::numeric_limits<int32_t>::max();
        for (auto value = begin; value != end; ++value)
            if (usable(*value)) minimum = std::min(minimum, *value);
        if (minimum == std::numeric_limits<int32_t>::max() || minimum == 0) continue;
        for (auto value = begin; value != end; ++value)
            if (usable(*value)) *value -= minimum;
    }
}

bool apply_gatk_output_allele_subset(const bcf_hdr_t* output_header,
                                     Record& record,
                                     Options& options,
                                     EmittedDeletions& upstream_deletions) {
    // Keep the historical native diagnostic profile lossless for assignment
    // and kernel regression. GenotypeGVCFs' production-compatible profile is
    // the path where GATK's final AF-based allele pruning is observable.
    if (!options.gatk_annotation_compatibility) return true;
    // GenotypingEngine asks isVcCoveredByDeletion(vc) while it builds the output
    // allele subset (:314) and records what the locus emitted only afterwards
    // (:178-179), so what is consulted here is the state of the EARLIER loci.
    // Both paths of this tool decide what a locus emits in the compute stage of
    // their pipeline, which runs on a single worker thread in traversal order
    // (runtime/pipeline.hpp), so that is where GATK's per-traversal state is
    // maintained.  This is why the '*' ownership answer cannot be computed by
    // the union/merge step, which runs ahead of every emission decision.
    const bool carries_spanning_deletion =
        std::find(record.alleles.begin(), record.alleles.end(), "*") !=
        record.alleles.end();
    if (carries_spanning_deletion) {
        const bool covered = upstream_deletions.covers(record.rid, record.pos);
        record.orphan_spanning_deletion = !covered;
        if (!covered) ++options.orphan_spanning_deletion_loci;
    } else {
        record.orphan_spanning_deletion = false;
    }
    if (record.allele_count < 2 || record.alleles.size() < 2 ||
        record.cohort_log10_p_allele_absent.size() !=
            static_cast<std::size_t>(record.allele_count)) {
        // No subset was computed, so every allele of the merged record is
        // emitted; a reference-only record contributes nothing to the state.
        upstream_deletions.record(record, record.alleles);
        // GenotypeGVCFsEngine.java:167 runs for every locus it regenotypes,
        // whether or not the standard-confidence subset pruned anything.
        apply_gatk_reverse_trim(output_header, record);
        return true;
    }
    const auto threshold = -0.1 * options.standard_confidence_for_calling;
    std::vector<std::string> output_alleles;
    output_alleles.reserve(record.alleles.size());
    output_alleles.push_back(record.alleles.front());
    std::size_t pruned = 0;
    for (std::size_t index = 1; index < record.alleles.size(); ++index) {
        const auto& allele = record.alleles[index];
        const bool spanning_deletion = allele == "*";
        const auto absent = record.cohort_log10_p_allele_absent[index];
        const bool plausible = std::isfinite(absent) && absent + 1.0e-10 < threshold;
        // GenotypingEngine.calculateOutputAlleleSubset() outputs an ALT only when
        // it is individually plausible AND is not a spurious spanning deletion:
        //
        //   toOutput = (isPlausible || forceKeepAllele(allele)
        //               || isNonRefWhichIsLoneAltAllele || forcedAlleles.contains(allele))
        //              && !isSpuriousSpanningDeletion;
        //
        // (GenotypingEngine.java:316, with isPlausible = passesThreshold() at
        // :312 and isSpuriousSpanningDeletion = isSpanningDeletion(allele) &&
        // !isVcCoveredByDeletion(vc) at :314, answered by the emitted-deletion
        // state consulted at the top of this function).  An owned '*' is
        // therefore NOT exempt from the standard-confidence threshold: a '*'
        // that a concrete deletion covers but whose allele count does not pass
        // is pruned like any other ALT.  The flag only ever marks records
        // carrying a '*'; it must not gate a concrete ALT.
        const bool owned = !spanning_deletion || !record.orphan_spanning_deletion;
        if (owned && plausible)
            output_alleles.push_back(allele);
        else
            ++pruned;
    }
    // GenotypingEngine.java:172-175 -- a locus whose only survived ALT is the
    // spanning deletion is refused outright unless the traversal emits all
    // active sites, i.e. unless --include-non-variant-sites is set:
    //
    //   if (! emitAllActiveSites() && outputAlternativeAlleles.alleles.size() == 1
    //           && Allele.SPAN_DEL.equals(outputAlternativeAlleles.alleles.get(0)))
    //       return null;
    //
    // The refusal sits BEFORE recordDeletions() at :178-179, so a refused locus
    // must not contribute to the emitted-deletion state either, and before any
    // site annotation.  An ORPHAN '*' never reaches this point: it is removed
    // from the subset as a spurious spanning deletion (:314), leaving the
    // REF-only case handled below.
    if (!options.include_non_variant_sites && output_alleles.size() == 2 &&
        output_alleles[1] == "*") {
        bcf_destroy(record.value);
        record.value = nullptr;
        return false;
    }
    if (std::getenv("FASTGATK_DEBUG_SUBSET") != nullptr) {
        std::fprintf(stderr, "[SUBSET] pos=%d nalleles=%zu pruned=%zu output=[", record.pos,
                     record.alleles.size(), pruned);
        for (const auto& allele : output_alleles) std::fprintf(stderr, "%s,", allele.c_str());
        std::fprintf(stderr, "] absent=[");
        for (const auto value : record.cohort_log10_p_allele_absent)
            std::fprintf(stderr, "%.6g,", value);
        std::fprintf(stderr, "] orphan=%d includeNV=%d\n",
                     record.orphan_spanning_deletion ? 1 : 0,
                     options.include_non_variant_sites ? 1 : 0);
    }
    if (pruned == 0) {
        upstream_deletions.record(record, output_alleles);
        // Nothing was pruned, but htsjdk still rebuilds the genotypes for the
        // output allele list, which normalises each PL row to its minimum.
        normalize_sample_pl(record);
        // The emitted allele list is the merged one, and the trim still applies.
        apply_gatk_reverse_trim(output_header, record);
        return true;
    }
    ++genotype_kernel_telemetry.output_allele_pruning_calls;
    genotype_kernel_telemetry.output_alleles_pruned += pruned;
    if (output_alleles.size() == 1) {
        // No ALT survived, so the emitted allele list is the reference alone and
        // recordDeletions() has nothing to record at this locus
        // (GenotypingEngine.java:349-355 skips a zero deletion size).
        if (!options.include_non_variant_sites) {
            // GATK returns null at :167-169 before :179, so an unemitted locus
            // must not contribute to the state either.
            bcf_destroy(record.value);
            record.value = nullptr;
            return false;
        }
        // Dense mode keeps the locus as GATK's REF-only no-call instead of
        // projecting it through the Number=G remap, which has no representation
        // for a single-allele output.
        materialize_gatk_monomorphic_ref_call(output_header, record);
        return true;
    }
    // GenotypingEngine.java:178-179: the deletions of the alleles this locus
    // actually emits are recorded *after* the ownership test above and before
    // the genotypes are subsetted.
    upstream_deletions.record(record, output_alleles);
    // Preserve force-no-call states inherited from an earlier merge, but do
    // not let this final structural projection alone decide the new call.
    // GenotypingEngine re-evaluates PREFER_PLS after it has projected the
    // Number=G row; the projected PLs, rather than the disappearing source
    // '*', decide whether a sample remains callable.
    const auto force_before_output_subset = record.force_no_call_samples;
    remap_record_to_allele_union(output_header, record, output_alleles);
    if (record.orphan_spanning_deletion && !record.pl.empty() && record.ploidy > 0) {
        const auto width = genotype_width(record.allele_count, record.ploidy);
        if (width > 0 && record.pl.size() % width == 0) {
            const auto sample_count = record.pl.size() / width;
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                auto begin = record.pl.begin() + static_cast<std::ptrdiff_t>(sample * width);
                auto end = begin + static_cast<std::ptrdiff_t>(width);
                int32_t minimum = std::numeric_limits<int32_t>::max();
                for (auto value = begin; value != end; ++value)
                    if (*value >= 0 && *value != bcf_int32_vector_end)
                        minimum = std::min(minimum, *value);
                if (minimum != std::numeric_limits<int32_t>::max())
                    for (auto value = begin; value != end; ++value)
                        if (*value >= 0 && *value != bcf_int32_vector_end) *value -= minimum;
            }
        }
    }
    // GATK's final allele subset is a GenotypingEngine PREFER_PLS pass, not
    // a source-GT projection.  Rebuild GT/GQ only after the Kokkos Number=G
    // projection and orphan-span PL normalization are complete.  Suppress
    // the temporary no-call introduced by this remap so a confident p3
    // 0/1/1 result is allowed through; retain pre-existing no-call state.
    record.force_no_call_samples = force_before_output_subset;
    record.force_no_call_samples.resize(record.output_samples.size(), 0);
    // Which samples take PREFER_PLS's second arm.  The projected row is
    // measured after the Number=G remap and after the orphan-span
    // normalization above; the test itself is min-shift invariant.
    const auto subset_width = record.ploidy > 0
        ? genotype_width(record.allele_count, record.ploidy) : 0;
    std::vector<std::uint8_t> pls_row_ignored(record.output_samples.size(), 0);
    if (subset_width > 0 && !record.pl.empty() &&
        record.pl.size() % subset_width == 0 &&
        !record.source_gt.empty() && !record.source_alleles.empty()) {
        const auto rows = record.pl.size() / subset_width;
        for (std::size_t sample = 0;
             sample < pls_row_ignored.size() && sample < rows; ++sample) {
            if (gatk_prefer_pls_row_is_uninformative(
                    record.pl, sample * subset_width,
                    (sample + 1) * subset_width))
                pls_row_ignored[sample] = 1;
        }
    }
    derive_gt_gq_from_pl(record);
    // PREFER_PLS's second arm: the projected PL row carries no usable signal,
    // so GATK projects the SOURCE call allele by allele.  bestMatchToOriginalGT
    // keeps every source allele that survived the subset, replaces an allele
    // that was pruned by the reference, and keeps a no-call copy a no-call
    // (GATKVariantContextUtils.java:333-338, definition at :397-403); the
    // branch assigns no GQ, so a source without FORMAT/GQ leaves the field
    // absent (AlleleSubsettingUtils.java:127-128 also gates the projected GQ
    // on the source genotype's own hasGQ()).
    // The stored merged call indexes the union allele list, so it must have the
    // same sample-major shape as the GT vector it is projected onto.
    const bool source_call_usable = !record.source_gt.empty() &&
        !record.source_alleles.empty() &&
        record.source_gt.size() == record.gt.size();
    for (std::size_t sample = 0; sample < pls_row_ignored.size(); ++sample) {
        if (pls_row_ignored[sample] == 0) continue;
        if (sample < force_before_output_subset.size() &&
            force_before_output_subset[sample] != 0)
            continue;
        if (!source_call_usable) break;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto index = sample * static_cast<std::size_t>(record.ploidy) +
                               static_cast<std::size_t>(copy);
            if (index >= record.gt.size()) break;
            const auto source = record.source_gt[index];
            if (source == bcf_int32_vector_end) {
                record.gt[index] = bcf_int32_vector_end;
                continue;
            }
            if (bcf_gt_is_missing(source)) {
                record.gt[index] = bcf_gt_missing;
                continue;
            }
            // bestMatchToOriginalGT() compares Allele objects, so match by
            // allele name: a source allele that survived the subset keeps its
            // identity, and every other source allele -- a pruned ALT, the
            // symbolic '*' or '<NON_REF>' -- becomes the reference.
            const auto old_index = bcf_gt_allele(source);
            int new_index = 0;
            if (old_index >= 0 &&
                old_index < static_cast<int>(record.source_alleles.size())) {
                const auto found = std::find(
                    record.alleles.begin(), record.alleles.end(),
                    record.source_alleles[static_cast<std::size_t>(old_index)]);
                if (found != record.alleles.end())
                    new_index = static_cast<int>(found - record.alleles.begin());
            }
            record.gt[index] = bcf_gt_is_phased(source)
                ? bcf_gt_phased(new_index) : bcf_gt_unphased(new_index);
        }
        if (!record.source_has_gq && sample < record.gq.size())
            record.gq[sample] = bcf_int32_missing;
    }
    if (std::any_of(pls_row_ignored.begin(), pls_row_ignored.end(),
                    [](const std::uint8_t value) { return value != 0; }) &&
        !record.gq.empty() &&
        std::all_of(record.gq.begin(), record.gq.end(),
                    [](const int32_t value) {
                        return value == bcf_int32_missing ||
                               value == bcf_int32_vector_end;
                    })) {
        // No genotype carries a GQ, so htsjdk does not emit the FORMAT key at
        // all.  HTSlib would otherwise print an empty GQ column.
        record.gq.clear();
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, "GQ") >= 0)
            (void)bcf_update_format_int32(output_header, record.value, "GQ",
                                          nullptr, 0);
    }
    // GATKVariantContextUtils.makeGenotypeCall additionally turns an
    // effectively uninformative hom-ref PREFER_PLS result into a no-call
    // (SUM_GL_THRESH_NOCALL = -0.1).  With integer PL output that boundary
    // is exactly the all-reference, zero-GQ case.  This is distinct from a
    // confident hom-ref call, and preserves the GATK orphan-* behaviour
    // while allowing a confident projected non-reference genotype above.
    for (std::size_t sample = 0; sample < record.output_samples.size(); ++sample) {
        // The no-call rule lives in PREFER_PLS's first arm (:351-352), which
        // requires an informative row; a sample that took the second arm is
        // decided by the projected source call instead.
        if (sample < pls_row_ignored.size() && pls_row_ignored[sample] != 0)
            continue;
        if (sample < force_before_output_subset.size() &&
            force_before_output_subset[sample] != 0)
            continue;
        if (sample >= record.gq.size() || record.gq[sample] != 0) continue;
        bool hom_ref = true;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto encoded = record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                                           static_cast<std::size_t>(copy)];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded) ||
                bcf_gt_allele(encoded) != 0) {
                hom_ref = false;
                break;
            }
        }
        if (!hom_ref) continue;
        record.force_no_call_samples[sample] = 1;
        for (int copy = 0; copy < record.ploidy; ++copy)
            record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                      static_cast<std::size_t>(copy)] = bcf_gt_missing;
    }
    // MLEAF is reported after that final genotype pass.  An orphan '*'
    // briefly leaves a partial source GT during projection; calculating the
    // denominator before the PREFER_PLS replacement would incorrectly use
    // two called copies for a triploid sample instead of three.
    if (record.orphan_spanning_deletion && !record.gt.empty() &&
        bcf_hdr_id2int(output_header, BCF_DT_ID, "MLEAC") >= 0) {
        int32_t* mle_ac = nullptr;
        int mle_ac_count = 0;
        const auto length = bcf_get_info_int32(output_header, record.value, "MLEAC",
                                               &mle_ac, &mle_ac_count);
        if (length > 0 && mle_ac_count >= record.allele_count - 1) {
            int32_t called_alleles = 0;
            for (const auto encoded : record.gt)
                if (encoded != bcf_int32_vector_end && !bcf_gt_is_missing(encoded))
                    ++called_alleles;
            std::vector<float> mle_af(static_cast<std::size_t>(record.allele_count - 1), 0.0F);
            for (int allele = 1; allele < record.allele_count; ++allele)
                mle_af[static_cast<std::size_t>(allele - 1)] =
                    called_alleles == 0 ? 0.0F : static_cast<float>(mle_ac[allele - 1]) /
                    static_cast<float>(called_alleles);
            (void)bcf_update_info_float(output_header, record.value, "MLEAF",
                                        mle_af.data(), static_cast<int>(mle_af.size()));
        }
        free(mle_ac);
    }
    // `remap_record_to_allele_union` works on staged Host vectors because the
    // global sample-major arrays are not known at source-shard decode time.
    // This helper runs after shard merge, so publish the rebuilt vectors at
    // the writer boundary.
    if (!record.gt.empty() && bcf_update_genotypes(
            output_header, record.value, record.gt.data(),
            static_cast<int>(record.gt.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish output GT subset");
    const auto publish_format = [&](const char* tag, const std::vector<int32_t>& values) {
        if (values.empty()) return;
        if (bcf_update_format_int32(output_header, record.value, tag, values.data(),
                                    static_cast<int>(values.size())) != 0)
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot publish output ") + tag + " subset");
    };
    publish_format("DP", record.dp);
    publish_format("GQ", record.gq);
    publish_format("RGQ", record.rgq);
    publish_format("MIN_DP", record.min_dp);
    publish_format("AD", record.ad);
    publish_format("SB", record.sb);
    normalize_sample_pl(record);
    publish_format("PL", record.pl);
    publish_format("PP", record.pp);
    // The vectors are no longer indexed by the published ALT list.  They are
    // only diagnostics and must not accidentally be consumed by a later
    // annotation pass after remapping.
    record.cohort_log10_p_allele_absent.clear();
    record.cohort_integer_allele_counts.clear();
    // GenotypeGVCFsEngine.java:167: the trim runs on the projected record, after
    // the deletions this locus emitted were recorded at :178-179 above and
    // before any site annotation (which the callers run after this helper).
    apply_gatk_reverse_trim(output_header, record);
    return true;
}

// Materialize a pure REF/<NON_REF> block for GATK's
// --include-non-variant-sites mode.  Reference-only sites have exactly one
// legal genotype row (all reference), so the PL vector is compacted to [0]
// while the source block's DP/AD/GQ values are retained.  Host-side record
// rewriting keeps this path independent of HTSlib FORMAT ordering; the common
// merge_sample_fields routine still performs deterministic sample union.
void materialize_reference_only(const bcf_hdr_t* output_header, Record& record,
                                bool gatk_annotation_compatibility = false) {
    if (!record.reference_block || record.alleles.empty() || record.ploidy <= 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid reference-only record");
    const auto sample_count = record.output_samples.size();
    const auto old_alleles = record.alleles.size();
    if (record.allele_count != static_cast<int>(old_alleles))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: invalid reference-only allele count");

    // Keep a source-width REF/<NON_REF> PL copy for the cohort posterior.  A
    // reference-only output has a compact one-cell PL row, so doing this
    // after compaction would lose the probability mass of the unseen ALT
    // genotype.  Only the canonical two-allele reference block is eligible;
    // malformed/partial rows remain fail-closed and retain the legacy path.
    const auto source_pl_width = genotype_width(static_cast<int>(old_alleles), record.ploidy);
    if (old_alleles >= 2 && source_pl_width > 0 &&
        record.pl.size() == sample_count * source_pl_width) {
        record.reference_confidence_pl = record.pl;
        record.reference_confidence_allele_count = static_cast<int>(old_alleles);
    } else {
        record.reference_confidence_pl.clear();
        record.reference_confidence_allele_count = 0;
    }

    // GenotypeGVCFs emits a bare no-call for an uncovered coordinate in
    // --include-non-variant-sites mode: GT=./. and no DP/GQ/AD/PL/MIN_DP
    // values.  The native path previously converted the same zero-depth
    // block into 0/0:0:0:0, which is numerically harmless but changes the
    // reference-confidence contract.  Detect an all-zero staged depth row
    // before compacting fields and clear every genotype annotation while
    // retaining the representable missing GT.
    bool uncovered = !record.dp.empty() &&
        std::all_of(record.dp.begin(), record.dp.end(), [](const int32_t value) {
            return value == bcf_int32_missing || value == bcf_int32_vector_end || value <= 0;
        });
    if (gatk_annotation_compatibility) {
        // Java's cleanup computes depth after MIN_DP->DP recovery.  Treat a
        // missing depth vector as zero depth (rather than manufacturing a
        // 0/0 call), and use MIN_DP whenever the source block supplies it.
        uncovered = sample_count > 0;
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            int32_t value = bcf_int32_missing;
            if (record.min_dp.size() == sample_count)
                value = record.min_dp[sample];
            else if (record.dp.size() == sample_count)
                value = record.dp[sample];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end && value > 0) {
                uncovered = false;
                break;
            }
        }
    }
    if (uncovered) {
        record.gt.assign(sample_count * static_cast<std::size_t>(record.ploidy),
                         bcf_gt_missing);
        record.dp.clear();
        record.gq.clear();
        record.rgq.clear();
        record.min_dp.clear();
        record.ad.clear();
        record.pl.clear();
        record.pp.clear();
        record.sb.clear();
        record.reference_confidence_pl.clear();
        record.reference_confidence_allele_count = 0;
        const auto reference = record.alleles.front();
        if (bcf_update_alleles_str(output_header, record.value, reference.c_str()) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize uncovered REF-only site");
        (void)bcf_update_info_int32(output_header, record.value, "DP", nullptr, 0);
        (void)bcf_update_info_int32(output_header, record.value, "MIN_DP", nullptr, 0);
        record.alleles = {reference};
        record.allele_count = 1;
        record.reference_block = false;
        return;
    }
    const auto old_ad_width = static_cast<std::size_t>(record.allele_count);
    if (!record.ad.empty() && old_ad_width != 0) {
        std::vector<int32_t> ref_ad(sample_count, bcf_int32_missing);
        for (std::size_t sample = 0; sample < sample_count; ++sample)
            ref_ad[sample] = record.ad[sample * old_ad_width];
        record.ad = std::move(ref_ad);
    }
    record.pl.assign(sample_count, 0);
    if (!record.pp.empty()) record.pp.assign(sample_count, 0);
    record.gt.assign(sample_count * static_cast<std::size_t>(record.ploidy),
                     bcf_gt_unphased(0));
    if (record.gq.size() != sample_count)
        record.gq.assign(sample_count, bcf_int32_missing);
    if (record.min_dp.size() != sample_count)
        record.min_dp.assign(sample_count, bcf_int32_missing);
    // GATK emits a no-call for an otherwise reference-only site when the
    // reference-confidence quality is exactly zero.  Preserve that small
    // but observable contract while retaining 0/0 for finite RGQ/GQ.
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        if (record.gq[sample] != 0) continue;
        for (int copy = 0; copy < record.ploidy; ++copy)
            record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                      static_cast<std::size_t>(copy)] = bcf_gt_missing;
    }
    if (gatk_annotation_compatibility) {
        // GenotypeGVCFsEngine.cleanupGenotypeAnnotations(createRefGTs=true)
        // moves a finite GQ to FORMAT/RGQ and removes GQ, PL, PP and MIN_DP
        // for a monomorphic reference call.  DP remains the recovered depth.
        // When MIN_DP is present, that recovered depth replaces the original
        // block DP (the block's span-wide DP may be larger than the minimum
        // depth at this coordinate).
        if (record.min_dp.size() == sample_count) {
            if (record.dp.size() != sample_count)
                record.dp.assign(sample_count, bcf_int32_missing);
            for (std::size_t sample = 0; sample < sample_count; ++sample) {
                const auto value = record.min_dp[sample];
                if (value != bcf_int32_missing && value != bcf_int32_vector_end)
                    record.dp[sample] = value;
            }
        }
        record.rgq = record.gq;
        record.gq.clear();
        record.pl.clear();
        record.pp.clear();
        record.min_dp.clear();
        if (record.dp.size() != sample_count)
            record.dp.assign(sample_count, bcf_int32_missing);
        if (record.rgq.size() != sample_count)
            record.rgq.assign(sample_count, bcf_int32_missing);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            const auto depth = record.dp[sample];
            const auto quality = record.rgq[sample];
            if (depth == bcf_int32_missing || depth == bcf_int32_vector_end || depth <= 0 ||
                quality == bcf_int32_missing || quality == bcf_int32_vector_end) {
                record.dp[sample] = bcf_int32_missing;
                record.rgq[sample] = bcf_int32_missing;
                for (int copy = 0; copy < record.ploidy; ++copy)
                    record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                              static_cast<std::size_t>(copy)] = bcf_gt_missing;
            }
        }
    }
    const auto reference = record.alleles.front();
    if (bcf_update_alleles_str(output_header, record.value, reference.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize REF-only alleles");
    record.alleles = {reference};
    record.allele_count = 1;
    // FORMAT arrays are written only after all source shards have been merged.
    // Updating source-sized arrays here would violate HTSlib's invariant that
    // n_values == n_samples * width once another shard adds samples.
    record.reference_block = false;
}

// A reference-confidence block can span a concrete variant emitted by a
// different sample.  GATK splits that block at the variant locus so the
// hom-ref sample participates in joint genotyping.  Do the split before the
// normal locus/ALT-union pass: residual segments remain blocks, while a
// one-base point record at each concrete variant is projected through the
// existing REF/<NON_REF> logic and merged with the variant record.
void split_reference_blocks_at_variants(const bcf_hdr_t* output_header,
                                        faidx_t* reference_index,
                                        std::vector<Record>& records) {
    if (records.empty()) return;
    std::vector<VariantSpan> spans;
    for (const auto& record : records) {
        if (record.reference_block) continue;
        const auto end = record_span_end(output_header, record);
        if (record.rid >= 0 && record.pos >= 0 && end > record.pos)
            spans.push_back(VariantSpan{record.rid, record.pos, end});
    }
    if (spans.empty()) return;
    std::sort(spans.begin(), spans.end(), [](const VariantSpan& left,
                                             const VariantSpan& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });

    std::vector<Record> rewritten;
    rewritten.reserve(records.size() + spans.size());
    for (auto& original : records) {
        if (!original.reference_block) {
            rewritten.push_back(std::move(original));
            continue;
        }
        const int block_begin = original.pos;
        const int block_end = record_span_end(output_header, original);
        std::vector<VariantSpan> overlaps;
        for (const auto& span : spans) {
            if (span.rid != original.rid) continue;
            if (span.end <= block_begin) continue;
            if (span.begin >= block_end) break;
            overlaps.push_back(VariantSpan{
                span.rid, std::max(block_begin, span.begin),
                std::min(block_end, span.end)});
        }
        if (overlaps.empty()) {
            rewritten.push_back(std::move(original));
            continue;
        }
        // A one-base block already has the exact concrete-site key and needs
        // no reference lookup or cloning (this is the common single-site
        // reference sample case).
        if (block_end == block_begin + 1 && overlaps.size() == 1 &&
            overlaps.front().begin == block_begin) {
            rewritten.push_back(std::move(original));
            continue;
        }
        if (reference_index == nullptr)
            throw std::runtime_error(
                "UNSUPPORTED_PARAMETER: cross-sample reference-block splitting requires an indexed -R reference");

        // Coalesce overlapping concrete spans so one long indel does not
        // manufacture multiple point records or leave a residual block in
        // the covered reference span.
        std::vector<VariantSpan> merged;
        for (const auto& span : overlaps) {
            if (span.end <= span.begin) continue;
            if (merged.empty() || span.begin > merged.back().end) {
                merged.push_back(span);
            } else {
                merged.back().end = std::max(merged.back().end, span.end);
            }
        }
        std::vector<std::pair<int, int>> segments;
        int cursor = block_begin;
        for (const auto& span : merged) {
            if (span.begin > cursor)
                segments.emplace_back(cursor, span.begin);
            // The point record is what joins this sample's REF/<NON_REF>
            // likelihood to the concrete variant at span.begin.
            if (span.begin >= block_begin && span.begin < block_end)
                segments.emplace_back(span.begin, span.begin + 1);
            cursor = std::max(cursor, span.end);
        }
        if (cursor < block_end) segments.emplace_back(cursor, block_end);
        if (segments.empty()) {
            bcf_destroy(original.value);
            original.value = nullptr;
            continue;
        }

        for (const auto [segment_begin, segment_end] : segments) {
            auto* copy = bcf_dup(original.value);
            if (copy == nullptr)
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot split reference block");
            copy->pos = segment_begin;
            // A segment that begins at the block's own start is the record's own
            // locus, so it keeps the record's REF (see the per-coordinate
            // expansion below for the measurement).
            const auto reference = segment_begin == block_begin
                ? original.alleles.front()
                : reference_base(reference_index, output_header, original.rid,
                                 segment_begin);
            const auto alleles = reference + ",<NON_REF>";
            if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
                bcf_destroy(copy);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot rewrite split reference block alleles");
            }
            if (segment_end - segment_begin > 1) {
                const int32_t end_value = segment_end;
                if (bcf_update_info_int32(output_header, copy, "END", &end_value, 1) != 0) {
                    bcf_destroy(copy);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write split reference block END");
                }
            } else {
                (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
            }
            bcf_unpack(copy, BCF_UN_STR);
            Record split = original;
            split.value = copy;
            split.pos = segment_begin;
            split.alleles.front() = reference;
            split.key = record_key(copy);
            split.reference_block = true;
            rewritten.push_back(std::move(split));
        }
        bcf_destroy(original.value);
        original.value = nullptr;
    }
    records = std::move(rewritten);
}

// GATK's locus traversal visits EVERY reference coordinate that an input record
// spans, not only the coordinates where a record starts
// (VariantLocusWalker.java:150-176), and --include-non-variant-sites then
// publishes such a locus.  Native's locus list is built from record starts plus
// the per-coordinate expansion of pure reference blocks, so positions covered
// by a SPANNING record were missing entirely.
//
// At such a locus the merger replaces the spanning event's deletion allele with
// the symbolic '*' (ReferenceConfidenceVariantContextMerger.java:150-151,
// :222-245) and projects the sample data onto [reference base, '*']
// (mergeRefConfidenceGenotypes(), :575-612).  For a source whose surviving ALT
// is that single deletion the projection is the identity map on allele indices,
// which is why this pass can express it as an ordinary allele-union remap:
// keeping the source REF and the called deletion allele, then renaming the
// envelope to the new single-base REF and '*'.
//
// Only positions with no record of their own are synthesized (GATK's
// GenotypeGVCFsEngine.java:339-354 keeps the record that starts at the locus),
// and only when the source's own genotype calls the deletion: a hom-reference
// spanning record contributes NO_CALL for that allele, which is the shape GATK
// publishes as a REF-only row instead -- a different, still unfixed shape.
int called_deletion_allele_index(const Record& record) {
    if (record.alleles.size() < 2 || !record.ploidy) return -1;
    const auto reference_length = record.alleles.front().size();
    for (std::size_t index = 1; index < record.alleles.size(); ++index) {
        const auto& allele = record.alleles[index];
        if (allele == "*" || (!allele.empty() && allele.front() == '<')) continue;
        if (allele.size() >= reference_length) continue;
        for (const auto encoded : record.gt) {
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            if (bcf_gt_allele(encoded) == static_cast<int>(index))
                return static_cast<int>(index);
        }
    }
    return -1;
}

void materialize_spanning_loci(const bcf_hdr_t* output_header,
                               faidx_t* reference_index,
                               std::vector<Record>& records,
                               const std::vector<Region>& regions) {
    if (records.empty()) return;
    std::set<std::pair<int, int>> occupied;
    for (const auto& record : records)
        if (record.rid >= 0 && record.pos >= 0) occupied.emplace(record.rid, record.pos);
    struct Coverage {
        int rid = -1;
        int begin = -1;
        int end = -1;
        std::size_t source = 0;
        int deletion = -1;
    };
    std::vector<Coverage> coverage;
    for (std::size_t index = 0; index < records.size(); ++index) {
        const auto& record = records[index];
        if (record.reference_block || record.rid < 0 || record.pos < 0) continue;
        const auto end = record_span_end(output_header, record);
        if (end <= record.pos + 1) continue;
        const auto deletion = called_deletion_allele_index(record);
        if (deletion < 0) continue;
        // A source locus outside the requested intervals is never emitted, so
        // GATK never records its deletion and the covered coordinate would be an
        // ORPHAN '*' there.  Restricting the pass to in-interval sources keeps
        // the emitted-deletion state consistent with the traversal.
        if (!regions.empty() &&
            !starts_in_regions(output_header, record.value, regions))
            continue;
        coverage.push_back(Coverage{record.rid, record.pos + 1, end, index, deletion});
    }
    if (coverage.empty()) return;
    if (reference_index == nullptr)
        throw std::runtime_error(
            "UNSUPPORTED_PARAMETER: spanning-locus materialization requires an indexed -R reference");
    std::sort(coverage.begin(), coverage.end(), [](const Coverage& left, const Coverage& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });

    std::vector<Record> synthesized;
    for (const auto& span : coverage) {
        const auto& source = records[span.source];
        if (source.alleles.size() < 2 || span.deletion < 1 ||
            static_cast<std::size_t>(span.deletion) >= source.alleles.size())
            continue;
        for (int position = span.begin; position < span.end; ++position) {
            if (!occupied.emplace(span.rid, position).second) continue;
            auto* copy = bcf_dup(source.value);
            if (copy == nullptr)
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot materialize a spanning locus");
            Record synthetic = source;
            synthetic.value = copy;
            // Project the sample data onto [source REF, the called deletion],
            // which is the identity map on allele indices whenever the source's
            // only surviving ALT is that deletion.
            const std::vector<std::string> projected{
                source.alleles.front(), source.alleles[static_cast<std::size_t>(span.deletion)]};
            remap_record_to_allele_union(output_header, synthetic, projected);
            const auto base = reference_base(reference_index, output_header, span.rid, position);
            const auto alleles = base + ",*";
            if (bcf_update_alleles_str(output_header, synthetic.value, alleles.c_str()) != 0) {
                bcf_destroy(synthetic.value);
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: cannot rewrite a materialized spanning locus");
            }
            (void)bcf_update_info_int32(output_header, synthetic.value, "END", nullptr, 0);
            synthetic.value->pos = position;
            bcf_unpack(synthetic.value, BCF_UN_STR);
            synthetic.pos = position;
            synthetic.alleles = {base, "*"};
            synthetic.allele_count = 2;
            synthetic.key = record_key(synthetic.value);
            synthetic.reference_block = false;
            synthetic.materialized_spanning_locus = true;
            // The cohort vectors describe the SOURCE allele list; the normal
            // compute stage recomputes them for the two-allele projection.
            synthetic.cohort_log10_p_allele_absent.clear();
            synthetic.cohort_integer_allele_counts.clear();
            synthetic.cohort_log10_p_no_variant = 0.0;
            synthetic.cohort_quality_available = false;
            synthesized.push_back(std::move(synthetic));
        }
    }
    if (synthesized.empty()) return;
    records.reserve(records.size() + synthesized.size());
    for (auto& record : synthesized) records.push_back(std::move(record));
}

void merge_sample_fields(const bcf_hdr_t* output_header,
                         std::vector<Record*>& group) {
    if (group.empty()) return;
    const int output_samples = output_header->n[BCF_DT_SAMPLE];
    if (output_samples <= 0) return;
    // A merged locus carries GQ when any contributing source record did, which
    // is the property GATK's cleanupGenotypeAnnotations() tests with
    // g.hasGQ() on the ReferenceConfidenceVariantContextMerger output.
    for (std::size_t index = 1; index < group.size(); ++index)
        group.front()->source_has_gq = group.front()->source_has_gq || group[index]->source_has_gq;
    const auto copy_fields = [&](const char* tag, int width,
                                 auto Record::*member, int32_t missing) {
        if (std::strcmp(tag, "GT") != 0 &&
            bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0)
            return;
        bool present = false;
        for (const auto* record : group) if (!(record->*member).empty()) { present = true; break; }
        if (!present) {
            // The staged bcf1_t can still carry source-shard FORMAT data even
            // when the Host representation intentionally cleared it (for
            // example an uncovered REF-only site).  Clear the tag at the
            // merged output boundary instead of leaving stale source fields.
            if (std::strcmp(tag, "GT") == 0) {
                (void)bcf_update_genotypes(output_header, group.front()->value, nullptr, 0);
            } else if (bcf_hdr_id2int(output_header, BCF_DT_ID, tag) >= 0) {
                (void)bcf_update_format_int32(output_header, group.front()->value, tag, nullptr, 0);
            }
            return;
        }
        std::vector<int32_t> merged(static_cast<std::size_t>(output_samples) * width, missing);
        for (const auto* record : group) {
            const auto& source = record->*member;
            if (source.empty()) continue;
            const int source_samples = static_cast<int>(record->output_samples.size());
            for (int sample = 0; sample < source_samples; ++sample) {
                const auto output_sample = record->output_samples[static_cast<std::size_t>(sample)];
                if (output_sample < 0 || output_sample >= output_samples) continue;
                for (int index = 0; index < width; ++index)
                    merged[static_cast<std::size_t>(output_sample) * width + index] =
                        source[static_cast<std::size_t>(sample) * width + index];
            }
        }
        // Keep the staged Host representation in sync with the HTSlib record.
        // Site annotations are computed from these vectors after shard/sample
        // merge; leaving the first shard's vector here would silently omit
        // samples contributed by later disjoint headers.
        group.front()->*member = merged;
        if (std::strcmp(tag, "GT") == 0) {
            if (bcf_update_genotypes(output_header, group.front()->value, merged.data(),
                                     static_cast<int>(merged.size())) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge GT fields");
        } else if (bcf_update_format_int32(output_header, group.front()->value, tag,
                                            merged.data(), static_cast<int>(merged.size())) != 0) {
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot merge ") + tag + " fields");
        }
    };
    copy_fields("GT", group.front()->ploidy > 0 ? group.front()->ploidy : 2,
                &Record::gt, bcf_gt_missing);
    // The stored merged source call is a Host-only staging vector with the same
    // sample-major shape as GT; it has no FORMAT tag of its own.  Merge it with
    // the identical sample projection so PREFER_PLS's bestMatchToOriginalGT()
    // fallback sees every sample's own source call.  All members of a group
    // share the union allele list, so the first non-empty one names it.
    {
        const int source_width = group.front()->ploidy > 0
            ? group.front()->ploidy : 2;
        const bool present = std::any_of(group.begin(), group.end(),
            [](const Record* record) { return !record->source_gt.empty(); });
        if (present && output_samples > 0) {
            std::vector<int32_t> merged_source(
                static_cast<std::size_t>(output_samples) *
                    static_cast<std::size_t>(source_width), bcf_gt_missing);
            for (const auto* record : group) {
                if (record->source_gt.empty()) continue;
                const int source_samples = static_cast<int>(record->output_samples.size());
                for (int sample = 0; sample < source_samples; ++sample) {
                    const auto output_sample = record->output_samples[static_cast<std::size_t>(sample)];
                    if (output_sample < 0 || output_sample >= output_samples) continue;
                    for (int index = 0; index < source_width; ++index) {
                        const auto source_index =
                            static_cast<std::size_t>(sample) * source_width +
                            static_cast<std::size_t>(index);
                        if (source_index >= record->source_gt.size()) continue;
                        merged_source[static_cast<std::size_t>(output_sample) * source_width +
                                      static_cast<std::size_t>(index)] =
                            record->source_gt[source_index];
                    }
                }
            }
            group.front()->source_gt = std::move(merged_source);
        } else if (!present) {
            group.front()->source_gt.clear();
        }
        const auto named = std::find_if(group.begin(), group.end(),
            [](const Record* record) { return !record->source_alleles.empty(); });
        if (named != group.end()) group.front()->source_alleles = (*named)->source_alleles;
        else group.front()->source_alleles.clear();
    }
    copy_fields("DP", 1, &Record::dp, bcf_int32_missing);
    copy_fields("GQ", 1, &Record::gq, bcf_int32_missing);
    copy_fields("RGQ", 1, &Record::rgq, bcf_int32_missing);
    copy_fields("MIN_DP", 1, &Record::min_dp, bcf_int32_missing);
    const int allele_count = std::max(1, group.front()->allele_count);
    copy_fields("AD", allele_count, &Record::ad, bcf_int32_missing);
    copy_fields("SB", 4, &Record::sb, bcf_int32_missing);
    const auto pl_width = genotype_width(allele_count, group.front()->ploidy);
    if (pl_width == 0) throw std::runtime_error("BAD_INPUT: unsupported genotype PL dimensions during merge");
    copy_fields("PL", static_cast<int>(pl_width), &Record::pl, bcf_int32_missing);
    copy_fields("PP", static_cast<int>(pl_width), &Record::pp, bcf_int32_missing);
    copy_fields("PS", 1, &Record::ps, bcf_int32_missing);

    const auto copy_string_fields = [&](const char* tag,
                                        std::vector<std::string> Record::*member) {
        if (bcf_hdr_id2int(output_header, BCF_DT_ID, tag) < 0) return;
        bool present = false;
        for (const auto* record : group) {
            if (!(record->*member).empty()) {
                present = true;
                break;
            }
        }
        if (!present) {
            (void)bcf_update_format(output_header, group.front()->value, tag,
                                    nullptr, 0, BCF_HT_STR);
            return;
        }
        std::vector<std::string> merged(static_cast<std::size_t>(output_samples), ".");
        for (const auto* record : group) {
            const auto& source = record->*member;
            for (std::size_t sample = 0; sample < source.size() &&
                 sample < record->output_samples.size(); ++sample) {
                const auto output_sample = record->output_samples[sample];
                if (output_sample < 0 || output_sample >= output_samples) continue;
                merged[static_cast<std::size_t>(output_sample)] = source[sample].empty()
                    ? "." : source[sample];
            }
        }
        group.front()->*member = merged;
        std::vector<const char*> pointers;
        pointers.reserve(merged.size());
        for (const auto& value : merged) pointers.push_back(value.c_str());
        if (bcf_update_format_string(output_header, group.front()->value, tag,
                                     pointers.data(), static_cast<int>(pointers.size())) != 0)
            throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot merge ") +
                                     tag + " FORMAT field");
    };
    copy_string_fields("PGT", &Record::pgt);
    copy_string_fields("PID", &Record::pid);

    std::vector<std::uint8_t> merged_force_no_call(
        static_cast<std::size_t>(output_samples), 0);
    bool merged_orphan_spanning_deletion = false;
    for (const auto* record : group) {
        merged_orphan_spanning_deletion =
            merged_orphan_spanning_deletion || record->orphan_spanning_deletion;
        for (std::size_t sample = 0; sample < record->force_no_call_samples.size() &&
             sample < record->output_samples.size(); ++sample) {
            const auto output_sample = record->output_samples[sample];
            if (output_sample >= 0 && output_sample < output_samples &&
                record->force_no_call_samples[sample] != 0)
                merged_force_no_call[static_cast<std::size_t>(output_sample)] = 1;
        }
    }
    group.front()->force_no_call_samples = std::move(merged_force_no_call);
    group.front()->orphan_spanning_deletion = merged_orphan_spanning_deletion;

    // Reference-only records have already compacted their writer-facing PL to
    // width one.  Merge the preserved source-width REF/<NON_REF> rows in
    // parallel so the cross-sample reference posterior still sees every
    // contributing sample after disjoint shard union.
    // A group may contain a source shard with no reference-confidence PL row
    // (for example an uncovered sample) before another shard that does carry
    // the canonical REF/<NON_REF> row.  Select the first valid source width
    // instead of making merge order decide whether cohort diagnostics survive.
    int reference_allele_count = 0;
    for (const auto* record : group) {
        if (record->reference_confidence_allele_count >= 2 &&
            !record->reference_confidence_pl.empty()) {
            reference_allele_count = record->reference_confidence_allele_count;
            break;
        }
    }
    const auto reference_pl_width = reference_allele_count >= 2
        ? genotype_width(reference_allele_count, group.front()->ploidy) : 0;
    if (reference_pl_width > 0) {
        bool present = false;
        for (const auto* record : group)
            if (record->reference_confidence_allele_count == reference_allele_count &&
                record->reference_confidence_pl.size() ==
                    record->output_samples.size() * reference_pl_width) {
                present = true;
                break;
            }
        if (present) {
            std::vector<int32_t> merged_reference(
                static_cast<std::size_t>(output_samples) * reference_pl_width,
                bcf_int32_missing);
            for (const auto* record : group) {
                if (record->reference_confidence_allele_count != reference_allele_count ||
                    record->reference_confidence_pl.size() !=
                        record->output_samples.size() * reference_pl_width)
                    continue;
                const int source_samples = static_cast<int>(record->output_samples.size());
                for (int sample = 0; sample < source_samples; ++sample) {
                    const auto output_sample = record->output_samples[static_cast<std::size_t>(sample)];
                    if (output_sample < 0 || output_sample >= output_samples) continue;
                    for (std::size_t index = 0; index < reference_pl_width; ++index)
                        merged_reference[static_cast<std::size_t>(output_sample) * reference_pl_width + index] =
                            record->reference_confidence_pl[static_cast<std::size_t>(sample) * reference_pl_width + index];
                }
            }
            group.front()->reference_confidence_pl = std::move(merged_reference);
            group.front()->reference_confidence_allele_count = reference_allele_count;
        }
    }
}

void update_site_annotations(const bcf_hdr_t* output_header, Record& record) {
    if (record.allele_count < 2 || record.gt.empty() || record.output_samples.empty()) return;
    if (record.ploidy <= 0 || record.gt.size() % static_cast<std::size_t>(record.ploidy) != 0)
        return;
    // After disjoint-shard merging, output_samples on the first staged
    // record intentionally remains its source-shard mapping.  The merged GT
    // vector is the authoritative global sample-major layout for annotations.
    const auto sample_count = record.gt.size() / static_cast<std::size_t>(record.ploidy);
    if (sample_count == 0) return;
    std::vector<int32_t> allele_indices(record.gt.size(), -1);
    for (std::size_t index = 0; index < record.gt.size(); ++index) {
        const auto encoded = record.gt[index];
        if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
        const int allele = bcf_gt_allele(encoded);
        if (allele >= 0 && allele < record.allele_count) allele_indices[index] = allele;
    }
    const auto counts = fastgatk::kernels::count_alleles_kokkos(
        allele_indices, sample_count, record.ploidy, record.allele_count);
    ++genotype_kernel_telemetry.allele_count_calls;
    genotype_kernel_telemetry.allele_count_prepare_seconds += counts.prepare_seconds;
    genotype_kernel_telemetry.allele_count_execute_seconds += counts.seconds;
    genotype_kernel_telemetry.allele_count_execution_space = counts.execution_space;
    const int32_t an = counts.an;
    std::vector<int32_t> ac(static_cast<std::size_t>(record.allele_count - 1), 0);
    for (std::size_t index = 1; index < counts.counts.size(); ++index)
        ac[index - 1] = counts.counts[index];
    if (bcf_update_info_int32(output_header, record.value, "AC", ac.data(),
                              static_cast<int>(ac.size())) != 0 ||
        bcf_update_info_int32(output_header, record.value, "AN", &an, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AC/AN annotations");
    std::vector<float> af(ac.size(), 0.0F);
    for (std::size_t index = 0; index < ac.size(); ++index)
        af[index] = an == 0 ? 0.0F : static_cast<float>(ac[index]) / static_cast<float>(an);
    if (bcf_update_info_float(output_header, record.value, "AF", af.data(),
                              static_cast<int>(af.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write AF annotations");
}

// GATK's standard INFO/DP is the cohort depth, not the depth copied from the
// first source VariantContext.  This distinction is visible for an all-sites
// reference block after disjoint multi-sample FORMAT merge: each sample keeps
// its own FORMAT/DP, while the input block's INFO/DP may still contain only
// the first sample's value.  Recompute the total from the merged sample-major
// vector at the final annotation boundary, including REF-only records where
// update_site_annotations() intentionally has no ALT to count.
void update_total_depth_annotation(const bcf_hdr_t* output_header,
                                   Record& record,
                                   const Options& options) {
    if (!options.gatk_annotation_compatibility || record.dp.empty() ||
        bcf_hdr_id2int(output_header, BCF_DT_ID, "DP") < 0)
        return;
    // The cohort-depth correction is specific to dense reference-only
    // materialization.  For concrete variant records GATK's INFO/DP is
    // produced by the variant annotation engine (and can intentionally
    // differ from the sum of FORMAT/DP, e.g. after likelihood/read filters),
    // so rewriting it here regresses the established variant oracle.
    if (record.allele_count != 1)
        return;
    std::int64_t total = 0;
    bool any = false;
    for (const auto value : record.dp) {
        if (value == bcf_int32_missing || value == bcf_int32_vector_end || value < 0)
            continue;
        any = true;
        total += value;
    }
    if (!any) {
        (void)bcf_update_info_int32(output_header, record.value, "DP", nullptr, 0);
        return;
    }
    const auto bounded = std::clamp<std::int64_t>(
        total, 0, std::numeric_limits<std::int32_t>::max());
    const auto depth = static_cast<std::int32_t>(bounded);
    if (bcf_update_info_int32(output_header, record.value, "DP", &depth, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write cohort INFO/DP");
}

// GATK's ExcessHet annotation is the right-sided exact Hardy-Weinberg test
// described by Wigginton, Cutler, and Abecasis.  It is intentionally kept as
// a small Host-side per-record calculation: the high-volume PL/AF work stays
// in the shared Kokkos kernels, while this irregular recurrence is bounded by
// the number of samples.  GATK only defines this annotation for diploid
// samples; fail closed for polyploid and multiallelic records rather than
// emitting a value with different semantics.
double exact_excess_het_phred(int hom_ref, int het, int hom_alt) {
    if (hom_ref < 0 || het < 0 || hom_alt < 0) return std::numeric_limits<double>::quiet_NaN();
    const int samples = hom_ref + het + hom_alt;
    if (samples == 0) return std::numeric_limits<double>::quiet_NaN();
    const int rare_copies = 2 * std::min(hom_ref, hom_alt) + het;
    int midpoint = rare_copies * (2 * samples - rare_copies) / (2 * samples);
    if ((midpoint & 1) != (rare_copies & 1)) ++midpoint;
    std::vector<double> probabilities(static_cast<std::size_t>(rare_copies) + 1, 0.0);
    probabilities[static_cast<std::size_t>(midpoint)] = 1.0;
    double probability_sum = 1.0;

    int current_hom_ref = (rare_copies - midpoint) / 2;
    int current_hom_alt = samples - current_hom_ref - midpoint;
    for (int current_het = midpoint; current_het > 1; current_het -= 2) {
        const auto value = probabilities[static_cast<std::size_t>(current_het)] *
            static_cast<double>(current_het) * static_cast<double>(current_het - 1) /
            (4.0 * static_cast<double>(current_hom_ref + 1) *
             static_cast<double>(current_hom_alt + 1));
        probabilities[static_cast<std::size_t>(current_het - 2)] = value;
        probability_sum += value;
        ++current_hom_ref;
        ++current_hom_alt;
    }

    current_hom_ref = (rare_copies - midpoint) / 2;
    current_hom_alt = samples - current_hom_ref - midpoint;
    for (int current_het = midpoint; current_het <= rare_copies - 2; current_het += 2) {
        const auto value = probabilities[static_cast<std::size_t>(current_het)] *
            4.0 * static_cast<double>(current_hom_ref) * static_cast<double>(current_hom_alt) /
            (static_cast<double>(current_het + 2) * static_cast<double>(current_het + 1));
        probabilities[static_cast<std::size_t>(current_het + 2)] = value;
        probability_sum += value;
        --current_hom_ref;
        --current_hom_alt;
    }
    if (!(probability_sum > 0.0) || !std::isfinite(probability_sum) ||
        het < 0 || het > rare_copies)
        return std::numeric_limits<double>::quiet_NaN();

    // Right tail: probability of observing the called number of heterozygotes
    // or more.  The small tolerance handles equal-probability states produced
    // by different recurrence paths without changing the exact tail.
    double tail = 0.0;
    for (int current_het = het; current_het <= rare_copies; current_het += 2)
        tail += probabilities[static_cast<std::size_t>(current_het)];
    // The recurrence visits only parity-compatible heterozygosity states; the
    // step of two above keeps the exact-test support explicit.
    if (!(tail > 0.0) || !std::isfinite(tail)) return std::numeric_limits<double>::quiet_NaN();
    const auto p_value = std::clamp(tail / probability_sum, 0.0, 1.0);
    if (!(p_value > 0.0)) return 999.0;
    return std::clamp(-10.0 * std::log10(p_value), 0.0, 999.0);
}

void update_excess_het_annotation(const bcf_hdr_t* output_header, Record& record) {
    // GATK's ExcessHet is a diploid cohort annotation.  Its genotype-count
    // reducer is not restricted to biallelic sites: 0/ALT genotypes enter the
    // heterozygous bucket, while ALT/ALT (including two different ALTs) enter
    // the homozygous-variant/non-reference bucket.  Restricting this helper
    // to allele_count==2 silently drops ExcessHet on complex multiallelic
    // records even though GenotypeGVCFs emits a valid 0.0000 value.
    if (record.allele_count < 2 || record.ploidy != 2 || record.gt.empty() ||
        record.gt.size() % 2 != 0)
        return;
    ++genotype_kernel_telemetry.excess_het_calls;
    int hom_ref = 0;
    int het = 0;
    int hom_alt = 0;
    const auto sample_count = record.gt.size() / 2;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto first = record.gt[sample * 2];
        const auto second = record.gt[sample * 2 + 1];
        if (first == bcf_int32_vector_end || second == bcf_int32_vector_end ||
            bcf_gt_is_missing(first) || bcf_gt_is_missing(second))
            continue;
        const int first_allele = bcf_gt_allele(first);
        const int second_allele = bcf_gt_allele(second);
        if (first_allele < 0 || second_allele < 0 ||
            first_allele >= record.allele_count || second_allele >= record.allele_count)
            return;
        if (first_allele == 0 && second_allele == 0) ++hom_ref;
        else if (first_allele == 0 || second_allele == 0) ++het;
        else ++hom_alt;
    }
    const auto phred = exact_excess_het_phred(hom_ref, het, hom_alt);
    if (!std::isfinite(phred)) return;
    const float value = static_cast<float>(phred);
    if (bcf_update_info_float(output_header, record.value, "ExcessHet", &value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write ExcessHet annotation");
}

// GATK's InbreedingCoeff is derived from the same diploid genotype-count
// reducer used by GenotypeUtils, but unlike ExcessHet it keeps the normalized
// likelihood mass (ROUND_GENOTYPE_COUNTS=false).  This is a deliberately
// bounded per-record reducer: the high-volume PL/AF work is still evaluated
// by the shared Kokkos cohort kernel, while this annotation preserves GATK's
// exact multiallelic "best alternate" projection and GQ-only hom-ref path.
void update_inbreeding_coeff_annotation(const bcf_hdr_t* output_header,
                                        Record& record) {
    constexpr std::size_t minimum_samples = 10;
    if (record.allele_count < 2 || record.ploidy != 2 || record.gt.empty() ||
        record.gt.size() % 2 != 0)
        return;

    const std::size_t sample_count = record.gt.size() / 2;
    if (sample_count < minimum_samples) return;

    const auto width = genotype_width(record.allele_count, record.ploidy);
    if (width < 3) return;
    if (!record.pl.empty() && record.pl.size() != sample_count * width) return;
    if (!record.gq.empty() && record.gq.size() != sample_count) return;

    double refs = 0.0;
    double hets = 0.0;
    double homs = 0.0;
    std::size_t usable_samples = 0;

    const auto sample_is_hom_ref = [&](const std::size_t sample) {
        const auto first = record.gt[sample * 2];
        const auto second = record.gt[sample * 2 + 1];
        return first != bcf_int32_vector_end && second != bcf_int32_vector_end &&
               !bcf_gt_is_missing(first) && !bcf_gt_is_missing(second) &&
               bcf_gt_allele(first) == 0 && bcf_gt_allele(second) == 0;
    };

    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const std::size_t offset = sample * width;
        bool has_likelihoods = false;
        if (!record.pl.empty()) {
            for (std::size_t genotype = 0; genotype < width; ++genotype) {
                const auto value = record.pl[offset + genotype];
                if (value >= 0 && value != bcf_int32_vector_end) {
                    has_likelihoods = true;
                    break;
                }
            }
        }

        if (!has_likelihoods) {
            // GenotypeUtils only synthesizes a contribution for a called
            // diploid hom-ref with GQ when PL is absent.  A GQ of zero is
            // intentionally split into three equal genotype masses.
            if (record.gq.empty() || record.gq[sample] < 0 ||
                !sample_is_hom_ref(sample))
                continue;
            ++usable_samples;
            const auto gq = record.gq[sample];
            if (gq == 0) {
                refs += 1.0 / 3.0;
                hets += 1.0 / 3.0;
                homs += 1.0 / 3.0;
            } else {
                const auto ref_probability =
                    1.0 - std::pow(10.0, -static_cast<double>(gq) / 10.0);
                refs += ref_probability;
                hets += 1.0 - ref_probability;
            }
            continue;
        }

        // Normalize PLs in linear space after subtracting the best PL.  This
        // is equivalent to MathUtils.normalizeFromLog10ToLinearSpace while
        // remaining stable for the large PL values emitted by gVCFs.
        int minimum_pl = std::numeric_limits<int>::max();
        for (std::size_t genotype = 0; genotype < width; ++genotype) {
            const auto value = record.pl[offset + genotype];
            if (value >= 0 && value < minimum_pl) minimum_pl = value;
        }
        if (minimum_pl == std::numeric_limits<int>::max()) continue;
        std::vector<double> normalized(width, 0.0);
        double normalizer = 0.0;
        for (std::size_t genotype = 0; genotype < width; ++genotype) {
            const auto value = record.pl[offset + genotype];
            if (value < 0 || value == bcf_int32_vector_end) continue;
            const auto likelihood = std::pow(
                10.0, -0.1 * static_cast<double>(value - minimum_pl));
            normalized[genotype] = std::isfinite(likelihood) ? likelihood : 0.0;
            normalizer += normalized[genotype];
        }
        if (!(normalizer > 0.0) || !std::isfinite(normalizer)) continue;
        for (auto& value : normalized) value /= normalizer;
        ++usable_samples;

        if (record.allele_count == 2) {
            refs += normalized[0];
            hets += normalized[1];
            homs += normalized[2];
            continue;
        }

        // For multiallelic sites GenotypeUtils first checks the maximum
        // likelihood genotype.  A maximum with no REF allele contributes all
        // mass to the homozygous-variant bucket; otherwise it projects the
        // full likelihood vector to the best biallelic REF/ALT pair.
        std::size_t max_index = 0;
        for (std::size_t genotype = 1; genotype < width; ++genotype)
            if (normalized[genotype] > normalized[max_index]) max_index = genotype;
        std::vector<int> max_alleles;
        unrank_genotype(max_index, record.allele_count, record.ploidy, max_alleles);
        if (max_alleles.size() != 2 ||
            (max_alleles[0] != 0 && max_alleles[1] != 0)) {
            homs += 1.0;
            continue;
        }

        std::size_t het_index = rank_genotype({0, 1}, record.allele_count);
        std::size_t var_index = rank_genotype({1, 1}, record.allele_count);
        double max_likelihood = normalized[het_index];
        for (int allele = 1; allele < record.allele_count; ++allele) {
            const auto candidate_het = rank_genotype({0, allele}, record.allele_count);
            const auto candidate_var = rank_genotype({allele, allele}, record.allele_count);
            if (candidate_het >= normalized.size() || candidate_var >= normalized.size())
                continue;
            if (normalized[candidate_het] > max_likelihood) {
                max_likelihood = normalized[candidate_het];
                het_index = candidate_het;
                var_index = candidate_var;
            }
        }
        const auto pair_sum = normalized[0] + normalized[het_index] + normalized[var_index];
        if (!(pair_sum > 0.0) || !std::isfinite(pair_sum)) continue;
        refs += normalized[0] / pair_sum;
        hets += normalized[het_index] / pair_sum;
        homs += normalized[var_index] / pair_sum;
    }

    if (usable_samples < minimum_samples) return;
    const auto denominator = 2.0 * (refs + hets + homs);
    if (!(denominator > 0.0) || !std::isfinite(denominator)) return;
    const auto p = (2.0 * refs + hets) / denominator;
    const auto q = 1.0 - p;
    const auto expected_hets = 2.0 * p * q * static_cast<double>(usable_samples);
    if (!(expected_hets > 0.0) || !std::isfinite(expected_hets)) return;
    const auto coefficient = 1.0 - hets / expected_hets;
    if (!std::isfinite(coefficient)) return;
    const float value = static_cast<float>(coefficient);
    if (bcf_update_info_float(output_header, record.value, "InbreedingCoeff", &value, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write InbreedingCoeff annotation");
    ++genotype_kernel_telemetry.inbreeding_coeff_calls;
}

// GATK's FisherStrand annotation uses a two-sided hypergeometric test and
// emits a phred-scaled p-value.  Keeping the calculation in Host code is
// intentional: this is a small per-record annotation, while all high-volume
// likelihood/count kernels remain Kokkos based.  The log-combination form
// avoids integer overflow and matches the Java FisherExactTest for the small
// tables normally carried by FORMAT/SB.
double fisher_two_sided_pvalue(int a, int b, int c, int d) {
    a = std::max(0, a); b = std::max(0, b);
    c = std::max(0, c); d = std::max(0, d);
    const int row0 = a + b;
    const int row1 = c + d;
    const int col0 = a + c;
    const int col1 = b + d;
    const int total = row0 + row1;
    if (total <= 0) return 1.0;
    const auto log_probability = [&](int x) {
        const int y = row0 - x;
        const int z = col0 - x;
        const int w = row1 - z;
        if (x < 0 || y < 0 || z < 0 || w < 0) return -std::numeric_limits<double>::infinity();
        return std::lgamma(static_cast<double>(row0) + 1.0) +
               std::lgamma(static_cast<double>(row1) + 1.0) +
               std::lgamma(static_cast<double>(col0) + 1.0) +
               std::lgamma(static_cast<double>(col1) + 1.0) -
               std::lgamma(static_cast<double>(total) + 1.0) -
               std::lgamma(static_cast<double>(x) + 1.0) -
               std::lgamma(static_cast<double>(y) + 1.0) -
               std::lgamma(static_cast<double>(z) + 1.0) -
               std::lgamma(static_cast<double>(w) + 1.0);
    };
    const auto observed = log_probability(a);
    const int lower = std::max(0, row0 - col1);
    const int upper = std::min(row0, col0);
    double probability = 0.0;
    for (int x = lower; x <= upper; ++x) {
        const auto value = log_probability(x);
        // The Java implementation includes tables no more likely than the
        // observed table.  A small tolerance stabilizes the boundary when
        // lgamma evaluates two mathematically equal integer tables.
        if (value <= observed + 1.0e-12)
            probability += std::exp(value);
    }
    return std::clamp(probability, 0.0, 1.0);
}

std::array<int, 4> merged_strand_bias_table(const bcf_hdr_t* output_header,
                                            const Record& record) {
    std::array<int, 4> table{0, 0, 0, 0};
    std::vector<int32_t> values = record.sb;
    if (values.empty()) {
        int32_t* raw = nullptr;
        int raw_count = 0;
        const auto length = bcf_get_format_int32(output_header, record.value, "SB",
                                                  &raw, &raw_count);
        if (length > 0 && raw_count > 0)
            values.assign(raw, raw + raw_count);
        free(raw);
    }
    if (values.size() < 4 || values.size() % 4 != 0) return table;
    for (std::size_t offset = 0; offset + 4 <= values.size(); offset += 4) {
        bool valid = true;
        for (int index = 0; index < 4; ++index) {
            const auto value = values[offset + static_cast<std::size_t>(index)];
            if (value == bcf_int32_missing || value == bcf_int32_vector_end || value < 0) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;
        // GATK only contributes a sample when the combined table has more
        // than ARRAY_DIM observations (MIN_COUNT=2).
        for (int index = 0; index < 4; ++index)
            table[static_cast<std::size_t>(index)] +=
                values[offset + static_cast<std::size_t>(index)];
    }
    return table;
}

void update_gatk_standard_annotations(const bcf_hdr_t* output_header, Record& record,
                                      GatkJavaRandom& random) {
    if (record.allele_count < 2 || record.gt.empty() || record.ploidy <= 0)
        return;

    // RMSMappingQuality finalizes the raw [sum(mapq^2), depth] pair produced
    // by HaplotypeCaller.  Keep RAW_MQandDP as an input contract but expose
    // the standard MQ annotation as GATK does.
    int32_t* raw_mq = nullptr;
    int raw_mq_count = 0;
    const auto raw_mq_length = bcf_get_info_int32(output_header, record.value,
                                                   "RAW_MQandDP", &raw_mq,
                                                   &raw_mq_count);
    if (raw_mq_length > 0 && raw_mq_count >= 2 && raw_mq[1] > 0) {
        const float mq = static_cast<float>(std::sqrt(
            static_cast<double>(raw_mq[0]) / static_cast<double>(raw_mq[1])));
        if (bcf_update_info_float(output_header, record.value, "MQ", &mq, 1) != 0) {
            free(raw_mq);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write MQ annotation");
        }
    }
    free(raw_mq);

    // QualByDepth uses AD depth from samples carrying a non-reference call.
    // If at least one sample has more than one ALT-supporting observation,
    // GATK switches to the restricted depth sum; otherwise the ordinary AD
    // depth is retained.  This reproduces the relevant GenotypeGVCFs path and
    // avoids using the all-sample INFO/DP for a multi-sample variant.
    const auto sample_count = record.gt.size() /
                               static_cast<std::size_t>(record.ploidy);
    const auto ad_width = static_cast<std::size_t>(record.allele_count);
    int depth = 0;
    int restricted_depth = 0;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        bool non_reference = false;
        for (int copy = 0; copy < record.ploidy; ++copy) {
            const auto encoded = record.gt[sample * static_cast<std::size_t>(record.ploidy) +
                                           static_cast<std::size_t>(copy)];
            if (encoded != bcf_int32_vector_end && !bcf_gt_is_missing(encoded) &&
                bcf_gt_allele(encoded) > 0) {
                non_reference = true;
                break;
            }
        }
        if (!non_reference) continue;
        int sample_depth = 0;
        int alt_depth = 0;
        if (record.ad.size() >= (sample + 1) * ad_width) {
            for (std::size_t allele = 0; allele < ad_width; ++allele) {
                const auto value = record.ad[sample * ad_width + allele];
                if (value == bcf_int32_missing || value < 0) continue;
                sample_depth += value;
                if (allele > 0) alt_depth += value;
            }
        }
        if (sample_depth > 0) {
            depth += sample_depth;
            if (alt_depth > 1) restricted_depth += sample_depth;
        } else if (sample < record.dp.size() && record.dp[sample] > 0) {
            depth += record.dp[sample];
        }
    }
    if (restricted_depth > 0) depth = restricted_depth;
    if (depth > 0 && std::isfinite(record.value->qual)) {
        // GATK's QualByDepth divides the DOUBLE it reads back from the record
        // (`-10.0 * vc.getLog10PError()`, QualByDepth.java:78+86), while this
        // tool's published QUAL is that double rounded to two decimals.  The two
        // differ only in the SIGN OF ZERO for a zero-confidence locus, and GATK
        // itself is not consistent there: measured in one and the same fixture,
        // the '*'-only locus at 3 comes out QD=-0.00 and the one at 4 comes out
        // QD=0.00 (the sign survives from a ~1e-16 round-off in the AF
        // calculator, whose direction is not reproducible without matching
        // GATK's arithmetic order).  A presentation-layer sign fix was tried and
        // reverted: it made the 3 case match and the 4 case diverge, which the
        // registered reverse-trim gate caught.  Keeping the rounded numerator is
        // therefore the measured-better choice until the AF kernel's near-zero
        // result can be aligned.
        // GATK's QualByDepth divides the double it reads back from the record
        // (-10.0 * vc.getLog10PError(), QualByDepth.java:78+86), and only the
        // sign of ZERO can differ from the rounded QUAL published here.  For a
        // REGENOTYPED locus that sign survives from a ~1e-16 round-off in the
        // AF calculator and is not reproducible without matching GATK's
        // arithmetic order (measured: QD=-0.00 for one '*' locus and QD=0.00 for
        // another inside a single fixture), so the rounded numerator is kept.
        // A locus this tool MATERIALIZED because a spanning record covers it is
        // different: it has no input record of its own, its confidence is
        // exactly zero, and pinned GATK publishes QD=-0.00 for every measured
        // instance of that shape (4/4), so -0.0 is reproduced there.
        // QualByDepth divides the UNROUNDED double it reads back from the record
        // (`-10.0 * vc.getLog10PError()`, QualByDepth.java:78+86), while the
        // published QUAL token is that double rounded to two decimals and stored
        // as float32.  Dividing the rounded value instead can land exactly on a
        // two-decimal tie: measured on GATK's own chr20 corpus at 20:10068160,
        // the rounded float gives 9.3349997202555333 (rendering as 9.34 after the
        // float32 round-trip) where GATK publishes 9.33 from 9.3345911138702871.
        // `call_confidence` IS that unrounded -10*log10PError double.
        // Magnitudes indistinguishable from zero are excluded: there the sign of
        // zero comes from a ~1e-16 round-off in the AF calculator that GATK does
        // not reproduce consistently (rounds 51/53), so the rounded QUAL is kept.
        const bool confidence_usable = record.call_confidence_available &&
            std::abs(record.call_confidence) > 1.0e-9;
        const double qd_numerator =
            (record.materialized_spanning_locus && record.value->qual == 0.0F)
                ? -0.0
                : (confidence_usable ? record.call_confidence
                                     : static_cast<double>(record.value->qual));
        const float qd = static_cast<float>(gatk_fix_high_qd(
            qd_numerator / static_cast<double>(depth), random));
        if (bcf_update_info_float(output_header, record.value, "QD", &qd, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write QD annotation");
    }

    const auto table = merged_strand_bias_table(output_header, record);
    const auto total = table[0] + table[1] + table[2] + table[3];
    if (total <= 0) return;
    // FisherStrand uses MIN_COUNT=2.  GATK still materializes the standard
    // FS field as zero for a lower-depth table, so preserve that observable
    // output contract while avoiding a spurious low-count p-value.
    const float fs = total <= 2 ? 0.0F : static_cast<float>(std::min(
        999.0, -10.0 * std::log10(std::max(fisher_two_sided_pvalue(
            table[0], table[1], table[2], table[3]), 1.0e-320))));
    if (bcf_update_info_float(output_header, record.value, "FS", &fs, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write FS annotation");

    const double t00 = static_cast<double>(table[0]) + 1.0;
    const double t01 = static_cast<double>(table[1]) + 1.0;
    const double t10 = static_cast<double>(table[2]) + 1.0;
    const double t11 = static_cast<double>(table[3]) + 1.0;
    const double ratio = (t00 / t01) * (t11 / t10) +
                         (t01 / t00) * (t10 / t11);
    const double ref_ratio = std::min(t00, t01) / std::max(t00, t01);
    const double alt_ratio = std::min(t10, t11) / std::max(t10, t11);
    const float sor = static_cast<float>(std::log(ratio) +
                                         std::log(ref_ratio) -
                                         std::log(alt_ratio));
    if (bcf_update_info_float(output_header, record.value, "SOR", &sor, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write SOR annotation");
}

void apply_gatk_annotation_compatibility(const bcf_hdr_t* output_header,
                                         Record& record,
                                         const Options& options) {
    if (!options.gatk_annotation_compatibility) return;
    // GATK never carries a source record's FILTER into the regenotyped call.
    // GenotypingEngine.calculateGenotypes() rebuilds the variant with
    //     new VariantContextBuilder(callSourceString(), vc.getContig(),
    //                               vc.getStart(), vc.getEnd(), outputAlleles)
    // (GenotypingEngine.java:181): that constructor receives only the source
    // string, the coordinates and the output allele list, so no filter state is
    // copied from the source VariantContext, `filtersWereApplied` stays false,
    // and htsjdk's VCFEncoder writes the unfiltered token `.` (its
    // addFilterString emits "." exactly when !vc.filtersWereApplied()).  The
    // only FILTER a GenotypeGVCFs record can carry is the engine's own LowQual
    // from `if (!passesCallThreshold(phredScaledConfidence))
    // builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME)` (:184-186) -- a
    // test on the RECOMPUTED confidence, never on the source filter.
    //
    // Inheriting the source filter is wrong in both name and presence.  HTSlib
    // keeps a single ID dictionary shared by FILTER/INFO/FORMAT and resolves
    // record.value->d.flt[] against the WRITING header, while the source record
    // was decoded with a header in which an undeclared source FILTER had just
    // been auto-registered by vcf_parse_filter() -- after this output header was
    // duplicated from the header-only first pass.  Measured on the dense-mode
    // non-PASS fixture: the REF-only row printed `RGQ` (the first id appended by
    // add_genotype_output_header_fields) where GATK writes `.`, and the same
    // row printed `GQ` once FORMAT/GQ was already declared in the input header.
    // Clearing the inherited filters also canonicalizes FILTER=PASS to `.`,
    // which is what htsjdk does for an unfiltered call.
    bcf_unpack(record.value, BCF_UN_FLT);
    if (record.value->d.n_flt > 0 &&
        bcf_update_filter(output_header, record.value, nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear inherited FILTER");
    // ... and then GATK's OWN decision, which is the only filter a
    // GenotypeGVCFs record can carry:
    //
    //   if ( ! passesCallThreshold(phredScaledConfidence) ) {
    //       builder.filter(GATKVCFConstants.LOW_QUAL_FILTER_NAME);
    //   }
    //
    // (GenotypingEngine.java:184-186, passesCallThreshold() = `conf >=
    // standardConfidenceForCalling` at :430-432, LOW_QUAL_FILTER_NAME =
    // "LowQual" at GATKVCFConstants.java:179).  It must be applied AFTER the
    // clear above, otherwise the clear would wipe it.
    //
    // The comparison uses the RAW recomputed confidence recorded by the engine
    // stage, never `record.value->qual`: the published token is that double
    // rounded to two decimals, so the measured weak-locus row (confidence in
    // [23.135, 23.1358), QUAL column 23.14) is filtered at
    // --standard-min-confidence-threshold-for-calling 23.1358 while its QUAL
    // already reads 23.14.  `call_confidence` is the PRE-posterior value,
    // because :184-186 precedes the --use-posteriors-to-calculate-qual update
    // at :192-199.  Records that never reached the genotyping engine (a
    // passthrough reference-confidence block) carry no confidence and stay
    // unfiltered, which is what GATK emits for them.  The filter is
    // site-level: `builder.filter()` takes no allele, GenotypeGVCFs never calls
    // AlleleFilterUtils.addAlleleAndSiteFilters() (the per-allele intersection
    // at AlleleFilterUtils.java:115-120 belongs to VariantFiltration/Mutect2),
    // and `--invalidate-previous-filters` is not a GenotypeGVCFs option at all
    // (measured: GATK exits 1 with "is not a recognized option").
    if (record.call_confidence_available &&
        !(record.call_confidence >= options.standard_confidence_for_calling)) {
        int low_qual = bcf_hdr_id2int(output_header, BCF_DT_ID, "LowQual");
        if (low_qual < 0)
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: cannot write LowQual FILTER: the id is not "
                "declared in the output header");
        if (bcf_update_filter(output_header, record.value, &low_qual, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write LowQual FILTER");
    }
    // RAW_MQandDP is the reducible input annotation used to derive MQ; GATK's
    // finalized GenotypeGVCFs records do not retain it.  RCQ/RCP are native
    // cross-sample diagnostics rather than standard GATK INFO fields.
    (void)bcf_update_info_int32(output_header, record.value, "RAW_MQandDP", nullptr, 0);
    (void)bcf_update_info_float(output_header, record.value, "RCQ", nullptr, 0);
    (void)bcf_update_info_float(output_header, record.value, "RCP", nullptr, 0);
    // FORMAT/SB is useful for debugging and is retained by the default native
    // profile, but GenotypeGVCFs final output emits only the standard genotype
    // fields.  Remove it after FS/SOR have consumed the table.
    if (bcf_update_format_int32(output_header, record.value, "SB", nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove native SB annotation");
    record.sb.clear();
    // GenotypingEngine's PREFER_PLS allele-subsetting path removes posterior
    // attributes before it builds the final Genotype.  In particular this is
    // observable when an input gVCF carries FORMAT/GP (or GATK's PHRED-scaled
    // FORMAT/PG): dropping <NON_REF> must not leave a stale Number=G vector
    // indexed by the pre-subset allele list.  Clear both fields after all
    // posterior/annotation consumers have run; the diagnostic native profile
    // intentionally keeps them available for focused posterior experiments.
    if (bcf_update_format_float(output_header, record.value, "GP", nullptr, 0) != 0 ||
        bcf_update_format_float(output_header, record.value, "PG", nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot remove stale posterior FORMAT fields");
}

std::vector<std::string> split_tab_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        fields.emplace_back(line.substr(begin, end == std::string::npos ?
                                                std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

std::string format_gatk_float_value(const std::string& key,
                                    const std::string& value) {
    if (value.empty() || value == ".") return value;
    int precision = -1;
    // htsjdk's VCFEncoder formats Float/Double values with the following
    // value-dependent rule (formatVCFDouble): values below one use three
    // decimals, values at least one use two, and negative values below
    // 0.01 select scientific notation.  GenotypeGVCFs carries the rank-sum
    // annotations as numeric attributes, so applying this rule is required
    // for raw-text parity (for example -3.891e+00 and 1.10), rather than
    // blindly forcing the three-decimal annotation format used by the
    // standalone RankSumTest producer.
    const bool rank_sum = key == "BaseQRankSum" || key == "MQRankSum" ||
                          key == "ReadPosRankSum";
    if (key == "AF" || key == "MLEAF") precision = 3;
    else if (key == "FS") precision = 3;
    else if (key == "ExcessHet" || key == "InbreedingCoeff") precision = 4;
    else if (key == "MQ" || key == "QD") precision = 2;
    else if (key == "SOR") precision = 3;
    if (precision < 0 && !rank_sum) return value;
    std::ostringstream output;
    std::size_t begin = 0;
    bool first = true;
    while (begin <= value.size()) {
        const auto end = value.find(',', begin);
        const auto token = value.substr(begin, end == std::string::npos ?
                                                std::string::npos : end - begin);
        if (!first) output << ',';
        first = false;
        try {
            auto number = std::stod(token);
            if (std::isfinite(number)) {
                // htsjdk's VCFEncoder formats the sign of a negative zero, and
                // the materialized spanning loci above are the only rows whose
                // QD can be -0.0 (GATK publishes -0.00 there).  Every other
                // annotation keeps the historical native normalisation.
                if (number == 0.0 && key != "QD") number = 0.0;
                int token_precision = precision;
                if (rank_sum) {
                    // Exact htsjdk VCFEncoder.formatVCFDouble ordering:
                    // any value below 0.01 is either scientific notation or
                    // the literal 0.00, including a zero rank-sum.  The
                    // latter is observable after GenotypeGVCFs even though
                    // HaplotypeCaller itself emits its raw zero as 0.000.
                    if (number < 0.01) {
                        if (std::abs(number) >= 1.0e-20) {
                            output << std::scientific << std::setprecision(3)
                                   << number;
                            token_precision = -2;
                        } else {
                            output << "0.00";
                            token_precision = -2;
                        }
                    } else if (number < 1.0) {
                        token_precision = 3;
                    } else {
                        token_precision = 2;
                    }
                }
                // HTSJDK's allele-frequency annotations use two decimals for
                // exact endpoints (0/1) and three for interior frequencies.
                if (token_precision >= 0 && (key == "AF" || key == "MLEAF")) {
                    token_precision = (number <= 0.0 || number >= 1.0) ? 2 : 3;
                }
                // HTSlib stores INFO floats as float32, while the Java
                // VariantContext retains the decimal source value as a
                // double.  `vcf_format1()` reconstructs the original text
                // (for example 1.045), but binary ties round-to-even in the
                // C++ stream and incorrectly produce 1.04.  Nudge only the
                // rank-sum decimal formatting across its representational
                // tie so it follows Java Formatter's 1.045 -> 1.05 rule.
                if (rank_sum && token_precision >= 0) {
                    const auto epsilon = std::max(1.0, std::abs(number)) * 1.0e-12;
                    number += std::copysign(epsilon, number == 0.0 ? 1.0 : number);
                }
                if (token_precision >= 0)
                    output << std::fixed << std::setprecision(token_precision) << number;
            }
            else
                output << token;
        } catch (...) {
            output << token;
        }
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return output.str();
}

std::string format_gatk_qual_value(const std::string& value) {
    if (value.empty() || value == ".") return value;
    // htsjdk's VCFEncoder.formatQualValue() is `String.format(Locale.US,
    // "%.2f", qual)` with a trailing ".00" removed, and Java's `%f` prints an
    // infinite double as the literal word "Infinity" (measured on the pinned
    // JDK 17: "Infinity" / "-Infinity" / "NaN").  HTSlib's kputd() renders the
    // same float through C's "%.g", i.e. "inf"/"-inf" (htslib kstring.c:38,
    // reached from vcf.c:4109), so translate the spelling.  A finite value
    // keeps the decimal path below.
    if (value == "inf") return "Infinity";
    if (value == "-inf") return "-Infinity";
    try {
        auto quality = std::stod(value);
        if (!std::isfinite(quality)) return value;
        if (quality == 0.0) quality = 0.0; // eliminate a possible negative zero
        std::ostringstream output;
        output << std::fixed << std::setprecision(2) << quality;
        auto text = output.str();
        // htsjdk VCFEncoder.formatQualValue() renders two decimals then
        // removes a trailing .00, unlike INFO float formatting.
        if (text.size() >= 3 && text.ends_with(".00"))
            text.resize(text.size() - 3);
        return text;
    } catch (...) {
        return value;
    }
}

std::string gatk_compatible_record_text(const std::string& formatted) {
    if (formatted.empty()) return formatted;
    std::string line = formatted;
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
    auto fields = split_tab_fields(line);
    if (fields.size() < 8) return formatted;
    fields[5] = format_gatk_qual_value(fields[5]);

    struct InfoValue { std::string key; std::string value; bool has_value = false; };
    std::vector<InfoValue> info;
    if (fields[7] != ".") {
        std::size_t begin = 0;
        while (begin <= fields[7].size()) {
            const auto end = fields[7].find(';', begin);
            const auto token = fields[7].substr(begin, end == std::string::npos ?
                                                        std::string::npos : end - begin);
            const auto equal = token.find('=');
            if (equal == std::string::npos)
                info.push_back(InfoValue{token, std::string{}, false});
            else
                info.push_back(InfoValue{token.substr(0, equal), token.substr(equal + 1), true});
            if (end == std::string::npos) break;
            begin = end + 1;
        }
    }
    const std::vector<std::string> standard_order{
        // Match the standard GenotypeGVCFs annotation order.  Rank-sum fields
        // are carried through when present; their full read-level recalculation
        // remains a separate GenotypeGVCFs annotation-engine task.
        "AC", "AF", "AN", "BaseQRankSum", "DP", "ExcessHet", "FS",
        "InbreedingCoeff", "MLEAC", "MLEAF", "MQ", "MQRankSum", "NDA", "QD",
        "ReadPosRankSum", "SOR"};
    std::vector<InfoValue> ordered;
    std::vector<bool> used(info.size(), false);
    for (const auto& key : standard_order) {
        for (std::size_t index = 0; index < info.size(); ++index) {
            if (!used[index] && info[index].key == key) {
                ordered.push_back(info[index]);
                used[index] = true;
                break;
            }
        }
    }
    for (std::size_t index = 0; index < info.size(); ++index)
        if (!used[index]) ordered.push_back(info[index]);
    if (ordered.empty()) {
        fields[7] = ".";
    } else {
        std::ostringstream info_text;
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            if (index != 0) info_text << ';';
            info_text << ordered[index].key;
            if (ordered[index].has_value) {
                info_text << '=' << format_gatk_float_value(ordered[index].key,
                                                             ordered[index].value);
            }
        }
        fields[7] = info_text.str();
    }

    // HTSJDK's VariantContextWriter emits the standard genotype fields in a
    // stable order independent of the order in which a source VCF declared
    // them.  Native HTSlib records retain the merged-header order instead;
    // normalize the textual compatibility profile here while keeping the
    // staged binary representation untouched.  This is important for direct
    // GATK replacement because row-level oracles and downstream tools often
    // compare the FORMAT/sample columns byte-for-byte.
    if (fields.size() >= 9 && fields[8] != "." && !fields[8].empty()) {
        std::vector<std::string> format_keys;
        std::size_t begin = 0;
        while (begin <= fields[8].size()) {
            const auto end = fields[8].find(':', begin);
            format_keys.emplace_back(fields[8].substr(
                begin, end == std::string::npos ? std::string::npos : end - begin));
            if (end == std::string::npos) break;
            begin = end + 1;
        }
        const std::vector<std::string> format_order{
            "GT", "AD", "DP", "GQ", "RGQ", "PGT", "PID", "PL", "PS",
            "PP", "GP", "PG"};
        std::vector<std::size_t> format_indices;
        std::vector<std::string> ordered_format_keys;
        std::vector<bool> format_used(format_keys.size(), false);
        for (const auto& key : format_order) {
            for (std::size_t index = 0; index < format_keys.size(); ++index) {
                if (!format_used[index] && format_keys[index] == key) {
                    ordered_format_keys.push_back(format_keys[index]);
                    format_indices.push_back(index);
                    format_used[index] = true;
                    break;
                }
            }
        }
        for (std::size_t index = 0; index < format_keys.size(); ++index) {
            if (!format_used[index]) {
                ordered_format_keys.push_back(format_keys[index]);
                format_indices.push_back(index);
            }
        }
        std::ostringstream format_text;
        for (std::size_t index = 0; index < ordered_format_keys.size(); ++index) {
            if (index != 0) format_text << ':';
            format_text << ordered_format_keys[index];
        }
        fields[8] = format_text.str();
        for (std::size_t sample = 9; sample < fields.size(); ++sample) {
            std::vector<std::string> values;
            begin = 0;
            while (begin <= fields[sample].size()) {
                const auto end = fields[sample].find(':', begin);
                values.emplace_back(fields[sample].substr(
                    begin, end == std::string::npos ? std::string::npos : end - begin));
                if (end == std::string::npos) break;
                begin = end + 1;
            }
            std::ostringstream sample_text;
            for (std::size_t index = 0; index < format_indices.size(); ++index) {
                if (index != 0) sample_text << ':';
                const auto source = format_indices[index];
                sample_text << (source < values.size() ? values[source] : ".");
            }
            fields[sample] = sample_text.str();
        }
    }
    std::ostringstream result;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) result << '\t';
        result << fields[index];
    }
    result << '\n';
    return result.str();
}

int write_vcf_text_line(htsFile* output, const std::string& line) {
    if (output->is_bgzf)
        return bgzf_write(output->fp.bgzf, line.data(), line.size()) ==
                       static_cast<ssize_t>(line.size()) ? 0 : -1;
    return hwrite(output->fp.hfile, line.data(), line.size()) ==
                   static_cast<ssize_t>(line.size()) ? 0 : -1;
}

// GATK's GenotypeGVCFs declares exactly this one filter line for every run of
// the tool.  GenotypeGVCFsEngine.setupVCFWriter() ends its header construction
// with `headerLines.add(GATKVCFHeaderLines.getFilterLine(
// GATKVCFConstants.LOW_QUAL_FILTER_NAME));` (GenotypeGVCFsEngine.java:416) --
// the last add() before the VCFHeader is built (:418-419) -- and that call is
// unguarded, so it is present whether or not the input declared any filter.
// The text is GATKVCFHeaderLines.java:89
// (`addFilterLine(new VCFFilterHeaderLine(LOW_QUAL_FILTER_NAME, "Low quality"))`)
// with the id from GATKVCFConstants.java:179 (LOW_QUAL_FILTER_NAME = "LowQual").
// It is the only filter the engine can apply to an output record
// (GenotypingEngine.java:184-186).
constexpr const char* kGatkLowQualFilterLine =
    "##FILTER=<ID=LowQual,Description=\"Low quality\">";

// GATK declares that filter unconditionally in every GenotypeGVCFs header
// (`headerLines.add(GATKVCFHeaderLines.getFilterLine(
// GATKVCFConstants.LOW_QUAL_FILTER_NAME))` is the last add() of
// GenotypeGVCFsEngine.setupVCFWriter(), :416, with the text from
// GATKVCFHeaderLines.java:89).  The declaration has to be in the bcf_hdr_t and
// not only in the header text written by gatk_compatible_header_text(),
// because GenotypingEngine.java:184-186 puts that id into a record's FILTER
// column and HTSlib resolves the id against this header.  Appending it here
// cannot change the emitted header text:
// gatk_compatible_header_text() re-collects the ##FILTER group from whatever
// the header contains, de-duplicates against this exact line, sorts the group
// and re-inserts it at the position the input's own filter lines occupied.
// Both output writers build their header through this call (the streaming
// setup via add_genotype_output_header_fields, the aggregate setup inline).
void ensure_gatk_low_qual_filter_declaration(bcf_hdr_t* output_header,
                                             const Options& options) {
    if (!options.gatk_annotation_compatibility) return;
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "LowQual") >= 0) return;
    bcf_hdr_append(output_header, kGatkLowQualFilterLine);
}

// HTSlib's synthetic PASS declaration, injected into every header it parses.
constexpr const char* kHtslibSyntheticPassFilterLine =
    "##FILTER=<ID=PASS,Description=\"All filters passed\">";

// htsjdk re-serializes every header line it parsed.  VCFHeaderLineTranslator
// .parseLine() builds a VCFCompoundHeaderLine out of the attributes of a
// `##INFO`/`##FORMAT`/`##FILTER`/`##ALT` line and the writer emits the line
// again out of that object, so a Description the input left unquoted comes back
// quoted: measured on the gate fixture, GATK writes
// `##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">` for an input
// that said `Description=Read depth`, and the same for `##ALT`,
// `##FORMAT=<ID=GT...>` and the rest.  The attribute order inside `<>` is NOT
// re-ordered -- htsjdk refuses an input that does not already use its
// (ID, Number, Type, Description) order ("Tag Description in wrong order (was
// #2, expected #4)", VCF4Parser.parseLine) -- so re-quoting is the whole of the
// re-serialization a valid input can need.  Only the measured kinds are
// touched; a line that is already quoted, or that has no Description, is
// returned unchanged.
std::string gatk_canonical_header_line(const std::string& line) {
    const bool compound =
        line.rfind("##INFO=<", 0) == 0 || line.rfind("##FORMAT=<", 0) == 0 ||
        line.rfind("##FILTER=<", 0) == 0 || line.rfind("##ALT=<", 0) == 0;
    if (!compound || line.size() < 3 || line.back() != '>') return line;
    static const std::string marker = "Description=";
    const auto at = line.find(marker);
    if (at == std::string::npos) return line;
    const auto value = at + marker.size();
    if (value >= line.size() - 1 || line[value] == '"') return line;
    return line.substr(0, value) + '"' +
           line.substr(value, line.size() - 1 - value) + "\">";
}

// GATK's own MLEAC/MLEAF lines (GATKVCFHeaderLines.java:148-149, added by
// GenotypeGVCFsEngine.setupVCFWriter() at :404-405).  Native's diagnostic
// (non-GATK) profile appends a short description of its own for these two keys;
// the GATK-compatibility profile must declare GATK's text instead, so these
// constants exist for the two places that append them.
constexpr const char* kGatkMleacInfoLine =
    "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">";
constexpr const char* kGatkMleafInfoLine =
    "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">";

// GATK's standard INFO/DP description (htsjdk VCFStandardHeaderLines, added by
// GenotypeGVCFsEngine.setupVCFWriter() at :407 -- "needed for gVCFs without DP
// tags").  It is deliberately DIFFERENT text from the htsjdk FORMAT/DP line
// ("Read depth"), which is why GATK's header carries two ##INFO=<ID=DP,...>
// lines when the input declares only the FORMAT/INFO pair this fixture has.
constexpr const char* kGatkStandardDpInfoLine =
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">";

// GATK's NDA line (GATKVCFHeaderLines.java:205), added by
// GenotypingEngine.getAppropriateVCFInfoHeaders() (GenotypingEngine.java:90-94)
// when --annotate-with-num-discovered-alleles is set.  Native's diagnostic
// profile appends a shorter description of its own.
constexpr const char* kGatkNdaInfoLine =
    "##INFO=<ID=NDA,Number=1,Type=Integer,Description=\"Number of alternate alleles discovered (but not necessarily genotyped) at this site\">";

// Every header line GATK 4.6.2.0 GenotypeGVCFs declares for itself, i.e. the
// lines of its output header that do not come from the input.  Measured byte for
// byte on the gate fixture of
// fastgatk-native/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py and cited
// to the source that produces each one:
//   * GenotypeGVCFsEngine.setupVCFWriter() adds
//     `annotationEngine.getVCFAnnotationDescriptions(false)` (:401) -- the
//     annotation lines, resolved through GATKVCFHeaderLines (or htsjdk's
//     VCFStandardHeaderLines for AD/DP, GATKVCFHeaderLines.java:18-43),
//     `genotypingEngine.getAppropriateVCFInfoHeaders()` (:402), MLEAC/MLEAF/RGQ
//     (:404-406), the standard INFO/DP line (:407) and finally the LowQual
//     filter (:416);
//   * the annotation set is the StandardAnnotation group
//     (GenotypeGVCFs.java:255-257), which is unconditional.
// GATK adds these to a Set<VCFHeaderLine> that de-duplicates by the whole line
// (ID + Number + Type + Description), not by ID, so a line is only absent from
// the output header when the input declared exactly the same text.  Measured:
// GATK's header carries BOTH `##INFO=<ID=DP,...Description="Read depth">` and
// `##INFO=<ID=DP,...Description="Approximate read depth; some reads may have
// been filtered">`, and both `##FORMAT=<ID=AD,...>` spellings, while the
// canonical MLEAC/MLEAF pair of a HaplotypeCaller gVCF appears exactly once.
const std::vector<std::string>& gatk_declared_info_lines() {
    static const std::vector<std::string> lines{
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
        "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">",
        "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">",
        "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">",
        "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">",
        "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">",
        "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">",
        "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">",
        "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">",
        "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">",
        "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">"};
    return lines;
}

// The FORMAT half of the same list: the standard AD line GATK's header carries
// because DepthPerAlleleBySample is a StandardAnnotation (its description comes
// from VariantAnnotation.java:22-33 -> GATKVCFHeaderLines.getFormatLine("AD",
// true) -> htsjdk's VCFStandardHeaderLines), and the RGQ line the writer adds
// itself (GenotypeGVCFsEngine.java:406 + GATKVCFHeaderLines.java:133).
const std::vector<std::string>& gatk_declared_format_lines() {
    static const std::vector<std::string> lines{
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">"};
    return lines;
}

// GATK's GenotypeGVCFs propagates the input header's own ##FILTER lines
// verbatim: the writer seeds its header set from the input's lines
// (`final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(
// inputVCFHeader.getMetaDataInInputOrder())`, GenotypeGVCFsEngine.java:395) and
// the only line it removes is the GVCF block band (`headerLines.removeIf(
// vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK))`, :398-399),
// so a `##FILTER=<ID=PASS,...>` the input declares is written back out with the
// input's own text.  GATK never synthesizes that line: neither the pinned
// gatk-package jar nor htsjdk 4.2.0 contains the string "All filters passed",
// and GATK's own PASS line (GATKVCFHeaderLines.java:90, "Site contains at least
// one allele that passes filters") is used only by VQSR.
//
// HTSlib cannot report the input's line.  bcf_hdr_parse() injects a synthetic
// `##FILTER=<ID=PASS,Description="All filters passed">` record as the second
// record of EVERY header it parses -- "The filter PASS must appear first in the
// dictionary" (third_party/htslib-build/htslib-src/vcf.c:1271-1276) -- and
// bcf_hdr_add_hrec()/bcf_hdr_register_hrec() then destroy any later ##FILTER
// record whose ID is already registered (vcf.c:1092-1112), so the parsed header
// carries HTSlib's canonical text whether or not the input declared PASS, and
// loses the input's own text when it did.  Measured: an input with no ##FILTER
// line at all still yields that line from bcf_hdr_format(), and an input
// declaring `Description="Some other PASS text"` yields the canonical text.
// Reading the input file's own text back is therefore the only way to tell
// "the input declared PASS" from "HTSlib injected it".
//
// Returns an empty string when the input has no readable text header (no input,
// stdin, BCF, an unreadable or GenomicsDB-expanded path), in which case the
// caller falls back to the previous behaviour of dropping the synthetic line.
std::string input_pass_filter_line(const std::string& path) {
    if (path.empty() || path == "-") return {};
    htsFile* input = hts_open(path.c_str(), "r");
    if (input == nullptr) return {};
    std::string found;
    // A BCF header is binary and its text form is HTSlib's reconstruction
    // anyway, so only a text-format input can answer this question.
    if (input->format.format == vcf) {
        kstring_t line{0, 0, nullptr};
        while (fastgatk::io::read_text_line(input, &line, path) >= 0) {
            if (line.l == 0) continue;
            const std::string text(line.s, line.l);
            if (text.rfind("#CHROM", 0) == 0) break;
            if (text.rfind("##FILTER=<ID=PASS,", 0) == 0 ||
                text.rfind("##FILTER=<ID=PASS>", 0) == 0) {
                found = text;
                break;
            }
        }
        free(line.s);
    }
    hts_close(input);
    return found;
}

// htsjdk writes a VCF header in SORTED order, so reproducing GATK's output
// header byte for byte means reproducing htsjdk's header sort as well as its
// content.  The comparator, read out of the pinned jar with
// `javap -p -c third_party/gatk-package/gatk-4.6.2.0/
// gatk-package-4.6.2.0-local.jar ...`, is:
//
//   htsjdk.variant.vcf.VCFHeaderLine.compareTo (every line type but contigs)
//     return toString().compareTo(o.toString());
//   htsjdk.variant.vcf.VCFContigHeaderLine.compareTo (overrides it)
//     if (o instanceof VCFContigHeaderLine)
//         return contigIndex.compareTo(((VCFContigHeaderLine) o).contigIndex);
//     return super.compareTo(o);            // == toString().compareTo(...)
//   htsjdk.variant.vcf.VCFHeader.getMetaDataInSortedOrder()
//     return makeGetMetaDataSet(new TreeSet<>(mMetaData));
//   htsjdk.variant.variantcontext.writer.VCFWriter.writeHeader(header, writer,
//                                                              version, name)
//     writer.write("##fileformat=" + version + "\n");
//     for (VCFHeaderLine line : header.getMetaDataInSortedOrder()) {
//         if (VCFHeaderVersion.isFormatString(line.getKey())) continue;
//         writer.write("##" + line + "\n");
//     }
//     writer.write("#" + HEADER_FIELDS + samples + "\n");
//
// So, precisely:
//   1. `##fileformat` is written first, out of the header version, and is NOT
//      part of the sort;
//   2. every other line is ordered by VCFHeaderLine.toString() -- the full line
//      text without the leading `##` (`key + "=" + value`) -- compared as Java
//      String.compareTo, i.e. by key first and then by the rest of the line
//      (ASCII byte order == UTF-16 code-unit order for every header line seen);
//   3. two contig lines compare by `contigIndex`, the index the line was given
//      when the INPUT header was parsed, i.e. the input's own contig
//      declaration order -- measured: an input declaring chr1, chr10, chr2
//      against a reference dictionary ordered chr2, chr10, chr1 comes back
//      chr1, chr10, chr2, and the reverse.  A contig line against a non-contig
//      line falls back to plain text, so the contig block keeps the input's
//      order but moves to where `##contig=...` sorts among the rest: after
//      every `##INFO=` and before `##source=`;
//   4. `#CHROM` is written last, after the sorted block.
// `##GATKCommandLine` and `##source` are ordinary lines in that sort (measured:
// GATKCommandLine between `##FORMAT=` and `##INFO=`, source last); native does
// not write a GATKCommandLine, and `##source` is already placed by the sort.
std::vector<std::string> gatk_htsjdk_sorted_header(std::vector<std::string> lines) {
    const std::size_t total = lines.size();
    std::vector<std::string> fileformat;
    std::vector<std::string> chrom;
    std::vector<std::string> contigs;
    std::vector<std::string> others;
    for (auto& line : lines) {
        if (line.rfind("##fileformat=", 0) == 0) {
            fileformat.push_back(std::move(line));
        } else if (line.rfind("#CHROM", 0) == 0) {
            chrom.push_back(std::move(line));
        } else if (line.rfind("##contig=", 0) == 0) {
            contigs.push_back(std::move(line));
        } else {
            others.push_back(std::move(line));
        }
    }
    std::sort(others.begin(), others.end());
    // Every contig line begins with this prefix and no non-contig line does, so
    // the first non-contig line that sorts after the prefix is exactly where the
    // contig block belongs (compareTo against a contig line is a text compare,
    // and all contig lines agree up to and including the prefix).
    const std::string contig_prefix = "##contig=";
    const auto at = std::lower_bound(others.begin(), others.end(), contig_prefix);
    std::vector<std::string> sorted;
    sorted.reserve(total);
    sorted.insert(sorted.end(), fileformat.begin(), fileformat.end());
    sorted.insert(sorted.end(), others.begin(), at);
    sorted.insert(sorted.end(), contigs.begin(), contigs.end());
    sorted.insert(sorted.end(), at, others.end());
    sorted.insert(sorted.end(), chrom.begin(), chrom.end());
    return sorted;
}

std::string gatk_compatible_header_text(const std::string& formatted,
                                        const std::string& input_pass_filter) {
    std::vector<std::string> lines;
    std::size_t begin = 0;
    while (begin <= formatted.size()) {
        const auto end = formatted.find('\n', begin);
        auto line = formatted.substr(begin, end == std::string::npos ?
                                             std::string::npos : end - begin);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) lines.push_back(std::move(line));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    std::vector<std::string> info_lines;
    std::vector<std::string> retained;
    bool inserted_info = false;
    bool inserted_source = false;
    for (const auto& raw_line : lines) {
        // htsjdk writes its header back out of parsed objects, so every
        // compound line it retained is re-quoted (see
        // gatk_canonical_header_line() above).  Everything below works on the
        // canonical spelling, which is what GATK's output carries.
        const auto line = gatk_canonical_header_line(raw_line);
        if (line.rfind("##fastgatk_genotype_gvcfs_status=", 0) == 0) continue;
        // GVCF block-band declarations are input-only metadata.  GATK's
        // GenotypeGVCFs writer does not propagate them to the materialized
        // VCF, even though HTSlib's duplicated input header would otherwise
        // retain them verbatim.
        if (line.rfind("##GVCFBlock", 0) == 0) continue;
        // HTSlib's own synthetic PASS declaration, not the input's: bcf_hdr_parse()
        // inserts `##FILTER=<ID=PASS,Description="All filters passed">` into every
        // header it parses (see input_pass_filter_line() below), so this line is
        // present here whether or not the input declared it.  The input's own
        // PASS line, when it has one, is re-added from the raw input text.
        if (line == kHtslibSyntheticPassFilterLine) continue;
        if (line.rfind("##INFO=<ID=RCQ,", 0) == 0 ||
            line.rfind("##INFO=<ID=RCP,", 0) == 0) continue;
        if (line.rfind("##INFO=", 0) == 0) {
            info_lines.push_back(line);
            continue;
        }
        if (line == "##source=HaplotypeCaller" && !inserted_source) {
            retained.emplace_back("##source=GenotypeGVCFs");
            inserted_source = true;
        }
        retained.push_back(line);
    }
    const std::vector<std::string> info_order{
        "AC", "AF", "AN", "BaseQRankSum", "DP", "END", "ExcessHet", "FS",
        "InbreedingCoeff", "MLEAC", "MLEAF", "MQ", "MQRankSum", "QD",
        "RAW_MQandDP", "ReadPosRankSum", "SOR"};
    // GATK's own declarations join the group the input's lines are already in;
    // a line is added only when its exact text is absent, because GATK's header
    // set de-duplicates by the whole line and not by ID (see
    // gatk_declared_info_lines() above).  This is what makes GATK's duplicate
    // ##INFO=<ID=DP,...>/##FORMAT=<ID=AD,...> pairs come out right: the input's
    // own spelling is a different line and is kept alongside the standard one.
    for (const auto& declared : gatk_declared_info_lines())
        if (std::find(info_lines.begin(), info_lines.end(), declared) == info_lines.end())
            info_lines.push_back(declared);
    std::stable_sort(info_lines.begin(), info_lines.end(), [&](const auto& left, const auto& right) {
        const auto id = [](const std::string& line) {
            const auto marker = line.find("ID=");
            if (marker == std::string::npos) return std::string{};
            const auto start = marker + 3;
            const auto end = line.find_first_of(",>", start);
            return line.substr(start, end == std::string::npos ? std::string::npos : end - start);
        };
        const auto left_id = id(left);
        const auto right_id = id(right);
        const auto left_it = std::find(info_order.begin(), info_order.end(), left_id);
        const auto right_it = std::find(info_order.begin(), info_order.end(), right_id);
        const auto left_rank = left_it == info_order.end() ? info_order.size() :
                               static_cast<std::size_t>(left_it - info_order.begin());
        const auto right_rank = right_it == info_order.end() ? info_order.size() :
                                static_cast<std::size_t>(right_it - info_order.begin());
        return left_rank < right_rank;
    });
    std::vector<std::string> output;
    output.reserve(retained.size() + info_lines.size() + 2);
    // The ##FILTER group.  htsjdk writes a VCFHeader in sorted order, so GATK's
    // group is ordered by the lines' own text (measured: LowQual precedes an
    // input `q10` declaration); adding GATK's unconditional LowQual declaration
    // and sorting the group reproduces that order for any input, while leaving
    // the group's position in the header where the input put its filter lines.
    std::vector<std::string> filter_lines;
    for (const auto& line : retained)
        if (line.rfind("##FILTER=", 0) == 0) filter_lines.push_back(line);
    // The input's own PASS declaration, when the input file really has one.  It
    // is read from the raw input text because HTSlib replaced it with its
    // synthetic record (input_pass_filter_line() above); GATK keeps it verbatim,
    // so it goes back into the group with the input's own Description and is
    // ordered with the rest by the sort below.  It is canonicalized like every
    // other retained line: htsjdk re-quotes the Description.
    const auto own_pass_filter = gatk_canonical_header_line(input_pass_filter);
    if (!own_pass_filter.empty() &&
        std::find(filter_lines.begin(), filter_lines.end(), own_pass_filter) ==
            filter_lines.end())
        filter_lines.emplace_back(own_pass_filter);
    if (std::find(filter_lines.begin(), filter_lines.end(), kGatkLowQualFilterLine) ==
        filter_lines.end())
        filter_lines.emplace_back(kGatkLowQualFilterLine);
    std::sort(filter_lines.begin(), filter_lines.end());
    std::vector<std::string> format_lines;
    for (const auto& line : retained)
        if (line.rfind("##FORMAT=", 0) == 0) format_lines.push_back(line);
    const auto format_id = [](const std::string& line) {
        const auto marker = line.find("ID=");
        if (marker == std::string::npos) return std::string{};
        const auto start = marker + 3;
        const auto end = line.find_first_of(",>", start);
        return line.substr(start, end == std::string::npos ? std::string::npos : end - start);
    };
    const auto has_rgq = std::any_of(format_lines.begin(), format_lines.end(),
                                     [&](const auto& line) { return format_id(line) == "RGQ"; });
    if (!has_rgq) {
        format_lines.emplace_back(
            "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">");
    }
    // The FORMAT half of GATK's own declarations (the standard AD line; RGQ is
    // already ensured above).  Same exact-text de-duplication rule as the INFO
    // group.
    for (const auto& declared : gatk_declared_format_lines())
        if (std::find(format_lines.begin(), format_lines.end(), declared) ==
            format_lines.end())
            format_lines.push_back(declared);
    std::stable_sort(format_lines.begin(), format_lines.end(), [&](const auto& left,
                                                                    const auto& right) {
        return format_id(left) < format_id(right);
    });
    bool inserted_filter = false;
    bool inserted_format = false;
    for (const auto& line : retained) {
        if (line.rfind("##FILTER=", 0) == 0) {
            if (!inserted_filter) {
                output.insert(output.end(), filter_lines.begin(), filter_lines.end());
                inserted_filter = true;
            }
            continue;
        }
        if (line.rfind("##FORMAT=", 0) == 0) {
            if (!inserted_format) {
                output.insert(output.end(), format_lines.begin(), format_lines.end());
                inserted_format = true;
            }
            continue;
        }
        if (line.rfind("#CHROM", 0) == 0 && !inserted_info) {
            output.insert(output.end(), info_lines.begin(), info_lines.end());
            inserted_info = true;
        }
        output.push_back(line);
    }
    if (!inserted_filter) {
        // The input declared no filter at all, so the group has no input
        // position to inherit: it goes where a ##FILTER line would have been in
        // the input's own ordering, immediately before the first ##FORMAT= line
        // (or before #CHROM when the input declared no FORMAT line either).
        const auto anchor = std::find_if(output.begin(), output.end(), [](const auto& line) {
            return line.rfind("##FORMAT=", 0) == 0 || line.rfind("#CHROM", 0) == 0;
        });
        output.insert(anchor, filter_lines.begin(), filter_lines.end());
    }
    if (!inserted_format) {
        const auto chrom = std::find_if(output.begin(), output.end(), [](const auto& line) {
            return line.rfind("#CHROM", 0) == 0;
        });
        output.insert(chrom, format_lines.begin(), format_lines.end());
    }
    if (!inserted_info) output.insert(output.end(), info_lines.begin(), info_lines.end());
    if (!inserted_source) {
        const auto chrom = std::find_if(output.begin(), output.end(), [](const auto& line) {
            return line.rfind("#CHROM", 0) == 0;
        });
        output.insert(chrom, "##source=GenotypeGVCFs");
    }
    // GATK's output header is htsjdk's SORTED header, not the input's order with
    // the tool's own lines appended, so the position the group insertions above
    // produced is not the position GATK writes: the whole header is re-sorted
    // exactly the way htsjdk's writer sorts it (see
    // gatk_htsjdk_sorted_header() above).  The insertions above still decide
    // which lines are present -- in particular the duplicate `##INFO=<ID=DP,...>`
    // and `##FORMAT=<ID=AD,...>` pairs, whose relative order at this point does
    // not matter because the sort orders them by their own text, as GATK's
    // `Description="Approximate read depth; ..."` precedes `Description="Read
    // depth"`.
    output = gatk_htsjdk_sorted_header(std::move(output));
    std::ostringstream result;
    for (const auto& line : output) result << line << '\n';
    return result.str();
}

// A lazy GenotypeGVCFs cursor owns one HTSlib stream and only the records that
// have not yet reached the current joint locus.  The aggregate implementation
// below remains the regression/oracle path; this cursor is deliberately kept
// separate so the streaming contract cannot accidentally grow a second full
// genome-sized vector.
struct GenotypeStreamCursor {
    GenotypeStreamCursor() = default;
    htsFile* input = nullptr;
    bcf_hdr_t* header = nullptr;
    bcf1_t* raw = nullptr;
    // When an interval selector and a compatible CSI/TBI sidecar are
    // available, stream mode advances this cursor through one indexed query
    // at a time instead of decoding the whole source file.  The iterator is
    // kept per source so the joint heap still owns at most one decoded locus
    // from each shard.
    std::unique_ptr<GenotypeTraversalIndex> traversal_index;
    hts_itr_t* iterator = nullptr;
    kstring_t indexed_line{0, 0, nullptr};
    std::size_t region_index = 0;
    std::vector<int> output_samples;
    std::deque<Record> pending;
    // Lazy state for a reference-confidence block that crosses one or more
    // concrete variant spans.  Only the original block and the next emitted
    // segment are resident; the span index itself is shared/read-only.
    std::optional<Record> split_block;
    std::size_t split_span_index = 0;
    int split_cursor = 0;
    bool dense_reference_block = false;
    bool have_last_coordinate = false;
    int last_rid = -1;
    int last_pos = -1;

    ~GenotypeStreamCursor() {
        if (iterator != nullptr) hts_itr_destroy(iterator);
        free(indexed_line.s);
        for (auto& record : pending) bcf_destroy(record.value);
        pending.clear();
        if (split_block.has_value()) bcf_destroy(split_block->value);
        traversal_index.reset();
        bcf_destroy(raw);
        if (header != nullptr) bcf_hdr_destroy(header);
        if (input != nullptr) bcf_close(input);
    }
    GenotypeStreamCursor(const GenotypeStreamCursor&) = delete;
    GenotypeStreamCursor& operator=(const GenotypeStreamCursor&) = delete;
};

// Read one raw record for a streaming cursor.  Indexed VCF text requires the
// small `tbx_itr_next` + `vcf_parse` bridge, while BCF/CSI can use
// `bcf_itr_next` directly.  Regions are normalized before cursors are built,
// so a record cannot be returned twice merely because two user selectors
// overlap.  The caller still applies `in_regions`: tabix/CSI index keys are
// record starts, whereas GVCF END can extend a reference block into a query.
bool read_next_genotype_cursor_record(GenotypeStreamCursor& cursor,
                                      const std::vector<Region>& regions,
                                      std::uint64_t& indexed_interval_queries) {
    if (cursor.traversal_index == nullptr) {
        // htslib signals end of input with -1 and an unparseable record with
        // -2 (vcf_parse's error return).  Collapsing the two into "no more
        // records" would silently truncate the traversal at the first
        // malformed record, so the parse failure is reported instead.
        const int status = bcf_read(cursor.input, cursor.header, cursor.raw);
        if (status < -1)
            throw std::runtime_error(
                "BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed");
        return status == 0;
    }

    while (true) {
        if (cursor.iterator != nullptr) {
            int status = -1;
            if (cursor.traversal_index->tabix != nullptr) {
                status = tbx_itr_next(cursor.input, cursor.traversal_index->tabix,
                                      cursor.iterator, &cursor.indexed_line);
                if (status >= 0) {
                    if (vcf_parse(&cursor.indexed_line, cursor.header, cursor.raw) < 0)
                        throw std::runtime_error(
                            "BAD_INPUT: indexed streaming GenotypeGVCFs VCF record parse failed");
                    return true;
                }
            } else {
                status = bcf_itr_next(cursor.input, cursor.iterator, cursor.raw);
                if (status >= 0) return true;
            }
            hts_itr_destroy(cursor.iterator);
            cursor.iterator = nullptr;
            if (status < -1)
                throw std::runtime_error(
                    "BAD_INPUT: indexed streaming GenotypeGVCFs traversal failed");
        }

        if (cursor.region_index >= regions.size()) return false;
        const auto& region = regions[cursor.region_index++];
        ++indexed_interval_queries;
        const int index_tid = cursor.traversal_index->tid_for(
            cursor.header, region.rid, region.contig);
        if (index_tid < 0) continue;
        if (cursor.traversal_index->tabix != nullptr) {
            cursor.iterator = tbx_itr_queryi(cursor.traversal_index->tabix, index_tid,
                                             region.begin, region.end);
        } else {
            cursor.iterator = bcf_itr_queryi(cursor.traversal_index->index, index_tid,
                                             region.begin, region.end);
        }
        // A valid index may still have no chunks for this contig/interval.
        // Move on to the next normalized selector without treating it as an
        // input error.
    }
}

struct GenotypeStreamHeapItem {
    Record record;
    std::size_t source = 0;
};

struct GenotypeStreamHeapCompare {
    bool operator()(const GenotypeStreamHeapItem& left,
                    const GenotypeStreamHeapItem& right) const {
        if (left.record.rid != right.record.rid)
            return left.record.rid > right.record.rid;
        if (left.record.pos != right.record.pos)
            return left.record.pos > right.record.pos;
        if (left.record.key != right.record.key)
            return left.record.key > right.record.key;
        return left.source > right.source;
    }
};

[[maybe_unused]] std::vector<Record> split_stream_reference_block(Record original,
                                                  const bcf_hdr_t* output_header,
                                                  faidx_t* reference_index,
                                                  const std::vector<VariantSpan>& spans) {
    const int block_begin = original.pos;
    const int block_end = record_span_end(output_header, original);
    std::vector<VariantSpan> overlaps;
    for (const auto& span : spans) {
        if (span.rid != original.rid) continue;
        if (span.end <= block_begin) continue;
        if (span.begin >= block_end) break;
        overlaps.push_back(VariantSpan{span.rid, std::max(block_begin, span.begin),
                                       std::min(block_end, span.end)});
    }
    if (overlaps.empty() || block_end <= block_begin + 1)
        return {std::move(original)};
    if (reference_index == nullptr)
        throw std::runtime_error(
            "UNSUPPORTED_PARAMETER: stream-by-locus cross-sample reference-block splitting requires an indexed -R reference");

    std::vector<VariantSpan> merged;
    for (const auto& span : overlaps) {
        if (span.end <= span.begin) continue;
        if (merged.empty() || span.begin > merged.back().end)
            merged.push_back(span);
        else
            merged.back().end = std::max(merged.back().end, span.end);
    }
    std::vector<std::pair<int, int>> segments;
    int cursor = block_begin;
    for (const auto& span : merged) {
        if (span.begin > cursor) segments.emplace_back(cursor, span.begin);
        if (span.begin >= block_begin && span.begin < block_end)
            segments.emplace_back(span.begin, span.begin + 1);
        cursor = std::max(cursor, span.end);
    }
    if (cursor < block_end) segments.emplace_back(cursor, block_end);
    if (segments.empty()) {
        bcf_destroy(original.value);
        original.value = nullptr;
        return {};
    }

    std::vector<Record> result;
    result.reserve(segments.size());
    for (const auto [segment_begin, segment_end] : segments) {
        auto* copy = bcf_dup(original.value);
        if (copy == nullptr)
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot split streaming reference block");
        copy->rid = original.rid;
        copy->pos = segment_begin;
        const auto reference = reference_base(const_cast<faidx_t*>(reference_index),
                                              output_header, original.rid, segment_begin);
        const auto alleles = reference + ",<NON_REF>";
        if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
            bcf_destroy(copy);
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: cannot rewrite streaming reference block alleles");
        }
        if (segment_end - segment_begin > 1) {
            const int32_t end_value = segment_end;
            if (bcf_update_info_int32(output_header, copy, "END", &end_value, 1) != 0) {
                bcf_destroy(copy);
                throw std::runtime_error(
                    "OUTPUT_CONTRACT_FAILURE: cannot write streaming reference block END");
            }
        } else {
            (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
        }
        bcf_unpack(copy, BCF_UN_STR);
        Record split = original;
        split.value = copy;
        split.rid = original.rid;
        split.pos = segment_begin;
        split.alleles.front() = reference;
        split.key = record_key(copy);
        split.reference_block = true;
        result.push_back(std::move(split));
    }
    bcf_destroy(original.value);
    original.value = nullptr;
    return result;
}

// Dense reference expansion is lazy in stream mode: keep the original block
// and clone one one-base record per call instead of allocating a vector for
// the entire END span.  The resulting point records enter the normal locus
// heap and therefore merge with concrete variants from other samples.
std::optional<Record> emit_dense_stream_reference_piece(
    GenotypeStreamCursor& cursor, const bcf_hdr_t* output_header,
    faidx_t* reference_index, bool gatk_annotation_compatibility) {
    if (!cursor.split_block.has_value() || !cursor.dense_reference_block)
        return std::nullopt;
    auto& original = *cursor.split_block;
    const int block_end = record_span_end(output_header, original);
    if (cursor.split_cursor >= block_end) {
        bcf_destroy(original.value);
        original.value = nullptr;
        cursor.split_block.reset();
        cursor.dense_reference_block = false;
        return std::nullopt;
    }
    const int position = cursor.split_cursor++;
    auto* copy = bcf_dup(original.value);
    if (copy == nullptr)
        throw std::runtime_error(
            "RESOURCE_EXHAUSTED: cannot clone dense streaming reference site");
    const auto reference = reference_base(reference_index, output_header,
                                          original.rid, position);
    copy->rid = original.rid;
    copy->pos = position;
    const auto alleles = reference + ",<NON_REF>";
    if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
        bcf_destroy(copy);
        throw std::runtime_error(
            "OUTPUT_CONTRACT_FAILURE: cannot rewrite dense streaming reference site");
    }
    (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
    (void)bcf_update_info_int32(output_header, copy, "AD", nullptr, 0);
    if (!original.dp.empty()) {
        const auto site_dp = !original.min_dp.empty() &&
                                     original.min_dp.front() != bcf_int32_missing &&
                                     original.min_dp.front() != bcf_int32_vector_end
                                 ? original.min_dp.front()
                                 : original.dp.front();
        if (site_dp == bcf_int32_missing || site_dp == bcf_int32_vector_end)
            (void)bcf_update_info_int32(output_header, copy, "DP", nullptr, 0);
        else
            (void)bcf_update_info_int32(output_header, copy, "DP", &site_dp, 1);
    }
    (void)bcf_update_info_int32(output_header, copy, "MIN_DP", nullptr, 0);
    bcf_unpack(copy, BCF_UN_STR);
    Record piece = original;
    piece.value = copy;
    piece.pos = position;
    piece.rid = original.rid;
    piece.alleles = {reference, "<NON_REF>"};
    piece.key = record_key(copy);
    materialize_reference_only(output_header, piece, gatk_annotation_compatibility);
    piece.key = record_key(copy);
    return piece;
}

void begin_lazy_stream_reference_split(GenotypeStreamCursor& cursor,
                                       Record original,
                                       const bcf_hdr_t* output_header,
                                       faidx_t* reference_index,
                                       const std::vector<VariantSpan>& spans,
                                       bool dense_reference_expansion = false) {
    const int block_begin = original.pos;
    const int block_end = record_span_end(output_header, original);
    if (dense_reference_expansion) {
        if (reference_index == nullptr)
            throw std::runtime_error(
                "UNSUPPORTED_PARAMETER: stream-by-locus --include-non-variant-sites requires an indexed -R reference");
        if (cursor.split_block.has_value())
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: nested streaming reference-block split");
        cursor.split_block.emplace(std::move(original));
        cursor.split_span_index = spans.size();
        cursor.split_cursor = block_begin;
        cursor.dense_reference_block = true;
        return;
    }
    if (block_end <= block_begin + 1) {
        cursor.pending.push_back(std::move(original));
        return;
    }
    auto first = std::lower_bound(
        spans.begin(), spans.end(), VariantSpan{original.rid, block_begin, block_begin},
        [](const VariantSpan& left, const VariantSpan& right) {
            if (left.rid != right.rid) return left.rid < right.rid;
            return left.begin < right.begin;
        });
    std::size_t index = static_cast<std::size_t>(first - spans.begin());
    if (index > 0 && spans[index - 1].rid == original.rid &&
        spans[index - 1].end > block_begin)
        --index;
    bool overlaps = false;
    for (std::size_t probe = index; probe < spans.size(); ++probe) {
        if (spans[probe].rid != original.rid || spans[probe].begin >= block_end) break;
        if (spans[probe].end > block_begin) { overlaps = true; break; }
    }
    if (!overlaps) {
        cursor.pending.push_back(std::move(original));
        return;
    }
    // The interior REF base of a split block cannot be recovered from the
    // block's start allele.  Require the same indexed FASTA contract as the
    // aggregate implementation instead of silently inventing N bases.
    if (reference_index == nullptr)
        throw std::runtime_error(
            "UNSUPPORTED_PARAMETER: stream-by-locus reference-block splitting requires an indexed -R reference");
    if (cursor.split_block.has_value())
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: nested streaming reference-block split");
    cursor.split_block.emplace(std::move(original));
    cursor.split_span_index = index;
    cursor.split_cursor = block_begin;
    cursor.dense_reference_block = false;
}

std::optional<Record> emit_lazy_stream_reference_piece(
    GenotypeStreamCursor& cursor, const bcf_hdr_t* output_header,
    faidx_t* reference_index, const std::vector<VariantSpan>& spans,
    bool gatk_annotation_compatibility) {
    if (cursor.dense_reference_block)
        return emit_dense_stream_reference_piece(cursor, output_header,
                                                 reference_index,
                                                 gatk_annotation_compatibility);
    if (!cursor.split_block.has_value()) return std::nullopt;
    auto& original = *cursor.split_block;
    const int block_begin = original.pos;
    const int block_end = record_span_end(output_header, original);
    while (cursor.split_span_index < spans.size()) {
        const auto& span = spans[cursor.split_span_index];
        if (span.rid != original.rid || span.begin >= block_end) break;
        if (span.end <= block_begin || span.end <= cursor.split_cursor) {
            ++cursor.split_span_index;
            continue;
        }
        const int span_begin = std::max(block_begin, span.begin);
        const int span_end = std::min(block_end, span.end);
        if (cursor.split_cursor < span_begin) {
            const int begin = cursor.split_cursor;
            cursor.split_cursor = span_begin;
            auto* copy = bcf_dup(original.value);
            if (!copy) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot clone streaming reference block");
            const auto reference = reference_base(reference_index, output_header, original.rid, begin);
            copy->rid = original.rid;
            copy->pos = begin;
            const auto alleles = reference + ",<NON_REF>";
            if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
                bcf_destroy(copy);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot rewrite streaming block segment");
            }
            const int32_t end_value = span_begin;
            if (span_begin - begin > 1)
                (void)bcf_update_info_int32(output_header, copy, "END", &end_value, 1);
            else
                (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
            bcf_unpack(copy, BCF_UN_STR);
            Record piece = original;
            piece.value = copy;
            piece.pos = begin;
            piece.rid = original.rid;
            piece.alleles.front() = reference;
            piece.key = record_key(copy);
            piece.reference_block = true;
            return piece;
        }
        const int begin = span_begin;
        cursor.split_cursor = std::max(cursor.split_cursor, span_end);
        ++cursor.split_span_index;
        auto* copy = bcf_dup(original.value);
        if (!copy) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot clone streaming reference point");
        const auto reference = reference_base(reference_index, output_header, original.rid, begin);
        copy->rid = original.rid;
        copy->pos = begin;
        const auto alleles = reference + ",<NON_REF>";
        if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
            bcf_destroy(copy);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot rewrite streaming reference point");
        }
        (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
        bcf_unpack(copy, BCF_UN_STR);
        Record piece = original;
        piece.value = copy;
        piece.pos = begin;
        piece.rid = original.rid;
        piece.alleles.front() = reference;
        piece.key = record_key(copy);
        piece.reference_block = true;
        return piece;
    }
    if (cursor.split_cursor < block_end) {
        const int begin = cursor.split_cursor;
        cursor.split_cursor = block_end;
        auto* copy = bcf_dup(original.value);
        if (!copy) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot clone trailing streaming block");
        const auto reference = reference_base(reference_index, output_header, original.rid, begin);
        copy->rid = original.rid;
        copy->pos = begin;
        const auto alleles = reference + ",<NON_REF>";
        if (bcf_update_alleles_str(output_header, copy, alleles.c_str()) != 0) {
            bcf_destroy(copy);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot rewrite trailing streaming block");
        }
        const int32_t end_value = block_end;
        if (block_end - begin > 1)
            (void)bcf_update_info_int32(output_header, copy, "END", &end_value, 1);
        else
            (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
        bcf_unpack(copy, BCF_UN_STR);
        Record piece = original;
        piece.value = copy;
        piece.pos = begin;
        piece.rid = original.rid;
        piece.alleles.front() = reference;
        piece.key = record_key(copy);
        piece.reference_block = true;
        return piece;
    }
    bcf_destroy(original.value);
    original.value = nullptr;
    cursor.split_block.reset();
    cursor.dense_reference_block = false;
    return std::nullopt;
}

// This is the stream-mode header contract shared by all input cursors.  It is
// intentionally the same field set as the aggregate path so a caller can
// switch modes without changing downstream VCF schemas.
void add_genotype_output_header_fields(bcf_hdr_t* output_header,
                                       const Options& options) {
    bcf_hdr_append(output_header,
                   "##fastgatk_genotype_gvcfs_status=reference-block-materialization");
    ensure_gatk_low_qual_filter_declaration(output_header, options);
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "GQ") < 0)
        bcf_hdr_append(output_header,
                       "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "RGQ") < 0)
        bcf_hdr_append(output_header,
                       "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">");
    if (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
        !options.gatk_annotation_compatibility) {
        if (!has_header_field(output_header, BCF_HL_FMT, "GP"))
            bcf_hdr_append(output_header,
                           "##FORMAT=<ID=GP,Number=G,Type=Float,Description=Genotype posterior in Phred Scale>");
        if (!has_header_field(output_header, BCF_HL_FMT, "PG"))
            bcf_hdr_append(output_header,
                           "##FORMAT=<ID=PG,Number=G,Type=Float,Description=Genotype priors in Phred Scale>");
    }
    // AC/AN/AF/DP/MLEAC/MLEAF/... are declared with GATK's own text in the
    // GATK-compatibility profile (see kGatkMleacInfoLine() and
    // kGatkStandardDpInfoLine() above); the diagnostic profile keeps the
    // shorter descriptions it has always written.
    const std::string dp_info_line = options.gatk_annotation_compatibility
        ? std::string(kGatkStandardDpInfoLine)
        : std::string("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">");
    const std::string mleac_info_line = options.gatk_annotation_compatibility
        ? std::string(kGatkMleacInfoLine)
        : std::string("##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=Maximum likelihood allele count>");
    const std::string mleaf_info_line = options.gatk_annotation_compatibility
        ? std::string(kGatkMleafInfoLine)
        : std::string("##INFO=<ID=MLEAF,Number=A,Type=Float,Description=Maximum likelihood allele frequency>");
    const std::vector<std::string> info_headers{
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">",
        dp_info_line,
        mleac_info_line,
        mleaf_info_line,
        "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">",
        "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">",
        "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">",
        "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">",
        "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">",
        "##INFO=<ID=RCQ,Number=1,Type=Float,Description=Joint cross-sample reference-confidence quality>",
        "##INFO=<ID=RCP,Number=1,Type=Float,Description=Joint probability that all samples are hom-reference>"};
    for (const auto& line : info_headers) {
        const auto marker = line.find("ID=");
        const auto end = line.find(',', marker == std::string::npos ? 0 : marker + 3);
        const auto id = marker == std::string::npos ? std::string{} :
            line.substr(marker + 3, end == std::string::npos ? std::string::npos : end - marker - 3);
        if (!id.empty() && bcf_hdr_id2int(output_header, BCF_DT_ID, id.c_str()) < 0)
            bcf_hdr_append(output_header, line.c_str());
    }
    if (options.annotate_with_num_discovered_alleles &&
        bcf_hdr_id2int(output_header, BCF_DT_ID, "NDA") < 0)
        bcf_hdr_append(output_header, options.gatk_annotation_compatibility
            ? kGatkNdaInfoLine
            : "##INFO=<ID=NDA,Number=1,Type=Integer,Description=Number of discovered alternate alleles>");
    if (bcf_hdr_sync(output_header) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync streaming VCF header");
}

void write_genotype_stream_manifest(const Options& options,
                                    const fastgatk::runtime::ResourceSnapshot& resources,
                                    const std::string& manifest_path,
                                    const std::string& index_path,
                                    std::uint64_t input_files,
                                    std::uint64_t input_records,
                                    std::uint64_t output_records,
                                    std::uint64_t skipped_reference_blocks,
                                    std::uint64_t interval_skipped,
                                    std::uint64_t indexed_inputs,
                                    std::uint64_t indexed_interval_queries,
                                    std::uint64_t stream_span_probe_records,
                                    int sample_count) {
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("cannot write manifest: " + manifest_path);
    manifest << "{\"schema_version\":1,\"tool\":\"GenotypeGVCFs\",\"implementation\":\"fastgatk-genotype-gvcf\",\"status\":\"contract-compatible\"," 
             << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
             << "\"compatibility\":{\"gatk_parameter_aliases\":true,\"reference_block_materialization\":true,\"nonref_pl_projection\":true,\"include_non_variant_sites\":"
             << (options.include_non_variant_sites ? "true" : "false")
             << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
             << ",\"interval_exclusion\":true"
             << ",\"output_interval_starts_in\":true"
             << ",\"format_compaction\":true,\"multi_allelic_materialization\":true,\"allele_union\":true,\"joint_genotype_from_pl\":true,\"kokkos_genotype_pl\":true,\"arbitrary_ploidy_gt_gq\":true,\"kokkos_allele_counts\":true,\"site_annotations\":true,\"excess_het_exact_diploid\":true,\"inbreeding_coeff_exact_diploid\":true,\"multi_sample_joint_merge\":true,\"cross_sample_reference_block_split\":true,\"cross_sample_reference_confidence\":true,\"genotype_prior_calculator_assumingHW\":true,\"allele_frequency_calculator_em\":true,\"cohort_af_calculator\":true,\"posterior_qual_opt_in\":true,\"output_allele_subset\":true,\"standard_confidence_allele_pruning\":true,\"max_alternate_alleles_likelihood_subset\":true,\"spanning_deletion_orphan_cleanup\":true,\"genotype_assignment_method\":true,\"posterior_genotype_assignment\":false,\"gatk_annotation_compatibility\":"
             << (options.gatk_annotation_compatibility ? "true" : "false")
             << ",\"num_discovered_alleles_annotation\":true"
             << ",\"stream_by_locus\":true,\"bounded_k_way_merge\":true,\"three_stage_pipeline\":true,\"vcf_index\":"
             << (index_path.empty() ? "false" : "true") << ",\"indexed_interval_traversal\":"
             << (indexed_inputs > 0 ? "true" : "false") << ",\"genomicsdb_bridge\":"
             << (genomicsdb_bridge_used ? "true" : "false") << ",\"bit_identical_to_gatk\":false},"
             << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
             << "\",\"kind\":\"vcf\",\"complete\":true}"
             << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":true}")
             << "],\"telemetry\":{\"resources\":" << resources.to_json()
             << ",\"input_files\":" << input_files << ",\"input_records\":" << input_records
             << ",\"output_records\":" << output_records << ",\"skipped_reference_blocks\":"
             << skipped_reference_blocks << ",\"intervals\":" << options.regions.size()
             << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
             << ",\"exclude_intervals\":" << options.excluded_regions.size()
             << ",\"only_output_calls_starting_in_intervals\":"
             << (options.only_output_calls_starting_in_intervals ? "true" : "false")
             << ",\"output_interval_skipped\":" << options.output_interval_skipped
             << ",\"include_non_variant_sites\":"
             << (options.include_non_variant_sites ? "true" : "false")
             << ",\"interval_skipped\":" << interval_skipped << ",\"indexed_inputs\":"
             << indexed_inputs << ",\"indexed_interval_queries\":" << indexed_interval_queries
             << ",\"stream_span_probe_records\":" << stream_span_probe_records
             << ",\"streamed_loci\":" << options.streamed_loci
             << ",\"streamed_peak_host_bytes\":" << options.streamed_peak_host_bytes
             << ",\"stream_max_inflight_records\":" << options.stream_max_inflight_records
             << ",\"genomicsdb_bridge\":" << (genomicsdb_bridge_used ? "true" : "false")
             << ",\"sample_count\":" << std::max(0, sample_count)
             << ",\"standard_confidence_for_calling\":" << options.standard_confidence_for_calling
             << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.calls
             << ",\"genotype_kernel_prepare_seconds\":" << genotype_kernel_telemetry.prepare_seconds
             << ",\"genotype_kernel_seconds\":" << genotype_kernel_telemetry.execute_seconds
             << ",\"genotype_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.execution_space) << "\""
             << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
             << ",\"pl_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.pl_remap_prepare_seconds
             << ",\"pl_remap_kernel_seconds\":" << genotype_kernel_telemetry.pl_remap_execute_seconds
             << ",\"pl_remap_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.pl_remap_execution_space) << "\""
             << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
             << ",\"allele_field_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
             << ",\"allele_field_remap_kernel_seconds\":" << genotype_kernel_telemetry.allele_field_remap_execute_seconds
             << ",\"allele_field_remap_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space) << "\""
             << ",\"allele_count_kernel_calls\":" << genotype_kernel_telemetry.allele_count_calls
             << ",\"allele_count_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_count_prepare_seconds
             << ",\"allele_count_kernel_seconds\":" << genotype_kernel_telemetry.allele_count_execute_seconds
             << ",\"allele_count_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.allele_count_execution_space) << "\""
             << ",\"posterior_kernel_calls\":" << genotype_kernel_telemetry.posterior_calls
             << ",\"posterior_kernel_seconds\":" << genotype_kernel_telemetry.posterior_execute_seconds
             << ",\"posterior_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.posterior_execution_space) << "\""
             << ",\"genotype_prior_kernel_calls\":" << genotype_kernel_telemetry.genotype_prior_calls
             << ",\"genotype_prior_kernel_prepare_seconds\":" << genotype_kernel_telemetry.genotype_prior_prepare_seconds
             << ",\"genotype_prior_kernel_seconds\":" << genotype_kernel_telemetry.genotype_prior_execute_seconds
             << ",\"genotype_prior_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.genotype_prior_execution_space) << "\""
             << ",\"cohort_af_kernel_calls\":" << genotype_kernel_telemetry.cohort_calls
             << ",\"cohort_af_kernel_prepare_seconds\":" << genotype_kernel_telemetry.cohort_prepare_seconds
             << ",\"cohort_af_kernel_seconds\":" << genotype_kernel_telemetry.cohort_execute_seconds
             << ",\"cohort_af_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.cohort_execution_space) << "\""
             << ",\"cohort_af_samples\":" << genotype_kernel_telemetry.cohort_samples
             << ",\"cohort_af_converged\":" << genotype_kernel_telemetry.cohort_converged
             << ",\"output_allele_pruning_calls\":" << genotype_kernel_telemetry.output_allele_pruning_calls
             << ",\"output_alleles_pruned\":" << genotype_kernel_telemetry.output_alleles_pruned
             << ",\"orphan_spanning_deletion_loci\":" << options.orphan_spanning_deletion_loci
             << ",\"max_alternate_alleles\":" << options.max_alternate_alleles
             << ",\"max_alt_pruning_calls\":" << genotype_kernel_telemetry.max_alt_pruning_calls
             << ",\"max_alt_alleles_pruned\":" << genotype_kernel_telemetry.max_alt_alleles_pruned
             << ",\"max_alt_score_kernel_calls\":" << genotype_kernel_telemetry.max_alt_score_kernel_calls
             << ",\"max_alt_score_kernel_prepare_seconds\":" << genotype_kernel_telemetry.max_alt_score_prepare_seconds
             << ",\"max_alt_score_kernel_seconds\":" << genotype_kernel_telemetry.max_alt_score_execute_seconds
             << ",\"max_alt_score_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.max_alt_score_execution_space) << "\""
             << ",\"annotate_with_num_discovered_alleles\":" << (options.annotate_with_num_discovered_alleles ? "true" : "false")
             << ",\"num_discovered_alleles_calls\":" << genotype_kernel_telemetry.num_discovered_alleles_calls
             << ",\"cross_sample_reference_kernel_calls\":" << genotype_kernel_telemetry.cross_sample_reference_calls
             << ",\"cross_sample_reference_kernel_prepare_seconds\":" << genotype_kernel_telemetry.cross_sample_reference_prepare_seconds
             << ",\"cross_sample_reference_kernel_seconds\":" << genotype_kernel_telemetry.cross_sample_reference_execute_seconds
             << ",\"cross_sample_reference_execution_space\":\"" << json_escape(genotype_kernel_telemetry.cross_sample_reference_execution_space) << "\""
             << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\",\"three_stage_pipeline\":true"
             << ",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
             << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
             << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
             << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
             << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
             << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
             << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
             << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
             << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
             << ",\"pipeline_stage_capacity_bytes\":" << options.pipeline_stage_capacity_bytes << "}}\n";
}

int run_streaming_genotype_gvcf(Options& options,
                                const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto input_paths = expand_genomicsdb_inputs(options.inputs, options.regions);
    bcf_hdr_t* output_header = nullptr;
    faidx_t* reference_index = nullptr;
    htsFile* output_file = nullptr;
    std::vector<Region> regions;
    std::vector<Region> excluded_regions;
    fastgatk::io::IntervalFileStats interval_file_stats;
    std::vector<std::string> sample_names;
    std::map<std::string, int> output_sample_indices;
    std::uint64_t input_records = 0;
    std::uint64_t span_probe_records = 0;
    std::uint64_t interval_skipped = 0;
    std::uint64_t skipped_reference_blocks = 0;
    std::uint64_t indexed_inputs = 0;
    std::uint64_t indexed_interval_queries = 0;
    int expected_sample_count = -1;
    std::vector<VariantSpan> variant_spans;

    try {
        if (!options.reference.empty()) reference_index = fai_load(options.reference.c_str());

        // First pass: merge dictionaries/sample columns and create the exact
        // output header before any cursor starts decoding records.  This pass
        // reads headers only and is independent of source record order.
        for (const auto& input_path : input_paths) {
            htsFile* input = bcf_open(input_path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) {
                bcf_close(input);
                throw std::runtime_error("BAD_INPUT: cannot read VCF header: " + input_path);
            }
            if (options.genotype_assignment_method == "USE_POSTERIORS_ANNOTATION" &&
                bcf_hdr_id2int(header, BCF_DT_ID, "PP") < 0) {
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error(
                    "BAD_INPUT: USE_POSTERIORS_ANNOTATION requires FORMAT/PP in every input header");
            }
            if (!output_header) {
                output_header = bcf_hdr_dup(header);
                if (!output_header) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate VCF header");
                }
                add_genotype_output_header_fields(output_header, options);
            }
            merge_contig_dictionary(output_header, header);
            std::vector<int> input_samples;
            input_samples.reserve(static_cast<std::size_t>(header->n[BCF_DT_SAMPLE]));
            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (name == nullptr || *name == '\0') {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("BAD_INPUT: GVCF header contains an invalid sample name");
                }
                auto found = output_sample_indices.find(name);
                if (found == output_sample_indices.end()) {
                    const int next = static_cast<int>(output_sample_indices.size());
                    output_sample_indices.emplace(name, next);
                    sample_names.emplace_back(name);
                    input_samples.push_back(next);
                } else {
                    input_samples.push_back(found->second);
                }
            }
            expected_sample_count = static_cast<int>(sample_names.size());
            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (bcf_hdr_id2int(output_header, BCF_DT_SAMPLE, name) < 0 &&
                    bcf_hdr_add_sample(output_header, name) != 0) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge VCF sample header");
                }
            }
            bcf_hdr_destroy(header);
            bcf_close(input);
        }
        if (!output_header) throw std::runtime_error("BAD_INPUT: no VCF/GVCF inputs");
        if (bcf_hdr_sync(output_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync merged VCF header");

        if (!input_paths.empty()) {
            htsFile* input = bcf_open(input_paths.front().c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen VCF/GVCF: " + input_paths.front());
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) {
                bcf_close(input);
                throw std::runtime_error("BAD_INPUT: cannot reread VCF header: " + input_paths.front());
            }
            bool first_selector = true;
            for (const auto& text : options.regions) {
                append_region_selector_with_rule(text, header, regions, interval_file_stats,
                                                  options.interval_set_rule, first_selector);
                first_selector = false;
            }
            for (const auto& text : options.excluded_regions)
                append_region_selector(text, header, excluded_regions, interval_file_stats);
            bcf_hdr_destroy(header);
            bcf_close(input);
            normalize_regions(regions);
            normalize_regions(excluded_regions);
        }

        // Probe concrete variant spans once.  Reference blocks are then split
        // lazily by each cursor only when they overlap a concrete variant from
        // another sample; no FORMAT/PL arrays are retained by this index.
        for (const auto& input_path : input_paths) {
            htsFile* input = bcf_open(input_path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            bcf_hdr_t* header = bcf_hdr_read(input);
            bcf1_t* record = bcf_init();
            if (!header || !record) {
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot initialize streaming span probe");
            }
            int probe_status = 0;
            while ((probe_status = bcf_read(input, header, record)) == 0) {
                ++span_probe_records;
                bcf_unpack(record, BCF_UN_ALL);
                if (!in_regions(header, record, regions, excluded_regions)) continue;
                bool concrete = false;
                for (int allele = 1; allele < record->n_allele; ++allele) {
                    const auto* alternate = record->d.allele[allele];
                    if (alternate == nullptr || std::strcmp(alternate, "<NON_REF>") == 0)
                        continue;
                    concrete = true;
                }
                if (!concrete) continue;
                const auto* name = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
                const int output_rid = name == nullptr ? -1 : bcf_hdr_name2id(output_header, name);
                if (output_rid < 0) continue;
                int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
                if (end <= record->pos) end = record->pos + 1;
                variant_spans.push_back(VariantSpan{output_rid, static_cast<int>(record->pos), end});
            }
            // A record htslib could not parse (-2) must not read as end of
            // input (-1); the probe would otherwise index a truncated file.
            if (probe_status < -1) {
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error(
                    "BAD_INPUT: streaming GenotypeGVCFs VCF record parse failed");
            }
            bcf_destroy(record);
            bcf_hdr_destroy(header);
            bcf_close(input);
        }
        std::sort(variant_spans.begin(), variant_spans.end(), [](const VariantSpan& left,
                                                                  const VariantSpan& right) {
            if (left.rid != right.rid) return left.rid < right.rid;
            if (left.begin != right.begin) return left.begin < right.begin;
            return left.end < right.end;
        });
        std::vector<VariantSpan> merged_spans;
        for (const auto& span : variant_spans) {
            if (merged_spans.empty() || span.rid != merged_spans.back().rid ||
                span.begin > merged_spans.back().end)
                merged_spans.push_back(span);
            else
                merged_spans.back().end = std::max(merged_spans.back().end, span.end);
        }
        variant_spans.swap(merged_spans);

        std::vector<std::unique_ptr<GenotypeStreamCursor>> cursors;
        cursors.reserve(input_paths.size());
        for (const auto& input_path : input_paths) {
            auto cursor = std::make_unique<GenotypeStreamCursor>();
            cursor->input = bcf_open(input_path.c_str(), "r");
            if (!cursor->input)
                throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            cursor->header = bcf_hdr_read(cursor->input);
            cursor->raw = bcf_init();
            if (!cursor->header || !cursor->raw)
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot initialize streaming cursor");
            if (!regions.empty()) {
                cursor->traversal_index = load_genotype_traversal_index(input_path);
                if (cursor->traversal_index != nullptr) ++indexed_inputs;
            }
            for (int sample = 0; sample < cursor->header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(cursor->header, BCF_DT_SAMPLE, sample);
                cursor->output_samples.push_back(output_sample_indices.at(name));
            }
            cursors.push_back(std::move(cursor));
        }

        // Decode/materialize one source record.  The lambda deliberately
        // mirrors the aggregate path's allele/FORMAT projection, but returns
        // immediately after filling a cursor's small pending deque.
        const auto check_cursor_order = [](GenotypeStreamCursor& cursor,
                                           const Record& record) {
            if (cursor.have_last_coordinate &&
                (record.rid < cursor.last_rid ||
                 (record.rid == cursor.last_rid && record.pos < cursor.last_pos)))
                throw std::runtime_error(
                    "BAD_INPUT: --stream-by-locus requires each input VCF/GVCF to be coordinate sorted");
            cursor.have_last_coordinate = true;
            cursor.last_rid = record.rid;
            cursor.last_pos = record.pos;
        };
        const auto fill_cursor = [&](GenotypeStreamCursor& cursor) -> bool {
            if (!cursor.pending.empty()) return true;
            if (auto piece = emit_lazy_stream_reference_piece(cursor, output_header,
                                                               reference_index, variant_spans,
                                                               options.gatk_annotation_compatibility)) {
                const auto piece_bytes = genotype_record_bytes(*piece);
                const auto block_bytes = cursor.split_block.has_value()
                    ? genotype_record_bytes(*cursor.split_block) : 0;
                options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                    options.streamed_peak_host_bytes, piece_bytes + block_bytes);
                cursor.pending.push_back(std::move(*piece));
                check_cursor_order(cursor, cursor.pending.back());
                return true;
            }
            while (read_next_genotype_cursor_record(cursor, regions,
                                                    indexed_interval_queries)) {
                ++input_records;
                bcf_unpack(cursor.raw, BCF_UN_ALL);
                if (!in_regions(cursor.header, cursor.raw, regions, excluded_regions)) {
                    ++interval_skipped;
                    continue;
                }
                std::vector<int> selected;
                int non_ref_index = -1;
                for (int allele = 1; allele < cursor.raw->n_allele; ++allele) {
                    if (cursor.raw->d.allele[allele] == nullptr) continue;
                    if (std::strcmp(cursor.raw->d.allele[allele], "<NON_REF>") == 0)
                        non_ref_index = allele;
                    else
                        selected.push_back(allele);
                }
                const bool reference_only = selected.empty() && non_ref_index >= 0;
                if (selected.empty() && non_ref_index < 0) {
                    ++skipped_reference_blocks;
                    continue;
                }
                if (non_ref_index >= 0) selected.push_back(non_ref_index);
                auto* materialized = bcf_dup(cursor.raw);
                if (!materialized)
                    throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_dup failed in streaming cursor");
                const int original_alleles = cursor.raw->n_allele;
                Record staged;
                extract_materialized_fields(cursor.header, cursor.raw, selected,
                                            original_alleles, cursor.output_samples, staged);
                staged.reference_block = reference_only;
                std::string alleles = cursor.raw->d.allele[0] == nullptr ? "N" : cursor.raw->d.allele[0];
                for (const auto allele : selected) {
                    alleles.push_back(',');
                    alleles += cursor.raw->d.allele[allele] == nullptr ? "N" :
                        cursor.raw->d.allele[allele];
                }
                if (bcf_update_alleles_str(output_header, materialized, alleles.c_str()) != 0) {
                    bcf_destroy(materialized);
                    throw std::runtime_error(
                        "OUTPUT_CONTRACT_FAILURE: cannot materialize streaming VCF alleles");
                }
                int32_t* info_ad = nullptr;
                int info_ad_count = 0;
                const auto info_ad_length = bcf_get_info_int32(output_header, materialized,
                                                               "AD", &info_ad, &info_ad_count);
                if (info_ad_length >= original_alleles && info_ad_count >= original_alleles) {
                    std::vector<int32_t> compact_info_ad(static_cast<std::size_t>(selected.size()) + 1,
                                                         bcf_int32_missing);
                    compact_info_ad[0] = info_ad[0];
                    for (std::size_t index = 0; index < selected.size(); ++index)
                        compact_info_ad[index + 1] = info_ad[selected[index]];
                    if (bcf_update_info_int32(output_header, materialized, "AD",
                                              compact_info_ad.data(),
                                              static_cast<int>(compact_info_ad.size())) != 0) {
                        free(info_ad);
                        bcf_destroy(materialized);
                        throw std::runtime_error(
                            "OUTPUT_CONTRACT_FAILURE: cannot compact streaming AD INFO field");
                    }
                }
                free(info_ad);
                bcf_unpack(materialized, BCF_UN_STR);
                const auto* contig_name = cursor.raw->rid >= 0
                    ? bcf_hdr_id2name(cursor.header, cursor.raw->rid) : nullptr;
                const int output_rid = contig_name == nullptr
                    ? -1 : bcf_hdr_name2id(output_header, contig_name);
                if (output_rid < 0) {
                    bcf_destroy(materialized);
                    throw std::runtime_error(
                        "OUTPUT_CONTRACT_FAILURE: cannot map streaming VCF contig");
                }
                materialized->rid = output_rid;
                staged.value = materialized;
                staged.rid = output_rid;
                staged.pos = static_cast<int>(materialized->pos);
                staged.alleles.reserve(selected.size() + 1);
                staged.alleles.emplace_back(cursor.raw->d.allele[0] == nullptr ? "N" :
                                            cursor.raw->d.allele[0]);
                for (const auto allele : selected)
                    staged.alleles.emplace_back(cursor.raw->d.allele[allele] == nullptr ? "N" :
                                                cursor.raw->d.allele[allele]);
                staged.key = record_key(materialized);
                if (reference_only) {
                    begin_lazy_stream_reference_split(cursor, std::move(staged), output_header,
                                                      reference_index, variant_spans,
                                                      options.include_non_variant_sites);
                    if (cursor.split_block.has_value())
                        options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                            options.streamed_peak_host_bytes,
                            genotype_record_bytes(*cursor.split_block));
                } else {
                    cursor.pending.push_back(std::move(staged));
                }
                if (cursor.pending.empty()) {
                    if (auto piece = emit_lazy_stream_reference_piece(cursor, output_header,
                                                                       reference_index, variant_spans,
                                                                       options.gatk_annotation_compatibility)) {
                        options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                            options.streamed_peak_host_bytes,
                            genotype_record_bytes(*piece) +
                                (cursor.split_block.has_value()
                                     ? genotype_record_bytes(*cursor.split_block) : 0));
                        cursor.pending.push_back(std::move(*piece));
                    }
                }
                if (!cursor.pending.empty()) {
                    check_cursor_order(cursor, cursor.pending.back());
                    return true;
                }
            }
            return false;
        };

        std::priority_queue<GenotypeStreamHeapItem,
                            std::vector<GenotypeStreamHeapItem>,
                            GenotypeStreamHeapCompare> heap;
        std::uint64_t queued_bytes = 0;
        const auto push_cursor = [&](const std::size_t source) {
            auto& cursor = *cursors[source];
            if (!fill_cursor(cursor)) return false;
            GenotypeStreamHeapItem item;
            item.record = std::move(cursor.pending.front());
            cursor.pending.pop_front();
            queued_bytes += genotype_record_bytes(item.record);
            heap.push(std::move(item));
            options.stream_max_inflight_records = std::max<std::uint64_t>(
                options.stream_max_inflight_records, heap.size());
            options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                options.streamed_peak_host_bytes, queued_bytes);
            return true;
        };
        for (std::size_t source = 0; source < cursors.size(); ++source)
            (void)push_cursor(source);

        output_file = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
        if (!output_file) throw std::runtime_error("cannot open output VCF: " + options.output);
        auto* output = output_file;
        if (options.gatk_annotation_compatibility) {
            kstring_t formatted_header{0, 0, nullptr};
            if (bcf_hdr_format(output_header, 0, &formatted_header) < 0 ||
                write_vcf_text_line(output, gatk_compatible_header_text(
                    std::string(formatted_header.s == nullptr ? "" : formatted_header.s,
                                formatted_header.l),
                    input_paths.empty()
                        ? std::string{}
                        : input_pass_filter_line(input_paths.front()))) != 0) {
                free(formatted_header.s);
                bcf_close(output);
                output_file = nullptr;
                throw std::runtime_error("cannot write GATK-compatible streaming VCF header");
            }
            free(formatted_header.s);
        } else if (bcf_hdr_write(output, output_header) != 0) {
            bcf_close(output);
            output_file = nullptr;
            throw std::runtime_error("cannot write streaming VCF header");
        }

        GatkJavaRandom gatk_random;
        std::uint64_t output_record_count = 0;
        // GenotypingEngine's per-traversal emitted-deletion state
        // (GenotypingEngine.java:52).  It is read and written only by the
        // single compute worker below, which visits loci in traversal order.
        EmittedDeletions upstream_deletions;
        const auto stage_capacity = std::max<std::uint64_t>(
            64ULL * 1024ULL * 1024ULL, resources.safe_memory_budget_bytes());
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<
            GenotypeDecoded, GenotypeComputed, GenotypeEncoded>;
        Pipeline pipeline(
            Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
            [&]() -> std::optional<GenotypeDecoded> {
                while (!heap.empty()) {
                    auto first = std::move(const_cast<GenotypeStreamHeapItem&>(heap.top()));
                    heap.pop();
                    queued_bytes -= genotype_record_bytes(first.record);
                    const auto locus_key = first.record.key;
                    std::vector<Record> group;
                    group.push_back(std::move(first.record));
                    (void)push_cursor(first.source);
                    while (!heap.empty() && heap.top().record.key == locus_key) {
                        auto next = std::move(const_cast<GenotypeStreamHeapItem&>(heap.top()));
                        heap.pop();
                        queued_bytes -= genotype_record_bytes(next.record);
                        group.push_back(std::move(next.record));
                        (void)push_cursor(next.source);
                    }
                    options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                        options.streamed_peak_host_bytes, queued_bytes);
                    options.stream_max_inflight_records = std::max<std::uint64_t>(
                        options.stream_max_inflight_records, heap.size() + group.size());
                    for (const auto& record : group)
                        options.streamed_peak_host_bytes = std::max<std::uint64_t>(
                            options.streamed_peak_host_bytes, queued_bytes + genotype_record_bytes(record));

                    std::vector<std::string> union_alleles{group.front().alleles.front()};
                    bool has_non_ref = false;
                    // A '*' is left in the merged allele union either way: GATK
                    // decides whether it is a spurious spanning deletion inside
                    // calculateOutputAlleleSubset() (:314), i.e. after the union
                    // exists and after the earlier loci have reported what they
                    // emitted.  See apply_gatk_output_allele_subset().
                    for (const auto& record : group) {
                        for (std::size_t allele = 1; allele < record.alleles.size(); ++allele) {
                            const auto& name = record.alleles[allele];
                            if (name == "<NON_REF>") { has_non_ref = true; continue; }
                            if (std::find(union_alleles.begin() + 1, union_alleles.end(), name) == union_alleles.end())
                                union_alleles.push_back(name);
                        }
                    }
                    if (has_non_ref) union_alleles.push_back("<NON_REF>");
                    std::vector<std::string> concrete;
                    for (const auto& allele : union_alleles)
                        if (allele != "<NON_REF>") concrete.push_back(allele);
                    if (concrete.size() < 2) {
                        if (!options.include_non_variant_sites) {
                            for (auto& record : group) bcf_destroy(record.value);
                            continue;
                        }
                        // Dense mode emits one REF-only record for every
                        // coordinate.  Each source is already a one-base
                        // materialized point, so compact it and merge
                        // disjoint sample columns exactly as the aggregate
                        // path does for a reference-only locus.
                        for (auto& record : group)
                            if (record.reference_block)
                                materialize_reference_only(output_header, record,
                                                           options.gatk_annotation_compatibility);
                        std::vector<Record*> merged{&group.front()};
                        std::map<int, bool> seen_samples;
                        for (const auto sample : group.front().output_samples)
                            seen_samples[sample] = true;
                        for (std::size_t index = 1; index < group.size(); ++index) {
                            bool overlaps = false;
                            for (const auto sample : group[index].output_samples)
                                if (seen_samples.contains(sample)) { overlaps = true; break; }
                            if (overlaps) {
                                bcf_destroy(group[index].value);
                                group[index].value = nullptr;
                                continue;
                            }
                            for (const auto sample : group[index].output_samples)
                                seen_samples[sample] = true;
                            merged.push_back(&group[index]);
                        }
                        merge_sample_fields(output_header, merged);
                        for (std::size_t index = 1; index < group.size(); ++index)
                            if (group[index].value != nullptr) {
                                bcf_destroy(group[index].value);
                                group[index].value = nullptr;
                            }
                        ++options.streamed_loci;
                        GenotypeDecoded decoded;
                        decoded.record = std::move(group.front());
                        return decoded;
                    }
                    for (auto& record : group) {
                        remap_record_to_allele_union(output_header, record, union_alleles);
                        // The merged genotype is the source call expressed in the
                        // union allele space; keep it before the PL-derived
                        // assignment below consumes it.
                        record.source_gt = record.gt;
                        record.source_alleles = record.alleles;
                        remove_non_ref_allele(output_header, record, concrete,
                                              options.genotype_assignment_method == "BEST_MATCH_TO_ORIGINAL");
                        if (options.genotype_assignment_method != "DO_NOT_ASSIGN_GENOTYPES" &&
                            options.genotype_assignment_method != "BEST_MATCH_TO_ORIGINAL" &&
                            options.genotype_assignment_method != "USE_POSTERIORS_ANNOTATION")
                            derive_gt_gq_from_pl(record);
                    }
                    std::vector<Record*> merged{&group.front()};
                    std::map<int, bool> seen_samples;
                    for (const auto sample : group.front().output_samples) seen_samples[sample] = true;
                    for (std::size_t index = 1; index < group.size(); ++index) {
                        bool overlaps = false;
                        for (const auto sample : group[index].output_samples)
                            if (seen_samples.contains(sample)) { overlaps = true; break; }
                        if (overlaps) { bcf_destroy(group[index].value); group[index].value = nullptr; continue; }
                        for (const auto sample : group[index].output_samples) seen_samples[sample] = true;
                        merged.push_back(&group[index]);
                    }
                    merge_sample_fields(output_header, merged);
                    for (std::size_t index = 1; index < group.size(); ++index)
                        if (group[index].value != nullptr) { bcf_destroy(group[index].value); group[index].value = nullptr; }
                    ++options.streamed_loci;
                    GenotypeDecoded decoded;
                    decoded.record = std::move(group.front());
                    return decoded;
                }
                return std::nullopt;
            },
            [&](GenotypeDecoded decoded) -> std::optional<GenotypeComputed> {
                auto record = std::move(decoded.record);
                if (gatk_skips_regenotyping(output_header, record, options)) {
                    // GenotypeGVCFsEngine.java:160 sends a variant record whose
                    // merged INFO/DP is not positive to the passthrough arm at
                    // :174, and :181 then rejects the locus because the depth
                    // test fails again; the default traversal therefore emits no
                    // record for it, even with a confident call in the
                    // genotypes.  No deletion is recorded either, because GATK
                    // never reaches recordDeletions() on this path.
                    bcf_destroy(record.value);
                    record.value = nullptr;
                    return std::nullopt;
                }
                std::vector<double> posterior_priors;
                annotate_num_discovered_alleles(output_header, record, options);
                (void)apply_gatk_max_alternate_alleles(output_header, record, options);
                update_cohort_af_annotations(output_header, record, options);
                if (!apply_gatk_output_allele_subset(output_header, record, options,
                                                     upstream_deletions))
                    return std::nullopt;
                if (record.finalized_monomorphic_ref) {
                    // The REF-only call is complete: GATK emits it directly from
                    // regenotypeVC, so no site or standard annotation runs.
                    apply_gatk_annotation_compatibility(output_header, record, options);
                    GenotypeComputed computed;
                    computed.record = std::move(record);
                    return computed;
                }
                if (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
                    !options.gatk_annotation_compatibility &&
                    record.allele_count >= 2 && record.ploidy > 0)
                    posterior_priors = build_log10_genotype_priors(
                        record.alleles, record.ploidy, options.snp_heterozygosity,
                        options.indel_heterozygosity, options.use_genotype_priors);
                apply_genotype_assignment(output_header, record, options, posterior_priors);
                normalize_haplotype_caller_phasing(output_header, record);
                update_site_annotations(output_header, record);
                update_total_depth_annotation(output_header, record, options);
                update_cross_sample_reference_confidence(output_header, record, options);
                update_site_quality_from_posteriors(record, options);
                update_gatk_standard_annotations(output_header, record, gatk_random);
                update_excess_het_annotation(output_header, record);
                update_inbreeding_coeff_annotation(output_header, record);
                apply_gatk_annotation_compatibility(output_header, record, options);
                GenotypeComputed computed;
                computed.record = std::move(record);
                return computed;
            },
            [&](GenotypeComputed computed) -> std::optional<GenotypeEncoded> {
                GenotypeEncoded encoded;
                encoded.record = std::move(computed.record);
                if (options.gatk_annotation_compatibility) {
                    kstring_t formatted{0, 0, nullptr};
                    if (vcf_format1(output_header, encoded.record.value, &formatted) != 0) {
                        free(formatted.s);
                        throw std::runtime_error("cannot format GATK-compatible streaming VCF record");
                    }
                    encoded.text = gatk_compatible_record_text(
                        std::string(formatted.s == nullptr ? "" : formatted.s, formatted.l));
                    free(formatted.s);
                }
                return encoded;
            },
            [&](GenotypeEncoded encoded) {
                if (options.only_output_calls_starting_in_intervals &&
                    !starts_in_regions(output_header, encoded.record.value, regions)) {
                    ++options.output_interval_skipped;
                    bcf_destroy(encoded.record.value);
                    encoded.record.value = nullptr;
                    return;
                }
                if (options.gatk_annotation_compatibility) {
                    if (write_vcf_text_line(output, encoded.text) != 0) {
                        bcf_destroy(encoded.record.value); encoded.record.value = nullptr;
                        throw std::runtime_error("cannot write GATK-compatible streaming VCF record");
                    }
                } else if (bcf_write(output, output_header, encoded.record.value) != 0) {
                    bcf_destroy(encoded.record.value); encoded.record.value = nullptr;
                    throw std::runtime_error("cannot write streaming VCF record");
                }
                ++output_record_count;
                bcf_destroy(encoded.record.value); encoded.record.value = nullptr;
            },
            [](const GenotypeDecoded& decoded) { return genotype_record_bytes(decoded.record); },
            [](const GenotypeComputed& computed) { return genotype_record_bytes(computed.record); },
            [](const GenotypeEncoded& encoded) {
                return genotype_record_bytes(encoded.record) + encoded.text.size();
            });
        const auto metrics = pipeline.run();
        options.pipeline_used = true;
        options.pipeline_decoded_items = metrics.decoded_items;
        options.pipeline_computed_items = metrics.computed_items;
        options.pipeline_encoded_items = metrics.encoded_items;
        options.pipeline_decoded_bytes = metrics.decoded_bytes;
        options.pipeline_computed_bytes = metrics.computed_bytes;
        options.pipeline_encoded_bytes = metrics.encoded_bytes;
        options.pipeline_peak_decoded_bytes = metrics.peak_decoded_bytes;
        options.pipeline_peak_computed_bytes = metrics.peak_computed_bytes;
        options.pipeline_peak_encoded_bytes = metrics.peak_encoded_bytes;
        options.pipeline_stage_capacity_bytes = stage_capacity;
        bcf_close(output);
        output_file = nullptr;
        std::string index_path;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0,
                                     &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        write_genotype_stream_manifest(options, resources, manifest_path, index_path,
                                       input_paths.size(), input_records, output_record_count,
                                       skipped_reference_blocks, interval_skipped, indexed_inputs,
                                       indexed_interval_queries, span_probe_records,
                                       expected_sample_count);
        std::cout << "{\"tool\":\"GenotypeGVCFs\",\"status\":\"contract-compatible\",\"stream_by_locus\":true,\"input_records\":"
                  << input_records << ",\"output_records\":" << output_record_count
                  << ",\"streamed_loci\":" << options.streamed_loci
                  << ",\"streamed_peak_host_bytes\":" << options.streamed_peak_host_bytes
                  << ",\"stream_max_inflight_records\":" << options.stream_max_inflight_records
                  << ",\"genomicsdb_bridge\":" << (genomicsdb_bridge_used ? "true" : "false")
                  << ",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
                  << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
                  << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items << "}\n";
        bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        if (output_file != nullptr) bcf_close(output_file);
        if (output_header) bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        throw;
    }
}

int run_tool(Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    if (options.stream_by_locus)
        return run_streaming_genotype_gvcf(options, resources);
    const auto kokkos_default_concurrency =
        Kokkos::DefaultExecutionSpace{}.concurrency();
    bcf_hdr_t* output_header = nullptr;
    faidx_t* reference_index = nullptr;
    std::vector<Record> records;
    std::uint64_t input_records = 0;
    std::uint64_t skipped_reference_blocks = 0;
    std::uint64_t interval_skipped = 0;
    std::uint64_t indexed_inputs = 0;
    std::uint64_t indexed_interval_queries = 0;
    std::uint64_t native_record_indexed_inputs = 0;
    std::uint64_t native_record_index_skipped_inputs = 0;
    fastgatk::io::IntervalFileStats interval_file_stats;
    std::vector<Region> regions;
    std::vector<Region> excluded_regions;
    bool intervals_parsed = false;
    int expected_sample_count = -1;
    std::vector<std::string> expected_sample_names;
    std::map<std::string, int> output_sample_indices;
    try {
        const auto input_paths = expand_genomicsdb_inputs(options.inputs, options.regions);
        NativeRecordIndex native_record_index;
        for (const auto& input : options.inputs) {
            if (input.rfind("gendb://", 0) != 0) continue;
            const auto workspace = std::filesystem::path(input.substr(std::strlen("gendb://")));
            const auto loaded = load_native_record_index(workspace);
            native_record_index.insert(loaded.begin(), loaded.end());
        }
        const bool has_genomicsdb_input = std::any_of(
            options.inputs.begin(), options.inputs.end(), [](const std::string& input) {
                return input.rfind("gendb://", 0) == 0;
            });
        if (!options.reference.empty()) {
            reference_index = fai_load(options.reference.c_str());
            if (!reference_index && options.include_non_variant_sites)
                throw std::runtime_error("BAD_INPUT: cannot load reference FASTA index for non-variant sites");
        }
        for (const auto& input_path : input_paths) {
            htsFile* input = bcf_open(input_path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/GVCF: " + input_path);
            bcf_hdr_t* header = bcf_hdr_read(input);
            if (!header) { bcf_close(input); throw std::runtime_error("BAD_INPUT: cannot read VCF header: " + input_path); }
            if (options.genotype_assignment_method == "USE_POSTERIORS_ANNOTATION" &&
                bcf_hdr_id2int(header, BCF_DT_ID, "PP") < 0) {
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error(
                    "BAD_INPUT: USE_POSTERIORS_ANNOTATION requires FORMAT/PP in every input header");
            }
            std::vector<int> output_samples;
            output_samples.reserve(static_cast<std::size_t>(header->n[BCF_DT_SAMPLE]));
            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (name == nullptr || *name == '\0') {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("BAD_INPUT: GVCF header contains an invalid sample name");
                }
                auto found = output_sample_indices.find(name);
                if (found == output_sample_indices.end()) {
                    const int next = static_cast<int>(output_sample_indices.size());
                    output_sample_indices.emplace(name, next);
                    expected_sample_names.emplace_back(name);
                    output_samples.push_back(next);
                } else {
                    output_samples.push_back(found->second);
                }
            }
            expected_sample_count = static_cast<int>(expected_sample_names.size());
            if (!output_header) {
                output_header = bcf_hdr_dup(header);
                if (!output_header) { bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate VCF header"); }
                bcf_hdr_append(output_header, "##fastgatk_genotype_gvcfs_status=reference-block-materialization");
                ensure_gatk_low_qual_filter_declaration(output_header, options);
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "GQ") < 0)
                    bcf_hdr_append(output_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "RGQ") < 0)
                    bcf_hdr_append(output_header, "##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description=\"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)\">");
                if (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
                    !options.gatk_annotation_compatibility) {
                    if (!has_header_field(output_header, BCF_HL_FMT, "GP"))
                        bcf_hdr_append(output_header, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=Genotype posterior in Phred Scale>");
                    if (!has_header_field(output_header, BCF_HL_FMT, "PG"))
                        bcf_hdr_append(output_header, "##FORMAT=<ID=PG,Number=G,Type=Float,Description=Genotype priors in Phred Scale>");
                }
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "DP") < 0)
                    bcf_hdr_append(output_header, options.gatk_annotation_compatibility
                        ? kGatkStandardDpInfoLine
                        : "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MLEAC") < 0)
                    bcf_hdr_append(output_header, options.gatk_annotation_compatibility
                        ? kGatkMleacInfoLine
                        : "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=Maximum likelihood allele count>");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MLEAF") < 0)
                    bcf_hdr_append(output_header, options.gatk_annotation_compatibility
                        ? kGatkMleafInfoLine
                        : "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=Maximum likelihood allele frequency>");
                if (options.annotate_with_num_discovered_alleles &&
                    bcf_hdr_id2int(output_header, BCF_DT_ID, "NDA") < 0)
                    bcf_hdr_append(output_header, options.gatk_annotation_compatibility
                        ? kGatkNdaInfoLine
                        : "##INFO=<ID=NDA,Number=1,Type=Integer,Description=Number of discovered alternate alleles>");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "FS") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "MQ") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "QD") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "SOR") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "ExcessHet") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "InbreedingCoeff") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "RCQ") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=RCQ,Number=1,Type=Float,Description=Joint cross-sample reference-confidence quality>");
                if (bcf_hdr_id2int(output_header, BCF_DT_ID, "RCP") < 0)
                    bcf_hdr_append(output_header, "##INFO=<ID=RCP,Number=1,Type=Float,Description=Joint probability that all samples are hom-reference>");
            }
            // Merge contig dictionaries before materializing records.  A
            // workspace may contain disjoint per-shard headers; output RID
            // values must be assigned from one global dictionary rather than
            // reusing header-local numeric IDs.
            merge_contig_dictionary(output_header, header);
            for (int sample = 0; sample < header->n[BCF_DT_SAMPLE]; ++sample) {
                const char* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
                if (bcf_hdr_id2int(output_header, BCF_DT_SAMPLE, name) < 0 &&
                    bcf_hdr_add_sample(output_header, name) != 0) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge VCF sample header");
                }
            }
            if (bcf_hdr_sync(output_header) != 0) {
                bcf_hdr_destroy(header);
                bcf_close(input);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync merged VCF header");
            }
            if (!intervals_parsed) {
                regions.reserve(options.regions.size());
                bool first_selector = true;
                for (const auto& text : options.regions) {
                    append_region_selector_with_rule(text, header, regions, interval_file_stats,
                                                      options.interval_set_rule, first_selector);
                    first_selector = false;
                }
                for (const auto& text : options.excluded_regions)
                    append_region_selector(text, header, excluded_regions, interval_file_stats);
                normalize_regions(regions);
                normalize_regions(excluded_regions);
                intervals_parsed = true;
            }
            bool skip_native_records = false;
            if (!regions.empty()) {
                const auto indexed = native_record_index.find(input_path);
                if (indexed != native_record_index.end()) {
                    ++native_record_indexed_inputs;
                    skip_native_records = !native_record_index_overlaps(indexed->second, header, regions);
                    if (skip_native_records) ++native_record_index_skipped_inputs;
                }
            }
            bcf1_t* record = bcf_init();
            if (!record) { bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed"); }
            const auto process_record = [&](bcf1_t* record) {
                ++input_records;
                bcf_unpack(record, BCF_UN_ALL);
                if (!in_regions(header, record, regions, excluded_regions)) {
                    ++interval_skipped;
                    return;
                }
                std::vector<int> selected;
                int non_ref_index = -1;
                for (int allele = 1; allele < record->n_allele; ++allele) {
                    if (record->d.allele[allele] == nullptr) continue;
                    if (std::strcmp(record->d.allele[allele], "<NON_REF>") == 0)
                        non_ref_index = allele;
                    else
                        selected.push_back(allele);
                }
                const bool reference_only = selected.empty() && non_ref_index >= 0;
                if (selected.empty() && non_ref_index < 0) { ++skipped_reference_blocks; return; }
                // Keep the symbolic allele in the staged record until all
                // shards have formed their concrete ALT union.  This lets a
                // pure reference block contribute its REF/<NON_REF>
                // likelihoods to a concrete variant emitted by another
                // sample.
                if (non_ref_index >= 0) selected.push_back(non_ref_index);
                const auto* duplicate = bcf_dup(record);
                if (!duplicate) { bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input); throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_dup failed"); }
                auto* materialized = const_cast<bcf1_t*>(duplicate);
                const int original_alleles = record->n_allele;
                Record staged;
                extract_materialized_fields(header, record, selected, original_alleles,
                                            output_samples, staged);
                staged.reference_block = reference_only;
                std::string alleles = record->d.allele[0] == nullptr ? "N" : record->d.allele[0];
                for (const auto allele : selected) {
                    alleles.push_back(',');
                    alleles += record->d.allele[allele] == nullptr ? "N" : record->d.allele[allele];
                }
                if (bcf_update_alleles_str(output_header, materialized, alleles.c_str()) != 0) {
                    bcf_destroy(materialized); bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot materialize VCF alleles");
                }
                int32_t* info_ad = nullptr;
                int info_ad_count = 0;
                const auto info_ad_length = bcf_get_info_int32(output_header, materialized, "AD", &info_ad, &info_ad_count);
                if (info_ad_length >= original_alleles && info_ad_count >= original_alleles) {
                    std::vector<int32_t> compact_info_ad(static_cast<std::size_t>(selected.size()) + 1,
                                                         bcf_int32_missing);
                    compact_info_ad[0] = info_ad[0];
                    for (std::size_t index = 0; index < selected.size(); ++index)
                        compact_info_ad[index + 1] = info_ad[selected[index]];
                    if (bcf_update_info_int32(output_header, materialized, "AD",
                                              compact_info_ad.data(), static_cast<int>(compact_info_ad.size())) != 0) {
                        free(info_ad); bcf_destroy(materialized); bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot compact AD INFO field");
                    }
                }
                free(info_ad);
                bcf_unpack(materialized, BCF_UN_STR);
                const auto* contig_name = record->rid >= 0
                    ? bcf_hdr_id2name(header, record->rid) : nullptr;
                const std::string staged_contig = contig_name == nullptr ? std::string{} :
                    std::string(contig_name);
                const auto output_rid = contig_name == nullptr
                    ? -1 : bcf_hdr_name2id(output_header, contig_name);
                if (output_rid < 0) {
                    bcf_destroy(materialized); bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot map VCF contig into merged header");
                }
                materialized->rid = output_rid;
                staged.value = materialized;
                staged.rid = output_rid;
                staged.pos = static_cast<int>(materialized->pos);
                staged.alleles.reserve(selected.size() + 1);
                staged.alleles.emplace_back(record->d.allele[0] == nullptr ? "N" : record->d.allele[0]);
                for (const auto allele : selected)
                    staged.alleles.emplace_back(record->d.allele[allele] == nullptr ? "N" : record->d.allele[allele]);
                staged.key = record_key(materialized);
                if (options.include_non_variant_sites && reference_only) {
                    // GATK expands a pure gVCF reference block into one
                    // REF-only output record per coordinate for
                    // --include-non-variant-sites.  Preserve the source
                    // sample/FORMAT arrays while cloning the small record
                    // envelope and removing END from each emitted site.
                    int32_t* end_value = nullptr;
                    int end_count = 0;
                    const auto end_length = bcf_get_info_int32(header, materialized,
                                                               "END", &end_value, &end_count);
                    const auto start = materialized->pos;
                    int end = start + 1;
                    if (end_length > 0 && end_count > 0 && end_value[0] > start + 1)
                        end = end_value[0];
                    free(end_value);
                    int32_t* info_min_dp = nullptr;
                    int info_min_dp_count = 0;
                    (void)bcf_get_info_int32(header, materialized, "MIN_DP",
                                             &info_min_dp, &info_min_dp_count);
                    int32_t* format_min_dp = nullptr;
                    int format_min_dp_count = 0;
                    (void)bcf_get_format_int32(header, materialized, "MIN_DP",
                                               &format_min_dp, &format_min_dp_count);
                    const bool has_min_dp = info_min_dp_count > 0 || format_min_dp_count > 0;
                    const auto min_dp_at = [&](const std::size_t sample) {
                        if (info_min_dp_count > 0)
                            return info_min_dp[std::min<std::size_t>(sample,
                                static_cast<std::size_t>(info_min_dp_count - 1))];
                        if (format_min_dp_count > 0)
                            return format_min_dp[std::min<std::size_t>(sample,
                                static_cast<std::size_t>(format_min_dp_count - 1))];
                        return staged.dp.empty() ? static_cast<int32_t>(0) : staged.dp[sample];
                    };
                    for (int position = start; position < end; ++position) {
                        auto* copy = bcf_dup(materialized);
                        if (!copy) {
                            bcf_destroy(materialized); bcf_destroy(record); bcf_hdr_destroy(header);
                            bcf_close(input);
                            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot expand reference block");
                        }
                        copy->pos = position;
                        // The merged REF of a locus is the REF of the record
                        // that STARTS there, and only an interior coordinate
                        // falls back to the reference base: the merger is handed
                        // `ref.getBase()` and uses it where no record starts
                        // (ReferenceConfidenceVariantContextMerger.merge()).
                        // Measured on a block whose REF is C over an all-A
                        // reference: GATK publishes C at the block's own
                        // coordinate and A at the interior ones; taking the FASTA
                        // base everywhere was the divergence behind 82 positions
                        // of GATK's own chr20 corpus.
                        std::string expanded_reference = staged.alleles.front();
                        if (reference_index && position != start) {
                            const char* contig = staged_contig.empty() ? nullptr : staged_contig.c_str();
                            int fetched_length = 0;
                            char* fetched = contig == nullptr ? nullptr :
                                faidx_fetch_seq(reference_index, contig, position, position,
                                                &fetched_length);
                            if (fetched != nullptr && fetched_length == 1)
                                expanded_reference.assign(1, fetched[0]);
                            free(fetched);
                        }
                        if (!expanded_reference.empty()) {
                            const auto expanded_alleles = expanded_reference + ",<NON_REF>";
                            (void)bcf_update_alleles_str(output_header, copy,
                                                         expanded_alleles.c_str());
                        }
                        (void)bcf_update_info_int32(output_header, copy, "END", nullptr, 0);
                        if (!staged.dp.empty()) {
                            const int32_t site_dp = has_min_dp ? min_dp_at(0) : staged.dp[0];
                            (void)bcf_update_info_int32(output_header, copy, "DP", &site_dp, 1);
                        }
                        (void)bcf_update_info_int32(output_header, copy, "MIN_DP", nullptr, 0);
                        (void)bcf_update_info_int32(output_header, copy, "AD", nullptr, 0);
                        Record expanded = staged;
                        expanded.value = copy;
                        expanded.pos = position;
                        expanded.alleles.front() = expanded_reference;
                        if (has_min_dp &&
                            (info_min_dp_count > 0 ||
                             format_min_dp_count >= static_cast<int>(expanded.dp.size())))
                            for (std::size_t sample = 0; sample < expanded.dp.size(); ++sample)
                                expanded.dp[sample] = min_dp_at(sample);
                        if (has_min_dp && expanded.min_dp.size() == expanded.dp.size())
                            for (std::size_t sample = 0; sample < expanded.min_dp.size(); ++sample)
                                expanded.min_dp[sample] = min_dp_at(sample);
                        expanded.key = record_key(copy);
                        records.push_back(std::move(expanded));
                    }
                    free(info_min_dp);
                    free(format_min_dp);
                    bcf_destroy(materialized);
                    return;
                }
                records.push_back(std::move(staged));
            };

            // Indexed traversal is the default whenever a caller supplies
            // intervals and the input has a compatible CSI/TBI index.  It
            // avoids decoding unrelated chromosomes/blocks and therefore
            // makes the Host byte budget meaningful for GenotypeGVCFs.  If an
            // index is absent or cannot be opened, retain the sequential
            // path and count records rejected by the interval predicate.
            bool used_index = false;
            if (skip_native_records) {
                // The workspace span index proved that this file has no
                // record overlapping the requested intervals.  Its header
                // and sample names were already consumed above, so skipping
                // the body preserves output columns while avoiding BGZF
                // decoding of an unrelated shard.
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(input);
                continue;
            }
            if (!regions.empty()) {
                const auto index = load_genotype_traversal_index(input_path);
                if (index != nullptr) {
                    used_index = true;
                    ++indexed_inputs;
                    for (const auto& region : regions) {
                        ++indexed_interval_queries;
                        const int index_tid = index->tid_for(header, region.rid, region.contig);
                        if (index_tid < 0) continue;
                        int status = 0;
                        if (index->tabix != nullptr) {
                            hts_itr_t* iterator = tbx_itr_queryi(index->tabix, index_tid,
                                                                 region.begin, region.end);
                            if (iterator == nullptr) continue;
                            kstring_t line{0, 0, nullptr};
                            while ((status = tbx_itr_next(input, index->tabix, iterator,
                                                          &line)) >= 0) {
                                if (vcf_parse(&line, header, record) < 0) {
                                    free(line.s);
                                    hts_itr_destroy(iterator);
                                    bcf_destroy(record);
                                    bcf_hdr_destroy(header);
                                    bcf_close(input);
                                    throw std::runtime_error(
                                        "BAD_INPUT: indexed GenotypeGVCFs VCF record parse failed");
                                }
                                process_record(record);
                            }
                            free(line.s);
                            hts_itr_destroy(iterator);
                        } else {
                            hts_itr_t* iterator = bcf_itr_queryi(index->index, index_tid,
                                                                 region.begin, region.end);
                            if (iterator == nullptr) continue;
                            while ((status = bcf_itr_next(input, iterator, record)) >= 0)
                                process_record(record);
                            hts_itr_destroy(iterator);
                        }
                        if (status < -1) {
                            bcf_destroy(record);
                            bcf_hdr_destroy(header);
                            bcf_close(input);
                            throw std::runtime_error("BAD_INPUT: indexed GenotypeGVCFs traversal failed");
                        }
                    }
                }
            }
            if (!used_index) {
                int read_status = 0;
                while ((read_status = bcf_read(input, header, record)) == 0)
                    process_record(record);
                // -1 is end of input; -2 is a record htslib could not parse.
                // Treating the parse failure as end of input would silently
                // discard this record and everything after it.
                if (read_status < -1) {
                    bcf_destroy(record);
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error(
                        "BAD_INPUT: GenotypeGVCFs input record parse failed");
                }
            }
            bcf_destroy(record);
            bcf_hdr_destroy(header);
            bcf_close(input);
        }
        // Materialize cross-sample reference-confidence points before sorting
        // and locus grouping.  This is intentionally reference-backed: a
        // block's REF base is not recoverable from its single start allele at
        // an interior split coordinate.
        split_reference_blocks_at_variants(output_header, reference_index, records);
        if (options.include_non_variant_sites) {
            materialize_spanning_loci(output_header, reference_index, records, regions);
            if (!regions.empty()) {
                // GATK traverses only the requested intervals, so the
                // per-coordinate dense record supply is clipped to them: a
                // reference block whose END reaches past the interval used to be
                // published in full (measured: 4809 in-block coordinates for a
                // 1000 bp window of GATK's own chr20 gVCF, against GATK's 1001).
                std::vector<Record> clipped;
                clipped.reserve(records.size());
                for (auto& record : records) {
                    if (record.value != nullptr && record.value->rid >= 0 &&
                        starts_in_regions(output_header, record.value, regions)) {
                        clipped.push_back(std::move(record));
                    } else {
                        bcf_destroy(record.value);
                        record.value = nullptr;
                    }
                }
                records = std::move(clipped);
            }
        }
        std::sort(records.begin(), records.end(), [](const Record& left, const Record& right) {
            if (left.rid != right.rid) return left.rid < right.rid;
            if (left.pos != right.pos) return left.pos < right.pos;
            return left.key < right.key;
        });
        std::vector<Record> unique;
        for (std::size_t begin = 0; begin < records.size();) {
            std::size_t end = begin + 1;
            while (end < records.size() && records[end].key == records[begin].key) ++end;
            // Build a deterministic concrete ALT union and keep the symbolic
            // <NON_REF> sentinel last.  Reference-only blocks can therefore
            // be remapped into a variant shard before their unseen-ALT
            // likelihoods are projected away.
            std::vector<std::string> union_alleles{records[begin].alleles.front()};
            bool has_non_ref = false;
            // See the streaming merge: the '*' stays in the union and its
            // ownership is decided where GATK decides it, in the compute stage
            // (apply_gatk_output_allele_subset), from the deletions earlier
            // loci actually emitted.
            for (std::size_t index = begin; index < end; ++index) {
                for (std::size_t allele = 1; allele < records[index].alleles.size(); ++allele) {
                    const auto& name = records[index].alleles[allele];
                    if (name == "<NON_REF>") {
                        has_non_ref = true;
                        continue;
                    }
                    if (std::find(union_alleles.begin() + 1, union_alleles.end(), name) == union_alleles.end())
                        union_alleles.push_back(name);
                }
            }
            if (has_non_ref) union_alleles.push_back("<NON_REF>");
            std::vector<std::string> concrete_alleles;
            concrete_alleles.reserve(union_alleles.size());
            for (const auto& allele : union_alleles)
                if (allele != "<NON_REF>") concrete_alleles.push_back(allele);
            if (concrete_alleles.size() < 2) {
                // A group containing only reference blocks has no concrete
                // variant to emit by default.  With GATK's explicit
                // --include-non-variant-sites flag, compact each source block
                // to a REF-only record and let the normal sample merge retain
                // disjoint samples.
                if (options.include_non_variant_sites) {
                    for (std::size_t index = begin; index < end; ++index)
                        materialize_reference_only(output_header, records[index],
                                                   options.gatk_annotation_compatibility);
                    unique.push_back(std::move(records[begin]));
                    std::vector<Record*> merged_group{&unique.back()};
                    std::map<int, bool> seen_samples;
                    for (const auto sample : unique.back().output_samples) seen_samples[sample] = true;
                    for (std::size_t index = begin + 1; index < end; ++index) {
                        bool overlaps = false;
                        for (const auto sample : records[index].output_samples)
                            if (seen_samples.contains(sample)) { overlaps = true; break; }
                        if (overlaps) {
                            bcf_destroy(records[index].value);
                            continue;
                        }
                        for (const auto sample : records[index].output_samples) seen_samples[sample] = true;
                        merged_group.push_back(&records[index]);
                    }
                    merge_sample_fields(output_header, merged_group);
                    for (std::size_t index = begin + 1; index < end; ++index) {
                        if (std::find(merged_group.begin(), merged_group.end(), &records[index]) == merged_group.end())
                            continue;
                        bcf_destroy(records[index].value);
                    }
                    begin = end;
                    continue;
                }
                // Otherwise consume and release it after preserving telemetry
                // for ignored reference blocks.
                for (std::size_t index = begin; index < end; ++index)
                    bcf_destroy(records[index].value);
                begin = end;
                continue;
            }
            for (std::size_t index = begin; index < end; ++index) {
                remap_record_to_allele_union(output_header, records[index], union_alleles);
                // See the streaming merge: keep the merged source call before the
                // PL-derived assignment consumes it.
                records[index].source_gt = records[index].gt;
                records[index].source_alleles = records[index].alleles;
                remove_non_ref_allele(
                    output_header, records[index], concrete_alleles,
                    options.genotype_assignment_method == "BEST_MATCH_TO_ORIGINAL");
                if (options.genotype_assignment_method != "DO_NOT_ASSIGN_GENOTYPES" &&
                    options.genotype_assignment_method != "BEST_MATCH_TO_ORIGINAL" &&
                    options.genotype_assignment_method != "USE_POSTERIORS_ANNOTATION")
                    derive_gt_gq_from_pl(records[index]);
            }
            unique.push_back(std::move(records[begin]));
            std::vector<Record*> merged_group{&unique.back()};
            std::map<int, bool> seen_samples;
            for (const auto sample : unique.back().output_samples) seen_samples[sample] = true;
            for (std::size_t index = begin + 1; index < end; ++index) {
                bool overlaps = false;
                for (const auto sample : records[index].output_samples)
                    if (seen_samples.contains(sample)) { overlaps = true; break; }
                if (overlaps) {
                    bcf_destroy(records[index].value);
                    continue;
                }
                for (const auto sample : records[index].output_samples) seen_samples[sample] = true;
                merged_group.push_back(&records[index]);
            }
            merge_sample_fields(output_header, merged_group);
            for (std::size_t index = begin + 1; index < end; ++index) {
                if (std::find(merged_group.begin(), merged_group.end(), &records[index]) == merged_group.end())
                    continue;
                bcf_destroy(records[index].value);
            }
            begin = end;
        }
        records.clear();
        records = std::move(unique);
        const char* mode = suffix(options.output, ".gz") ? "wz" : "w";
        htsFile* output = bcf_open(options.output.c_str(), mode);
        if (!output) throw std::runtime_error("cannot open output VCF: " + options.output);
        if (options.gatk_annotation_compatibility) {
            kstring_t formatted_header{0, 0, nullptr};
            if (bcf_hdr_format(output_header, 0, &formatted_header) < 0 ||
                write_vcf_text_line(output, gatk_compatible_header_text(
                    std::string(formatted_header.s == nullptr ? "" : formatted_header.s,
                                formatted_header.l),
                    input_paths.empty()
                        ? std::string{}
                        : input_pass_filter_line(input_paths.front()))) != 0) {
                free(formatted_header.s);
                bcf_close(output);
                throw std::runtime_error("cannot write GATK-compatible VCF header");
            }
            free(formatted_header.s);
        } else if (bcf_hdr_write(output, output_header) != 0) {
            bcf_close(output);
            throw std::runtime_error("cannot write output VCF header");
        }
        // Match the process-wide deterministic GATK Utils RNG used by
        // QualByDepth.  The state is shared by the single compute worker, so
        // each high-QD annotation consumes the same Gaussian sequence as Java
        // while the caller thread remains the only VCF writer.
        GatkJavaRandom gatk_random;
        std::uint64_t output_record_count = 0;
        std::vector<Record> pipeline_records = std::move(records);
        std::size_t decode_index = 0;
        // GenotypingEngine's per-traversal emitted-deletion state
        // (GenotypingEngine.java:52).  It is read and written only by the
        // single compute worker below, which visits loci in traversal order.
        EmittedDeletions upstream_deletions;
        const auto stage_capacity = std::max<std::uint64_t>(
            64ULL * 1024ULL * 1024ULL, resources.safe_memory_budget_bytes());
        using Pipeline = fastgatk::runtime::ThreeStagePipeline<
            GenotypeDecoded, GenotypeComputed, GenotypeEncoded>;
        Pipeline pipeline(
            Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
            [&]() -> std::optional<GenotypeDecoded> {
                if (decode_index >= pipeline_records.size()) return std::nullopt;
                GenotypeDecoded decoded;
                auto& source = pipeline_records[decode_index++];
                decoded.record = std::move(source);
                // Record owns a raw bcf1_t*.  Its implicit move constructor
                // cannot null the source pointer, so clear it explicitly to
                // prevent the post-pipeline cleanup from double-freeing the
                // value already released by the sink.
                source.value = nullptr;
                return decoded;
            },
            [&](GenotypeDecoded decoded) -> std::optional<GenotypeComputed> {
                auto record = std::move(decoded.record);
                if (gatk_skips_regenotyping(output_header, record, options)) {
                    // GenotypeGVCFsEngine.java:160 sends a variant record whose
                    // merged INFO/DP is not positive to the passthrough arm at
                    // :174, and :181 then rejects the locus because the depth
                    // test fails again; the default traversal therefore emits no
                    // record for it, even with a confident call in the
                    // genotypes.  No deletion is recorded either, because GATK
                    // never reaches recordDeletions() on this path.
                    bcf_destroy(record.value);
                    record.value = nullptr;
                    return std::nullopt;
                }
                std::vector<double> posterior_priors;
                annotate_num_discovered_alleles(output_header, record, options);
                (void)apply_gatk_max_alternate_alleles(output_header, record, options);
                update_cohort_af_annotations(output_header, record, options);
                if (!apply_gatk_output_allele_subset(output_header, record, options,
                                                     upstream_deletions))
                    return std::nullopt;
                if (record.finalized_monomorphic_ref) {
                    // The REF-only call is complete: GATK emits it directly from
                    // regenotypeVC, so no site or standard annotation runs.
                    apply_gatk_annotation_compatibility(output_header, record, options);
                    GenotypeComputed computed;
                    computed.record = std::move(record);
                    return computed;
                }
                if (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
                    !options.gatk_annotation_compatibility &&
                    record.allele_count >= 2 && record.ploidy > 0)
                    posterior_priors = build_log10_genotype_priors(
                        record.alleles, record.ploidy, options.snp_heterozygosity,
                        options.indel_heterozygosity, options.use_genotype_priors);
                apply_genotype_assignment(output_header, record, options, posterior_priors);
                normalize_haplotype_caller_phasing(output_header, record);
                update_site_annotations(output_header, record);
                update_total_depth_annotation(output_header, record, options);
                update_cross_sample_reference_confidence(output_header, record, options);
                update_site_quality_from_posteriors(record, options);
                update_gatk_standard_annotations(output_header, record, gatk_random);
                update_excess_het_annotation(output_header, record);
                update_inbreeding_coeff_annotation(output_header, record);
                apply_gatk_annotation_compatibility(output_header, record, options);
                GenotypeComputed computed;
                computed.record = std::move(record);
                return computed;
            },
            [&](GenotypeComputed computed) -> std::optional<GenotypeEncoded> {
                GenotypeEncoded encoded;
                encoded.record = std::move(computed.record);
                if (options.gatk_annotation_compatibility) {
                    kstring_t formatted{0, 0, nullptr};
                    if (vcf_format1(output_header, encoded.record.value, &formatted) != 0) {
                        free(formatted.s);
                        throw std::runtime_error("cannot format GATK-compatible VCF record");
                    }
                    encoded.text = gatk_compatible_record_text(
                        std::string(formatted.s == nullptr ? "" : formatted.s, formatted.l));
                    free(formatted.s);
                }
                return encoded;
            },
            [&](GenotypeEncoded encoded) {
                if (options.only_output_calls_starting_in_intervals &&
                    !starts_in_regions(output_header, encoded.record.value, regions)) {
                    ++options.output_interval_skipped;
                    bcf_destroy(encoded.record.value);
                    encoded.record.value = nullptr;
                    return;
                }
                if (options.gatk_annotation_compatibility) {
                    if (write_vcf_text_line(output, encoded.text) != 0) {
                        bcf_destroy(encoded.record.value);
                        encoded.record.value = nullptr;
                        throw std::runtime_error("cannot write GATK-compatible VCF record");
                    }
                } else if (bcf_write(output, output_header, encoded.record.value) != 0) {
                    bcf_destroy(encoded.record.value);
                    encoded.record.value = nullptr;
                    throw std::runtime_error("cannot write output VCF record");
                }
                bcf_destroy(encoded.record.value);
                encoded.record.value = nullptr;
                ++output_record_count;
            },
            [](const GenotypeDecoded& decoded) {
                return genotype_record_bytes(decoded.record);
            },
            [](const GenotypeComputed& computed) {
                return genotype_record_bytes(computed.record);
            },
            [](const GenotypeEncoded& encoded) {
                return genotype_record_bytes(encoded.record) +
                       static_cast<std::uint64_t>(encoded.text.size());
            });
        const auto pipeline_metrics = pipeline.run();
        destroy_records(pipeline_records);
        options.pipeline_used = true;
        options.pipeline_decoded_items = pipeline_metrics.decoded_items;
        options.pipeline_computed_items = pipeline_metrics.computed_items;
        options.pipeline_encoded_items = pipeline_metrics.encoded_items;
        options.pipeline_decoded_bytes = pipeline_metrics.decoded_bytes;
        options.pipeline_computed_bytes = pipeline_metrics.computed_bytes;
        options.pipeline_encoded_bytes = pipeline_metrics.encoded_bytes;
        options.pipeline_peak_decoded_bytes = pipeline_metrics.peak_decoded_bytes;
        options.pipeline_peak_computed_bytes = pipeline_metrics.peak_computed_bytes;
        options.pipeline_peak_encoded_bytes = pipeline_metrics.peak_encoded_bytes;
        options.pipeline_stage_capacity_bytes = stage_capacity;
        bcf_close(output);
        std::string index_path;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VCF index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const bool output_complete = file_complete(options.output) &&
            (index_path.empty() || file_complete(index_path));
        if (!output_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: primary output or index is incomplete");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"GenotypeGVCFs\",\"implementation\":\"fastgatk-genotype-gvcf\",\"status\":\"contract-compatible\"," 
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"gatk_parameter_aliases\":true,\"reference_block_materialization\":true,\"nonref_pl_projection\":true,\"include_non_variant_sites\":"
                 << (options.include_non_variant_sites ? "true" : "false") << ",\"format_compaction\":true,\"multi_allelic_materialization\":true,\"allele_union\":true,\"joint_genotype_from_pl\":true,\"kokkos_genotype_pl\":true,\"arbitrary_ploidy_gt_gq\":true,\"kokkos_allele_counts\":true,\"site_annotations\":true,\"excess_het_exact_diploid\":true,\"inbreeding_coeff_exact_diploid\":true,\"multi_sample_joint_merge\":true,\"cross_sample_reference_block_split\":true,\"cross_sample_reference_confidence\":true,\"genotype_prior_calculator_assumingHW\":true,\"allele_frequency_calculator_em\":true,\"gq_only_homref_af_approximation\":true,\"posterior_qual_opt_in\":true,\"output_allele_subset\":true,\"standard_confidence_allele_pruning\":true,\"max_alternate_alleles_likelihood_subset\":true,\"cohort_af_calculator\":true,\"spanning_deletion_nonvariant_set\":true,\"spanning_deletion_orphan_cleanup\":true,\"spanning_deletion_posterior_qual\":true,\"genotype_assignment_method\":true,\"posterior_genotype_assignment\":"
                 << (options.genotype_assignment_method == "USE_POSTERIOR_PROBABILITIES" &&
                     !options.gatk_annotation_compatibility ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"set_to_no_call_no_annotations\":true"
                 << ",\"all_format_annotation_cleanup\":true"
                 << ",\"best_match_to_original\":true"
                 << ",\"use_posteriors_annotation\":"
                 << (options.genotype_assignment_method == "USE_POSTERIORS_ANNOTATION" ? "true" : "false")
                 << ",\"gatk_annotation_compatibility\":"
                 << (options.gatk_annotation_compatibility ? "true" : "false")
                 << ",\"interval_exclusion\":true"
                 << ",\"output_interval_starts_in\":true"
                 << ",\"num_discovered_alleles_annotation\":true"
                 << ",\"genomicsdb_workspace\":"
                 << (has_genomicsdb_input ? "true" : "false") << ",\"genomicsdb_bridge\":"
                 << (genomicsdb_bridge_used ? "true" : "false") << ",\"vcf_index\":"
                 << (index_path.empty() ? "false" : "true") << ",\"interval_subset\":"
                 << (!options.regions.empty() ? "true" : "false") << ",\"indexed_interval_traversal\":true"
                 << ",\"native_record_index\":" << (!native_record_index.empty() ? "true" : "false")
                 << ",\"native_record_indexed_inputs\":" << native_record_indexed_inputs
                 << ",\"native_record_index_skipped_inputs\":" << native_record_index_skipped_inputs
                 << ",\"kokkos_default_concurrency\":" << kokkos_default_concurrency
                 << ",\"three_stage_pipeline\":" << (options.pipeline_used ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":"
                 << (output_complete ? "true" : "false") << "}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":" +
                     (output_complete ? "true" : "false") + "}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << input_paths.size() << ",\"input_records\":" << input_records << ",\"output_records\":" << output_record_count
                 << ",\"genomicsdb_workspace_input_index\":" << (has_genomicsdb_input ? "true" : "false")
                 << ",\"genomicsdb_bridge\":" << (genomicsdb_bridge_used ? "true" : "false")
                 << ",\"skipped_reference_blocks\":" << skipped_reference_blocks
                 << ",\"intervals\":" << options.regions.size()
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"exclude_intervals\":" << options.excluded_regions.size()
                 << ",\"only_output_calls_starting_in_intervals\":"
                 << (options.only_output_calls_starting_in_intervals ? "true" : "false")
                 << ",\"output_interval_skipped\":" << options.output_interval_skipped
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"indexed_inputs\":" << indexed_inputs
                 << ",\"indexed_interval_queries\":" << indexed_interval_queries
                 << ",\"native_record_indexed_inputs\":" << native_record_indexed_inputs
                 << ",\"native_record_index_skipped_inputs\":" << native_record_index_skipped_inputs
                 << ",\"kokkos_default_concurrency\":" << kokkos_default_concurrency
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records
                 << ",\"sample_count\":" << std::max(0, expected_sample_count)
                 << ",\"snp_heterozygosity\":" << options.snp_heterozygosity
                 << ",\"indel_heterozygosity\":" << options.indel_heterozygosity
                 << ",\"heterozygosity_stdev\":" << options.heterozygosity_stdev
                 << ",\"standard_confidence_for_calling\":" << options.standard_confidence_for_calling
                 << ",\"use_genotype_priors\":" << (options.use_genotype_priors ? "true" : "false")
                 << ",\"use_new_qual_calculator\":" << (options.use_new_qual_calculator ? "true" : "false")
                 << ",\"use_posteriors_to_calculate_qual\":"
                 << (options.use_posteriors_to_calculate_qual ? "true" : "false")
                 << ",\"genotype_assignment_method\":\""
                 << json_escape(options.genotype_assignment_method) << "\""
                 << ",\"include_non_variant_sites\":"
                 << (options.include_non_variant_sites ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"genotype_kernel_calls\":" << genotype_kernel_telemetry.calls
                 << ",\"genotype_kernel_prepare_seconds\":" << genotype_kernel_telemetry.prepare_seconds
                 << ",\"genotype_kernel_seconds\":" << genotype_kernel_telemetry.execute_seconds
                 << ",\"genotype_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.execution_space) << "\""
                 << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                 << ",\"pl_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.pl_remap_prepare_seconds
                 << ",\"pl_remap_kernel_seconds\":" << genotype_kernel_telemetry.pl_remap_execute_seconds
                 << ",\"pl_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.pl_remap_execution_space) << "\""
                 << ",\"allele_count_kernel_calls\":" << genotype_kernel_telemetry.allele_count_calls
                 << ",\"allele_count_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_count_prepare_seconds
                 << ",\"allele_count_kernel_seconds\":" << genotype_kernel_telemetry.allele_count_execute_seconds
                 << ",\"allele_count_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_count_execution_space) << "\""
                 << ",\"posterior_kernel_calls\":" << genotype_kernel_telemetry.posterior_calls
                 << ",\"posterior_kernel_prepare_seconds\":" << genotype_kernel_telemetry.posterior_prepare_seconds
                 << ",\"posterior_kernel_seconds\":" << genotype_kernel_telemetry.posterior_execute_seconds
                 << ",\"posterior_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.posterior_execution_space) << "\""
                 << ",\"posterior_samples\":" << genotype_kernel_telemetry.posterior_samples
                 << ",\"cross_sample_reference_kernel_calls\":" << genotype_kernel_telemetry.cross_sample_reference_calls
                 << ",\"cross_sample_reference_kernel_prepare_seconds\":" << genotype_kernel_telemetry.cross_sample_reference_prepare_seconds
                 << ",\"cross_sample_reference_kernel_seconds\":" << genotype_kernel_telemetry.cross_sample_reference_execute_seconds
                 << ",\"cross_sample_reference_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.cross_sample_reference_execution_space) << "\""
                 << ",\"cross_sample_reference_samples\":" << genotype_kernel_telemetry.cross_sample_reference_samples
                 << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
                 << ",\"allele_field_remap_kernel_prepare_seconds\":" << genotype_kernel_telemetry.allele_field_remap_prepare_seconds
                 << ",\"allele_field_remap_kernel_seconds\":" << genotype_kernel_telemetry.allele_field_remap_execute_seconds
                 << ",\"allele_field_remap_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.allele_field_remap_execution_space) << "\""
                 << ",\"cohort_af_kernel_calls\":" << genotype_kernel_telemetry.cohort_calls
                 << ",\"cohort_af_kernel_prepare_seconds\":" << genotype_kernel_telemetry.cohort_prepare_seconds
                 << ",\"cohort_af_kernel_seconds\":" << genotype_kernel_telemetry.cohort_execute_seconds
                 << ",\"cohort_af_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.cohort_execution_space) << "\""
                 << ",\"cohort_af_samples\":" << genotype_kernel_telemetry.cohort_samples
                 << ",\"cohort_af_approximate_gq_samples\":" << genotype_kernel_telemetry.cohort_approximate_gq_samples
                 << ",\"cohort_af_iterations\":" << genotype_kernel_telemetry.cohort_iterations
                 << ",\"cohort_af_converged\":" << genotype_kernel_telemetry.cohort_converged
                 << ",\"output_allele_pruning_calls\":" << genotype_kernel_telemetry.output_allele_pruning_calls
                 << ",\"output_alleles_pruned\":" << genotype_kernel_telemetry.output_alleles_pruned
                 << ",\"orphan_spanning_deletion_loci\":" << options.orphan_spanning_deletion_loci
                 << ",\"max_alternate_alleles\":" << options.max_alternate_alleles
                 << ",\"max_alt_pruning_calls\":" << genotype_kernel_telemetry.max_alt_pruning_calls
                 << ",\"max_alt_alleles_pruned\":" << genotype_kernel_telemetry.max_alt_alleles_pruned
                 << ",\"max_alt_score_kernel_calls\":" << genotype_kernel_telemetry.max_alt_score_kernel_calls
                 << ",\"max_alt_score_kernel_prepare_seconds\":" << genotype_kernel_telemetry.max_alt_score_prepare_seconds
                 << ",\"max_alt_score_kernel_seconds\":" << genotype_kernel_telemetry.max_alt_score_execute_seconds
                 << ",\"max_alt_score_kernel_execution_space\":\"" << json_escape(genotype_kernel_telemetry.max_alt_score_execution_space) << "\""
                 << ",\"annotate_with_num_discovered_alleles\":" << (options.annotate_with_num_discovered_alleles ? "true" : "false")
                 << ",\"num_discovered_alleles_calls\":" << genotype_kernel_telemetry.num_discovered_alleles_calls
                 << ",\"excess_het_calls\":" << genotype_kernel_telemetry.excess_het_calls
                 << ",\"inbreeding_coeff_calls\":" << genotype_kernel_telemetry.inbreeding_coeff_calls
                 << ",\"posterior_assignment_kernel_calls\":" << genotype_kernel_telemetry.posterior_assignment_calls
                 << ",\"posterior_assignment_kernel_prepare_seconds\":" << genotype_kernel_telemetry.posterior_assignment_prepare_seconds
                 << ",\"posterior_assignment_kernel_seconds\":" << genotype_kernel_telemetry.posterior_assignment_execute_seconds
                 << ",\"posterior_assignment_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.posterior_assignment_execution_space) << "\""
                 << ",\"posterior_assignment_samples\":" << genotype_kernel_telemetry.posterior_assignment_samples
                 << ",\"posterior_annotation_assignment_kernel_calls\":"
                 << genotype_kernel_telemetry.posterior_annotation_assignment_calls
                 << ",\"posterior_annotation_assignment_prepare_seconds\":"
                 << genotype_kernel_telemetry.posterior_annotation_assignment_prepare_seconds
                 << ",\"posterior_annotation_assignment_kernel_seconds\":"
                 << genotype_kernel_telemetry.posterior_annotation_assignment_execute_seconds
                 << ",\"posterior_annotation_assignment_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.posterior_annotation_assignment_execution_space) << "\""
                 << ",\"posterior_annotation_assignment_samples\":"
                 << genotype_kernel_telemetry.posterior_annotation_assignment_samples
                 << ",\"genotype_prior_kernel_calls\":"
                 << genotype_kernel_telemetry.genotype_prior_calls
                 << ",\"genotype_prior_kernel_prepare_seconds\":"
                 << genotype_kernel_telemetry.genotype_prior_prepare_seconds
                 << ",\"genotype_prior_kernel_seconds\":"
                 << genotype_kernel_telemetry.genotype_prior_execute_seconds
                 << ",\"genotype_prior_kernel_execution_space\":\""
                 << json_escape(genotype_kernel_telemetry.genotype_prior_execution_space) << "\""
                 << ",\"no_annotation_format_fields\":"
                 << genotype_kernel_telemetry.no_annotation_format_fields
                 << ",\"pipeline_lifecycle\":\""
                 << (options.pipeline_used ? "Host decode->bounded queue->Kokkos compute->encode->sink" : "not-used")
                 << "\",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
                 << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
                 << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
                 << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
                 << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
                 << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
                 << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
                 << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
                 << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
                 << ",\"pipeline_stage_capacity_bytes\":" << options.pipeline_stage_capacity_bytes << "}}\n";
        std::cout << "{\"tool\":\"GenotypeGVCFs\",\"status\":\"contract-compatible\",\"input_records\":"
                  << input_records << ",\"output_records\":" << output_record_count << ",\"skipped_reference_blocks\":"
                  << skipped_reference_blocks << ",\"interval_skipped\":" << interval_skipped
                  << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"sample_count\":" << std::max(0, expected_sample_count)
                  << ",\"pl_remap_kernel_calls\":" << genotype_kernel_telemetry.pl_remap_calls
                  << ",\"allele_field_remap_kernel_calls\":" << genotype_kernel_telemetry.allele_field_remap_calls
                  << ",\"cohort_af\":"
                  << (options.use_new_qual_calculator ? "true" : "false")
                  << ",\"cohort_af_approximate_gq_samples\":"
                  << genotype_kernel_telemetry.cohort_approximate_gq_samples
                  << ",\"excess_het_calls\":" << genotype_kernel_telemetry.excess_het_calls
                  << ",\"inbreeding_coeff_calls\":" << genotype_kernel_telemetry.inbreeding_coeff_calls
                  << ",\"posterior_qual\":"
                  << (options.use_posteriors_to_calculate_qual ? "true" : "false")
                  << ",\"posterior_annotation_assignment_kernel_calls\":"
                  << genotype_kernel_telemetry.posterior_annotation_assignment_calls
                  << ",\"genotype_prior_kernel_calls\":"
                  << genotype_kernel_telemetry.genotype_prior_calls
                  << ",\"no_annotation_format_fields\":"
                  << genotype_kernel_telemetry.no_annotation_format_fields
                  << ",\"pipeline_lifecycle\":\""
                  << (options.pipeline_used ? "Host decode->bounded queue->Kokkos compute->encode->sink" : "not-used")
                  << "\",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
                  << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
                  << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
                  << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
                  << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
                  << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
                  << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
                  << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
                  << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
                  << ",\"pipeline_stage_capacity_bytes\":" << options.pipeline_stage_capacity_bytes
                  << ",\"genotype_assignment_method\":\""
                  << json_escape(options.genotype_assignment_method) << "\""
                  << ",\"gatk_annotation_compatibility\":"
                  << (options.gatk_annotation_compatibility ? "true" : "false")
                  << ",\"genomicsdb_workspace\":" << (has_genomicsdb_input ? "true" : "false")
                  << ",\"genomicsdb_bridge\":" << (genomicsdb_bridge_used ? "true" : "false") << "}\n";
        destroy_records(records);
        bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        return 0;
    } catch (...) {
        destroy_records(records);
        if (output_header) bcf_hdr_destroy(output_header);
        if (reference_index) fai_destroy(reference_index);
        throw;
    }
}

#else
int run_tool(Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for GenotypeGVCFs");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        auto options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        // GenotypeGVCFs has no GATK-facing compute-thread flag, but it still
        // must respect the worker envelope in a SLURM task.  Kokkos otherwise
        // falls back to the host-wide OpenMP team when OMP_NUM_THREADS is not
        // exported by the scheduler.  Prefer an explicit OMP value when one
        // is present, otherwise cap the team to SLURM_CPUS_PER_TASK; both are
        // bounded again by ResourceSnapshot::effective_threads().
        std::size_t requested_threads = resources.slurm_threads;
        if (const char* value = std::getenv("OMP_NUM_THREADS"); value != nullptr && *value != '\0') {
            char* end = nullptr;
            const auto parsed = std::strtoull(value, &end, 10);
            if (end != value && *end == '\0' && parsed > 0)
                requested_threads = static_cast<std::size_t>(parsed);
        }
        Kokkos::InitializationSettings settings;
        if (requested_threads != 0)
            settings.set_num_threads(resources.effective_threads(requested_threads));
        Kokkos::initialize(settings);
        initialized = true;
        const auto result = run_tool(options, resources);
        Kokkos::finalize();
        return result;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "fastgatk-genotype-gvcf: " << error.what() << '\n';
        return 2;
    }
}
