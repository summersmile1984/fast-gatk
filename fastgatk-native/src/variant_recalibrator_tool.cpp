#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "vqsr_fastmath.hpp"
#include "vqsr_jama_lu.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "optional_boolean.hpp"

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct ResourceSpec {
    std::string label;
    std::string path;
    bool training = false;
    bool truth = false;
    bool known = false;
    // GATK treats a resource prior as a phred-scaled prior quality.  The
    // default in VariantDataManager is Q=2, even when a resource does not
    // explicitly carry a prior label.
    double prior = 2.0;
};

struct Options {
    std::string input;
    std::string output;
    std::string output_model;
    std::string input_model;
    std::string tranches_output;
    std::string manifest;
    std::string rscript_output;
    std::string mode = "SNP";
    std::string reference;
    std::vector<ResourceSpec> resources;
    std::vector<std::string> annotations;
    bool annotations_explicit = false;
    std::vector<double> truth_sensitivity;
    bool output_tranches_for_scatter = false;
    std::vector<double> vqslod_tranches;
    std::size_t sample_every = 1;
    // Match VariantRecalibratorArgumentCollection defaults.  Tests and
    // reproducibility-sensitive callers can still request a single Gaussian
    // explicitly; the native default must not silently change GATK's model
    // complexity.
    std::size_t max_gaussians = 8;
    std::size_t max_negative_gaussians = 2;
    // GATK separates retry attempts from VBEM iterations.  A retry rebuilds
    // a model after a failed initialization, while --max-iterations controls
    // the expectation/maximization loop itself.
    std::size_t max_attempts = 1;
    std::size_t max_iterations = 150;
    std::size_t kmeans_iterations = 100;
    std::size_t maximum_training_variants = 2500000;
    std::size_t minimum_bad_variants = 1000;
    double bad_lod_score_cutoff = -5.0;
    // GATK's MultivariateGaussian always stores and evaluates the full
    // covariance matrix; keep this as the native default.  The command-line
    // switch remains an explicit compatibility spelling for older manifests.
    bool full_covariance = true;
    bool allele_specific = false;
    bool sites_only_vcf_output = false;
    bool create_output_variant_index = true;
    double shrinkage = 1.0;
    double dirichlet = 0.001;
    double prior_counts = 20.0;
    double standard_deviation_threshold = 10.0;
    // GATK's VariantDataManager applies deterministic jitter while decoding
    // annotations.  Keep the two public controls even when the default
    // (MQ_CAP=0, MQ_JITTER=0.05) is used so a native invocation can be replayed
    // with the same command line and random stream as Java GATK.
    int mq_cap = 0;
    double mq_jitter = 0.05;
    int threads = 1;
};

// GATK's VQSR missing-annotation path uses Utils.getRandomGenerator(), which
// is a java.util.Random seeded with GATK_RANDOM_SEED.  Keep the 48-bit LCG and
// cached Box-Muller spare value here so missing-value normalization and the
// twenty-draw marginalization path consume the same sequence on every
// execution space.  This is intentionally Host-only; the stochastic control
// flow must not enter a device kernel.
class JavaRandom {
public:
    explicit JavaRandom(std::uint64_t seed)
        : seed_((seed ^ multiplier_) & mask_) {}

    std::uint32_t next(int bits) {
        seed_ = (seed_ * multiplier_ + addend_) & mask_;
        return static_cast<std::uint32_t>(seed_ >> (48 - bits));
    }

    double next_double() {
        const auto high = static_cast<std::uint64_t>(next(26));
        const auto low = static_cast<std::uint64_t>(next(27));
        return static_cast<double>((high << 27) + low) /
               static_cast<double>(std::uint64_t{1} << 53);
    }

    // java.util.Random.nextBoolean() consumes exactly one 32-bit LCG draw.
    // Keeping this operation explicit matters for VQSR model initialization:
    // MultivariateGaussian.initializeRandomSigma uses one boolean after each
    // lower-triangular covariance draw, and any shortcut here would shift the
    // random stream used by every subsequent Gaussian.
    bool next_boolean() {
        return next(1) != 0;
    }

    // java.util.Random.nextInt(bound), including its rejection step.  This
    // is used by Collections.shuffle when GATK downsamples large training
    // sets, so the native sampler consumes the same 48-bit LCG stream.
    int next_int(int bound) {
        if (bound <= 0) throw std::invalid_argument("random bound must be positive");
        if ((bound & -bound) == bound)
            return static_cast<int>((static_cast<std::int64_t>(bound) * next(31)) >> 31);
        int bits = 0;
        int value = 0;
        do {
            bits = static_cast<int>(next(31));
            value = bits % bound;
        } while (bits - value + (bound - 1) < 0);
        return value;
    }

    double next_gaussian() {
        if (have_spare_) {
            have_spare_ = false;
            return spare_;
        }
        double u = 0.0;
        double v = 0.0;
        double radius = 0.0;
        do {
            u = 2.0 * next_double() - 1.0;
            v = 2.0 * next_double() - 1.0;
            radius = u * u + v * v;
        } while (radius >= 1.0 || radius == 0.0);
        const auto multiplier = std::sqrt(-2.0 * std::log(radius) / radius);
        spare_ = v * multiplier;
        have_spare_ = true;
        return u * multiplier;
    }

private:
    static constexpr std::uint64_t multiplier_ = 0x5DEECE66DULL;
    static constexpr std::uint64_t addend_ = 0xBULL;
    static constexpr std::uint64_t mask_ = (std::uint64_t{1} << 48) - 1;
    std::uint64_t seed_ = 0;
    bool have_spare_ = false;
    double spare_ = 0.0;
};

struct KernelTelemetry {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

struct Entry {
    std::string key;
    std::string type;
    std::vector<double> values;
    std::vector<std::uint8_t> missing;
    bool training = false;
    bool negative_training = false;
    bool truth = false;
    bool known = false;
    bool scoreable = false;
    bool bad_candidate = false;
    // VariantDataManager marks normalized rows beyond STD_THRESHOLD before
    // selecting positive/anti-training data.  The rows remain in the output
    // VCF, but must not participate in either fitted model or worst-row
    // selection.
    bool failing_std_threshold = false;
    bool is_snp = false;
    bool is_transition = false;
    // Stored in the same log10-odds space that GATK adds to the contrastive
    // positive-vs-negative model score.
    double prior_lod = 0.0;
    double score = std::numeric_limits<double>::quiet_NaN();
    // Java VariantDatum.worstAnnotation: the annotation dimension whose
    // per-dimension log10 likelihood separates the positive and negative
    // models least (VariantRecalibratorEngine.calculateWorstPerformingAnnotation).
    // -1 means every dimension was null, which GATK writes as culprit=NULL.
    int worst_annotation = -1;
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

bool complete_file(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() && path.compare(path.size() - value.size(), value.size(), value) == 0;
}

void parse_resource_attributes(ResourceSpec& resource, const std::string& attributes) {
    std::size_t begin = 0;
    while (begin <= attributes.size()) {
        const auto end = attributes.find(',', begin);
        const auto token = attributes.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
        const auto equals = token.find('=');
        const auto name = token.substr(0, equals);
        const auto value = equals == std::string::npos ? std::string{} : token.substr(equals + 1);
        if (name == "training") resource.training = value == "true";
        else if (name == "truth") resource.truth = value == "true";
        else if (name == "known") resource.known = value == "true";
        else if (name == "prior") {
            resource.prior = std::stod(value);
            if (!std::isfinite(resource.prior) || resource.prior < 0.0)
                throw std::invalid_argument("resource prior must be a finite non-negative phred quality");
        }
        if (end == std::string::npos) break;
        begin = end + 1;
    }
}

// QualityUtils.qualToProb(Q) followed by log10(p/(1-p)).  This is deliberately
// kept at the scoring boundary: the Gaussian models are trained only on the
// annotation vectors, while resource priors shift the final contrastive LOD.
double prior_quality_to_log10_odds(double quality) {
    if (!std::isfinite(quality) || quality < 0.0)
        throw std::invalid_argument("resource prior must be a finite non-negative phred quality");
    const double error_probability = std::pow(10.0, -quality / 10.0);
    const double true_probability = 1.0 - error_probability;
    if (!(true_probability > 0.0)) return -std::numeric_limits<double>::infinity();
    if (!(error_probability > 0.0)) return std::numeric_limits<double>::infinity();
    return std::log10(true_probability) - std::log10(error_probability);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-variant-recalibrator (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE                 input VCF/BCF\n"
                         "      --resource[:labels] FILE       training/truth/known VCF (repeatable)\n"
                         "  -an, --use-annotation FIELD        numeric INFO annotation (repeatable)\n"
                         "      --mode SNP|INDEL|BOTH           recalibration mode\n"
                         "  -O, --output FILE                  recal VCF with INFO/VQSLOD\n"
                         "      --tranches-file FILE           output tranches CSV\n"
                         "      --truth-sensitivity-tranche F   target sensitivity percent\n"
                         "      --output-tranches-for-scatter   write VQSLOD tranche slices\n"
                         "      --vqslod-tranche F              VQSLOD slice (repeatable)\n"
                         "      --sample-every-Nth-variant N    deterministic sampling\n"
                         "      --maximum-training-variants N   cap training rows (default 2500000)\n"
                         "      --minimum-bad-variants N        minimum fallback bad rows (default 1000)\n"
                         "      --bad-lod-score-cutoff F        select bad rows at/below LOD (default -5)\n"
                         "      --max-gaussians N               full-covariance GMM components (default 8)\n"
                         "      --max-negative-gaussians N     maximum negative-model components\n"
                         "      --max-attempts N                model-build retries (default 1)\n"
                         "      --max-iterations N              VBEM iterations (default 150)\n"
                         "      --k-means-iterations N          deterministic initialization iterations (default 100)\n"
                         "      --shrinkage F                   VBEM mean shrinkage (default 1)\n"
                         "      --dirichlet F                   VBEM mixture prior (default 0.001)\n"
                         "      --prior-counts F                VBEM covariance prior counts (default 20)\n"
                         "      --standard-deviation-threshold F outlier training threshold\n"
                         "      --mq-cap N                      MQ logit cap/jitter transform (default 0)\n"
                         "      --mq-jitter F                   MQ/AS_MQ jitter scale (default 0.05)\n"
                         "      --full-covariance               use regularized full-covariance GMM\n"
                         "      --output-model FILE             write reusable GATKReport model\n"
                         "      --input-model FILE              reuse a serialized GATKReport model\n"
                         "      --sites-only-vcf-output BOOL   omit FORMAT/sample columns in VCF output\n"
                         "      --create-output-variant-index  optional VCF index (default true)\n"
                         "      --threads N                     Kokkos execution threads\n"
                         "      --output-manifest FILE          OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-V" || has_option(argument, "--variant")) {
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        } else if (argument.rfind("--resource:", 0) == 0) {
            if (index + 1 >= argc) throw std::invalid_argument("missing value for --resource");
            ResourceSpec resource;
            resource.label = argument.substr(std::string("--resource:").size());
            parse_resource_attributes(resource, resource.label);
            resource.path = argv[++index];
            options.resources.push_back(std::move(resource));
        } else if (argument == "--resource" || has_option(argument, "--resource")) {
            ResourceSpec resource;
            resource.label = "resource";
            resource.training = true;
            resource.path = require_value(index, argc, argv, argument, "--resource");
            options.resources.push_back(std::move(resource));
        } else if (argument == "-an" || has_option(argument, "--use-annotation")) {
            options.annotations_explicit = true;
            options.annotations.push_back(require_value(index, argc, argv, argument, "--use-annotation", "-an"));
        } else if (has_option(argument, "--truth-sensitivity-tranche") || has_option(argument, "--tranche")) {
            const auto name = argument.rfind("--tranche", 0) == 0 ? "--tranche" : "--truth-sensitivity-tranche";
            options.truth_sensitivity.push_back(std::stod(require_value(index, argc, argv, argument, name)));
        } else if (has_option(argument, "--mode")) {
            options.mode = require_value(index, argc, argv, argument, "--mode");
            if (options.mode != "SNP" && options.mode != "INDEL" && options.mode != "BOTH")
                throw std::invalid_argument("--mode must be SNP, INDEL, or BOTH");
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--tranches-file")) {
            options.tranches_output = require_value(index, argc, argv, argument, "--tranches-file");
        } else if (has_option(argument, "--vqslod-tranche")) {
            const auto value = std::stod(require_value(index, argc, argv, argument, "--vqslod-tranche"));
            if (!std::isfinite(value)) throw std::invalid_argument("vqslod-tranche must be finite");
            options.vqslod_tranches.push_back(value);
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (has_option(argument, "--rscript-file")) {
            options.rscript_output = require_value(index, argc, argv, argument, "--rscript-file");
        } else if (has_option(argument, "--output-model")) {
            options.output_model = require_value(index, argc, argv, argument, "--output-model");
        } else if (has_option(argument, "--input-model")) {
            options.input_model = require_value(index, argc, argv, argument, "--input-model");
        } else if (has_option(argument, "--sample-every-Nth-variant")) {
            options.sample_every = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--sample-every-Nth-variant")));
            if (options.sample_every == 0) throw std::invalid_argument("sample-every-Nth-variant must be positive");
        } else if (has_option(argument, "--maximum-training-variants") ||
                   has_option(argument, "--max-num-training-data")) {
            const auto name = has_option(argument, "--max-num-training-data")
                ? "--max-num-training-data" : "--maximum-training-variants";
            options.maximum_training_variants = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, name)));
            if (options.maximum_training_variants == 0)
                throw std::invalid_argument("maximum-training-variants must be positive");
        } else if (has_option(argument, "--max-gaussians")) {
            options.max_gaussians = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--max-gaussians")));
            if (options.max_gaussians == 0) throw std::invalid_argument("max-gaussians must be positive");
        } else if (has_option(argument, "--minimum-bad-variants")) {
            options.minimum_bad_variants = static_cast<std::size_t>(std::stoull(
                require_value(index, argc, argv, argument, "--minimum-bad-variants")));
            if (options.minimum_bad_variants == 0)
                throw std::invalid_argument("minimum-bad-variants must be positive");
        } else if (has_option(argument, "--bad-lod-score-cutoff") ||
                   has_option(argument, "--bad-lod-cutoff")) {
            const auto name = has_option(argument, "--bad-lod-cutoff")
                ? "--bad-lod-cutoff" : "--bad-lod-score-cutoff";
            options.bad_lod_score_cutoff = std::stod(require_value(index, argc, argv, argument, name));
            if (!std::isfinite(options.bad_lod_score_cutoff))
                throw std::invalid_argument("bad-lod-score-cutoff must be finite");
        } else if (has_option(argument, "--max-negative-gaussians")) {
            options.max_negative_gaussians = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--max-negative-gaussians")));
            if (options.max_negative_gaussians == 0)
                throw std::invalid_argument("max-negative-gaussians must be positive");
        } else if (has_option(argument, "--max-attempts")) {
            options.max_attempts = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--max-attempts")));
            if (options.max_attempts == 0) throw std::invalid_argument("max-attempts must be positive");
        } else if (has_option(argument, "--max-iterations")) {
            options.max_iterations = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--max-iterations")));
            if (options.max_iterations == 0) throw std::invalid_argument("max-iterations must be positive");
        } else if (has_option(argument, "--k-means-iterations")) {
            options.kmeans_iterations = static_cast<std::size_t>(std::stoull(require_value(
                index, argc, argv, argument, "--k-means-iterations")));
            if (options.kmeans_iterations == 0)
                throw std::invalid_argument("k-means-iterations must be positive");
        } else if (has_option(argument, "--shrinkage")) {
            options.shrinkage = std::stod(require_value(index, argc, argv, argument, "--shrinkage"));
            if (!std::isfinite(options.shrinkage) || options.shrinkage < 0.0)
                throw std::invalid_argument("shrinkage must be finite and non-negative");
        } else if (has_option(argument, "--dirichlet")) {
            options.dirichlet = std::stod(require_value(index, argc, argv, argument, "--dirichlet"));
            if (!std::isfinite(options.dirichlet) || options.dirichlet <= 0.0)
                throw std::invalid_argument("dirichlet must be finite and positive");
        } else if (has_option(argument, "--prior-counts")) {
            options.prior_counts = std::stod(require_value(index, argc, argv, argument, "--prior-counts"));
            if (!std::isfinite(options.prior_counts) || options.prior_counts < 0.0)
                throw std::invalid_argument("prior-counts must be finite and non-negative");
        } else if (has_option(argument, "--standard-deviation-threshold") || has_option(argument, "--std")) {
            const auto name = has_option(argument, "--std") ? "--std" : "--standard-deviation-threshold";
            options.standard_deviation_threshold = std::stod(require_value(index, argc, argv, argument, name));
            if (!std::isfinite(options.standard_deviation_threshold) ||
                options.standard_deviation_threshold <= 0.0)
                throw std::invalid_argument("standard-deviation-threshold must be finite and positive");
        } else if (has_option(argument, "--mq-cap") || has_option(argument, "--mq-cap-for-logit-jitter-transform")) {
            const auto name = argument.rfind("--mq-cap-for-logit-jitter-transform", 0) == 0
                ? "--mq-cap-for-logit-jitter-transform" : "--mq-cap";
            options.mq_cap = std::stoi(require_value(index, argc, argv, argument, name));
            if (options.mq_cap < 0) throw std::invalid_argument("mq-cap must be non-negative");
        } else if (has_option(argument, "--mq-jitter")) {
            options.mq_jitter = std::stod(require_value(index, argc, argv, argument, "--mq-jitter"));
            if (!std::isfinite(options.mq_jitter) || options.mq_jitter < 0.0)
                throw std::invalid_argument("mq-jitter must be finite and non-negative");
        } else if (argument == "--full-covariance" ||
                   argument.rfind("--full-covariance=", 0) == 0) {
            options.full_covariance = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--full-covariance");
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
            if (options.threads < 1) throw std::invalid_argument("threads must be positive");
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "--AS" || argument.rfind("--AS=", 0) == 0 ||
                   argument == "--use-allele-specific-annotations" ||
                   argument.rfind("--use-allele-specific-annotations=", 0) == 0) {
            const char* name = argument.rfind("--use-allele-specific-annotations", 0) == 0
                ? "--use-allele-specific-annotations" : "--AS";
            options.allele_specific = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (argument == "--output-tranches-for-scatter" ||
                   argument.rfind("--output-tranches-for-scatter=", 0) == 0) {
            options.output_tranches_for_scatter = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--output-tranches-for-scatter");
        } else if (argument == "--create-output-variant-index" ||
                   argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_output_variant_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--sites-only-vcf-output" ||
                   argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--sites-only-vcf-output");
        } else if (argument == "--quiet" || argument.rfind("--quiet=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--quiet");
        } else if (argument == "--ignore-all-filters" ||
                   argument.rfind("--ignore-all-filters=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--ignore-all-filters");
        } else if (argument == "--dont-run-rscript" ||
                   argument.rfind("--dont-run-rscript=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--dont-run-rscript");
        } else if (argument == "--trust-all-polymorphic" ||
                   argument.rfind("--trust-all-polymorphic=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--trust-all-polymorphic");
        } else if (argument == "--disable-sequence-dictionary-validation" ||
                   argument.rfind("--disable-sequence-dictionary-validation=", 0) == 0) {
            (void)fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--disable-sequence-dictionary-validation");
        } else if (has_option(argument, "--ignore-filter") || has_option(argument, "--aggregate") ||
                   has_option(argument, "--target-titv")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument, argument.substr(0, argument.find('=' )).c_str());
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.tranches_output.empty()) throw std::invalid_argument("--tranches-file is required");
    if (options.annotations.empty()) options.annotations = {"QD", "MQ", "FS", "SOR"};
    if (options.truth_sensitivity.empty()) options.truth_sensitivity = {100.0, 99.9, 99.0, 90.0};
    if (options.vqslod_tranches.empty()) {
        // Keep the exact resolution used by VariantRecalibrator's hidden
        // VQSLOD_TRANCHES list: coarse above 5, fine around zero, and coarse
        // again below -5.  Explicit --vqslod-tranche values are intended for
        // tests and scatter plans with a smaller transfer footprint.
        for (double value = 10.0; value > 5.0; value -= 0.1)
            options.vqslod_tranches.push_back(value);
        for (double value = 5.0; value > -5.0; value -= 0.01)
            options.vqslod_tranches.push_back(value);
        for (double value = -5.0; value > -10.0; value -= 0.1)
            options.vqslod_tranches.push_back(value);
    }
    for (const auto target : options.truth_sensitivity)
        if (!std::isfinite(target) || target < 0.0 || target > 100.0)
            throw std::invalid_argument("truth-sensitivity-tranche must be in [0,100]");
    return options;
}

#if FASTGATK_HAS_HTSLIB

bool equals_ignore_case(std::string left, std::string right) {
    if (left.size() != right.size()) return false;
    std::transform(left.begin(), left.end(), left.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    std::transform(right.begin(), right.end(), right.begin(),
                   [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
    return left == right;
}

// Mirror VariantDataManager.decodeAnnotation's deterministic jitter boundary.
// The Java implementation compares values with a 0.01 tolerance rather than
// requiring an exact decimal token, so use the same inclusive comparison here.
double decode_annotation_jitter(std::string_view annotation, double value,
                                JavaRandom* random, int mq_cap, double mq_jitter) {
    if (random == nullptr || !std::isfinite(value)) return value;
    const auto near = [](double lhs, double rhs) { return std::abs(lhs - rhs) <= 0.01; };
    const auto is = [&](const char* name) { return equals_ignore_case(std::string(annotation), name); };
    if ((is("HaplotypeScore") || is("FS") || is("AS_FS") || is("AS_FilterStatus") ||
         is("InbreedingCoeff")) && near(value, 0.0)) {
        return value + 0.01 * random->next_gaussian();
    }
    // SOR is encoded as ln(2) at its lower bound by GATK's annotation engine.
    if ((is("SOR") || is("AS_SOR") || is("AS_StrandOddsRatio")) &&
        near(value, 0.6931472)) {
        return value + 0.01 * random->next_gaussian();
    }
    if (is("MQ")) {
        if (mq_cap > 0) {
            const double transformed = std::log((value + 0.01) /
                                                 (static_cast<double>(mq_cap) + 0.01 - value));
            const double cap_transformed = std::log((static_cast<double>(mq_cap) + 0.01) /
                                                    0.01);
            return near(transformed, cap_transformed)
                ? transformed + mq_jitter * random->next_gaussian() : transformed;
        }
        if (near(value, static_cast<double>(mq_cap)))
            return value + mq_jitter * random->next_gaussian();
    }
    if (is("AS_MQ") || is("AS_RMS_MAPPING_QUALITY"))
        return value + mq_jitter * random->next_gaussian();
    return value;
}

std::string variant_type(bcf1_t* record) {
    const auto types = bcf_get_variant_types(record);
    const bool snp = (types & VCF_SNP) != 0;
    const bool mnp = (types & VCF_MNP) != 0;
    const bool indel = (types & VCF_INDEL) != 0;
    if ((snp || mnp) && indel) return "MIXED";
    if (indel) return "INDEL";
    if (snp) return "SNP";
    if (mnp) return "MNP";
    bcf_unpack(record, BCF_UN_STR);
    if (record->d.allele != nullptr)
        for (int allele = 1; allele < record->n_allele; ++allele)
            if (record->d.allele[allele] != nullptr && record->d.allele[allele][0] == '<')
                return "SYMBOLIC";
    return "OTHER";
}

bool is_transition_base(char reference, char alternate) {
    const char ref = static_cast<char>(std::toupper(static_cast<unsigned char>(reference)));
    const char alt = static_cast<char>(std::toupper(static_cast<unsigned char>(alternate)));
    return (ref == 'A' && alt == 'G') || (ref == 'G' && alt == 'A') ||
           (ref == 'C' && alt == 'T') || (ref == 'T' && alt == 'C');
}

bool mode_matches_variant(const std::string& type, const std::string& mode) {
    // VariantDataManager.checkVariationClass treats SNP mode as SNP or MNP,
    // and INDEL mode as indel, mixed, structural/symbolic.  Comparing the
    // collapsed HTSlib type string directly would silently drop MNPs and
    // mixed/symbolic records from otherwise valid GATK runs.
    if (mode == "BOTH") return true;
    if (mode == "SNP") return type == "SNP" || type == "MNP";
    if (mode == "INDEL") return type == "INDEL" || type == "MIXED" || type == "SYMBOLIC";
    return false;
}

std::vector<std::string> record_keys(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    std::vector<std::string> keys;
    const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
    if (!contig || !record->d.allele || record->n_allele < 2) return keys;
    for (int alt = 1; alt < record->n_allele; ++alt) {
        if (!record->d.allele[alt]) continue;
        keys.push_back(std::string(contig) + ":" + std::to_string(record->pos + 1) + ":" +
                       record->d.allele[0] + ":" + record->d.allele[alt]);
    }
    return keys;
}

bool resource_record_is_valid(const bcf_hdr_t* header, bcf1_t* record) {
    if (!header || !record) return false;
    bcf_unpack(record, BCF_UN_ALL);
    // VariantDataManager.isValidVariant requires an unfiltered resource
    // record.  HTSlib keeps PASS as an empty FILTER vector (and may expose a
    // literal PASS id in hand-written VCFs), so accept both representations
    // while rejecting every other filter label.
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name && std::string(name) != "PASS" && std::string(name) != ".") return false;
    }
    if (!record->d.allele || record->n_allele < 2 || !record->d.allele[0]) return false;
    for (int alt = 1; alt < record->n_allele; ++alt) {
        const char* alternate = record->d.allele[alt];
        if (!alternate || alternate[0] == '\0' || std::string(alternate) == "<NON_REF>") return false;
    }
    return true;
}

std::unordered_set<std::string> read_resource_keys(const std::string& path) {
    htsFile* input = bcf_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open resource VCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    bcf1_t* record = bcf_init();
    if (!header || !record) throw std::runtime_error("BAD_INPUT: cannot read resource header");
    std::unordered_set<std::string> keys;
    while (bcf_read(input, header, record) == 0) {
        if (resource_record_is_valid(header, record))
            for (const auto& key : record_keys(header, record)) keys.insert(key);
        bcf_clear(record);
    }
    bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
    if (keys.empty()) throw std::runtime_error("BAD_INPUT: resource contains no variant keys: " + path);
    return keys;
}

bool annotation_values(const bcf_hdr_t* header, bcf1_t* record,
                       const std::vector<std::string>& names, std::vector<double>& values,
                       std::vector<std::uint8_t>& missing,
                       int allele_index = 0, bool allele_specific = false,
                       JavaRandom* random = nullptr, int mq_cap = 0,
                       double mq_jitter = 0.05) {
    values.clear();
    missing.clear();
    for (const auto& name : names) {
        float* buffer = nullptr;
        int count = 0;
        const int length = bcf_get_info_float(header, record, name.c_str(), &buffer, &count);
        int selected = allele_specific ? allele_index : 0;
        if (allele_specific) {
            const int field_id = bcf_hdr_id2int(header, BCF_DT_ID, name.c_str());
            const auto number = field_id >= 0 ? bcf_hdr_id2length(header, BCF_HL_INFO, field_id) : BCF_VL_FIXED;
            if (number != BCF_VL_A && number != BCF_VL_R) {
                free(buffer);
                values.clear();
                missing.clear();
                return false;
            }
            if (number == BCF_VL_R) ++selected; // Number=R includes the REF value.
        }
        const bool missing_value = length <= selected || buffer == nullptr || !std::isfinite(buffer[selected]);
        values.push_back(missing_value ? 0.0 : decode_annotation_jitter(
            name, static_cast<double>(buffer[selected]), random, mq_cap, mq_jitter));
        missing.push_back(missing_value ? 1 : 0);
        free(buffer);
    }
    return values.size() == names.size();
}

std::unordered_map<std::string, std::unordered_set<std::string>> build_resource_sets(
    const std::vector<ResourceSpec>& resources) {
    std::unordered_map<std::string, std::unordered_set<std::string>> sets;
    for (const auto& resource : resources)
        sets.emplace(resource.path, read_resource_keys(resource.path));
    return sets;
}

#endif

struct GaussianMixture {
    std::vector<double> weights;
    std::vector<std::vector<double>> means;
    std::vector<std::vector<double>> variances;
    std::vector<std::vector<double>> covariances;
    std::vector<std::vector<double>> precision;
    std::vector<double> log_determinants;
    // Java MultivariateGaussian.precomputeDenominatorForEvaluation caches
    // log10(pow(2*PI,-N/2)) + log10(pow(det,-0.5)) per component for final
    // (non-variational) evaluation; keep that exact log10 composition here.
    std::vector<double> eval_log10_denominators;
    // Variational-Bayes state used only during training.  The serialized
    // GATKReport intentionally stores the final Gaussian parameters; these
    // posterior hyperparameters are rebuilt for a fresh training run and are
    // not needed by the final scoring kernel.
    std::vector<double> hyper_a;
    std::vector<double> hyper_b;
    std::vector<double> hyper_lambda;
    std::vector<double> variational_log_denom;
    std::vector<std::vector<double>> variational_precision;
    bool variational_ready = false;
    std::size_t count = 0;
    bool full_covariance = false;
    // Training and scoring policies are explicit so a serialized model can
    // distinguish model construction from the Kokkos scoring path.
    std::string training_execution_space;
    std::string training_execution_policy;
    std::string scoring_execution_space;
    std::string scoring_execution_policy;
};

// Apache Commons Math 3.5's Gamma.digamma implementation (used by GATK VQSR)
// is ported bit-exactly in vqsr_fastmath.hpp (including the FastMath.log
// table path it delegates to for arguments >= 49); gatk_digamma_exact keeps
// the recursive recurrence association and the inv/252.0 division.  The
// gatk_digamma name is kept as the single call-site spelling.
double gatk_digamma(double value) {
    return fastgatk::vqsr_fastmath::gatk_digamma_exact(value);
}

struct ModelBundle {
    GaussianMixture good;
    GaussianMixture bad;
    std::vector<std::string> annotations;
    std::vector<double> annotation_means;
    std::vector<double> annotation_stdevs;
};

// Jama LU inverse + determinant wrapper.  GATK 4.6.2.0's MultivariateGaussian
// uses Jama.Matrix.inverse()/det() (Jama 1.0.3), so every VQSR precision and
// determinant value must come from the bit-exact port in vqsr_jama_lu.hpp
// rather than from a Gauss-Jordan elimination whose pivoting and reduction
// order differs from Java.
bool invert_jama_lu(const std::vector<double>& matrix, std::size_t dimensions,
                    std::vector<double>& inverse, double& determinant) {
    return fastgatk::vqsr_jama::jama_lu_inverse_det(matrix, dimensions, inverse, determinant);
}

void refresh_precision(GaussianMixture& mixture) {
    const auto dimensions = mixture.means.empty() ? 0 : mixture.means.front().size();
    mixture.precision.assign(mixture.weights.size(), std::vector<double>(dimensions * dimensions, 0.0));
    mixture.log_determinants.assign(mixture.weights.size(), 0.0);
    mixture.eval_log10_denominators.assign(mixture.weights.size(), 0.0);
    // Java precomputeDenominatorForEvaluation: cachedDenomLog10 =
    // log10(pow(2*PI,-N/2)) + log10(pow(det,-0.5)), with Math.pow on the raw
    // LU product determinant.
    const double pi = 3.14159265358979323846;
    const double gaussian_constant_log10 =
        std::log10(std::pow(2.0 * pi, -static_cast<double>(dimensions) / 2.0));
    for (std::size_t component = 0; component < mixture.weights.size(); ++component) {
        std::vector<double> covariance(dimensions * dimensions, 0.0);
        if (mixture.full_covariance && component < mixture.covariances.size()) {
            covariance = mixture.covariances[component];
        } else {
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                covariance[dimension * dimensions + dimension] =
                    std::max(1.0e-6, mixture.variances[component][dimension]);
        }
        std::vector<double> inverse;
        double log_determinant = 0.0;
        double raw_determinant = 0.0;
        bool valid = false;
        {
            if (invert_jama_lu(covariance, dimensions, inverse, raw_determinant) &&
                std::isfinite(raw_determinant) && raw_determinant > 0.0) {
                // Java logs sigma.det() AFTER the LU product (Math.log(det)),
                // not as a sum of per-pivot logs.
                log_determinant = std::log(raw_determinant);
                valid = true;
            }
        }
        if (!valid) {
            covariance.assign(dimensions * dimensions, 0.0);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                covariance[dimension * dimensions + dimension] =
                    std::max(1.0e-4, mixture.variances[component][dimension]);
            if (!invert_jama_lu(covariance, dimensions, inverse, raw_determinant) ||
                !std::isfinite(raw_determinant) || !(raw_determinant > 0.0))
                throw std::runtime_error("MODEL_FAILURE: covariance matrix is singular");
            log_determinant = std::log(raw_determinant);
            if (mixture.full_covariance) mixture.covariances[component] = covariance;
        }
        mixture.precision[component] = std::move(inverse);
        mixture.log_determinants[component] = log_determinant;
        mixture.eval_log10_denominators[component] =
            gaussian_constant_log10 +
            std::log10(std::pow(raw_determinant, -0.5));
    }
    // Any covariance refresh invalidates the Normal-Wishart cache.  Training
    // explicitly rebuilds it after each M-step; serialized/final scoring uses
    // the ordinary Gaussian precision path.
    mixture.variational_ready = false;
}

void refresh_variational_cache(GaussianMixture& mixture) {
    const auto components = mixture.weights.size();
    const auto dimensions = mixture.means.empty() ? 0 : mixture.means.front().size();
    if (components == 0 || dimensions == 0 || mixture.precision.size() != components) {
        mixture.variational_ready = false;
        return;
    }
    if (mixture.hyper_a.size() != components || mixture.hyper_b.size() != components ||
        mixture.hyper_lambda.size() != components) {
        mixture.variational_ready = false;
        return;
    }
    const double minimum_a = std::max(1.0, static_cast<double>(dimensions) + 1.0e-9);
    const double total_lambda = std::accumulate(
        mixture.hyper_lambda.begin(), mixture.hyper_lambda.end(), 0.0);
    if (!(total_lambda > 0.0) || !std::isfinite(total_lambda)) {
        mixture.variational_ready = false;
        return;
    }
    mixture.variational_precision.assign(
        components, std::vector<double>(dimensions * dimensions, 0.0));
    mixture.variational_log_denom.assign(components, 0.0);
    const double log_two = std::log(2.0);
    // Java precomputeDenominatorForVariationalBayes divides EACH natural-log
    // term by Math.log(10.0) separately and then sums in log10 space:
    // cachedDenomLog10 = (pi/log10) + (lambda/log10) + (beta/log10).
    const double log_ten = std::log(10.0);
    for (std::size_t component = 0; component < components; ++component) {
        const double a = std::max(minimum_a, mixture.hyper_a[component]);
        const double b = std::max(1.0e-12, mixture.hyper_b[component]);
        const double lambda = std::max(1.0e-12, mixture.hyper_lambda[component]);
        double wishart_log_term = 0.0;
        for (std::size_t dimension = 1; dimension <= dimensions; ++dimension)
            wishart_log_term += gatk_digamma(
                (a + 1.0 - static_cast<double>(dimension)) / 2.0);
        const double lambda_term =
            0.5 * (wishart_log_term - mixture.log_determinants[component] +
                   static_cast<double>(dimensions) * log_two);
        const double pi_term = gatk_digamma(lambda) - gatk_digamma(total_lambda);
        const double beta_term =
            -static_cast<double>(dimensions) / (2.0 * b);
        const double denom_log10 =
            (pi_term / log_ten) + (lambda_term / log_ten) + (beta_term / log_ten);
        if (!std::isfinite(denom_log10)) {
            mixture.variational_ready = false;
            return;
        }
        mixture.variational_log_denom[component] = denom_log10;
        for (std::size_t index = 0; index < dimensions * dimensions; ++index)
            mixture.variational_precision[component][index] =
                a * mixture.precision[component][index];
    }
    mixture.variational_ready = true;
}

std::vector<std::string> split_whitespace(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> fields;
    for (std::string field; input >> field;) fields.push_back(std::move(field));
    return fields;
}

struct ReportTable {
    std::string name;
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;
};

template <typename Fn>
void for_each_report_line(const std::string& path, Fn&& callback) {
#if FASTGATK_HAS_HTSLIB
    if (suffix(path, ".gz")) {
        htsFile* input = hts_open(path.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open compressed VQSR model report: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (hts_getline(input, '\n', &line) >= 0)
                callback(line.s == nullptr ? std::string{} : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(input);
            throw;
        }
        const auto close_status = hts_close(input);
        free(line.s);
        if (close_status != 0)
            throw std::runtime_error("BAD_INPUT: cannot read compressed VQSR model report: " + path);
        return;
    }
#else
    if (suffix(path, ".gz"))
        throw std::runtime_error("BACKEND_UNAVAILABLE: compressed VQSR model reports require HTSlib");
#endif
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open VQSR model report: " + path);
    std::string line;
    while (std::getline(input, line)) callback(line);
    if (input.bad()) throw std::runtime_error("BAD_INPUT: cannot read VQSR model report: " + path);
}

std::unordered_map<std::string, ReportTable> read_gatk_report(const std::string& path) {
    std::unordered_map<std::string, ReportTable> tables;
    ReportTable current;
    bool in_table = false;
    for_each_report_line(path, [&](const std::string& line) {
        if (line.rfind("#:GATKTable:", 0) == 0) {
            const auto payload = line.substr(std::string("#:GATKTable:").size());
            // The metadata line has a numeric first token.  The following
            // title line has the table name and description.
            if (!payload.empty() && std::isdigit(static_cast<unsigned char>(payload.front()))) return;
            const auto separator = payload.find(':');
            const auto name = payload.substr(0, separator);
            if (name.empty()) return;
            if (in_table) tables[current.name] = std::move(current);
            current = ReportTable{};
            current.name = name;
            in_table = true;
            return;
        }
        if (!in_table || line.empty() || line[0] == '#') return;
        const auto fields = split_whitespace(line);
        if (fields.empty()) return;
        if (current.columns.empty()) current.columns = fields;
        else current.rows.push_back(fields);
    });
    if (in_table) tables[current.name] = std::move(current);
    if (tables.empty()) throw std::runtime_error("BAD_INPUT: VQSR model has no GATKReport tables");
    return tables;
}

std::size_t column_index(const ReportTable& table, const std::string& name) {
    for (std::size_t index = 0; index < table.columns.size(); ++index)
        if (table.columns[index] == name) return index;
    throw std::runtime_error("BAD_INPUT: VQSR model table " + table.name +
                             " is missing column " + name);
}

double report_double(const std::string& value) {
    if (value == "-Inf" || value == "-Infinity" || value == "-inf")
        return -std::numeric_limits<double>::infinity();
    if (value == "Inf" || value == "Infinity" || value == "inf")
        return std::numeric_limits<double>::infinity();
    const auto parsed = std::stod(value);
    if (std::isnan(parsed)) throw std::runtime_error("BAD_INPUT: NaN in VQSR model report");
    return parsed;
}

std::string format_report_double(double value) {
    if (std::isinf(value)) return value < 0.0 ? "-Infinity" : "Infinity";
    if (!std::isfinite(value)) return "NaN";
    std::ostringstream output;
    output << std::scientific << std::setprecision(16) << value;
    return output.str();
}

GaussianMixture read_report_mixture(
    const std::unordered_map<std::string, ReportTable>& tables,
    const std::string& prefix, const std::vector<std::string>& serialized_annotations) {
    const auto& means = tables.at(prefix + "ModelMeans");
    const auto& covariances = tables.at(prefix + "ModelCovariances");
    const auto& pmix = tables.at(prefix == "Positive" ? "GoodGaussianPMix" : "BadGaussianPMix");
    const auto gaussian_column = column_index(means, "Gaussian");
    const std::size_t dimensions = serialized_annotations.size();
    GaussianMixture mixture;
    mixture.full_covariance = true;
    const auto components = means.rows.size();
    if (components == 0 || dimensions == 0)
        throw std::runtime_error("BAD_INPUT: empty VQSR Gaussian model");
    mixture.weights.assign(components, 0.0);
    mixture.means.assign(components, std::vector<double>(dimensions, 0.0));
    mixture.variances.assign(components, std::vector<double>(dimensions, 1.0));
    mixture.covariances.assign(components, std::vector<double>(dimensions * dimensions, 0.0));
    const std::vector<std::size_t> mean_columns = [&] {
        std::vector<std::size_t> indices;
        indices.reserve(dimensions);
        for (const auto& annotation : serialized_annotations)
            indices.push_back(column_index(means, annotation));
        return indices;
    }();
    for (std::size_t row = 0; row < components; ++row) {
        if (means.rows[row].size() < means.columns.size())
            throw std::runtime_error("BAD_INPUT: malformed VQSR means row");
        const auto gaussian = static_cast<std::size_t>(std::stoull(means.rows[row][gaussian_column]));
        if (gaussian >= components) throw std::runtime_error("BAD_INPUT: non-contiguous VQSR Gaussian id");
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            mixture.means[gaussian][dimension] = report_double(means.rows[row][mean_columns[dimension]]);
    }
    const auto pmix_gaussian = column_index(pmix, "Gaussian");
    const auto pmix_value = column_index(pmix, "pMixLog10");
    for (const auto& row : pmix.rows) {
        if (row.size() <= std::max(pmix_gaussian, pmix_value))
            throw std::runtime_error("BAD_INPUT: malformed VQSR PMix row");
        const auto gaussian = static_cast<std::size_t>(std::stoull(row[pmix_gaussian]));
        if (gaussian >= components) throw std::runtime_error("BAD_INPUT: VQSR PMix Gaussian id out of range");
        mixture.weights[gaussian] = std::pow(10.0, report_double(row[pmix_value]));
    }
    double weight_sum = std::accumulate(mixture.weights.begin(), mixture.weights.end(), 0.0);
    if (!(weight_sum > 0.0) || !std::isfinite(weight_sum))
        throw std::runtime_error("BAD_INPUT: VQSR PMix weights are not finite");
    for (auto& weight : mixture.weights) weight /= weight_sum;
    const auto covariance_gaussian = column_index(covariances, "Gaussian");
    const auto covariance_annotation = column_index(covariances, "Annotation");
    for (const auto& row : covariances.rows) {
        if (row.size() < covariances.columns.size())
            throw std::runtime_error("BAD_INPUT: malformed VQSR covariance row");
        const auto gaussian = static_cast<std::size_t>(std::stoull(row[covariance_gaussian]));
        if (gaussian >= components) throw std::runtime_error("BAD_INPUT: VQSR covariance Gaussian id out of range");
        const auto annotation_iterator = std::find(serialized_annotations.begin(), serialized_annotations.end(),
                                                    row[covariance_annotation]);
        if (annotation_iterator == serialized_annotations.end())
            throw std::runtime_error("BAD_INPUT: VQSR covariance annotation is unknown");
        const auto left = static_cast<std::size_t>(annotation_iterator - serialized_annotations.begin());
        for (std::size_t right = 0; right < dimensions; ++right) {
            const auto column = column_index(covariances, serialized_annotations[right]);
            mixture.covariances[gaussian][left * dimensions + right] = report_double(row[column]);
        }
    }
    for (std::size_t component = 0; component < components; ++component)
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            mixture.variances[component][dimension] =
                std::max(1.0e-6, mixture.covariances[component][dimension * dimensions + dimension]);
    mixture.count = components;
    refresh_precision(mixture);
    return mixture;
}

ModelBundle read_model_report(const std::string& path,
                              const std::vector<std::string>& requested_annotations,
                              bool annotations_explicit,
                              bool reorder_to_requested) {
    const auto tables = read_gatk_report(path);
    const auto& means = tables.at("AnnotationMeans");
    const auto& stdevs = tables.at("AnnotationStdevs");
    const auto annotation_column = column_index(means, "Annotation");
    const auto mean_column = column_index(means, "Mean");
    const auto stdev_annotation = column_index(stdevs, "Annotation");
    const auto stdev_column = column_index(stdevs, "StandardDeviation");
    ModelBundle bundle;
    for (const auto& row : means.rows) {
        if (row.size() <= std::max(annotation_column, mean_column))
            throw std::runtime_error("BAD_INPUT: malformed VQSR AnnotationMeans row");
        bundle.annotations.push_back(row[annotation_column]);
        bundle.annotation_means.push_back(report_double(row[mean_column]));
    }
    if (bundle.annotations.empty()) throw std::runtime_error("BAD_INPUT: VQSR model has no annotations");
    bundle.annotation_stdevs.assign(bundle.annotations.size(), 1.0);
    for (const auto& row : stdevs.rows) {
        if (row.size() <= std::max(stdev_annotation, stdev_column)) continue;
        const auto iterator = std::find(bundle.annotations.begin(), bundle.annotations.end(), row[stdev_annotation]);
        if (iterator != bundle.annotations.end())
            bundle.annotation_stdevs[iterator - bundle.annotations.begin()] =
                std::max(1.0e-12, report_double(row[stdev_column]));
    }
    if (annotations_explicit) {
        if (requested_annotations.size() != bundle.annotations.size())
            throw std::runtime_error("BAD_INPUT: requested VQSR annotations do not match input model");
        for (const auto& annotation : requested_annotations)
            if (std::find(bundle.annotations.begin(), bundle.annotations.end(), annotation) == bundle.annotations.end())
                throw std::runtime_error("BAD_INPUT: requested VQSR annotation is absent from input model: " + annotation);
    }
    // Reorder the serialized model into command-line annotation order.  GATK
    // permits a different -an order when it can be mapped by name.
    // When loading a serialized GATK model, the model's AnnotationMeans order
    // is authoritative.  Barclay still validates explicit -an names, but it
    // does not rewrite the persisted model dimensions to the caller's order;
    // doing so would change the serialized model and covariance coordinates.
    const auto desired = (annotations_explicit && reorder_to_requested)
        ? requested_annotations : bundle.annotations;
    std::vector<std::size_t> permutation;
    for (const auto& annotation : desired)
        permutation.push_back(static_cast<std::size_t>(std::find(bundle.annotations.begin(), bundle.annotations.end(), annotation) - bundle.annotations.begin()));
    auto reorder = [&](GaussianMixture& mixture) {
        for (auto& mean : mixture.means) {
            std::vector<double> old = mean;
            for (std::size_t i = 0; i < permutation.size(); ++i) mean[i] = old[permutation[i]];
        }
        for (auto& covariance : mixture.covariances) {
            const auto old = covariance;
            const auto dimensions = permutation.size();
            for (std::size_t i = 0; i < dimensions; ++i)
                for (std::size_t j = 0; j < dimensions; ++j)
                    covariance[i * dimensions + j] = old[permutation[i] * dimensions + permutation[j]];
        }
        mixture.variances.assign(mixture.weights.size(), std::vector<double>(permutation.size(), 1.0));
        for (std::size_t component = 0; component < mixture.weights.size(); ++component)
            for (std::size_t dimension = 0; dimension < permutation.size(); ++dimension)
                mixture.variances[component][dimension] =
                    std::max(1.0e-6, mixture.covariances[component][dimension * permutation.size() + dimension]);
        refresh_precision(mixture);
    };
    bundle.good = read_report_mixture(tables, "Positive", bundle.annotations);
    bundle.bad = read_report_mixture(tables, "Negative", bundle.annotations);
    if (desired != bundle.annotations) {
        reorder(bundle.good);
        reorder(bundle.bad);
        std::vector<double> means, stdevs_reordered;
        for (const auto index : permutation) {
            means.push_back(bundle.annotation_means[index]);
            stdevs_reordered.push_back(bundle.annotation_stdevs[index]);
        }
        bundle.annotation_means = std::move(means);
        bundle.annotation_stdevs = std::move(stdevs_reordered);
        bundle.annotations = desired;
    }
    return bundle;
}

void write_model_table(std::ostream& output, const std::string& name,
                       const std::string& description, std::size_t rows,
                       const std::vector<std::string>& columns,
                       const std::vector<std::vector<std::string>>& values) {
    output << "#:GATKTable:" << columns.size() << ":" << rows << ":";
    for (std::size_t i = 0; i < columns.size(); ++i) {
        const bool string_column = i == 0 || columns[i] == "Annotation";
        output << (string_column ? "%s:" : "%.16E:");
    }
    output << ";\n#:GATKTable:" << name << ":" << description << "\n";
    for (const auto& column : columns) output << std::left << std::setw(28) << column;
    output << '\n';
    for (const auto& row : values) {
        for (const auto& value : row) output << std::left << std::setw(28) << value;
        output << '\n';
    }
    output << '\n';
}

void write_model_report(const std::string& path, const GaussianMixture& good,
                        const GaussianMixture& bad, const std::vector<std::string>& annotations,
                        const std::vector<double>* annotation_means = nullptr,
                        const std::vector<double>* annotation_stdevs = nullptr) {
    std::ostringstream model;
    auto& output = model;
    const auto dimensions = annotations.size();
    output << "#:GATKReport.v1.1:8\n";
    std::vector<std::vector<std::string>> rows;
    for (std::size_t index = 0; index < annotations.size(); ++index)
        rows.push_back({annotations[index], format_report_double(
            annotation_means && index < annotation_means->size() ? (*annotation_means)[index] : 0.0)});
    write_model_table(output, "AnnotationMeans", "Mean for each annotation, used to normalize data",
                      rows.size(), {"Annotation", "Mean"}, rows);
    rows.clear();
    for (std::size_t index = 0; index < annotations.size(); ++index)
        rows.push_back({annotations[index], format_report_double(
            annotation_stdevs && index < annotation_stdevs->size() ? (*annotation_stdevs)[index] : 1.0)});
    write_model_table(output, "AnnotationStdevs", "Standard deviation for each annotation, used to normalize data",
                      rows.size(), {"Annotation", "StandardDeviation"}, rows);
    auto write_pmix = [&](const GaussianMixture& mixture, const std::string& name) {
        rows.clear();
        for (std::size_t component = 0; component < mixture.weights.size(); ++component)
            rows.push_back({std::to_string(component), format_report_double(std::log10(std::max(1.0e-300, mixture.weights[component]))) });
        write_model_table(output, name, "Pmixture log 10 used to evaluate model", rows.size(),
                          {"Gaussian", "pMixLog10"}, rows);
    };
    write_pmix(good, "GoodGaussianPMix");
    auto write_means = [&](const GaussianMixture& mixture, const std::string& name) {
        rows.clear();
        for (std::size_t component = 0; component < mixture.means.size(); ++component) {
            std::vector<std::string> row{std::to_string(component)};
            for (const auto value : mixture.means[component]) row.push_back(format_report_double(value));
            rows.push_back(std::move(row));
        }
        std::vector<std::string> columns{"Gaussian"};
        columns.insert(columns.end(), annotations.begin(), annotations.end());
        write_model_table(output, name, "Vector of annotation values to describe the normalized mean for each Gaussian", rows.size(), columns, rows);
    };
    auto write_covariance = [&](const GaussianMixture& mixture, const std::string& name) {
        rows.clear();
        for (std::size_t component = 0; component < mixture.weights.size(); ++component)
            for (std::size_t left = 0; left < dimensions; ++left) {
                std::vector<std::string> row{std::to_string(component), annotations[left]};
                for (std::size_t right = 0; right < dimensions; ++right)
                    row.push_back(format_report_double(mixture.covariances[component][left * dimensions + right]));
                rows.push_back(std::move(row));
            }
        std::vector<std::string> columns{"Gaussian", "Annotation"};
        columns.insert(columns.end(), annotations.begin(), annotations.end());
        write_model_table(output, name, "Matrix to describe the normalized covariance for each Gaussian", rows.size(), columns, rows);
    };
    write_means(good, "PositiveModelMeans");
    write_covariance(good, "PositiveModelCovariances");
    write_pmix(bad, "BadGaussianPMix");
    write_means(bad, "NegativeModelMeans");
    write_covariance(bad, "NegativeModelCovariances");

    const auto text = model.str();
    if (suffix(path, ".gz")) {
#if FASTGATK_HAS_ZLIB
        gzFile compressed = gzopen(path.c_str(), "wb");
        if (!compressed)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed VQSR model: " + path);
        std::size_t offset = 0;
        while (offset < text.size()) {
            const auto chunk = std::min<std::size_t>(text.size() - offset, 1U << 20);
            if (gzwrite(compressed, text.data() + offset, static_cast<unsigned>(chunk)) !=
                static_cast<int>(chunk)) {
                gzclose(compressed);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed VQSR model: " + path);
            }
            offset += chunk;
        }
        if (gzclose(compressed) != Z_OK)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize compressed VQSR model: " + path);
#else
        throw std::runtime_error("BACKEND_UNAVAILABLE: compressed VQSR model reports require zlib");
#endif
    } else {
        std::ofstream output_file(path);
        if (!output_file)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VQSR model: " + path);
        output_file << text;
        if (!output_file)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write VQSR model: " + path);
    }
    if (!complete_file(path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete VQSR model: " + path);
}

void responsibility_step_kokkos(const std::vector<std::vector<double>>& points,
                                const GaussianMixture& mixture,
                                std::vector<double>& responsibilities,
                                KernelTelemetry* telemetry = nullptr) {
    const auto rows = points.size();
    const auto components = mixture.weights.size();
    const auto dimensions = mixture.means.empty() ? 0 : mixture.means.front().size();
    responsibilities.assign(rows * components, 0.0);
    if (rows == 0 || components == 0 || dimensions == 0) return;
    Kokkos::View<double*> device_points("variant_recalibrator_em_points", rows * dimensions);
    Kokkos::View<double**> device_means("variant_recalibrator_em_means", components, dimensions);
    Kokkos::View<double***> device_precision("variant_recalibrator_em_precision",
                                              components, dimensions, dimensions);
    Kokkos::View<double*> device_weights("variant_recalibrator_em_weights", components);
    Kokkos::View<double*> device_logdet("variant_recalibrator_em_logdet", components);
    Kokkos::View<double*> device_responsibilities("variant_recalibrator_em_responsibilities",
                                                  rows * components);
    // Java evaluateDatumLog10 accumulates a cross-product vector per datum
    // (crossProdTmp[iii] += (x[jjj]-mu[jjj]) * Inv[jjj][iii]) before the final
    // inner product; the scratch view preserves that association order.
    Kokkos::View<double*> device_cross("variant_recalibrator_em_cross",
                                       rows * dimensions);
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host_batch("variant-recalibrator-estep-v1");
    host_batch.records = rows;
    host_batch.bytes = rows * dimensions * sizeof(double) + rows * components * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("variant-recalibrator-estep");
    kernel_plan.begin_prepare(host_batch);
    auto host_points = Kokkos::create_mirror_view(device_points);
    auto host_means = Kokkos::create_mirror_view(device_means);
    auto host_precision = Kokkos::create_mirror_view(device_precision);
    auto host_weights = Kokkos::create_mirror_view(device_weights);
    auto host_logdet = Kokkos::create_mirror_view(device_logdet);
    const bool variational = mixture.variational_ready &&
        mixture.variational_precision.size() == components &&
        mixture.variational_log_denom.size() == components;
    for (std::size_t row = 0; row < rows; ++row)
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            host_points(row * dimensions + dimension) = points[row][dimension];
    for (std::size_t component = 0; component < components; ++component) {
        host_weights(component) = variational ? 1.0 :
            std::max(1.0e-12, mixture.weights[component]);
        // The regular path stores log(det Sigma) and evaluates in natural
        // log.  The variational path stores Java's cachedDenomLog10 (each
        // Normal-Wishart/Dirichlet term divided by log(10) separately) and
        // evaluates log10-likelihoods, matching precomputeDenominatorFor-
        // VariationalBayes/evaluateDatumLog10 bit for bit.
        host_logdet(component) = variational
            ? mixture.variational_log_denom[component]
            : mixture.log_determinants[component];
        for (std::size_t left = 0; left < dimensions; ++left) {
            host_means(component, left) = mixture.means[component][left];
            for (std::size_t right = 0; right < dimensions; ++right)
                host_precision(component, left, right) =
                    (variational ? mixture.variational_precision[component][left * dimensions + right]
                                 : mixture.precision[component][left * dimensions + right]);
        }
    }
    Kokkos::deep_copy(device_points, host_points);
    Kokkos::deep_copy(device_means, host_means);
    Kokkos::deep_copy(device_precision, host_precision);
    Kokkos::deep_copy(device_weights, host_weights);
    Kokkos::deep_copy(device_logdet, host_logdet);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(rows);
    device_batch.bind("points", device_points);
    device_batch.bind("means", device_means);
    device_batch.bind("precision", device_precision);
    device_batch.bind("weights", device_weights);
    device_batch.bind("logdet", device_logdet);
    device_batch.bind("responsibilities", device_responsibilities);
    device_batch.bind("cross", device_cross);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    kernel_plan.begin_execute();
    const double log_ten = std::log(10.0);
    Kokkos::parallel_for("variant_recalibrator_em_estep",
                         Kokkos::RangePolicy<ExecSpace>(0, rows),
                         KOKKOS_LAMBDA(const std::size_t row) {
                             auto datum_eval = [&](const std::size_t component,
                                                   double& quadratic) {
                                 // Zero the Java crossProdTmp accumulator.
                                 for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                                     device_cross(row * dimensions + dimension) = 0.0;
                                 for (std::size_t first = 0; first < dimensions; ++first)
                                     for (std::size_t second = 0; second < dimensions; ++second)
                                         device_cross(row * dimensions + first) +=
                                             (device_points(row * dimensions + second) -
                                              device_means(component, second)) *
                                             device_precision(component, second, first);
                                 quadratic = 0.0;
                                 for (std::size_t first = 0; first < dimensions; ++first)
                                     quadratic += device_cross(row * dimensions + first) *
                                         (device_points(row * dimensions + first) -
                                          device_means(component, first));
                             };
                             double maximum = -1.0e300;
                             for (std::size_t component = 0; component < components; ++component) {
                                 double quadratic = 0.0;
                                 datum_eval(component, quadratic);
                                 double logp;
                                 if (variational) {
                                     // Java evaluateDatumLog10: log10 pdf value
                                     // is ((-0.5*q)/log(10)) + cachedDenomLog10.
                                     logp = ((-0.5 * quadratic) / log_ten) +
                                         device_logdet(component);
                                 } else {
                                     logp = log(device_weights(component)) -
                                         0.5 * (quadratic + device_logdet(component));
                                 }
                                 maximum = maximum > logp ? maximum : logp;
                             }
                             double normalizer = 0.0;
                             for (std::size_t component = 0; component < components; ++component) {
                                 double quadratic = 0.0;
                                 datum_eval(component, quadratic);
                                 double logp;
                                 if (variational) {
                                     logp = ((-0.5 * quadratic) / log_ten) +
                                         device_logdet(component);
                                 } else {
                                     logp = log(device_weights(component)) -
                                         0.5 * (quadratic + device_logdet(component));
                                 }
                                 // Java normalizeLog10DeleteMePlease turns each
                                 // log10 value into 10^(x-max) and normalizes by
                                 // the sum; the natural-log branch keeps exp().
                                 const auto responsibility = variational
                                     ? pow(10.0, logp - maximum)
                                     : exp(logp - maximum);
                                 device_responsibilities(row * components + component) = responsibility;
                                 normalizer += responsibility;
                             }
                             for (std::size_t component = 0; component < components; ++component)
                                 device_responsibilities(row * components + component) /=
                                     normalizer > 1.0e-12 ? normalizer : 1.0e-12;
                         });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) {
        ++telemetry->batches;
        telemetry->observations += rows;
        telemetry->prepare_seconds += kernel_plan.telemetry().prepare_seconds;
        telemetry->execute_seconds += kernel_plan.telemetry().execute_seconds;
    }
    auto host_responsibilities = Kokkos::create_mirror_view(device_responsibilities);
    Kokkos::deep_copy(host_responsibilities, device_responsibilities);
    for (std::size_t index = 0; index < responsibilities.size(); ++index)
        responsibilities[index] = host_responsibilities(index);
}

struct EmMeanStatistics {
    std::vector<double> weights;
    std::vector<double> weighted_sums;
    std::string execution_space;
    std::string execution_policy;
};

// Accumulate the E-step responsibilities with one local table per TeamPolicy
// league entry.  Each team uses a single lane for its bounded tile and walks
// rows/components in a fixed order; this avoids floating-point atomics whose
// order would vary with the backend.  The Host merge is ordered by league and
// then by component/dimension, preserving Strict reproducibility.
EmMeanStatistics em_mean_statistics_kokkos(
    const std::vector<std::vector<double>>& points,
    const std::vector<double>& responsibilities,
    std::size_t components, std::size_t dimensions,
    double weight_base = 0.0,
    KernelTelemetry* telemetry = nullptr) {
    EmMeanStatistics result;
    result.weights.assign(components, 0.0);
    result.weighted_sums.assign(components * dimensions, 0.0);
    result.execution_space = Kokkos::DefaultExecutionSpace::name();
    result.execution_policy = "TeamPolicy";
    if (points.empty() || components == 0 || dimensions == 0) return result;
    if (responsibilities.size() != points.size() * components)
        throw std::invalid_argument("VQSR EM responsibilities have an invalid shape");
    for (const auto& point : points)
        if (point.size() != dimensions)
            throw std::invalid_argument("VQSR EM points have an invalid shape");

    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using MemorySpace = typename ExecSpace::memory_space;
    constexpr std::size_t observations_per_team = 256;
    const auto observation_count = points.size();
    const auto league_size = (observation_count + observations_per_team - 1) /
                             observations_per_team;
    Kokkos::View<double*, MemorySpace> device_points(
        "variant_recalibrator_mstep_points", points.size() * dimensions);
    Kokkos::View<double*, MemorySpace> device_responsibilities(
        "variant_recalibrator_mstep_responsibilities", responsibilities.size());
    Kokkos::View<double**, Kokkos::LayoutRight, MemorySpace> team_weights(
        "variant_recalibrator_mstep_team_weights", league_size, components);
    Kokkos::View<double***, Kokkos::LayoutRight, MemorySpace> team_sums(
        "variant_recalibrator_mstep_team_sums", league_size, components, dimensions);
    fastgatk::core::HostBatch host_batch("variant-recalibrator-mstep-v1");
    host_batch.records = observation_count;
    host_batch.bytes = points.size() * dimensions * sizeof(double) +
                       responsibilities.size() * sizeof(double);
    fastgatk::core::KernelPlan<ExecSpace> kernel_plan("variant-recalibrator-mstep");
    kernel_plan.begin_prepare(host_batch);
    auto host_points = Kokkos::create_mirror_view(device_points);
    auto host_responsibilities = Kokkos::create_mirror_view(device_responsibilities);
    for (std::size_t row = 0; row < points.size(); ++row) {
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            host_points(row * dimensions + dimension) = points[row][dimension];
        for (std::size_t component = 0; component < components; ++component)
            host_responsibilities(row * components + component) =
                responsibilities[row * components + component];
    }
    Kokkos::deep_copy(device_points, host_points);
    Kokkos::deep_copy(device_responsibilities, host_responsibilities);
    Kokkos::deep_copy(team_weights, 0.0);
    Kokkos::deep_copy(team_sums, 0.0);
    // Java MultivariateGaussian.maximizeGaussian initializes sumProb to 1e-10
    // BEFORE accumulating responsibilities per datum, so the first added term
    // is the base and every later add accumulates onto it.  Seeding the first
    // league's accumulator with the base reproduces that exact association
    // when a single league covers the observations; for larger inputs the
    // league merge order remains league-ordered (a documented boundary on the
    // multi-league raw-bit guarantee).
    if (weight_base != 0.0 && league_size > 0) {
        auto host_init = Kokkos::create_mirror_view(team_weights);
        Kokkos::deep_copy(host_init, team_weights);
        for (std::size_t component = 0; component < components; ++component)
            host_init(0, component) = weight_base;
        Kokkos::deep_copy(team_weights, host_init);
    }
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(observation_count);
    device_batch.bind("points", device_points);
    device_batch.bind("responsibilities", device_responsibilities);
    device_batch.bind("team_weights", team_weights);
    device_batch.bind("team_sums", team_sums);
    ExecSpace().fence();
    kernel_plan.end_prepare(device_batch);
    using TeamMember = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
    const Kokkos::TeamPolicy<ExecSpace> policy(league_size, Kokkos::AUTO());
    kernel_plan.begin_execute();
    Kokkos::parallel_for("variant_recalibrator_em_mstep", policy,
        KOKKOS_LAMBDA(const TeamMember& team) {
            const auto league = static_cast<std::size_t>(team.league_rank());
            const auto begin = league * observations_per_team;
            const auto end = (begin + observations_per_team < observation_count)
                ? begin + observations_per_team : observation_count;
            Kokkos::single(Kokkos::PerTeam(team), [&] {
                for (std::size_t row = begin; row < end; ++row)
                    for (std::size_t component = 0; component < components; ++component) {
                        const auto responsibility = device_responsibilities(
                            row * components + component);
                        team_weights(league, component) += responsibility;
                        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                            team_sums(league, component, dimension) +=
                                responsibility * device_points(row * dimensions + dimension);
                    }
            });
        });
    ExecSpace().fence();
    kernel_plan.end_execute();
    if (telemetry != nullptr) {
        ++telemetry->batches;
        telemetry->observations += observation_count;
        telemetry->prepare_seconds += kernel_plan.telemetry().prepare_seconds;
        telemetry->execute_seconds += kernel_plan.telemetry().execute_seconds;
    }
    const auto host_team_weights = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), team_weights);
    const auto host_team_sums = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), team_sums);
    for (std::size_t league = 0; league < league_size; ++league) {
        for (std::size_t component = 0; component < components; ++component) {
            result.weights[component] += host_team_weights(league, component);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                result.weighted_sums[component * dimensions + dimension] +=
                    host_team_sums(league, component, dimension);
        }
    }
    return result;
}

std::vector<std::vector<double>> selected_points(const std::vector<Entry>& entries,
                                                 bool training, std::size_t dimensions,
                                                 double standard_deviation_threshold,
                                                 std::size_t maximum_points,
                                                 JavaRandom* random) {
    std::vector<std::vector<double>> points;
    for (const auto& entry : entries) {
        if (!entry.scoreable || entry.failing_std_threshold || (training && !entry.training) ||
            (!training && !entry.bad_candidate && !entry.negative_training)) continue;
        if (training && standard_deviation_threshold > 0.0 &&
            std::any_of(entry.values.begin(), entry.values.end(), [&](const double value) {
                return std::isfinite(value) && std::abs(value) > standard_deviation_threshold;
            })) continue;
        points.push_back(entry.values);
    }
    // GATK bounds the amount of training data to keep resident with
    // Collections.shuffle(trainingData, Utils.getRandomGenerator()).  Keep
    // that operation Host-side and consume the same Java LCG stream; device
    // kernels still receive a compact, deterministic point matrix.
    if (maximum_points > 0 && points.size() > maximum_points) {
        std::vector<std::size_t> order(points.size());
        std::iota(order.begin(), order.end(), 0);
        if (random != nullptr) {
            for (std::size_t index = order.size(); index > 1; --index) {
                const auto swap_index = static_cast<std::size_t>(random->next_int(
                    static_cast<int>(index)));
                std::swap(order[index - 1], order[swap_index]);
            }
        }
        std::vector<std::vector<double>> sampled;
        sampled.reserve(maximum_points);
        for (std::size_t rank = 0; rank < maximum_points; ++rank)
            sampled.push_back(std::move(points[order[rank]]));
        points.swap(sampled);
    }
    (void)dimensions;
    return points;
}

GaussianMixture fit_mixture(const std::vector<Entry>& entries, bool training,
                            std::size_t dimensions, std::size_t requested_components,
                            std::size_t max_attempts, bool full_covariance,
                            double shrinkage, double dirichlet, double prior_counts,
                            double standard_deviation_threshold,
                            std::size_t kmeans_iterations,
                            std::size_t maximum_training_variants,
                            JavaRandom* random = nullptr,
                            KernelTelemetry* telemetry = nullptr) {
    auto points = selected_points(entries, training, dimensions, standard_deviation_threshold,
                                  maximum_training_variants, random);
    GaussianMixture mixture;
    mixture.count = points.size();
    mixture.full_covariance = full_covariance;
    mixture.training_execution_space = Kokkos::DefaultExecutionSpace::name();
    mixture.training_execution_policy = "TeamPolicy";
    if (points.empty()) return mixture;
    const auto components = std::min(requested_components, points.size());
    mixture.weights.assign(components, 1.0 / static_cast<double>(components));
    mixture.means.assign(components, std::vector<double>(dimensions, 0.0));
    mixture.variances.assign(components, std::vector<double>(dimensions, 1.0));
    mixture.covariances.assign(components, std::vector<double>(dimensions * dimensions, 0.0));

    // Match GaussianMixtureModel.initializeRandomModel.  GATK first consumes
    // one Java-Random draw per annotation for every Gaussian, then runs
    // K-means.  The previous quantile seed was deterministic but it described
    // a different model and could change both the selected local optimum and
    // every downstream VQSLOD for max-gaussians > 1.
    for (std::size_t component = 0; component < components; ++component)
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            mixture.means[component][dimension] = random == nullptr
                ? 0.0
                : -4.0 + 8.0 * random->next_double();

    // GATK uses K-means before VBEM.  Keep assignment and reduction order
    // stable across execution spaces while preserving Java's random restart
    // for an empty cluster.
    for (std::size_t iteration = 0; iteration < kmeans_iterations; ++iteration) {
        std::vector<std::size_t> assignments(points.size(), 0);
        std::vector<std::vector<double>> sums(components, std::vector<double>(dimensions, 0.0));
        std::vector<std::size_t> counts(components, 0);
        for (std::size_t row = 0; row < points.size(); ++row) {
            std::size_t best = 0;
            double best_distance = std::numeric_limits<double>::infinity();
            for (std::size_t component = 0; component < components; ++component) {
                double distance = 0.0;
                for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                    const auto delta = points[row][dimension] - mixture.means[component][dimension];
                    distance += delta * delta;
                }
                if (distance < best_distance || (distance == best_distance && component < best)) {
                    best_distance = distance;
                    best = component;
                }
            }
            assignments[row] = best;
            ++counts[best];
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                sums[best][dimension] += points[row][dimension];
        }
        for (std::size_t component = 0; component < components; ++component) {
            if (counts[component] == 0) {
                for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                    mixture.means[component][dimension] = random == nullptr
                        ? 0.0
                        : -4.0 + 8.0 * random->next_double();
                continue;
            }
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                mixture.means[component][dimension] = sums[component][dimension] /
                    static_cast<double>(counts[component]);
        }
    }

    // initializeRandomSigma builds a random lower-triangular matrix L and
    // stores L*L^T.  Its diagonal remains strictly positive, so this is SPD
    // before the first variational E-step exactly as in the Java model.
    for (std::size_t component = 0; component < components; ++component) {
        std::vector<double> lower(dimensions * dimensions, 0.0);
        // Java MultivariateGaussian.initializeRandomSigma walks the column
        // index first and then the rows below that diagonal (iii outer,
        // jjj=iii..N inner).  The matrix is lower triangular, but walking it
        // row-first consumes the JavaRandom stream in a different order and
        // changes multi-Gaussian local optima even when every later formula
        // is otherwise equivalent.
        for (std::size_t column = 0; column < dimensions; ++column)
            for (std::size_t row = column; row < dimensions; ++row) {
                auto value = random == nullptr
                    ? 1.0
                    : 0.55 + 1.25 * random->next_double();
                if (random != nullptr && random->next_boolean()) value *= -1.0;
                lower[row * dimensions + column] = value;
            }
        for (std::size_t row = 0; row < dimensions; ++row)
            for (std::size_t column = 0; column < dimensions; ++column) {
                double value = 0.0;
                for (std::size_t inner = 0; inner < dimensions; ++inner)
                    value += lower[row * dimensions + inner] * lower[column * dimensions + inner];
                mixture.covariances[component][row * dimensions + column] = value;
                if (row == column) mixture.variances[component][row] = value;
            }
    }
    mixture.hyper_a.assign(components, std::max(1.0, prior_counts));
    mixture.hyper_b.assign(components, std::max(1.0e-12, shrinkage));
    mixture.hyper_lambda.assign(components, std::max(1.0e-12, dirichlet));
    refresh_precision(mixture);
    refresh_variational_cache(mixture);

    std::vector<double> responsibilities(points.size() * components, 0.0);
    // Java's GaussianMixtureModel starts with an expectation step, then each
    // iteration performs M-step -> pMix normalization -> expectation.  Keep
    // the same state transition here; in particular convergence is based on
    // the L1 change in mixture coefficients, not a mean-coordinate delta.
    responsibility_step_kokkos(points, mixture, responsibilities, telemetry);
    for (std::size_t attempt = 0; attempt < max_attempts; ++attempt) {
        const auto previous_weights = mixture.weights;
        const auto statistics = em_mean_statistics_kokkos(
            points, responsibilities, components, dimensions, 1.0e-10, telemetry);
        const double empirical_sigma = 1.0 / 200.0;
        const double total_probability = std::max(1.0e-10,
            std::accumulate(statistics.weights.begin(), statistics.weights.end(), 0.0));
        for (std::size_t component = 0; component < components; ++component) {
            // MultivariateGaussian.maximizeGaussian uses the raw posterior
            // mass for pMix normalization.  The Dirichlet term affects the
            // variational denominator through hyperParameter_lambda; it is
            // not added to the serialized mixture coefficient.  The 1e-10
            // base is already part of the seeded accumulation (matching
            // Java's sumProb = 1E-10 before the per-datum adds).
            const double weight = statistics.weights[component];
            mixture.weights[component] = weight / total_probability;
            // Java MultivariateGaussian.maximizeGaussian updates these three
            // Normal-Wishart/Dirichlet posterior hyperparameters before the
            // next variational E-step.
            mixture.hyper_a[component] = weight + prior_counts;
            mixture.hyper_b[component] = weight + shrinkage;
            mixture.hyper_lambda[component] = weight + dirichlet;
            std::vector<double> raw_mean(dimensions, 0.0);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                raw_mean[dimension] = statistics.weighted_sums[component * dimensions + dimension] / weight;
            const double shrinkage_factor = (shrinkage * weight) / (shrinkage + weight);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                // VariantRecalibrator's empirical mean is zero because the
                // data were normalized before model fitting.
                mixture.means[component][dimension] =
                    (weight * raw_mean[dimension]) / (weight + shrinkage);

            auto& covariance = mixture.covariances[component];
            std::fill(covariance.begin(), covariance.end(), 0.0);
            for (std::size_t row = 0; row < points.size(); ++row) {
                const auto responsibility = responsibilities[row * components + component];
                for (std::size_t left = 0; left < dimensions; ++left)
                    for (std::size_t right = 0; right < dimensions; ++right) {
                        // Java's pVarSigma.set(iii, jjj, prob*(...)*(...) + regCovar)
                        // adds the per-datum regularizer to the product BEFORE
                        // the running plusEquals; keep the same association.
                        covariance[left * dimensions + right] += responsibility *
                            (points[row][left] - raw_mean[left]) *
                            (points[row][right] - raw_mean[right]) +
                            (left == right ? 1.0e-6 : 0.0);
                    }
            }
            for (std::size_t left = 0; left < dimensions; ++left)
                for (std::size_t right = 0; right < dimensions; ++right) {
                    covariance[left * dimensions + right] += empirical_sigma *
                        (left == right ? 1.0 : 0.0);
                    // GATK's Normal-Wishart prior contributes one shrinkage
                    // factor to the outer product.  Applying the factor to
                    // both deltas squares the prior strength and changes
                    // multi-Gaussian responsibilities materially.
                    covariance[left * dimensions + right] += shrinkage_factor *
                        raw_mean[left] * raw_mean[right];
                    covariance[left * dimensions + right] = std::max(
                        left == right ? 1.0e-6 : -std::numeric_limits<double>::infinity(),
                        covariance[left * dimensions + right]);
                }
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                mixture.variances[component][dimension] = std::max(
                    1.0e-6, covariance[dimension * dimensions + dimension]);
        }
        refresh_precision(mixture);
        refresh_variational_cache(mixture);
        // Complete this iteration's E-step before checking convergence, as
        // Java does.  The resulting responsibilities are also the state fed
        // to evaluateFinalModelParameters after the loop.
        responsibility_step_kokkos(points, mixture, responsibilities, telemetry);
        double mixture_change = 0.0;
        for (std::size_t component = 0; component < components; ++component)
            mixture_change += std::abs(mixture.weights[component] -
                (component < previous_weights.size() ? previous_weights[component] : 0.0));
        if (attempt + 1 > 2 && mixture_change < 2.0e-3) break;
    }
    // GATK's VBEM loop already completes an E-step after every M-step.  Its
    // evaluateFinalModelParameters() consumes that last responsibility table
    // directly; it does *not* run another E-step after convergence.  Keeping
    // the final table here is important for multi-component models: an extra
    // E-step can harden a narrow component and materially change its final
    // covariance and downstream VQSLOD, while a one-Gaussian fixture masks the
    // error completely.
    const auto final_statistics = em_mean_statistics_kokkos(
        points, responsibilities, components, dimensions, 0.0, telemetry);
    // Java's evaluateFinalModelParameters() resets sumProb to 0.0 (no 1e-10
    // base) and divides mu by the raw posterior mass; normalizePMixtureLog10
    // then works on those raw masses.  Keep the same un-clamped values.
    const double final_total = std::max(1.0e-10,
        std::accumulate(final_statistics.weights.begin(),
                        final_statistics.weights.end(), 0.0));
    for (std::size_t component = 0; component < components; ++component) {
        // Java's evaluateFinalModelParameters() divides by the raw posterior
        // mass, but a numerically empty component must not poison downstream
        // scoring with NaN means/covariances; the 1e-10 floor keeps the
        // jitter/zero-variance fixtures well-formed while staying far below
        // any mass a real component carries.
        const double weight = std::max(1.0e-10, final_statistics.weights[component]);
        mixture.weights[component] = weight / final_total;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
            mixture.means[component][dimension] =
                final_statistics.weighted_sums[component * dimensions + dimension] / weight;
        if (full_covariance) {
            auto& covariance = mixture.covariances[component];
            std::fill(covariance.begin(), covariance.end(), 0.0);
            for (std::size_t row = 0; row < points.size(); ++row) {
                const auto responsibility = responsibilities[row * components + component];
                for (std::size_t left = 0; left < dimensions; ++left)
                    for (std::size_t right = 0; right < dimensions; ++right) {
                        // Java: pVarSigma.set(iii, jjj, prob*(...)*(...) + regCovar)
                        // then sigma.plusEquals; reg is one term of each datum.
                        covariance[left * dimensions + right] += responsibility *
                            (points[row][left] - mixture.means[component][left]) *
                            (points[row][right] - mixture.means[component][right]) +
                            (left == right ? 1.0e-6 : 0.0);
                    }
            }
            // Java: sigma.timesEquals(1.0 / sumProb) multiplies every element
            // by one reciprocal rather than dividing element-wise.
            const double inverse_weight = 1.0 / weight;
            for (std::size_t left = 0; left < dimensions; ++left)
                for (std::size_t right = 0; right < dimensions; ++right)
                    covariance[left * dimensions + right] *= inverse_weight;
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                mixture.variances[component][dimension] = std::max(
                    1.0e-6, covariance[dimension * dimensions + dimension]);
        } else {
            auto& covariance = mixture.covariances[component];
            std::fill(covariance.begin(), covariance.end(), 0.0);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                double variance = 0.0;
                for (std::size_t row = 0; row < points.size(); ++row) {
                    const auto delta = points[row][dimension] - mixture.means[component][dimension];
                    variance += responsibilities[row * components + component] * delta * delta;
                    variance += 1.0e-6;
                }
                mixture.variances[component][dimension] = std::max(
                    1.0e-6, (variance + 1.0 / 200.0) / weight);
                covariance[dimension * dimensions + dimension] =
                    mixture.variances[component][dimension];
            }
        }
    }
    refresh_precision(mixture);
    // Final scoring follows GATK's non-variational Gaussian evaluator.
    mixture.variational_ready = false;
    mixture.training_execution_policy = "VBEM-TeamPolicy-KMeans";
    return mixture;
}

// Java GaussianMixtureModel.evaluateDatum for a complete annotation vector:
// per component log10 term pMixLog10 + (-0.5*q)/log(10) + cachedDenomLog10
// (with the two-phase crossProdTmp quadratic and precomputeDenominator-
// ForEvaluation composition), combined by log10sumLog10.  This is the same
// quantity the score kernel computes, kept here for the host-only
// select-worst and missing-annotation paths so every LOD uses one numeric
// definition (matching VariantRecalibratorEngine's evaluateDatum).
double mixture_log10_probability(const GaussianMixture& mixture,
                                 const std::vector<double>& values) {
    const auto dimensions = mixture.means.empty() ? 0 : mixture.means.front().size();
    if (dimensions == 0 || values.size() != dimensions || mixture.weights.empty())
        return -std::numeric_limits<double>::infinity();
    if (mixture.eval_log10_denominators.size() != mixture.weights.size())
        return -std::numeric_limits<double>::infinity();
    const auto components = mixture.weights.size();
    const double log_ten = std::log(10.0);
    std::vector<double> cross(dimensions, 0.0);
    std::vector<double> terms(components, -std::numeric_limits<double>::infinity());
    for (std::size_t component = 0; component < components; ++component) {
        if (component >= mixture.precision.size() ||
            component >= mixture.means.size() ||
            !(mixture.weights[component] > 0.0))
            continue;
        // Java evaluateDatumLog10's crossProdTmp: tmp[iii] +=
        // (x[jjj]-mu[jjj]) * Inv[jjj][iii], then the inner product.
        std::fill(cross.begin(), cross.end(), 0.0);
        for (std::size_t first = 0; first < dimensions; ++first)
            for (std::size_t second = 0; second < dimensions; ++second)
                cross[first] += (values[second] - mixture.means[component][second]) *
                    mixture.precision[component][second * dimensions + first];
        double quadratic = 0.0;
        for (std::size_t first = 0; first < dimensions; ++first)
            quadratic += cross[first] * (values[first] - mixture.means[component][first]);
        terms[component] = std::log10(mixture.weights[component]) +
            ((-0.5 * quadratic) / log_ten) +
            mixture.eval_log10_denominators[component];
    }
    if (components == 1) return terms[0];
    std::size_t maximum_index = 0;
    for (std::size_t component = 1; component < components; ++component)
        if (terms[component] > terms[maximum_index]) maximum_index = component;
    const double maximum = terms[maximum_index];
    if (!std::isfinite(maximum)) return maximum;
    // log10sumLog10: sum = 1.0 + sum_{i != max} 10^(x_i - max).
    double sum = 1.0;
    for (std::size_t component = 0; component < components; ++component) {
        if (component == maximum_index || !std::isfinite(terms[component])) continue;
        sum += std::pow(10.0, terms[component] - maximum);
    }
    return maximum + std::log10(sum);
}

// Java MathUtils.normalDistributionLog10(mean, sd, x).  GATK's VQSR call site
// (GaussianMixtureModel.evaluateDatumInOneDimension) passes the covariance
// diagonal gaussian.sigma.get(iii, iii) in the "sd" slot, i.e. a variance; the
// formula below -- a = -log10(sd * sqrt(2*pi)); b = -(x-mean)^2 / (2*sd^2) /
// log(10) -- reproduces that call verbatim, including the naming quirk.
double gatk_normal_distribution_log10(double mean, double sd, double x) {
    const double root_two_pi = std::sqrt(2.0 * 3.14159265358979323846);
    const double a = -1.0 * std::log10(sd * root_two_pi);
    const double b = -1.0 * (((x - mean) * (x - mean)) / (2.0 * sd * sd)) / std::log(10.0);
    return a + b;
}

// Java GaussianMixtureModel.evaluateDatumInOneDimension: -- the single-
// annotation log10 mixture likelihood that
// VariantRecalibratorEngine.calculateWorstPerformingAnnotation compares
// between the positive and the negative model in order to pick the `culprit`
// annotation.  Every mixture component contributes
// pMixtureLog10 + normalDistributionLog10(mu[iii], sigma[iii][iii], x[iii]),
// combined by log10sumLog10 (NaN-tolerant, as in the Java helper).
double mixture_log10_probability_one_dimension(const GaussianMixture& mixture,
                                               double value, std::size_t dimension) {
    const auto dimensions = mixture.means.empty() ? 0 : mixture.means.front().size();
    if (dimensions == 0 || dimension >= dimensions || mixture.weights.empty())
        return -std::numeric_limits<double>::infinity();
    const auto components = mixture.weights.size();
    std::vector<double> terms(components, -std::numeric_limits<double>::infinity());
    for (std::size_t component = 0; component < components; ++component) {
        // Java: pVarInGaussianLog10 = gaussian.pMixtureLog10, and the normal
        // term is added only when that log10 weight is not -infinity.
        const double mixture_log10 = std::log10(mixture.weights[component]);
        terms[component] = mixture_log10;
        if (mixture_log10 == -std::numeric_limits<double>::infinity()) continue;
        if (component >= mixture.means.size() ||
            mixture.means[component].size() != dimensions) continue;
        double sigma = 0.0;
        if (mixture.full_covariance && component < mixture.covariances.size() &&
            mixture.covariances[component].size() == dimensions * dimensions)
            sigma = mixture.covariances[component][dimension * dimensions + dimension];
        else if (component < mixture.variances.size() &&
                 dimension < mixture.variances[component].size())
            sigma = mixture.variances[component][dimension];
        terms[component] += gatk_normal_distribution_log10(
            mixture.means[component][dimension], sigma, value);
    }
    if (components == 1) return terms[0];
    // Java nanTolerantLog10SumLog10: a NaN term makes the whole result NaN,
    // and `NaN < minProb` is false, so that dimension is skipped by the caller.
    for (const auto term : terms)
        if (std::isnan(term)) return std::numeric_limits<double>::quiet_NaN();
    std::size_t maximum_index = 0;
    for (std::size_t component = 1; component < components; ++component)
        if (terms[component] > terms[maximum_index]) maximum_index = component;
    const double maximum = terms[maximum_index];
    if (!std::isfinite(maximum)) return maximum;
    double sum = 1.0;
    for (std::size_t component = 0; component < components; ++component) {
        if (component == maximum_index || !std::isfinite(terms[component])) continue;
        sum += std::pow(10.0, terms[component] - maximum);
    }
    return maximum + std::log10(sum);
}

// Java VariantRecalibratorEngine.calculateWorstPerformingAnnotation: keep the
// dimension with the smallest good-minus-bad one-dimension log10 likelihood.
// The comparison is Java's strict `prob < minProb` (minProb starts at
// Double.MAX_VALUE), so ties keep the first dimension in GATK's information
// order and a null (missing) dimension is skipped entirely.
void calculate_worst_annotation(std::vector<Entry>& entries, const GaussianMixture& good,
                                const GaussianMixture& bad) {
    for (auto& entry : entries) {
        entry.worst_annotation = -1;
        double minimum_probability = std::numeric_limits<double>::max();
        for (std::size_t dimension = 0; dimension < entry.values.size(); ++dimension) {
            if (dimension < entry.missing.size() && entry.missing[dimension] != 0) continue;
            const double good_probability =
                mixture_log10_probability_one_dimension(good, entry.values[dimension], dimension);
            const double bad_probability =
                mixture_log10_probability_one_dimension(bad, entry.values[dimension], dimension);
            const double probability = good_probability - bad_probability;
            if (probability < minimum_probability) {
                minimum_probability = probability;
                entry.worst_annotation = static_cast<int>(dimension);
            }
        }
    }
}

struct BadVariantSelection {
    std::size_t cutoff_selected = 0;
    std::size_t fallback_selected = 0;
};

BadVariantSelection select_worst_variants(std::vector<Entry>& entries,
                                          const GaussianMixture& good,
                                          double bad_lod_score_cutoff,
                                          std::size_t minimum_bad_variants) {
    struct ScoredIndex { std::size_t index = 0; double lod = 0.0; };
    std::vector<ScoredIndex> candidates;
    candidates.reserve(entries.size());
    for (std::size_t index = 0; index < entries.size(); ++index) {
        auto& entry = entries[index];
        entry.bad_candidate = false;
        // GATK's selectWorstVariants walks the complete data list.  Positive
        // training rows are eligible when their good-model LOD is below the
        // cutoff (and are then marked NEGATIVE_TRAIN_SITE as well); do not
        // pre-filter them here.
        if (!entry.scoreable || entry.failing_std_threshold || entry.negative_training ||
            std::any_of(entry.missing.begin(), entry.missing.end(),
                        [](std::uint8_t value) { return value != 0; }))
            continue;
        const auto log_probability = mixture_log10_probability(good, entry.values);
        if (!std::isfinite(log_probability)) continue;
        const auto lod = log_probability;
        candidates.push_back(ScoredIndex{index, lod});
        if (lod <= bad_lod_score_cutoff) entry.bad_candidate = true;
    }
    BadVariantSelection selection;
    for (const auto& candidate : candidates)
        if (entries[candidate.index].bad_candidate) ++selection.cutoff_selected;
    // A tiny fixture often has no datum below the production cutoff.  GATK's
    // command then warns and may fail its negative model; selecting the
    // deterministic worst rows keeps the native workflow usable while still
    // honoring the cutoff whenever it yields enough data.
    if (selection.cutoff_selected == 0 && !candidates.empty()) {
        std::stable_sort(candidates.begin(), candidates.end(), [](const ScoredIndex& left,
                                                                  const ScoredIndex& right) {
            if (left.lod != right.lod) return left.lod < right.lod;
            return left.index < right.index;
        });
        const auto fallback_count = std::min(minimum_bad_variants, candidates.size());
        for (std::size_t rank = 0; rank < fallback_count; ++rank) {
            entries[candidates[rank].index].bad_candidate = true;
            ++selection.fallback_selected;
        }
    }
    return selection;
}

double marginalized_log10_probability(const GaussianMixture& mixture,
                                      const std::vector<double>& values,
                                      const std::vector<std::uint8_t>& missing,
                                      JavaRandom& random) {
    if (values.empty() || missing.size() != values.size())
        return -std::numeric_limits<double>::infinity();
    std::vector<double> sampled = values;
    double probability_sum = 0.0;
    std::size_t draws = 0;
    // This mirrors GaussianMixtureModel.evaluateDatumMarginalized: each
    // missing dimension receives twenty independent N(0,1) draws and the
    // resulting mixture probabilities are averaged.  The values of other
    // missing dimensions remain at their normalized jitter values, matching
    // the Java loop's in-place update order.
    for (std::size_t dimension = 0; dimension < values.size(); ++dimension) {
        if (missing[dimension] == 0) continue;
        for (int draw = 0; draw < 20; ++draw) {
            sampled[dimension] = random.next_gaussian();
            const auto log_probability = mixture_log10_probability(mixture, sampled);
            if (std::isfinite(log_probability)) probability_sum += std::pow(10.0, log_probability);
            ++draws;
        }
    }
    if (draws == 0 || !(probability_sum > 0.0))
        return -std::numeric_limits<double>::infinity();
    return std::log10(probability_sum / static_cast<double>(draws));
}

void score_entries(std::vector<Entry>& entries, const GaussianMixture& good,
                   const GaussianMixture& bad, int threads,
                   KernelTelemetry* telemetry = nullptr,
                   JavaRandom* random = nullptr) {
    std::vector<std::size_t> indices;
    for (std::size_t index = 0; index < entries.size(); ++index)
        if (entries[index].scoreable) indices.push_back(index);
    if (indices.empty()) return;
    const std::size_t dimensions = good.means.empty() ? 0 : good.means.front().size();
    const std::size_t good_components = good.weights.size();
    const std::size_t bad_components = bad.weights.size();
    std::vector<double> input(indices.size() * dimensions);
    std::vector<std::uint8_t> input_missing(indices.size() * dimensions, 0);
    for (std::size_t row = 0; row < indices.size(); ++row)
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
            input[row * dimensions + dimension] = entries[indices[row]].values[dimension];
            if (dimension < entries[indices[row]].missing.size())
                input_missing[row * dimensions + dimension] = entries[indices[row]].missing[dimension];
        }
    std::vector<double> final_scores(indices.size(), std::numeric_limits<double>::quiet_NaN());
    {
        Kokkos::View<double**> device_input("variant_recalibrator_input", indices.size(), dimensions);
        Kokkos::View<std::uint8_t**> device_missing("variant_recalibrator_missing", indices.size(), dimensions);
        Kokkos::View<double**> device_good_mean("variant_recalibrator_good_mean", good_components, dimensions);
        Kokkos::View<double**> device_bad_mean("variant_recalibrator_bad_mean", bad_components, dimensions);
        Kokkos::View<double***> device_good_precision("variant_recalibrator_good_precision",
                                                       good_components, dimensions, dimensions);
        Kokkos::View<double***> device_bad_precision("variant_recalibrator_bad_precision",
                                                      bad_components, dimensions, dimensions);
        Kokkos::View<double*> device_good_logdet("variant_recalibrator_good_logdet", good_components);
        Kokkos::View<double*> device_bad_logdet("variant_recalibrator_bad_logdet", bad_components);
        Kokkos::View<double*> device_good_weight("variant_recalibrator_good_weight", good_components);
        Kokkos::View<double*> device_bad_weight("variant_recalibrator_bad_weight", bad_components);
        Kokkos::View<double*> device_scores("variant_recalibrator_scores", indices.size());
        Kokkos::View<double*> device_cross("variant_recalibrator_score_cross",
                                           indices.size() * dimensions);
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host_batch("variant-recalibrator-score-v1");
        host_batch.records = indices.size();
        host_batch.bytes = input.size() * sizeof(double) + input_missing.size() * sizeof(std::uint8_t);
        fastgatk::core::KernelPlan<ExecSpace> kernel_plan("variant-recalibrator-score");
        kernel_plan.begin_prepare(host_batch);
        auto host_input = Kokkos::create_mirror_view(device_input);
        auto host_missing = Kokkos::create_mirror_view(device_missing);
        auto host_good_mean = Kokkos::create_mirror_view(device_good_mean);
        auto host_bad_mean = Kokkos::create_mirror_view(device_bad_mean);
        auto host_good_precision = Kokkos::create_mirror_view(device_good_precision);
        auto host_bad_precision = Kokkos::create_mirror_view(device_bad_precision);
        auto host_good_logdet = Kokkos::create_mirror_view(device_good_logdet);
        auto host_bad_logdet = Kokkos::create_mirror_view(device_bad_logdet);
        auto host_good_weight = Kokkos::create_mirror_view(device_good_weight);
        auto host_bad_weight = Kokkos::create_mirror_view(device_bad_weight);
        for (std::size_t row = 0; row < indices.size(); ++row)
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                host_input(row, dimension) = input[row * dimensions + dimension];
                host_missing(row, dimension) = input_missing[row * dimensions + dimension];
            }
        for (std::size_t component = 0; component < good_components; ++component) {
            host_good_weight(component) = good.weights[component];
            host_good_logdet(component) = good.eval_log10_denominators[component];
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                host_good_mean(component, dimension) = good.means[component][dimension];
                for (std::size_t other = 0; other < dimensions; ++other)
                    host_good_precision(component, dimension, other) =
                        good.precision[component][dimension * dimensions + other];
            }
        }
        for (std::size_t component = 0; component < bad_components; ++component) {
            host_bad_weight(component) = bad.weights[component];
            host_bad_logdet(component) = bad.eval_log10_denominators[component];
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                host_bad_mean(component, dimension) = bad.means[component][dimension];
                for (std::size_t other = 0; other < dimensions; ++other)
                    host_bad_precision(component, dimension, other) =
                        bad.precision[component][dimension * dimensions + other];
            }
        }
        Kokkos::deep_copy(device_input, host_input);
        Kokkos::deep_copy(device_missing, host_missing);
        Kokkos::deep_copy(device_good_mean, host_good_mean);
        Kokkos::deep_copy(device_bad_mean, host_bad_mean);
        Kokkos::deep_copy(device_good_precision, host_good_precision);
        Kokkos::deep_copy(device_bad_precision, host_bad_precision);
        Kokkos::deep_copy(device_good_logdet, host_good_logdet);
        Kokkos::deep_copy(device_bad_logdet, host_bad_logdet);
        Kokkos::deep_copy(device_good_weight, host_good_weight);
        Kokkos::deep_copy(device_bad_weight, host_bad_weight);
        fastgatk::core::DeviceBatch<ExecSpace> device_batch(indices.size());
        device_batch.bind("input", device_input);
        device_batch.bind("missing", device_missing);
        device_batch.bind("good_mean", device_good_mean);
        device_batch.bind("bad_mean", device_bad_mean);
        device_batch.bind("good_precision", device_good_precision);
        device_batch.bind("bad_precision", device_bad_precision);
        device_batch.bind("good_logdet", device_good_logdet);
        device_batch.bind("bad_logdet", device_bad_logdet);
        device_batch.bind("good_weight", device_good_weight);
        device_batch.bind("bad_weight", device_bad_weight);
        device_batch.bind("scores", device_scores);
        device_batch.bind("score_cross", device_cross);
        ExecSpace().fence();
        kernel_plan.end_prepare(device_batch);
        kernel_plan.begin_execute();
        const double score_log_ten = std::log(10.0);
        Kokkos::parallel_for("variant_recalibrator_gaussian_score",
            Kokkos::RangePolicy<ExecSpace>(0, indices.size()), KOKKOS_LAMBDA(const std::size_t row) {
                // Java's non-variational evaluator builds each log10 Gaussian
                // term as pMixtureLog10 + (-0.5*q)/log(10) + cachedDenomLog10
                // (precomputeDenominatorForEvaluation) and combines terms with
                // log10sumLog10: max + log10(1 + sum_{i!=max} 10^(x_i - max)).
                auto component_log10 = [&](const std::size_t component,
                                           const Kokkos::View<double*>& weight,
                                           const Kokkos::View<double*>& logdet,
                                           const Kokkos::View<double**>& mean,
                                           const Kokkos::View<double***>& precision,
                                           double& quadratic) {
                    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
                        device_cross(row * dimensions + dimension) = 0.0;
                    for (std::size_t first = 0; first < dimensions; ++first)
                        for (std::size_t second = 0; second < dimensions; ++second)
                            if (device_missing(row, first) == 0 &&
                                device_missing(row, second) == 0)
                                device_cross(row * dimensions + first) +=
                                    (device_input(row, second) - mean(component, second)) *
                                    precision(component, second, first);
                    quadratic = 0.0;
                    for (std::size_t first = 0; first < dimensions; ++first)
                        if (device_missing(row, first) == 0)
                            quadratic += device_cross(row * dimensions + first) *
                                (device_input(row, first) - mean(component, first));
                    return log10(weight(component)) +
                        ((-0.5 * quadratic) / score_log_ten) + logdet(component);
                };
                double good_values = -1.0e300;
                double bad_values = -1.0e300;
                std::size_t good_max_index = 0;
                std::size_t bad_max_index = 0;
                for (std::size_t component = 0; component < good_components; ++component) {
                    double quadratic = 0.0;
                    const double value = component_log10(
                        component, device_good_weight, device_good_logdet,
                        device_good_mean, device_good_precision, quadratic);
                    if (value > good_values) {
                        good_values = value;
                        good_max_index = component;
                    }
                }
                for (std::size_t component = 0; component < bad_components; ++component) {
                    double quadratic = 0.0;
                    const double value = component_log10(
                        component, device_bad_weight, device_bad_logdet,
                        device_bad_mean, device_bad_precision, quadratic);
                    if (value > bad_values) {
                        bad_values = value;
                        bad_max_index = component;
                    }
                }
                double good_sum = 1.0;
                double bad_sum = 1.0;
                for (std::size_t component = 0; component < good_components; ++component) {
                    if (component == good_max_index) continue;
                    double quadratic = 0.0;
                    const double value = component_log10(
                        component, device_good_weight, device_good_logdet,
                        device_good_mean, device_good_precision, quadratic);
                    good_sum += pow(10.0, value - good_values);
                }
                for (std::size_t component = 0; component < bad_components; ++component) {
                    if (component == bad_max_index) continue;
                    double quadratic = 0.0;
                    const double value = component_log10(
                        component, device_bad_weight, device_bad_logdet,
                        device_bad_mean, device_bad_precision, quadratic);
                    bad_sum += pow(10.0, value - bad_values);
                }
                device_scores(row) =
                    (good_values + log10(good_sum)) -
                    (bad_values + log10(bad_sum));
            });
        ExecSpace().fence();
        kernel_plan.end_execute();
        if (telemetry != nullptr) {
            ++telemetry->batches;
            telemetry->observations += indices.size();
            telemetry->prepare_seconds += kernel_plan.telemetry().prepare_seconds;
            telemetry->execute_seconds += kernel_plan.telemetry().execute_seconds;
        }
        auto host_scores = Kokkos::create_mirror_view(device_scores);
        Kokkos::deep_copy(host_scores, device_scores);
        for (std::size_t row = 0; row < indices.size(); ++row)
            final_scores[row] = host_scores(row);
    }
    // The regular scorer handles complete vectors in one Kokkos batch.  GATK
    // intentionally uses a seeded Monte-Carlo marginalization for a missing
    // annotation, however, and the RNG is global/ordered.  Re-evaluate only
    // those rows on Host in the same two-pass order as
    // VariantRecalibratorEngine (all good-model rows, then all bad-model
    // rows), leaving complete rows on the device path.
    if (random != nullptr) {
        std::vector<double> marginalized_good(indices.size(),
            std::numeric_limits<double>::quiet_NaN());
        std::vector<double> marginalized_bad(indices.size(),
            std::numeric_limits<double>::quiet_NaN());
        for (std::size_t row = 0; row < indices.size(); ++row) {
            const auto& entry = entries[indices[row]];
            if (std::none_of(entry.missing.begin(), entry.missing.end(),
                             [](std::uint8_t value) { return value != 0; })) continue;
            marginalized_good[row] = marginalized_log10_probability(
                good, entry.values, entry.missing, *random);
        }
        for (std::size_t row = 0; row < indices.size(); ++row) {
            const auto& entry = entries[indices[row]];
            if (std::none_of(entry.missing.begin(), entry.missing.end(),
                             [](std::uint8_t value) { return value != 0; })) continue;
            marginalized_bad[row] = marginalized_log10_probability(
                bad, entry.values, entry.missing, *random);
        }
        for (std::size_t row = 0; row < indices.size(); ++row) {
            if (std::isfinite(marginalized_good[row]) && std::isfinite(marginalized_bad[row]))
                final_scores[row] = marginalized_good[row] - marginalized_bad[row];
        }
    }
    (void)threads;
    for (std::size_t row = 0; row < indices.size(); ++row) entries[indices[row]].score = final_scores[row];
}

struct TrancheRow {
    double requested = 0.0;
    double min_lod = -std::numeric_limits<double>::infinity();
    std::int64_t known = 0;
    std::int64_t novel = 0;
    std::int64_t known_ti = 0;
    std::int64_t known_tv = 0;
    std::int64_t novel_ti = 0;
    std::int64_t novel_tv = 0;
    std::int64_t accessible_truth_sites = 0;
    std::int64_t calls_at_truth_sites = 0;
};

double titv_ratio(std::int64_t transitions, std::int64_t transversions) {
    return static_cast<double>(transitions) /
           static_cast<double>(std::max<std::int64_t>(transversions, 1));
}

void write_tranches(const std::string& path, const std::vector<Entry>& entries,
                    const std::vector<double>& targets, const std::string& mode,
                    bool scatter, const std::vector<double>& vqslod_targets) {
    std::vector<const Entry*> scored;
    std::int64_t accessible_truth_sites = 0;
    for (const auto& entry : entries) {
        if (!entry.scoreable || !std::isfinite(entry.score)) continue;
        scored.push_back(&entry);
        if (entry.truth) ++accessible_truth_sites;
    }
    std::stable_sort(scored.begin(), scored.end(), [](const Entry* lhs, const Entry* rhs) {
        return lhs->score > rhs->score;
    });
    std::vector<double> truth_scores;
    truth_scores.reserve(static_cast<std::size_t>(accessible_truth_sites));
    for (const auto* entry : scored)
        if (entry->truth) truth_scores.push_back(entry->score);

    auto build_row = [&](double requested, bool use_vqslod) {
        TrancheRow row;
        row.requested = requested;
        bool has_variants = false;
        if (use_vqslod) {
            const auto iterator = std::find_if(scored.begin(), scored.end(),
                [&](const Entry* entry) { return entry->score >= requested; });
            has_variants = iterator != scored.end();
            if (has_variants) row.min_lod = (*iterator)->score;
            else if (!scored.empty()) row.min_lod = scored.back()->score;
        } else {
            row.min_lod = truth_scores.empty() ? 0.0 : truth_scores.back();
            if (!truth_scores.empty()) {
                const auto required = static_cast<std::size_t>(std::ceil(requested / 100.0 * truth_scores.size()));
                const auto index = required >= truth_scores.size()
                    ? truth_scores.size() - 1 : truth_scores.size() - required;
                row.min_lod = truth_scores[index];
                has_variants = true;
            }
        }
        row.accessible_truth_sites = accessible_truth_sites;
        for (const auto* entry : scored) {
            const bool include = use_vqslod
                ? (has_variants && entry->score >= row.min_lod)
                : (entry->score >= row.min_lod);
            if (!include) continue;
            if (entry->known) {
                ++row.known;
                if (entry->is_snp) (entry->is_transition ? ++row.known_ti : ++row.known_tv);
            } else {
                ++row.novel;
                if (entry->is_snp) (entry->is_transition ? ++row.novel_ti : ++row.novel_tv);
            }
            if (entry->truth) ++row.calls_at_truth_sites;
        }
        // Java's VQSLODTranche.emptyTranche reports zero known/novel counts
        // but still records calls at the lowest observed LOD.  Preserve that
        // edge case when a requested scatter slice is above all observations.
        if (use_vqslod && !has_variants) {
            row.known = row.novel = row.known_ti = row.known_tv = row.novel_ti = row.novel_tv = 0;
            row.calls_at_truth_sites = 0;
            for (const auto* entry : scored)
                if (entry->truth && entry->score >= row.min_lod) ++row.calls_at_truth_sites;
        }
        return row;
    };

    std::vector<TrancheRow> rows;
    const auto& requested = scatter ? vqslod_targets : targets;
    rows.reserve(requested.size());
    for (const auto value : requested) rows.push_back(build_row(value, scatter));
    if (scatter && rows.size() > 1) {
        std::stable_sort(rows.begin(), rows.end(), [](const TrancheRow& lhs, const TrancheRow& rhs) {
            return lhs.calls_at_truth_sites < rhs.calls_at_truth_sites;
        });
    }

    std::ofstream output(path);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write tranches: " + path);
    if (scatter) {
        output << "# Variant quality score tranches file\n# Version number 6\n"
                  "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";
    } else {
        output << "# Variant quality score tranches file\n# Version number 5\n"
                  "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";
    }
    output << std::fixed << std::setprecision(2);
    double previous = 0.0;
    for (const auto& row : rows) {
        const double sensitivity = row.accessible_truth_sites == 0 ? 0.0
            : static_cast<double>(row.calls_at_truth_sites) / row.accessible_truth_sites;
        const double known_titv = titv_ratio(row.known_ti, row.known_tv);
        const double novel_titv = titv_ratio(row.novel_ti, row.novel_tv);
        output << row.requested << ',' << row.known << ',' << row.novel << ','
               << std::setprecision(4) << known_titv << ',' << novel_titv << ','
               << row.min_lod << ",VQSRTranche" << mode
               << std::setprecision(2) << previous << "to" << row.requested << ',' << mode << ','
               << row.accessible_truth_sites << ',' << row.calls_at_truth_sites << ','
               << std::setprecision(4) << sensitivity << '\n';
        previous = row.requested;
        output << std::setprecision(2);
    }
    output.close();
    if (!output || !complete_file(path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete tranches");
}

#if FASTGATK_HAS_HTSLIB

std::vector<Entry> read_entries(const Options& options,
                                const std::unordered_map<std::string, std::unordered_set<std::string>>& resource_sets,
                                std::uint64_t& input_records,
                                JavaRandom* random = nullptr) {
    htsFile* input = bcf_open(options.input.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open input VCF: " + options.input);
    bcf_hdr_t* header = bcf_hdr_read(input);
    bcf1_t* record = bcf_init();
    if (!header || !record) throw std::runtime_error("BAD_INPUT: cannot read input VCF header");
    std::vector<Entry> entries;
    while (bcf_read(input, header, record) == 0) {
        ++input_records;
        const auto type = variant_type(record);
        const bool mode_match = mode_matches_variant(type, options.mode);
        const auto keys = record_keys(header, record);
        bool record_training = false, record_truth = false, record_known = false;
        int positive_flag_count = 0;
        int negative_flag_count = 0;
        const bool input_positive_label = bcf_get_info_flag(
            header, record, "POSITIVE_TRAIN_SITE", nullptr, &positive_flag_count) > 0;
        const bool input_negative_label = bcf_get_info_flag(
            header, record, "NEGATIVE_TRAIN_SITE", nullptr, &negative_flag_count) > 0;
        double record_prior = 2.0;
        for (const auto& resource : options.resources) {
            const auto& keys_for_resource = resource_sets.at(resource.path);
            bool overlaps = false;
            for (const auto& key : keys) if (keys_for_resource.count(key)) { overlaps = true; break; }
            record_training = record_training || (resource.training && overlaps);
            record_truth = record_truth || (resource.truth && overlaps);
            record_known = record_known || (resource.known && overlaps);
            if (overlaps) record_prior = std::max(record_prior, resource.prior);
        }
        const std::size_t alt_count = options.allele_specific ? keys.size() : std::min<std::size_t>(1, keys.size());
        for (std::size_t alt = 0; alt < alt_count; ++alt) {
            std::vector<double> values;
            std::vector<std::uint8_t> missing;
            // VariantRecalibrator's hidden --sample-every-Nth-variant is
            // zero-based at the traversal boundary: the first input record
            // is retained, then every Nth record thereafter.  input_records
            // is incremented before decoding, so subtract one here to keep
            // Java's counter phase (and the sampled recalibration-table
            // membership) exact.
            const bool sampled = input_records > 0 &&
                ((input_records - 1) % options.sample_every == 0);
            const bool scoreable = sampled && mode_match && annotation_values(
                header, record, options.annotations, values, missing,
                static_cast<int>(alt), options.allele_specific, random,
                options.mq_cap, options.mq_jitter);
            Entry entry;
            entry.key = keys[alt];
            entry.type = type;
            entry.values = std::move(values);
            entry.missing = std::move(missing);
            entry.scoreable = scoreable;
            if (record->n_allele > static_cast<int>(alt + 1) && record->d.allele != nullptr) {
                const std::string reference_allele(record->d.allele[0]);
                const std::string alternate_allele(record->d.allele[alt + 1]);
                entry.is_snp = reference_allele.size() == 1 && alternate_allele.size() == 1;
                entry.is_transition = entry.is_snp &&
                    is_transition_base(reference_allele.front(), alternate_allele.front());
            }
            entry.prior_lod = prior_quality_to_log10_odds(record_prior);
            if (options.allele_specific) {
                double allele_prior = 2.0;
                for (const auto& resource : options.resources) {
                    const auto& keys_for_resource = resource_sets.at(resource.path);
                    const bool overlaps = keys_for_resource.count(entry.key) != 0;
                    entry.training = entry.training || (resource.training && overlaps);
                    entry.truth = entry.truth || (resource.truth && overlaps);
                    entry.known = entry.known || (resource.known && overlaps);
                    if (overlaps) allele_prior = std::max(allele_prior, resource.prior);
                }
                entry.prior_lod = prior_quality_to_log10_odds(allele_prior);
                entry.training = entry.training || input_positive_label;
            } else {
                entry.training = record_training || input_positive_label;
                entry.negative_training = input_negative_label;
                entry.truth = record_truth;
                entry.known = record_known;
            }
            if (options.allele_specific) entry.negative_training = input_negative_label;
            if (options.resources.empty()) entry.training = entry.truth = true;
            entries.push_back(std::move(entry));
        }
        bcf_clear(record);
    }
    bcf_destroy(record); bcf_hdr_destroy(header); bcf_close(input);
    return entries;
}

void normalize_entries(std::vector<Entry>& entries, const ModelBundle& model,
                       JavaRandom* random = nullptr) {
    if (model.annotation_means.size() != model.annotation_stdevs.size())
        throw std::runtime_error("BAD_INPUT: VQSR model normalization vectors have different sizes");
    for (auto& entry : entries) {
        if (!entry.scoreable) continue;
        if (entry.values.size() != model.annotation_means.size()) {
            entry.scoreable = false;
            entry.values.clear();
            continue;
        }
        for (std::size_t dimension = 0; dimension < entry.values.size(); ++dimension) {
            const bool missing = dimension < entry.missing.size() && entry.missing[dimension] != 0;
            entry.values[dimension] = missing
                ? (random == nullptr ? 0.0 : 0.1 * random->next_gaussian())
                : (entry.values[dimension] - model.annotation_means[dimension]) /
                  std::max(1.0e-12, model.annotation_stdevs[dimension]);
        }
    }
}

// GATK's Gaussian model is fit in a normalized annotation space and the
// normalization vectors are serialized alongside the mixture.  Keep this
// computation deterministic and independent of the execution backend so a
// model produced on Host can be replayed by the Kokkos scorer on another
// backend.
void compute_annotation_normalization(const std::vector<Entry>& entries,
                                      std::size_t dimensions,
                                      std::vector<double>& means,
                                      std::vector<double>& stdevs) {
    means.assign(dimensions, 0.0);
    stdevs.assign(dimensions, 1.0);
    std::vector<std::size_t> counts(dimensions, 0);
    for (const auto& entry : entries) {
        // VariantDataManager computes means/stdevs from positive training
        // sites, not from the full evaluation callset.  Keeping this subset
        // stable is important for model reuse and for GATK-compatible scores.
        if (!entry.scoreable || !entry.training || entry.values.size() != dimensions) continue;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
            const bool missing = dimension < entry.missing.size() && entry.missing[dimension] != 0;
            if (missing) continue;
            means[dimension] += entry.values[dimension];
            ++counts[dimension];
        }
    }
    if (std::all_of(counts.begin(), counts.end(), [](std::size_t count) { return count == 0; }))
        throw std::runtime_error("BAD_INPUT: cannot normalize an empty VQSR training set");
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension)
        if (counts[dimension] != 0) means[dimension] /= static_cast<double>(counts[dimension]);
    std::vector<double> variances(dimensions, 0.0);
    for (const auto& entry : entries) {
        if (!entry.scoreable || !entry.training || entry.values.size() != dimensions) continue;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
            const bool missing = dimension < entry.missing.size() && entry.missing[dimension] != 0;
            if (missing || counts[dimension] == 0) continue;
            const auto delta = entry.values[dimension] - means[dimension];
            variances[dimension] += delta * delta;
        }
    }
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
        // A constant annotation is valid; use unit scale rather than
        // manufacturing an infinite normalized value.  GATK does not accept
        // this as a usable model, however: VariantDataManager aborts the
        // normalization pass when any training annotation has standard
        // deviation < 1e-5.  Silently substituting unit scale here used to
        // let a native run publish a model whose likelihood geometry differed
        // from Java (and could make the subsequent covariance singular).
        const auto variance = counts[dimension] == 0
            ? 0.0 : variances[dimension] / static_cast<double>(counts[dimension]);
        if (counts[dimension] != 0 &&
            (!std::isfinite(variance) || std::sqrt(std::max(0.0, variance)) < 1.0e-5))
            throw std::runtime_error("BAD_INPUT: Found annotations with zero variance. They must be excluded before proceeding.");
        stdevs[dimension] = variance > 1.0e-24 && std::isfinite(variance)
            ? std::sqrt(variance) : 1.0;
    }
}

// VariantDataManager.normalizeData reorders dimensions after normalization so
// that model initialization is independent of the command-line annotation
// order.  The ordering key is the absolute shift between the training mean
// and the non-training mean (descending); Java's stable Collections.sort
// preserves the original order for ties.  This is not cosmetic: the order
// controls the JavaRandom stream used by random covariance initialization and
// therefore changes multi-Gaussian VBEM models.
std::vector<std::size_t> gatk_annotation_order(const std::vector<Entry>& entries,
                                                const std::vector<double>& means,
                                                std::size_t dimensions) {
    std::vector<double> non_training_sum(dimensions, 0.0);
    std::vector<std::size_t> non_training_count(dimensions, 0);
    for (const auto& entry : entries) {
        if (!entry.scoreable || entry.training || entry.values.size() != dimensions) continue;
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
            const bool missing = dimension < entry.missing.size() && entry.missing[dimension] != 0;
            if (missing) continue;
            non_training_sum[dimension] += entry.values[dimension];
            ++non_training_count[dimension];
        }
    }
    std::vector<std::size_t> order(dimensions);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        const auto left_mean = non_training_count[left] == 0
            ? std::numeric_limits<double>::quiet_NaN()
            : non_training_sum[left] / static_cast<double>(non_training_count[left]);
        const auto right_mean = non_training_count[right] == 0
            ? std::numeric_limits<double>::quiet_NaN()
            : non_training_sum[right] / static_cast<double>(non_training_count[right]);
        const auto left_key = std::isfinite(left_mean) && left < means.size()
            ? -std::abs(means[left] - left_mean) : std::numeric_limits<double>::quiet_NaN();
        const auto right_key = std::isfinite(right_mean) && right < means.size()
            ? -std::abs(means[right] - right_mean) : std::numeric_limits<double>::quiet_NaN();
        // java.lang.Double.compare orders NaN above all finite values.  Since
        // the key is ascending, finite shifts precede missing non-training
        // means, while equal values remain stable.
        if (std::isnan(left_key) || std::isnan(right_key))
            return std::isnan(right_key) && !std::isnan(left_key);
        return left_key < right_key;
    });
    return order;
}

void reorder_annotation_dimensions(std::vector<Entry>& entries,
                                   std::vector<std::string>& annotations,
                                   std::vector<double>& means,
                                   std::vector<double>& stdevs,
                                   const std::vector<std::size_t>& order) {
    if (order.size() != annotations.size() || means.size() != order.size() ||
        stdevs.size() != order.size())
        throw std::runtime_error("BAD_INPUT: invalid VQSR annotation reorder");
    auto old_annotations = annotations;
    auto old_means = means;
    auto old_stdevs = stdevs;
    for (std::size_t dimension = 0; dimension < order.size(); ++dimension) {
        annotations[dimension] = old_annotations[order[dimension]];
        means[dimension] = old_means[order[dimension]];
        stdevs[dimension] = old_stdevs[order[dimension]];
    }
    for (auto& entry : entries) {
        if (entry.values.size() != order.size() || entry.missing.size() != order.size()) continue;
        auto old_values = entry.values;
        auto old_missing = entry.missing;
        for (std::size_t dimension = 0; dimension < order.size(); ++dimension) {
            entry.values[dimension] = old_values[order[dimension]];
            entry.missing[dimension] = old_missing[order[dimension]];
        }
    }
}

void normalize_entries_in_place(std::vector<Entry>& entries,
                                const std::vector<double>& means,
                                const std::vector<double>& stdevs,
                                JavaRandom* random = nullptr) {
    if (means.size() != stdevs.size())
        throw std::runtime_error("BAD_INPUT: VQSR normalization vectors have different sizes");
    for (auto& entry : entries) {
        if (!entry.scoreable) continue;
        if (entry.values.size() != means.size()) {
            entry.scoreable = false;
            entry.values.clear();
            continue;
        }
        for (std::size_t dimension = 0; dimension < entry.values.size(); ++dimension) {
            const bool missing = dimension < entry.missing.size() && entry.missing[dimension] != 0;
            entry.values[dimension] = missing
                ? (random == nullptr ? 0.0 : 0.1 * random->next_gaussian())
                : (entry.values[dimension] - means[dimension]) /
                  std::max(1.0e-12, stdevs[dimension]);
        }
    }
}

void mark_failing_std_threshold(std::vector<Entry>& entries, double threshold) {
    for (auto& entry : entries) {
        entry.failing_std_threshold = false;
        if (!entry.scoreable || !(threshold > 0.0)) continue;
        entry.failing_std_threshold = std::any_of(
            entry.values.begin(), entry.values.end(), [&](const double value) {
                return std::isfinite(value) && std::abs(value) > threshold;
            });
    }
}

std::string write_recal_vcf(const Options& options, const std::unordered_map<std::string, double>& scores,
                     const std::unordered_map<std::string, bool>& positive_training,
                     const std::unordered_map<std::string, bool>& negative_training,
                     const std::unordered_map<std::string, std::string>& culprits,
                     const std::string& model, std::uint64_t input_records,
                     std::uint64_t& scored_records, bool allele_specific) {
    htsFile* input = bcf_open(options.input.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot reopen input VCF");
    bcf_hdr_t* input_header = bcf_hdr_read(input);
    bcf_hdr_t* output_header = input_header ? bcf_hdr_dup(input_header) : nullptr;
    bcf1_t* record = bcf_init();
    if (!input_header || !output_header || !record) throw std::runtime_error("BAD_INPUT: cannot duplicate VCF header");
    // GATK's VariantDataManager always writes a site-only recalibration
    // artifact: scalar records use the dummy N/<VQSR> allele pair and never
    // carry the input FORMAT/sample payload.  Keep the sample mask active
    // while records are scored, then strip FORMAT only at this writer
    // boundary so input annotations/model decisions are unaffected.
    if (bcf_hdr_set_samples(output_header, nullptr, 0) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot configure sites-only VCF header");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "VQSLOD") < 0)
        bcf_hdr_append(output_header, (std::string("##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=fastgatk ") + model + " VQSR score>").c_str());
    if (allele_specific && bcf_hdr_id2int(output_header, BCF_DT_ID, "AS_VQSLOD") < 0)
        bcf_hdr_append(output_header, (std::string("##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=fastgatk ") + model + " allele-specific VQSR score>").c_str());
    if (allele_specific && bcf_hdr_id2int(output_header, BCF_DT_ID, "AS_culprit") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=AS_culprit,Number=A,Type=String,Description=fastgatk allele-specific VQSR culprit>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "CULPRIT") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=CULPRIT,Number=1,Type=String,Description=fastgatk VQSR model>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "culprit") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=culprit,Number=1,Type=String,Description=GATK VQSR model annotation>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "END") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=END,Number=1,Type=Integer,Description=Stop position of the variant>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "POSITIVE_TRAIN_SITE") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the positive training set of good variants>");
    if (bcf_hdr_id2int(output_header, BCF_DT_ID, "NEGATIVE_TRAIN_SITE") < 0)
        bcf_hdr_append(output_header, "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=This variant was used to build the negative training set of bad variants>");
    bcf_hdr_append(output_header, "##source=fastgatk-variant-recalibrator");
    bcf_hdr_append(output_header, (std::string("##fastgatk_variant_recalibrator_status=prototype-") + model).c_str());
    // GATK writes an annotation NAME here (VariantDataManager.java:485), so the
    // value comes from the per-datum worst-annotation selection computed in
    // main(); the model name is used only for native provenance headers.
    // "NULL" is GATK's spelling for a datum whose dimensions were all null.
    const auto culprit_value = [&culprits](const std::string& key) -> std::string {
        const auto iterator = culprits.find(key);
        return iterator == culprits.end() ? std::string("NULL") : iterator->second;
    };
    if (bcf_hdr_sync(output_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync VCF header");
    htsFile* output = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
    if (!output || bcf_hdr_write(output, output_header) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write recal VCF header");
    std::uint64_t output_record_index = 0;
    while (bcf_read(input, input_header, record) == 0) {
        // GATK's sample-every-Nth-variant mode downsamples the recalibration
        // table itself, not merely the score computation.  Keep the same
        // zero-based phase as read_entries/apply: unsampled input records do
        // not appear in the output VCF at all.
        if (output_record_index++ % options.sample_every != 0) {
            bcf_clear(record);
            continue;
        }
        const auto keys = record_keys(input_header, record);
        // VariantRecalibrator annotates the training membership in the recal
        // VCF.  Remove any stale flags inherited from the input before
        // rebuilding them from the exact allele/resource key that produced
        // the score; this keeps a reused input deterministic.
        bcf_update_info_flag(output_header, record, "POSITIVE_TRAIN_SITE", nullptr, 0);
        bcf_update_info_flag(output_header, record, "NEGATIVE_TRAIN_SITE", nullptr, 0);
        bool positive_label = false;
        bool negative_label = false;
        if (allele_specific) {
            std::vector<float> allele_scores;
            bool any_score = false;
            for (const auto& key : keys) {
                const auto iterator = scores.find(key);
                const float score = iterator == scores.end()
                    ? std::numeric_limits<float>::quiet_NaN() : iterator->second;
                any_score = any_score || iterator != scores.end();
                allele_scores.push_back(score);
            }
            if (any_score && !allele_scores.empty()) {
                if (bcf_update_info_float(output_header, record, "AS_VQSLOD",
                                          allele_scores.data(), static_cast<int>(allele_scores.size())) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update AS_VQSLOD");
                std::ostringstream culprit_list;
                for (std::size_t index = 0; index < allele_scores.size(); ++index) {
                    if (index) culprit_list << ',';
                    culprit_list << culprit_value(keys[index]);
                }
                if (bcf_update_info_string(output_header, record, "AS_culprit", culprit_list.str().c_str()) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update AS_culprit");
                bcf_update_info_string(output_header, record, "CULPRIT", culprit_value(keys.front()).c_str());
                ++scored_records;
            }
        } else {
            const auto iterator = keys.empty() ? scores.end() : scores.find(keys.front());
            if (iterator != scores.end()) {
                // GATK VariantDataManager stores the VQSLOD INFO attribute as
                // String.format("%.4f", datum.lod) TEXT of the full-precision
                // double; write the same text so raw-bit score parity survives
                // serialization.
                char vqslod_buffer[64];
                std::snprintf(vqslod_buffer, sizeof(vqslod_buffer), "%.4f", iterator->second);
                if (bcf_update_info_string(output_header, record, "VQSLOD", vqslod_buffer) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update VQSLOD");
                if (bcf_update_info_string(output_header, record, "culprit",
                                           culprit_value(keys.front()).c_str()) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update culprit");
                ++scored_records;
            }
        }
        for (const auto& key : keys) {
            const auto positive = positive_training.find(key);
            const auto negative = negative_training.find(key);
            positive_label = positive_label || (positive != positive_training.end() && positive->second);
            negative_label = negative_label || (negative != negative_training.end() && negative->second);
        }
        if (positive_label && bcf_update_info_flag(output_header, record,
                                                   "POSITIVE_TRAIN_SITE", nullptr, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write POSITIVE_TRAIN_SITE");
        if (negative_label && bcf_update_info_flag(output_header, record,
                                                   "NEGATIVE_TRAIN_SITE", nullptr, 1) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write NEGATIVE_TRAIN_SITE");
        if (!allele_specific) {
            // VariantDataManager::writeOutRecalibrationTable intentionally
            // replaces the original callset allele with a symbolic marker;
            // ApplyVQSR joins this record back by END/coordinate.  Remove
            // stale callset annotations before writing the recal artifact.
            bcf_unpack(record, BCF_UN_INFO);
            std::vector<std::string> stale_info;
            stale_info.reserve(static_cast<std::size_t>(record->n_info));
            for (int info_index = 0; info_index < record->n_info; ++info_index) {
                const char* id = bcf_hdr_int2id(output_header, BCF_DT_ID,
                                                record->d.info[info_index].key);
                if (id != nullptr && std::strcmp(id, "END") != 0 &&
                    std::strcmp(id, "VQSLOD") != 0 && std::strcmp(id, "culprit") != 0 &&
                    std::strcmp(id, "CULPRIT") != 0 &&
                    std::strcmp(id, "POSITIVE_TRAIN_SITE") != 0 &&
                    std::strcmp(id, "NEGATIVE_TRAIN_SITE") != 0)
                    stale_info.emplace_back(id);
            }
            for (const auto& id : stale_info)
                if (bcf_update_info(output_header, record, id.c_str(), nullptr, 0, BCF_HT_STR) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear recalibration INFO/" + id);
            if (bcf_update_info(output_header, record, "CULPRIT", nullptr, 0, BCF_HT_STR) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear legacy CULPRIT");
            const int32_t end_value = static_cast<int32_t>(record->pos + 1);
            if (bcf_update_info_int32(output_header, record, "END", &end_value, 1) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write recalibration END");
            if (bcf_update_alleles_str(output_header, record, "N,<VQSR>") != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write recalibration dummy alleles");
            const std::uint32_t missing_qual_bits = bcf_float_missing;
            std::memcpy(&record->qual, &missing_qual_bits, sizeof(record->qual));
            if (bcf_update_filter(output_header, record, nullptr, 0) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear recalibration FILTER");
        }
        // Filtering/scoring above consumed the full input record.  Only at
        // the writer boundary do we remove FORMAT/sample payload, matching
        // GATK's sites-only output contract without changing annotations or
        // model decisions.
        if (bcf_subset_format(output_header, record) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sites-only VCF record");
        if (bcf_translate(output_header, input_header, record) != 0 || bcf_write(output, output_header, record) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write recal VCF record");
        bcf_clear(record);
    }
    bcf_destroy(record); bcf_close(output); bcf_hdr_destroy(output_header); bcf_hdr_destroy(input_header); bcf_close(input);
    if (!complete_file(options.output)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete recal VCF");
    std::string index_path;
    if (options.create_output_variant_index && options.output != "-") {
        if (suffix(options.output, ".gz")) {
            index_path = options.output + ".tbi";
            if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0 || !complete_file(index_path))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build recal VCF index");
        } else {
            index_path = options.output + ".idx";
            fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
        }
    }
    (void)input_records;
    return index_path;
}

#endif

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        auto options = parse(argc, argv);
#if !FASTGATK_HAS_HTSLIB
        (void)options;
        throw std::runtime_error("BACKEND_UNAVAILABLE: VariantRecalibrator requires HTSlib");
#else
        const bool input_model_requested = !options.input_model.empty();
        ModelBundle serialized_model;
        if (input_model_requested) {
            const auto requested = options.annotations_explicit ? options.annotations : std::vector<std::string>{};
            serialized_model = read_model_report(options.input_model, requested,
                                                 options.annotations_explicit,
                                                 /*reorder_to_requested=*/false);
            options.annotations = serialized_model.annotations;
        }
        const auto resource_sets = build_resource_sets(options.resources);
        std::uint64_t input_records = 0;
        JavaRandom gatk_random(47382911ULL);
        auto entries = read_entries(options, resource_sets, input_records, &gatk_random);
        std::size_t scoreable = 0;
        std::size_t missing_annotation_records = 0;
        for (const auto& entry : entries) {
            if (!entry.scoreable) continue;
            ++scoreable;
            if (std::any_of(entry.missing.begin(), entry.missing.end(),
                            [](std::uint8_t value) { return value != 0; }))
                ++missing_annotation_records;
        }
        if (scoreable == 0 && !input_model_requested) {
            if (options.allele_specific)
                throw std::runtime_error("BACKEND_UNAVAILABLE: --AS requires every requested annotation to be Number=A");
            throw std::runtime_error("BAD_INPUT: no records have all requested annotations");
        }
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(static_cast<unsigned>(options.threads));
        Kokkos::initialize(settings);
        initialized = true;
        KernelTelemetry kernel_telemetry;
        GaussianMixture good;
        GaussianMixture bad;
        bool model_loaded = false;
        std::size_t bad_cutoff_selected = 0;
        std::size_t bad_fallback_selected = 0;
        std::vector<double> normalization_means;
        std::vector<double> normalization_stdevs;
        if (input_model_requested) {
            normalize_entries(entries, serialized_model, &gatk_random);
            mark_failing_std_threshold(entries, options.standard_deviation_threshold);
            scoreable = 0;
            missing_annotation_records = 0;
            for (const auto& entry : entries) {
                if (!entry.scoreable) continue;
                ++scoreable;
                if (std::any_of(entry.missing.begin(), entry.missing.end(),
                                [](std::uint8_t value) { return value != 0; }))
                    ++missing_annotation_records;
            }
            if (scoreable == 0)
                throw std::runtime_error("BAD_INPUT: no input variants match the serialized VQSR model annotations");
            good = serialized_model.good;
            bad = serialized_model.bad;
            normalization_means = serialized_model.annotation_means;
            normalization_stdevs = serialized_model.annotation_stdevs;
            model_loaded = true;
        } else {
            compute_annotation_normalization(entries, options.annotations.size(),
                                             normalization_means, normalization_stdevs);
            const auto annotation_order = gatk_annotation_order(
                entries, normalization_means, options.annotations.size());
            reorder_annotation_dimensions(entries, options.annotations,
                                          normalization_means, normalization_stdevs,
                                          annotation_order);
            normalize_entries_in_place(entries, normalization_means, normalization_stdevs,
                                       &gatk_random);
            mark_failing_std_threshold(entries, options.standard_deviation_threshold);
            // GATK's --max-attempts is a retry budget around model
            // construction, not an alias for --max-iterations.  A failed
            // initialization consumes the Java random stream and retries
            // with the next state; a successful first attempt is therefore
            // byte-for-byte unchanged from the previous native path.
            auto fit_with_retries = [&](bool training, std::size_t requested_components) {
                std::exception_ptr failure;
                for (std::size_t attempt = 0; attempt < options.max_attempts; ++attempt) {
                    try {
                        return fit_mixture(entries, training, options.annotations.size(),
                                           requested_components, options.max_iterations,
                                           options.full_covariance, options.shrinkage,
                                           options.dirichlet, options.prior_counts,
                                           options.standard_deviation_threshold,
                                           options.kmeans_iterations, options.maximum_training_variants,
                                           &gatk_random, &kernel_telemetry);
                    } catch (...) {
                        failure = std::current_exception();
                    }
                }
                if (failure) std::rethrow_exception(failure);
                throw std::runtime_error("MODEL_FAILURE: no model-build attempts were available");
            };
            good = fit_with_retries(true, options.max_gaussians);
            const auto bad_selection = select_worst_variants(
                entries, good, options.bad_lod_score_cutoff, options.minimum_bad_variants);
            bad_cutoff_selected = bad_selection.cutoff_selected;
            bad_fallback_selected = bad_selection.fallback_selected;
            bad = fit_with_retries(false,
                                   std::min(options.max_gaussians, options.max_negative_gaussians));
            if (good.count == 0)
                throw std::runtime_error("BAD_INPUT: training resources do not overlap scoreable input variants");
            if (bad.count == 0) {
                bad = good;
                for (auto& component : bad.variances)
                    for (auto& variance : component) variance *= 4.0;
                for (auto& covariance : bad.covariances)
                    for (auto& value : covariance) value *= 4.0;
                refresh_precision(bad);
            }
        }
        if (good.training_execution_space.empty())
            good.training_execution_space = Kokkos::DefaultExecutionSpace::name();
        if (bad.training_execution_space.empty())
            bad.training_execution_space = Kokkos::DefaultExecutionSpace::name();
        if (good.training_execution_policy.empty()) good.training_execution_policy = "SerializedModel";
        if (bad.training_execution_policy.empty()) bad.training_execution_policy = "SerializedModel";
        good.scoring_execution_space = Kokkos::DefaultExecutionSpace::name();
        bad.scoring_execution_space = Kokkos::DefaultExecutionSpace::name();
        good.scoring_execution_policy = "RangePolicy";
        bad.scoring_execution_policy = "RangePolicy";
        score_entries(entries, good, bad, options.threads, &kernel_telemetry, &gatk_random);
        // GATK's contrastive evaluation adds the per-datum prior log-odds to
        // the positive-minus-negative model score.  It is not part of the
        // serialized Gaussian model, so it must be re-applied on both fresh
        // and --input-model runs.
        double minimum_prior_lod = std::numeric_limits<double>::infinity();
        double maximum_prior_lod = -std::numeric_limits<double>::infinity();
        for (auto& entry : entries) {
            if (!entry.scoreable || !std::isfinite(entry.score)) continue;
            entry.score += entry.prior_lod;
            if (std::isfinite(entry.prior_lod)) {
                minimum_prior_lod = std::min(minimum_prior_lod, entry.prior_lod);
                maximum_prior_lod = std::max(maximum_prior_lod, entry.prior_lod);
            }
        }
        // GATK VariantDataManager writes the NAME of the per-datum worst
        // annotation dimension as `culprit` (VariantDataManager.java:485),
        // never a model/provenance name.  The dimension itself comes from
        // VariantRecalibratorEngine.calculateWorstPerformingAnnotation, which
        // runs on every datum after the contrastive evaluation.
        calculate_worst_annotation(entries, good, bad);
        std::unordered_map<std::string, std::string> culprits;
        for (const auto& entry : entries) {
            if (entry.key.empty()) continue;
            const bool named = entry.worst_annotation >= 0 &&
                static_cast<std::size_t>(entry.worst_annotation) < options.annotations.size();
            culprits[entry.key] = named
                ? options.annotations[static_cast<std::size_t>(entry.worst_annotation)]
                : "NULL";
        }
        const std::string model = model_loaded ? "serialized-gmm" : (options.full_covariance
            ? "full-covariance-gmm"
            : (options.max_gaussians > 1 ? "diagonal-gmm" : "diagonal-gaussian"));
        std::unordered_map<std::string, double> scores;
        std::unordered_map<std::string, bool> positive_training;
        std::unordered_map<std::string, bool> negative_training;
        for (const auto& entry : entries)
            if (!entry.key.empty()) {
                if (std::isfinite(entry.score)) scores[entry.key] = entry.score;
                // Keep labels independent of scoreability: GATK records the
                // resource membership in the recal VCF whenever the allele
                // is part of a training set, even if a requested annotation
                // was missing and therefore did not receive VQSLOD.
                positive_training[entry.key] = positive_training[entry.key] || entry.training;
                // select_worst_variants marks the rows chosen for the
                // negative model as bad_candidate.  GATK writes those
                // anti-training rows back as NEGATIVE_TRAIN_SITE in the
                // recalibration table; retaining only an input-side flag
                // silently lost that provenance.
                negative_training[entry.key] = negative_training[entry.key] ||
                    entry.negative_training || entry.bad_candidate;
            }
        std::uint64_t scored_records = 0;
        const auto vcf_index = write_recal_vcf(options, scores, positive_training,
                                               negative_training, culprits, model, input_records,
                                               scored_records, options.allele_specific);
        write_tranches(options.tranches_output, entries, options.truth_sensitivity, options.mode,
                       options.output_tranches_for_scatter, options.vqslod_tranches);
        if (!options.output_model.empty())
            write_model_report(options.output_model, good, bad, options.annotations,
                               &normalization_means, &normalization_stdevs);
        if (!options.rscript_output.empty()) {
            std::ofstream script(options.rscript_output);
            if (!script) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write rscript file");
            script << "# fastgatk VariantRecalibrator prototype; model=" << model << "; annotations=";
            for (std::size_t index = 0; index < options.annotations.size(); ++index) {
                if (index) script << ',';
                script << options.annotations[index];
            }
            script << '\n';
            script.close();
            if (!complete_file(options.rscript_output)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete rscript");
        }
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest");
        manifest << "{\"schema_version\":1,\"tool\":\"VariantRecalibrator\",\"implementation\":\"fastgatk-variant-recalibrator\","
                 << "\"status\":\"prototype\",\"model\":\"" << model << "\",\"components\":"
                 << good.weights.size() << ",\"max_attempts\":" << options.max_attempts
                 << ",\"em_iterations\":" << options.max_iterations
                 << ",\"sample_every_nth_variant\":" << options.sample_every
                 << ",\"maximum_training_variants\":" << options.maximum_training_variants
                 << ",\"minimum_bad_variants\":" << options.minimum_bad_variants
                 << ",\"bad_lod_score_cutoff\":" << options.bad_lod_score_cutoff
                 << ",\"bad_cutoff_selected\":" << bad_cutoff_selected
                 << ",\"bad_fallback_selected\":" << bad_fallback_selected
                 << ",\"scatter_tranches\":" << (options.output_tranches_for_scatter ? "true" : "false")
                 << ",\"vqslod_tranches\":" << options.vqslod_tranches.size()
                 << ",\"max_negative_gaussians\":" << options.max_negative_gaussians
                 << ",\"kmeans_iterations\":" << options.kmeans_iterations
                 << ",\"shrinkage\":" << std::setprecision(17) << options.shrinkage
                 << ",\"dirichlet\":" << options.dirichlet
                 << ",\"prior_counts\":" << options.prior_counts
                 << ",\"standard_deviation_threshold\":" << options.standard_deviation_threshold
                 << ",\"annotation_jitter\":true,\"gatk_random_seed\":47382911"
                 << ",\"mq_cap\":" << options.mq_cap
                 << ",\"mq_jitter\":" << options.mq_jitter
                 << ",\"variational_normal_wishart\":" << (model_loaded ? "false" : "true")
                 << ",\"full_covariance\":" << (options.full_covariance ? "true" : "false")
                 << ",\"allele_specific\":" << (options.allele_specific ? "true" : "false")
                 << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"mode\":\""
                 << json_escape(options.mode) << "\",\"input_records\":" << input_records
                 << ",\"scoreable_records\":" << scoreable << ",\"scored_records\":" << scored_records
                 << ",\"missing_annotation_records\":" << missing_annotation_records
                 << ",\"annotations\":[";
        for (std::size_t index = 0; index < options.annotations.size(); ++index) {
            if (index) manifest << ',';
            manifest << '"' << json_escape(options.annotations[index]) << '"';
        }
        manifest << "],\"training_records\":" << good.count << ",\"background_records\":" << bad.count
                 << ",\"resources\":" << options.resources.size() << ",\"output\":\""
                 << json_escape(options.output) << "\",\"tranches\":\"" << json_escape(options.tranches_output)
                 << "\",\"input_model\":\"" << json_escape(options.input_model)
                 << "\",\"output_model\":\"" << json_escape(options.output_model)
                 << "\",\"model_loaded\":" << (model_loaded ? "true" : "false")
                 << ",\"vcf_index\":\"" << json_escape(vcf_index)
                 << "\",\"vcf_index_enabled\":" << (!vcf_index.empty() ? "true" : "false")
                 << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                 << ",\"training_execution_space\":\""
                 << json_escape(good.training_execution_space)
                 << "\",\"training_execution_policy\":\""
                 << json_escape(good.training_execution_policy)
                 << "\",\"scoring_execution_space\":\""
                 << json_escape(good.scoring_execution_space)
                 << "\",\"scoring_execution_policy\":\""
                 << json_escape(good.scoring_execution_policy) << "\""
                 << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"kernel_execution_policy\":\"TeamPolicy+RangePolicy\""
                 << ",\"kernel_batches\":" << kernel_telemetry.batches
                 << ",\"kernel_observations\":" << kernel_telemetry.observations
                 << ",\"kernel_prepare_seconds\":" << kernel_telemetry.prepare_seconds
                 << ",\"kernel_execute_seconds\":" << kernel_telemetry.execute_seconds
                 << ",\"resource_prior\":{\"default_quality\":2.0,\"log10_odds_min\":"
                 << std::setprecision(17) << (std::isfinite(minimum_prior_lod) ? minimum_prior_lod : 0.0)
                 << ",\"log10_odds_max\":"
                 << (std::isfinite(maximum_prior_lod) ? maximum_prior_lod : 0.0) << "}"
                 << ",\"normalization\":{\"means\":[";
        for (std::size_t index = 0; index < normalization_means.size(); ++index) {
            if (index) manifest << ',';
            manifest << std::setprecision(17) << normalization_means[index];
        }
        manifest << "],\"stdevs\":[";
        for (std::size_t index = 0; index < normalization_stdevs.size(); ++index) {
            if (index) manifest << ',';
            manifest << std::setprecision(17) << normalization_stdevs[index];
        }
        manifest << "]}"
                 << ",\"bit_identical_to_gatk\":false}\n";
        manifest.close();
        if (!manifest || !complete_file(manifest_path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete manifest");
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"tool\":\"VariantRecalibrator\",\"status\":\"prototype\",\"model\":\"" << model
                  << "\",\"components\":" << good.weights.size()
                  << ",\"max_negative_gaussians\":" << options.max_negative_gaussians
                  << ",\"kmeans_iterations\":" << options.kmeans_iterations
                  << ",\"variational_normal_wishart\":" << (model_loaded ? "false" : "true")
                  << ",\"full_covariance\":" << (options.full_covariance ? "true" : "false")
                  << ",\"allele_specific\":" << (options.allele_specific ? "true" : "false")
                  << ",\"sites_only_vcf_output\":" << (options.sites_only_vcf_output ? "true" : "false")
                  << ",\"scatter_tranches\":"
                  << (options.output_tranches_for_scatter ? "true" : "false")
                  << ",\"vqslod_tranches\":" << options.vqslod_tranches.size()
                  << ",\"input_records\":"
                  << input_records << ",\"scoreable_records\":" << scoreable << ",\"scored_records\":" << scored_records
                  << ",\"missing_annotation_records\":" << missing_annotation_records
                  << ",\"annotation_jitter\":true,\"gatk_random_seed\":47382911"
                  << ",\"mq_cap\":" << options.mq_cap
                  << ",\"mq_jitter\":" << options.mq_jitter
                  << ",\"input_model\":" << (model_loaded ? "true" : "false")
                  << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                  << ",\"training_execution_policy\":\""
                  << good.training_execution_policy
                  << "\",\"scoring_execution_policy\":\""
                  << good.scoring_execution_policy << "\"}\n";
        return 0;
#endif
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << error.what() << '\n';
        return 2;
    }
}
