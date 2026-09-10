#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::vector<std::string> evals;
    std::vector<std::string> comps;
    std::string output;
   std::string manifest;
   std::string reference;
    std::string pedigree;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> modules;
    std::vector<std::string> stratifications;
    bool list = false;
    bool no_standard_modules = false;
    bool no_standard_stratifications = false;
    bool strict_allele_match = false;
   bool gatk_report = false;
    // GATK's default treats genotyped AC=0 sites as reference for
    // polymorphism-tracking evaluators. --keep-ac0 promotes them back to
    // their site-level SNP/INDEL type, including per-sample/family strata.
    bool keep_ac0 = false;
    double mendelian_qual_threshold = 50.0;
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

bool has_option_value(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-variant-eval (GATK-compatible native prototype)\n"
                         "  -eval, --eval FILE              evaluation VCF/BCF (repeatable)\n"
                         "  -comp, --comp FILE              comparison VCF/BCF (repeatable)\n"
                         "  -O, --output FILE               VariantEval report\n"
                        "  -R, --reference FILE            reference compatibility option\n"
                         "      --pedigree, -ped FILE       pedigree for MendelianViolationEvaluator\n"
                         "      --mendelian-violation-qual-threshold, -mvq FLOAT\n"
                         "                                  minimum trio GQ (default 50)\n"
                         "  -L, --intervals REGION           interval selector (repeatable)\n"
                         "      --interval-set-rule RULE     UNION (default) or INTERSECTION\n"
                         "  -EV, --eval-module NAME          evaluator module (repeatable)\n"
                         "  -S, --select NAME                stratification (repeatable)\n"
                         "  -ST, --stratification-module NAME GATK stratification alias\n"
                         "  -ls, --list                      list supported modules\n"
                         "      --keep-ac0, -keep-ac0         retain genotyped AC=0 sites\n"
                         "      --gatk-report                emit GATKReport v1.1-compatible tables\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-ls" || argument == "--list") options.list = true;
        else if (argument == "--gatk-report" || argument == "--gatk-compatible-report")
            options.gatk_report = true;
        else if (argument == "-eval" || has_option_value(argument, "--eval"))
            options.evals.push_back(require_value(index, argc, argv, argument, "--eval", "-eval"));
        else if (argument == "-comp" || has_option_value(argument, "--comp"))
            options.comps.push_back(require_value(index, argc, argv, argument, "--comp", "-comp"));
        else if (argument == "-O" || has_option_value(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
       else if (argument == "-R" || has_option_value(argument, "--reference"))
           options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "--pedigree" || argument == "--ped" || has_option_value(argument, "--pedigree") || has_option_value(argument, "--ped"))
            options.pedigree = require_value(index, argc, argv, argument,
                argument.rfind("--pedigree", 0) == 0 ? "--pedigree" : "--ped", argument == "--ped" ? "-ped" : "--pedigree");
        else if (argument == "-mvq" || argument == "--mendelian-violation-qual-threshold" ||
                 has_option_value(argument, "--mendelian-violation-qual-threshold") ||
                 has_option_value(argument, "-mvq")) {
            const auto inline_short = option_value(argument, "-mvq");
            options.mendelian_qual_threshold = inline_short.empty()
                ? std::stod(require_value(index, argc, argv, argument,
                    "--mendelian-violation-qual-threshold", "-mvq"))
                : std::stod(inline_short);
        }
        else if (argument == "-L" || has_option_value(argument, "--intervals") ||
                 has_option_value(argument, "--interval") || has_option_value(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                has_option_value(argument, "--region") ? "--region" :
                has_option_value(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (has_option_value(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-EV" || has_option_value(argument, "--eval-module"))
            options.modules.push_back(require_value(index, argc, argv, argument, "--eval-module", "-EV"));
        else if (argument == "-S" || has_option_value(argument, "--select") ||
                 argument == "-ST" || has_option_value(argument, "--stratification-module")) {
            const char* name = (argument == "-ST" || has_option_value(argument, "--stratification-module"))
                ? "--stratification-module" : "--select";
            options.stratifications.push_back(require_value(
                index, argc, argv, argument, name, argument == "-ST" ? "-ST" : "-S"));
        }
        else if (argument == "-no-ev" || argument == "--do-not-use-all-standard-modules")
            options.no_standard_modules = true;
        else if (argument == "-no-st" || argument == "--do-not-use-all-standard-stratifications")
            options.no_standard_stratifications = true;
        else if (argument == "-strict" || argument == "--require-strict-allele-match")
            options.strict_allele_match = true;
        else if (argument == "--keep-ac0" || argument == "-keep-ac0" ||
                 argument.rfind("--keep-ac0=", 0) == 0 ||
                 argument.rfind("-keep-ac0=", 0) == 0) {
            const char* name = argument.rfind("-keep-ac0", 0) == 0
                ? "-keep-ac0" : "--keep-ac0";
            options.keep_ac0 = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        }
        else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Accepted compatibility flags.  The native evaluator is deterministic
            // and does not need a dictionary validation pass.
        } else if (has_option_value(argument, "--java-options") ||
                   has_option_value(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else if (has_option_value(argument, "--output-manifest") ||
                   has_option_value(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--known-name" || argument == "--gold-standard" ||
                   argument == "--sample" ||
                   argument == "--stratification-module" || argument == "--strat-intervals" ||
                   argument == "--num-samples" || argument == "--sample-ploidy") {
            // These controls are accepted only to preserve a predictable
            // compatibility surface; unsupported semantics are reported in the
            // manifest instead of being silently treated as native.
            if (index + 1 >= argc) throw std::invalid_argument("missing value for " + argument);
            ++index;
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (!options.list && options.evals.empty()) throw std::invalid_argument("-eval/--eval is required");
    if (!options.list && options.output.empty()) throw std::invalid_argument("-O/--output is required");
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct GenotypeCounts {
    std::uint64_t called = 0;
    std::uint64_t no_call = 0;
    std::uint64_t hom_ref = 0;
    std::uint64_t het = 0;
    std::uint64_t hom_var = 0;
    std::uint64_t chromosomes = 0;
    std::uint64_t alt_chromosomes = 0;
    std::uint64_t alt_carrier_samples = 0;
    std::uint64_t called_not_filtered = 0;
    std::uint64_t no_call_or_filtered = 0;
    std::uint64_t filtered = 0;
    bool singleton = false;
    bool polymorphic = false;
    std::vector<int> classes;
    std::vector<int> gq;
    std::vector<int> alt_chromosomes_per_sample;
    // Called allele indices per sample, retained for evaluators (for example
    // ThetaVariantEvaluator) that need pairwise allele differences rather
    // than only the hom-ref/het/hom-var class.
    std::vector<std::vector<int>> alleles;
};

struct Record {
    std::string contig;
    int rid = -1;
    int position = 0;
    int end_exclusive = 0;
    std::string ref;
    std::vector<std::string> alts;
    bool filtered = false;
    bool symbolic = false;
    bool snp = false;
    bool mnp = false;
    bool insertion = false;
    bool deletion = false;
    bool complex_indel = false;
    bool mixed = false;
    bool variant = false;
    bool transition = false;
    GenotypeCounts genotypes;
    std::string filter = "PASS";
    std::string variant_type = "NO_VARIATION";
    double allele_frequency = 0.0;
    bool has_allele_count = false;
    std::int64_t allele_count = 0;
    std::vector<std::string> sample_names;
};

// GATK's CountVariants evaluator treats a record with FORMAT/GT differently
// from a sites-only record.  For genotyped input, an ALT allele present only
// in the site representation is not a variant locus unless at least one
// sample has a called ALT genotype; all-hom-ref and all-no-call records are
// reference loci.  Sites-only records retain the site-level VCF semantics.
bool has_effective_variant_genotype(const Record& record, bool keep_ac0 = false) {
    return record.genotypes.classes.empty() || keep_ac0
        ? record.variant
        : record.genotypes.polymorphic;
}

// CountVariants is evaluated once per sample when the GATK Sample
// stratification is enabled.  Keep this separate from the site-level
// StratumMetrics so a sample's hom-ref/no-call genotype does not get mistaken
// for a site-level variant merely because another sample carries the ALT.
struct SampleCountMetrics {
    std::uint64_t called_loci = 0;
    std::uint64_t reference_loci = 0;
    std::uint64_t variant_loci = 0;
    std::uint64_t snps = 0;
    std::uint64_t mnps = 0;
    std::uint64_t insertions = 0;
    std::uint64_t deletions = 0;
    std::uint64_t complex = 0;
    std::uint64_t symbolic = 0;
    std::uint64_t mixed = 0;
    std::uint64_t no_calls = 0;
    std::uint64_t het = 0;
    std::uint64_t hom_ref = 0;
    std::uint64_t hom_var = 0;
    std::uint64_t singletons = 0;
};

struct ReadStats {
    std::size_t records = 0;
    fastgatk::io::IntervalFileStats interval_files;
};

struct HeaderShape {
    std::vector<std::string> contigs;
    std::vector<std::int64_t> contig_lengths;
    std::vector<std::string> samples;
};

struct Trio {
    std::string family;
    std::string mother;
    std::string father;
    std::string child;
};

struct MendelianMetrics {
    std::uint64_t n_variants = 0;
    std::uint64_t n_skipped = 0;
    std::uint64_t n_fam_called = 0;
    std::uint64_t n_var_fam_called = 0;
    std::uint64_t n_low_qual = 0;
    std::uint64_t n_no_call = 0;
    std::uint64_t n_loci_violations = 0;
    std::uint64_t n_violations = 0;
    std::uint64_t mv_ref_ref_var = 0;
    std::uint64_t mv_ref_ref_het = 0;
    std::uint64_t mv_ref_het_var = 0;
    std::uint64_t mv_ref_var_var = 0;
    std::uint64_t mv_ref_var_ref = 0;
    std::uint64_t mv_var_het_ref = 0;
    std::uint64_t mv_var_var_ref = 0;
    std::uint64_t mv_var_var_het = 0;
    std::uint64_t hom_ref_hom_ref_hom_ref = 0;
    std::uint64_t het_het_het = 0;
    std::uint64_t het_het_hom_ref = 0;
    std::uint64_t het_het_hom_var = 0;
    std::uint64_t hom_var_hom_var_hom_var = 0;
    std::uint64_t hom_ref_hom_var_het = 0;
    std::uint64_t het_het_inherited_ref = 0;
    std::uint64_t het_het_inherited_var = 0;
    std::uint64_t hom_ref_het_inherited_ref = 0;
    std::uint64_t hom_ref_het_inherited_var = 0;
    std::uint64_t hom_var_het_inherited_ref = 0;
    std::uint64_t hom_var_het_inherited_var = 0;
};

enum class ValidationSiteStatus : int {
    NoCall = 0,
    Filtered = 1,
    Mono = 2,
    Poly = 3,
};

struct ValidationReportMetrics {
    static constexpr std::size_t status_count = 4;
    std::array<std::uint64_t, status_count * status_count> counts{};

    std::uint64_t get(ValidationSiteStatus comp, ValidationSiteStatus eval) const {
        return counts[static_cast<std::size_t>(comp) * status_count +
                      static_cast<std::size_t>(eval)];
    }
};

struct Metrics {
    std::uint64_t processed_loci = 0;
    std::uint64_t called_loci = 0;
    std::uint64_t reference_loci = 0;
    std::uint64_t variant_loci = 0;
    std::uint64_t snps = 0;
    std::uint64_t mnps = 0;
    std::uint64_t insertions = 0;
    std::uint64_t deletions = 0;
    std::uint64_t complex = 0;
    std::uint64_t symbolic = 0;
    std::uint64_t mixed = 0;
    std::uint64_t no_calls = 0;
    std::uint64_t hom_ref = 0;
    std::uint64_t het = 0;
    std::uint64_t hom_var = 0;
    std::uint64_t ti = 0;
    std::uint64_t tv = 0;
    std::uint64_t comp_records = 0;
    std::uint64_t comp_ti = 0;
    std::uint64_t comp_tv = 0;
    std::uint64_t eval_comp_overlap = 0;
    std::uint64_t strict_allele_mismatch = 0;
    std::uint64_t missing_comp_snps = 0;
    std::uint64_t compared_eval_genotypes = 0;
    std::uint64_t concordant_genotypes = 0;
    std::map<std::string, std::uint64_t> genotype_concordance;
    // Standard GATK IndelSummary/MultiallelicSummary aggregates.  The
    // existing CountVariants site counters remain unchanged; these fields
    // deliberately count both sites and concrete alternate alleles where the
    // Java evaluator does so.
    std::uint64_t snp_sites = 0;
    std::uint64_t multiallelic_snp_sites = 0;
    std::uint64_t indel_sites = 0;
    std::uint64_t multiallelic_indel_sites = 0;
    std::uint64_t indel_alleles = 0;
    std::uint64_t indel_insertions = 0;
    std::uint64_t indel_deletions = 0;
    std::uint64_t large_insertions = 0;
    std::uint64_t large_deletions = 0;
    std::uint64_t indel_singletons = 0;
    std::uint64_t snp_singletons = 0;
    std::uint64_t multiallelic_snp_ti = 0;
    std::uint64_t multiallelic_snp_tv = 0;
    std::uint64_t known_multiallelic_snp_partial = 0;
    std::uint64_t known_multiallelic_snp_complete = 0;
    std::map<int, std::uint64_t> indel_length_histogram;
    std::uint64_t indel_histogram_total = 0;
    std::uint64_t genotype_called_not_filtered = 0;
    std::uint64_t genotype_no_call_or_filtered = 0;
    std::uint64_t filtered_sites = 0;
    // VariantAFEvaluator (GATK's diploid SNP allele-fraction summary).
    double variant_af_sum = 0.0;
    std::uint64_t variant_af_called_sites = 0;
    std::uint64_t variant_af_het_sites = 0;
    std::uint64_t variant_af_hom_var_sites = 0;
    std::uint64_t variant_af_hom_ref_sites = 0;
    double theta_total_het = 0.0;
    double theta_total_avg_diffs = 0.0;
    double theta_region_num_sites = 0.0;
    std::uint64_t theta_num_sites = 0;
    MendelianMetrics mendelian;
    ValidationReportMetrics validation_report;
};

struct StratumMetrics {
    std::uint64_t loci = 0;
    std::uint64_t variants = 0;
    std::uint64_t snps = 0;
    std::uint64_t mnps = 0;
    std::uint64_t insertions = 0;
    std::uint64_t deletions = 0;
    std::uint64_t complex = 0;
    std::uint64_t symbolic = 0;
    std::uint64_t mixed = 0;
    std::uint64_t ti = 0;
    std::uint64_t tv = 0;
};

bool is_base_snp(const std::string& ref, const std::string& alt) {
    return ref.size() == 1 && alt.size() == 1 && alt[0] != '<' && alt[0] != '*';
}

bool is_transition_base(char ref, char alt) {
    ref = static_cast<char>(std::toupper(static_cast<unsigned char>(ref)));
    alt = static_cast<char>(std::toupper(static_cast<unsigned char>(alt)));
    return (ref == 'A' && alt == 'G') || (ref == 'G' && alt == 'A') ||
           (ref == 'C' && alt == 'T') || (ref == 'T' && alt == 'C');
}

bool is_variant_key_equal(const Record& lhs, const Record& rhs, bool strict) {
    if (lhs.contig != rhs.contig || lhs.position != rhs.position) return false;
    if (!strict) return true;
    return lhs.ref == rhs.ref && lhs.alts == rhs.alts;
}

std::string record_key(const Record& record) {
    std::ostringstream key;
    key << record.contig << ':' << record.position << ':' << record.ref << ':';
    for (const auto& alt : record.alts) key << alt << ',';
    return key.str();
}

bool record_matches_comparison(const Record& eval,
                               const std::vector<Record>& comparisons,
                               bool strict_allele_match) {
    return std::any_of(comparisons.begin(), comparisons.end(), [&](const auto& comp) {
        return is_variant_key_equal(eval, comp, strict_allele_match);
    });
}

GenotypeCounts read_genotypes(const bcf_hdr_t* header, bcf1_t* record) {
    GenotypeCounts counts;
    if (bcf_hdr_nsamples(header) == 0) return counts;
    int32_t* raw = nullptr;
    int raw_count = 0;
    const int length = bcf_get_genotypes(header, record, &raw, &raw_count);
    if (length <= 0) {
        free(raw);
        return counts;
    }
    const int samples = bcf_hdr_nsamples(header);
    if (samples <= 0 || length % samples != 0)
        throw std::runtime_error("BAD_INPUT: malformed FORMAT/GT cardinality");
    const int ploidy = length / samples;
    if (ploidy <= 0)
        throw std::runtime_error("BAD_INPUT: FORMAT/GT has zero ploidy");
    counts.classes.assign(static_cast<std::size_t>(samples), 0);
    counts.gq.assign(static_cast<std::size_t>(samples), -1);
    counts.alt_chromosomes_per_sample.assign(static_cast<std::size_t>(samples), 0);
    counts.alleles.assign(static_cast<std::size_t>(samples), {});
    int32_t* gq_raw = nullptr;
    int gq_raw_count = 0;
    const int gq_length = bcf_get_format_int32(header, record, "GQ", &gq_raw, &gq_raw_count);
    if (gq_length > 0 && gq_raw != nullptr) {
        const auto available = std::min(gq_length, samples);
        for (int sample = 0; sample < available; ++sample) {
            const auto value = gq_raw[sample];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end)
                counts.gq[static_cast<std::size_t>(sample)] = value;
        }
    }
    free(gq_raw);
    char** filter_values = nullptr;
    int filter_count = 0;
    const int filter_length = bcf_get_format_string(header, record, "FT", &filter_values, &filter_count);
    for (int sample = 0; sample < samples; ++sample) {
        bool missing = false;
        bool any_alt = false;
        int first_allele = -1;
        bool all_same = true;
        bool genotype_filtered = false;
        if (filter_length > sample && filter_values != nullptr && filter_values[sample] != nullptr) {
            const std::string filter(filter_values[sample]);
            genotype_filtered = !filter.empty() && filter != "." && filter != "PASS";
        }
        for (int allele = 0; allele < ploidy; ++allele) {
            const int32_t value = raw[sample * ploidy + allele];
            if (bcf_gt_is_missing(value) || value == bcf_int32_vector_end) {
                missing = true;
                continue;
            }
            const int index = bcf_gt_allele(value);
            counts.alleles[static_cast<std::size_t>(sample)].push_back(index);
            if (first_allele < 0) first_allele = index;
            else if (first_allele != index) all_same = false;
            if (index > 0) any_alt = true;
            if (index > 0)
                ++counts.alt_chromosomes_per_sample[static_cast<std::size_t>(sample)];
        }
        if (missing || first_allele < 0) {
            ++counts.no_call;
            ++counts.no_call_or_filtered;
            if (genotype_filtered) ++counts.filtered;
            counts.classes[static_cast<std::size_t>(sample)] = 0;
            continue;
        }
        if (genotype_filtered) {
            ++counts.filtered;
            ++counts.no_call_or_filtered;
        } else {
            ++counts.called_not_filtered;
        }
        ++counts.called;
        counts.chromosomes += static_cast<std::uint64_t>(ploidy);
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto value = raw[sample * ploidy + allele];
            if (!bcf_gt_is_missing(value) && value != bcf_int32_vector_end &&
                bcf_gt_allele(value) > 0)
                ++counts.alt_chromosomes;
        }
        if (any_alt) {
            counts.polymorphic = true;
            ++counts.alt_carrier_samples;
        }
        if (!any_alt) {
            ++counts.hom_ref;
            counts.classes[static_cast<std::size_t>(sample)] = 1;
        } else if (all_same) {
            ++counts.hom_var;
            counts.classes[static_cast<std::size_t>(sample)] = 3;
        } else {
            ++counts.het;
            counts.classes[static_cast<std::size_t>(sample)] = 2;
        }
    }
    if (filter_values != nullptr) {
        if (filter_values[0] != nullptr) free(filter_values[0]);
        free(filter_values);
    }
    counts.singleton = counts.alt_carrier_samples == 1;
    free(raw);
    return counts;
}

Record decode_record(const bcf_hdr_t* header, bcf1_t* raw) {
    if (raw == nullptr || raw->pos < 0 || raw->n_allele < 1 || raw->d.allele == nullptr ||
        raw->d.allele[0] == nullptr || raw->d.allele[0][0] == '\0')
        throw std::runtime_error("BAD_INPUT: malformed variant record coordinates/REF");
    Record record;
    record.rid = raw->rid;
    record.contig = raw->rid >= 0 ? bcf_hdr_id2name(header, raw->rid) : "*";
    record.position = raw->pos + 1;
    record.end_exclusive = static_cast<int>(raw->pos + std::max<hts_pos_t>(1, raw->rlen));
    int32_t* end_values = nullptr;
    int end_count = 0;
    const auto end_length = bcf_get_info_int32(header, raw, "END", &end_values, &end_count);
    if (end_length > 0 && end_count > 0 && end_values[0] > 0) {
        if (end_values[0] < record.position)
            throw std::runtime_error("BAD_INPUT: INFO/END precedes POS");
        record.end_exclusive = std::max(record.end_exclusive, end_values[0]);
    }
    free(end_values);
    if (raw->n_allele > 0) record.ref = raw->d.allele[0];
    for (int index = 1; index < raw->n_allele; ++index) {
        if (raw->d.allele[index] == nullptr || raw->d.allele[index][0] == '\0')
            throw std::runtime_error("BAD_INPUT: malformed variant record ALT");
        record.alts.emplace_back(raw->d.allele[index]);
    }
    record.sample_names.reserve(static_cast<std::size_t>(bcf_hdr_nsamples(header)));
    for (int sample = 0; sample < bcf_hdr_nsamples(header); ++sample) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, sample);
        record.sample_names.emplace_back(name == nullptr ? std::string{} : std::string(name));
    }
    record.filter = "PASS";
    if (raw->d.n_flt > 0 && raw->d.flt != nullptr) {
        std::ostringstream filter;
        for (int index = 0; index < raw->d.n_flt; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, raw->d.flt[index]);
            if (name == nullptr || std::string(name) == "PASS" || std::string(name) == ".") continue;
            if (filter.tellp() > 0) filter << ';';
            filter << name;
        }
        if (filter.tellp() > 0) record.filter = filter.str();
    }
    record.filtered = record.filter != "PASS";
    record.genotypes = read_genotypes(header, raw);

    bool has_snp = false;
    bool has_mnp = false;
    bool has_insertion = false;
    bool has_deletion = false;
    bool has_other_indel = false;
    for (const auto& alt : record.alts) {
        if (alt.empty() || alt[0] == '<' || alt[0] == '*') {
            record.symbolic = true;
            continue;
        }
        if (is_base_snp(record.ref, alt)) {
            has_snp = true;
            if (is_transition_base(record.ref[0], alt[0])) record.transition = true;
        } else if (record.ref.size() == alt.size() && record.ref.size() > 1) {
            has_mnp = true;
        } else if (alt.size() > record.ref.size()) {
            has_insertion = true;
        } else if (alt.size() < record.ref.size()) {
            has_deletion = true;
        } else {
            has_other_indel = true;
        }
    }
    const int concrete_classes = static_cast<int>(has_snp) + static_cast<int>(has_mnp) +
        static_cast<int>(has_insertion || has_deletion || has_other_indel) + static_cast<int>(record.symbolic);
    record.mixed = concrete_classes > 1;
    record.snp = has_snp && !record.mixed;
    record.mnp = has_mnp && !record.mixed;
    record.insertion = has_insertion && !has_deletion && !has_other_indel && !record.mixed;
    record.deletion = has_deletion && !has_insertion && !has_other_indel && !record.mixed;
    record.complex_indel = has_other_indel || (has_insertion && has_deletion);
    record.variant = !record.alts.empty() && concrete_classes > 0;
    if (!record.variant) record.variant_type = "NO_VARIATION";
    else if (record.mixed) record.variant_type = "MIXED";
    else if (record.symbolic) record.variant_type = "SYMBOLIC";
    else if (record.snp) record.variant_type = "SNP";
    else if (record.mnp) record.variant_type = "MNP";
    else if (record.insertion) record.variant_type = "INSERTION";
    else if (record.deletion) record.variant_type = "DELETION";
    else record.variant_type = "COMPLEX_INDEL";
    if (record.genotypes.chromosomes != 0)
        record.allele_frequency = static_cast<double>(record.genotypes.alt_chromosomes) /
                                  static_cast<double>(record.genotypes.chromosomes);
    else {
        float* allele_frequencies = nullptr;
        int allele_frequency_count = 0;
        const auto length = bcf_get_info_float(header, raw, "AF",
                                               &allele_frequencies, &allele_frequency_count);
        if (length > 0 && allele_frequency_count > 0 &&
            !bcf_float_is_missing(allele_frequencies[0]) &&
            !bcf_float_is_vector_end(allele_frequencies[0]) &&
            std::isfinite(allele_frequencies[0]))
            record.allele_frequency = std::clamp(
                static_cast<double>(allele_frequencies[0]), 0.0, 1.0);
        free(allele_frequencies);
    }
    int32_t* allele_counts = nullptr;
    int allele_count_capacity = 0;
    const auto allele_count_length = bcf_get_info_int32(
        header, raw, "AC", &allele_counts, &allele_count_capacity);
    if (allele_count_length > 0 && allele_counts != nullptr) {
        record.has_allele_count = true;
        for (int index = 0; index < allele_count_length; ++index) {
            const auto value = allele_counts[index];
            if (value != bcf_int32_missing && value != bcf_int32_vector_end)
                record.allele_count += value;
        }
    }
    free(allele_counts);
    return record;
}

ValidationSiteStatus validation_site_status(
    const Record* record, const std::vector<std::string>* eval_samples = nullptr) {
    if (record == nullptr) return ValidationSiteStatus::NoCall;
    if (record->filtered) return ValidationSiteStatus::Filtered;

    if (!record->genotypes.classes.empty()) {
        bool polymorphic = false;
        bool used_subset = false;
        if (eval_samples != nullptr && !eval_samples->empty()) {
            used_subset = std::all_of(eval_samples->begin(), eval_samples->end(),
                [&](const std::string& sample) {
                    return std::find(record->sample_names.begin(), record->sample_names.end(), sample) !=
                           record->sample_names.end();
                });
            if (used_subset) {
                for (const auto& sample : *eval_samples) {
                    const auto found = std::find(
                        record->sample_names.begin(), record->sample_names.end(), sample);
                    const auto index = static_cast<std::size_t>(
                        std::distance(record->sample_names.begin(), found));
                    const auto genotype_class = record->genotypes.classes[index];
                    if (genotype_class == 2 || genotype_class == 3) {
                        polymorphic = true;
                        break;
                    }
                }
            }
        }
        if (!used_subset) polymorphic = record->genotypes.polymorphic;
        return polymorphic ? ValidationSiteStatus::Poly : ValidationSiteStatus::Mono;
    }

    // ValidationReport.calcSiteStatus treats a multi-allelic sites-only
    // record as polymorphic regardless of AC, then consults AC for a
    // biallelic record.  Without either GT or AC, GATK deliberately treats a
    // present record as called/polymorphic.
    if (record->alts.size() > 1) return ValidationSiteStatus::Poly;
    if (record->has_allele_count)
        return record->allele_count > 0 ? ValidationSiteStatus::Poly
                                        : ValidationSiteStatus::Mono;
    return ValidationSiteStatus::Poly;
}

std::string validation_locus_key(const Record& record) {
    return record.contig + '\t' + std::to_string(record.position);
}

ValidationReportMetrics compute_validation_report(
    const std::vector<Record>& eval_records,
    const std::vector<Record>& comp_records,
    const std::vector<std::string>& eval_samples) {
    std::map<std::string, const Record*> eval_by_locus;
    for (const auto& record : eval_records) {
        // Without the Filter stratifier, VariantEvalEngine.bindVariantContexts
        // removes filtered records before ValidationReport.update2 sees them.
        // This bounded evaluator deliberately mirrors that default path.
        if (record.filtered) continue;
        const auto [iterator, inserted] = eval_by_locus.emplace(
            validation_locus_key(record), &record);
        (void)iterator;
        if (!inserted)
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: ValidationReport requires at most one eval record per locus");
    }

    std::vector<int> encoded_counts;
    encoded_counts.reserve(comp_records.size());
    for (const auto& comp : comp_records) {
        if (comp.filtered) continue;
        const auto found = eval_by_locus.find(validation_locus_key(comp));
        const Record* eval = found == eval_by_locus.end() ? nullptr : found->second;
        const auto comp_status = validation_site_status(&comp, &eval_samples);
        const auto eval_status = validation_site_status(eval);
        encoded_counts.push_back(static_cast<int>(comp_status) *
            static_cast<int>(ValidationReportMetrics::status_count) +
            static_cast<int>(eval_status));
    }

    Kokkos::View<int*> encoded("variant_eval_validation_status", encoded_counts.size());
    auto encoded_host = Kokkos::create_mirror_view(encoded);
    for (std::size_t index = 0; index < encoded_counts.size(); ++index)
        encoded_host(index) = encoded_counts[index];
    Kokkos::deep_copy(encoded, encoded_host);

    Kokkos::View<std::uint64_t*> counts(
        "variant_eval_validation_counts",
        ValidationReportMetrics::status_count * ValidationReportMetrics::status_count);
    Kokkos::deep_copy(counts, std::uint64_t{0});
    Kokkos::parallel_for(
        "variant_eval_validation_report",
        Kokkos::RangePolicy<>(0, encoded_counts.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            Kokkos::atomic_fetch_add(&counts(encoded(index)), std::uint64_t{1});
        });
    Kokkos::fence("variant_eval_validation_report_complete");

    const auto counts_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counts);
    ValidationReportMetrics metrics;
    for (std::size_t index = 0; index < metrics.counts.size(); ++index)
        metrics.counts[index] = counts_host(index);
    return metrics;
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_eval_intervals(std::vector<fastgatk::io::IndexedInterval>& intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& left, const auto& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<fastgatk::io::IndexedInterval> merged;
    merged.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (interval.end <= interval.begin) continue;
        if (!merged.empty() && merged.back().rid == interval.rid && interval.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    intervals.swap(merged);
}

void append_eval_interval_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<fastgatk::io::IndexedInterval>& intervals,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<fastgatk::io::IndexedInterval> incoming;
    fastgatk::io::append_interval_selector(selector, header, incoming, stats);
    normalize_eval_intervals(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        intervals.insert(intervals.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<fastgatk::io::IndexedInterval> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < intervals.size() && right < incoming.size()) {
        if (intervals[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < intervals[left].rid) { ++right; continue; }
        const int begin = std::max(intervals[left].begin, incoming[right].begin);
        const int end = std::min(intervals[left].end, incoming[right].end);
        if (begin < end) {
            auto overlap = intervals[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(std::move(overlap));
        }
        if (intervals[left].end < incoming[right].end) ++left;
        else ++right;
    }
    intervals.swap(intersection);
}

bool in_intervals(const Record& record, const std::vector<fastgatk::io::IndexedInterval>& intervals) {
    if (intervals.empty()) return true;
    return std::any_of(intervals.begin(), intervals.end(), [&](const auto& interval) {
        return interval.rid == record.rid && record.end_exclusive > interval.begin &&
               record.position - 1 < interval.end;
    });
}

HeaderShape header_shape(const bcf_hdr_t* header, const std::string& path) {
    HeaderShape shape;
    if (header == nullptr) throw std::runtime_error("BAD_INPUT: null variant header: " + path);
    shape.contigs.reserve(static_cast<std::size_t>(header->n[BCF_DT_CTG]));
    shape.contig_lengths.reserve(static_cast<std::size_t>(header->n[BCF_DT_CTG]));
    for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
        const auto* name = bcf_hdr_id2name(header, rid);
        if (name == nullptr || *name == '\0')
            throw std::runtime_error("BAD_INPUT: unnamed contig in variant header: " + path);
        const auto* entry = header->id[BCF_DT_CTG][rid].val;
        const auto length = entry == nullptr ? 0 : entry->info[0];
        shape.contigs.emplace_back(name);
        shape.contig_lengths.push_back(length);
    }
    const int samples = bcf_hdr_nsamples(header);
    shape.samples.reserve(static_cast<std::size_t>(std::max(samples, 0)));
    std::set<std::string> unique_samples;
    for (int index = 0; index < samples; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
        if (name == nullptr || *name == '\0' || !unique_samples.insert(name).second)
            throw std::runtime_error("BAD_INPUT: duplicate/empty sample in variant header: " + path);
        shape.samples.emplace_back(name);
    }
    return shape;
}

void validate_header_shape(const bcf_hdr_t* header, const HeaderShape& expected,
                           const std::string& path, bool require_samples) {
    const auto actual = header_shape(header, path);
    if (actual.contigs != expected.contigs || actual.contig_lengths != expected.contig_lengths)
        throw std::runtime_error("BAD_INPUT: incompatible contig dictionary: " + path);
    if (require_samples && actual.samples != expected.samples)
        throw std::runtime_error("BAD_INPUT: incompatible sample columns: " + path);
}

std::vector<Record> load_records(const std::string& path,
                                 const std::vector<fastgatk::io::IndexedInterval>& intervals,
                                 ReadStats& stats,
                                 const HeaderShape* expected_shape = nullptr,
                                 bool require_samples = false) {
    htsFile* file = bcf_open(path.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: cannot open variant input: " + path);
    bcf_hdr_t* header = bcf_hdr_read(file);
    if (!header) {
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: cannot read variant header: " + path);
    }
    if (expected_shape != nullptr)
        validate_header_shape(header, *expected_shape, path, require_samples);
    std::vector<Record> records;
    bcf1_t* raw = bcf_init();
    int read_status = 0;
    while ((read_status = bcf_read(file, header, raw)) >= 0) {
        bcf_unpack(raw, BCF_UN_ALL);
        Record record = decode_record(header, raw);
        if (in_intervals(record, intervals)) records.push_back(std::move(record));
    }
    // HTSlib uses -1 for a clean EOF and negative values below -1 for a
    // malformed/truncated stream.  Treating every negative value as EOF
    // would let a partial VCF/BCF produce a plausible but incomplete report.
    if (read_status < -1) {
        bcf_destroy(raw);
        bcf_hdr_destroy(header);
        bcf_close(file);
        throw std::runtime_error("BAD_INPUT: failed reading variant input: " + path);
    }
    bcf_destroy(raw);
    bcf_hdr_destroy(header);
    bcf_close(file);
    stats.records += records.size();
    return records;
}

std::uint64_t evaluation_territory(const bcf_hdr_t* header,
                                   const std::vector<fastgatk::io::IndexedInterval>& intervals) {
    std::vector<fastgatk::io::IndexedInterval> covered = intervals;
    if (covered.empty()) {
        for (int rid = 0; rid < header->n[BCF_DT_CTG]; ++rid) {
            const auto* value = header->id[BCF_DT_CTG][rid].val;
            if (value != nullptr && value->info[0] > 0) {
                const auto length = std::min<std::uint64_t>(
                    value->info[0], static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
                covered.push_back({rid, 0, static_cast<int>(length), {}});
            }
        }
    } else {
        for (auto& interval : covered) {
            if (interval.end != std::numeric_limits<int>::max() ||
                interval.rid < 0 || interval.rid >= header->n[BCF_DT_CTG]) continue;
            const auto* value = header->id[BCF_DT_CTG][interval.rid].val;
            if (value != nullptr && value->info[0] > 0)
                interval.end = static_cast<int>(std::min<std::uint64_t>(
                    value->info[0], static_cast<std::uint64_t>(std::numeric_limits<int>::max())));
        }
    }
    std::sort(covered.begin(), covered.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.rid != rhs.rid ? lhs.rid < rhs.rid :
               (lhs.begin != rhs.begin ? lhs.begin < rhs.begin : lhs.end < rhs.end);
    });
    std::uint64_t total = 0;
    int current_rid = -1;
    int current_begin = 0;
    int current_end = 0;
    for (const auto& interval : covered) {
        if (interval.rid != current_rid || interval.begin > current_end) {
            if (current_rid >= 0 && current_end > current_begin)
                total += static_cast<std::uint64_t>(current_end - current_begin);
            current_rid = interval.rid;
            current_begin = interval.begin;
            current_end = interval.end;
        } else {
            current_end = std::max(current_end, interval.end);
        }
    }
    if (current_rid >= 0 && current_end > current_begin)
        total += static_cast<std::uint64_t>(current_end - current_begin);
    return total;
}

void accumulate_record(const Record& record, Metrics& metrics, bool keep_ac0 = false) {
    ++metrics.processed_loci;
    // CountVariants' nCalledLoci is a record/locus counter, not a count of
    // samples with a non-missing genotype.  It therefore includes records
    // containing only no-call genotypes and sites-only records.
    ++metrics.called_loci;
    if (record.filtered) ++metrics.filtered_sites;
    metrics.no_calls += record.genotypes.no_call;
    metrics.genotype_called_not_filtered += record.genotypes.called_not_filtered;
    metrics.genotype_no_call_or_filtered += record.genotypes.no_call_or_filtered;
    metrics.hom_ref += record.genotypes.hom_ref;
    metrics.het += record.genotypes.het;
    metrics.hom_var += record.genotypes.hom_var;
    const bool effective_variant = has_effective_variant_genotype(record, keep_ac0);
    if (!effective_variant) {
        ++metrics.reference_loci;
    } else if (record.snp) { ++metrics.variant_loci; ++metrics.snps; }
    else if (record.mnp) { ++metrics.variant_loci; ++metrics.mnps; }
    else if (record.insertion) { ++metrics.variant_loci; ++metrics.insertions; }
    else if (record.deletion) { ++metrics.variant_loci; ++metrics.deletions; }
    else if (record.complex_indel) { ++metrics.variant_loci; ++metrics.complex; }
    // CountVariants reports symbolic events in nSymbolic but, matching the
    // GATK evaluator's switch, does not add them to nVariantLoci.
    else if (record.symbolic) ++metrics.symbolic;
    else if (record.mixed) { ++metrics.variant_loci; ++metrics.mixed; }
    if (effective_variant && record.snp) {
        if (record.transition) ++metrics.ti;
        else ++metrics.tv;
    }

    // The standard IndelSummary and MultiallelicSummary evaluators count
    // sites separately from concrete alternate alleles.  Keep this update
    // alongside the existing site counters so mixed/symbolic records are not
    // accidentally promoted to indels.
    if (effective_variant && record.snp) {
        ++metrics.snp_sites;
        if (record.genotypes.alt_chromosomes == 1) ++metrics.snp_singletons;
        if (record.alts.size() > 1) {
            ++metrics.multiallelic_snp_sites;
            for (const auto& alt : record.alts) {
                if (!is_base_snp(record.ref, alt)) continue;
                if (is_transition_base(record.ref[0], alt[0])) ++metrics.multiallelic_snp_ti;
                else ++metrics.multiallelic_snp_tv;
            }
        }
    }
    if (effective_variant && (record.insertion || record.deletion || record.complex_indel)) {
        ++metrics.indel_sites;
        if (record.alts.size() > 1) ++metrics.multiallelic_indel_sites;
        for (const auto& alt : record.alts) {
            if (alt.empty() || alt[0] == '<' || alt[0] == '*') continue;
            const auto delta = static_cast<long long>(alt.size()) - static_cast<long long>(record.ref.size());
            if (delta == 0) continue;
            ++metrics.indel_alleles;
            if (record.genotypes.alt_chromosomes == 1) ++metrics.indel_singletons;
            if (delta > 0) {
                ++metrics.indel_insertions;
                if (delta > 10) ++metrics.large_insertions;
            } else {
                ++metrics.indel_deletions;
                if (-delta > 10) ++metrics.large_deletions;
            }
            if (!record.complex_indel && std::abs(delta) <= 10) {
                ++metrics.indel_length_histogram[static_cast<int>(delta)];
                ++metrics.indel_histogram_total;
            }
        }
    }
}

// Port the bounded, diploid part of GATK VariantAFEvaluator.  GATK ignores
// AC=0 SNP sites by default; records with FORMAT/GT contribute one observation
// per called sample, while sites-only records contribute their INFO/AF value.
void accumulate_variant_af(const Record& record, Metrics& metrics, bool keep_ac0 = false) {
    if (!record.snp || (!keep_ac0 && !record.genotypes.classes.empty() &&
                        !record.genotypes.polymorphic)) return;
    if (!record.genotypes.classes.empty()) {
        if (!record.genotypes.polymorphic) return;
        for (const auto genotype : record.genotypes.classes) {
            if (genotype == 0) continue; // no-call
            ++metrics.variant_af_called_sites;
            if (genotype == 1) {
                ++metrics.variant_af_hom_ref_sites;
                continue;
            }
            if (genotype == 2) {
                metrics.variant_af_sum += 0.5;
                ++metrics.variant_af_het_sites;
            } else if (genotype == 3) {
                metrics.variant_af_sum += 1.0;
                ++metrics.variant_af_hom_var_sites;
            }
        }
        return;
    }
    ++metrics.variant_af_called_sites;
    metrics.variant_af_sum += record.allele_frequency;
}

// Port the site-level portion of GATK ThetaVariantEvaluator.  The evaluator
// ignores AC=0 SNPs, averages heterozygosity over called samples, and uses all
// called allele copies for the pairwise-difference estimate.
void accumulate_theta(const Record& record, Metrics& metrics, bool keep_ac0 = false) {
    if (!record.snp || record.genotypes.classes.empty() ||
        (!keep_ac0 && !record.genotypes.polymorphic))
        return;
    const auto sample_count = record.genotypes.classes.size();
    std::uint64_t called_genotypes = 0;
    std::uint64_t het_genotypes = 0;
    std::map<int, std::uint64_t> allele_counts;
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
        const auto genotype = record.genotypes.classes[sample];
        if (genotype == 0) continue;
        ++called_genotypes;
        if (genotype == 2) ++het_genotypes;
        if (sample >= record.genotypes.alleles.size()) continue;
        for (const auto allele : record.genotypes.alleles[sample])
            if (allele >= 0) ++allele_counts[allele];
    }
    if (called_genotypes == 0) return;
    ++metrics.theta_num_sites;
    metrics.theta_total_het += static_cast<double>(het_genotypes) /
                               static_cast<double>(called_genotypes);
    double harmonic = 0.0;
    for (std::size_t index = 1; index <= sample_count; ++index)
        harmonic += 1.0 / static_cast<double>(index);
    if (harmonic > 0.0) metrics.theta_region_num_sites += 1.0 / harmonic;

    std::uint64_t total_alleles = 0;
    std::uint64_t same_pairs = 0;
    for (const auto& [allele, count] : allele_counts) {
        (void)allele;
        total_alleles += count;
        same_pairs += count * (count - 1) / 2;
    }
    const auto total_pairs = total_alleles * (total_alleles - 1) / 2;
    if (total_pairs > 0)
        metrics.theta_total_avg_diffs += static_cast<double>(total_pairs - same_pairs) /
                                         static_cast<double>(total_pairs);
}

std::string decimal(double value, int precision) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(precision) << value;
    return output.str();
}

std::vector<Trio> read_pedigree(const std::string& path) {
    std::vector<Trio> trios;
    if (path.empty()) return trios;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open pedigree: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::vector<std::string> values;
        std::string value;
        while (fields >> value) values.push_back(value);
        if (values.size() < 4 || values[0].empty() || values[1].empty()) continue;
        const auto& family = values[0];
        const auto& child = values[1];
        const auto& father = values[2];
        const auto& mother = values[3];
        if (father == "0" || father == "." || mother == "0" || mother == ".")
            continue;
        trios.push_back({family, mother, father, child});
    }
    return trios;
}

int sample_index(const Record& record, const std::string& name) {
    const auto found = std::find(record.sample_names.begin(), record.sample_names.end(), name);
    return found == record.sample_names.end()
        ? -1 : static_cast<int>(std::distance(record.sample_names.begin(), found));
}

int trio_genotype(const Record& record, int sample) {
    return sample < 0 || static_cast<std::size_t>(sample) >= record.genotypes.classes.size()
        ? 0 : record.genotypes.classes[static_cast<std::size_t>(sample)];
}

int trio_gq(const Record& record, int sample) {
    return sample < 0 || static_cast<std::size_t>(sample) >= record.genotypes.gq.size()
        ? -1 : record.genotypes.gq[static_cast<std::size_t>(sample)];
}

bool mendelian_violation(int mother, int father, int child) {
    if (child == 0) return false;
    if (mother == 0)
        return (father == 1 && child == 3) || (father == 3 && child == 1);
    if (father == 0)
        return (mother == 1 && child == 3) || (mother == 3 && child == 1);
    if (mother == 1 && father == 1) return child != 1;
    if ((mother == 1 && father == 2) || (mother == 2 && father == 1)) return child == 3;
    if ((mother == 1 && father == 3) || (mother == 3 && father == 1)) return child != 2;
    if (mother == 2 && father == 2) return false;
    if ((mother == 2 && father == 3) || (mother == 3 && father == 2)) return child == 1;
    if (mother == 3 && father == 3) return child != 3;
    return false;
}

void accumulate_mendelian(const Record& record, const std::vector<Trio>& trios,
                           double min_gq, MendelianMetrics& metrics) {
    if (record.alts.size() != 1 || record.genotypes.classes.empty()) return;
    if (trios.empty()) {
        ++metrics.n_skipped;
        return;
    }
    std::uint64_t called = 0;
    std::uint64_t low_qual = 0;
    std::uint64_t no_call = 0;
    std::uint64_t var_called = 0;
    std::uint64_t violations = 0;
    for (const auto& trio : trios) {
        const int mother_index = sample_index(record, trio.mother);
        const int father_index = sample_index(record, trio.father);
        const int child_index = sample_index(record, trio.child);
        if (mother_index < 0 || father_index < 0 || child_index < 0) continue;
        const int mother = trio_genotype(record, mother_index);
        const int father = trio_genotype(record, father_index);
        const int child = trio_genotype(record, child_index);
        if ((mother == 0 && father == 0) || child == 0) {
            ++no_call;
            continue;
        }
        if (min_gq > 0.0 &&
            (trio_gq(record, mother_index) < min_gq ||
             trio_gq(record, father_index) < min_gq ||
             trio_gq(record, child_index) < min_gq)) {
            ++low_qual;
            continue;
        }
        ++called;
        if (mother != 1 || father != 1 || child != 1) ++var_called;
        if (mendelian_violation(mother, father, child)) {
            ++violations;
            if (mother == 1 && father == 1 && child == 3) ++metrics.mv_ref_ref_var;
            if (mother == 1 && father == 1 && child == 2) ++metrics.mv_ref_ref_het;
            if (((mother == 1 && father == 2) || (mother == 2 && father == 1)) && child == 3)
                ++metrics.mv_ref_het_var;
            if (((mother == 1 && father == 3) || (mother == 3 && father == 1)) && child == 3)
                ++metrics.mv_ref_var_var;
            if (((mother == 1 && father == 3) || (mother == 3 && father == 1)) && child == 1)
                ++metrics.mv_ref_var_ref;
            if (((mother == 2 && father == 3) || (mother == 3 && father == 2)) && child == 1)
                ++metrics.mv_var_het_ref;
            if (mother == 3 && father == 3 && child == 1) ++metrics.mv_var_var_ref;
            if (mother == 3 && father == 3 && child == 2) ++metrics.mv_var_var_het;
        }
        if (mother == 1 && father == 1 && child == 1) ++metrics.hom_ref_hom_ref_hom_ref;
        if (mother == 2 && father == 2 && child == 2) ++metrics.het_het_het;
        if (mother == 2 && father == 2 && child == 1) ++metrics.het_het_hom_ref;
        if (mother == 2 && father == 2 && child == 3) ++metrics.het_het_hom_var;
        if (mother == 3 && father == 3 && child == 3) ++metrics.hom_var_hom_var_hom_var;
        if (((mother == 1 && father == 3) || (mother == 3 && father == 1)) && child == 2)
            ++metrics.hom_ref_hom_var_het;
        if (mother == 2 && father == 2) {
            if (child == 2) {
                ++metrics.het_het_inherited_ref;
                ++metrics.het_het_inherited_var;
            } else if (child == 1) {
                metrics.het_het_inherited_ref += 2;
            } else if (child == 3) {
                metrics.het_het_inherited_var += 2;
            }
        }
        if (((mother == 1 && father == 2) || (mother == 2 && father == 1))) {
            if (child == 1) ++metrics.hom_ref_het_inherited_ref;
            if (child == 2) ++metrics.hom_ref_het_inherited_var;
        }
        if (((mother == 3 && father == 2) || (mother == 2 && father == 3))) {
            if (child == 2) ++metrics.hom_var_het_inherited_ref;
            if (child == 3) ++metrics.hom_var_het_inherited_var;
        }
    }
    if (called == 0 && low_qual == 0 && no_call == 0) {
        ++metrics.n_skipped;
        return;
    }
    ++metrics.n_variants;
    metrics.n_fam_called += called;
    metrics.n_var_fam_called += var_called;
    metrics.n_low_qual += low_qual;
    metrics.n_no_call += no_call;
    if (violations > 0) {
        ++metrics.n_loci_violations;
        metrics.n_violations += violations;
    }
}

void accumulate_sample_record(const Record& record, std::size_t sample,
                              SampleCountMetrics& metrics, bool keep_ac0 = false) {
    ++metrics.called_loci;
    const int genotype = sample < record.genotypes.classes.size()
        ? record.genotypes.classes[sample] : 0;
    const bool promote_ac0 = keep_ac0 && record.variant;
    const auto count_variant_type = [&]() {
        ++metrics.variant_loci;
        if ((record.snp || record.mnp) &&
            sample < record.genotypes.alt_chromosomes_per_sample.size() &&
            record.genotypes.alt_chromosomes_per_sample[sample] == 1)
            ++metrics.singletons;
        if (record.snp) ++metrics.snps;
        else if (record.mnp) ++metrics.mnps;
        else if (record.insertion) ++metrics.insertions;
        else if (record.deletion) ++metrics.deletions;
        else if (record.complex_indel) ++metrics.complex;
        else if (record.symbolic) ++metrics.symbolic;
        else if (record.mixed) ++metrics.mixed;
    };
    if (genotype == 0) {
        ++metrics.no_calls;
        if (promote_ac0) count_variant_type();
        else ++metrics.reference_loci;
        return;
    }
    if (genotype == 1) {
        ++metrics.hom_ref;
        if (promote_ac0) count_variant_type();
        else ++metrics.reference_loci;
        return;
    }
    count_variant_type();
    if (genotype == 2) ++metrics.het;
    else if (genotype == 3) ++metrics.hom_var;
}

void accumulate_family_record(const Record& record,
                              const std::vector<std::string>& members,
                              SampleCountMetrics& metrics, bool keep_ac0 = false) {
    ++metrics.called_loci;
    std::size_t present = 0;
    std::size_t alt_carriers = 0;
    bool any_alt = false;
    for (const auto& member : members) {
        const auto index = sample_index(record, member);
        if (index < 0) continue;
        ++present;
        const int genotype = trio_genotype(record, index);
        if (genotype == 0) {
            ++metrics.no_calls;
        } else if (genotype == 1) {
            ++metrics.hom_ref;
        } else {
            any_alt = true;
            ++alt_carriers;
            if (genotype == 2) ++metrics.het;
            else if (genotype == 3) ++metrics.hom_var;
        }
    }
    if (present == 0) return;
    // Family stratification subsets genotypes before CountVariants sees the
    // record, so an all-hom-ref/all-no-call family is AC=0 as well.
    if (!any_alt && keep_ac0 && record.variant) any_alt = true;
    if (any_alt) {
        ++metrics.variant_loci;
        if ((record.snp || record.mnp) && alt_carriers == 1) ++metrics.singletons;
        if (record.snp) ++metrics.snps;
        else if (record.mnp) ++metrics.mnps;
        else if (record.insertion) ++metrics.insertions;
        else if (record.deletion) ++metrics.deletions;
        else if (record.complex_indel) ++metrics.complex;
        else if (record.symbolic) ++metrics.symbolic;
        else if (record.mixed) ++metrics.mixed;
    } else {
        ++metrics.reference_loci;
    }
}

void accumulate_stratum(const Record& record, StratumMetrics& metrics,
                        bool keep_ac0 = false) {
    ++metrics.loci;
    const bool effective_variant = has_effective_variant_genotype(record, keep_ac0);
    if (effective_variant && record.variant && !record.symbolic) ++metrics.variants;
    if (effective_variant && record.snp) ++metrics.snps;
    if (effective_variant && record.mnp) ++metrics.mnps;
    if (effective_variant && record.insertion) ++metrics.insertions;
    if (effective_variant && record.deletion) ++metrics.deletions;
    if (effective_variant && record.complex_indel) ++metrics.complex;
    if (effective_variant && record.symbolic) ++metrics.symbolic;
    if (effective_variant && record.mixed) ++metrics.mixed;
    if (effective_variant && record.snp) {
        if (record.transition) ++metrics.ti;
        else ++metrics.tv;
    }
}

std::string allele_frequency_bin(double frequency) {
    if (frequency <= 0.0) return "0";
    if (frequency < 0.01) return "0-0.01";
    if (frequency < 0.05) return "0.01-0.05";
    if (frequency < 0.10) return "0.05-0.10";
    if (frequency < 0.50) return "0.10-0.50";
    if (frequency < 1.0) return "0.50-1.00";
    return "1.00";
}

void write_stratification_table(
    std::ostream& output, const std::string& table,
    const std::map<std::string, StratumMetrics>& strata) {
    output << "##table=" << table << '\n';
    output << "Stratum\tnProcessedLoci\tnVariantLoci\tnSNPs\tnMNPs\tnInsertions\tnDeletions\tnComplex\tnSymbolic\tnMixed\tnTi\tnTv\n";
    for (const auto& [name, metrics] : strata) {
        output << name << '\t' << metrics.loci << '\t' << metrics.variants << '\t'
               << metrics.snps << '\t' << metrics.mnps << '\t' << metrics.insertions << '\t'
               << metrics.deletions << '\t' << metrics.complex << '\t' << metrics.symbolic << '\t'
               << metrics.mixed << '\t' << metrics.ti << '\t' << metrics.tv << '\n';
    }
    output << '\n';
}

bool contains_name(const std::vector<std::string>& values, const std::string& name) {
    return std::find(values.begin(), values.end(), name) != values.end();
}

std::string ratio(std::uint64_t numerator, std::uint64_t denominator) {
    if (denominator == 0) return "0.00000000";
    std::ostringstream output;
    output << std::fixed << std::setprecision(8)
           << static_cast<double>(numerator) / static_cast<double>(denominator);
    return output.str();
}

void write_metric_table(std::ostream& output, const std::string& table,
                        const std::vector<std::pair<std::string, std::string>>& rows) {
    output << "##table=" << table << '\n';
    output << "Metric\tValue\n";
    for (const auto& row : rows) output << row.first << '\t' << row.second << '\n';
    output << '\n';
}

void write_sample_count_table(std::ostream& output,
                              const std::map<std::string, SampleCountMetrics>& samples,
                              std::uint64_t territory) {
    output << "##table=CountVariantsBySample\n"
           << "Sample\tnProcessedLoci\tnCalledLoci\tnReferenceLoci\tnVariantLoci"
              "\tnSNPs\tnMNPs\tnInsertions\tnDeletions\tnComplex\tnSymbolic\tnMixed"
              "\tnNoCalls\tnHets\tnHomRef\tnHomVar\tnSingletons\n";
    for (const auto& [name, sample] : samples) {
        output << name << '\t' << territory << '\t' << sample.called_loci << '\t'
               << sample.reference_loci << '\t' << sample.variant_loci << '\t'
               << sample.snps << '\t' << sample.mnps << '\t' << sample.insertions << '\t'
               << sample.deletions << '\t' << sample.complex << '\t' << sample.symbolic << '\t'
               << sample.mixed << '\t' << sample.no_calls << '\t' << sample.het << '\t'
               << sample.hom_ref << '\t' << sample.hom_var << '\t' << sample.singletons << '\n';
    }
    output << '\n';
}

void write_family_count_table(std::ostream& output,
                              const std::map<std::string, SampleCountMetrics>& families,
                              std::uint64_t territory) {
    output << "##table=CountVariantsByFamily\n"
           << "Family\tnProcessedLoci\tnCalledLoci\tnReferenceLoci\tnVariantLoci"
              "\tnSNPs\tnMNPs\tnInsertions\tnDeletions\tnComplex\tnSymbolic\tnMixed"
              "\tnNoCalls\tnHets\tnHomRef\tnHomVar\tnSingletons\n";
    for (const auto& [name, family] : families) {
        output << name << '\t' << territory << '\t' << family.called_loci << '\t'
               << family.reference_loci << '\t' << family.variant_loci << '\t'
               << family.snps << '\t' << family.mnps << '\t' << family.insertions << '\t'
               << family.deletions << '\t' << family.complex << '\t' << family.symbolic << '\t'
               << family.mixed << '\t' << family.no_calls << '\t' << family.het << '\t'
               << family.hom_ref << '\t' << family.hom_var << '\t' << family.singletons << '\n';
    }
    output << '\n';
}

void write_indel_histogram_table(std::ostream& output,
                                 const std::map<int, std::uint64_t>& histogram,
                                 const std::uint64_t total) {
    output << "##table=IndelLengthHistogram\nLength\tFreq\n";
    output << std::fixed << std::setprecision(2);
    for (int length = -10; length <= 10; ++length) {
        if (length == 0) continue;
        const auto found = histogram.find(length);
        const auto count = found == histogram.end() ? 0ULL : found->second;
        const double frequency = total == 0 ? 0.0 : static_cast<double>(count) / static_cast<double>(total);
        output << length << '\t' << frequency << '\n';
    }
    output << '\n';
}

std::string format_gatk_cell(const std::string& value, const std::string& format) {
    if (format == "%s") return value;
    if (value == "NaN" || value == "Infinity" || value == "-Infinity") return value;
    if (format == "%d") {
        std::size_t consumed = 0;
        try {
            const auto integer = std::stoll(value, &consumed);
            if (consumed != value.size()) return value;
            return std::to_string(integer);
        } catch (const std::exception&) {
            return value;
        }
    }
    if (format.size() < 4 || format.front() != '%' ||
        (format.back() != 'f' && format.back() != 'e'))
        throw std::invalid_argument("unsupported GATKReport format: " + format);
    const auto dot = format.find('.');
    const auto precision_end = format.find(format.back(), dot == std::string::npos ? 0 : dot);
    if (dot == std::string::npos || precision_end == std::string::npos)
        throw std::invalid_argument("malformed GATKReport format: " + format);
    int precision = 0;
    try {
        precision = std::stoi(format.substr(dot + 1, precision_end - dot - 1));
    } catch (const std::exception&) {
        throw std::invalid_argument("malformed GATKReport format: " + format);
    }
    std::size_t consumed = 0;
    double number = 0.0;
    try {
        number = std::stod(value, &consumed);
        if (consumed != value.size()) return value;
    } catch (const std::exception&) {
        return value;
    }
    std::ostringstream formatted;
    if (format.back() == 'e') formatted << std::scientific;
    else formatted << std::fixed;
    formatted << std::setprecision(precision) << number;
    return formatted.str();
}

void write_gatk_table(std::ostream& output, const std::string& name,
                      const std::string& description,
                      const std::vector<std::string>& columns,
                      const std::vector<std::vector<std::string>>& rows,
                      const std::vector<std::string>& formats = {}) {
    const std::vector<std::string> default_formats(columns.size(), "%s");
    const auto& cell_formats = formats.empty() ? default_formats : formats;
    if (cell_formats.size() != columns.size())
        throw std::invalid_argument("GATKReport format/column size mismatch for " + name);
    std::vector<std::vector<std::string>> formatted_rows;
    formatted_rows.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.size() != columns.size())
            throw std::invalid_argument("GATKReport row/column size mismatch for " + name);
        std::vector<std::string> formatted;
        formatted.reserve(row.size());
        for (std::size_t index = 0; index < row.size(); ++index)
            formatted.push_back(format_gatk_cell(row[index], cell_formats[index]));
        formatted_rows.push_back(std::move(formatted));
    }

    std::vector<std::size_t> widths(columns.size(), 0);
    for (std::size_t index = 0; index < columns.size(); ++index) {
        widths[index] = columns[index].size();
        for (const auto& row : formatted_rows)
            widths[index] = std::max(widths[index], row[index].size());
    }
    output << "#:GATKTable:" << columns.size() << ':' << rows.size() << ':';
    for (const auto& format : cell_formats) output << format << ':';
    output << ";\n#:GATKTable:" << name << ':' << description << '\n';
    for (std::size_t index = 0; index < columns.size(); ++index) {
        if (index) output << "  ";
        output << std::left << std::setw(static_cast<int>(widths[index])) << columns[index];
    }
    output << '\n';
    for (const auto& row : formatted_rows) {
        for (std::size_t index = 0; index < row.size(); ++index) {
            if (index) output << "  ";
            if (cell_formats[index] == "%s") output << std::left;
            else output << std::right;
            output << std::setw(static_cast<int>(widths[index])) << row[index];
        }
        output << '\n';
    }
    output << '\n';
}

void write_mendelian_gatk_table(std::ostream& output, const std::string& comp_name,
                                const MendelianMetrics& m) {
    const std::vector<std::string> columns{
        "MendelianViolationEvaluator", "CompFeatureInput", "EvalFeatureInput",
        "nVariants", "nSkipped", "nFamCalled", "nVarFamCalled", "nLowQual", "nNoCall",
        "nLociViolations", "nViolations", "mvRefRef_Var", "mvRefRef_Het",
        "mvRefHet_Var", "mvRefVar_Var", "mvRefVar_Ref", "mvVarHet_Ref",
        "mvVarVar_Ref", "mvVarVar_Het", "HomRefHomRef_HomRef", "HetHet_Het",
        "HetHet_HomRef", "HetHet_HomVar", "HomVarHomVar_HomVar", "HomRefHomVAR_Het",
        "HetHet_inheritedRef", "HetHet_inheritedVar", "HomRefHet_inheritedRef",
        "HomRefHet_inheritedVar", "HomVarHet_inheritedRef", "HomVarHet_inheritedVar"};
    const std::vector<std::string> row{
        "MendelianViolationEvaluator", comp_name, "eval",
        std::to_string(m.n_variants), std::to_string(m.n_skipped),
        std::to_string(m.n_fam_called), std::to_string(m.n_var_fam_called),
        std::to_string(m.n_low_qual), std::to_string(m.n_no_call),
        std::to_string(m.n_loci_violations), std::to_string(m.n_violations),
        std::to_string(m.mv_ref_ref_var), std::to_string(m.mv_ref_ref_het),
        std::to_string(m.mv_ref_het_var), std::to_string(m.mv_ref_var_var),
        std::to_string(m.mv_ref_var_ref), std::to_string(m.mv_var_het_ref),
        std::to_string(m.mv_var_var_ref), std::to_string(m.mv_var_var_het),
        std::to_string(m.hom_ref_hom_ref_hom_ref), std::to_string(m.het_het_het),
        std::to_string(m.het_het_hom_ref), std::to_string(m.het_het_hom_var),
        std::to_string(m.hom_var_hom_var_hom_var), std::to_string(m.hom_ref_hom_var_het),
        std::to_string(m.het_het_inherited_ref), std::to_string(m.het_het_inherited_var),
        std::to_string(m.hom_ref_het_inherited_ref), std::to_string(m.hom_ref_het_inherited_var),
        std::to_string(m.hom_var_het_inherited_ref), std::to_string(m.hom_var_het_inherited_var)};
    write_gatk_table(output, "MendelianViolationEvaluator",
                     "Mendelian Violation Evaluator", columns, {row},
                     {"%s", "%s", "%s",
                      "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d",
                      "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d",
                      "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d",
                      "%d"});
}

struct ValidationReportSummary {
    std::uint64_t n_comp = 0;
    std::uint64_t tp = 0;
    std::uint64_t fp = 0;
    std::uint64_t fn = 0;
    std::uint64_t tn = 0;
    std::uint64_t comp_mono_eval_no_call = 0;
    std::uint64_t comp_mono_eval_filtered = 0;
    std::uint64_t comp_mono_eval_mono = 0;
    std::uint64_t comp_mono_eval_poly = 0;
    std::uint64_t comp_poly_eval_no_call = 0;
    std::uint64_t comp_poly_eval_filtered = 0;
    std::uint64_t comp_poly_eval_mono = 0;
    std::uint64_t comp_poly_eval_poly = 0;
    std::uint64_t comp_filtered = 0;
};

ValidationReportSummary summarize_validation_report(const ValidationReportMetrics& metrics) {
    ValidationReportSummary summary;
    using Status = ValidationSiteStatus;
    summary.comp_mono_eval_no_call = metrics.get(Status::Mono, Status::NoCall);
    summary.comp_mono_eval_filtered = metrics.get(Status::Mono, Status::Filtered);
    summary.comp_mono_eval_mono = metrics.get(Status::Mono, Status::Mono);
    summary.comp_mono_eval_poly = metrics.get(Status::Mono, Status::Poly);
    summary.comp_poly_eval_no_call = metrics.get(Status::Poly, Status::NoCall);
    summary.comp_poly_eval_filtered = metrics.get(Status::Poly, Status::Filtered);
    summary.comp_poly_eval_mono = metrics.get(Status::Poly, Status::Mono);
    summary.comp_poly_eval_poly = metrics.get(Status::Poly, Status::Poly);
    for (std::size_t eval = 0; eval < ValidationReportMetrics::status_count; ++eval)
        summary.comp_filtered += metrics.counts[
            static_cast<std::size_t>(Status::Filtered) * ValidationReportMetrics::status_count + eval];
    for (const auto value : metrics.counts) summary.n_comp += value;
    summary.tp = summary.comp_poly_eval_poly;
    summary.fn = summary.comp_poly_eval_no_call + summary.comp_poly_eval_filtered +
                 summary.comp_poly_eval_mono;
    summary.fp = summary.comp_mono_eval_poly;
    summary.tn = summary.comp_mono_eval_no_call + summary.comp_mono_eval_filtered +
                 summary.comp_mono_eval_mono;
    if (summary.n_comp != summary.tp + summary.fn + summary.fp + summary.tn +
                          summary.comp_filtered)
        throw std::runtime_error("INTERNAL_ERROR: ValidationReport confusion matrix is inconsistent");
    return summary;
}

std::string validation_percent(std::uint64_t numerator, std::uint64_t denominator) {
    if (denominator == 0) return "NaN";
    return std::to_string(100.0 * static_cast<double>(numerator) /
                          static_cast<double>(denominator));
}

void write_validation_gatk_table(std::ostream& output, const std::string& comp_name,
                                 const ValidationReportMetrics& metrics) {
    const auto summary = summarize_validation_report(metrics);
    const auto specificity = summary.tn + summary.fp == 0
        ? std::string("100")
        : validation_percent(summary.tn, summary.tn + summary.fp);
    write_gatk_table(
        output, "ValidationReport",
        "Assess site accuracy and sensitivity of callset against follow-up validation assay",
        {"ValidationReport", "CompFeatureInput", "EvalFeatureInput", "nComp", "TP", "FP",
         "FN", "TN", "sensitivity", "specificity", "PPV", "FDR",
         "CompMonoEvalNoCall", "CompMonoEvalFiltered", "CompMonoEvalMono",
         "CompMonoEvalPoly", "CompPolyEvalNoCall", "CompPolyEvalFiltered",
         "CompPolyEvalMono", "CompPolyEvalPoly", "CompFiltered", "nDifferentAlleleSites"},
        {{"ValidationReport", comp_name, "eval", std::to_string(summary.n_comp),
          std::to_string(summary.tp), std::to_string(summary.fp), std::to_string(summary.fn),
          std::to_string(summary.tn), validation_percent(summary.tp, summary.tp + summary.fn),
          specificity, validation_percent(summary.tp, summary.tp + summary.fp),
          validation_percent(summary.fp, summary.fp + summary.tp),
          std::to_string(summary.comp_mono_eval_no_call),
          std::to_string(summary.comp_mono_eval_filtered),
          std::to_string(summary.comp_mono_eval_mono),
          std::to_string(summary.comp_mono_eval_poly),
          std::to_string(summary.comp_poly_eval_no_call),
          std::to_string(summary.comp_poly_eval_filtered),
          std::to_string(summary.comp_poly_eval_mono),
          std::to_string(summary.comp_poly_eval_poly),
          std::to_string(summary.comp_filtered), "0"}},
        {"%s", "%s", "%s", "%d", "%d", "%d", "%d", "%d", "%.2f", "%.2f",
         "%.2f", "%.2f", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d",
         "%d", "%d"});
}

void write_gatk_report(std::ostream& output, const Options& options,
                       const Metrics& metrics,
                       const std::map<std::string, StratumMetrics>& contig_strata,
                       const std::map<std::string, StratumMetrics>& filter_strata,
                       const std::map<std::string, StratumMetrics>& type_strata,
                       const std::map<std::string, StratumMetrics>& allele_frequency_strata,
                       const std::map<std::string, StratumMetrics>& novelty_strata,
                       const std::map<std::string, StratumMetrics>& sample_strata,
                       const std::map<std::string, SampleCountMetrics>& sample_counts,
                       const std::map<std::string, SampleCountMetrics>& family_counts,
                       const std::map<std::string, Metrics>& novelty_counts,
                       std::uint64_t territory,
                       std::size_t eval_records) {
    // GATKReport uses the literal feature-source label "none" when no
    // comparison track is supplied; keeping this value exact helps reports
    // serve as a drop-in oracle instead of requiring a formatting shim.
    const std::string comp_name = options.comps.empty() ? "none" : "comp";
    const std::string eval_name = "eval";
    const auto want_module = [&](const std::string& name) {
        return !options.no_standard_modules || contains_name(options.modules, name);
    };
    const auto want_stratification = [&](const std::string& name) {
        return !options.no_standard_stratifications || contains_name(options.stratifications, name);
    };
    const auto report_territory = territory == 0 ? static_cast<std::uint64_t>(eval_records) : territory;
    const bool emit_count_variants = want_module("CountVariants");
    const bool emit_family_counts = emit_count_variants &&
        want_stratification("Family") && !family_counts.empty();
    const bool emit_novelty_counts = emit_count_variants &&
        options.no_standard_stratifications &&
        contains_name(options.stratifications, "Novelty") &&
        !novelty_counts.empty() && !emit_family_counts &&
        !(want_stratification("Sample") && !sample_counts.empty());
    const std::vector<std::string> count_columns{
        "CountVariants", "CompFeatureInput", "EvalFeatureInput", "nProcessedLoci",
        "nCalledLoci", "nRefLoci", "nVariantLoci", "variantRate", "variantRatePerBp",
        "nSNPs", "nMNPs", "nInsertions", "nDeletions", "nComplex", "nSymbolic",
        "nMixed", "nNoCalls", "nHets", "nHomRef", "nHomVar", "nSingletons",
        "nHomDerived", "heterozygosity", "heterozygosityPerBp", "hetHomRatio",
        "indelRate", "indelRatePerBp", "insertionDeletionRatio"};
    const std::vector<std::string> count_formats{
        "%s", "%s", "%s", "%d", "%d", "%d", "%d", "%.8f", "%.8f",
        "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d", "%d",
        "%d", "%d", "%d", "%.2e", "%.2f", "%.2f", "%.2e", "%.2f", "%.2f"};
    const auto inverse_ratio_decimal = [](std::uint64_t numerator,
                                          std::uint64_t denominator) {
        if (numerator == 0 || denominator == 0) return std::string("0");
        // CountVariants' reciprocal fields are produced by GATK's integer
        // RatioMetric helper before the report's floating format is applied.
        return std::to_string(denominator / numerator);
    };
    const auto het_hom_ratio = [](std::uint64_t het, std::uint64_t hom) {
        // GATK's RatioMetric reports one for a het-only row rather than a
        // divide-by-zero zero; an entirely empty row remains zero.
        if (hom == 0) return het == 0 ? std::string("0") : std::string("1");
        return std::to_string(static_cast<double>(het) / static_cast<double>(hom));
    };
    const auto gatk_titv_ratio = [](std::uint64_t ti, std::uint64_t tv) {
        if (tv == 0) return ti == 0 ? std::string("0") : std::to_string(ti);
        return std::to_string(static_cast<double>(ti) / static_cast<double>(tv));
    };
    const auto gatk_insertion_deletion_ratio = [](std::uint64_t insertions,
                                                   std::uint64_t deletions) {
        // GATK's RatioMetric divides by max(denominator, 1), so an
        // insertion-only stratum reports its insertion count rather than 0.
        if (deletions == 0) return std::to_string(insertions);
        return std::to_string(static_cast<double>(insertions) /
                              static_cast<double>(deletions));
    };
    const auto count_row_for = [&](const Metrics& value, std::uint64_t denominator,
                                   const std::string& stratum = std::string()) {
        // CountVariants' nSingletons is populated only for SNP/MNP records;
        // IndelSummary keeps its separate indel singleton counter.
        const auto value_singletons = value.snp_singletons;
        const auto value_indels = value.insertions + value.deletions + value.complex;
        auto row = std::vector<std::string>{
            "CountVariants", comp_name, eval_name, std::to_string(denominator),
            std::to_string(value.called_loci), std::to_string(value.reference_loci),
            std::to_string(value.variant_loci), ratio(value.variant_loci, denominator),
            inverse_ratio_decimal(value.variant_loci, denominator), std::to_string(value.snps),
            std::to_string(value.mnps), std::to_string(value.insertions),
            std::to_string(value.deletions), std::to_string(value.complex),
            std::to_string(value.symbolic), std::to_string(value.mixed),
            std::to_string(value.no_calls), std::to_string(value.het),
            std::to_string(value.hom_ref), std::to_string(value.hom_var),
            std::to_string(value_singletons), "0",
            ratio(value.het, denominator), inverse_ratio_decimal(value.het, denominator),
            het_hom_ratio(value.het, value.hom_var), ratio(value_indels, denominator),
            inverse_ratio_decimal(value_indels, denominator),
            gatk_insertion_deletion_ratio(value.insertions, value.deletions)};
        if (!stratum.empty()) row.insert(row.begin() + 3, stratum);
        return row;
    };
    const auto count_row = count_row_for(metrics, report_territory);
    std::size_t table_count = emit_count_variants ? 1 : 0;
    if (contains_name(options.modules, "TiTvVariantEvaluator") || !options.no_standard_modules) ++table_count;
    if (contains_name(options.modules, "VariantAFEvaluator")) ++table_count;
    if (contains_name(options.modules, "ThetaVariantEvaluator")) ++table_count;
    if (contains_name(options.modules, "MendelianViolationEvaluator")) ++table_count;
    if (contains_name(options.modules, "VariantSummary") || !options.no_standard_modules) ++table_count;
    if (contains_name(options.modules, "IndelSummary") || !options.no_standard_modules) ++table_count;
    if (contains_name(options.modules, "MultiallelicSummary") || !options.no_standard_modules) ++table_count;
    if (contains_name(options.modules, "GenotypeFilterSummary")) ++table_count;
    if (contains_name(options.modules, "PrintMissingComp")) ++table_count;
    if (contains_name(options.modules, "ValidationReport")) ++table_count;
    if (want_module("Contig")) ++table_count;
    if (want_module("Filter")) ++table_count;
    if (want_module("VariantType")) ++table_count;
    if (want_module("AlleleFrequency")) ++table_count;
    if (want_stratification("Novelty") && !novelty_strata.empty() && !emit_novelty_counts)
        ++table_count;
    if (want_module("Sample")) ++table_count;
    if (emit_family_counts && want_stratification("Sample") && !sample_counts.empty())
        ++table_count;
    output << "#:GATKReport.v1.1:" << table_count << '\n';
    if (emit_novelty_counts) {
        auto novelty_columns = count_columns;
        auto novelty_formats = count_formats;
        novelty_columns.insert(novelty_columns.begin() + 3, "Novelty");
        novelty_formats.insert(novelty_formats.begin() + 3, "%s");
        std::vector<std::vector<std::string>> rows;
        for (const auto& [name, value] : novelty_counts) {
            auto row = count_row_for(value, report_territory, name);
            rows.push_back(std::move(row));
        }
        write_gatk_table(output, "CountVariants", "Counts different classes of variants in the sample",
                         novelty_columns, rows, novelty_formats);
    } else if (emit_count_variants && want_stratification("Sample") && !sample_counts.empty()) {
        auto sample_columns = count_columns;
        sample_columns.insert(sample_columns.begin() + 3, "Sample");
        std::vector<std::vector<std::string>> rows;
        rows.reserve(sample_counts.size() + 1);
        for (const auto& [name, sample] : sample_counts) {
            Metrics sample_value;
            sample_value.called_loci = sample.called_loci;
            sample_value.reference_loci = sample.reference_loci;
            sample_value.variant_loci = sample.variant_loci;
            sample_value.snps = sample.snps;
            sample_value.mnps = sample.mnps;
            sample_value.insertions = sample.insertions;
            sample_value.deletions = sample.deletions;
            sample_value.complex = sample.complex;
            sample_value.symbolic = sample.symbolic;
            sample_value.mixed = sample.mixed;
            sample_value.no_calls = sample.no_calls;
            sample_value.het = sample.het;
            sample_value.hom_ref = sample.hom_ref;
            sample_value.hom_var = sample.hom_var;
            sample_value.snp_singletons = sample.singletons;
            rows.push_back(count_row_for(sample_value, report_territory, name));
        }
        // Keep the aggregate row available to downstream consumers by using
        // the synthetic ALL_SAMPLE_NAME state used by GATK's Sample stratum.
        auto aggregate = count_row;
        aggregate.insert(aggregate.begin() + 3, "all");
        rows.push_back(std::move(aggregate));
        write_gatk_table(output, "CountVariants", "Counts different classes of variants in the sample",
                         sample_columns, rows,
                         [&]() { auto f = count_formats; f.insert(f.begin() + 3, "%s"); return f; }());
    } else if (emit_count_variants && !emit_family_counts) {
        write_gatk_table(output, "CountVariants", "Counts different classes of variants in the sample",
                         count_columns, std::vector<std::vector<std::string>>{count_row}, count_formats);
    }
    if (emit_family_counts) {
        auto family_columns = count_columns;
        family_columns.insert(family_columns.begin() + 3, "Family");
        std::vector<std::vector<std::string>> rows;
        for (const auto& [name, family] : family_counts) {
            Metrics family_value;
            family_value.called_loci = family.called_loci;
            family_value.reference_loci = family.reference_loci;
            family_value.variant_loci = family.variant_loci;
            family_value.snps = family.snps;
            family_value.mnps = family.mnps;
            family_value.insertions = family.insertions;
            family_value.deletions = family.deletions;
            family_value.complex = family.complex;
            family_value.symbolic = family.symbolic;
            family_value.mixed = family.mixed;
            family_value.no_calls = family.no_calls;
            family_value.het = family.het;
            family_value.hom_ref = family.hom_ref;
            family_value.hom_var = family.hom_var;
            family_value.snp_singletons = family.singletons;
            rows.push_back(count_row_for(family_value, report_territory, name));
        }
        auto family_formats = count_formats;
        family_formats.insert(family_formats.begin() + 3, "%s");
        write_gatk_table(output, "CountVariants", "Counts different classes of variants in the sample",
                         family_columns, rows, family_formats);
    }
    if (contains_name(options.modules, "VariantAFEvaluator")) {
        const auto average = metrics.variant_af_called_sites == 0
            ? 0.0 : metrics.variant_af_sum /
                static_cast<double>(metrics.variant_af_called_sites);
        write_gatk_table(output, "VariantAFEvaluator",
                         "Computes different estimates of theta based on variant sites and genotypes",
                         {"VariantAFEvaluator", "CompFeatureInput", "EvalFeatureInput",
                          "avgVarAF", "totalCalledSites", "totalHetSites",
                          "totalHomVarSites", "totalHomRefSites"},
                         {{"VariantAFEvaluator", comp_name, eval_name,
                           decimal(average, 8),
                           std::to_string(metrics.variant_af_called_sites),
                           std::to_string(metrics.variant_af_het_sites),
                           std::to_string(metrics.variant_af_hom_var_sites),
                           std::to_string(metrics.variant_af_hom_ref_sites)}},
                         {"%s", "%s", "%s", "%.8f", "%d", "%d", "%d", "%d"});
    }
    if (contains_name(options.modules, "ThetaVariantEvaluator")) {
        const auto average_het = metrics.theta_num_sites == 0
            ? 0.0 : metrics.theta_total_het /
                static_cast<double>(metrics.theta_num_sites);
        const auto average_diffs = metrics.theta_num_sites == 0
            ? 0.0 : metrics.theta_total_avg_diffs /
                static_cast<double>(metrics.theta_num_sites);
        write_gatk_table(output, "ThetaVariantEvaluator",
                         "Computes different estimates of theta based on variant sites and genotypes",
                         {"ThetaVariantEvaluator", "CompFeatureInput", "EvalFeatureInput",
                          "avgHet", "avgAvgDiffs", "totalHet", "totalAvgDiffs",
                          "thetaRegionNumSites"},
                         {{"ThetaVariantEvaluator", comp_name, eval_name,
                           decimal(average_het, 8), decimal(average_diffs, 8),
                           decimal(metrics.theta_total_het, 8),
                           decimal(metrics.theta_total_avg_diffs, 8),
                           decimal(metrics.theta_region_num_sites, 8)}},
                         {"%s", "%s", "%s", "%.8f", "%.8f", "%.8f", "%.8f", "%.8f"});
    }
    if (contains_name(options.modules, "MendelianViolationEvaluator"))
        write_mendelian_gatk_table(output, comp_name, metrics.mendelian);
    if (contains_name(options.modules, "TiTvVariantEvaluator") || !options.no_standard_modules) {
        write_gatk_table(output, "TiTvVariantEvaluator", "Ti/Tv Variant Evaluator",
                         {"TiTvVariantEvaluator", "CompFeatureInput", "EvalFeatureInput",
                          "nTi", "nTv", "tiTvRatio", "nTiInComp", "nTvInComp",
                          "TiTvRatioStandard", "nTiDerived", "nTvDerived", "tiTvDerivedRatio"},
                         {{"TiTvVariantEvaluator", comp_name, eval_name,
                           std::to_string(metrics.ti), std::to_string(metrics.tv),
                           gatk_titv_ratio(metrics.ti, metrics.tv), std::to_string(metrics.comp_ti),
                           std::to_string(metrics.comp_tv), gatk_titv_ratio(metrics.comp_ti, metrics.comp_tv),
                           "0", "0", "0"}},
                         {"%s", "%s", "%s", "%d", "%d", "%.2f", "%d", "%d",
                          "%.2f", "%d", "%d", "%.2f"});
    }
    if (contains_name(options.modules, "VariantSummary") || !options.no_standard_modules) {
        write_gatk_table(output, "VariantSummary", "Summary of variant records",
                         {"VariantSummary", "CompFeatureInput", "EvalFeatureInput",
                          "nProcessedLoci", "nVariantLoci", "nFilteredLoci", "nNoCalls",
                          "nCalledNotFiltered", "nNoCallOrFiltered"},
                         {{"VariantSummary", comp_name, eval_name,
                           std::to_string(metrics.processed_loci), std::to_string(metrics.variant_loci),
                           std::to_string(metrics.filtered_sites), std::to_string(metrics.no_calls),
                           std::to_string(metrics.genotype_called_not_filtered),
                           std::to_string(metrics.genotype_no_call_or_filtered)}});
    }
    if (contains_name(options.modules, "IndelSummary") || !options.no_standard_modules) {
        write_gatk_table(output, "IndelSummary", "Indel summary",
                         {"IndelSummary", "CompFeatureInput", "EvalFeatureInput",
                          "n_indels", "n_indel_sites", "n_insertions", "n_deletions",
                          "n_large_insertions", "n_large_deletions"},
                         {{"IndelSummary", comp_name, eval_name,
                           std::to_string(metrics.indel_alleles), std::to_string(metrics.indel_sites),
                           std::to_string(metrics.indel_insertions), std::to_string(metrics.indel_deletions),
                           std::to_string(metrics.large_insertions), std::to_string(metrics.large_deletions)}});
    }
    if (contains_name(options.modules, "MultiallelicSummary") || !options.no_standard_modules) {
        write_gatk_table(output, "MultiallelicSummary", "Evaluation summary for multi-allelic variants",
                         {"MultiallelicSummary", "CompFeatureInput", "EvalFeatureInput",
                          "nProcessedLoci", "nSNPs", "nMultiSNPs", "nIndels", "nMultiIndels",
                          "knownSNPsPartial", "knownSNPsComplete"},
                         {{"MultiallelicSummary", comp_name, eval_name,
                           std::to_string(metrics.processed_loci), std::to_string(metrics.snp_sites),
                           std::to_string(metrics.multiallelic_snp_sites), std::to_string(metrics.indel_sites),
                           std::to_string(metrics.multiallelic_indel_sites),
                           std::to_string(metrics.known_multiallelic_snp_partial),
                           std::to_string(metrics.known_multiallelic_snp_complete)}});
    }
    if (contains_name(options.modules, "GenotypeFilterSummary")) {
        write_gatk_table(output, "GenotypeFilterSummary", "Genotype filter summary",
                         {"GenotypeFilterSummary", "CompFeatureInput", "EvalFeatureInput",
                          "nCalledNotFiltered", "nNoCallOrFiltered"},
                         {{"GenotypeFilterSummary", comp_name, eval_name,
                           std::to_string(metrics.genotype_called_not_filtered),
                           std::to_string(metrics.genotype_no_call_or_filtered)}});
    }
    if (contains_name(options.modules, "PrintMissingComp")) {
        write_gatk_table(output, "PrintMissingComp", "Count of comparison SNPs missing from eval",
                         {"PrintMissingComp", "CompFeatureInput", "EvalFeatureInput", "nMissing"},
                         {{"PrintMissingComp", comp_name, eval_name,
                           std::to_string(metrics.missing_comp_snps)}});
    }
    if (contains_name(options.modules, "ValidationReport"))
        write_validation_gatk_table(output, comp_name, metrics.validation_report);
    const auto write_strata = [&](const std::string& name,
                                  const std::map<std::string, StratumMetrics>& strata) {
        std::vector<std::vector<std::string>> rows;
        for (const auto& [stratum, value] : strata)
            rows.push_back({name, comp_name, stratum, eval_name, std::to_string(value.loci),
                            std::to_string(value.variants), std::to_string(value.snps),
                            std::to_string(value.mnps), std::to_string(value.insertions),
                            std::to_string(value.deletions), std::to_string(value.complex),
                            std::to_string(value.symbolic), std::to_string(value.mixed),
                            std::to_string(value.ti), std::to_string(value.tv)});
        write_gatk_table(output, name, "VariantEval stratification " + name,
                         {name, "CompFeatureInput", name, "EvalFeatureInput", "nProcessedLoci",
                          "nVariantLoci", "nSNPs", "nMNPs", "nInsertions", "nDeletions",
                          "nComplex", "nSymbolic", "nMixed", "nTi", "nTv"}, rows);
    };
    if (want_module("Contig")) write_strata("Contig", contig_strata);
    if (want_module("Filter")) write_strata("Filter", filter_strata);
    if (want_module("VariantType")) write_strata("VariantType", type_strata);
    if (want_module("AlleleFrequency")) write_strata("AlleleFrequency", allele_frequency_strata);
    if (want_module("Sample")) write_strata("Sample", sample_strata);
    if (want_stratification("Novelty") && !novelty_strata.empty() && !emit_novelty_counts)
        write_strata("Novelty", novelty_strata);
    (void)eval_records;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot&) {
    if (options.list) {
        std::cout << "Available stratification modules:\n"
                     "  Contig\n  Filter\n  VariantType\n  AlleleFrequency\n  Novelty\n  Sample\n  Family\n\n"
                     "Available evaluation modules:\n"
                     "  CountVariants*\n  TiTvVariantEvaluator*\n  VariantAFEvaluator\n  ThetaVariantEvaluator\n  MendelianViolationEvaluator\n  VariantSummary*\n  IndelSummary*\n  MultiallelicSummary*\n  IndelLengthHistogram*\n  GenotypeFilterSummary\n  PrintMissingComp\n  ValidationReport\n  CompOverlap*\n  GenotypeConcordance\n";
        return 0;
    }

    // The eval header is the authoritative contig dictionary for selectors.
    htsFile* header_file = bcf_open(options.evals.front().c_str(), "r");
    if (!header_file) throw std::runtime_error("BAD_INPUT: cannot open evaluation input");
    bcf_hdr_t* eval_header = bcf_hdr_read(header_file);
    if (!eval_header) {
        bcf_close(header_file);
        throw std::runtime_error("BAD_INPUT: cannot read evaluation header");
    }
    const auto eval_shape = header_shape(eval_header, options.evals.front());
    std::vector<fastgatk::io::IndexedInterval> intervals;
    fastgatk::io::IntervalFileStats interval_stats;
    bool first_interval_selector = true;
    for (const auto& selector : options.regions) {
        append_eval_interval_selector_with_rule(
            selector, eval_header, intervals, interval_stats,
            options.interval_set_rule, first_interval_selector);
        first_interval_selector = false;
    }
    normalize_eval_intervals(intervals);
    const auto territory = evaluation_territory(eval_header, intervals);
    bcf_hdr_destroy(eval_header);
    bcf_close(header_file);

    ReadStats eval_stats;
    std::vector<Record> eval_records;
    for (const auto& input : options.evals) {
        auto records = load_records(input, intervals, eval_stats, &eval_shape, true);
        eval_records.insert(eval_records.end(), records.begin(), records.end());
    }
    ReadStats comp_stats;
    std::vector<Record> comp_records;
    const bool require_comparison_samples = !options.comps.empty() &&
        (!options.no_standard_modules || contains_name(options.modules, "GenotypeConcordance"));
    for (const auto& input : options.comps) {
        auto records = load_records(input, intervals, comp_stats, &eval_shape,
                                    require_comparison_samples);
        comp_records.insert(comp_records.end(), records.begin(), records.end());
    }
    const auto trios = read_pedigree(options.pedigree);

    Metrics metrics;
    for (const auto& record : eval_records) {
        accumulate_record(record, metrics, options.keep_ac0);
        accumulate_variant_af(record, metrics, options.keep_ac0);
        accumulate_theta(record, metrics, options.keep_ac0);
        accumulate_mendelian(record, trios, options.mendelian_qual_threshold, metrics.mendelian);
    }
    metrics.comp_records = comp_records.size();
    for (const auto& record : comp_records) {
        if (!record.snp || !record.genotypes.polymorphic) continue;
        if (record.transition) ++metrics.comp_ti;
        else ++metrics.comp_tv;
    }
    if (contains_name(options.modules, "ValidationReport")) {
        if (options.evals.size() != 1 || options.comps.size() > 1)
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: bounded ValidationReport supports one eval and at most one comp track");
        if (!options.no_standard_stratifications || !options.stratifications.empty())
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: bounded ValidationReport requires -no-st without stratification modules");
        metrics.validation_report = compute_validation_report(
            eval_records, comp_records, eval_shape.samples);
    }

    std::map<std::string, StratumMetrics> contig_strata;
    std::map<std::string, StratumMetrics> filter_strata;
    std::map<std::string, StratumMetrics> type_strata;
    std::map<std::string, StratumMetrics> allele_frequency_strata;
    std::map<std::string, StratumMetrics> novelty_strata;
    std::map<std::string, StratumMetrics> sample_strata;
    std::map<std::string, StratumMetrics> family_strata;
    std::map<std::string, SampleCountMetrics> sample_counts;
    std::map<std::string, std::vector<std::string>> family_members;
    for (const auto& trio : trios) {
        auto& members = family_members[trio.family];
        for (const auto& name : {trio.mother, trio.father, trio.child})
            if (std::find(members.begin(), members.end(), name) == members.end())
                members.push_back(name);
    }
    if (contains_name(options.stratifications, "Family") && !family_members.empty()) {
        std::vector<std::string> all_members;
        for (const auto& [family, members] : family_members) {
            (void)family;
            all_members.insert(all_members.end(), members.begin(), members.end());
        }
        std::sort(all_members.begin(), all_members.end());
        all_members.erase(std::unique(all_members.begin(), all_members.end()), all_members.end());
        family_members["all"] = std::move(all_members);
    }
    std::map<std::string, SampleCountMetrics> family_counts;
    std::map<std::string, Metrics> novelty_counts;
    if (contains_name(options.stratifications, "Novelty")) {
        for (const auto& name : {std::string("all"), std::string("known"), std::string("novel")}) {
            novelty_strata.emplace(name, StratumMetrics{});
            novelty_counts.emplace(name, Metrics{});
        }
    }
    for (const auto& record : eval_records) {
        accumulate_stratum(record, contig_strata[record.contig], options.keep_ac0);
        accumulate_stratum(record, filter_strata[record.filter], options.keep_ac0);
        const auto effective_type = has_effective_variant_genotype(record, options.keep_ac0)
            ? record.variant_type : std::string("NO_VARIATION");
        accumulate_stratum(record, type_strata[effective_type], options.keep_ac0);
        accumulate_stratum(record, allele_frequency_strata[allele_frequency_bin(record.allele_frequency)],
                           options.keep_ac0);
        const bool known = record_matches_comparison(record, comp_records, options.strict_allele_match);
        accumulate_stratum(record, novelty_strata["all"], options.keep_ac0);
        accumulate_stratum(record, novelty_strata[known ? "known" : "novel"], options.keep_ac0);
        accumulate_record(record, novelty_counts["all"], options.keep_ac0);
        accumulate_record(record, novelty_counts[known ? "known" : "novel"], options.keep_ac0);
        accumulate_stratum(record, sample_strata["ALL_SAMPLES"], options.keep_ac0);
        for (std::size_t sample = 0; sample < record.sample_names.size(); ++sample) {
            const auto& name = record.sample_names[sample];
            if (name.empty()) continue;
            accumulate_stratum(record, sample_strata[name], options.keep_ac0);
            accumulate_sample_record(record, sample, sample_counts[name], options.keep_ac0);
        }
        for (const auto& [family, members] : family_members) {
            accumulate_stratum(record, family_strata[family], options.keep_ac0);
            accumulate_family_record(record, members, family_counts[family], options.keep_ac0);
        }
    }

    std::unordered_set<std::string> comp_keys;
    comp_keys.reserve(comp_records.size() * 2 + 1);
    for (const auto& record : comp_records) comp_keys.insert(record_key(record));
    for (const auto& eval : eval_records) {
        const auto exact = comp_keys.find(record_key(eval)) != comp_keys.end();
        bool overlap = exact;
        if (!overlap && !options.strict_allele_match) {
            overlap = std::any_of(comp_records.begin(), comp_records.end(), [&](const auto& comp) {
                return is_variant_key_equal(eval, comp, false);
            });
        }
        if (overlap) ++metrics.eval_comp_overlap;
        if (!exact && overlap) ++metrics.strict_allele_mismatch;
        if (eval.snp && eval.alts.size() > 1) {
            const auto comp = std::find_if(comp_records.begin(), comp_records.end(), [&](const auto& candidate) {
                return is_variant_key_equal(eval, candidate, false);
            });
            if (comp != comp_records.end()) {
                std::size_t known = 0;
                for (const auto& alt : eval.alts)
                    if (std::find(comp->alts.begin(), comp->alts.end(), alt) != comp->alts.end()) ++known;
                if (known == eval.alts.size()) ++metrics.known_multiallelic_snp_complete;
                else if (known > 0) ++metrics.known_multiallelic_snp_partial;
            }
        }
    }

    // PrintMissingComp is comparison-driven: count PASS comparison SNPs for
    // which no SNP evaluation record exists at the same locus (or with the
    // same alleles when strict matching is requested).
    if (contains_name(options.modules, "PrintMissingComp")) {
        for (const auto& comp : comp_records) {
            if (comp.filtered || !comp.snp) continue;
            const auto found = std::any_of(eval_records.begin(), eval_records.end(),
                [&](const auto& eval) {
                    return eval.snp && is_variant_key_equal(eval, comp, options.strict_allele_match);
                });
            if (!found) ++metrics.missing_comp_snps;
        }
    }

    const auto want_module = [&](const std::string& name) {
        return !options.no_standard_modules || contains_name(options.modules, name);
    };
    const bool need_genotype_concordance = !comp_records.empty() && want_module("GenotypeConcordance");
    if (need_genotype_concordance) {
        for (const auto& eval : eval_records) {
            const auto found = std::find_if(comp_records.begin(), comp_records.end(), [&](const auto& comp) {
                return is_variant_key_equal(eval, comp, options.strict_allele_match);
            });
            if (found == comp_records.end()) continue;
            const auto count = std::min(eval.genotypes.classes.size(), found->genotypes.classes.size());
            for (std::size_t sample = 0; sample < count; ++sample) {
                ++metrics.compared_eval_genotypes;
                const auto eval_class = eval.genotypes.classes[sample];
                const auto comp_class = found->genotypes.classes[sample];
                const auto class_name = [](int value) {
                    switch (value) {
                        case 1: return std::string("HOM_REF");
                        case 2: return std::string("HET");
                        case 3: return std::string("HOM_VAR");
                        default: return std::string("NO_CALL");
                    }
                };
                if (eval_class == comp_class) {
                    ++metrics.concordant_genotypes;
                    ++metrics.genotype_concordance["CONCORDANT"];
                } else {
                    ++metrics.genotype_concordance["DISCORDANT"];
                }
                ++metrics.genotype_concordance[class_name(eval_class) + "_EVAL"];
                ++metrics.genotype_concordance[class_name(comp_class) + "_COMP"];
            }
        }
    }

    const auto unsupported_module = std::any_of(options.modules.begin(), options.modules.end(), [](const auto& name) {
        return name != "CountVariants" && name != "TiTvVariantEvaluator" &&
               name != "VariantAFEvaluator" &&
               name != "ThetaVariantEvaluator" &&
               name != "MendelianViolationEvaluator" &&
               name != "VariantSummary" &&
               name != "PrintMissingComp" &&
               name != "IndelSummary" && name != "MultiallelicSummary" &&
               name != "IndelLengthHistogram" &&
               name != "GenotypeFilterSummary" &&
               name != "ValidationReport" &&
               name != "CompOverlap" && name != "GenotypeConcordance";
    });
    const auto unsupported_stratification = std::any_of(
        options.stratifications.begin(), options.stratifications.end(), [](const auto& name) {
            return name != "Contig" && name != "Filter" && name != "VariantType" &&
                   name != "AlleleFrequency" && name != "Novelty" &&
                   name != "Sample" && name != "Family";
        });
    const bool validation_report_active = contains_name(options.modules, "ValidationReport");
    const auto validation_summary = summarize_validation_report(metrics.validation_report);

    std::ofstream output(options.output);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open report: " + options.output);
    const bool need_genotype_filter_summary = contains_name(options.modules, "GenotypeFilterSummary");
    const bool sample_stratification_requested =
        !options.no_standard_stratifications || contains_name(options.stratifications, "Sample");
    const bool family_stratification_active =
        (!options.no_standard_stratifications || contains_name(options.stratifications, "Family")) &&
        !family_strata.empty();
    const auto stratification_table_count = [&]() {
        if (!options.no_standard_stratifications)
            return static_cast<std::size_t>(5 + (family_stratification_active ? 1 : 0));
        std::size_t count = 0;
        if (contains_name(options.stratifications, "Contig") && !contig_strata.empty()) ++count;
        if (contains_name(options.stratifications, "Filter") && !filter_strata.empty()) ++count;
        if (contains_name(options.stratifications, "VariantType") && !type_strata.empty()) ++count;
        if (contains_name(options.stratifications, "AlleleFrequency") && !allele_frequency_strata.empty()) ++count;
        if (contains_name(options.stratifications, "Novelty") && !novelty_counts.empty()) ++count;
        if (contains_name(options.stratifications, "Sample") && !sample_strata.empty()) ++count;
        if (contains_name(options.stratifications, "Family") && !family_strata.empty()) ++count;
        return count;
    };
    const auto novelty_record_count = [&](const std::string& name) {
        const auto found = novelty_counts.find(name);
        return found == novelty_counts.end() ? std::size_t{0} : found->second.processed_loci;
    };
    if (options.gatk_report) {
        write_gatk_report(output, options, metrics, contig_strata, filter_strata, type_strata,
                          allele_frequency_strata, novelty_strata, sample_strata, sample_counts,
                          family_counts, novelty_counts,
                          territory == 0 ? eval_records.size() : territory, eval_records.size());
    } else {
    output << "# fastgatk VariantEval native prototype\n"
           << "# status=prototype bit_identical_to_gatk=false\n";
    if (want_module("CountVariants")) write_metric_table(output, "CountVariants", {
        {"nProcessedLoci", std::to_string(metrics.processed_loci)},
        {"nCalledLoci", std::to_string(metrics.called_loci)},
        {"nReferenceLoci", std::to_string(metrics.reference_loci)},
        {"nVariantLoci", std::to_string(metrics.variant_loci)},
        {"variantRate", ratio(metrics.variant_loci, metrics.processed_loci)},
        {"nSNPs", std::to_string(metrics.snps)},
        {"nMNPs", std::to_string(metrics.mnps)},
        {"nInsertions", std::to_string(metrics.insertions)},
        {"nDeletions", std::to_string(metrics.deletions)},
        {"nComplex", std::to_string(metrics.complex)},
        {"nSymbolic", std::to_string(metrics.symbolic)},
        {"nMixed", std::to_string(metrics.mixed)},
        {"nNoCalls", std::to_string(metrics.no_calls)},
        {"nHomRef", std::to_string(metrics.hom_ref)},
        {"nHets", std::to_string(metrics.het)},
        {"nHomVar", std::to_string(metrics.hom_var)},
        {"heterozygosity", ratio(metrics.het, metrics.processed_loci)},
    });
    if (contains_name(options.modules, "VariantAFEvaluator")) {
        const auto average = metrics.variant_af_called_sites == 0
            ? 0.0 : metrics.variant_af_sum /
                static_cast<double>(metrics.variant_af_called_sites);
        write_metric_table(output, "VariantAFEvaluator", {
            {"avgVarAF", decimal(average, 8)},
            {"totalCalledSites", std::to_string(metrics.variant_af_called_sites)},
            {"totalHetSites", std::to_string(metrics.variant_af_het_sites)},
            {"totalHomVarSites", std::to_string(metrics.variant_af_hom_var_sites)},
            {"totalHomRefSites", std::to_string(metrics.variant_af_hom_ref_sites)},
        });
    }
    if (contains_name(options.modules, "ThetaVariantEvaluator")) {
        const auto average_het = metrics.theta_num_sites == 0
            ? 0.0 : metrics.theta_total_het /
                static_cast<double>(metrics.theta_num_sites);
        const auto average_diffs = metrics.theta_num_sites == 0
            ? 0.0 : metrics.theta_total_avg_diffs /
                static_cast<double>(metrics.theta_num_sites);
        write_metric_table(output, "ThetaVariantEvaluator", {
            {"avgHet", decimal(average_het, 8)},
            {"avgAvgDiffs", decimal(average_diffs, 8)},
            {"totalHet", decimal(metrics.theta_total_het, 8)},
            {"totalAvgDiffs", decimal(metrics.theta_total_avg_diffs, 8)},
            {"thetaRegionNumSites", decimal(metrics.theta_region_num_sites, 8)},
        });
    }
    if (contains_name(options.modules, "MendelianViolationEvaluator")) {
        const auto& m = metrics.mendelian;
        write_metric_table(output, "MendelianViolationEvaluator", {
            {"nVariants", std::to_string(m.n_variants)},
            {"nSkipped", std::to_string(m.n_skipped)},
            {"nFamCalled", std::to_string(m.n_fam_called)},
            {"nVarFamCalled", std::to_string(m.n_var_fam_called)},
            {"nLowQual", std::to_string(m.n_low_qual)},
            {"nNoCall", std::to_string(m.n_no_call)},
            {"nLociViolations", std::to_string(m.n_loci_violations)},
            {"nViolations", std::to_string(m.n_violations)},
            {"mvRefRef_Var", std::to_string(m.mv_ref_ref_var)},
            {"mvRefRef_Het", std::to_string(m.mv_ref_ref_het)},
            {"mvRefHet_Var", std::to_string(m.mv_ref_het_var)},
            {"mvRefVar_Var", std::to_string(m.mv_ref_var_var)},
            {"mvRefVar_Ref", std::to_string(m.mv_ref_var_ref)},
            {"mvVarHet_Ref", std::to_string(m.mv_var_het_ref)},
            {"mvVarVar_Ref", std::to_string(m.mv_var_var_ref)},
            {"mvVarVar_Het", std::to_string(m.mv_var_var_het)},
            {"HomRefHomRef_HomRef", std::to_string(m.hom_ref_hom_ref_hom_ref)},
            {"HetHet_Het", std::to_string(m.het_het_het)},
            {"HetHet_HomRef", std::to_string(m.het_het_hom_ref)},
            {"HetHet_HomVar", std::to_string(m.het_het_hom_var)},
            {"HomVarHomVar_HomVar", std::to_string(m.hom_var_hom_var_hom_var)},
            {"HomRefHomVAR_Het", std::to_string(m.hom_ref_hom_var_het)},
            {"HetHet_inheritedRef", std::to_string(m.het_het_inherited_ref)},
            {"HetHet_inheritedVar", std::to_string(m.het_het_inherited_var)},
            {"HomRefHet_inheritedRef", std::to_string(m.hom_ref_het_inherited_ref)},
            {"HomRefHet_inheritedVar", std::to_string(m.hom_ref_het_inherited_var)},
            {"HomVarHet_inheritedRef", std::to_string(m.hom_var_het_inherited_ref)},
            {"HomVarHet_inheritedVar", std::to_string(m.hom_var_het_inherited_var)},
        });
    }
    if (sample_stratification_requested && !sample_counts.empty())
        write_sample_count_table(output, sample_counts,
                                 territory == 0 ? eval_records.size() : territory);
    if (family_stratification_active && !family_counts.empty())
        write_family_count_table(output, family_counts,
                                 territory == 0 ? eval_records.size() : territory);
    if (want_module("TiTvVariantEvaluator")) write_metric_table(output, "TiTvVariantEvaluator", {
        {"nTi", std::to_string(metrics.ti)},
        {"nTv", std::to_string(metrics.tv)},
        {"tiTvRatio", ratio(metrics.ti, metrics.tv)},
    });
    if (want_module("VariantSummary")) write_metric_table(output, "VariantSummary", {
        {"nProcessedLoci", std::to_string(metrics.processed_loci)},
        {"nCalledLoci", std::to_string(metrics.called_loci)},
        {"nReferenceLoci", std::to_string(metrics.reference_loci)},
        {"nVariantLoci", std::to_string(metrics.variant_loci)},
        {"nFilteredLoci", std::to_string(metrics.filtered_sites)},
        {"nNoCallGenotypes", std::to_string(metrics.no_calls)},
        {"nCalledNotFilteredGenotypes", std::to_string(metrics.genotype_called_not_filtered)},
        {"nNoCallOrFilteredGenotypes", std::to_string(metrics.genotype_no_call_or_filtered)},
        {"nSNPs", std::to_string(metrics.snps)},
        {"nMNPs", std::to_string(metrics.mnps)},
        {"nInsertions", std::to_string(metrics.insertions)},
        {"nDeletions", std::to_string(metrics.deletions)},
        {"nComplex", std::to_string(metrics.complex)},
        {"nSymbolic", std::to_string(metrics.symbolic)},
        {"nMixed", std::to_string(metrics.mixed)},
    });
    if (contains_name(options.modules, "PrintMissingComp")) write_metric_table(
        output, "PrintMissingComp", {
            {"nMissing", std::to_string(metrics.missing_comp_snps)},
        });
    if (contains_name(options.modules, "ValidationReport")) {
        const auto validation = summarize_validation_report(metrics.validation_report);
        write_metric_table(output, "ValidationReport", {
            {"nComp", std::to_string(validation.n_comp)},
            {"TP", std::to_string(validation.tp)},
            {"FP", std::to_string(validation.fp)},
            {"FN", std::to_string(validation.fn)},
            {"TN", std::to_string(validation.tn)},
            {"CompMonoEvalNoCall", std::to_string(validation.comp_mono_eval_no_call)},
            {"CompMonoEvalFiltered", std::to_string(validation.comp_mono_eval_filtered)},
            {"CompMonoEvalMono", std::to_string(validation.comp_mono_eval_mono)},
            {"CompMonoEvalPoly", std::to_string(validation.comp_mono_eval_poly)},
            {"CompPolyEvalNoCall", std::to_string(validation.comp_poly_eval_no_call)},
            {"CompPolyEvalFiltered", std::to_string(validation.comp_poly_eval_filtered)},
            {"CompPolyEvalMono", std::to_string(validation.comp_poly_eval_mono)},
            {"CompPolyEvalPoly", std::to_string(validation.comp_poly_eval_poly)},
            {"CompFiltered", std::to_string(validation.comp_filtered)},
        });
    }
    if (want_module("IndelSummary")) {
        const auto insertion_deletion_ratio = ratio(metrics.indel_insertions, metrics.indel_deletions);
        write_metric_table(output, "IndelSummary", {
            {"n_SNPs", std::to_string(metrics.snp_sites)},
            {"n_singleton_SNPs", std::to_string(metrics.snp_singletons)},
            {"n_indels", std::to_string(metrics.indel_alleles)},
            {"n_singleton_indels", std::to_string(metrics.indel_singletons)},
            {"n_indel_sites", std::to_string(metrics.indel_sites)},
            {"n_multiallelic_indel_sites", std::to_string(metrics.multiallelic_indel_sites)},
            {"n_insertions", std::to_string(metrics.indel_insertions)},
            {"n_deletions", std::to_string(metrics.indel_deletions)},
            {"insertion_to_deletion_ratio", insertion_deletion_ratio},
            {"n_large_insertions", std::to_string(metrics.large_insertions)},
            {"n_large_deletions", std::to_string(metrics.large_deletions)},
            {"percent_of_sites_with_more_than_2_alleles", ratio(metrics.multiallelic_indel_sites, metrics.indel_sites)},
        });
    }
    if (want_module("IndelLengthHistogram"))
        write_indel_histogram_table(output, metrics.indel_length_histogram, metrics.indel_histogram_total);
    if (need_genotype_filter_summary) {
        write_metric_table(output, "GenotypeFilterSummary", {
            {"nCalledNotFiltered", std::to_string(metrics.genotype_called_not_filtered)},
            {"nNoCallOrFiltered", std::to_string(metrics.genotype_no_call_or_filtered)},
        });
    }
    if (want_module("MultiallelicSummary")) {
        write_metric_table(output, "MultiallelicSummary", {
            {"nProcessedLoci", std::to_string(metrics.processed_loci)},
            {"nSNPs", std::to_string(metrics.snp_sites)},
            {"nMultiSNPs", std::to_string(metrics.multiallelic_snp_sites)},
            {"processedMultiSnpRatio", ratio(metrics.multiallelic_snp_sites, metrics.processed_loci)},
            {"variantMultiSnpRatio", ratio(metrics.multiallelic_snp_sites, metrics.snp_sites)},
            {"nIndels", std::to_string(metrics.indel_sites)},
            {"nMultiIndels", std::to_string(metrics.multiallelic_indel_sites)},
            {"processedMultiIndelRatio", ratio(metrics.multiallelic_indel_sites, metrics.processed_loci)},
            {"variantMultiIndelRatio", ratio(metrics.multiallelic_indel_sites, metrics.indel_sites)},
            {"nTi", std::to_string(metrics.multiallelic_snp_ti)},
            {"nTv", std::to_string(metrics.multiallelic_snp_tv)},
            {"TiTvRatio", ratio(metrics.multiallelic_snp_ti, metrics.multiallelic_snp_tv)},
            {"knownSNPsPartial", std::to_string(metrics.known_multiallelic_snp_partial)},
            {"knownSNPsComplete", std::to_string(metrics.known_multiallelic_snp_complete)},
        });
    }
    if (!options.comps.empty() && want_module("CompOverlap")) {
        write_metric_table(output, "CompOverlap", {
            {"nEvalRecords", std::to_string(eval_records.size())},
            {"nCompRecords", std::to_string(metrics.comp_records)},
            {"nEvalCompOverlap", std::to_string(metrics.eval_comp_overlap)},
            {"nStrictAlleleMismatches", std::to_string(metrics.strict_allele_mismatch)},
            {"evalCompOverlapRate", ratio(metrics.eval_comp_overlap, eval_records.size())},
        });
    }
    const auto want_stratification = [&](const std::string& name) {
        return !options.no_standard_stratifications || contains_name(options.stratifications, name);
    };
    if (want_stratification("Contig"))
        write_stratification_table(output, "Contig", contig_strata);
    if (want_stratification("Filter"))
        write_stratification_table(output, "Filter", filter_strata);
    if (want_stratification("VariantType"))
        write_stratification_table(output, "VariantType", type_strata);
    if (want_stratification("AlleleFrequency"))
        write_stratification_table(output, "AlleleFrequency", allele_frequency_strata);
    if (want_stratification("Novelty") && !novelty_strata.empty())
        write_stratification_table(output, "Novelty", novelty_strata);
    if (want_stratification("Sample"))
        write_stratification_table(output, "Sample", sample_strata);
    if (want_stratification("Family") && !family_strata.empty())
        write_stratification_table(output, "Family", family_strata);
    if (need_genotype_concordance) {
        write_metric_table(output, "GenotypeConcordance", {
            {"nComparedGenotypes", std::to_string(metrics.compared_eval_genotypes)},
            {"nConcordantGenotypes", std::to_string(metrics.concordant_genotypes)},
            {"nDiscordantGenotypes", std::to_string(
                metrics.compared_eval_genotypes - metrics.concordant_genotypes)},
            {"genotypeConcordanceRate", ratio(metrics.concordant_genotypes,
                                               metrics.compared_eval_genotypes)},
            {"nEvalHomRef", std::to_string(metrics.genotype_concordance["HOM_REF_EVAL"])},
            {"nEvalHet", std::to_string(metrics.genotype_concordance["HET_EVAL"])},
            {"nEvalHomVar", std::to_string(metrics.genotype_concordance["HOM_VAR_EVAL"])},
            {"nEvalNoCall", std::to_string(metrics.genotype_concordance["NO_CALL_EVAL"])},
        });
    }
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing report");

    if (!options.manifest.empty()) {
        std::ofstream manifest(options.manifest);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open manifest");
        manifest << "{\"schema_version\":1,\"tool\":\"VariantEval\","
                 << "\"implementation\":\"fastgatk-variant-eval\",\"status\":\"prototype\"," 
                 << "\"compatibility\":{\"count_variants\":true,\"ti_tv\":true,"
                 << "\"variant_af_evaluator\":" << (contains_name(options.modules, "VariantAFEvaluator") ? "true" : "false") << ","
                 << "\"theta_variant_evaluator\":" << (contains_name(options.modules, "ThetaVariantEvaluator") ? "true" : "false") << ","
                 << "\"mendelian_violation_evaluator\":" << (contains_name(options.modules, "MendelianViolationEvaluator") ? "true" : "false") << ","
                 << "\"validation_report\":" << (validation_report_active ? "true" : "false") << ","
                 << "\"pedigree_trios\":" << trios.size() << ","
                 << "\"mendelian_qual_threshold\":" << decimal(options.mendelian_qual_threshold, 8) << ","
                 << "\"gatk_report\":" << (options.gatk_report ? "true" : "false") << ","
                 << "\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\","
                 << "\"variant_summary\":" << (want_module("VariantSummary") ? "true" : "false")
                 << ","
                 << "\"print_missing_comp\":" << (contains_name(options.modules, "PrintMissingComp") ? "true" : "false")
                 << ","
                 << "\"indel_summary\":" << (want_module("IndelSummary") ? "true" : "false")
                 << ",\"multiallelic_summary\":" << (want_module("MultiallelicSummary") ? "true" : "false")
                 << ",\"indel_length_histogram\":" << (want_module("IndelLengthHistogram") ? "true" : "false")
                 << ",\"genotype_filter_summary\":" << (need_genotype_filter_summary ? "true" : "false")
                 << ",\"comparison_overlap\":" << (!options.comps.empty() ? "true" : "false")
                 << ",\"keep_ac0\":" << (options.keep_ac0 ? "true" : "false")
                 << ",\"standard_stratifications\":true"
                 << ",\"novelty_stratification\":" << (!novelty_counts.empty() ? "true" : "false")
                 << ",\"sample_stratification\":" << (sample_stratification_requested && !sample_counts.empty() ? "true" : "false")
                 << ",\"family_stratification\":" << (family_stratification_active ? "true" : "false")
                 << ",\"genotype_concordance\":" << (need_genotype_concordance ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":false},\"telemetry\":{"
                 << "\"eval_records\":" << eval_records.size()
                 << ",\"comp_records\":" << comp_records.size()
                 << ",\"eval_comp_overlap\":" << metrics.eval_comp_overlap
                 << ",\"strict_allele_mismatches\":" << metrics.strict_allele_mismatch
                 << ",\"compared_genotypes\":" << metrics.compared_eval_genotypes
                 << ",\"concordant_genotypes\":" << metrics.concordant_genotypes
                 << ",\"indel_histogram_total\":" << metrics.indel_histogram_total
                 << ",\"genotype_called_not_filtered\":" << metrics.genotype_called_not_filtered
                 << ",\"genotype_no_call_or_filtered\":" << metrics.genotype_no_call_or_filtered
                 << ",\"filtered_sites\":" << metrics.filtered_sites
                 << ",\"keep_ac0\":" << (options.keep_ac0 ? "true" : "false")
                 << ",\"missing_comp_snps\":" << metrics.missing_comp_snps
                 << ",\"variant_af_called_sites\":" << metrics.variant_af_called_sites
                 << ",\"variant_af_het_sites\":" << metrics.variant_af_het_sites
                 << ",\"variant_af_hom_var_sites\":" << metrics.variant_af_hom_var_sites
                 << ",\"variant_af_hom_ref_sites\":" << metrics.variant_af_hom_ref_sites
                 << ",\"theta_num_sites\":" << metrics.theta_num_sites
                 << ",\"theta_total_het\":" << decimal(metrics.theta_total_het, 8)
                 << ",\"theta_total_avg_diffs\":" << decimal(metrics.theta_total_avg_diffs, 8)
                 << ",\"mendelian_n_variants\":" << metrics.mendelian.n_variants
                 << ",\"mendelian_n_skipped\":" << metrics.mendelian.n_skipped
                 << ",\"mendelian_n_violations\":" << metrics.mendelian.n_violations
                 << ",\"validation_n_comp\":" << validation_summary.n_comp
                 << ",\"validation_tp\":" << validation_summary.tp
                 << ",\"validation_fp\":" << validation_summary.fp
                 << ",\"validation_fn\":" << validation_summary.fn
                 << ",\"validation_tn\":" << validation_summary.tn
                 << ",\"validation_report_execution_space\":\""
                 << (validation_report_active ? Kokkos::DefaultExecutionSpace::name() : "none") << "\""
                 << ",\"novelty_stratification\":" << (!novelty_counts.empty() ? "true" : "false")
                 << ",\"novelty_known_records\":" << novelty_record_count("known")
                 << ",\"novelty_novel_records\":" << novelty_record_count("novel")
                 << ",\"stratification_tables\":" << stratification_table_count()
                 << ",\"family_names\":" << family_counts.size()
                 << ",\"sample_names\":" << sample_counts.size()
                 << ",\"territory_bases\":" << (territory == 0 ? eval_records.size() : territory)
                 << ",\"interval_list_inputs\":" << interval_stats.files
                 << ",\"interval_list_records\":" << interval_stats.records
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"execution_space\":\"Host\"},\"fallback\":{"
                 << "\"unsupported_modules\":" << (unsupported_module ? "true" : "false")
                 << ",\"unsupported_stratifications\":" << (unsupported_stratification ? "true" : "false")
                 << ",\"cloud_inputs\":false}}\n";
    }

    std::cout << "{\"tool\":\"VariantEval\",\"status\":\"prototype\","
              << "\"eval_records\":" << eval_records.size()
              << ",\"comp_records\":" << comp_records.size()
              << ",\"ti\":" << metrics.ti << ",\"tv\":" << metrics.tv
              << ",\"overlap\":" << metrics.eval_comp_overlap
              << ",\"compared_genotypes\":" << metrics.compared_eval_genotypes
              << ",\"concordant_genotypes\":" << metrics.concordant_genotypes
              << ",\"sample_names\":" << sample_counts.size()
              << ",\"territory_bases\":" << (territory == 0 ? eval_records.size() : territory)
              << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
              << ",\"stratification_tables\":" << stratification_table_count()
              << ",\"novelty_stratification\":" << (!novelty_counts.empty() ? "true" : "false") << "}\n";
    return 0;
}

#else

int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for VariantEval");
}

#endif

}  // namespace

int main(int argc, char** argv) {
    bool kokkos_initialized = false;
    try {
        const auto options = parse(argc, argv);
        Kokkos::initialize();
        kokkos_initialized = true;
        const auto result = run_tool(options, fastgatk::runtime::ResourceSnapshot::probe());
        Kokkos::finalize();
        kokkos_initialized = false;
        return result;
    } catch (const std::exception& error) {
        if (kokkos_initialized) Kokkos::finalize();
        std::cerr << "fastgatk-variant-eval: " << error.what() << '\n';
        return 2;
    }
}
