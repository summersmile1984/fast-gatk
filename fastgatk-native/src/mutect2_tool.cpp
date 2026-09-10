#include "fastgatk/calling/pipeline.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/kernels/pairhmm_kokkos.hpp"
#include "fastgatk/kernels/read_filter.hpp"
#include "fastgatk/kernels/somatic.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/somatic/f1r2_archive.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <atomic>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <sys/stat.h>

#if FASTGATK_HAS_HTSLIB
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

namespace {

enum class ReferenceConfidenceMode : std::uint8_t {
    None,
    Gvcf,
    BpResolution,
};

struct Options {
    std::string tumor;
    // GATK permits repeatable -I/--input shards.  `tumor` remains the first
    // path for the legacy single-input/streaming code; aggregate paths use
    // this ordered vector and merge them at the Host boundary.
    std::vector<std::string> tumor_inputs;
    std::string tumor_sample;
    std::string normal;
    // GATK permits --normal-sample/-normal to be supplied zero or more times.
    // Keep the historical singular field as the first-name compatibility
    // view, while the vector drives per-sample normal output and evidence.
    std::vector<std::string> normal_samples;
    std::string normal_sample;
    std::string reference;
    // GATK Mutect2 feature inputs.  These belong to the Host feature-query
    // boundary: assembly and all PairHMM/somatic numerical kernels remain
    // unchanged, while the resolved per-locus values are applied when calls
    // are selected and serialized.
    std::string germline_resource;
    std::string panel_of_normals;
    bool genotype_pon_sites = false;
    bool genotype_germline_sites = false;
    double max_population_allele_frequency = 0.01;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::string output;
    // GATK AssemblyRegionWalker debug contract.  The file is an IGV line
    // track containing both active and inactive raw profile segments; graph
    // and PairHMM continue to consume the padded active windows internally.
    std::string assembly_region_out;
    // GATK's common sites-only switch suppresses FORMAT/sample columns while
    // retaining the complete site INFO contract.  Keep the flag at the
    // writer boundary; calling and filtering still consume the same evidence.
    bool sites_only_vcf_output = false;
    // Mutect2 clips terminal inverted-repeat (ITR) artifacts by default
    // after its common read filters.  This is a Host read transformation:
    // all graph construction and PairHMM work remains in the Kokkos path.
    bool ignore_itr_artifacts = false;
    // GATK's advanced --independent-mates option supplies each read itself
    // as the SomaticGenotypingEngine grouping key instead of read name.
    bool independent_mates = false;
    // Mutect2 reference confidence is a first-class somatic mode.  Unlike
    // HaplotypeCaller gVCF, GVCF compression is by TLOD bands and BP mode
    // writes every reference context independently.
    ReferenceConfidenceMode reference_confidence_mode = ReferenceConfidenceMode::None;
    std::vector<double> gvcf_lod_bands{-2.5, -2.0, -1.5, -1.0,
                                       -0.5, 0.0, 0.5, 1.0};
    bool gvcf_lod_bands_explicit = false;
    // GATK Mutect2 --force-active keeps ActivityProfile boundaries but forces
    // every segment through AssemblyRegion graph/PairHMM processing.
    bool force_active = false;
    // GATK's common engine option controls whether a compressed VCF receives
    // a sibling Tabix index.  Keep indexing at the final writer boundary so
    // disabling it never changes caller/posterior/F1R2 computation.
    bool create_output_variant_index = true;
    // GATKTool's common VCF provenance switch.  When disabled it removes
    // both the source and command-line header records, while leaving all
    // biological/schema fields untouched.
    bool add_output_vcf_command_line = true;
    std::string stats;
    // Mutect2Engine records the number of pileup loci with this many reads in
    // its companion ``.stats`` GATK table.  This is deliberately independent
    // from --min-depth, which is a native calling/assembly control rather
    // than Mutect2's callable-site threshold.
    std::uint32_t callable_depth = 10;
    std::string f1r2;
    std::string manifest;
    std::size_t interval_file_inputs = 0;
    std::size_t interval_file_records = 0;
    std::size_t batch_records = 4096;
    std::size_t requested_batch_records = 4096;
    std::size_t effective_batch_records = 0;
    std::uint64_t adaptive_batch_reductions = 0;
    // Zero keeps the historical aggregate path.  A positive value enables
    // indexed core/halo traversal with this many reference bases per tile.
    std::size_t stream_region_size = 0;
    bool stream_indexed = false;
    std::uint64_t streamed_regions = 0;
    std::uint64_t streamed_peak_host_bytes = 0;
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
    // Mutect2's native-pair-hmm-threads default is four in GATK 4.6.2.0.
    // Keep the Kokkos launch width aligned with that public contract when no
    // explicit --threads/--native-pair-hmm-threads override is supplied.
    int threads = 4;
    // PairHMMNativeArgumentCollection defaults this public Mutect2 flag to
    // false.  Keep an explicit tool-level value instead of inheriting the
    // reusable calling pipeline's strict-double default: the latter is a
    // native implementation choice, while this option selects GATK's public
    // float/default versus double recurrence contract.
    bool native_pair_hmm_use_double_precision = false;
    double contamination = 0.0;
    // GATK's --minimum-allele-fraction changes the ALT Dirichlet
    // pseudocount used by SomaticGenotypingEngine.  Keep the default zero
    // (flat [1,1] prior) and pass explicit values through the Kokkos kernel.
    double minimum_allele_fraction = 0.0;
    double somatic_prior = 1.0e-4;
    double germline_prior = 1.0e-3;
    double artifact_prior = 1.0e-3;
    double population_allele_frequency = 5.0e-8;
    // GATK keeps the activity-region seed threshold and final emission
    // threshold separate.  The activity kernel consumes the former using
    // Mutect2's quality-aware pileup likelihood model; final emission still
    // applies the stricter threshold until the full somatic posterior path is
    // release-aligned.
    double initial_tumor_lod = 2.0;
    double tumor_lod_to_emit = 3.0;
    // GATK M2ArgumentCollection normal-lod is a log10 normal-evidence
    // threshold used by SomaticGenotypingEngine to suppress sites whose ALT
    // is sufficiently supported in the matched normal.  The native
    // SomaticLikelihoods result is already expressed in log10, so preserve
    // the option in that same domain at the output emission boundary.
    double normal_log10_odds = 2.2;
    bool population_allele_frequency_explicit = false;
    bool initial_tumor_lod_explicit = false;
    bool tumor_lod_to_emit_explicit = false;
    bool mitochondria_mode = false;
    bool graph_pruning_log_odds_threshold_explicit = false;
    bool graph_pruning_seeding_log_odds_threshold_explicit = false;
    // ReadThreadingAssembler accepts --kmer-size repeatedly.  Preserve the
    // release default until the first explicit occurrence, at which point
    // the caller-provided ordered list replaces it.
    bool graph_kmer_sizes_explicit = false;
    fastgatk::calling::Options calling;
    std::vector<std::string> read_filters;
    std::vector<std::string> disabled_read_filters;
};

// Mutect2 switches its read-to-haplotype AlleleLikelihoods to natural-log
// space before realignment, marginalization, and every BestAllele annotation.
// BestAllele::isInformative() still compares against the Java constant 0.2,
// so the equivalent margin for the native log10 PairHMM matrices is 0.2 / ln(10).
// Keep the likelihood matrix in log10 space for the Kokkos kernels; only the
// Host-side predicate changes units at this GATK API boundary.
constexpr double kGatkMutectLog10InformativeMargin = 0.08685889638065036;
// searchBestAllele() applies this earlier threshold in natural-log space
// while resolving close calls using REF's higher default priority.  Converted
// back to the native PairHMM domain it is exactly 0.2 log10.
constexpr double kGatkMutectLog10ReferenceTieMargin = 0.2;

std::optional<std::size_t> gatk_mutect_best_allele(
    const std::vector<double>& values, const double missing) {
    std::size_t best = 0;
    double best_value = -std::numeric_limits<double>::infinity();
    double second_value = -std::numeric_limits<double>::infinity();
    std::size_t observed = 0;
    for (std::size_t allele = 0; allele < values.size(); ++allele) {
        const auto value = values[allele];
        if (!std::isfinite(value) || value <= missing) continue;
        ++observed;
        // GATK's scan is stable: REF is first, and equal likelihoods retain
        // the earlier allele.  Do not move this order into a kernel because
        // it represents Java allele-list semantics rather than arithmetic.
        if (value > best_value) {
            second_value = best_value;
            best_value = value;
            best = allele;
        } else if (value > second_value) {
            second_value = value;
        }
    }
    if (observed < 2 || !std::isfinite(best_value)) return std::nullopt;

    // AlleleLikelihoods.bestAllelesBreakingTies() uses REF's priority when
    // it lies strictly within 0.2 log10 of a better ALT.  It then reports a
    // negative REF confidence, which is not informative.  Testing only the
    // final top-two margin would wrongly retain such ALT reads.
    if (best != 0 && !values.empty() && std::isfinite(values[0]) && values[0] > missing &&
        best_value - values[0] < kGatkMutectLog10ReferenceTieMargin)
        return std::nullopt;
    if (!(best_value - second_value > kGatkMutectLog10InformativeMargin))
        return std::nullopt;
    return best;
}

double effective_tumor_lod_emit(const Options& options) {
    if (options.reference_confidence_mode != ReferenceConfidenceMode::None)
        return -std::numeric_limits<double>::infinity();
    return std::max(options.initial_tumor_lod, options.tumor_lod_to_emit);
}

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name != nullptr && argument == short_name)) &&
        index + 1 < argc) return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

const char* interval_set_rule_name(const fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

std::vector<std::string> selected_normal_samples(const Options& options) {
    if (!options.normal_samples.empty()) return options.normal_samples;
    if (!options.normal_sample.empty()) return {options.normal_sample};
    return {};
}

// Keep feature resources in the C++ Host layer.  GATK's FeatureContext is a
// coordinate lookup over VariantContexts; resource VCF decoding must happen
// before a small per-locus value can be handed to the caller/output boundary.
// Neither the assembly graph nor Kokkos PairHMM/somatic arithmetic owns this
// metadata.
struct MutectFeatureVariant {
    std::string reference;
    std::vector<std::string> alternates;
    // An empty vector represents a missing or malformed AF annotation.  This
    // mirrors SomaticGenotypingEngine, which retains the configured fallback
    // when AF is absent or has the wrong Number=A cardinality.
    std::vector<double> allele_frequencies;
};

using MutectFeatureLocus = std::pair<std::string, std::int64_t>;
using MutectFeatureTable = std::map<MutectFeatureLocus,
                                    std::vector<MutectFeatureVariant>>;

struct MutectFeatureResources {
    MutectFeatureTable germline;
    MutectFeatureTable panel_of_normals;
};

#if FASTGATK_HAS_HTSLIB
MutectFeatureTable load_mutect_feature_table(const std::string& path,
                                             const bool read_allele_frequencies) {
    MutectFeatureTable result;
    if (path.empty()) return result;
    htsFile* input = bcf_open(path.c_str(), "r");
    if (input == nullptr)
        throw std::runtime_error("BAD_INPUT: cannot open Mutect2 feature VCF/BCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (header == nullptr) {
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read Mutect2 feature VCF/BCF header: " + path);
    }
    bcf1_t* record = bcf_init();
    if (record == nullptr) {
        bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate Mutect2 feature record");
    }
    while (bcf_read(input, header, record) == 0) {
        bcf_unpack(record, BCF_UN_STR);
        if (record->rid < 0 || record->n_allele < 1 || record->d.allele == nullptr) continue;
        const char* contig = bcf_hdr_id2name(header, record->rid);
        if (contig == nullptr) continue;
        MutectFeatureVariant value;
        value.reference = record->d.allele[0] == nullptr ? "" : record->d.allele[0];
        for (int allele = 1; allele < record->n_allele; ++allele)
            value.alternates.emplace_back(record->d.allele[allele] == nullptr
                                              ? "" : record->d.allele[allele]);
        if (read_allele_frequencies) {
            float* af_values = nullptr;
            int af_capacity = 0;
            const int af_count = bcf_get_info_float(header, record, "AF",
                                                    &af_values, &af_capacity);
            if (af_count == static_cast<int>(value.alternates.size())) {
                value.allele_frequencies.reserve(static_cast<std::size_t>(af_count));
                bool valid = true;
                for (int index = 0; index < af_count; ++index) {
                    if (bcf_float_is_missing(af_values[index]) ||
                        bcf_float_is_vector_end(af_values[index]) ||
                        !std::isfinite(af_values[index])) {
                        valid = false;
                        break;
                    }
                    value.allele_frequencies.push_back(af_values[index]);
                }
                if (!valid) value.allele_frequencies.clear();
            }
            std::free(af_values);
        }
        result[MutectFeatureLocus{contig, record->pos}].push_back(std::move(value));
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(input);
    return result;
}

MutectFeatureResources load_mutect_feature_resources(const Options& options) {
    return MutectFeatureResources{
        load_mutect_feature_table(options.germline_resource, true),
        load_mutect_feature_table(options.panel_of_normals, false)};
}
#endif

bool starts_with(const std::string& value, const std::string& prefix) {
    return value.size() >= prefix.size() &&
        std::equal(prefix.begin(), prefix.end(), value.begin());
}

std::optional<double> resource_allele_frequency(
    const MutectFeatureTable& resources, const std::string& contig,
    const std::int64_t position, const std::string& reference,
    const std::string& alternate) {
    const auto found = resources.find(MutectFeatureLocus{contig, position});
    if (found == resources.end()) return std::nullopt;
    std::optional<double> result;
    for (const auto& resource : found->second) {
        if (resource.allele_frequencies.size() != resource.alternates.size()) continue;
        if (!starts_with(reference, resource.reference) &&
            !starts_with(resource.reference, reference)) {
            throw std::runtime_error(
                "BAD_INPUT: Mutect2 germline resource has incompatible reference allele at " +
                contig + ":" + std::to_string(position + 1));
        }
        const auto& common_reference = reference.size() >= resource.reference.size()
            ? reference : resource.reference;
        const auto extend = [&](const std::string& allele, const std::string& source_reference) {
            return allele + common_reference.substr(source_reference.size());
        };
        const auto normalized_alternate = extend(alternate, reference);
        for (std::size_t index = 0; index < resource.alternates.size(); ++index) {
            if (normalized_alternate == extend(resource.alternates[index], resource.reference))
                result = resource.allele_frequencies[index];
        }
    }
    return result;
}

bool resource_has_locus(const MutectFeatureTable& resources, const std::string& contig,
                        const std::int64_t position) {
    return resources.find(MutectFeatureLocus{contig, position}) != resources.end();
}

std::uint8_t activity_base_index(const char base) {
    switch (base) {
        case 'A': case 'a': return 0U;
        case 'C': case 'c': return 1U;
        case 'G': case 'g': return 2U;
        case 'T': case 't': return 3U;
        default: return 4U;
    }
}

// Mutect2Engine.isActive() applies germline-resource and PON predicates
// before an AssemblyRegion is built. Resolve VCF records into the small
// caller-independent facts the shared Host pipeline needs: a common indel,
// common germline substitution reference bases, and PON membership. The
// pipeline derives the dominant tumor PileupQualBuffer bucket from reads and
// passes only the final byte-per-locus veto to Kokkos.
void configure_mutect2_activity_feature_masks(
    fastgatk::calling::Options& calling_options, const Options& options,
    const fastgatk::io::HeaderSummary& header,
    const MutectFeatureResources& feature_resources, const bool has_normal) {
    using Mask = fastgatk::calling::SomaticActivityFeatureMask;
    using MaskKey = std::pair<std::int32_t, std::int32_t>;
    calling_options.activity_genotype_germline_sites = options.genotype_germline_sites;
    calling_options.activity_feature_masks.clear();
    std::map<std::string, std::int32_t> tid_by_contig;
    for (std::size_t tid = 0; tid < header.contigs.size(); ++tid) {
        if (tid <= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
            tid_by_contig.emplace(header.contigs[tid], static_cast<std::int32_t>(tid));
    }
    std::map<MaskKey, Mask> masks;
    const auto mask_for = [&](const std::string& contig, const std::int64_t position)
        -> Mask* {
        const auto tid = tid_by_contig.find(contig);
        if (tid == tid_by_contig.end() || position < 0 ||
            position > std::numeric_limits<std::int32_t>::max())
            return nullptr;
        const MaskKey key{tid->second, static_cast<std::int32_t>(position)};
        auto [found, inserted] = masks.try_emplace(key);
        if (inserted) {
            found->second.tid = key.first;
            found->second.position = key.second;
        }
        return &found->second;
    };
    // GATK chooses the normal-pileup early veto whenever a matched normal is
    // present. The germline feature predicate is tumor-only and disabled by
    // --genotype-germline-sites.
    if (!has_normal && !options.genotype_germline_sites) {
        for (const auto& [locus, records] : feature_resources.germline) {
            for (const auto& record : records) {
                if (record.allele_frequencies.size() != record.alternates.size()) continue;
                for (std::size_t allele = 0; allele < record.alternates.size(); ++allele) {
                    if (record.allele_frequencies[allele] <
                        options.max_population_allele_frequency)
                        continue;
                    auto* mask = mask_for(locus.first, locus.second);
                    if (mask == nullptr) continue;
                    if (record.alternates[allele].size() != record.reference.size()) {
                        mask->common_germline_indel = true;
                    } else if (!record.reference.empty()) {
                        const auto reference_base = activity_base_index(record.reference.front());
                        if (reference_base < 4U)
                            mask->common_germline_substitution_reference_bases |=
                                static_cast<std::uint8_t>(1U << reference_base);
                    }
                }
            }
        }
    }
    if (!options.genotype_pon_sites) {
        for (const auto& [locus, records] : feature_resources.panel_of_normals) {
            if (records.empty()) continue;
            if (auto* mask = mask_for(locus.first, locus.second); mask != nullptr)
                mask->panel_of_normals = true;
        }
    }
    calling_options.activity_feature_masks.reserve(masks.size());
    for (const auto& [key, mask] : masks) {
        (void)key;
        calling_options.activity_feature_masks.push_back(mask);
    }
}

// M2ArgumentCollection defines every input-header sample that is not named by
// --normal-sample as a tumor.  --tumor-sample is retained as the legacy
// default/metadata selector, but is not a request to discard other tumor
// samples: SomaticGenotypingEngine combines all of them for TLOD and writes a
// FORMAT column for each one.
std::vector<std::string> selected_tumor_samples(
    const Options& options, const fastgatk::io::HeaderSummary& header) {
    const auto normal_samples = selected_normal_samples(options);
    std::vector<std::string> tumors;
    tumors.reserve(header.samples.size());
    for (const auto& sample : header.samples) {
        if (std::find(normal_samples.begin(), normal_samples.end(), sample) == normal_samples.end())
            tumors.push_back(sample);
    }
    std::sort(tumors.begin(), tumors.end());
    tumors.erase(std::unique(tumors.begin(), tumors.end()), tumors.end());
    // GATK also accepts a deliberately self-matched invocation in which an
    // explicit --tumor-sample is repeated as --normal-sample (the compact
    // contract uses it to exercise normal-lod behavior on one BAM).  In that
    // compatibility case the explicit tumor role is authoritative; do not
    // turn the input matrix empty merely because both role selectors name it.
    if (tumors.empty() && !options.tumor_sample.empty() &&
        std::find(header.samples.begin(), header.samples.end(), options.tumor_sample) !=
            header.samples.end())
        tumors.push_back(options.tumor_sample);
    return tumors;
}

// These classes are representable by the shared Kokkos ordinal mask.  Keep
// unknown custom filters fail-closed; the dispatcher can then preserve argv
// and invoke the configured GATK fallback explicitly.
void enable_read_filter(fastgatk::calling::Options& calling,
                        const std::string& filter_name) {
    if (filter_name == "MappingQualityReadFilter") {
        // Uses the shared --minimum-mapping-quality setting.
    } else if (filter_name == "MappingQualityAvailableReadFilter") {
        calling.exclude_mapping_quality_unavailable = true;
    } else if (filter_name == "MappingQualityNotZeroReadFilter") {
        calling.exclude_mapping_quality_zero = true;
    } else if (filter_name == "NotDuplicateReadFilter") {
        calling.exclude_duplicates = true;
    } else if (filter_name == "MappedReadFilter") {
        calling.exclude_unmapped = true;
    } else if (filter_name == "NotSecondaryAlignmentReadFilter") {
        calling.exclude_secondary = true;
    } else if (filter_name == "NotSupplementaryAlignmentReadFilter") {
        calling.exclude_supplementary = true;
    } else if (filter_name == "PassesVendorQualityCheckReadFilter") {
        calling.exclude_qcfail = true;
    } else if (filter_name == "GoodCigarReadFilter") {
        calling.require_good_cigar = true;
    } else if (filter_name == "WellformedReadFilter") {
        calling.require_good_cigar = true;
        calling.require_no_n_cigar = true;
        calling.require_read_group = true;
    } else if (filter_name == "NonZeroReferenceLengthAlignmentReadFilter") {
        calling.require_nonzero_reference_span = true;
    } else if (filter_name == "NonChimericOriginalAlignmentReadFilter") {
        calling.require_non_chimeric_original_alignment = true;
    } else if (filter_name == "ReadLengthReadFilter") {
        calling.require_read_length = true;
    } else {
        throw std::runtime_error("UNSUPPORTED_PARAMETER: --read-filter " + filter_name);
    }
}

void disable_read_filter(fastgatk::calling::Options& calling,
                         const std::string& filter_name) {
    if (filter_name == "NotDuplicateReadFilter") {
        calling.exclude_duplicates = false;
    } else if (filter_name == "MappingQualityAvailableReadFilter") {
        calling.exclude_mapping_quality_unavailable = false;
    } else if (filter_name == "MappingQualityNotZeroReadFilter") {
        calling.exclude_mapping_quality_zero = false;
    } else if (filter_name == "MappingQualityReadFilter") {
        calling.minimum_mapping_quality = 0;
    } else if (filter_name == "MappedReadFilter") {
        calling.exclude_unmapped = false;
    } else if (filter_name == "NotSecondaryAlignmentReadFilter") {
        calling.exclude_secondary = false;
    } else if (filter_name == "NotSupplementaryAlignmentReadFilter") {
        calling.exclude_supplementary = false;
    } else if (filter_name == "PassesVendorQualityCheckReadFilter") {
        calling.exclude_qcfail = false;
    } else if (filter_name == "GoodCigarReadFilter") {
        calling.require_good_cigar = false;
    } else if (filter_name == "WellformedReadFilter") {
        calling.require_good_cigar = false;
        calling.require_no_n_cigar = false;
        calling.require_read_group = false;
    } else if (filter_name == "NonZeroReferenceLengthAlignmentReadFilter") {
        calling.require_nonzero_reference_span = false;
    } else if (filter_name == "NonChimericOriginalAlignmentReadFilter") {
        calling.require_non_chimeric_original_alignment = false;
    } else if (filter_name == "ReadLengthReadFilter") {
        calling.require_read_length = false;
    } else {
        throw std::runtime_error("UNSUPPORTED_PARAMETER: --disable-read-filter " + filter_name);
    }
}

void clear_tool_default_read_filters(fastgatk::calling::Options& calling) {
    calling.minimum_mapping_quality = 0;
    calling.exclude_mapping_quality_unavailable = false;
    calling.exclude_mapping_quality_zero = false;
    calling.exclude_duplicates = false;
    calling.exclude_unmapped = false;
    calling.exclude_secondary = false;
    calling.exclude_supplementary = false;
    calling.exclude_qcfail = false;
    calling.require_good_cigar = false;
    calling.require_nonzero_reference_span = false;
    calling.require_non_chimeric_original_alignment = false;
    calling.require_no_n_cigar = false;
    calling.require_read_group = false;
    calling.require_read_length = false;
}

void resolve_read_filter_options(Options& options) {
    if (!options.calling.disable_tool_default_read_filters) return;
    const auto configured_mapping_quality = options.calling.minimum_mapping_quality;
    clear_tool_default_read_filters(options.calling);
    for (const auto& filter_name : options.read_filters) {
        if (filter_name == "MappingQualityReadFilter")
            options.calling.minimum_mapping_quality = configured_mapping_quality;
        enable_read_filter(options.calling, filter_name);
    }
    for (const auto& filter_name : options.disabled_read_filters)
        disable_read_filter(options.calling, filter_name);
    options.calling.apply_read_filters = true;
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

// JSON does not admit IEEE NaN or infinity tokens.  PairHMM diagnostics can
// legitimately retain an unscored candidate at a narrow interval boundary;
// keep that absence explicit in native telemetry rather than emitting an
// invalid sidecar that downstream JSON readers cannot consume.
std::string json_number_or_null(const double value) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream text;
    text << std::setprecision(17) << value;
    return text.str();
}

std::string normal_sample_list_json(const Options& options) {
    auto names = selected_normal_samples(options);
    if (names.empty() && !options.normal.empty()) names.push_back("NORMAL");
    std::ostringstream text;
    text << '[';
    for (std::size_t index = 0; index < names.size(); ++index) {
        if (index != 0) text << ',';
        text << '"' << json_escape(names[index]) << '"';
    }
    text << ']';
    return text.str();
}

Options parse(int argc, char** argv) {
    Options options;
    // ReadThreadingAssembler's release defaults are the retry list [10,25],
    // not a single k=10 attempt.  Keep the scalar field for source
    // compatibility and pass the complete list to the shared graph API.
    options.calling.graph_kmer_size = 10;
    options.calling.graph_kmer_sizes = {10, 25};
    // M2ArgumentCollection overrides HC's MNP default to one base of phase
    // distance.  Explicit --max-mnp-distance/--mnp-dist can still set zero.
    options.calling.max_mnp_distance = 1;
    // Match Mutect2's StandardMutect2ReadFilters. GoodCigar,
    // NonZeroReferenceLengthAlignment, NonChimericOriginalAlignment and
    // ReadLength(30, infinity) are
    // implicit GATK defaults (not
    // optional quality gates), so enable them at the CLI boundary while
    // keeping the shared library's permissive defaults for direct callers.
    options.calling.require_good_cigar = true;
    options.calling.require_nonzero_reference_span = true;
    options.calling.require_non_chimeric_original_alignment = true;
    options.calling.require_no_n_cigar = true;
    options.calling.require_read_group = true;
    options.calling.require_read_length = true;
    options.calling.min_read_length = 30;
    options.calling.max_read_length = std::numeric_limits<std::uint32_t>::max();
    options.calling.exclude_mapping_quality_unavailable = true;
    options.calling.exclude_mapping_quality_zero = true;
    // MutectReadThreadingAssemblerArgumentCollection enables adaptive
    // pruning by default; --disable-adaptive-pruning switches back to the
    // fixed min-pruning path. HC keeps its separate GATK default (disabled).
    options.calling.graph_use_adaptive_pruning = true;
    // Mutect2 owns its somatic posterior thresholding downstream; the shared
    // HC call-confidence gate must not discard low-VAF candidates here.
    options.calling.standard_confidence_for_calling = 0.0;
    options.calling.somatic_mode = true;
    // GATK Mutect2 overrides AssemblyRegionWalker's regular positional
    // sampler with MutectDownsampler.  Keep this tool policy explicit as it
    // crosses the Host read-filter/assembly boundary.
    options.calling.mutect2_downsampling = true;
    // Mutect2's PairHMMLikelihoodCalculationEngine removes poorly-modelled
    // evidence after its normalized PairHMM matrix is built.  Candidate
    // emission remains a later SomaticGenotypingEngine decision, but this
    // read-level gate is still part of the likelihood/evidence contract and
    // must remain enabled here (as it is in the shared Options default).
    options.calling.pairhmm_filter_poorly_modeled_reads = true;
    options.calling.activity_quality_aware_somatic = true;
    options.calling.activity_initial_tumor_log10_odds = options.initial_tumor_lod;
    options.calling.activity_pcr_snv_quality = 40;
    options.calling.activity_multiple_substitution_quality_correction = 5;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-mutect2 (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 tumor BAM/CRAM (repeatable for aggregate shards)\n"
                         "      --tumor-sample NAME          tumor sample name (default: first header sample)\n"
                         "      --normal-input FILE          normal BAM/CRAM (optional)\n"
                         "      -normal, --normal-sample NAME normal sample (repeatable)\n"
                         "  -R, --reference FILE             reference FASTA\n"
                         "  -L, --intervals REGION           contig:start-end\n"
                         "      -isr, --interval-set-rule R  UNION (default) or INTERSECTION\n"
                         "  -O, --output FILE                somatic VCF/VCF.GZ\n"
                         "      -ERC, --emit-ref-confidence MODE  NONE, GVCF, or BP_RESOLUTION (default NONE)\n"
                         "      -LODB, --gvcf-lod-band F     repeatable exclusive upper TLOD band bound\n"
                         "      --assembly-region-out FILE  IGV activity/profile regions\n"
                         "      --sites-only-vcf-output B   omit FORMAT/sample columns (default false)\n"
                         "      --ignore-itr-artifacts[=BOOL]  disable default ITR palindrome artifact clipping\n"
                         "      --independent-mates[=BOOL]  let paired reads support haplotypes independently\n"
                         "      --force-active[=BOOL]       force all activity segments active (default false)\n"
                         "      --create-output-variant-index[=B]  write .tbi for compressed VCF (default true)\n"
                         "      --add-output-vcf-command-line[=BOOL]  write VCF source/provenance headers (default true)\n"
                         "      --stats FILE                 Mutect2 stats sidecar\n"
                         "      --callable-depth N           minimum pileup depth recorded in Mutect2 .stats (default 10)\n"
                         "      --f1r2-tar-gz FILE            F1R2 sidecar (.tar.gz standard; .tsv legacy)\n"
                         "      --stream-by-region N         indexed core/halo tile size in bases\n"
                         "      --native-pair-hmm-threads N  Kokkos PairHMM thread count (default 4; --threads alias)\n"
                         "      --native-pair-hmm-use-double-precision[=BOOL]  use double PairHMM recurrence (default false)\n"
                         "      --contamination F             contamination fraction [0,1)\n"
                         "      --minimum-allele-fraction F  soft ALT fraction prior [0,1) (-min-AF alias)\n"
                         "      --somatic-prior F             somatic posterior prior\n"
                         "      --germline-prior F            germline posterior prior\n"
                         "      --artifact-prior F            orientation-artifact prior\n"
                         "      --germline-resource FILE     population AF VCF/BCF feature input\n"
                         "      -pon, --panel-of-normals FILE panel-of-normals VCF/BCF feature input\n"
                         "      --genotype-pon-sites[=BOOL] emit PON sites and annotate PON\n"
                         "      --genotype-germline-sites[=BOOL] emit high-AF tumor-only sites\n"
                         "      -max-af, --max-population-af F  tumor-only germline activity AF ceiling (default 0.01)\n"
                         "      --af-of-alleles-not-in-resource F  population AF fallback (5e-8 tumor-only, 1e-6 tumor-normal; -default-af alias)\n"
                         "      --initial-tumor-lod F         active-region seed TLOD (default 2; -init-lod alias)\n"
                         "      --tumor-lod-to-emit F         minimum TLOD for emitted candidate (default 3; -emit-lod alias)\n"
                         "      --normal-lod F                matched-normal ALT log10-odds emission threshold (default 2.2)\n"
                         "      --pcr-snv-qual Q             PCR SNV quality cap for overlapping fragments (default 40)\n"
                         "      --pcr-indel-qual Q           PCR indel quality cap for overlapping fragments (default 40)\n"
                         "      --base-qual-correction-factor Q  substitution base-quality correction (default 5)\n"
                         "      --mitochondria-mode[=BOOL]   GATK mitochondrial defaults (optional boolean)\n"
                         "      --active-probability-threshold P  AssemblyRegion activity threshold\n"
                         "      --assembly-region-padding N  AssemblyRegion halo/padding\n"
                         "      --min-assembly-region-size N  minimum AssemblyRegion size (default 50)\n"
                         "      --max-assembly-region-size N  hard AssemblyRegion length cap\n"
                         "      --max-prob-propagation-distance N  active-locus grouping distance\n"
                         "      --minimum-mapping-quality Q  GATK read-filter MAPQ floor (default 20)\n"
                         "      --min-read-length N          ReadLengthReadFilter lower bound\n"
                         "      --max-read-length N          ReadLengthReadFilter upper bound\n"
                         "      --base-quality-score-threshold Q  PairHMM base-quality floor (default 18)\n"
                         "      --min-base-quality-score Q   minimum base quality for calling (default 10)\n"
                         "      --disable-cap-base-qualities-to-map-quality  disable PairHMM MAPQ cap\n"
                         "      --phred-scaled-global-read-mismapping-rate Q  per-read likelihood cap (default 45; negative disables)\n"
                         "      --pcr-indel-model M       NONE, HOSTILE, AGGRESSIVE, or CONSERVATIVE (default)\n"
                         "      --flow-assembly-collapse-hmer-size N  cap flow haplotype homopolymers (-1=auto, 0=off)\n"
                         "      --flow-assembly-collapse-partial-mode  stop HMER restoration at reference mismatch\n"
                         "      --include-duplicates         keep duplicate reads (debug/compatibility)\n"
                         "      --dont-use-soft-clipped-bases  exclude soft-clipped bases from assembly evidence\n"
                         "      --max-reads-per-locus N     positional reservoir downsampling cap (0=unbounded)\n"
                         "      --max-reads-per-alignment-start N  GATK positional downsampling cap (default 50)\n"
                         "      --downsampling-stride N     Mutect2 downsampling pool width (default 1; -stride alias)\n"
                         "      --max-suspicious-reads-per-alignment-start N  reject stride after N MAPQ<=50 reads (0=off)\n"
                         "      --downsampling-seed N      Java-compatible reservoir seed\n"
                         "      --allele-informative-reads-overlap-margin N  retain-evidence overlap margin (default 2)\n"
                         "      --kmer-size K                repeatable graph k-mer size (positive uint32, default 10,25)\n"
                         "      --dont-increase-kmer-sizes-for-cycles  disable GATK cyclic graph k-mer retries\n"
                         "      --min-kmer-count N           graph read-kmer support floor\n"
                         "      --min-pruning N              GATK graph branch pruning support floor\n"
                         "      --num-pruning-samples N      samples that must pass min-pruning (default 1)\n"
                         "      --adaptive-pruning           enable GATK-shaped likelihood chain pruning\n"
                         "      --adaptive-pruning-initial-error-rate P  adaptive pruning error prior\n"
                         "      --pruning-lod-threshold D    adaptive pruning log10 threshold\n"
                         "      --pruning-seeding-lod-threshold D  adaptive seeding log10 threshold\n"
                         "      --max-unpruned-variants N   adaptive pruning variant cap\n"
                         "      --disable-adaptive-pruning  use fixed min-pruning only\n"
                         "      --linked-de-bruijn-graph     retain linked de-Bruijn topology (experimental)\n"
                         "      --disable-artificial-haplotype-recovery  disable linked-graph recovery\n"
                         "      --enable-legacy-graph-cycle-detection  use legacy cycle policy\n"
                         "      --min-dangling-branch-length N  minimum retained dangling branch length (default 4)\n"
                         "      --min-dangling-matching-bases N  exact-prefix gate; -1 uses GATK legacy recovery (default -1)\n"
                         "      --allow-non-unique-kmers-in-ref  allow reference k-mer repeats (GATK default: reject/retry)\n"
                         "      --recover-all-dangling-branches  recover connected dangling graph forks\n"
                         "      --error-correct-reads        correct low-frequency assembly k-mers (off by default)\n"
                         "      --error-correction-log-odds D  hidden GATK pileup correction threshold (default -inf)\n"
                         "      --kmer-length-for-read-error-correction K  correction k-mer length (default 25)\n"
                         "      --min-observations-for-kmer-to-be-solid N  solid correction k-mer floor (default 20)\n"
                         "      --max-haplotype-paths N      bounded graph path count\n"
                         "      --max-num-haplotypes-in-population N  GATK alias for graph path cap\n"
                         "      --max-haplotype-depth N      bounded graph path depth\n"
                         "      --max-haplotype-combination-alleles N  complete local haplotype set cap (1..16)\n"
                         "      --max-mnp-distance N         merge phased graph SNPs into an MNP (default 1)\n"
                         "      --mnp-dist N                 short alias for --max-mnp-distance\n";
            std::exit(0);
        } else if (is_option(argument, "--input") || argument == "-I") {
            const auto input = require_value(index, argc, argv, argument, "--input", "-I");
            if (options.tumor.empty()) options.tumor = input;
            options.tumor_inputs.push_back(input);
        }
        else if (is_option(argument, "--tumor-sample"))
            options.tumor_sample = require_value(index, argc, argv, argument, "--tumor-sample");
        else if (is_option(argument, "--normal-input"))
            options.normal = require_value(index, argc, argv, argument, "--normal-input");
        else if (argument == "-normal" || is_option(argument, "--normal") ||
                 is_option(argument, "--normal-sample"))
        {
            const auto sample = require_value(
                index, argc, argv, argument,
                argument.rfind("--normal-sample", 0) == 0 ? "--normal-sample" : "--normal",
                "-normal");
            if (std::find(options.normal_samples.begin(), options.normal_samples.end(), sample) ==
                options.normal_samples.end())
                options.normal_samples.push_back(sample);
            if (options.normal_sample.empty()) options.normal_sample = sample;
        }
        else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--intervals") || is_option(argument, "--region") || argument == "-L")
            options.regions.push_back(require_value(index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                argument.rfind("--region", 0) == 0 ? "--region" : "--intervals", "-L"));
        else if (argument == "-isr" || is_option(argument, "--interval-set-rule")) {
            auto rule = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char value) {
                return static_cast<char>(std::toupper(value));
            });
            if (rule == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (rule == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument(
                "BAD_INPUT: --interval-set-rule must be UNION or INTERSECTION");
        }
        else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--assembly-region-out"))
            options.assembly_region_out = require_value(
                index, argc, argv, argument, "--assembly-region-out");
        else if (argument == "--force-active" ||
                 argument.rfind("--force-active=", 0) == 0)
            options.force_active = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--force-active");
        else if (is_option(argument, "--sites-only-vcf-output")) {
            const auto value = require_value(index, argc, argv, argument,
                                             "--sites-only-vcf-output");
            if (value == "true" || value == "TRUE" || value == "1")
                options.sites_only_vcf_output = true;
            else if (value == "false" || value == "FALSE" || value == "0")
                options.sites_only_vcf_output = false;
            else
                throw std::invalid_argument(
                    "--sites-only-vcf-output must be true or false");
        }
        else if (argument == "-ERC" || argument.rfind("-ERC=", 0) == 0 ||
                 is_option(argument, "--emit-ref-confidence")) {
            auto mode = argument.rfind("-ERC=", 0) == 0
                ? argument.substr(5)
                : require_value(index, argc, argv, argument,
                                "--emit-ref-confidence", "-ERC");
            std::transform(mode.begin(), mode.end(), mode.begin(),
                           [](const unsigned char value) {
                               return static_cast<char>(std::toupper(value));
                           });
            if (mode == "NONE")
                options.reference_confidence_mode = ReferenceConfidenceMode::None;
            else if (mode == "GVCF")
                options.reference_confidence_mode = ReferenceConfidenceMode::Gvcf;
            else if (mode == "BP_RESOLUTION")
                options.reference_confidence_mode = ReferenceConfidenceMode::BpResolution;
            else
                throw std::invalid_argument(
                    "--emit-ref-confidence must be NONE, GVCF, or BP_RESOLUTION");
        }
        else if (argument == "-LODB" || argument.rfind("-LODB=", 0) == 0 ||
                 is_option(argument, "--gvcf-lod-band")) {
            const auto value = argument.rfind("-LODB=", 0) == 0
                ? argument.substr(6)
                : require_value(index, argc, argv, argument,
                                "--gvcf-lod-band", "-LODB");
            if (!options.gvcf_lod_bands_explicit) {
                options.gvcf_lod_bands.clear();
                options.gvcf_lod_bands_explicit = true;
            }
            options.gvcf_lod_bands.push_back(std::stod(value));
        }
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0)
            options.create_output_variant_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        else if (argument == "--add-output-vcf-command-line" ||
                 argument.rfind("--add-output-vcf-command-line=", 0) == 0)
            options.add_output_vcf_command_line = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--add-output-vcf-command-line");
        else if (is_option(argument, "--stats"))
            options.stats = require_value(index, argc, argv, argument, "--stats");
        else if (is_option(argument, "--callable-depth")) {
            const auto value = require_value(index, argc, argv, argument, "--callable-depth");
            const auto parsed = std::stoll(value);
            if (parsed < 1 || parsed > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--callable-depth must be a positive 32-bit integer");
            options.callable_depth = static_cast<std::uint32_t>(parsed);
        }
        else if (is_option(argument, "--f1r2-tar-gz"))
            options.f1r2 = require_value(index, argc, argv, argument, "--f1r2-tar-gz");
        else if (is_option(argument, "--contamination") || is_option(argument, "--contamination-fraction"))
            options.contamination = std::stod(require_value(
                index, argc, argv, argument,
                argument.rfind("--contamination-fraction", 0) == 0
                    ? "--contamination-fraction" : "--contamination"));
        else if (is_option(argument, "--minimum-allele-fraction") ||
                 argument == "-min-AF" || argument.rfind("-min-AF=", 0) == 0) {
            if (argument.rfind("-min-AF=", 0) == 0)
                options.minimum_allele_fraction = std::stod(argument.substr(8));
            else
                options.minimum_allele_fraction = std::stod(require_value(
                    index, argc, argv, argument, "--minimum-allele-fraction", "-min-AF"));
        }
        else if (is_option(argument, "--somatic-prior"))
            options.somatic_prior = std::stod(require_value(index, argc, argv, argument, "--somatic-prior"));
        else if (is_option(argument, "--germline-prior"))
            options.germline_prior = std::stod(require_value(index, argc, argv, argument, "--germline-prior"));
        else if (is_option(argument, "--artifact-prior"))
            options.artifact_prior = std::stod(require_value(index, argc, argv, argument, "--artifact-prior"));
        else if (is_option(argument, "--germline-resource"))
            options.germline_resource = require_value(
                index, argc, argv, argument, "--germline-resource");
        else if (argument == "-pon" || is_option(argument, "--panel-of-normals"))
            options.panel_of_normals = require_value(
                index, argc, argv, argument, "--panel-of-normals", "-pon");
        else if (argument == "--genotype-pon-sites" ||
                 argument.rfind("--genotype-pon-sites=", 0) == 0)
            options.genotype_pon_sites = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--genotype-pon-sites");
        else if (argument == "--genotype-germline-sites" ||
                 argument.rfind("--genotype-germline-sites=", 0) == 0)
            options.genotype_germline_sites = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--genotype-germline-sites");
        else if (argument == "-max-af" || is_option(argument, "--max-population-af"))
            options.max_population_allele_frequency = std::stod(require_value(
                index, argc, argv, argument, "--max-population-af", "-max-af"));
        else if (is_option(argument, "--af-of-alleles-not-in-resource") ||
                 is_option(argument, "--population-af") ||
                 is_option(argument, "--default-af") || argument == "-default-af") {
            const auto name = argument == "-default-af" ||
                              argument.rfind("--default-af", 0) == 0
                                  ? "--default-af"
                                  : (argument.rfind("--population-af", 0) == 0
                                         ? "--population-af" : "--af-of-alleles-not-in-resource");
            options.population_allele_frequency = std::stod(require_value(
                index, argc, argv, argument,
                name, "-default-af"));
            options.population_allele_frequency_explicit = true;
        }
        else if (is_option(argument, "--initial-tumor-lod") || argument == "-init-lod") {
            options.initial_tumor_lod = std::stod(require_value(
                index, argc, argv, argument, "--initial-tumor-lod", "-init-lod"));
            options.initial_tumor_lod_explicit = true;
        }
        else if (is_option(argument, "--tumor-lod-to-emit") || argument == "-emit-lod") {
            options.tumor_lod_to_emit = std::stod(require_value(
                index, argc, argv, argument, "--tumor-lod-to-emit", "-emit-lod"));
            options.tumor_lod_to_emit_explicit = true;
        }
        else if (is_option(argument, "--normal-lod"))
            options.normal_log10_odds = std::stod(require_value(
                index, argc, argv, argument, "--normal-lod"));
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (is_option(argument, "--batch-records"))
            options.batch_records = std::stoull(require_value(index, argc, argv, argument, "--batch-records"));
        else if (is_option(argument, "--stream-by-region")) {
            options.stream_region_size = std::stoull(require_value(
                index, argc, argv, argument, "--stream-by-region"));
            if (options.stream_region_size == 0)
                throw std::invalid_argument("--stream-by-region must be positive");
        }
        else if (is_option(argument, "--threads") || is_option(argument, "--native-pair-hmm-threads"))
            options.threads = std::stoi(require_value(index, argc, argv, argument,
                argument.rfind("--native-pair-hmm-threads", 0) == 0
                    ? "--native-pair-hmm-threads" : "--threads"));
        else if (argument == "--native-pair-hmm-use-double-precision" ||
                 argument.rfind("--native-pair-hmm-use-double-precision=", 0) == 0)
            options.native_pair_hmm_use_double_precision =
                fastgatk::native::parse_optional_boolean(
                    index, argc, argv, argument,
                    "--native-pair-hmm-use-double-precision");
        else if (is_option(argument, "--min-depth"))
            options.calling.min_depth = std::stoull(require_value(index, argc, argv, argument, "--min-depth"));
        else if (is_option(argument, "--min-alt-support"))
            options.calling.min_alt_support = std::stoull(require_value(index, argc, argv, argument, "--min-alt-support"));
        else if (is_option(argument, "--active-probability-threshold")) {
            options.calling.active_probability_threshold = std::stod(require_value(
                index, argc, argv, argument, "--active-probability-threshold"));
            if (!std::isfinite(options.calling.active_probability_threshold) ||
                options.calling.active_probability_threshold < 0.0 ||
                options.calling.active_probability_threshold > 1.0)
                throw std::invalid_argument(
                    "--active-probability-threshold must be finite and in [0,1]");
        }
        else if (is_option(argument, "--assembly-region-padding"))
            options.calling.assembly_region_padding = static_cast<std::uint32_t>(std::stoul(
                require_value(index, argc, argv, argument, "--assembly-region-padding")));
        else if (is_option(argument, "--min-assembly-region-size"))
            options.calling.min_assembly_region_size = static_cast<std::uint32_t>(std::stoul(
                require_value(index, argc, argv, argument, "--min-assembly-region-size")));
        else if (is_option(argument, "--max-assembly-region-size"))
            options.calling.max_assembly_region_size = static_cast<std::uint32_t>(std::stoul(
                require_value(index, argc, argv, argument, "--max-assembly-region-size")));
        else if (is_option(argument, "--max-prob-propagation-distance"))
            options.calling.max_probability_propagation_distance = static_cast<std::uint32_t>(std::stoul(
                require_value(index, argc, argv, argument, "--max-prob-propagation-distance")));
        else if (is_option(argument, "--minimum-mapping-quality")) {
            const auto value = std::stoul(require_value(index, argc, argv, argument,
                                                        "--minimum-mapping-quality"));
            if (value > 255) throw std::invalid_argument("--minimum-mapping-quality must be <=255");
            options.calling.minimum_mapping_quality = static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--min-read-length")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--min-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-read-length is out of range");
            options.calling.min_read_length = static_cast<std::uint32_t>(value);
            options.calling.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        }
        else if (is_option(argument, "--max-read-length")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--max-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-read-length is out of range");
            options.calling.max_read_length = static_cast<std::uint32_t>(value);
            options.calling.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        }
        else if (is_option(argument, "--base-quality-score-threshold")) {
            const auto value = std::stoul(require_value(index, argc, argv, argument,
                                                         "--base-quality-score-threshold"));
            if (value < 6 || value > 255)
                throw std::invalid_argument("--base-quality-score-threshold must be in [6,255]");
            options.calling.pairhmm_base_quality_score_threshold =
                static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--pcr-snv-qual")) {
            const auto value = std::stoul(require_value(index, argc, argv, argument,
                                                         "--pcr-snv-qual"));
            if (value > 255)
                throw std::invalid_argument("--pcr-snv-qual must be <=255");
            options.calling.activity_pcr_snv_quality = static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--pcr-indel-qual")) {
            const auto value = std::stoul(require_value(index, argc, argv, argument,
                                                         "--pcr-indel-qual"));
            if (value > 255)
                throw std::invalid_argument("--pcr-indel-qual must be <=255");
            options.calling.overlapping_pcr_indel_quality =
                static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--base-qual-correction-factor")) {
            const auto value = std::stoul(require_value(index, argc, argv, argument,
                                                         "--base-qual-correction-factor"));
            if (value > 255)
                throw std::invalid_argument("--base-qual-correction-factor must be <=255");
            options.calling.activity_multiple_substitution_quality_correction =
                static_cast<std::uint8_t>(value);
        }
        else if (argument == "-mbq" || is_option(argument, "--min-base-quality-score") ||
                 is_option(argument, "--min-base-quality")) {
            const auto name = argument == "-mbq" ||
                              argument.rfind("--min-base-quality-score", 0) == 0
                                  ? "--min-base-quality-score" : "--min-base-quality";
            const auto value = std::stoul(require_value(index, argc, argv, argument, name, "-mbq"));
            if (value > 255)
                throw std::invalid_argument("--min-base-quality-score must be <=255");
            options.calling.min_base_quality = static_cast<std::uint8_t>(value);
        }
        else if (argument == "--disable-cap-base-qualities-to-map-quality")
            options.calling.pairhmm_disable_cap_base_qualities_to_mapq = true;
        else if (is_option(argument, "--phred-scaled-global-read-mismapping-rate")) {
            options.calling.phred_scaled_global_read_mismapping_rate = std::stod(
                require_value(index, argc, argv, argument,
                              "--phred-scaled-global-read-mismapping-rate"));
            if (!std::isfinite(options.calling.phred_scaled_global_read_mismapping_rate))
                throw std::invalid_argument(
                    "--phred-scaled-global-read-mismapping-rate must be finite");
        }
        else if (is_option(argument, "--pcr-indel-model")) {
            auto model = require_value(index, argc, argv, argument, "--pcr-indel-model");
            std::transform(model.begin(), model.end(), model.begin(),
                           [](const unsigned char value) { return static_cast<char>(std::toupper(value)); });
            if (model == "NONE") {
                options.calling.pairhmm_conservative_indel_model = false;
                options.calling.pairhmm_pcr_indel_model = "NONE";
                options.calling.pairhmm_pcr_error_rate_factor = 0.0;
            } else if (model == "HOSTILE" || model == "AGGRESSIVE" || model == "CONSERVATIVE") {
                options.calling.pairhmm_conservative_indel_model = true;
                options.calling.pairhmm_pcr_indel_model = model;
                options.calling.pairhmm_pcr_error_rate_factor =
                    model == "HOSTILE" ? 1.0 : (model == "AGGRESSIVE" ? 2.0 : 3.0);
            } else throw std::invalid_argument(
                "--pcr-indel-model must be NONE, HOSTILE, AGGRESSIVE, or CONSERVATIVE");
        }
        else if (is_option(argument, "--flow-assembly-collapse-hmer-size")) {
            const auto value = std::stoll(require_value(index, argc, argv, argument,
                                                         "--flow-assembly-collapse-hmer-size"));
            if (value < -1 || value > 255)
                throw std::invalid_argument(
                    "--flow-assembly-collapse-hmer-size must be in [-1,255]");
            options.calling.flow_assembly_collapse_hmer_size =
                static_cast<std::int32_t>(value);
        }
        else if (argument == "--flow-assembly-collapse-partial-mode")
            options.calling.flow_assembly_collapse_partial_mode = true;
        else if (is_option(argument, "--expected-mismatch-rate-for-read-disqualification")) {
            options.calling.pairhmm_expected_error_rate_per_base = std::stod(
                require_value(index, argc, argv, argument,
                              "--expected-mismatch-rate-for-read-disqualification"));
            if (!std::isfinite(options.calling.pairhmm_expected_error_rate_per_base) ||
                options.calling.pairhmm_expected_error_rate_per_base < 0.0)
                throw std::invalid_argument(
                    "--expected-mismatch-rate-for-read-disqualification must be finite and non-negative");
        }
        else if (argument == "--include-duplicates")
            options.calling.exclude_duplicates = false;
        else if (argument == "--dont-use-soft-clipped-bases" ||
                 argument.rfind("--dont-use-soft-clipped-bases=", 0) == 0)
            options.calling.use_soft_clipped_bases = !fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--dont-use-soft-clipped-bases");
        else if (argument == "--ignore-itr-artifacts" ||
                 argument.rfind("--ignore-itr-artifacts=", 0) == 0)
            options.ignore_itr_artifacts = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--ignore-itr-artifacts");
        else if (argument == "--independent-mates" ||
                 argument.rfind("--independent-mates=", 0) == 0) {
            options.independent_mates = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--independent-mates");
            options.calling.somatic_independent_mates = options.independent_mates;
        }
        else if (argument == "--do-not-correct-overlapping-quality" ||
                 argument.rfind("--do-not-correct-overlapping-quality=", 0) == 0 ||
                 argument == "--do-not-correct-overlapping-base-qualities" ||
                 argument.rfind("--do-not-correct-overlapping-base-qualities=", 0) == 0) {
            const char* name = argument.rfind("--do-not-correct-overlapping-base-qualities", 0) == 0
                ? "--do-not-correct-overlapping-base-qualities"
                : "--do-not-correct-overlapping-quality";
            options.calling.do_not_correct_overlapping_base_qualities =
                fastgatk::native::parse_optional_boolean(index, argc, argv, argument, name);
        }
        else if (is_option(argument, "--max-reads-per-locus") ||
                 is_option(argument, "--max-reads-per-alignment-start")) {
            const auto name = argument.rfind("--max-reads-per-alignment-start", 0) == 0
                ? "--max-reads-per-alignment-start" : "--max-reads-per-locus";
            options.calling.max_reads_per_locus = std::stoull(
                require_value(index, argc, argv, argument, name));
        }
        else if (argument == "-stride" || is_option(argument, "--downsampling-stride")) {
            const auto name = argument == "-stride" ? "-stride" : "--downsampling-stride";
            const auto value = std::stoull(require_value(index, argc, argv, argument, name, "-stride"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--downsampling-stride must be a positive uint32");
            options.calling.mutect2_downsampling_stride = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--max-suspicious-reads-per-alignment-start")) {
            const auto value = std::stoll(require_value(
                index, argc, argv, argument, "--max-suspicious-reads-per-alignment-start"));
            if (value < std::numeric_limits<std::int32_t>::min() ||
                value > std::numeric_limits<std::int32_t>::max())
                throw std::invalid_argument(
                    "--max-suspicious-reads-per-alignment-start must fit int32");
            options.calling.mutect2_max_suspicious_reads_per_alignment_start =
                static_cast<std::int32_t>(value);
        }
        else if (is_option(argument, "--downsampling-seed"))
            options.calling.downsampling_seed = std::stoull(
                require_value(index, argc, argv, argument, "--downsampling-seed"));
        else if (is_option(argument, "--allele-informative-reads-overlap-margin")) {
            const auto value = std::stoull(require_value(
                index, argc, argv, argument, "--allele-informative-reads-overlap-margin"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument(
                    "--allele-informative-reads-overlap-margin must fit uint32");
            options.calling.informative_read_overlap_margin =
                static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--kmer-size")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--kmer-size"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--kmer-size must be a positive uint32");
            options.calling.graph_kmer_size = static_cast<std::uint32_t>(value);
            if (!options.graph_kmer_sizes_explicit) {
                options.calling.graph_kmer_sizes.clear();
                options.graph_kmer_sizes_explicit = true;
            }
            options.calling.graph_kmer_sizes.push_back(static_cast<std::uint32_t>(value));
        }
        else if (argument == "--dont-increase-kmer-sizes-for-cycles" ||
                 argument.rfind("--dont-increase-kmer-sizes-for-cycles=", 0) == 0)
            options.calling.graph_dont_increase_kmer_sizes_for_cycles =
                fastgatk::native::parse_optional_boolean(
                    index, argc, argv, argument, "--dont-increase-kmer-sizes-for-cycles");
        else if (is_option(argument, "--min-kmer-count")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--min-kmer-count"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-kmer-count must be positive and fit uint32");
            options.calling.graph_min_kmer_count = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-pruning")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--min-pruning"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-pruning must fit uint32");
            options.calling.graph_min_pruning = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--num-pruning-samples")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--num-pruning-samples"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--num-pruning-samples must be positive and fit uint32");
            options.calling.graph_num_pruning_samples = static_cast<std::uint32_t>(value);
        }
        else if (argument == "--adaptive-pruning")
            options.calling.graph_use_adaptive_pruning = true;
        else if (argument == "--disable-adaptive-pruning" ||
                 argument.rfind("--disable-adaptive-pruning=", 0) == 0)
            options.calling.graph_use_adaptive_pruning = !fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--disable-adaptive-pruning");
        else if (argument == "--linked-de-bruijn-graph")
            options.calling.graph_linked_de_bruijn = true;
        else if (argument == "--disable-artificial-haplotype-recovery")
            options.calling.graph_disable_artificial_haplotype_recovery = true;
        else if (argument == "--enable-legacy-graph-cycle-detection")
            options.calling.graph_enable_legacy_cycle_detection = true;
        else if (is_option(argument, "--adaptive-pruning-initial-error-rate"))
            options.calling.graph_initial_error_rate_for_pruning = std::stod(
                require_value(index, argc, argv, argument, "--adaptive-pruning-initial-error-rate"));
        else if (is_option(argument, "--pruning-lod-threshold")) {
            options.calling.graph_pruning_log_odds_threshold = std::stod(
                require_value(index, argc, argv, argument, "--pruning-lod-threshold"));
            options.graph_pruning_log_odds_threshold_explicit = true;
        }
        else if (is_option(argument, "--pruning-seeding-lod-threshold")) {
            options.calling.graph_pruning_seeding_log_odds_threshold = std::stod(
                require_value(index, argc, argv, argument, "--pruning-seeding-lod-threshold"));
            options.graph_pruning_seeding_log_odds_threshold_explicit = true;
        }
        else if (is_option(argument, "--max-unpruned-variants")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--max-unpruned-variants"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-unpruned-variants must be positive and fit uint32");
            options.calling.graph_max_unpruned_variants = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-dangling-branch-length")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument,
                                                          "--min-dangling-branch-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-dangling-branch-length must fit uint32");
            options.calling.graph_min_dangling_branch_length = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-dangling-matching-bases")) {
            const auto value = std::stoll(require_value(index, argc, argv, argument,
                                                         "--min-dangling-matching-bases"));
            if (value < -1 || value > std::numeric_limits<std::int32_t>::max())
                throw std::invalid_argument("--min-dangling-matching-bases must be -1 or a non-negative int32");
            options.calling.graph_min_dangling_matching_bases = static_cast<std::int32_t>(value);
        }
        else if (argument == "--allow-non-unique-kmers-in-ref")
            options.calling.graph_allow_non_unique_kmers_in_ref = true;
        else if (argument == "--recover-all-dangling-branches" ||
                 argument.rfind("--recover-all-dangling-branches=", 0) == 0)
            options.calling.graph_recover_all_dangling_branches =
                fastgatk::native::parse_optional_boolean(
                    index, argc, argv, argument, "--recover-all-dangling-branches");
        else if (argument == "--error-correct-reads" ||
                 argument.rfind("--error-correct-reads=", 0) == 0)
            options.calling.error_correct_reads = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--error-correct-reads");
        else if (is_option(argument, "--error-correction-log-odds"))
            options.calling.pileup_error_correction_log_odds = std::stod(
                require_value(index, argc, argv, argument, "--error-correction-log-odds"));
        else if (is_option(argument, "--kmer-length-for-read-error-correction")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument,
                                                          "--kmer-length-for-read-error-correction"));
            if (value == 0 || value > 31)
                throw std::invalid_argument("--kmer-length-for-read-error-correction must be in [1,31]");
            options.calling.error_correction_kmer_length = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-observations-for-kmer-to-be-solid")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument,
                                                          "--min-observations-for-kmer-to-be-solid"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-observations-for-kmer-to-be-solid must be positive and fit uint32");
            options.calling.error_correction_min_solid_observations = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--max-haplotype-paths") ||
                 is_option(argument, "--max-num-haplotypes-in-population")) {
            const auto name = is_option(argument, "--max-num-haplotypes-in-population")
                ? "--max-num-haplotypes-in-population" : "--max-haplotype-paths";
            const auto value = std::stoull(require_value(index, argc, argv, argument, name));
            if (value == 0) throw std::invalid_argument("--max-haplotype-paths must be positive");
            options.calling.graph_max_paths = value;
        }
        else if (is_option(argument, "--max-haplotype-depth")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument, "--max-haplotype-depth"));
            if (value == 0) throw std::invalid_argument("--max-haplotype-depth must be positive");
            options.calling.graph_max_depth = value;
        }
        else if (is_option(argument, "--max-haplotype-combination-alleles")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument,
                                                          "--max-haplotype-combination-alleles"));
            if (value == 0 || value > 16)
                throw std::invalid_argument("--max-haplotype-combination-alleles must be in [1,16]");
            options.calling.max_haplotype_combination_alleles = value;
        }
        else if (is_option(argument, "--max-mnp-distance") || is_option(argument, "--mnp-dist")) {
            const auto value = std::stoull(require_value(index, argc, argv, argument,
                                                          argument.rfind("--mnp-dist", 0) == 0
                                                              ? "--mnp-dist" : "--max-mnp-distance"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-mnp-distance must fit uint32");
            options.calling.max_mnp_distance = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--read-filter")) {
            const auto filter_name = require_value(index, argc, argv, argument, "--read-filter");
            enable_read_filter(options.calling, filter_name);
            options.read_filters.push_back(filter_name);
        } else if (is_option(argument, "--disable-read-filter")) {
            const auto filter_name = require_value(index, argc, argv, argument,
                                                    "--disable-read-filter");
            options.disabled_read_filters.push_back(filter_name);
            disable_read_filter(options.calling, filter_name);
        } else if (argument == "--disable-tool-default-read-filters" ||
                   argument.rfind("--disable-tool-default-read-filters=", 0) == 0)
            options.calling.disable_tool_default_read_filters =
                fastgatk::native::parse_optional_boolean(
                    index, argc, argv, argument, "--disable-tool-default-read-filters");
        else if (argument == "--mitochondria-mode" ||
                 argument.rfind("--mitochondria-mode=", 0) == 0)
            options.mitochondria_mode = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--mitochondria-mode");
        else
            throw std::invalid_argument("unknown option: " + argument);
    }
    // M2ArgumentCollection selects the default AF from the presence of normal
    // samples: tumor-only uses 5e-8, tumor-normal uses 1e-6.  Apply this
    // after parsing sample selectors so an explicit value (or mitochondrial
    // mode) always wins.
    if (!options.mitochondria_mode && !options.population_allele_frequency_explicit &&
        (!options.normal.empty() || !options.normal_sample.empty() ||
         !options.normal_samples.empty()))
        options.population_allele_frequency = 1.0e-6;
    if (options.mitochondria_mode) {
        // M2ArgumentCollection.getMitochondriaModeNameValuePairs().  Explicit
        // caller values win; otherwise apply the mode's sensitivity defaults.
        if (!options.population_allele_frequency_explicit)
            options.population_allele_frequency = 4.0e-3;
        if (!options.initial_tumor_lod_explicit)
            options.initial_tumor_lod = 0.0;
        if (!options.tumor_lod_to_emit_explicit)
            options.tumor_lod_to_emit = 0.0;
        options.calling.graph_recover_all_dangling_branches = true;
        if (!options.graph_pruning_log_odds_threshold_explicit)
            options.calling.graph_pruning_log_odds_threshold = -4.0 * std::log(10.0);
        // M2ArgumentCollection.getMitochondriaModeNameValuePairs() overrides
        // only --pruning-lod-threshold.  The stricter default seeding gate
        // (ln(10^4)) is intentionally preserved: lowering it admits a
        // different bounded subgraph in AdaptiveChainPruner.
    }
    if (options.tumor.empty()) throw std::invalid_argument("--input is required");
    if (options.output.empty()) throw std::invalid_argument("--output is required");
    if (options.batch_records == 0 || options.threads < 1)
        throw std::invalid_argument("invalid batch/thread count");
    if (!std::isfinite(options.contamination) || options.contamination < 0.0 ||
        options.contamination >= 1.0 || !std::isfinite(options.somatic_prior) ||
        !std::isfinite(options.germline_prior) || !std::isfinite(options.artifact_prior) ||
        !std::isfinite(options.initial_tumor_lod) || options.initial_tumor_lod < 0.0 ||
        !std::isfinite(options.tumor_lod_to_emit) || options.tumor_lod_to_emit < 0.0 ||
        !std::isfinite(options.normal_log10_odds) ||
        !std::isfinite(options.minimum_allele_fraction) ||
        options.minimum_allele_fraction < 0.0 || options.minimum_allele_fraction >= 1.0 ||
        options.somatic_prior <= 0.0 || options.germline_prior <= 0.0 || options.artifact_prior <= 0.0)
        throw std::invalid_argument("somatic posterior contamination/priors are invalid");
    if (!std::isfinite(options.population_allele_frequency) ||
        options.population_allele_frequency <= 0.0 ||
        options.population_allele_frequency >= 1.0)
        throw std::invalid_argument("population allele frequency must be finite and in (0,1)");
    if (!std::isfinite(options.max_population_allele_frequency) ||
        options.max_population_allele_frequency < 0.0 ||
        options.max_population_allele_frequency > 1.0)
        throw std::invalid_argument("--max-population-af must be finite and in [0,1]");
    if (options.reference_confidence_mode == ReferenceConfidenceMode::Gvcf) {
        if (options.gvcf_lod_bands.empty())
            throw std::invalid_argument("--gvcf-lod-band requires at least one bound in GVCF mode");
        for (std::size_t index = 0; index < options.gvcf_lod_bands.size(); ++index) {
            const auto value = options.gvcf_lod_bands[index];
            if (!std::isfinite(value) ||
                (index != 0U && !(value > options.gvcf_lod_bands[index - 1U])))
                throw std::invalid_argument(
                    "--gvcf-lod-band values must be finite and strictly increasing");
        }
    }
    // Mutect2Engine's reference-confidence call delegates to
    // ReferenceConfidenceModel, whose 4.6.2 contract accepts exactly one
    // sample likelihood matrix.  A matched-normal invocation reaches that
    // call with tumor plus normal matrices and GATK rejects it rather than
    // emitting a mixed-sample gVCF.  Fail at the C++ Host option boundary as
    // well: the Kokkos PairHMM/somatic kernels must not fabricate output for
    // an input combination the source tool cannot represent.
    if (options.reference_confidence_mode != ReferenceConfidenceMode::None &&
        (!options.normal.empty() || !selected_normal_samples(options).empty()))
        throw std::invalid_argument(
            "BAD_INPUT: Mutect2 reference confidence requires exactly one sample; "
            "matched normals are unsupported by GATK 4.6.2");
    // Match HTSJDK's normalized repeatable argument list before the Kokkos
    // graph retry policy receives it.  graph_kmer_size intentionally keeps
    // the final explicit value for legacy telemetry consumers.
    std::sort(options.calling.graph_kmer_sizes.begin(),
              options.calling.graph_kmer_sizes.end());
    options.calling.graph_kmer_sizes.erase(
        std::unique(options.calling.graph_kmer_sizes.begin(),
                    options.calling.graph_kmer_sizes.end()),
        options.calling.graph_kmer_sizes.end());
    const bool emit_reference_confidence =
        options.reference_confidence_mode != ReferenceConfidenceMode::None;
    options.calling.emit_reference_confidence = emit_reference_confidence;
    options.calling.emit_reference_confidence_bp_resolution =
        options.reference_confidence_mode == ReferenceConfidenceMode::BpResolution;
    // M2ArgumentCollection uses the negative-infinity reference-confidence
    // emission/initial LOD so low-scoring EventMap alleles remain available
    // to the SomaticGVCF writer.  The Kokkos activity profile consumes that
    // same source threshold; normal VCF mode retains the user/default value.
    options.calling.activity_initial_tumor_log10_odds = emit_reference_confidence
        ? -std::numeric_limits<double>::infinity() : options.initial_tumor_lod;
    options.calling.force_active = options.force_active;
    // GATK's PairHMMNativeArgumentCollection exposes this Mutect2 option
    // with a default of false.  The calling pipeline also serves strict HC
    // tests and defaults to double there, so select the Mutect2 contract at
    // the C++ Host option boundary rather than changing the shared kernel.
    options.calling.pairhmm_use_double_precision =
        options.native_pair_hmm_use_double_precision;
    resolve_read_filter_options(options);
    return options;
}

std::vector<std::string> load_reference(const std::string& path,
                                        const fastgatk::io::HeaderSummary& header) {
    std::vector<std::string> sequences(header.contigs.size());
    if (path.empty()) return sequences;
#if FASTGATK_HAS_HTSLIB
    // GATK's ReferenceFileSource is random-access based and therefore reads
    // bgzip FASTA through its .fai/.gzi indexes.  An ifstream sees compressed
    // bytes as one non-FASTA line, which silently left every sequence empty
    // for production references such as hs37d5.fa.gz and made PairHMM take
    // the no-reference fallback.  Fetch only BAM-header contigs: this keeps
    // the host-side contiguous reference contract while avoiding a whole
    // genome materialization for a targeted BAM.
    // Avoid creating an FAI behind the caller's back for an unindexed plain
    // FASTA.  Indexed production references take this random-access path;
    // compact smoke fixtures continue through the non-mutating parser.
    std::error_code fai_error;
    const bool has_fai = std::filesystem::exists(path + ".fai", fai_error);
    if (fai_error)
        throw std::runtime_error("BACKEND_UNAVAILABLE: cannot inspect reference FAI: " + path + ".fai");
    if (has_fai) {
        faidx_t* fai = fai_load(path.c_str());
        if (fai == nullptr)
            throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA has an unreadable .fai: " + path);
        try {
            for (std::size_t index = 0; index < header.contigs.size(); ++index) {
                const auto& contig = header.contigs[index];
                const auto length = faidx_seq_len64(fai, contig.c_str());
                // Keep the indexed loader compatible with bounded traversals
                // that use a subset reference.  The walker consumes only the
                // contigs it visits; a missing unrelated header contig must
                // not turn that valid traversal into an eager-load failure.
                if (length < 0) continue;
                if (length == 0) continue;
                hts_pos_t fetched = 0;
                char* bases = faidx_fetch_seq64(fai, contig.c_str(), 0, length - 1, &fetched);
                if (bases == nullptr || fetched != length) {
                    free(bases);
                    throw std::runtime_error(
                        "BAD_INPUT: failed to fetch reference contig: " + contig);
                }
                sequences[index].assign(bases, static_cast<std::size_t>(fetched));
                free(bases);
            }
            fai_destroy(fai);
            return sequences;
        } catch (...) {
            fai_destroy(fai);
            throw;
        }
    }
#endif
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open reference FASTA: " + path);
    std::map<std::string, std::size_t> ids;
    for (std::size_t index = 0; index < header.contigs.size(); ++index)
        ids.emplace(header.contigs[index], index);
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

void append_batch(fastgatk::io::ReadBatch& all, const fastgatk::io::ReadBatch& batch) {
    if (!batch.indel_quality_layout_valid())
        throw std::runtime_error("BAD_INPUT: malformed BI/BD indel quality metadata in HTS batch");
    const auto base_offset = static_cast<std::uint32_t>(all.bases.size());
    const auto record_offset = all.records();
    const bool has_cigar = batch.cigar_offsets.size() == batch.records() + 1;
    const bool has_names = batch.name_offsets.size() == batch.records() + 1;
    const bool has_rg = batch.read_group_offsets.size() == batch.records() + 1;
    const bool has_oa = batch.original_alignment_offsets.size() == batch.records() + 1;
    const bool has_xm = batch.mate_contig_offsets.size() == batch.records() + 1;
    const bool has_oa_presence = batch.original_alignment_present.size() == batch.records();
    const bool has_xm_presence = batch.mate_contig_present.size() == batch.records();
    if ((!batch.original_alignment_present.empty() && !has_oa_presence) ||
        (!batch.mate_contig_present.empty() && !has_xm_presence))
        throw std::runtime_error("BAD_INPUT: malformed OA/XM presence metadata in HTS batch");
    const bool has_insertion_qualities =
        batch.insertion_quality_offsets.size() == batch.records() + 1;
    const bool has_deletion_qualities =
        batch.deletion_quality_offsets.size() == batch.records() + 1;
    const bool has_flow_tp = batch.flow_tp_offsets.size() == batch.records() + 1;
    const bool has_flow_t0 = batch.flow_t0_offsets.size() == batch.records() + 1;
    const bool has_flow_order = batch.flow_order_offsets.size() == batch.records() + 1;
    if (has_cigar && all.cigar_offsets.empty()) all.cigar_offsets.push_back(0);
    if (has_names && all.name_offsets.empty()) all.name_offsets.push_back(0);
    if (has_rg && all.read_group_offsets.empty()) all.read_group_offsets.push_back(0);
    // OA/XM are optional per record too.  Preserve absent tags as empty
    // spans/presence=false when merging independently produced BAMs, instead
    // of requiring every input header to expose the same optional-tag set.
    if (has_oa && all.original_alignment_offsets.empty())
        all.original_alignment_offsets.assign(record_offset + 1, 0);
    if (has_xm && all.mate_contig_offsets.empty())
        all.mate_contig_offsets.assign(record_offset + 1, 0);
    const auto materialize_presence = [record_offset](auto& presence,
                                                        const auto& offsets) {
        if (!presence.empty()) return;
        presence.reserve(record_offset);
        for (std::size_t record = 0; record < record_offset; ++record) {
            const bool present = offsets.size() == record_offset + 1 &&
                offsets[record + 1] > offsets[record];
            presence.push_back(present ? 1U : 0U);
        }
    };
    if (has_oa_presence) materialize_presence(
        all.original_alignment_present, all.original_alignment_offsets);
    if (has_xm_presence) materialize_presence(
        all.mate_contig_present, all.mate_contig_offsets);
    // BI/BD tags are optional per read, not per input.  A run can therefore
    // legitimately combine (for example) an untagged tumor BAM with a
    // tagged normal BAM.  Materialize empty spans for the earlier records so
    // the flat Host batch preserves HTSlib's per-record absence semantics;
    // PairHMM then applies GATK's default Q45 only to those empty spans.
    if (has_insertion_qualities && all.insertion_quality_offsets.empty())
        all.insertion_quality_offsets.assign(record_offset + 1, 0);
    if (has_deletion_qualities && all.deletion_quality_offsets.empty())
        all.deletion_quality_offsets.assign(record_offset + 1, 0);
    if (has_flow_tp && all.flow_tp_offsets.empty()) all.flow_tp_offsets.push_back(0);
    if (has_flow_t0 && all.flow_t0_offsets.empty()) all.flow_t0_offsets.push_back(0);
    if (has_flow_order && all.flow_order_offsets.empty()) all.flow_order_offsets.push_back(0);
    const bool output_has_insertion_qualities =
        !all.insertion_quality_offsets.empty();
    const bool output_has_deletion_qualities =
        !all.deletion_quality_offsets.empty();
    const bool output_has_oa = !all.original_alignment_offsets.empty();
    const bool output_has_xm = !all.mate_contig_offsets.empty();
    const bool output_has_oa_presence =
        has_oa_presence || !all.original_alignment_present.empty();
    const bool output_has_xm_presence =
        has_xm_presence || !all.mate_contig_present.empty();
    auto append_string = [](auto& payload, auto& offsets, const auto& source,
                            const auto& source_offsets, std::size_t record) {
        const auto begin = source_offsets[record];
        const auto end = source_offsets[record + 1];
        payload.insert(payload.end(), source.begin() + begin, source.begin() + end);
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    };
    for (std::size_t record = 0; record < batch.records(); ++record) {
        const auto begin = batch.offsets[record];
        const auto end = batch.offsets[record + 1];
        all.bases.insert(all.bases.end(), batch.bases.begin() + begin, batch.bases.begin() + end);
        all.qualities.insert(all.qualities.end(), batch.qualities.begin() + begin, batch.qualities.begin() + end);
        all.offsets.push_back(base_offset + static_cast<std::uint32_t>(all.bases.size() - base_offset));
        all.positions.push_back(batch.positions[record]);
        all.tids.push_back(batch.tids[record]);
        all.mapq.push_back(batch.mapq[record]);
        const auto output_record = record_offset + record;
        if (output_record > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("RESOURCE_EXHAUSTED: calling record ordinals exceed uint32 range");
        // This creates a new outer input collection, so its source ordinal is
        // the destination record rather than an ordinal inherited from a
        // temporary decode/sort batch.
        all.source_records.push_back(static_cast<std::uint32_t>(output_record));
        all.flags.push_back(batch.flags.empty() ? 0 : batch.flags[record]);
        all.mate_tids.push_back(batch.mate_tids.size() == batch.records()
                                    ? batch.mate_tids[record] : -1);
        all.mate_positions.push_back(batch.mate_positions.size() == batch.records()
                                         ? batch.mate_positions[record] : -1);
        all.template_lengths.push_back(batch.template_lengths.size() == batch.records()
                                           ? batch.template_lengths[record] : 0);
        if (has_cigar) {
            const auto cigar_begin = batch.cigar_offsets[record];
            const auto cigar_end = batch.cigar_offsets[record + 1];
            all.cigar_ops.insert(all.cigar_ops.end(), batch.cigar_ops.begin() + cigar_begin,
                                 batch.cigar_ops.begin() + cigar_end);
            all.cigar_offsets.push_back(static_cast<std::uint32_t>(all.cigar_ops.size()));
        }
        if (has_names) append_string(all.names, all.name_offsets, batch.names, batch.name_offsets, record);
        if (has_rg) append_string(all.read_groups, all.read_group_offsets,
                                  batch.read_groups, batch.read_group_offsets, record);
        if (output_has_oa) {
            if (has_oa) append_string(all.original_alignments, all.original_alignment_offsets,
                                      batch.original_alignments, batch.original_alignment_offsets, record);
            else
                all.original_alignment_offsets.push_back(
                    static_cast<std::uint32_t>(all.original_alignments.size()));
        }
        if (output_has_xm) {
            if (has_xm) append_string(all.mate_contigs, all.mate_contig_offsets,
                                      batch.mate_contigs, batch.mate_contig_offsets, record);
            else
                all.mate_contig_offsets.push_back(
                    static_cast<std::uint32_t>(all.mate_contigs.size()));
        }
        if (output_has_oa_presence) {
            const auto present = has_oa_presence
                ? batch.original_alignment_present[record]
                : (has_oa && batch.original_alignment_offsets[record + 1] >
                    batch.original_alignment_offsets[record]);
            all.original_alignment_present.push_back(present ? 1U : 0U);
        }
        if (output_has_xm_presence) {
            const auto present = has_xm_presence
                ? batch.mate_contig_present[record]
                : (has_xm && batch.mate_contig_offsets[record + 1] >
                    batch.mate_contig_offsets[record]);
            all.mate_contig_present.push_back(present ? 1U : 0U);
        }
        if (output_has_insertion_qualities) {
            if (has_insertion_qualities)
                append_string(all.insertion_qualities, all.insertion_quality_offsets,
                              batch.insertion_qualities, batch.insertion_quality_offsets, record);
            else
                all.insertion_quality_offsets.push_back(
                    static_cast<std::uint32_t>(all.insertion_qualities.size()));
        }
        if (output_has_deletion_qualities) {
            if (has_deletion_qualities)
                append_string(all.deletion_qualities, all.deletion_quality_offsets,
                              batch.deletion_qualities, batch.deletion_quality_offsets, record);
            else
                all.deletion_quality_offsets.push_back(
                    static_cast<std::uint32_t>(all.deletion_qualities.size()));
        }
        if (has_flow_tp) {
            const auto flow_begin = batch.flow_tp_offsets[record];
            const auto flow_end = batch.flow_tp_offsets[record + 1];
            all.flow_tp.insert(all.flow_tp.end(), batch.flow_tp.begin() + flow_begin,
                               batch.flow_tp.begin() + flow_end);
            all.flow_tp_offsets.push_back(static_cast<std::uint32_t>(all.flow_tp.size()));
        }
        if (has_flow_t0) {
            const auto flow_begin = batch.flow_t0_offsets[record];
            const auto flow_end = batch.flow_t0_offsets[record + 1];
            all.flow_t0_phred.insert(all.flow_t0_phred.end(), batch.flow_t0_phred.begin() + flow_begin,
                                     batch.flow_t0_phred.begin() + flow_end);
            all.flow_t0_offsets.push_back(static_cast<std::uint32_t>(all.flow_t0_phred.size()));
        }
        if (has_flow_order) {
            const auto flow_begin = batch.flow_order_offsets[record];
            const auto flow_end = batch.flow_order_offsets[record + 1];
            all.flow_orders.insert(all.flow_orders.end(), batch.flow_orders.begin() + flow_begin,
                                   batch.flow_orders.begin() + flow_end);
            all.flow_order_offsets.push_back(static_cast<std::uint32_t>(all.flow_orders.size()));
        }
        all.flow_max_hmer.push_back(batch.flow_max_hmer.size() == batch.records()
                                        ? batch.flow_max_hmer[record] : 0);
        if (batch.sample_ids.size() == batch.records())
            all.sample_ids.push_back(batch.sample_ids[record]);
    }
}

// HtsReader deliberately exposes source-local @RG/SM names rather than a
// process-global ordinal.  Once Mutect2 has resolved its sample roles, assign
// a stable Host identity before batches are merged: pair cleanup, assembly
// fragment de-duplication, and likelihood grouping must never coalesce reads
// merely because two tumor samples reuse a read name or read-group ID.
void label_batch_sample(fastgatk::io::ReadBatch& batch, const std::int32_t sample_id) {
    batch.sample_ids.assign(batch.records(), sample_id);
}

// Materialize one HTS record without retaining the source reader.  The
// region-streaming path uses this to perform a bounded coordinate merge across
// repeatable -I shards.  Keeping the copy at the ReadBatch boundary preserves
// CIGAR, BI/BD, flow tags and mate metadata exactly as in the ordinary
// aggregate path; no bam1_t pointer escapes the Host reader.
void append_record(fastgatk::io::ReadBatch& destination,
                   const fastgatk::io::ReadBatch& source,
                   const std::size_t index) {
    const auto count = source.records();
    if (index >= count || source.offsets.size() != count + 1)
        throw std::runtime_error("BAD_INPUT: malformed record offsets in HTS batch");
    const auto begin = source.offsets[index];
    const auto end = source.offsets[index + 1];
    if (begin > end || end > source.bases.size() || end > source.qualities.size())
        throw std::runtime_error("BAD_INPUT: malformed base/quality span in HTS batch");
    fastgatk::io::ReadBatch one;
    one.offsets = {0, static_cast<std::uint32_t>(end - begin)};
    one.bases.assign(source.bases.begin() + begin, source.bases.begin() + end);
    one.qualities.assign(source.qualities.begin() + begin, source.qualities.begin() + end);
    one.positions = {source.positions.size() == count ? source.positions[index] : -1};
    one.tids = {source.tids.size() == count ? source.tids[index] : -1};
    one.mapq = {static_cast<std::uint8_t>(source.mapq.size() == count ? source.mapq[index] : 0)};

    const auto copy_span = [&](const auto& source_offsets, const auto& source_values,
                               auto& destination_offsets, auto& destination_values) {
        if (source_offsets.size() != count + 1) return;
        const auto span_begin = source_offsets[index];
        const auto span_end = source_offsets[index + 1];
        if (span_begin > span_end || span_end > source_values.size())
            throw std::runtime_error("BAD_INPUT: malformed per-record metadata span");
        destination_offsets = {0, static_cast<std::uint32_t>(span_end - span_begin)};
        destination_values.assign(source_values.begin() + span_begin,
                                  source_values.begin() + span_end);
    };
    copy_span(source.cigar_offsets, source.cigar_ops, one.cigar_offsets, one.cigar_ops);
    copy_span(source.name_offsets, source.names, one.name_offsets, one.names);
    copy_span(source.read_group_offsets, source.read_groups,
              one.read_group_offsets, one.read_groups);
    copy_span(source.original_alignment_offsets, source.original_alignments,
              one.original_alignment_offsets, one.original_alignments);
    if (source.original_alignment_present.size() == count)
        one.original_alignment_present = {source.original_alignment_present[index]};
    copy_span(source.mate_contig_offsets, source.mate_contigs,
              one.mate_contig_offsets, one.mate_contigs);
    if (source.mate_contig_present.size() == count)
        one.mate_contig_present = {source.mate_contig_present[index]};
    copy_span(source.insertion_quality_offsets, source.insertion_qualities,
              one.insertion_quality_offsets, one.insertion_qualities);
    copy_span(source.deletion_quality_offsets, source.deletion_qualities,
              one.deletion_quality_offsets, one.deletion_qualities);
    copy_span(source.flow_tp_offsets, source.flow_tp, one.flow_tp_offsets, one.flow_tp);
    copy_span(source.flow_t0_offsets, source.flow_t0_phred,
              one.flow_t0_offsets, one.flow_t0_phred);
    copy_span(source.flow_order_offsets, source.flow_orders,
              one.flow_order_offsets, one.flow_orders);
    if (source.flow_max_hmer.size() == count) one.flow_max_hmer = {source.flow_max_hmer[index]};
    if (source.flags.size() == count) one.flags = {source.flags[index]};
    if (source.mate_tids.size() == count) one.mate_tids = {source.mate_tids[index]};
    if (source.mate_positions.size() == count) one.mate_positions = {source.mate_positions[index]};
    if (source.template_lengths.size() == count) one.template_lengths = {source.template_lengths[index]};
    if (source.sample_ids.size() == count) one.sample_ids = {source.sample_ids[index]};
    append_batch(destination, one);
}

// Match Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(): unless
// --ignore-itr-artifacts is set, GATK applies
// PalindromeArtifactClipReadTransformer with a five-base minimum palindrome.
// This stays at the C++ Host boundary because it changes a decoded GATKRead;
// the compact ReadBatch then proceeds unchanged into the Kokkos assembly and
// PairHMM kernels.
constexpr std::size_t kMutect2MinPalindromeSize = 5;
constexpr double kMutect2MinItrMatchingFraction = 0.9;
constexpr std::uint16_t kSamFlagPaired = 0x1U;
constexpr std::uint16_t kSamFlagProperPair = 0x2U;
constexpr std::uint16_t kSamFlagUnmapped = 0x4U;
constexpr std::uint16_t kSamFlagMateUnmapped = 0x8U;
constexpr std::uint16_t kSamFlagReverse = 0x10U;
constexpr std::uint16_t kSamFlagMateReverse = 0x20U;

std::uint8_t gatk_itr_complement(const std::uint8_t base) {
    // BaseUtils.getComplement() is deliberately stricter than a generic IUPAC
    // complement: the Mutect2 transformer accepts A/C/G/T/N reference bases
    // and throws for any other reference symbol.
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        case 'N': case 'n': return 'N';
        default:
            throw std::runtime_error(
                "BAD_INPUT: ITR palindrome clipping encountered an unsupported reference base");
    }
}

bool itr_hard_clip_terminal(fastgatk::io::ReadBatch& read, const bool leading,
                            const std::size_t clip_count) {
    if (read.records() != 1 || read.offsets.size() != 2 ||
        read.offsets.front() > read.offsets.back() ||
        read.offsets.back() > read.bases.size() || read.offsets.back() > read.qualities.size())
        throw std::runtime_error("BAD_INPUT: malformed ITR clipping read batch");
    const auto read_length = static_cast<std::size_t>(read.offsets.back());
    if (clip_count == 0 || clip_count > read_length)
        throw std::runtime_error("BAD_INPUT: ITR clipping span is outside the read");

    // ClippingOp.applyHardClipBases() returns ReadUtils.emptyRead() for a
    // whole-read clip.  Keep the same safe, unmapped Host representation so
    // downstream read predicates reject it naturally.
    if (clip_count == read_length) {
        read.bases.clear();
        read.qualities.clear();
        read.offsets = {0, 0};
        read.cigar_ops.clear();
        read.cigar_offsets = {0, 0};
        read.insertion_qualities.clear();
        read.deletion_qualities.clear();
        if (!read.insertion_quality_offsets.empty()) read.insertion_quality_offsets = {0, 0};
        if (!read.deletion_quality_offsets.empty()) read.deletion_quality_offsets = {0, 0};
        if (read.flags.size() == 1) read.flags.front() |= kSamFlagUnmapped;
        if (read.mapq.size() == 1) read.mapq.front() = 0;
        return true;
    }

    if (!read.cigar_record_layout_valid(0))
        throw std::runtime_error("BAD_INPUT: malformed CIGAR during ITR clipping");
    const auto cigar_begin = static_cast<std::size_t>(read.cigar_offsets.front());
    const auto cigar_end = static_cast<std::size_t>(read.cigar_offsets.back());
    if (cigar_begin == cigar_end)
        throw std::runtime_error("BAD_INPUT: empty CIGAR during ITR clipping");
    const auto terminal_index = leading ? cigar_begin : cigar_end - 1;
    const auto terminal = fastgatk::io::CigarOp::unpack(read.cigar_ops[terminal_index]);
    if (!terminal.valid() || static_cast<std::size_t>(terminal.length) != clip_count ||
        (terminal.code != fastgatk::io::CigarOpCode::SoftClip &&
         terminal.code != fastgatk::io::CigarOpCode::Insertion))
        throw std::runtime_error("BAD_INPUT: ITR clip no longer matches terminal S/I CIGAR");

    const auto new_length = read_length - clip_count;
    const auto trim_indel_qualities = [&](auto& offsets, auto& qualities) {
        if (offsets.empty()) {
            offsets = {0, static_cast<std::uint32_t>(new_length)};
            qualities.assign(new_length, 45U);
            return;
        }
        if (offsets.size() != 2 || offsets.front() > offsets.back() ||
            offsets.back() > qualities.size())
            throw std::runtime_error("BAD_INPUT: malformed BI/BD metadata during ITR clipping");
        const auto span_length = static_cast<std::size_t>(offsets.back() - offsets.front());
        if (span_length == 0) {
            // ReadClipper materializes the missing companion BI/BD tag as
            // GATK's Q45 default whenever either indel-quality tag exists.
            qualities.assign(new_length, 45U);
        } else {
            if (span_length != read_length)
                throw std::runtime_error("BAD_INPUT: BI/BD length differs from read during ITR clipping");
            if (leading)
                qualities.erase(qualities.begin(), qualities.begin() + clip_count);
            else
                qualities.erase(qualities.begin() + new_length, qualities.end());
        }
        offsets = {0, static_cast<std::uint32_t>(qualities.size())};
    };
    const bool has_insertion_qualities = read.insertion_quality_offsets.size() == 2 &&
        read.insertion_quality_offsets.back() > read.insertion_quality_offsets.front();
    const bool has_deletion_qualities = read.deletion_quality_offsets.size() == 2 &&
        read.deletion_quality_offsets.back() > read.deletion_quality_offsets.front();

    if (leading) {
        read.bases.erase(read.bases.begin(), read.bases.begin() + clip_count);
        read.qualities.erase(read.qualities.begin(), read.qualities.begin() + clip_count);
    } else {
        read.bases.erase(read.bases.begin() + new_length, read.bases.end());
        read.qualities.erase(read.qualities.begin() + new_length, read.qualities.end());
    }
    read.offsets = {0, static_cast<std::uint32_t>(new_length)};
    read.cigar_ops[terminal_index] = fastgatk::io::CigarOp{
        terminal.length, fastgatk::io::CigarOpCode::HardClip}.pack();
    if (has_insertion_qualities || has_deletion_qualities) {
        trim_indel_qualities(read.insertion_quality_offsets, read.insertion_qualities);
        trim_indel_qualities(read.deletion_quality_offsets, read.deletion_qualities);
    }
    return true;
}

bool clip_mutect2_itr_artifact_record(
    fastgatk::io::ReadBatch& read, const std::vector<std::string>& references) {
    if (read.records() != 1 || read.offsets.size() != 2 || read.flags.size() != 1 ||
        read.positions.size() != 1 || read.tids.size() != 1 || read.mate_positions.size() != 1 ||
        read.template_lengths.size() != 1 || !read.cigar_record_layout_valid(0))
        return false;
    const auto flags = read.flags.front();
    if ((flags & kSamFlagProperPair) == 0 || (flags & kSamFlagPaired) == 0 ||
        (flags & (kSamFlagUnmapped | kSamFlagMateUnmapped)) != 0 ||
        ((flags & kSamFlagReverse) != 0) == ((flags & kSamFlagMateReverse) != 0))
        return false;
    const auto tid = read.tids.front();
    const auto position = read.positions.front();
    const auto mate_position = read.mate_positions.front();
    const auto template_length = read.template_lengths.front();
    if (tid < 0 || static_cast<std::size_t>(tid) >= references.size() || position < 0 ||
        mate_position < 0 || template_length == 0)
        return false;
    const auto read_length = static_cast<std::size_t>(read.offsets.back() - read.offsets.front());
    const auto cigar_begin = static_cast<std::size_t>(read.cigar_offsets.front());
    const auto cigar_end = static_cast<std::size_t>(read.cigar_offsets.back());
    if (cigar_begin == cigar_end || cigar_end > read.cigar_ops.size()) return false;
    std::int64_t reference_span = 0;
    for (std::size_t cigar = cigar_begin; cigar < cigar_end; ++cigar) {
        const auto operation = fastgatk::io::CigarOp::unpack(read.cigar_ops[cigar]);
        if (!operation.valid()) return false;
        if (operation.consumes_reference()) reference_span += operation.length;
    }
    const auto read_start = static_cast<std::int64_t>(position) + 1;
    const auto read_end = static_cast<std::int64_t>(position) + reference_span;
    const auto mate_start = static_cast<std::int64_t>(mate_position) + 1;
    const bool reverse = (flags & kSamFlagReverse) != 0;
    // ReadUtils.hasWellDefinedFragmentSize().  Keep this before computing an
    // adaptor boundary, including its intentionally permissive handling of a
    // malformed cross-contig nonzero TLEN just as the Java implementation.
    if ((reverse && read_end <= mate_start) ||
        (!reverse && read_start > mate_start + static_cast<std::int64_t>(template_length)))
        return false;
    const bool upstream_of_mate = template_length > 0;
    const auto terminal_index = upstream_of_mate ? cigar_begin : cigar_end - 1;
    const auto terminal = fastgatk::io::CigarOp::unpack(read.cigar_ops[terminal_index]);
    if (!terminal.valid() || (terminal.code != fastgatk::io::CigarOpCode::SoftClip &&
                              terminal.code != fastgatk::io::CigarOpCode::Insertion))
        return false;
    const auto potential_artifact_bases = static_cast<std::size_t>(terminal.length);
    const auto bases_to_compare = std::min(
        potential_artifact_bases + kMutect2MinPalindromeSize, read_length);
    if (bases_to_compare == 0) return false;
    const auto absolute_template_length = template_length < 0
        ? -static_cast<std::int64_t>(template_length)
        : static_cast<std::int64_t>(template_length);
    const auto adaptor_boundary = reverse ? mate_start - 1 : read_start + absolute_template_length;
    const auto ref_start = upstream_of_mate
        ? adaptor_boundary - static_cast<std::int64_t>(bases_to_compare)
        : adaptor_boundary + 1;
    const auto ref_end = upstream_of_mate ? adaptor_boundary - 1
                                           : adaptor_boundary + static_cast<std::int64_t>(bases_to_compare);
    const auto& reference = references[static_cast<std::size_t>(tid)];
    if (ref_start < 1 || ref_end > static_cast<std::int64_t>(reference.size()) ||
        (upstream_of_mate && ref_start < read_start) ||
        (!upstream_of_mate && read_end < ref_end))
        return false;
    std::size_t matching_bases = 0;
    auto read_index = upstream_of_mate ? bases_to_compare - 1 : read_length - 1;
    for (std::size_t offset = 0; offset < bases_to_compare; ++offset) {
        const auto reference_index = static_cast<std::size_t>(ref_start - 1) + offset;
        if (gatk_itr_complement(static_cast<std::uint8_t>(reference[reference_index])) ==
            read.bases[read_index])
            ++matching_bases;
        --read_index;
    }
    if (static_cast<double>(matching_bases) / static_cast<double>(bases_to_compare) <
        kMutect2MinItrMatchingFraction)
        return false;
    return itr_hard_clip_terminal(read, upstream_of_mate, potential_artifact_bases);
}

std::uint64_t apply_mutect2_itr_artifact_clipping(
    fastgatk::io::ReadBatch& reads, const std::vector<std::string>& references) {
    if (reads.records() == 0) return 0;
    fastgatk::io::ReadBatch transformed;
    transformed.offsets.push_back(0);
    std::uint64_t clipped = 0;
    for (std::size_t record = 0; record < reads.records(); ++record) {
        fastgatk::io::ReadBatch one;
        one.offsets.push_back(0);
        append_record(one, reads, record);
        if (clip_mutect2_itr_artifact_record(one, references)) {
            ++clipped;
        }
        append_batch(transformed, one);
    }
    // The transformer changes a read in place in GATK; it must not replace
    // the outer input identity assembled by the Kokkos read-filter selection.
    if (reads.source_records.size() == reads.records())
        transformed.source_records = reads.source_records;
    if (reads.source_input_ids.size() == reads.records())
        transformed.source_input_ids = reads.source_input_ids;
    reads = std::move(transformed);
    return clipped;
}

void configure_mutect2_itr_artifact_transform(
    Options& options, const std::vector<std::string>& references) {
    options.calling.post_read_filter_transform = {};
    if (options.ignore_itr_artifacts) return;
    options.calling.post_read_filter_transform = [&references](fastgatk::io::ReadBatch& reads) {
        return static_cast<std::size_t>(apply_mutect2_itr_artifact_clipping(reads, references));
    };
}

// HTSlib's indexed iterators are coordinate ordered per file.  GATK merges
// multiple ReadsDataSource iterators with the same comparator; the Host
// staging path materializes one bounded tile and performs the equivalent
// stable merge, using input order as the deterministic tie-break.
void sort_batch_by_coordinate(fastgatk::io::ReadBatch& batch) {
    const auto count = batch.records();
    if (count < 2) return;
    std::vector<std::size_t> order(count);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](const auto lhs, const auto rhs) {
        const auto lhs_tid = lhs < batch.tids.size() ? batch.tids[lhs] : -1;
        const auto rhs_tid = rhs < batch.tids.size() ? batch.tids[rhs] : -1;
        const auto lhs_pos = lhs < batch.positions.size() ? batch.positions[lhs] : -1;
        const auto rhs_pos = rhs < batch.positions.size() ? batch.positions[rhs] : -1;
        return std::tie(lhs_tid, lhs_pos) < std::tie(rhs_tid, rhs_pos);
    });
    fastgatk::io::ReadBatch sorted;
    sorted.offsets.push_back(0);
    for (const auto index : order) append_record(sorted, batch, index);
    batch = std::move(sorted);
}

// Mutect2's matched-normal likelihood pass does not discover new active
// regions.  Mutect2Engine has already made that Host-side traversal on the
// joint tumor/normal pileup, and each normal read is subsequently scored only
// in the AssemblyRegion that owns its shared graph and EventMap.  Keep just
// those read windows before invoking the sample-local run: otherwise its
// preliminary Host activity/assembly preparation walks an entire normal BAM
// only to replace the result with the supplied AssemblyRegions below.
fastgatk::io::ReadBatch select_reads_for_assembly_regions(
    const fastgatk::io::ReadBatch& reads,
    const std::vector<fastgatk::calling::AssemblyRegionAssembly>& assemblies) {
    fastgatk::io::ReadBatch selected;
    selected.offsets.push_back(0);
    for (std::size_t record = 0; record < reads.records(); ++record) {
        if (record >= reads.tids.size() || record >= reads.positions.size() ||
            reads.positions[record] < 0)
            continue;
        auto reference_end = fastgatk::io::reference_end(reads, record);
        if (reference_end < 0) {
            if (record + 1 >= reads.offsets.size()) continue;
            reference_end = static_cast<std::int64_t>(reads.positions[record]) +
                static_cast<std::int64_t>(reads.offsets[record + 1] - reads.offsets[record]);
        }
        const auto reference_start = static_cast<std::int64_t>(reads.positions[record]);
        const auto overlaps = std::any_of(assemblies.begin(), assemblies.end(),
            [&](const auto& assembly) {
                const auto& region = assembly.region;
                return region.tid == reads.tids[record] &&
                    reference_start <= static_cast<std::int64_t>(region.end) &&
                    reference_end > static_cast<std::int64_t>(region.start);
            });
        if (overlaps) append_record(selected, reads, record);
    }
    return selected;
}

fastgatk::io::ReadBatch read_all(const std::string& path, const std::string& reference,
                                 const std::vector<std::string>& regions, std::size_t batch_records,
                                 fastgatk::io::HeaderSummary& header,
                                 std::uint64_t safe_memory_budget,
                                 const std::string& selected_sample = {},
                                 bool default_to_first_sample = false,
                                 std::string* effective_sample = nullptr,
                                 std::size_t* interval_file_inputs = nullptr,
                                 std::size_t* interval_file_records = nullptr,
                                 std::uint64_t* adaptive_batch_reductions = nullptr,
                                 std::size_t* effective_batch_records = nullptr,
                                 std::vector<fastgatk::io::HtsInterval>* requested_intervals = nullptr,
                                 bool allow_empty = false,
                                 fastgatk::io::HtsIntervalSetRule interval_set_rule =
                                     fastgatk::io::HtsIntervalSetRule::Union) {
    fastgatk::runtime::ResourceBudget budget;
    budget.host_target_bytes = safe_memory_budget;
    const auto controller = fastgatk::runtime::AdaptiveController{};
    const auto initial_limits = controller.initial(
        budget, fastgatk::runtime::WorkEstimate{
            static_cast<std::uint32_t>(std::min<std::size_t>(
                batch_records, std::numeric_limits<std::uint32_t>::max())),
            1, 2048, 1024, 2048});
    const auto initial_batch_records = std::min<std::size_t>(batch_records, initial_limits.max_reads);
    // Sample selection belongs at the HTSlib Host boundary.  Filtering after
    // decode would still mix evidence from unrelated @RG/@SM records into
    // Mutect2's fragment grouping and PairHMM matrices.  Mutect2's tumor
    // default is the first header sample, so resolve that before the real
    // reader is opened.  GATK rejects a tumor input with no @RG/SM sample
    // (ReferenceConfidenceModel requires a non-empty sample set); never fall
    // back to a synthetic TUMOR name after decoding all reads.
    std::string resolved_sample = selected_sample;
    if (resolved_sample.empty() && default_to_first_sample) {
        fastgatk::io::HtsReader probe(path, reference, regions, initial_batch_records,
                                      {}, {}, 0, 0, interval_set_rule);
        if (probe.header().samples.empty())
            throw std::runtime_error(
                "BAD_INPUT: tumor input header has no @RG SM sample (GATK requires a sample)");
        resolved_sample = probe.header().samples.front();
    }
    fastgatk::io::HtsReader reader(path, reference, regions, initial_batch_records,
                                   resolved_sample, {}, 0, 0, interval_set_rule);
    if (effective_sample != nullptr) *effective_sample = resolved_sample;
    auto batch_limits = initial_limits;
    batch_limits.max_reads = static_cast<std::uint32_t>(std::min<std::size_t>(
        reader.batch_records(), std::numeric_limits<std::uint32_t>::max()));
    if (effective_batch_records != nullptr) *effective_batch_records = reader.batch_records();
    if (interval_file_inputs != nullptr) *interval_file_inputs = reader.interval_file_inputs();
    if (interval_file_records != nullptr) *interval_file_records = reader.interval_file_records();
    if (requested_intervals != nullptr) *requested_intervals = reader.requested_intervals();
    header = reader.header();
    fastgatk::io::ReadBatch all, batch;
    all.offsets.push_back(0);
    while (reader.next(batch)) {
        const auto batch_bytes = batch.bytes();
        const auto accumulated_bytes = all.bytes();
        const auto estimated_bytes = accumulated_bytes >
                std::numeric_limits<std::uint64_t>::max() - batch_bytes
            ? std::numeric_limits<std::uint64_t>::max()
            : accumulated_bytes + batch_bytes;
        if (safe_memory_budget != 0 && estimated_bytes > safe_memory_budget)
            throw std::runtime_error("RESOURCE_EXHAUSTED: input staging exceeds safe memory budget");
        append_batch(all, batch);
        fastgatk::runtime::RuntimeTelemetry runtime;
        runtime.host_bytes = batch_bytes;
        runtime.inflight_bytes = batch_bytes;
        runtime.compute_queue_empty = true;
        const auto pressure = safe_memory_budget != 0 && batch_bytes > safe_memory_budget
            ? fastgatk::runtime::Pressure::HostMemory
            : fastgatk::runtime::Pressure::Normal;
        const auto next_limits = controller.next(runtime, batch_limits, pressure);
        if (next_limits.max_reads < reader.batch_records()) {
            reader.set_batch_records(next_limits.max_reads);
            if (adaptive_batch_reductions != nullptr) ++*adaptive_batch_reductions;
        }
        batch_limits = next_limits;
        if (effective_batch_records != nullptr) *effective_batch_records = reader.batch_records();
    }
    if (all.records() == 0 && !allow_empty)
        throw std::runtime_error("BAD_INPUT: no reads selected");
    return all;
}

void merge_input_headers(fastgatk::io::HeaderSummary& destination,
                         const fastgatk::io::HeaderSummary& source) {
    if (destination.contigs.empty()) {
        destination = source;
        return;
    }
    if (destination.contigs != source.contigs ||
        destination.contig_lengths != source.contig_lengths)
        throw std::runtime_error(
            "BAD_INPUT: repeated -I inputs have incompatible sequence dictionaries");
    if (destination.contig_assemblies.size() != destination.contigs.size())
        destination.contig_assemblies.resize(destination.contigs.size());
    for (std::size_t index = 0; index < source.contig_assemblies.size() &&
                                  index < destination.contig_assemblies.size(); ++index) {
        if (destination.contig_assemblies[index].empty())
            destination.contig_assemblies[index] = source.contig_assemblies[index];
        else if (!source.contig_assemblies[index].empty() &&
                 destination.contig_assemblies[index] != source.contig_assemblies[index])
            throw std::runtime_error(
                "BAD_INPUT: repeated -I inputs have conflicting @SQ AS assembly tags");
    }
    for (const auto& sample : source.samples) {
        if (std::find(destination.samples.begin(), destination.samples.end(), sample) ==
            destination.samples.end())
            destination.samples.push_back(sample);
    }
    for (const auto& [read_group, sample] : source.read_group_samples) {
        const auto existing = destination.read_group_samples.find(read_group);
        if (existing == destination.read_group_samples.end()) {
            destination.read_group_samples[read_group] = sample;
            continue;
        }
        if (existing->second == sample) continue;
        // @RG IDs are scoped to the header of their originating read input.
        // GATK's ReadsPathDataSource accepts multiple -I files that reuse an
        // ID for different SM values and resolves each record against that
        // input's header.  read_all_many() does the same before batches are
        // appended, so rejecting this otherwise valid layout here made the
        // merged *metadata* view stricter than GATK.  HeaderSummary has no
        // per-input RG namespace because Mutect2 never uses this merged map
        // to classify a decoded record; retain the first mapping solely for
        // legacy metadata consumers and keep source-local maps authoritative.
    }
}

fastgatk::io::ReadBatch read_all_many(
    const std::vector<std::string>& paths, const std::string& reference,
    const std::vector<std::string>& regions, std::size_t batch_records,
    fastgatk::io::HeaderSummary& header, std::uint64_t safe_memory_budget,
    const std::string& selected_sample = {}, bool default_to_first_sample = false,
    std::string* effective_sample = nullptr,
    std::size_t* interval_file_inputs = nullptr,
    std::size_t* interval_file_records = nullptr,
    std::uint64_t* adaptive_batch_reductions = nullptr,
    std::size_t* effective_batch_records = nullptr,
    std::vector<fastgatk::io::HtsInterval>* requested_intervals = nullptr,
    bool allow_empty = false,
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union,
    // Keep the public `regions` selectors as the source of interval metadata
    // and output-domain ownership, while allowing AssemblyRegion callers to
    // fetch a wider, Host-only read window.  GATK does this before graph
    // construction: `-L` is the core call domain, not a request to clip every
    // read used to assemble an event at its boundary.
    const std::vector<std::string>* read_regions = nullptr) {
    if (paths.empty()) throw std::invalid_argument("BAD_INPUT: at least one -I input is required");
    if (batch_records == 0) throw std::invalid_argument("HTS batch size must be positive");

    fastgatk::io::HeaderSummary merged_header;
    std::vector<fastgatk::io::HeaderSummary> input_headers;
    input_headers.reserve(paths.size());
    std::vector<std::vector<fastgatk::io::HtsInterval>> input_intervals;
    input_intervals.reserve(paths.size());
    std::vector<std::size_t> input_interval_file_inputs;
    std::vector<std::size_t> input_interval_file_records;
    input_interval_file_inputs.reserve(paths.size());
    input_interval_file_records.reserve(paths.size());
    for (const auto& path : paths) {
        fastgatk::io::HtsReader probe(path, reference, regions, batch_records,
                                      {}, {}, 0, 0, interval_set_rule);
        input_headers.push_back(probe.header());
        input_intervals.push_back(probe.requested_intervals());
        input_interval_file_inputs.push_back(probe.interval_file_inputs());
        input_interval_file_records.push_back(probe.interval_file_records());
        merge_input_headers(merged_header, input_headers.back());
    }

    std::string resolved_sample = selected_sample;
    if (resolved_sample.empty() && default_to_first_sample) {
        if (merged_header.samples.empty())
            throw std::runtime_error(
                "BAD_INPUT: tumor input header has no @RG SM sample (GATK requires a sample)");
        resolved_sample = merged_header.samples.front();
    }
    std::size_t selected_input_count = 0;
    for (const auto& input_header : input_headers) {
        if (resolved_sample.empty() ||
            std::find(input_header.samples.begin(), input_header.samples.end(), resolved_sample) !=
                input_header.samples.end())
            ++selected_input_count;
    }
    if (!resolved_sample.empty() && selected_input_count == 0)
        throw std::runtime_error(
            "BAD_INPUT: requested sample is absent from tumor input headers: " + resolved_sample);
    if (effective_sample != nullptr) *effective_sample = resolved_sample;
    if (interval_file_inputs != nullptr) *interval_file_inputs = 0;
    if (interval_file_records != nullptr) *interval_file_records = 0;
    if (adaptive_batch_reductions != nullptr) *adaptive_batch_reductions = 0;
    if (effective_batch_records != nullptr) *effective_batch_records = 0;
    if (requested_intervals != nullptr)
        *requested_intervals = input_intervals.empty() ? std::vector<fastgatk::io::HtsInterval>{}
                                                       : input_intervals.front();

    fastgatk::io::ReadBatch all;
    all.offsets.push_back(0);
    std::uint64_t accumulated_bytes = 0;
    const auto& effective_read_regions = read_regions == nullptr ? regions : *read_regions;
    for (std::size_t index = 0; index < paths.size(); ++index) {
        if (!resolved_sample.empty() &&
            std::find(input_headers[index].samples.begin(), input_headers[index].samples.end(),
                      resolved_sample) == input_headers[index].samples.end())
            continue;
        std::uint64_t file_reductions = 0;
        std::size_t file_effective_batch = 0;
        auto batch = read_all(paths[index], reference, effective_read_regions, batch_records,
                              input_headers[index], safe_memory_budget, resolved_sample,
                              false, nullptr, nullptr, nullptr, &file_reductions,
                              &file_effective_batch, nullptr,
                              true, interval_set_rule);
        accumulated_bytes = accumulated_bytes >
                std::numeric_limits<std::uint64_t>::max() - batch.bytes()
            ? std::numeric_limits<std::uint64_t>::max()
            : accumulated_bytes + batch.bytes();
        if (safe_memory_budget != 0 && accumulated_bytes > safe_memory_budget)
            throw std::runtime_error(
                "RESOURCE_EXHAUSTED: repeated -I input staging exceeds safe memory budget");
        append_batch(all, batch);
        if (interval_file_inputs != nullptr)
            *interval_file_inputs += input_interval_file_inputs[index];
        if (interval_file_records != nullptr)
            *interval_file_records += input_interval_file_records[index];
        if (adaptive_batch_reductions != nullptr) *adaptive_batch_reductions += file_reductions;
        if (effective_batch_records != nullptr) {
            *effective_batch_records = *effective_batch_records == 0
                ? file_effective_batch : std::min(*effective_batch_records, file_effective_batch);
        }
    }
    header = std::move(merged_header);
    if (all.records() == 0 && !allow_empty)
        throw std::runtime_error("BAD_INPUT: no reads selected");
    return all;
}

std::vector<fastgatk::io::HtsInterval> expanded_assembly_fetch_intervals(
    const fastgatk::io::HeaderSummary& header,
    const std::vector<fastgatk::io::HtsInterval>& core_intervals,
    const std::uint32_t padding) {
    // An absent -L means an unbounded traversal, which is already the
    // HtsReader empty-region contract.  Do not manufacture a contig list in
    // that case: it would turn an intentional whole-input traversal into an
    // indexed query with subtly different unmapped/read-order behavior.
    if (core_intervals.empty()) return {};

    const auto context = static_cast<std::int64_t>(padding);
    std::vector<fastgatk::io::HtsInterval> expanded;
    expanded.reserve(core_intervals.size());
    for (const auto& core : core_intervals) {
        if (core.tid < 0 || static_cast<std::size_t>(core.tid) >= header.contigs.size() ||
            static_cast<std::size_t>(core.tid) >= header.contig_lengths.size())
            continue;
        const auto contig_length = header.contig_lengths[static_cast<std::size_t>(core.tid)];
        const auto start = std::max<std::int64_t>(0, core.start - context);
        const auto end = std::min<std::int64_t>(contig_length, core.end + context);
        if (end > start) expanded.push_back(fastgatk::io::HtsInterval{core.tid, start, end});
    }
    std::sort(expanded.begin(), expanded.end(), [](const auto& left, const auto& right) {
        return std::tie(left.tid, left.start, left.end) <
               std::tie(right.tid, right.start, right.end);
    });
    std::vector<fastgatk::io::HtsInterval> merged;
    for (const auto& interval : expanded) {
        if (!merged.empty() && merged.back().tid == interval.tid &&
            interval.start <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }

    return merged;
}

std::vector<std::string> assembly_fetch_regions(
    const fastgatk::io::HeaderSummary& header,
    const std::vector<fastgatk::io::HtsInterval>& core_intervals,
    const std::uint32_t padding) {
    const auto intervals = expanded_assembly_fetch_intervals(
        header, core_intervals, padding);
    std::vector<std::string> regions;
    regions.reserve(intervals.size());
    for (const auto& interval : intervals) {
        const auto& contig = header.contigs[static_cast<std::size_t>(interval.tid)];
        regions.push_back(contig + ":" + std::to_string(interval.start + 1) +
                          "-" + std::to_string(interval.end));
    }
    return regions;
}

std::string allele_ref(const fastgatk::calling::AssemblyCandidate& candidate) {
    return candidate.reference_allele.empty()
        ? std::string(1, static_cast<char>(candidate.reference)) : candidate.reference_allele;
}

std::string allele_alt(const fastgatk::calling::AssemblyCandidate& candidate) {
    return candidate.alternate_allele.empty()
        ? std::string(1, static_cast<char>(candidate.alternate)) : candidate.alternate_allele;
}

std::size_t candidate_index(const fastgatk::calling::Result& result,
                            const fastgatk::calling::AssemblyCandidate& candidate) {
    std::size_t selected = std::numeric_limits<std::size_t>::max();
    int selected_priority = -1;
    const bool event_map_ownership_available =
        result.likelihood_candidate_event_map_owned.size() == result.candidates.size();
    for (std::size_t index = 0; index < result.candidates.size(); ++index) {
        const auto& item = result.candidates[index];
        if (item.tid != candidate.tid || item.position != candidate.position ||
            allele_ref(item) != allele_ref(candidate) || allele_alt(item) != allele_alt(candidate))
            continue;
        // The local Host candidate inventory may contain both a raw pileup
        // description and the EventMap-materialized description of the same
        // allele.  Every source-facing consumer (TLOD, annotations and VCF
        // writing) must resolve that identity to the EventMap row when it is
        // available; choosing the first raw duplicate loses the PairHMM
        // context and produces incorrect F1R2/F2R1 despite correct TLOD.
        int priority = 0;
        if (event_map_ownership_available &&
            result.likelihood_candidate_event_map_owned[index] != 0)
            priority += 16;
        if (item.graph_derived) priority += 4;
        if (item.forced_by_tumor == candidate.forced_by_tumor) priority += 2;
        if (item.somatic_low_quality_assembly_recovery ==
            candidate.somatic_low_quality_assembly_recovery)
            ++priority;
        if (priority > selected_priority) {
            selected_priority = priority;
            selected = index;
        }
    }
    return selected;
}

std::vector<fastgatk::calling::ForcedAllele> tumor_forced_alleles(
    const fastgatk::calling::Result& tumor,
    const fastgatk::kernels::SomaticLikelihoodResult* tumor_somatic,
    const double tumor_lod_emit) {
    std::vector<fastgatk::calling::ForcedAllele> forced;
    // This is the tumor's pre-normal somatic emission set: it includes every
    // candidate above Mutect2's tumor LOD boundary, but excludes uncalled
    // assembly branches. GATK's normal likelihood matrix needs the former,
    // while injecting all raw graph candidates needlessly creates a large
    // normal-only haplotype cross-product.
    if (tumor_somatic == nullptr) return forced;
    forced.reserve(tumor.candidates.size());
    for (std::size_t index = 0; index < tumor.candidates.size(); ++index) {
        const auto& candidate = tumor.candidates[index];
        if (index >= tumor_somatic->tlod.size() ||
            !std::isfinite(tumor_somatic->tlod[index]) ||
            tumor_somatic->tlod[index] <= tumor_lod_emit)
            continue;
        const auto reference = allele_ref(candidate);
        const auto alternate = allele_alt(candidate);
        if (candidate.tid < 0 || candidate.position < 0 || reference.empty() ||
            alternate.empty() || reference == alternate)
            continue;
        forced.push_back(fastgatk::calling::ForcedAllele{
            candidate.tid, candidate.position, reference, alternate,
            candidate.somatic_low_quality_assembly_recovery,
            candidate.linked_insertion_reference_position,
            candidate.linked_insertion_bases,
            candidate.reference_count,
            candidate.alternate_count,
            candidate.graph_derived});
    }
    return forced;
}

// An AssemblyResult EventMap is owned by the complete source-to-sink paths
// materialized in its graph.  Its CIGAR-derived per-haplotype map is the
// authority for that ownership: a preceding insertion/deletion shifts the
// raw path string, so indexing it directly at a reference coordinate can
// reject a valid downstream SNP.  Keep this Host predicate purely as an
// EventMap lookup; it neither synthesizes a haplotype nor performs any
// likelihood work.
bool graph_materializes_event(const fastgatk::calling::AssemblyRegionAssembly& assembly,
                              const fastgatk::calling::ForcedAllele& event) {
    if (event.tid < 0 || event.position < 0 || event.reference.empty() ||
        event.alternate.empty())
        return false;
    const auto same_event = [&](const fastgatk::calling::ForcedAllele& owner) {
        return owner.tid == event.tid && owner.position == event.position &&
            owner.reference == event.reference && owner.alternate == event.alternate;
    };
    const auto& graph = assembly.graph;
    // Every current AssemblyResult snapshot carries these maps.  Preserve a
    // narrow sequence fallback only for older callers that supplied a graph
    // without provenance, rather than letting it override an explicit map.
    if (!assembly.haplotype_event_maps.empty()) {
        for (std::size_t path = 0; path < assembly.haplotype_event_maps.size(); ++path) {
            if (path >= graph.haplotype_path_tids.size() ||
                path >= graph.haplotype_path_has_non_reference_edge.size() ||
                graph.haplotype_path_tids[path] != event.tid ||
                graph.haplotype_path_has_non_reference_edge[path] == 0)
                continue;
            if (std::any_of(assembly.haplotype_event_maps[path].events.begin(),
                            assembly.haplotype_event_maps[path].events.end(), same_event))
                return true;
        }
        return false;
    }
    for (std::size_t path = 0; path < graph.haplotype_path_sequences.size(); ++path) {
        if (path >= graph.haplotype_path_tids.size() ||
            path >= graph.haplotype_path_starts.size() ||
            path >= graph.haplotype_path_has_non_reference_edge.size() ||
            graph.haplotype_path_tids[path] != event.tid ||
            graph.haplotype_path_has_non_reference_edge[path] == 0)
            continue;
        const auto path_start = graph.haplotype_path_starts[path];
        if (event.position < path_start) continue;
        const auto offset = static_cast<std::size_t>(event.position - path_start);
        const auto& sequence = graph.haplotype_path_sequences[path];
        if (offset + event.alternate.size() <= sequence.size() &&
            sequence.compare(offset, event.alternate.size(), event.alternate) == 0 &&
            (offset + event.reference.size() > sequence.size() ||
             sequence.compare(offset, event.reference.size(), event.reference) != 0))
            return true;
    }
    return false;
}

// Preserve the joint AssemblyResult EventMap through per-sample replay.  The
// preliminary tumor TLOD is not an EventMap pruning criterion: GATK keeps a
// graph-materialized allele, scores it with the final sample-local PairHMM
// matrix, and only then applies the somatic emission gate.  Retain the narrow
// tumor recovery fallback for compact paths whose coordinate projection is
// not yet available, but never promote arbitrary pileup candidates.
void retain_eventmap_projection_and_tumor_recovery(
    std::vector<fastgatk::calling::AssemblyRegionAssembly>& assemblies,
    const fastgatk::calling::Result& tumor,
    const fastgatk::kernels::SomaticLikelihoodResult* tumor_somatic,
    const double tumor_lod_emit) {
    using Key = std::tuple<std::int32_t, std::int32_t, std::string, std::string>;
    std::set<Key> tumor_recovery;
    const auto recovered_events = tumor_forced_alleles(
        tumor, tumor_somatic, tumor_lod_emit);
    // TLOD is an emission calculation, not evidence that a raw pileup row
    // belongs to the joint AssemblyResult EventMap.
    // Keep only the explicit low-BQ compound recovery here: it is the one
    // deliberately documented exception whose upstream insertion branch was
    // assembled but whose downstream SNP is not yet projected by the compact
    // graph.  Admitting every TLOD-passing simple pileup candidate bypasses
    // GATK's EventMap pruning and, on HCC1143 chr20, turns the unassembled
    // 361521/361522 substitutions into two native-only calls.
    for (const auto& event : recovered_events) {
        if (!event.somatic_low_quality_assembly_recovery) continue;
        tumor_recovery.emplace(event.tid, event.position, event.reference, event.alternate);
    }
    for (auto& assembly : assemblies) {
        const auto before = assembly.event_map.size();
        assembly.event_map.erase(std::remove_if(assembly.event_map.begin(), assembly.event_map.end(),
            [&](const auto& event) {
                const Key key{event.tid, event.position, event.reference, event.alternate};
                const auto graph_owned = graph_materializes_event(assembly, event);
                const auto tumor_owned = tumor_recovery.find(key) != tumor_recovery.end();
                if (std::getenv("FASTGATK_DEBUG_EVENTMAP_LIFECYCLE") != nullptr) {
                    std::cerr << "[FASTGATK_EVENTMAP_RETAIN] region="
                              << assembly.region.tid << ':' << assembly.region.start
                              << '-' << assembly.region.end
                              << " tid=" << event.tid
                              << " position=" << event.position
                              << " ref=" << event.reference
                              << " alt=" << event.alternate
                              << " graph_owned=" << (graph_owned ? 1 : 0)
                              << " tumor_owned=" << (tumor_owned ? 1 : 0) << '\n';
                }
                return !graph_owned && !tumor_owned;
            }), assembly.event_map.end());
        if (std::getenv("FASTGATK_DEBUG_EVENTMAP_LIFECYCLE") != nullptr) {
            std::cerr << "[FASTGATK_EVENTMAP_RETAIN_SUMMARY] region="
                      << assembly.region.tid << ':' << assembly.region.start
                      << '-' << assembly.region.end
                      << " before=" << before
                      << " after=" << assembly.event_map.size() << '\n';
        }
    }
    // A compound low-BQ recovery can be a valid tumor EventMap allele even
    // when the compact joint graph did not project its downstream SNP edge.
    // The historical one-graph path recovered this exact already-emittable
    // allele from the tumor result.  Preserve that behavior at the only
    // source-aligned scope available here: its owning AssemblyRegion.  Do
    // not generalize this to ordinary TLOD-passing pileup rows, which would
    // bypass EventMap pruning and reintroduce cross-region combinations.
    for (const auto& event : recovered_events) {
        if (!event.somatic_low_quality_assembly_recovery) continue;
        const auto already_present = [&](const auto& assembly) {
            return std::any_of(assembly.event_map.begin(), assembly.event_map.end(),
                [&](const auto& existing) {
                    return existing.tid == event.tid && existing.position == event.position &&
                        existing.reference == event.reference && existing.alternate == event.alternate;
                });
        };
        auto owner = std::find_if(assemblies.begin(), assemblies.end(), [&](const auto& assembly) {
            const auto& region = assembly.region;
            return region.tid == event.tid && event.position >= region.active_start &&
                event.position <= region.active_end;
        });
        if (owner == assemblies.end()) {
            owner = std::find_if(assemblies.begin(), assemblies.end(), [&](const auto& assembly) {
                const auto& region = assembly.region;
                return region.tid == event.tid && event.position >= region.start &&
                    event.position <= region.end;
            });
        }
        if (owner == assemblies.end() && event.linked_insertion_reference_position >= 0) {
            // The recovered SNP is projected downstream of a concrete
            // insertion carried by the same assembled haplotype.  Its active
            // point can sit beyond the compact activity island even though
            // the insertion that owns the graph lies inside it.
            owner = std::find_if(assemblies.begin(), assemblies.end(), [&](const auto& assembly) {
                const auto& region = assembly.region;
                return region.tid == event.tid &&
                    event.linked_insertion_reference_position >= region.start &&
                    event.linked_insertion_reference_position <= region.end;
            });
        }
        if (owner == assemblies.end()) {
            // AssemblyRegionTrimmer uses 75 bp indel genotyping padding.
            // When a low-BQ SNP is the downstream projection of a dangling
            // insertion, native's compact activity profile may terminate at
            // the insertion itself.  Assign only within that source padding
            // to the nearest local AssemblyResult; more distant candidates
            // remain excluded rather than becoming free-standing pileup rows.
            constexpr std::int32_t kIndelGenotypingPadding = 75;
            const auto distance = [&](const auto& assembly) {
                const auto& region = assembly.region;
                if (region.tid != event.tid) return std::numeric_limits<std::int32_t>::max();
                if (event.position < region.start) return region.start - event.position;
                if (event.position > region.end) return event.position - region.end;
                return 0;
            };
            const auto nearest = std::min_element(assemblies.begin(), assemblies.end(),
                [&](const auto& left, const auto& right) {
                    return distance(left) < distance(right);
                });
            if (nearest != assemblies.end() && distance(*nearest) <= kIndelGenotypingPadding)
                owner = nearest;
        }
        if (owner != assemblies.end() && !already_present(*owner)) {
            // Keep this explicit EventMap event inside the replay/output
            // envelope.  Its active-core ownership remains the upstream
            // insertion region, while the ordinary halo merge resolves the
            // downstream SNP deterministically.
            owner->region.start = std::min(owner->region.start, event.position);
            owner->region.end = std::max(owner->region.end, event.position);
            owner->event_map.push_back(event);
        }
    }
}

// GATK's AssemblyResultSet exposes variation events by projecting each
// retained haplotype back onto the reference (EventMap.fromHaplotype), not by
// replaying every intermediate pileup/dangling candidate.  `graph_derived`
// is the native equivalent of that projected EventMap membership: it is set
// only after the retained graph path has passed the reference-connected
// SW/CIGAR projection in calling::run().
//
// Mutect2's optional PileupBasedAlleles union is disabled by default in GATK
// 4.6.  Until every compact-graph path projects correctly, the narrow
// read-supported-SNP fallback below covers only missing EventMap projections.
// A final TLOD-gated fallback remains quarantined below for paths native has
// not yet materialized at all; it is deliberately not treated as EventMap
// truth and will be removed once those projections are complete.
std::vector<fastgatk::calling::ForcedAllele> event_map_forced_alleles(
    const fastgatk::calling::Result& joint_assembly,
    const fastgatk::calling::Result& tumor,
    const fastgatk::kernels::SomaticLikelihoodResult* tumor_somatic,
    const double tumor_lod_emit) {
    using AlleleKey = std::tuple<std::int32_t, std::int32_t, std::string, std::string>;
    std::vector<fastgatk::calling::ForcedAllele> forced;
    std::set<AlleleKey> inserted;
    const auto append = [&](const fastgatk::calling::AssemblyCandidate& candidate) {
        const auto reference = allele_ref(candidate);
        const auto alternate = allele_alt(candidate);
        if (candidate.tid < 0 || candidate.position < 0 || reference.empty() ||
            alternate.empty() || reference == alternate)
            return;
        const AlleleKey key{candidate.tid, candidate.position, reference, alternate};
        if (!inserted.insert(key).second) return;
        forced.push_back(fastgatk::calling::ForcedAllele{
            candidate.tid, candidate.position, reference, alternate,
            candidate.somatic_low_quality_assembly_recovery,
            candidate.linked_insertion_reference_position,
            candidate.linked_insertion_bases,
            candidate.reference_count,
            candidate.alternate_count});
    };
    const auto tumor_recovery = tumor_forced_alleles(
        tumor, tumor_somatic, tumor_lod_emit);
    std::set<AlleleKey> tumor_recovery_keys;
    for (const auto& candidate : tumor_recovery) {
        tumor_recovery_keys.emplace(candidate.tid, candidate.position,
                                    candidate.reference, candidate.alternate);
    }
    forced.reserve(joint_assembly.candidates.size() + tumor.candidates.size());
    const auto belongs_to_joint_assembly_region =
        [&](const fastgatk::calling::AssemblyCandidate& candidate) {
            return std::any_of(joint_assembly.calling_regions.begin(),
                               joint_assembly.calling_regions.end(),
                               [&](const auto& region) {
                return candidate.tid == region.tid &&
                    candidate.position >= region.start &&
                    candidate.position <= region.end;
            });
        };
    for (const auto& candidate : joint_assembly.candidates) {
        const auto reference = allele_ref(candidate);
        const auto alternate = allele_alt(candidate);
        const bool supported_simple_substitution =
            reference.size() == 1 && alternate.size() == 1 &&
            candidate.alternate_count > 0;
        // The graph EventMap projection is authoritative.  Native's compact
        // graph still misses a subset of concrete SNP paths (the chr17
        // contract has three such retained paths), so carry a strictly
        // read-supported single-base assembly candidate as an interim
        // projection fallback.  In particular, do not admit the zero-support
        // long/dangling bookkeeping rows that are present in Result::candidates.
        // AssemblyResultSet.trimTo() exposes events from the active
        // AssemblyRegion only.  A compact graph path that SW anchors outside
        // every joint region is an untrimmed traversal artifact, not an
        // EventMap allele for the matched-normal likelihood matrix.
        const AlleleKey key{candidate.tid, candidate.position, reference, alternate};
        if (belongs_to_joint_assembly_region(candidate) &&
            (candidate.graph_derived || supported_simple_substitution) &&
            tumor_recovery_keys.find(key) != tumor_recovery_keys.end())
            append(candidate);
    }
    for (const auto& candidate : tumor_recovery) {
        const AlleleKey key{candidate.tid, candidate.position,
                            candidate.reference, candidate.alternate};
        if (inserted.insert(key).second) forced.push_back(candidate);
    }
    return forced;
}

fastgatk::kernels::KmerGraphResult assembly_graph(
    const fastgatk::calling::Result& assembly) {
    fastgatk::kernels::KmerGraphResult result;
    result.used = assembly.graph_used;
    result.kmer_size = assembly.graph_kmer_size_selected;
    result.kmer_iterations = assembly.graph_kmer_iterations;
    result.has_non_reference_cycles = assembly.graph_has_non_reference_cycles;
    result.haplotype_path_count = assembly.graph_haplotype_paths;
    result.haplotype_path_sequences = assembly.graph_haplotype_sequences;
    result.haplotype_path_tids = assembly.graph_haplotype_path_tids;
    result.haplotype_path_starts = assembly.graph_haplotype_path_starts;
    result.haplotype_path_ends = assembly.graph_haplotype_path_ends;
    result.haplotype_path_support = assembly.graph_haplotype_path_support;
    result.haplotype_path_scores = assembly.graph_haplotype_path_scores;
    result.haplotype_path_has_non_reference_edge =
        assembly.graph_haplotype_path_has_non_reference_edge;
    result.haplotype_path_alt_read_starts = assembly.graph_haplotype_path_alt_read_starts;
    result.haplotype_path_alt_read_ends = assembly.graph_haplotype_path_alt_read_ends;
    return result;
}

// Return whether one tumor ALT remains eligible under GATK's --normal-lod
// boundary.  Matched-normal scoring receives every tumor allele before its
// PairHMM pass, so a candidate row must be present whenever normal evidence
// exists; do not restrict the lookup to normal *emitted* calls.
bool normal_lod_allows_call(
    const double threshold,
    const fastgatk::calling::GenotypeCall& tumor_call,
    const fastgatk::calling::Result* normal,
    const fastgatk::kernels::SomaticLikelihoodResult* normal_somatic) {
    if (normal == nullptr) return true;
    const auto index = candidate_index(*normal, tumor_call.candidate);
    if (index == std::numeric_limits<std::size_t>::max()) return true;
    double nlod = index < normal->likelihoods.size()
        ? normal->likelihoods[index].hom_ref - normal->likelihoods[index].het
        : -std::numeric_limits<double>::infinity();
    if (normal_somatic != nullptr && index < normal_somatic->normal_log10_odds.size())
        nlod = normal_somatic->normal_log10_odds[index];
    return std::isfinite(nlod) && nlod > threshold;
}

// The all-low-BQ compound-haplotype recovery is confirmed by the local
// assembly and is ultimately emitted by SomaticGenotypingEngine's TLOD, not
// by the generic diploid hard-call gate used by the reusable pipeline.  Keep
// that exceptional route explicit: ordinary candidates still need a concrete
// pipeline GenotypeCall before reaching Mutect2 output.
void materialize_low_quality_assembly_recovery_calls(
    fastgatk::calling::Result& tumor,
    const fastgatk::kernels::SomaticLikelihoodResult* somatic,
    const double tumor_lod_emit) {
    if (somatic == nullptr) return;
    for (std::size_t index = 0; index < tumor.candidates.size(); ++index) {
        const auto& candidate = tumor.candidates[index];
        if (!candidate.somatic_low_quality_assembly_recovery ||
            index >= somatic->tlod.size() || !std::isfinite(somatic->tlod[index]) ||
            somatic->tlod[index] <= tumor_lod_emit)
            continue;
        const auto already_called = std::any_of(tumor.calls.begin(), tumor.calls.end(),
            [&](const auto& call) {
                return candidate_index(tumor, call.candidate) == index;
            });
        if (already_called) continue;
        fastgatk::calling::GenotypeCall call;
        call.candidate = candidate;
        if (index < tumor.likelihoods.size()) call.likelihoods = tumor.likelihoods[index];
        call.genotype = 1;
        call.qual = std::max(0.0, somatic->tlod[index] * 10.0);
        tumor.calls.push_back(std::move(call));
    }
}

// In reference-confidence mode Mutect2's negative-infinity emission threshold
// is applied to every EventMap allele, not only to candidates that survived
// the reusable diploid hard-call gate.  The latter gate remains correct for
// ordinary VCF mode, but it would otherwise turn low-TLOD graph alleles into
// reference blocks and erase the somatic gVCF's allelic context.  Keep the
// decision and the construction on the Host: EventMap provenance and
// GenotypeCall ownership are Host state, while the candidate likelihoods and
// TLODs remain the existing Kokkos results.
void materialize_reference_confidence_eventmap_calls(
    fastgatk::calling::Result& tumor,
    const fastgatk::kernels::SomaticLikelihoodResult* somatic,
    const bool emit_reference_confidence) {
    if (!emit_reference_confidence || somatic == nullptr) return;
    const bool has_explicit_event_map_ownership =
        tumor.likelihood_candidate_event_map_owned.size() == tumor.candidates.size();
    for (std::size_t index = 0; index < tumor.candidates.size(); ++index) {
        const auto& candidate = tumor.candidates[index];
        const bool event_map_owned = has_explicit_event_map_ownership
            ? tumor.likelihood_candidate_event_map_owned[index] != 0U
            : candidate.graph_derived;
        if (!event_map_owned || index >= somatic->tlod.size() ||
            !std::isfinite(somatic->tlod[index]))
            continue;
        const auto belongs_to_active_event_map = [&] {
            // A flattened Result contains EventMap alleles from every
            // AssemblyRegion plus its padded halos.  GATK
            // AssemblyResultSet.trimTo() only hands callMutations an allele
            // from the region whose *active core* owns it.  Testing merely
            // against the union of global intervals incorrectly promotes a
            // right-halo event from the preceding region into the next core.
            // Resolve the ownership in the local Result before emitting the
            // reference-confidence VariantContext.
            if (!tumor.assembly_region_likelihood_results.empty()) {
                for (const auto& part_pointer : tumor.assembly_region_likelihood_results) {
                    if (part_pointer == nullptr) continue;
                    const auto& part = *part_pointer;
                    const bool part_has_explicit_ownership =
                        part.likelihood_candidate_event_map_owned.size() == part.candidates.size();
                    for (std::size_t local = 0; local < part.candidates.size(); ++local) {
                        const auto& local_candidate = part.candidates[local];
                        if (local_candidate.tid != candidate.tid ||
                            local_candidate.position != candidate.position ||
                            allele_ref(local_candidate) != allele_ref(candidate) ||
                            allele_alt(local_candidate) != allele_alt(candidate))
                            continue;
                        const bool local_event_map_owned = part_has_explicit_ownership
                            ? part.likelihood_candidate_event_map_owned[local] != 0U
                            : local_candidate.graph_derived;
                        if (!local_event_map_owned) continue;
                        for (const auto& region : part.calling_regions) {
                            if (local_candidate.tid == region.tid &&
                                local_candidate.position >= region.active_start &&
                                local_candidate.position <= region.active_end)
                                return true;
                        }
                    }
                }
                return false;
            }
            if (tumor.calling_regions.empty()) return true;
            return std::any_of(tumor.calling_regions.begin(), tumor.calling_regions.end(),
                [&](const auto& region) {
                    return candidate.tid == region.tid &&
                        candidate.position >= region.active_start &&
                        candidate.position <= region.active_end;
                });
        };
        if (!belongs_to_active_event_map()) continue;
        const auto already_called = std::any_of(tumor.calls.begin(), tumor.calls.end(),
            [&](const auto& call) { return candidate_index(tumor, call.candidate) == index; });
        if (already_called) continue;
        fastgatk::calling::GenotypeCall call;
        call.candidate = candidate;
        if (index < tumor.likelihoods.size()) call.likelihoods = tumor.likelihoods[index];
        call.genotype = 1;
        call.qual = std::max(0.0, somatic->tlod[index] * 10.0);
        tumor.calls.push_back(std::move(call));
    }
}

struct GroupedSomaticEvidence {
    std::size_t read_count = 0;
    bool grouped = false;
    // The established dense vectors remain only for the legacy ungrouped
    // fallback below.  Production Mutect2 PairHMM output is candidate-major
    // CSR so wide windows do not materialize candidate x fragment matrices
    // merely to hand values back to a Kokkos reduction.
    bool sparse_rows = false;
    const std::vector<std::uint32_t>* sparse_candidate_offsets = nullptr;
    const std::vector<fastgatk::calling::SomaticGroupLikelihood>* sparse_likelihood_rows = nullptr;
    std::vector<double> reference;
    std::vector<double> alternate;
    std::vector<double> non_reference;

    double reference_at(const std::size_t candidate, const std::size_t group) const {
        constexpr double kMissing = -1.0e299;
        if (!sparse_rows) {
            const auto index = candidate * read_count + group;
            return index < reference.size() ? reference[index] : kMissing;
        }
        if (sparse_candidate_offsets == nullptr || sparse_likelihood_rows == nullptr ||
            candidate + 1 >= sparse_candidate_offsets->size()) return kMissing;
        const auto begin = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate]);
        const auto end = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate + 1]);
        if (begin > end || end > sparse_likelihood_rows->size()) return kMissing;
        const auto found = std::lower_bound(sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(begin),
            sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end), static_cast<std::uint32_t>(group),
            [](const auto& row, const auto queried_group) { return row.group < queried_group; });
        return found != sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end) &&
                found->group == group ? found->reference : kMissing;
    }

    double alternate_at(const std::size_t candidate, const std::size_t group) const {
        constexpr double kMissing = -1.0e299;
        if (!sparse_rows) {
            const auto index = candidate * read_count + group;
            return index < alternate.size() ? alternate[index] : kMissing;
        }
        if (sparse_candidate_offsets == nullptr || sparse_likelihood_rows == nullptr ||
            candidate + 1 >= sparse_candidate_offsets->size()) return kMissing;
        const auto begin = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate]);
        const auto end = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate + 1]);
        if (begin > end || end > sparse_likelihood_rows->size()) return kMissing;
        const auto found = std::lower_bound(sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(begin),
            sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end), static_cast<std::uint32_t>(group),
            [](const auto& row, const auto queried_group) { return row.group < queried_group; });
        return found != sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end) &&
                found->group == group ? found->alternate : kMissing;
    }

    double spanning_deletion_at(const std::size_t candidate, const std::size_t group) const {
        constexpr double kMissing = -1.0e299;
        // The legacy dense fallback predates GATK's spanning-deletion
        // allele mapper and has no `*` row. Production PairHMM output is the
        // sparse EventMap representation above.
        if (!sparse_rows) return kMissing;
        if (sparse_candidate_offsets == nullptr || sparse_likelihood_rows == nullptr ||
            candidate + 1 >= sparse_candidate_offsets->size()) return kMissing;
        const auto begin = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate]);
        const auto end = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate + 1]);
        if (begin > end || end > sparse_likelihood_rows->size()) return kMissing;
        const auto found = std::lower_bound(sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(begin),
            sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end), static_cast<std::uint32_t>(group),
            [](const auto& row, const auto queried_group) { return row.group < queried_group; });
        return found != sparse_likelihood_rows->begin() + static_cast<std::ptrdiff_t>(end) &&
                found->group == group ? found->spanning_deletion : kMissing;
    }

    std::vector<std::uint32_t> groups_for_candidate(const std::size_t candidate) const {
        std::vector<std::uint32_t> groups;
        if (!sparse_rows) {
            groups.resize(read_count);
            for (std::size_t group = 0; group < read_count; ++group)
                groups[group] = static_cast<std::uint32_t>(group);
            return groups;
        }
        if (sparse_candidate_offsets == nullptr || sparse_likelihood_rows == nullptr ||
            candidate + 1 >= sparse_candidate_offsets->size()) return groups;
        const auto begin = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate]);
        const auto end = static_cast<std::size_t>((*sparse_candidate_offsets)[candidate + 1]);
        if (begin > end || end > sparse_likelihood_rows->size()) return groups;
        groups.reserve(end - begin);
        for (std::size_t row = begin; row < end; ++row)
            groups.push_back((*sparse_likelihood_rows)[row].group);
        return groups;
    }
};

// GATK's SomaticGenotypingEngine calls AlleleLikelihoods.groupEvidence before
// it evaluates tumor/normal likelihoods.  The grouping key is GATKRead::name
// and the per-allele log likelihoods of all reads in one fragment are summed.
// Result keeps the post-filter read-name vector, so this Host reduction can
// reproduce the same boundary without passing HTSlib objects to Kokkos.
std::optional<GroupedSomaticEvidence> group_somatic_evidence(
    const fastgatk::calling::Result& result) {
    if (result.candidates.empty() || result.reference_read_likelihoods.size() != result.candidates.size() ||
        result.allele_read_likelihoods.size() != result.candidates.size())
        return std::nullopt;
    const std::size_t read_count = result.reference_read_likelihoods.empty()
        ? 0 : result.reference_read_likelihoods.front().size();
    if (result.somatic_likelihoods_grouped &&
        result.somatic_group_representative_sources.size() != 0 &&
        result.somatic_candidate_group_offsets.size() == result.candidates.size() + 1 &&
        result.somatic_candidate_group_offsets.front() == 0 &&
        result.somatic_candidate_group_offsets.back() == result.somatic_group_likelihoods.size()) {
        for (std::size_t candidate = 0; candidate < result.candidates.size(); ++candidate) {
            const auto begin = static_cast<std::size_t>(result.somatic_candidate_group_offsets[candidate]);
            const auto end = static_cast<std::size_t>(result.somatic_candidate_group_offsets[candidate + 1]);
            if (begin > end || end > result.somatic_group_likelihoods.size())
                return std::nullopt;
            for (std::size_t row = begin + 1; row < end; ++row)
                if (result.somatic_group_likelihoods[row - 1].group >=
                    result.somatic_group_likelihoods[row].group)
                    return std::nullopt;
        }
        GroupedSomaticEvidence grouped;
        grouped.read_count = result.somatic_group_representative_sources.size();
        grouped.grouped = true;
        grouped.sparse_rows = true;
        grouped.sparse_candidate_offsets = &result.somatic_candidate_group_offsets;
        grouped.sparse_likelihood_rows = &result.somatic_group_likelihoods;
        return grouped;
    }
    if (read_count == 0) return std::nullopt;

    std::vector<std::size_t> group_for_read(read_count, 0);
    std::map<std::string, std::size_t> group_indices;
    const bool names_available = result.likelihood_read_names.size() == read_count;
    std::size_t group_count = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        // Empty/missing names fail closed to singleton evidence. This is
        // equivalent to Fragment::createAndAvoidFailure for legacy batches
        // where a stable fragment key is unavailable.
        const auto& name = names_available ? result.likelihood_read_names[read] : std::string{};
        if (name.empty()) {
            group_for_read[read] = group_count++;
            continue;
        }
        const auto [iter, inserted] = group_indices.emplace(name, group_count);
        if (inserted) ++group_count;
        group_for_read[read] = iter->second;
    }
    if (group_count == 0) return std::nullopt;

    constexpr double kMissing = -1.0e299;
    GroupedSomaticEvidence grouped;
    grouped.read_count = group_count;
    grouped.grouped = group_count < read_count;
    grouped.reference.assign(result.candidates.size() * group_count, kMissing);
    grouped.alternate.assign(result.candidates.size() * group_count, kMissing);
    const bool has_non_reference = result.non_ref_read_likelihoods.size() == result.candidates.size();
    if (has_non_reference)
        grouped.non_reference.assign(result.candidates.size() * group_count, kMissing);
    std::vector<std::uint8_t> reference_seen(grouped.reference.size(), 0);
    std::vector<std::uint8_t> alternate_seen(grouped.alternate.size(), 0);
    std::vector<std::uint8_t> non_reference_seen(grouped.non_reference.size(), 0);
    const bool coordinates_available = result.likelihood_read_tids.size() == read_count &&
        result.likelihood_read_starts.size() == read_count &&
        result.likelihood_read_ends.size() == read_count;
    for (std::size_t candidate = 0; candidate < result.candidates.size(); ++candidate) {
        if (result.reference_read_likelihoods[candidate].size() != read_count ||
            result.allele_read_likelihoods[candidate].size() != read_count)
            return std::nullopt;
        const bool candidate_has_non_reference = has_non_reference &&
            result.non_ref_read_likelihoods[candidate].size() == read_count;
        // GATK groups reads into Fragments first and then retains a fragment
        // when its span overlaps the merged variant expanded by
        // informativeReadOverlapMargin (default 2). A fragment may therefore
        // contribute both mates even when only one mate overlaps the event.
        // Mark groups before summing likelihood rows so the native matrix has
        // the same candidate-specific evidence boundary.
        std::vector<std::uint8_t> group_overlaps(group_count,
            coordinates_available ? 0 : 1);
        if (coordinates_available) {
            const auto reference_length = allele_ref(result.candidates[candidate]).size();
            const auto event_start = std::max<std::int64_t>(
                0, static_cast<std::int64_t>(result.candidates[candidate].position) - 2);
            const auto event_end = std::max<std::int64_t>(
                event_start + 1,
                static_cast<std::int64_t>(result.candidates[candidate].position) +
                static_cast<std::int64_t>(reference_length == 0 ? 1 : reference_length) + 2);
            for (std::size_t read = 0; read < read_count; ++read) {
                const auto tid = result.likelihood_read_tids[read];
                const auto start = result.likelihood_read_starts[read];
                const auto end = result.likelihood_read_ends[read];
                if (tid == result.candidates[candidate].tid && start >= 0 && end > start &&
                    start < event_end && end > event_start)
                    group_overlaps[group_for_read[read]] = 1;
            }
        }
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto group = group_for_read[read];
            if (group >= group_overlaps.size() || group_overlaps[group] == 0)
                continue;
            const auto index = candidate * group_count + group;
            const auto ref = result.reference_read_likelihoods[candidate][read];
            const auto alt = result.allele_read_likelihoods[candidate][read];
            if (std::isfinite(ref) && ref > kMissing) {
                if (reference_seen[index] == 0) {
                    grouped.reference[index] = ref;
                    reference_seen[index] = 1;
                } else {
                    grouped.reference[index] += ref;
                }
            }
            if (std::isfinite(alt) && alt > kMissing) {
                if (alternate_seen[index] == 0) {
                    grouped.alternate[index] = alt;
                    alternate_seen[index] = 1;
                } else {
                    grouped.alternate[index] += alt;
                }
            }
            if (candidate_has_non_reference) {
                const auto non_reference = result.non_ref_read_likelihoods[candidate][read];
                if (std::isfinite(non_reference) && non_reference > kMissing) {
                    if (non_reference_seen[index] == 0) {
                        grouped.non_reference[index] = non_reference;
                        non_reference_seen[index] = 1;
                    } else {
                        grouped.non_reference[index] += non_reference;
                    }
                }
            }
        }
        // A somatic candidate may have been restored from the Host pileup
        // boundary after the HC-oriented read×haplotype CIGAR gate rejected
        // every ALT owner.  Preserve its aggregate evidence as one stable
        // fragment so Mutect2 still computes a finite TLOD instead of silently
        // dropping the candidate from the somatic model.
        const auto fallback = candidate < result.likelihoods.size()
            ? result.likelihoods[candidate] : fastgatk::calling::Likelihoods{};
        // Keep the historical pileup fallback, but attach it to the first
        // fragment that actually survives the candidate overlap boundary;
        // attaching it unconditionally to group zero can silently erase a
        // low-support candidate when the first BAM read belongs elsewhere.
        const auto active_group = std::find(group_overlaps.begin(), group_overlaps.end(), 1U);
        if (active_group != group_overlaps.end()) {
            const auto group = static_cast<std::size_t>(active_group - group_overlaps.begin());
            const auto first = candidate * group_count + group;
            if (reference_seen[first] == 0 && std::isfinite(fallback.hom_ref)) {
                grouped.reference[first] = fallback.hom_ref;
                reference_seen[first] = 1;
            }
            if (alternate_seen[first] == 0 && std::isfinite(fallback.hom_alt)) {
                grouped.alternate[first] = fallback.hom_alt;
                alternate_seen[first] = 1;
            }
        }
    }
    return grouped;
}

struct GroupedSomaticDepth {
    std::vector<std::uint32_t> depth;
    std::vector<std::uint32_t> reference;
    std::vector<std::uint32_t> alternate;
};

// FragmentDepthPerAlleleBySample computes FORMAT/FAD after
// SomaticGenotypingEngine has grouped evidence, marginalised to concrete
// alleles, and retained only the informative-overlap fragments.  It then uses
// the ordinary BestAllele rule (strict natural-log margin 0.2, equivalent to
// the log10 margin above), rather than the
// fractional effective counts used transiently by addGenotypes() before the
// annotation engine overwrites FORMAT/AD and FORMAT/DP with read-level values.
std::optional<GroupedSomaticDepth> grouped_somatic_fragment_depth(
    const fastgatk::calling::Result& result,
    const std::vector<fastgatk::calling::AssemblyCandidate>* emitted_candidates = nullptr) {
    const auto grouped = group_somatic_evidence(result);
    if (!grouped.has_value()) return std::nullopt;
    const auto candidate_count = result.candidates.size();
    const auto group_count = grouped->read_count;
    if (candidate_count == 0 || group_count == 0) return std::nullopt;
    GroupedSomaticDepth output;
    output.depth.assign(candidate_count, 0U);
    output.reference.assign(candidate_count, 0U);
    output.alternate.assign(candidate_count, 0U);

    // FAD is calculated from the same EventMap-marginalized fragment allele
    // list as callMutations' other BestAllele annotations.  Result retains
    // raw candidate siblings for Host bookkeeping, but they must not become
    // additional fragment alleles merely because they share a locus.
    const bool has_pairhmm_event_map_ownership =
        result.likelihood_candidate_event_map_owned.size() == candidate_count;
    const auto event_map_owns_candidate = [&](const std::size_t candidate) {
        if (candidate >= candidate_count) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position &&
                    event.reference == allele_ref(item) && event.alternate == allele_alt(item);
            });
        });
    };
    const auto event_map_has_locus = [&](const std::size_t candidate) {
        if (candidate >= candidate_count) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position;
            });
        });
    };
    const auto belongs_to_fragment_allele_list = [&](const std::size_t candidate) {
        const auto& item = result.candidates[candidate];
        const bool recovery_without_event_map_locus =
            item.somatic_low_quality_assembly_recovery &&
            (result.graph_haplotype_event_maps.empty() || !event_map_has_locus(candidate));
        if (recovery_without_event_map_locus) return true;
        if (has_pairhmm_event_map_ownership)
            return result.likelihood_candidate_event_map_owned[candidate] != 0;
        if (!result.graph_haplotype_event_maps.empty())
            return event_map_owns_candidate(candidate);
        return true;
    };

    using SomaticGroupKey = std::tuple<std::int32_t, std::int32_t, std::string>;
    std::map<SomaticGroupKey, std::vector<std::size_t>> groups;
    if (emitted_candidates != nullptr) {
        // AnnotationEngine receives the post-emission VariantContext.  Do
        // not let a sibling EventMap ALT that failed emission compete in FAD.
        for (const auto& emitted : *emitted_candidates) {
            const auto candidate = candidate_index(result, emitted);
            if (candidate == std::numeric_limits<std::size_t>::max() ||
                candidate >= candidate_count)
                return std::nullopt;
            const auto& item = result.candidates[candidate];
            groups[SomaticGroupKey{item.tid, item.position, allele_ref(item)}].push_back(candidate);
        }
    } else {
        for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
            if (!belongs_to_fragment_allele_list(candidate)) continue;
            const auto& item = result.candidates[candidate];
            groups[SomaticGroupKey{item.tid, item.position, allele_ref(item)}].push_back(candidate);
        }
    }
    const auto saturating_increment = [](std::uint32_t& value) {
        if (value != std::numeric_limits<std::uint32_t>::max()) ++value;
    };
    constexpr double kMissing = -1.0e299;
    for (const auto& [key, indices] : groups) {
        (void)key;
        if (indices.empty()) continue;
        std::set<std::uint32_t> active_groups;
        for (const auto candidate : indices) {
            for (const auto group : grouped->groups_for_candidate(candidate))
                active_groups.insert(group);
        }
        for (const auto group : active_groups) {
            std::vector<double> values(indices.size() + 1, kMissing);
            for (const auto candidate : indices)
                values[0] = std::max(values[0], grouped->reference_at(candidate, group));
            for (std::size_t alternate = 0; alternate < indices.size(); ++alternate)
                values[alternate + 1] = grouped->alternate_at(indices[alternate], group);
            const auto best = gatk_mutect_best_allele(values, kMissing);
            if (!best.has_value()) continue;
            for (const auto candidate : indices) {
                saturating_increment(output.depth[candidate]);
                if (*best == 0) saturating_increment(output.reference[candidate]);
            }
            if (*best != 0 && *best - 1 < indices.size())
                saturating_increment(output.alternate[indices[*best - 1]]);
        }
    }
    return output;
}

void apply_grouped_somatic_depth(fastgatk::calling::Result& result) {
    // SomaticGenotypingEngine annotates each AssemblyRegion before the walker
    // merges its VariantContexts.  The flattened Result deliberately omits a
    // global read/fragment likelihood matrix because its read ordinals are
    // region-local.  Mirror the existing TLOD and annotation reducers: keep
    // the Kokkos-produced regional fragment likelihoods in their owning
    // Result, annotate there, then project those scalar FORMAT counts onto
    // the merged candidate list, preferring the active-region owner over a
    // padded halo.
    if (!result.assembly_region_likelihood_results.empty()) {
        const auto count = result.candidates.size();
        std::vector<std::uint32_t> depth(count, 0U);
        std::vector<std::uint32_t> reference(count, 0U);
        std::vector<std::uint32_t> alternate(count, 0U);
        std::vector<int> ownership(count, -1);
        for (const auto& part_pointer : result.assembly_region_likelihood_results) {
            if (!part_pointer) continue;
            auto& part = *part_pointer;
            apply_grouped_somatic_depth(part);
            for (std::size_t local_index = 0; local_index < part.candidates.size(); ++local_index) {
                const auto merged_index = candidate_index(result, part.candidates[local_index]);
                if (merged_index >= count) continue;
                int priority = 1;
                if (part.calling_regions.size() == 1) {
                    const auto& region = part.calling_regions.front();
                    const auto& candidate = part.candidates[local_index];
                    if (candidate.tid == region.tid &&
                        candidate.position >= region.active_start &&
                        candidate.position <= region.active_end)
                        priority = 2;
                }
                if (priority < ownership[merged_index]) continue;
                ownership[merged_index] = priority;
                if (local_index < part.somatic_fragment_depth.size())
                    depth[merged_index] = part.somatic_fragment_depth[local_index];
                if (local_index < part.somatic_fragment_reference_count.size())
                    reference[merged_index] = part.somatic_fragment_reference_count[local_index];
                if (local_index < part.somatic_fragment_alternate_count.size())
                    alternate[merged_index] = part.somatic_fragment_alternate_count[local_index];
            }
        }
        result.somatic_fragment_depth = std::move(depth);
        result.somatic_fragment_reference_count = std::move(reference);
        result.somatic_fragment_alternate_count = std::move(alternate);
        return;
    }
    const auto grouped = grouped_somatic_fragment_depth(result);
    if (!grouped.has_value() || grouped->depth.size() != result.candidates.size() ||
        grouped->reference.size() != result.candidates.size() ||
        grouped->alternate.size() != result.candidates.size())
        return;
    result.somatic_fragment_depth = grouped->depth;
    result.somatic_fragment_reference_count = grouped->reference;
    result.somatic_fragment_alternate_count = grouped->alternate;
}

struct EmittedSomaticFragmentDepth {
    std::uint32_t depth = 0U;
    // Allele-list order: shared REF followed by the emitted concrete ALTs.
    std::vector<std::uint32_t> allele_depths;
};

std::optional<EmittedSomaticFragmentDepth> emitted_somatic_fragment_depth_local(
    const fastgatk::calling::Result& result,
    const std::vector<fastgatk::calling::AssemblyCandidate>& emitted_candidates) {
    if (emitted_candidates.empty()) return std::nullopt;
    const auto grouped = grouped_somatic_fragment_depth(result, &emitted_candidates);
    if (!grouped.has_value()) return std::nullopt;
    EmittedSomaticFragmentDepth output;
    output.allele_depths.assign(emitted_candidates.size() + 1, 0U);
    for (std::size_t alternate = 0; alternate < emitted_candidates.size(); ++alternate) {
        const auto candidate = candidate_index(result, emitted_candidates[alternate]);
        if (candidate == std::numeric_limits<std::size_t>::max() ||
            candidate >= grouped->depth.size() || candidate >= grouped->reference.size() ||
            candidate >= grouped->alternate.size())
            return std::nullopt;
        if (alternate == 0) {
            output.depth = grouped->depth[candidate];
            output.allele_depths.front() = grouped->reference[candidate];
        }
        output.allele_depths[alternate + 1] = grouped->alternate[candidate];
    }
    return output;
}

std::optional<EmittedSomaticFragmentDepth> emitted_somatic_fragment_depth(
    const fastgatk::calling::Result& result,
    const std::vector<fastgatk::calling::AssemblyCandidate>& emitted_candidates) {
    if (result.assembly_region_likelihood_results.empty())
        return emitted_somatic_fragment_depth_local(result, emitted_candidates);

    const auto& first = emitted_candidates.front();
    const fastgatk::calling::Result* owner = nullptr;
    int owner_priority = -1;
    for (const auto& part_pointer : result.assembly_region_likelihood_results) {
        if (!part_pointer) continue;
        const auto& part = *part_pointer;
        const auto has_all = std::all_of(emitted_candidates.begin(), emitted_candidates.end(),
            [&](const auto& emitted) {
                return candidate_index(part, emitted) != std::numeric_limits<std::size_t>::max();
            });
        if (!has_all) continue;
        int priority = 1;
        for (const auto& region : part.calling_regions) {
            if (region.tid == first.tid && first.position >= region.active_start &&
                first.position <= region.active_end) {
                priority = 2;
                break;
            }
        }
        if (priority > owner_priority) {
            owner = &part;
            owner_priority = priority;
        }
    }
    return owner == nullptr ? std::nullopt
                            : emitted_somatic_fragment_depth_local(*owner, emitted_candidates);
}

// Mutect2 runs the annotation engine on the retained read likelihoods, not on
// the original pileup.  FORMAT/SB is a read-level BestAllele table, whereas
// FORMAT/F1R2 and FORMAT/F2R1 are fragment-level tables whose orientation is
// taken from the first read of the fragment.  Keep those two granularities
// separate here, just as StrandBiasBySample and OrientationBiasReadCounts do
// in GATK.
struct PerAlleleAnnotationValues {
    std::vector<std::uint8_t> base_qualities;
    std::vector<std::uint8_t> mapping_qualities;
    std::vector<std::int32_t> fragment_lengths;
    std::vector<std::int32_t> read_positions;
};

struct SomaticAnnotationCounts {
    std::vector<std::uint32_t> retained_read_depth;
    std::vector<std::uint32_t> reference_forward;
    std::vector<std::uint32_t> reference_reverse;
    std::vector<std::uint32_t> alternate_forward;
    std::vector<std::uint32_t> alternate_reverse;
    std::vector<std::uint32_t> fragment_f1r2_reference;
    std::vector<std::uint32_t> fragment_f2r1_reference;
    std::vector<std::uint32_t> fragment_f1r2_alternate;
    std::vector<std::uint32_t> fragment_f2r1_alternate;
    std::vector<PerAlleleAnnotationValues> reference_annotation_values;
    std::vector<PerAlleleAnnotationValues> alternate_annotation_values;
};

std::string read_name_for_record(const fastgatk::io::ReadBatch& reads,
                                 const std::size_t record) {
    if (reads.name_offsets.size() != reads.records() + 1 || record >= reads.records())
        return {};
    const auto begin = reads.name_offsets[record];
    const auto end = reads.name_offsets[record + 1];
    if (begin > end || end > reads.names.size()) return {};
    return std::string(reads.names.begin() + static_cast<std::ptrdiff_t>(begin),
                       reads.names.begin() + static_cast<std::ptrdiff_t>(end));
}

const fastgatk::calling::CandidateReadRealignment* candidate_read_realignment(
    const fastgatk::calling::Result& result, const std::size_t candidate,
    const std::size_t likelihood_record) {
    if (candidate >= result.candidates.size() ||
        candidate >= result.likelihood_candidate_read_context_ordinals.size() ||
        likelihood_record > std::numeric_limits<std::uint32_t>::max())
        return nullptr;
    const auto context_ordinal = result.likelihood_candidate_read_context_ordinals[candidate];
    if (context_ordinal == std::numeric_limits<std::uint32_t>::max()) return nullptr;
    const auto& contexts = result.likelihood_candidate_read_realignments;
    const auto iter = std::lower_bound(
        contexts.begin(), contexts.end(), std::tuple<std::uint32_t, std::uint32_t>{
            context_ordinal, static_cast<std::uint32_t>(likelihood_record)},
        [](const auto& context, const auto& key) {
            return std::tie(context.context_ordinal, context.source_record) < key;
        });
    if (iter == contexts.end() || iter->context_ordinal != context_ordinal ||
        iter->source_record != likelihood_record)
        return nullptr;
    return &*iter;
}

bool realigned_interval_for_candidate(
    const fastgatk::calling::Result& result, const std::size_t candidate,
    const std::size_t likelihood_record, std::int32_t& tid, std::int64_t& start,
    std::int64_t& end) {
    const auto* candidate_context = candidate_read_realignment(
        result, candidate, likelihood_record);
    if (candidate_context != nullptr && candidate_context->tid >= 0 &&
        candidate_context->start >= 0 && candidate_context->end > candidate_context->start) {
        tid = candidate_context->tid;
        start = candidate_context->start;
        end = candidate_context->end;
        return true;
    }
    // A read admitted to this candidate's AssemblyRegion must not borrow a
    // selected-haplotype projection from an overlapping region. If this
    // context could not produce a projection, GATK leaves the original read
    // alignment in place for its annotation lookup.
    if (candidate_context != nullptr) return false;
    if (likelihood_record < result.likelihood_read_realigned_tids.size() &&
        likelihood_record < result.likelihood_read_realigned_starts.size() &&
        likelihood_record < result.likelihood_read_realigned_ends.size() &&
        result.likelihood_read_realigned_tids[likelihood_record] >= 0 &&
        result.likelihood_read_realigned_starts[likelihood_record] >= 0 &&
        result.likelihood_read_realigned_ends[likelihood_record] >
            result.likelihood_read_realigned_starts[likelihood_record]) {
        tid = result.likelihood_read_realigned_tids[likelihood_record];
        start = result.likelihood_read_realigned_starts[likelihood_record];
        end = result.likelihood_read_realigned_ends[likelihood_record];
        return true;
    }
    return false;
}

// Coverage receives the precomputed read-likelihood evidence cardinality,
// after PairHMM's poorly-modelled-read filter but before the BestAllele
// informativeness reduction that produces FORMAT/DP and AD.  A flattened
// result delegates to the active-core AssemblyRegion that owns this EventMap
// allele; its public candidate scalar has already preserved that count.
std::uint32_t coverage_depth_for_candidate(const fastgatk::calling::Result& result,
                                           const std::size_t candidate) {
    if (candidate >= result.candidates.size()) return 0U;
    // A partitioned Mutect2 traversal intentionally flattens the public
    // Result after each AssemblyRegion has produced its own PairHMM matrix.
    // The flattened Result has no single read x allele matrix by design, so
    // resolve INFO/DP back to the AssemblyRegion that owns this EventMap
    // allele. This is the same active-core ownership rule used for the
    // regional somatic-posterior reduction; choosing the outer max-depth
    // candidate here loses the original AlleleLikelihoods evidence count.
    if (!result.assembly_region_likelihood_results.empty()) {
        const auto& event = result.candidates[candidate];
        std::uint32_t selected_depth = event.variant_depth;
        int selected_priority = -1;
        for (std::size_t part_index = 0;
             part_index < result.assembly_region_likelihood_results.size(); ++part_index) {
            const auto& part = result.assembly_region_likelihood_results[part_index];
            if (!part) continue;
            const auto local_candidate = candidate_index(*part, event);
            if (local_candidate == std::numeric_limits<std::size_t>::max()) continue;
            int priority = 1;
            for (const auto& region : part->calling_regions) {
                if (region.tid == event.tid && event.position >= region.active_start &&
                    event.position <= region.active_end) {
                    priority = 2;
                    break;
                }
            }
            if (priority < selected_priority) continue;
            const auto local_depth = coverage_depth_for_candidate(*part, local_candidate);
            if (priority > selected_priority ||
                (priority == selected_priority && selected_priority < 0)) {
                selected_depth = local_depth;
                selected_priority = priority;
            }
        }
        return selected_depth;
    }
    return result.candidates[candidate].variant_depth;
}

std::optional<std::uint32_t> read_offset_at_candidate(
    const fastgatk::calling::Result& result, const std::size_t candidate_index,
    const std::size_t likelihood_record,
    const fastgatk::io::ReadBatch& reads, const std::size_t record,
    const fastgatk::calling::AssemblyCandidate& candidate,
    const bool allow_deleted_reference_coordinate) {
    if (record >= reads.records() || record >= reads.offsets.size() - 1 ||
        record >= reads.tids.size() || reads.tids[record] != candidate.tid)
        return std::nullopt;
    // BaseQuality and ReadPosition first reject a VariantContext start that
    // lies outside the read's realigned CIGAR interval.  The sparse
    // projection can contain an endpoint observation for a terminal
    // insertion/deletion; that coordinate is not a read base at this SNV
    // position and must not leak into either annotation.  Intervals here are
    // zero-based half-open, whereas GATK's read end is one-based inclusive.
    std::int32_t realigned_tid = candidate.tid;
    std::int64_t realigned_start = -1;
    std::int64_t realigned_end = -1;
    if (realigned_interval_for_candidate(result, candidate_index, likelihood_record,
                                         realigned_tid, realigned_start, realigned_end) &&
        (realigned_tid != candidate.tid || candidate.position < realigned_start ||
         candidate.position >= realigned_end))
        return std::nullopt;
    const auto* candidate_context = candidate_read_realignment(
        result, candidate_index, likelihood_record);
    const auto candidate_realigned_count = std::min({
        result.likelihood_candidate_read_realigned_base_contexts.size(),
        result.likelihood_candidate_read_realigned_base_records.size(),
        result.likelihood_candidate_read_realigned_base_offsets.size(),
        result.likelihood_candidate_read_realigned_base_positions.size()});
    for (std::size_t index = 0; index < candidate_realigned_count; ++index) {
        if (candidate_context == nullptr ||
            result.likelihood_candidate_read_realigned_base_contexts[index] !=
                candidate_context->context_ordinal ||
            result.likelihood_candidate_read_realigned_base_records[index] != likelihood_record ||
            result.likelihood_candidate_read_realigned_base_positions[index] != candidate.position)
            continue;
        const bool concrete =
            index >= result.likelihood_candidate_read_realigned_base_is_concrete.size() ||
            result.likelihood_candidate_read_realigned_base_is_concrete[index] != 0;
        if (!allow_deleted_reference_coordinate && !concrete) continue;
        const auto begin = reads.offsets[record];
        const auto offset = result.likelihood_candidate_read_realigned_base_offsets[index];
        if (begin <= reads.qualities.size() && offset < reads.qualities.size() - begin)
            return offset;
    }
    // A candidate-scoped realignment is the exact GATKRead used by
    // BaseQualityRankSumTest and ReadPosition. If that realigned CIGAR has
    // no base at the locus (for example, it clips a terminal base that the
    // original BAM CIGAR still covers), GATK returns an empty lookup. Do not
    // resurrect the original-CIGAR base in that case. Only a failed/missing
    // candidate realignment retains the original read alignment.
    if (candidate_context != nullptr && candidate_context->tid >= 0 &&
        candidate_context->start >= 0 && candidate_context->end > candidate_context->start)
        return std::nullopt;
    const auto begin = reads.offsets[record];
    const auto end = reads.offsets[record + 1];
    if (begin > end || end > reads.qualities.size()) return std::nullopt;
    for (std::size_t offset = begin; offset < end; ++offset) {
        fastgatk::io::ReadProjection projection;
        if (fastgatk::io::project_read_offset(reads, record, offset - begin, projection) &&
            projection.reference_position == candidate.position)
            return static_cast<std::uint32_t>(offset - begin);
    }
    return std::nullopt;
}

std::optional<std::uint8_t> base_quality_at_candidate(
    const fastgatk::calling::Result& result, const std::size_t candidate_index,
    const std::size_t likelihood_record,
    const fastgatk::io::ReadBatch& reads, const std::size_t record,
    const fastgatk::calling::AssemblyCandidate& candidate) {
    const auto offset = read_offset_at_candidate(
        result, candidate_index, likelihood_record, reads, record, candidate, false);
    if (!offset.has_value() || record + 1 >= reads.offsets.size()) return std::nullopt;
    const auto begin = reads.offsets[record];
    if (begin > reads.qualities.size() || *offset >= reads.qualities.size() - begin)
        return std::nullopt;
    // PairHMM rows are formed from the overlap-corrected Host batch.  GATK's
    // Mutect2 annotations observe those same callRegion() qualities after
    // best-haplotype realignment, whereas `record` is only the immutable
    // outer-BAM provenance used to recover flags/MAPQ.  Prefer the compact
    // row payload whenever it is present; retain the raw path for Results
    // created by older/direct callers.
    if (likelihood_record + 1 < result.likelihood_read_quality_offsets.size()) {
        const auto quality_begin = result.likelihood_read_quality_offsets[likelihood_record];
        const auto quality_end = result.likelihood_read_quality_offsets[likelihood_record + 1];
        if (quality_begin <= quality_end && quality_end <= result.likelihood_read_qualities.size() &&
            *offset < quality_end - quality_begin)
            return result.likelihood_read_qualities[quality_begin + *offset];
    }
    return reads.qualities[begin + *offset];
}

std::optional<SomaticAnnotationCounts> somatic_annotation_counts_local(
    const fastgatk::calling::Result& result,
    const fastgatk::io::ReadBatch& reads,
    const std::vector<fastgatk::calling::AssemblyCandidate>* emitted_candidates = nullptr) {
    const auto grouped = group_somatic_evidence(result);
    if (!grouped.has_value()) return std::nullopt;
    const bool debug_annotations =
        std::getenv("FASTGATK_DEBUG_SOMATIC_ANNOTATIONS") != nullptr;
    const auto candidate_count = result.candidates.size();
    const auto read_count = result.likelihood_read_names.size();
    const auto fragment_count = grouped->read_count;
    if (candidate_count == 0 || read_count == 0 || fragment_count == 0 ||
        result.reference_read_likelihoods.size() != candidate_count ||
        result.allele_read_likelihoods.size() != candidate_count)
        return std::nullopt;

    SomaticAnnotationCounts output;
    output.retained_read_depth.assign(candidate_count, 0U);
    output.reference_forward.assign(candidate_count, 0U);
    output.reference_reverse.assign(candidate_count, 0U);
    output.alternate_forward.assign(candidate_count, 0U);
    output.alternate_reverse.assign(candidate_count, 0U);
    output.fragment_f1r2_reference.assign(candidate_count, 0U);
    output.fragment_f2r1_reference.assign(candidate_count, 0U);
    output.fragment_f1r2_alternate.assign(candidate_count, 0U);
    output.fragment_f2r1_alternate.assign(candidate_count, 0U);
    output.reference_annotation_values.resize(candidate_count);
    output.alternate_annotation_values.resize(candidate_count);

    // callMutations annotates the EventMap-marginalized AlleleLikelihoods,
    // not the broader Host candidate inventory.  The latter retains raw
    // CIGAR/pileup observations for diagnostics and recovery bookkeeping;
    // when one of those raw siblings shares a likelihood context with a
    // concrete EventMap allele, context equality alone is not sufficient to
    // admit it to BestAllele/F1R2.  Preserve the same PairHMM-time ownership
    // boundary used by the somatic evidence reducer below.
    const bool has_pairhmm_event_map_ownership =
        result.likelihood_candidate_event_map_owned.size() == candidate_count;
    const auto event_map_owns_candidate = [&](const std::size_t candidate) {
        if (candidate >= candidate_count) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position &&
                    event.reference == allele_ref(item) && event.alternate == allele_alt(item);
            });
        });
    };
    const auto event_map_has_locus = [&](const std::size_t candidate) {
        if (candidate >= candidate_count) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position;
            });
        });
    };
    const auto belongs_to_annotation_allele_list = [&](const std::size_t candidate) {
        const auto& item = result.candidates[candidate];
        const bool recovery_without_event_map_locus =
            item.somatic_low_quality_assembly_recovery &&
            (result.graph_haplotype_event_maps.empty() || !event_map_has_locus(candidate));
        if (recovery_without_event_map_locus) return true;
        if (has_pairhmm_event_map_ownership)
            return result.likelihood_candidate_event_map_owned[candidate] != 0;
        if (!result.graph_haplotype_event_maps.empty())
            return event_map_owns_candidate(candidate);
        return true;
    };

    using FragmentKey = std::tuple<std::int32_t, std::int32_t, std::string>;
    std::map<FragmentKey, std::vector<std::size_t>> candidate_groups;
    if (emitted_candidates != nullptr) {
        // Depth/strand/rank annotations are applied to the final
        // VariantContext after SomaticGenotypingEngine has removed
        // non-emitted EventMap siblings.  Restrict the BestAllele table to
        // that concrete final list before any of its consumers run.
        for (const auto& emitted : *emitted_candidates) {
            const auto candidate = candidate_index(result, emitted);
            if (candidate == std::numeric_limits<std::size_t>::max() ||
                candidate >= candidate_count)
                return std::nullopt;
            const auto& item = result.candidates[candidate];
            candidate_groups[FragmentKey{item.tid, item.position, allele_ref(item)}].push_back(candidate);
        }
    } else {
        for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
            if (!belongs_to_annotation_allele_list(candidate)) continue;
            const auto& item = result.candidates[candidate];
            candidate_groups[FragmentKey{item.tid, item.position, allele_ref(item)}].push_back(candidate);
        }
    }

    // The PairHMM pre-grouping path preserves one source-record ordinal per
    // fragment. Use it as the authority for reattaching individual likelihood
    // rows to fragments: an earlier all-read name scan can include reads that
    // were filtered before PairHMM and therefore no longer exist in GATK's
    // annotation evidence collection.
    std::vector<std::size_t> fragment_for_read(
        read_count, std::numeric_limits<std::size_t>::max());
    const auto& representative = result.somatic_group_representative_sources;
    if (representative.size() == fragment_count) {
        std::map<std::string, std::size_t> fragment_by_name;
        for (std::size_t fragment = 0; fragment < fragment_count; ++fragment) {
            const auto source = static_cast<std::size_t>(representative[fragment]);
            if (source >= read_count || result.likelihood_read_names[source].empty()) continue;
            fragment_by_name.emplace(result.likelihood_read_names[source], fragment);
        }
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto found = fragment_by_name.find(result.likelihood_read_names[read]);
            if (found != fragment_by_name.end()) fragment_for_read[read] = found->second;
        }
    } else {
        std::map<std::string, std::size_t> fragments_by_name;
        std::size_t discovered_fragments = 0;
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto& name = result.likelihood_read_names[read];
            if (name.empty()) {
                fragment_for_read[read] = discovered_fragments++;
                continue;
            }
            const auto [iter, inserted] = fragments_by_name.emplace(name, discovered_fragments);
            if (inserted) ++discovered_fragments;
            fragment_for_read[read] = iter->second;
        }
        if (discovered_fragments != fragment_count) return std::nullopt;
    }

    // Result coordinates are post-read-filter, as are the likelihood rows.
    // Match them back to the decoded read only to obtain immutable SAM flags,
    // MAPQ and base quality for the source annotation predicates.
    using ReadKey = std::tuple<std::string, std::int32_t, std::int64_t, std::int64_t>;
    std::map<ReadKey, std::size_t> input_record_by_key;
    for (std::size_t record = 0; record < reads.records(); ++record) {
        const auto tid = record < reads.tids.size() ? reads.tids[record] : -1;
        const auto start = record < reads.positions.size() ?
            static_cast<std::int64_t>(reads.positions[record]) : -1;
        input_record_by_key.emplace(ReadKey{read_name_for_record(reads, record), tid, start,
                                             fastgatk::io::reference_end(reads, record)}, record);
    }
    std::vector<std::optional<std::size_t>> input_record_for_read(read_count);
    const bool source_records_available =
        result.likelihood_read_source_records.size() == read_count;
    const bool coordinates_available = result.likelihood_read_tids.size() == read_count &&
        result.likelihood_read_starts.size() == read_count &&
        result.likelihood_read_ends.size() == read_count;
    for (std::size_t read = 0; read < read_count; ++read) {
        // A PairHMM result created from an AssemblyRegion carries the outer
        // decoded-BAM ordinal explicitly.  This is authoritative: a
        // name/start/end map is not injective for duplicate or supplementary
        // records and can silently flip strand/orientation annotations.
        if (source_records_available) {
            const auto record = static_cast<std::size_t>(
                result.likelihood_read_source_records[read]);
            if (record < reads.records()) input_record_for_read[read] = record;
            continue;
        }
        // Compatibility path for Result objects produced by older callers.
        if (coordinates_available) {
            const auto found = input_record_by_key.find(ReadKey{
                result.likelihood_read_names[read], result.likelihood_read_tids[read],
                result.likelihood_read_starts[read], result.likelihood_read_ends[read]});
            if (found != input_record_by_key.end()) input_record_for_read[read] = found->second;
        }
        if (!input_record_for_read[read].has_value() && read < reads.records())
            input_record_for_read[read] = read;
    }

    constexpr double kMissing = -1.0e299;
    const auto best_allele = [&](const std::vector<double>& values)
        -> std::optional<std::size_t> {
        return gatk_mutect_best_allele(values, kMissing);
    };

    for (const auto& [key, candidate_indices] : candidate_groups) {
        (void)key;
        if (candidate_indices.empty()) continue;
        // `Result::candidates` retains raw pileup observations in addition
        // to the EventMap alleles that entered the AssemblyRegion's
        // likelihood collection.  At a multiallelic locus, an unmodeled raw
        // sibling has no context ordinal.  It must not make the modeled
        // siblings (and every subsequent annotation) invalid: GATK builds
        // BestAllele from the EventMap allele list only.
        std::vector<std::size_t> indices(candidate_indices);
        const bool candidate_contexts_available =
            result.likelihood_candidate_read_context_ordinals.size() == candidate_count;
        if (candidate_contexts_available) {
            constexpr auto kMissingContext = std::numeric_limits<std::uint32_t>::max();
            auto context = kMissingContext;
            for (const auto candidate : indices) {
                const auto candidate_context =
                    result.likelihood_candidate_read_context_ordinals[candidate];
                if (candidate_context != kMissingContext) {
                    context = candidate_context;
                    break;
                }
            }
            if (context == kMissingContext) continue;
            indices.erase(std::remove_if(indices.begin(), indices.end(),
                [&](const auto candidate) {
                    return result.likelihood_candidate_read_context_ordinals[candidate] != context;
                }), indices.end());
            if (indices.empty()) continue;
        }
        // The sparse PairHMM result can contain the same read name in more
        // than one AssemblyRegion context.  Reconstruct the Fragment lookup
        // for this callMutations() likelihood collection, not globally by
        // name: OrientationBiasReadCounts consumes that exact Fragment and
        // its first read.  The legacy name-only table above remains only for
        // Results produced before context provenance was available.
        auto candidate_fragment_for_read = fragment_for_read;
        const auto& fragment_contexts = result.somatic_group_context_ordinals;
        if (fragment_contexts.size() == fragment_count &&
            indices.front() < result.likelihood_candidate_read_context_ordinals.size()) {
            const auto context = result.likelihood_candidate_read_context_ordinals[indices.front()];
            if (context != std::numeric_limits<std::uint32_t>::max()) {
                for (const auto candidate : indices) {
                    if (candidate >= result.likelihood_candidate_read_context_ordinals.size() ||
                        result.likelihood_candidate_read_context_ordinals[candidate] != context)
                        return std::nullopt;
                }
                using ScopedFragmentKey = std::tuple<std::uint32_t, std::uint32_t, std::string>;
                std::map<ScopedFragmentKey, std::size_t> fragments_in_context;
                for (std::size_t fragment = 0; fragment < fragment_count; ++fragment) {
                    if (fragment_contexts[fragment] != context) continue;
                    const auto source = static_cast<std::size_t>(representative[fragment]);
                    if (source >= read_count || result.likelihood_read_names[source].empty()) continue;
                    const auto sample = source < reads.sample_ids.size()
                        ? reads.sample_ids[source] : 0U;
                    fragments_in_context.emplace(ScopedFragmentKey{
                        context, sample, result.likelihood_read_names[source]}, fragment);
                }
                candidate_fragment_for_read.assign(
                    read_count, std::numeric_limits<std::size_t>::max());
                for (std::size_t read = 0; read < read_count; ++read) {
                    if (!input_record_for_read[read].has_value()) continue;
                    const auto record = *input_record_for_read[read];
                    const auto sample = record < reads.sample_ids.size()
                        ? reads.sample_ids[record] : 0U;
                    const auto found = fragments_in_context.find(ScopedFragmentKey{
                        context, sample, result.likelihood_read_names[read]});
                    if (found != fragments_in_context.end())
                        candidate_fragment_for_read[read] = found->second;
                }
            }
        }
        std::vector<std::uint8_t> fragment_overlaps(fragment_count, 0U);
        if (representative.size() != fragment_count || !coordinates_available) return std::nullopt;
        const auto reference_length = allele_ref(result.candidates[indices.front()]).size();
        const auto event_start = std::max<std::int64_t>(
            0, static_cast<std::int64_t>(result.candidates[indices.front()].position) - 2);
        const auto event_end = std::max<std::int64_t>(event_start + 1,
            static_cast<std::int64_t>(result.candidates[indices.front()].position) +
            static_cast<std::int64_t>(reference_length == 0 ? 1 : reference_length) + 2);
        for (std::size_t read = 0; read < read_count; ++read) {
            const auto fragment = candidate_fragment_for_read[read];
            if (fragment >= fragment_overlaps.size()) continue;
            std::int32_t tid = result.likelihood_read_tids[read];
            std::int64_t start = result.likelihood_read_starts[read];
            std::int64_t end = result.likelihood_read_ends[read];
            realigned_interval_for_candidate(
                result, indices.front(), read, tid, start, end);
            if (tid == result.candidates[indices.front()].tid && start >= 0 && end > start &&
                start < event_end && end > event_start)
                fragment_overlaps[fragment] = 1U;
        }

        // FragmentDepthPerAlleleBySample and OrientationBiasReadCounts share
        // this fragment BestAllele reduction. The latter additionally gates
        // its first read by usable MAPQ and a base quality of at least 20.
        std::set<std::uint32_t> active_fragments;
        for (const auto candidate : indices) {
            for (const auto fragment : grouped->groups_for_candidate(candidate))
                active_fragments.insert(fragment);
        }
        for (const auto fragment_id : active_fragments) {
            const auto fragment = static_cast<std::size_t>(fragment_id);
            if (fragment >= fragment_count) continue;
            if (fragment_overlaps[fragment] == 0) continue;
            // Annotations receive the final VariantContext allele list:
            // its reference plus the emitted concrete ALTs.  A temporary
            // spanning-deletion state from marginalization is not an allele
            // of that VariantContext, so it cannot compete in this fragment
            // BestAllele/F1R2/F2R1 reduction.
            std::vector<double> values(indices.size() + 1, kMissing);
            // The locus has one shared REF allele.  Raw pileup siblings can
            // remain in Result even when only another concrete ALT entered
            // the EventMap, so the first candidate is not necessarily the
            // one whose read/fragment row carries the REF marginalization.
            // GATK's allele matrix has exactly one REF row: reconstruct it
            // as the maximum available REF likelihood across this locus.
            for (const auto candidate : indices)
                values[0] = std::max(values[0], grouped->reference_at(candidate, fragment));
            for (std::size_t alternate = 0; alternate < indices.size(); ++alternate)
                values[alternate + 1] = grouped->alternate_at(indices[alternate], fragment);
            const auto best = best_allele(values);
            if (!best.has_value() || fragment >= representative.size()) continue;
            // OrientationBiasReadCounts calls Fragment.getFirstRead(), then
            // applies its MAPQ/base-quality gates to that exact GATKRead.
            // Fragment is built in AssemblyRegion order, so its retained
            // representative already carries the source-defined first-read
            // identity. Re-selecting a mate from the local likelihood rows
            // can invert F1R2/F2R1 in overlapping regions.
            const std::size_t orientation_source =
                static_cast<std::size_t>(representative[fragment]);
            if (orientation_source >= input_record_for_read.size() ||
                !input_record_for_read[orientation_source].has_value())
                continue;
            const auto orientation_record = *input_record_for_read[orientation_source];
            const auto mapq = orientation_record < reads.mapq.size()
                ? reads.mapq[orientation_record] : 0U;
            const auto base_quality = base_quality_at_candidate(
                result, indices.front(), orientation_source, reads, orientation_record,
                result.candidates[indices.front()]);
            if (debug_annotations) {
                const auto flags = orientation_record < reads.flags.size()
                    ? reads.flags[orientation_record] : 0U;
                std::int32_t realigned_tid = orientation_source < result.likelihood_read_tids.size()
                    ? result.likelihood_read_tids[orientation_source] : -1;
                std::int64_t realigned_start = orientation_source < result.likelihood_read_starts.size()
                    ? result.likelihood_read_starts[orientation_source] : -1;
                std::int64_t realigned_end = orientation_source < result.likelihood_read_ends.size()
                    ? result.likelihood_read_ends[orientation_source] : -1;
                const bool has_candidate_realignment = realigned_interval_for_candidate(
                    result, indices.front(), orientation_source, realigned_tid,
                    realigned_start, realigned_end);
                std::cerr << "[FASTGATK_SOMATIC_FRAGMENT]"
                          << " position=" << result.candidates[indices.front()].position
                          << " fragment=" << fragment
                          << " source=" << orientation_source
                          << " record=" << orientation_record
                          << " name=" << read_name_for_record(reads, orientation_record)
                          << " ref=" << values[0]
                          << " alt=" << (values.size() > 1 ? values[1] : kMissing)
                          << " best=" << *best
                          << " flags=" << flags
                          << " mapq=" << static_cast<unsigned>(mapq)
                          << " baseq=" << (base_quality.has_value()
                              ? std::to_string(*base_quality) : std::string{"missing"})
                          << " overlaps=" << static_cast<unsigned>(fragment_overlaps[fragment])
                          << " realigned=" << (has_candidate_realignment ? 1 : 0)
                          << " realigned_tid=" << realigned_tid
                          << " realigned_start=" << realigned_start
                          << " realigned_end=" << realigned_end
                          << '\n';
            }
            if (mapq == 0U || mapq == 255U || !base_quality.has_value() || *base_quality < 20U)
                continue;
            const auto flags = orientation_record < reads.flags.size()
                ? reads.flags[orientation_record] : 0U;
            const bool reverse = (flags & 0x10U) != 0;
            const bool first_of_pair = (flags & 0x40U) != 0;
            const bool f2r1 = reverse == first_of_pair;
            if (*best == 0) {
                for (const auto candidate : indices) {
                    ++(f2r1 ? output.fragment_f2r1_reference[candidate]
                            : output.fragment_f1r2_reference[candidate]);
                }
            } else if (*best - 1 < indices.size()) {
                const auto candidate = indices[*best - 1];
                ++(f2r1 ? output.fragment_f2r1_alternate[candidate]
                        : output.fragment_f1r2_alternate[candidate]);
            }
        }

        // StrandBiasBySample counts individual retained reads. Coverage uses
        // the same retained population but does not require an informative
        // BestAllele, so retain its count before the strict natural-log 0.2
        // margin (converted above for native log10 likelihoods).
        for (std::size_t read = 0; read < read_count; ++read) {
            // Unlike FAD/F1R2/F2R1, read likelihood annotations do not keep
            // a non-overlapping mate merely because its fragment overlaps.
            // SomaticGenotypingEngine calls retainEvidence() separately on
            // logReadAlleleLikelihoods with the expanded event interval.
            std::int32_t tid = result.likelihood_read_tids[read];
            std::int64_t start = result.likelihood_read_starts[read];
            std::int64_t end = result.likelihood_read_ends[read];
            realigned_interval_for_candidate(
                result, indices.front(), read, tid, start, end);
            if (tid != result.candidates[indices.front()].tid || start < 0 || end <= start ||
                start >= event_end || end <= event_start)
                continue;
            if (read < result.annotation_read_qualified.size() &&
                result.annotation_read_qualified[read] == 0)
                continue;
            // SomaticGenotypingEngine maps `*` while constructing its local
            // EventMap matrix, then drops it when it trims to the REF plus
            // emitted concrete ALT alleles for read annotations. Keep the
            // Kokkos marginalisation's true REF row, but do not resurrect
            // `*` as a third BestAllele here.
            std::vector<double> values(indices.size() + 1, kMissing);
            // A concrete ALT that never entered EventMap can have an empty
            // copied row even though a sibling EventMap allele has the same
            // locus-wide REF likelihood.  Use the one shared REF row's
            // maximum, as AlleleLikelihoods.marginalize() does, rather than
            // making the arbitrary first raw sibling the REF owner.
            std::size_t reference_likelihood_candidate = indices.front();
            for (const auto candidate : indices) {
                if (result.reference_read_likelihoods[candidate].size() != read_count)
                    return std::nullopt;
                const auto reference = result.reference_read_likelihoods[candidate][read];
                if (reference > values[0]) {
                    values[0] = reference;
                    reference_likelihood_candidate = candidate;
                }
            }
            bool has_likelihood = std::isfinite(values[0]) && values[0] > kMissing;
            for (std::size_t alternate = 0; alternate < indices.size(); ++alternate) {
                if (result.allele_read_likelihoods[indices[alternate]].size() != read_count)
                    return std::nullopt;
                values[alternate + 1] = result.allele_read_likelihoods[indices[alternate]][read];
                has_likelihood = has_likelihood ||
                    (std::isfinite(values[alternate + 1]) && values[alternate + 1] > kMissing);
            }
            if (!has_likelihood) continue;
            for (const auto candidate : indices) ++output.retained_read_depth[candidate];
            const auto best = best_allele(values);
            if (debug_annotations) {
                const auto record = read < input_record_for_read.size() &&
                    input_record_for_read[read].has_value()
                    ? std::to_string(*input_record_for_read[read]) : std::string{"missing"};
                std::cerr << "[FASTGATK_SOMATIC_READ]"
                          << " position=" << result.candidates[indices.front()].position
                          << " source=" << read
                          << " record=" << record
                          << " name=" << (input_record_for_read[read].has_value()
                              ? read_name_for_record(reads, *input_record_for_read[read])
                              : std::string{})
                          << " tid=" << tid
                          << " start=" << start
                          << " end=" << end
                          << " ref=" << values[0]
                          << " alt=" << (values.size() > 1 ? values[1] : kMissing)
                          << " best=" << (best.has_value() ? std::to_string(*best)
                              : std::string{"uninformative"})
                          << " qualified=" << (read < result.annotation_read_qualified.size()
                              ? static_cast<unsigned>(result.annotation_read_qualified[read]) : 1U)
                          << '\n';
            }
            if (!best.has_value() || read >= input_record_for_read.size() ||
                !input_record_for_read[read].has_value())
                continue;
            const auto flags = *input_record_for_read[read] < reads.flags.size()
                ? reads.flags[*input_record_for_read[read]] : 0U;
            const bool reverse = (flags & 0x10U) != 0;
            if (*best == 0) {
                for (const auto candidate : indices) {
                    ++(reverse ? output.reference_reverse[candidate]
                                : output.reference_forward[candidate]);
                }
            } else if (*best - 1 < indices.size()) {
                const auto candidate = indices[*best - 1];
                ++(reverse ? output.alternate_reverse[candidate]
                            : output.alternate_forward[candidate]);
            }

            // PerAlleleAnnotation only consumes informative reads whose
            // mapping quality is available. Unlike strand bias, it then
            // gathers a metric-specific value (base quality/read position
            // may be absent even when MAPQ and fragment length are valid).
            const auto record = *input_record_for_read[read];
            const auto mapq = record < reads.mapq.size() ? reads.mapq[record] : 0U;
            if (mapq == 0U || mapq == 255U) continue;
            std::vector<PerAlleleAnnotationValues*> annotation_values;
            std::size_t annotation_candidate = reference_likelihood_candidate;
            if (*best == 0) {
                // The one shared REF AlleleLikelihoods row is visible from
                // every concrete output ALT at the locus. Materialize its
                // annotation values for each candidate so a later VCF call
                // whose only emitted ALT is not the first raw pileup sibling
                // still observes the same REF distribution.
                annotation_values.reserve(indices.size());
                for (const auto candidate : indices)
                    annotation_values.push_back(&output.reference_annotation_values[candidate]);
            } else if (*best - 1 < indices.size()) {
                annotation_candidate = indices[*best - 1];
                annotation_values.push_back(&output.alternate_annotation_values[annotation_candidate]);
            }
            if (annotation_values.empty()) continue;
            const auto base_quality = base_quality_at_candidate(
                result, annotation_candidate, read, reads, record,
                result.candidates[annotation_candidate]);
            std::optional<std::int32_t> read_position;
            if (const auto read_offset = read_offset_at_candidate(
                    result, annotation_candidate, read, reads, record,
                    result.candidates[annotation_candidate], true);
                read_offset.has_value() && record + 1 < reads.offsets.size()) {
                const auto begin = reads.offsets[record];
                const auto end = reads.offsets[record + 1];
                if (end > begin && *read_offset < end - begin) {
                    const auto left = *read_offset;
                    const auto right = static_cast<std::uint32_t>(end - begin - 1) - left;
                    read_position = static_cast<std::int32_t>(std::min(left, right));
                }
            }
            for (auto* annotation_values_for_allele : annotation_values) {
                annotation_values_for_allele->mapping_qualities.push_back(mapq);
                if (record < reads.template_lengths.size()) {
                    const auto length = static_cast<std::int64_t>(reads.template_lengths[record]);
                    annotation_values_for_allele->fragment_lengths.push_back(static_cast<std::int32_t>(
                        std::min<std::int64_t>(std::numeric_limits<std::int32_t>::max(),
                                               std::llabs(length))));
                }
                if (base_quality.has_value())
                    annotation_values_for_allele->base_qualities.push_back(*base_quality);
                if (read_position.has_value())
                    annotation_values_for_allele->read_positions.push_back(*read_position);
            }
        }
    }
    return output;
}

std::optional<SomaticAnnotationCounts> somatic_annotation_counts(
    const fastgatk::calling::Result& result,
    const fastgatk::io::ReadBatch& reads,
    const std::vector<fastgatk::calling::AssemblyCandidate>* emitted_candidates = nullptr) {
    if (result.assembly_region_likelihood_results.empty())
        return somatic_annotation_counts_local(result, reads, emitted_candidates);

    if (emitted_candidates != nullptr) {
        if (emitted_candidates->empty()) return std::nullopt;
        // The output locus is owned by one active AssemblyRegion.  Resolve it
        // once, calculate the final-allele annotation table in its native
        // likelihood-row coordinate, then project just those counts back to
        // the flattened Result used by the VCF writer.
        const auto& first = emitted_candidates->front();
        const fastgatk::calling::Result* owner = nullptr;
        int owner_priority = -1;
        for (const auto& part_pointer : result.assembly_region_likelihood_results) {
            if (!part_pointer) continue;
            const auto& part = *part_pointer;
            const auto has_all = std::all_of(emitted_candidates->begin(), emitted_candidates->end(),
                [&](const auto& emitted) {
                    return candidate_index(part, emitted) !=
                        std::numeric_limits<std::size_t>::max();
                });
            if (!has_all) continue;
            int priority = 1;
            for (const auto& region : part.calling_regions) {
                if (region.tid == first.tid && first.position >= region.active_start &&
                    first.position <= region.active_end) {
                    priority = 2;
                    break;
                }
            }
            if (priority > owner_priority) {
                owner = &part;
                owner_priority = priority;
            }
        }
        if (owner == nullptr) return std::nullopt;
        const auto local = somatic_annotation_counts_local(*owner, reads, emitted_candidates);
        if (!local.has_value()) return std::nullopt;
        const auto count = result.candidates.size();
        SomaticAnnotationCounts output;
        output.retained_read_depth.assign(count, 0U);
        output.reference_forward.assign(count, 0U);
        output.reference_reverse.assign(count, 0U);
        output.alternate_forward.assign(count, 0U);
        output.alternate_reverse.assign(count, 0U);
        output.fragment_f1r2_reference.assign(count, 0U);
        output.fragment_f2r1_reference.assign(count, 0U);
        output.fragment_f1r2_alternate.assign(count, 0U);
        output.fragment_f2r1_alternate.assign(count, 0U);
        output.reference_annotation_values.resize(count);
        output.alternate_annotation_values.resize(count);
        for (const auto& emitted : *emitted_candidates) {
            const auto local_index = candidate_index(*owner, emitted);
            const auto merged_index = candidate_index(result, emitted);
            if (local_index == std::numeric_limits<std::size_t>::max() ||
                merged_index == std::numeric_limits<std::size_t>::max() ||
                local_index >= local->retained_read_depth.size() || merged_index >= count)
                return std::nullopt;
            const auto copy_count = [&](auto& destination, const auto& source) {
                if (local_index < source.size()) destination[merged_index] = source[local_index];
            };
            copy_count(output.retained_read_depth, local->retained_read_depth);
            copy_count(output.reference_forward, local->reference_forward);
            copy_count(output.reference_reverse, local->reference_reverse);
            copy_count(output.alternate_forward, local->alternate_forward);
            copy_count(output.alternate_reverse, local->alternate_reverse);
            copy_count(output.fragment_f1r2_reference, local->fragment_f1r2_reference);
            copy_count(output.fragment_f2r1_reference, local->fragment_f2r1_reference);
            copy_count(output.fragment_f1r2_alternate, local->fragment_f1r2_alternate);
            copy_count(output.fragment_f2r1_alternate, local->fragment_f2r1_alternate);
            if (local_index < local->reference_annotation_values.size())
                output.reference_annotation_values[merged_index] =
                    local->reference_annotation_values[local_index];
            if (local_index < local->alternate_annotation_values.size())
                output.alternate_annotation_values[merged_index] =
                    local->alternate_annotation_values[local_index];
        }
        return output;
    }

    // AssemblyRegion likelihood rows have region-local read ordinals.  GATK
    // evaluates annotations before CallAssemblyRegions merges its output, so
    // reconstruct the flattened annotation vector by taking the candidate's
    // owning active-region result.  A halo result is retained only when the
    // owning region did not materialize that exact allele.
    const auto count = result.candidates.size();
    SomaticAnnotationCounts output;
    output.retained_read_depth.assign(count, 0U);
    output.reference_forward.assign(count, 0U);
    output.reference_reverse.assign(count, 0U);
    output.alternate_forward.assign(count, 0U);
    output.alternate_reverse.assign(count, 0U);
    output.fragment_f1r2_reference.assign(count, 0U);
    output.fragment_f2r1_reference.assign(count, 0U);
    output.fragment_f1r2_alternate.assign(count, 0U);
    output.fragment_f2r1_alternate.assign(count, 0U);
    output.reference_annotation_values.resize(count);
    output.alternate_annotation_values.resize(count);
    std::vector<int> ownership(count, -1);
    bool recovered = false;

    const auto find_merged_candidate = [&](const auto& candidate) {
        return std::find_if(result.candidates.begin(), result.candidates.end(),
            [&](const auto& merged) {
                return merged.tid == candidate.tid &&
                    merged.position == candidate.position &&
                    allele_ref(merged) == allele_ref(candidate) &&
                    allele_alt(merged) == allele_alt(candidate);
            });
    };
    for (const auto& part_pointer : result.assembly_region_likelihood_results) {
        if (!part_pointer) continue;
        const auto& part = *part_pointer;
        const auto local = somatic_annotation_counts(part, reads);
        if (!local.has_value()) continue;
        for (std::size_t local_index = 0; local_index < part.candidates.size(); ++local_index) {
            const auto merged = find_merged_candidate(part.candidates[local_index]);
            if (merged == result.candidates.end()) continue;
            const auto merged_index = static_cast<std::size_t>(merged - result.candidates.begin());
            int priority = 1;
            if (part.calling_regions.size() == 1) {
                const auto& region = part.calling_regions.front();
                if (part.candidates[local_index].tid == region.tid &&
                    part.candidates[local_index].position >= region.active_start &&
                    part.candidates[local_index].position <= region.active_end)
                    priority = 2;
            }
            if (priority < ownership[merged_index]) continue;
            ownership[merged_index] = priority;
            recovered = true;
            const auto copy_count = [&](auto& destination, const auto& source) {
                if (local_index < source.size()) destination[merged_index] = source[local_index];
            };
            copy_count(output.retained_read_depth, local->retained_read_depth);
            copy_count(output.reference_forward, local->reference_forward);
            copy_count(output.reference_reverse, local->reference_reverse);
            copy_count(output.alternate_forward, local->alternate_forward);
            copy_count(output.alternate_reverse, local->alternate_reverse);
            copy_count(output.fragment_f1r2_reference, local->fragment_f1r2_reference);
            copy_count(output.fragment_f2r1_reference, local->fragment_f2r1_reference);
            copy_count(output.fragment_f1r2_alternate, local->fragment_f1r2_alternate);
            copy_count(output.fragment_f2r1_alternate, local->fragment_f2r1_alternate);
            if (local_index < local->reference_annotation_values.size())
                output.reference_annotation_values[merged_index] =
                    local->reference_annotation_values[local_index];
            if (local_index < local->alternate_annotation_values.size())
                output.alternate_annotation_values[merged_index] =
                    local->alternate_annotation_values[local_index];
        }
    }
    return recovered ? std::optional<SomaticAnnotationCounts>{std::move(output)} : std::nullopt;
}

// DepthPerAlleleBySample is applied by GATK's annotation engine after
// SomaticGenotypingEngine has removed non-emitted EventMap alleles from the
// VariantContext.  Re-evaluate the read-level BestAllele table against that
// exact final REF + ALT list.  This is deliberately separate from FAD: the
// latter remains the fragment-level annotation calculated above.
struct EmittedSomaticReadDepth {
    std::uint32_t depth = 0U;
    // Allele-list order: shared REF followed by the caller's emitted ALTs.
    std::vector<std::uint32_t> allele_depths;
};

std::optional<EmittedSomaticReadDepth> emitted_somatic_read_depth_local(
    const fastgatk::calling::Result& result,
    const fastgatk::io::ReadBatch& reads,
    const std::vector<fastgatk::calling::AssemblyCandidate>& emitted_candidates) {
    (void)reads;
    if (emitted_candidates.empty()) return std::nullopt;
    const auto read_count = result.likelihood_read_names.size();
    if (read_count == 0 || result.likelihood_read_tids.size() != read_count ||
        result.likelihood_read_starts.size() != read_count ||
        result.likelihood_read_ends.size() != read_count)
        return std::nullopt;

    std::vector<std::size_t> indices;
    indices.reserve(emitted_candidates.size());
    for (const auto& emitted : emitted_candidates) {
        const auto index = candidate_index(result, emitted);
        if (index == std::numeric_limits<std::size_t>::max() ||
            index >= result.candidates.size() ||
            index >= result.reference_read_likelihoods.size() ||
            index >= result.allele_read_likelihoods.size() ||
            result.reference_read_likelihoods[index].size() != read_count ||
            result.allele_read_likelihoods[index].size() != read_count)
            return std::nullopt;
        indices.push_back(index);
    }
    const auto& first = result.candidates[indices.front()];
    for (const auto index : indices) {
        const auto& candidate = result.candidates[index];
        if (candidate.tid != first.tid || candidate.position != first.position ||
            allele_ref(candidate) != allele_ref(first))
            return std::nullopt;
    }

    constexpr double kMissing = -1.0e299;
    const auto reference_length = allele_ref(first).size();
    const auto event_start = std::max<std::int64_t>(
        0, static_cast<std::int64_t>(first.position) - 2);
    const auto event_end = std::max<std::int64_t>(event_start + 1,
        static_cast<std::int64_t>(first.position) +
        static_cast<std::int64_t>(reference_length == 0 ? 1 : reference_length) + 2);

    if (read_count > std::numeric_limits<std::uint32_t>::max()) return std::nullopt;
    std::vector<fastgatk::pairhmm::ReadAlleleBestRequest> requests;
    requests.reserve(read_count * (indices.size() + 1));
    std::vector<double> reference_by_row;
    reference_by_row.reserve(read_count);
    std::size_t row_count = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        std::int32_t tid = result.likelihood_read_tids[read];
        std::int64_t start = result.likelihood_read_starts[read];
        std::int64_t end = result.likelihood_read_ends[read];
        realigned_interval_for_candidate(result, indices.front(), read, tid, start, end);
        if (tid != first.tid || start < 0 || end <= start ||
            start >= event_end || end <= event_start)
            continue;
        if (read < result.annotation_read_qualified.size() &&
            result.annotation_read_qualified[read] == 0)
            continue;

        double reference = kMissing;
        for (const auto index : indices)
            reference = std::max(reference, result.reference_read_likelihoods[index][read]);
        if (std::isfinite(reference) && reference > kMissing) {
            requests.push_back(fastgatk::pairhmm::ReadAlleleBestRequest{
                static_cast<std::uint32_t>(row_count), 0U, reference});
        }
        reference_by_row.push_back(reference);
        for (std::size_t alternate = 0; alternate < indices.size(); ++alternate) {
            const auto likelihood = result.allele_read_likelihoods[indices[alternate]][read];
            if (std::isfinite(likelihood) && likelihood > kMissing) {
                requests.push_back(fastgatk::pairhmm::ReadAlleleBestRequest{
                    static_cast<std::uint32_t>(row_count),
                    static_cast<std::uint32_t>(alternate + 1), likelihood});
            }
        }
        ++row_count;
    }
    if (row_count == 0 || requests.empty()) return std::nullopt;

    const auto reduction = fastgatk::pairhmm::reduce_read_allele_best_kokkos(requests, row_count);
    if (reduction.best_second_allele_by_row.size() != row_count * 2 ||
        reduction.best_second_likelihood_by_row.size() != row_count * 2)
        return std::nullopt;

    EmittedSomaticReadDepth output;
    output.allele_depths.assign(indices.size() + 1, 0U);
    constexpr auto kMissingAllele = std::numeric_limits<std::uint32_t>::max();
    for (std::size_t row = 0; row < row_count; ++row) {
        const auto best = reduction.best_second_allele_by_row[row * 2];
        const auto best_value = reduction.best_second_likelihood_by_row[row * 2];
        const auto second_value = reduction.best_second_likelihood_by_row[row * 2 + 1];
        const auto margin = best_value - second_value;
        const bool supported = !(best == kMissingAllele ||
            best >= output.allele_depths.size() ||
            !std::isfinite(best_value) ||
            (best != 0U && row < reference_by_row.size() &&
             std::isfinite(reference_by_row[row]) && reference_by_row[row] > kMissing &&
             best_value - reference_by_row[row] < kGatkMutectLog10ReferenceTieMargin) ||
            !(margin > kGatkMutectLog10InformativeMargin));
        if (!supported)
            continue;
        if (output.depth != std::numeric_limits<std::uint32_t>::max()) ++output.depth;
        if (output.allele_depths[best] != std::numeric_limits<std::uint32_t>::max())
            ++output.allele_depths[best];
    }
    return output;
}

std::optional<EmittedSomaticReadDepth> emitted_somatic_read_depth(
    const fastgatk::calling::Result& result,
    const fastgatk::io::ReadBatch& reads,
    const std::vector<fastgatk::calling::AssemblyCandidate>& emitted_candidates) {
    if (result.assembly_region_likelihood_results.empty())
        return emitted_somatic_read_depth_local(result, reads, emitted_candidates);

    // A flattened walker result has per-AssemblyRegion likelihood row
    // ordinals.  All ALTs at one emitted VCF locus are owned by one active
    // region; resolve that owner before invoking the Kokkos row reduction.
    const auto& first = emitted_candidates.front();
    const fastgatk::calling::Result* owner = nullptr;
    int owner_priority = -1;
    for (const auto& part_pointer : result.assembly_region_likelihood_results) {
        if (!part_pointer) continue;
        const auto& part = *part_pointer;
        bool has_all = true;
        for (const auto& emitted : emitted_candidates) {
            if (candidate_index(part, emitted) == std::numeric_limits<std::size_t>::max()) {
                has_all = false;
                break;
            }
        }
        if (!has_all) continue;
        int priority = 1;
        for (const auto& region : part.calling_regions) {
            if (region.tid == first.tid && first.position >= region.active_start &&
                first.position <= region.active_end) {
                priority = 2;
                break;
            }
        }
        if (priority > owner_priority) {
            owner = &part;
            owner_priority = priority;
        }
    }
    return owner == nullptr ? std::nullopt
                            : emitted_somatic_read_depth_local(*owner, reads, emitted_candidates);
}

std::optional<fastgatk::kernels::SomaticLikelihoodResult> calculate_somatic(
    const fastgatk::calling::Result& result,
    const double minimum_allele_fraction = 0.0,
    const bool include_non_reference_allele = false) {
    // Mutect2Engine evaluates SomaticGenotypingEngine inside callRegion(),
    // before the walker merges VariantContexts from neighboring
    // AssemblyRegions.  A flattened Result intentionally has no one global
    // read-likelihood matrix (read ordinals are region-local), so reduce the
    // completed regional matrices back onto the flattened candidate list
    // rather than accidentally treating the merged object as unscored.
    if (!result.assembly_region_likelihood_results.empty()) {
        const bool debug_somatic_region_reduce =
            std::getenv("FASTGATK_DEBUG_REGION_SCHEDULING") != nullptr;
        fastgatk::kernels::SomaticLikelihoodResult output;
        const auto count = result.candidates.size();
        const auto missing = -std::numeric_limits<double>::infinity();
        output.tlod.assign(count, missing);
        output.non_reference_tlod.assign(count, missing);
        output.non_reference_allele_fraction.assign(count, 0.0);
        output.normal_log10_odds.assign(count, missing);
        output.best_allele_fraction.assign(count, 0.0);
        output.reference_log10_likelihood.assign(count, missing);
        output.best_log10_likelihood.assign(count, missing);
        output.informative_reads.assign(count, 0);
        std::vector<int> ownership(count, -1);
        // AssemblyRegionWalker merges one VariantContext per exact EventMap
        // allele.  The former linear scan for every local candidate made the
        // Host-only regional reducer O(regions * local_alleles * all_alleles)
        // on real whole-contig traversals.  It could keep an otherwise
        // completed HCC1143 run in the reducer for hours.  Build the same
        // identity lookup once; `emplace` deliberately keeps the first entry,
        // which is exactly what the previous find_if selected if a caller
        // supplied duplicate compatibility candidates.
        using MergedAlleleKey =
            std::tuple<std::int32_t, std::int32_t, std::string, std::string>;
        const auto allele_key = [](const auto& candidate) {
            return MergedAlleleKey{candidate.tid, candidate.position,
                                   allele_ref(candidate), allele_alt(candidate)};
        };
        std::map<MergedAlleleKey, std::size_t> merged_candidate_index;
        for (std::size_t index = 0; index < result.candidates.size(); ++index)
            merged_candidate_index.emplace(allele_key(result.candidates[index]), index);
        if (debug_somatic_region_reduce) {
            std::cerr << "[FASTGATK_SOMATIC_REGION_REDUCE] stage=begin regions="
                      << result.assembly_region_likelihood_results.size()
                      << " candidates=" << result.candidates.size() << '\n';
        }
        for (const auto& part_pointer : result.assembly_region_likelihood_results) {
            if (!part_pointer) continue;
            const auto& part = *part_pointer;
            const auto local = calculate_somatic(
                part, minimum_allele_fraction, include_non_reference_allele);
            if (!local.has_value()) continue;
            output.evidence_groups += local->evidence_groups;
            output.prepare_seconds += local->prepare_seconds;
            output.seconds += local->seconds;
            if (output.execution_space.empty())
                output.execution_space = local->execution_space;
            for (std::size_t local_index = 0; local_index < part.candidates.size(); ++local_index) {
                if (local_index >= local->tlod.size()) continue;
                const auto merged = merged_candidate_index.find(
                    allele_key(part.candidates[local_index]));
                if (merged == merged_candidate_index.end()) continue;
                const auto merged_index = merged->second;
                // Overlapping AssemblyRegion halos can contain the same
                // event.  GATK gives its active core ownership; a halo value
                // is only a fallback if the owning region did not retain the
                // exact allele.
                int priority = 1;
                if (part.calling_regions.size() == 1) {
                    const auto& region = part.calling_regions.front();
                    if (part.candidates[local_index].tid == region.tid &&
                        part.candidates[local_index].position >= region.active_start &&
                        part.candidates[local_index].position <= region.active_end)
                        priority = 2;
                }
                if (priority < ownership[merged_index]) continue;
                ownership[merged_index] = priority;
                output.tlod[merged_index] = local->tlod[local_index];
                if (local_index < local->non_reference_tlod.size())
                    output.non_reference_tlod[merged_index] =
                        local->non_reference_tlod[local_index];
                if (local_index < local->non_reference_allele_fraction.size())
                    output.non_reference_allele_fraction[merged_index] =
                        local->non_reference_allele_fraction[local_index];
                if (local_index < local->normal_log10_odds.size())
                    output.normal_log10_odds[merged_index] =
                        local->normal_log10_odds[local_index];
                if (local_index < local->best_allele_fraction.size())
                    output.best_allele_fraction[merged_index] =
                        local->best_allele_fraction[local_index];
                if (local_index < local->reference_log10_likelihood.size())
                    output.reference_log10_likelihood[merged_index] =
                        local->reference_log10_likelihood[local_index];
                if (local_index < local->best_log10_likelihood.size())
                    output.best_log10_likelihood[merged_index] =
                        local->best_log10_likelihood[local_index];
                if (local_index < local->informative_reads.size())
                    output.informative_reads[merged_index] = local->informative_reads[local_index];
            }
        }
        if (debug_somatic_region_reduce) {
            std::cerr << "[FASTGATK_SOMATIC_REGION_REDUCE] stage=end regions="
                      << result.assembly_region_likelihood_results.size()
                      << " candidates=" << result.candidates.size() << '\n';
        }
        return output;
    }
    const auto grouped = group_somatic_evidence(result);
    if (!grouped.has_value()) return std::nullopt;
    fastgatk::kernels::SomaticLikelihoodResult output;
    if (grouped->sparse_rows) {
        if (grouped->sparse_candidate_offsets == nullptr ||
            grouped->sparse_likelihood_rows == nullptr)
            return std::nullopt;
        std::vector<double> sparse_reference;
        std::vector<double> sparse_alternate;
        sparse_reference.reserve(grouped->sparse_likelihood_rows->size());
        sparse_alternate.reserve(grouped->sparse_likelihood_rows->size());
        for (const auto& row : *grouped->sparse_likelihood_rows) {
            sparse_reference.push_back(row.reference);
            sparse_alternate.push_back(row.alternate);
        }
        output = fastgatk::kernels::calculate_somatic_likelihood_sparse_kokkos(
            sparse_reference, sparse_alternate, *grouped->sparse_candidate_offsets,
            result.candidates.size(), grouped->read_count, 101, minimum_allele_fraction);
    } else {
        output = fastgatk::kernels::calculate_somatic_likelihood_kokkos(
            grouped->reference, grouped->alternate, result.candidates.size(),
            grouped->read_count, 101, minimum_allele_fraction);
    }
    if (include_non_reference_allele)
        output.non_reference_tlod.assign(result.candidates.size(),
            -std::numeric_limits<double>::infinity());
    if (include_non_reference_allele)
        output.non_reference_allele_fraction.assign(result.candidates.size(), 0.0);

    // SomaticGenotypingEngine evaluates one merged allele list per locus and
    // obtains each TLOD as logEvidence(all alleles) -
    // logEvidence(all alleles except this ALT).  The historical native path
    // evaluated every candidate as an independent biallelic model, which is
    // numerically correct for singleton loci but overstates evidence when
    // sibling ALTs compete for the same read.  Keep the established matrix
    // layout and overwrite loci with two or more concrete ALTs, or a GATK
    // symbolic spanning-deletion (`*`) allele.
    using SomaticGroupKey = std::tuple<std::int32_t, std::int32_t, std::string>;
    std::map<SomaticGroupKey, std::vector<std::size_t>> groups;
    // The somatic likelihood matrix is defined by AssemblyResult's retained
    // EventMap.  Result also carries raw CIGAR candidates for annotations;
    // when provenance exists, those raw siblings must not become additional
    // alleles in all-alleles-vs-without-ALT TLOD calculations.
    const auto event_map_owns_candidate = [&](const std::size_t candidate) {
        if (candidate >= result.candidates.size()) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position &&
                    event.reference == allele_ref(item) && event.alternate == allele_alt(item);
            });
        });
    };
    const auto event_map_has_locus = [&](const std::size_t candidate) {
        if (candidate >= result.candidates.size()) return false;
        const auto& item = result.candidates[candidate];
        return std::any_of(result.graph_haplotype_event_maps.begin(),
                           result.graph_haplotype_event_maps.end(),
                           [&](const auto& event_map) {
            return std::any_of(event_map.events.begin(), event_map.events.end(),
                               [&](const auto& event) {
                return event.tid == item.tid && event.position == item.position;
            });
        });
    };
    const bool has_pairhmm_event_map_ownership =
        result.likelihood_candidate_event_map_owned.size() == result.candidates.size();
    for (std::size_t candidate = 0; candidate < result.candidates.size(); ++candidate) {
        // The flattened native candidate list deliberately retains raw
        // pileup/CIGAR observations for annotations and future local assembly,
        // even when no retained AssemblyRegion haplotype realizes the allele.
        // SomaticGenotypingEngine, however, builds its allele matrix from the
        // EventMap of the assembled haplotypes.  Do not let an unmodeled raw
        // sibling turn a concrete biallelic EventMap allele into a synthetic
        // multiallelic model: it has no likelihood row and therefore is not an
        // allele in GATK's matrix.  Concrete zero-probability alleles still
        // have a sparse row and remain included.
        const auto& item = result.candidates[candidate];
        // The compound low-BQ recovery is valid only when the local EventMap
        // has no allele at this locus.  Once it does, GATK's disabled
        // pileup-allele path cannot introduce a raw competing ALT.
        const bool explicit_recovery = item.forced_by_tumor ||
            (item.somatic_low_quality_assembly_recovery &&
             (result.graph_haplotype_event_maps.empty() ||
              !event_map_has_locus(candidate)));
        if (!explicit_recovery) {
            // A regional Result can be rebuilt from an AssemblyResult
            // snapshot without retaining its complete EventMap objects.  In
            // that normal path, use PairHMM-time ownership; the direct map
            // lookup remains the compatibility fallback for older producers.
            if (has_pairhmm_event_map_ownership &&
                result.likelihood_candidate_event_map_owned[candidate] == 0)
                continue;
            if (!has_pairhmm_event_map_ownership &&
                !result.graph_haplotype_event_maps.empty() &&
                !event_map_owns_candidate(candidate))
                continue;
        }
        if (grouped->groups_for_candidate(candidate).empty()) continue;
        groups[SomaticGroupKey{item.tid, item.position, allele_ref(item)}].push_back(candidate);
    }
    for (const auto& [key, indices] : groups) {
        (void)key;
        std::set<std::uint32_t> active_groups;
        for (const auto candidate : indices) {
            for (const auto group : grouped->groups_for_candidate(candidate))
                active_groups.insert(group);
        }
        if (active_groups.empty()) continue;
        bool has_spanning_deletion = false;
        for (const auto group : active_groups) {
            for (const auto candidate : indices) {
                const auto spanning = grouped->spanning_deletion_at(candidate, group);
                if (std::isfinite(spanning) && spanning > -1.0e299) {
                    has_spanning_deletion = true;
                    break;
                }
            }
            if (has_spanning_deletion) break;
        }
        if (indices.size() < 2 && !has_spanning_deletion && !include_non_reference_allele)
            continue;
        const auto allele_count = indices.size() + 1 +
            static_cast<std::size_t>(has_spanning_deletion);
        std::vector<double> matrix;
        matrix.reserve(allele_count * active_groups.size());
        const auto append_row = [&](const std::size_t candidate, const bool reference) {
            for (const auto group : active_groups) {
                matrix.push_back(reference ? grouped->reference_at(candidate, group) :
                                             grouped->alternate_at(candidate, group));
            }
        };
        // REF and symbolic `*` are each one shared merged allele. Raw
        // candidate siblings may not carry the sparse EventMap row, so do
        // not let their arbitrary sort order turn a present row into
        // missing evidence.
        for (const auto group : active_groups) {
            double reference = -1.0e299;
            for (const auto candidate : indices)
                reference = std::max(reference, grouped->reference_at(candidate, group));
            matrix.push_back(reference);
        }
        for (const auto candidate : indices) append_row(candidate, false);
        if (has_spanning_deletion) {
            // GATK adds exactly one symbolic `*` allele to the merged locus.
            for (const auto group : active_groups) {
                double spanning_deletion = -1.0e299;
                for (const auto candidate : indices)
                    spanning_deletion = std::max(spanning_deletion,
                        grouped->spanning_deletion_at(candidate, group));
                matrix.push_back(spanning_deletion);
            }
        }
        if (indices.size() >= 2 || has_spanning_deletion) {
            const auto multi = fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                matrix, allele_count, active_groups.size(), minimum_allele_fraction);
            for (std::size_t alt = 0; alt < indices.size(); ++alt) {
                const auto candidate = indices[alt];
                if (alt >= multi.tlod.size()) continue;
                output.tlod[candidate] = multi.tlod[alt];
                // Keep the biallelic normal LOD.  The native PairHMM aggregate
                // materializes each candidate's REF row as "all haplotypes that
                // do not carry this ALT"; it therefore includes sibling ALT
                // haplotypes.  That is suitable for the all-alleles-vs-without-
                // this-ALT TLOD calculation below, but not for GATK's diploid
                // REF/ALT normal comparison.  Reusing one sibling's REF row in
                // the multiallelic normal calculation turned a demonstrably
                // germline site (NLOD -35.13) into +5.92.
                output.best_allele_fraction[candidate] = multi.best_allele_fraction[alt];
                output.reference_log10_likelihood[candidate] = multi.reference_log10_likelihood[alt];
                output.best_log10_likelihood[candidate] = multi.best_log10_likelihood[alt];
                output.informative_reads[candidate] = multi.informative_reads[alt];
            }
            output.prepare_seconds += multi.prepare_seconds;
            output.seconds += multi.seconds;
        }
        if (include_non_reference_allele) {
            // AlleleLikelihoods.addNonReferenceAllele assigns <NON_REF> the
            // median likelihood strictly below each fragment's best concrete
            // allele likelihood.  If no concrete row is lower, the Java
            // implementation uses the shared best value (there are always at
            // least REF and one concrete ALT at an EventMap locus).  Construct
            // that Host-owned row from the existing PairHMM matrix, then keep
            // the all-alleles-vs-without-allele numerical reduction inside the
            // same Kokkos multiallelic kernel as concrete TLOD.
            std::vector<double> with_non_reference = matrix;
            std::vector<double> non_reference_row(active_groups.size(),
                -std::numeric_limits<double>::infinity());
            for (std::size_t group_index = 0; group_index < active_groups.size(); ++group_index) {
                double best = -std::numeric_limits<double>::infinity();
                for (std::size_t allele = 0; allele < allele_count; ++allele) {
                    const auto likelihood = matrix[allele * active_groups.size() + group_index];
                    if (std::isfinite(likelihood) && likelihood > -1.0e299)
                        best = std::max(best, likelihood);
                }
                if (!std::isfinite(best) || best <= -1.0e299) continue;
                std::vector<double> qualified;
                qualified.reserve(allele_count);
                for (std::size_t allele = 0; allele < allele_count; ++allele) {
                    const auto likelihood = matrix[allele * active_groups.size() + group_index];
                    if (std::isfinite(likelihood) && likelihood > -1.0e299 && likelihood < best)
                        qualified.push_back(likelihood);
                }
                if (qualified.empty()) {
                    non_reference_row[group_index] = best;
                } else {
                    std::sort(qualified.begin(), qualified.end());
                    const auto middle = qualified.size() / 2U;
                    non_reference_row[group_index] = (qualified.size() & 1U) != 0U
                        ? qualified[middle]
                        : (qualified[middle - 1U] + qualified[middle]) / 2.0;
                }
            }
            with_non_reference.insert(with_non_reference.end(),
                non_reference_row.begin(), non_reference_row.end());
            const auto with_non_reference_likelihood =
                fastgatk::kernels::calculate_somatic_multiallelic_likelihood_kokkos(
                    with_non_reference, allele_count + 1U, active_groups.size(),
                    minimum_allele_fraction);
            if (!with_non_reference_likelihood.tlod.empty()) {
                const auto non_reference_lod = with_non_reference_likelihood.tlod.back();
                const auto non_reference_fraction =
                    with_non_reference_likelihood.best_allele_fraction.empty() ? 0.0 :
                    with_non_reference_likelihood.best_allele_fraction.back();
                for (std::size_t alt = 0; alt < indices.size(); ++alt) {
                    const auto candidate = indices[alt];
                    // In reference-confidence mode GATK calculates concrete
                    // ALT TLOD and posterior AF after adding <NON_REF>, not
                    // from the earlier biallelic convenience matrix.  The
                    // shared Kokkos result is therefore the source-authority
                    // value for both concrete and symbolic writer fields.
                    if (alt < with_non_reference_likelihood.tlod.size())
                        output.tlod[candidate] = with_non_reference_likelihood.tlod[alt];
                    if (alt < with_non_reference_likelihood.best_allele_fraction.size())
                        output.best_allele_fraction[candidate] =
                            with_non_reference_likelihood.best_allele_fraction[alt];
                    output.non_reference_tlod[candidate] = non_reference_lod;
                    output.non_reference_allele_fraction[candidate] = non_reference_fraction;
                }
            }
            output.prepare_seconds += with_non_reference_likelihood.prepare_seconds;
            output.seconds += with_non_reference_likelihood.seconds;
        }
    }
    // Keep the numerical path completely unchanged, but expose the exact
    // candidate-by-fragment matrix when explicitly requested.  This is the
    // matrix that is handed to the Kokkos SomaticLikelihoods reduction after
    // Host-side fragment grouping/marginalization, and is therefore the
    // appropriate audit boundary for a source-level transcription of GATK's
    // SomaticLikelihoodsEngine.  Do not print a synthesized dense matrix:
    // absent sparse rows represent GATK's missing-allele state and must stay
    // absent in the verifier as well.
    if (std::getenv("FASTGATK_DEBUG_TLOD") != nullptr) {
        std::ostringstream debug;
        debug << std::setprecision(17);
        for (std::size_t candidate = 0; candidate < result.candidates.size(); ++candidate) {
            const auto fragment_groups = grouped->groups_for_candidate(candidate);
            if (fragment_groups.empty()) continue;
            const auto& allele = result.candidates[candidate];
            debug << "[FASTGATK_DEBUG_TLOD: candidate] pos=" << allele.position
                  << " ref=" << allele_ref(allele)
                  << " alt=" << allele_alt(allele) << '\n';
            for (const auto group : fragment_groups) {
                debug << "  g=" << group
                      << " ref=" << grouped->reference_at(candidate, group)
                      << " alt=" << grouped->alternate_at(candidate, group)
                      << " span=" << grouped->spanning_deletion_at(candidate, group);
                if (group < result.somatic_group_representative_sources.size())
                    debug << " source=" << result.somatic_group_representative_sources[group];
                if (group < result.somatic_group_context_ordinals.size())
                    debug << " context=" << result.somatic_group_context_ordinals[group];
                debug << '\n';
            }
            if (candidate < output.tlod.size() && std::isfinite(output.tlod[candidate])) {
                debug << "[FASTGATK_TLOD_FULL] pos=" << allele.position
                      << " ref=" << allele_ref(allele)
                      << " alt=" << allele_alt(allele)
                      << " tlod=" << output.tlod[candidate] << '\n';
            }
        }
        std::cerr << debug.str();
    }
    return output;
}

struct OrientationCounts {
    std::uint32_t f1r2 = 0;
    std::uint32_t r1f2 = 0;
};

// Return the observed base aligned to one reference coordinate.  This is the
// same M/=/X-only projection contract as io::project_read_offset(), but the
// caller already knows the one coordinate it needs.  Walking a record's CIGAR
// once avoids the former O(read_length x cigar_ops) sequence of defensive
// public-API projections in the full Mutect2 posterior pass.
std::optional<std::uint8_t> observed_base_at_reference_position(
    const fastgatk::io::ReadBatch& reads, const std::size_t record,
    const std::int32_t reference_position) {
    if (record >= reads.records() || record >= reads.positions.size() ||
        reads.positions[record] < 0 || record + 1 >= reads.offsets.size())
        return std::nullopt;
    const auto begin = static_cast<std::size_t>(reads.offsets[record]);
    const auto end = static_cast<std::size_t>(reads.offsets[record + 1]);
    if (begin > end || end > reads.bases.size()) return std::nullopt;
    const auto read_length = end - begin;
    const auto target = static_cast<std::int64_t>(reference_position);
    std::int64_t reference_cursor = reads.positions[record];
    std::size_t read_cursor = 0;
    if (!reads.has_cigar()) {
        if (target < reference_cursor) return std::nullopt;
        const auto offset = static_cast<std::uint64_t>(target - reference_cursor);
        if (offset >= read_length) return std::nullopt;
        return reads.bases[begin + static_cast<std::size_t>(offset)];
    }
    if (!reads.cigar_record_layout_valid(record)) return std::nullopt;
    const auto cigar_begin = static_cast<std::size_t>(reads.cigar_offsets[record]);
    const auto cigar_end = static_cast<std::size_t>(reads.cigar_offsets[record + 1]);
    for (std::size_t index = cigar_begin; index < cigar_end; ++index) {
        const auto operation = fastgatk::io::CigarOp::unpack(reads.cigar_ops[index]);
        if (!operation.valid()) return std::nullopt;
        const auto length = static_cast<std::size_t>(operation.length);
        if (operation.projects_base()) {
            if (target >= reference_cursor &&
                static_cast<std::uint64_t>(target - reference_cursor) < length) {
                const auto offset = read_cursor + static_cast<std::size_t>(
                    target - reference_cursor);
                if (offset >= read_length) return std::nullopt;
                return reads.bases[begin + offset];
            }
            if (length > read_length - std::min(read_cursor, read_length))
                return std::nullopt;
            read_cursor += length;
            reference_cursor += static_cast<std::int64_t>(operation.length);
        } else {
            if (operation.consumes_read()) {
                if (length > read_length - std::min(read_cursor, read_length))
                    return std::nullopt;
                read_cursor += length;
            }
            if (operation.consumes_reference())
                reference_cursor += static_cast<std::int64_t>(operation.length);
        }
    }
    return std::nullopt;
}

// Mutect2's full orientation-bias model is intentionally still a fallback,
// but emitting deterministic strand-support counts makes the native sidecar
// useful to FilterMutectCalls and to downstream artifact diagnostics.
OrientationCounts orientation_counts(const fastgatk::io::ReadBatch& reads,
                                     const fastgatk::calling::AssemblyCandidate& candidate) {
    const auto alternate = allele_alt(candidate);
    if (alternate.size() != 1) return {};
    const auto expected = static_cast<char>(std::toupper(static_cast<unsigned char>(alternate.front())));
    OrientationCounts counts;
    for (std::size_t record = 0; record < reads.records(); ++record) {
        if (record >= reads.tids.size() || reads.tids[record] != candidate.tid) continue;
        const auto observed_base = observed_base_at_reference_position(
            reads, record, candidate.position);
        if (!observed_base.has_value()) continue;
        const auto observed = static_cast<char>(std::toupper(
            static_cast<unsigned char>(*observed_base)));
        if (observed != expected) continue;
        const bool reverse = !reads.flags.empty() && (reads.flags[record] & 0x10U) != 0;
        if (reverse) ++counts.r1f2;
        else ++counts.f1r2;
    }
    return counts;
}

// OrientationBiasReadCounts is evaluated from the local
// AlleleLikelihoods collection, not by rescanning a walker's complete BAM
// for every AssemblyRegion EventMap allele.  In a partitioned Result the
// flattened public candidate is only a merge view, so first resolve the one
// active-core local collection that owns it.  This is both the GATK evidence
// boundary and avoids a full chromosome scan for every local posterior.
OrientationCounts orientation_counts(
    const fastgatk::calling::Result& result,
    const fastgatk::io::ReadBatch& reads, const std::size_t candidate_index) {
    if (candidate_index >= result.candidates.size()) return {};
    const auto& candidate = result.candidates[candidate_index];
    if (!result.assembly_region_likelihood_results.empty()) {
        const fastgatk::calling::Result* owner = nullptr;
        std::size_t owner_index = std::numeric_limits<std::size_t>::max();
        int owner_priority = -1;
        for (const auto& part_pointer : result.assembly_region_likelihood_results) {
            if (!part_pointer) continue;
            const auto& part = *part_pointer;
            for (std::size_t local_index = 0; local_index < part.candidates.size(); ++local_index) {
                const auto& local = part.candidates[local_index];
                if (local.tid != candidate.tid || local.position != candidate.position ||
                    allele_ref(local) != allele_ref(candidate) ||
                    allele_alt(local) != allele_alt(candidate))
                    continue;
                int priority = 1;
                if (part.calling_regions.size() == 1) {
                    const auto& region = part.calling_regions.front();
                    if (local.tid == region.tid && local.position >= region.active_start &&
                        local.position <= region.active_end)
                        priority = 2;
                }
                if (priority > owner_priority) {
                    owner = &part;
                    owner_index = local_index;
                    owner_priority = priority;
                }
            }
        }
        return owner == nullptr ? OrientationCounts{}
                                : orientation_counts(*owner, reads, owner_index);
    }
    const auto read_count = result.likelihood_read_names.size();
    if (read_count == 0 || result.likelihood_read_tids.size() != read_count ||
        result.likelihood_read_starts.size() != read_count ||
        result.likelihood_read_ends.size() != read_count ||
        candidate_index >= result.reference_read_likelihoods.size() ||
        candidate_index >= result.allele_read_likelihoods.size() ||
        result.reference_read_likelihoods[candidate_index].size() != read_count ||
        result.allele_read_likelihoods[candidate_index].size() != read_count)
        return {};

    // An AssemblyRegion child preserves its outer ReadBatch ordinal in this
    // vector.  That is the same identity carried by GATK's local
    // AlleleLikelihoods collection, and it lets the normal full-contig path
    // resolve its handful of local reads directly.  Rebuilding a name/start/
    // end map over the complete BAM for every candidate was semantically
    // redundant and made the final posterior reduction quadratic in the
    // number of regions.  Validate that the ordinal still describes this
    // input batch before using it; the older keyed fallback below remains for
    // direct API results whose provenance is not relative to `reads`.
    bool direct_source_records =
        result.likelihood_read_source_records.size() == read_count;
    for (std::size_t read = 0; direct_source_records && read < read_count; ++read) {
        const auto record = static_cast<std::size_t>(
            result.likelihood_read_source_records[read]);
        direct_source_records = record < reads.records() && record < reads.tids.size() &&
            record < reads.positions.size() &&
            reads.tids[record] == result.likelihood_read_tids[read] &&
            static_cast<std::int64_t>(reads.positions[record]) ==
                result.likelihood_read_starts[read] &&
            fastgatk::io::reference_end(reads, record) ==
                result.likelihood_read_ends[read];
    }
    using ReadKey = std::tuple<std::string, std::int32_t, std::int64_t, std::int64_t>;
    std::map<ReadKey, std::size_t> input_record_by_key;
    if (!direct_source_records) {
        for (std::size_t record = 0; record < reads.records(); ++record) {
            const auto tid = record < reads.tids.size() ? reads.tids[record] : -1;
            const auto start = record < reads.positions.size()
                ? static_cast<std::int64_t>(reads.positions[record]) : -1;
            input_record_by_key.emplace(ReadKey{read_name_for_record(reads, record), tid, start,
                                                 fastgatk::io::reference_end(reads, record)}, record);
        }
    }

    constexpr double kMissingLikelihood = -1.0e299;
    const auto alternate = allele_alt(candidate);
    if (alternate.size() != 1) return {};
    const auto expected = static_cast<char>(std::toupper(
        static_cast<unsigned char>(alternate.front())));
    const auto event_start = std::max<std::int64_t>(
        0, static_cast<std::int64_t>(candidate.position) - 2);
    const auto event_end = std::max<std::int64_t>(event_start + 1,
        static_cast<std::int64_t>(candidate.position) +
        static_cast<std::int64_t>(std::max<std::size_t>(1, allele_ref(candidate).size())) + 2);
    OrientationCounts counts;
    for (std::size_t read = 0; read < read_count; ++read) {
        if (read < result.annotation_read_qualified.size() &&
            result.annotation_read_qualified[read] == 0)
            continue;
        const auto reference = result.reference_read_likelihoods[candidate_index][read];
        const auto alternative = result.allele_read_likelihoods[candidate_index][read];
        if ((!std::isfinite(reference) || reference <= kMissingLikelihood) &&
            (!std::isfinite(alternative) || alternative <= kMissingLikelihood))
            continue;
        std::size_t input_record = std::numeric_limits<std::size_t>::max();
        if (direct_source_records) {
            input_record = result.likelihood_read_source_records[read];
        } else {
            const auto input = input_record_by_key.find(ReadKey{
                result.likelihood_read_names[read], result.likelihood_read_tids[read],
                result.likelihood_read_starts[read], result.likelihood_read_ends[read]});
            if (input == input_record_by_key.end()) continue;
            input_record = input->second;
        }
        std::int32_t tid = result.likelihood_read_tids[read];
        std::int64_t start = result.likelihood_read_starts[read];
        std::int64_t end = result.likelihood_read_ends[read];
        realigned_interval_for_candidate(result, candidate_index, read, tid, start, end);
        if (tid != candidate.tid || start < 0 || end <= start ||
            start >= event_end || end <= event_start)
            continue;
        const auto observed_base = observed_base_at_reference_position(
            reads, input_record, candidate.position);
        if (!observed_base.has_value() || static_cast<char>(std::toupper(
                static_cast<unsigned char>(*observed_base))) != expected)
            continue;
        const bool reverse = input_record < reads.flags.size() &&
            (reads.flags[input_record] & 0x10U) != 0;
        if (reverse) ++counts.r1f2;
        else ++counts.f1r2;
    }
    return counts;
}

template <typename T>
int median_integer(std::vector<T> values) {
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    const auto middle = values.size() / 2;
    if (values.size() % 2 != 0) return static_cast<int>(values[middle]);
    return static_cast<int>(std::llround((static_cast<double>(values[middle - 1]) +
                                         static_cast<double>(values[middle])) / 2.0));
}

void append_per_allele_annotation_values(PerAlleleAnnotationValues& destination,
                                         const PerAlleleAnnotationValues& source) {
    destination.base_qualities.insert(destination.base_qualities.end(),
                                      source.base_qualities.begin(), source.base_qualities.end());
    destination.mapping_qualities.insert(destination.mapping_qualities.end(),
                                         source.mapping_qualities.begin(), source.mapping_qualities.end());
    destination.fragment_lengths.insert(destination.fragment_lengths.end(),
                                        source.fragment_lengths.begin(), source.fragment_lengths.end());
    destination.read_positions.insert(destination.read_positions.end(),
                                      source.read_positions.begin(), source.read_positions.end());
}

int per_allele_annotation_statistic(const PerAlleleAnnotationValues& values, const char metric) {
    if (metric == 'b') return median_integer(values.base_qualities);
    if (metric == 'f') return median_integer(values.fragment_lengths);
    if (metric == 'm') return values.mapping_qualities.empty()
        ? 60 : median_integer(values.mapping_qualities);
    return values.read_positions.empty() ? 50 : median_integer(values.read_positions);
}

std::string per_allele_annotation_values(const PerAlleleAnnotationValues& reference,
                                         const std::vector<PerAlleleAnnotationValues>& alternates,
                                         const char metric) {
    std::ostringstream output;
    output << per_allele_annotation_statistic(reference, metric);
    for (const auto& alternate : alternates)
        output << ',' << per_allele_annotation_statistic(alternate, metric);
    return output.str();
}

std::string alternate_per_allele_annotation_values(
    const std::vector<PerAlleleAnnotationValues>& alternates, const char metric) {
    std::ostringstream output;
    for (std::size_t index = 0; index < alternates.size(); ++index) {
        if (index != 0) output << ',';
        output << per_allele_annotation_statistic(alternates[index], metric);
    }
    return output.str();
}

int base_index_for_f1r2(const char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

std::vector<std::string> f1r2_contexts() {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    std::vector<std::string> contexts;
    contexts.reserve(64);
    for (const auto left : bases)
        for (const auto middle : bases)
            for (const auto right : bases)
                contexts.emplace_back(std::string{left, middle, right});
    return contexts;
}

OrientationCounts orientation_counts_for_base(
    const fastgatk::io::ReadBatch& reads,
    const fastgatk::calling::AssemblyCandidate& candidate,
    const char expected_base) {
    OrientationCounts counts;
    const auto expected = static_cast<char>(std::toupper(static_cast<unsigned char>(expected_base)));
    for (std::size_t record = 0; record < reads.records(); ++record) {
        if (record >= reads.tids.size() || reads.tids[record] != candidate.tid) continue;
        const auto begin = reads.offsets[record];
        const auto end = reads.offsets[record + 1];
        for (std::size_t offset = begin; offset < end; ++offset) {
            fastgatk::io::ReadProjection projection;
            if (!fastgatk::io::project_read_offset(reads, record, offset - begin, projection) ||
                projection.reference_position != candidate.position) continue;
            const auto observed = static_cast<char>(std::toupper(
                static_cast<unsigned char>(reads.bases[offset])));
            if (observed != expected) continue;
            const bool reverse = !reads.flags.empty() && (reads.flags[record] & 0x10U) != 0;
            if (reverse) ++counts.r1f2;
            else ++counts.f1r2;
            break;
        }
    }
    return counts;
}

std::uint32_t depth_at_position(const fastgatk::io::ReadBatch& reads,
                                const std::int32_t tid,
                                const std::int32_t position) {
    std::uint32_t depth = 0;
    for (std::size_t record = 0; record < reads.records(); ++record) {
        if (record >= reads.tids.size() || reads.tids[record] != tid) continue;
        const auto begin = reads.offsets[record];
        const auto end = reads.offsets[record + 1];
        for (std::size_t offset = begin; offset < end; ++offset) {
            fastgatk::io::ReadProjection projection;
            if (fastgatk::io::project_read_offset(reads, record, offset - begin, projection) &&
                projection.reference_position == position) {
                ++depth;
                break;
            }
        }
    }
    return depth;
}

std::optional<std::string> f1r2_context(
    const std::vector<std::string>& references,
    const fastgatk::calling::AssemblyCandidate& candidate) {
    if (candidate.tid < 0 || static_cast<std::size_t>(candidate.tid) >= references.size() ||
        candidate.position < 1) return std::nullopt;
    const auto& reference = references[static_cast<std::size_t>(candidate.tid)];
    const auto position = static_cast<std::size_t>(candidate.position);
    if (position + 1 >= reference.size()) return std::nullopt;
    std::string context = reference.substr(position - 1, 3);
    for (auto& base : context)
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    if (std::any_of(context.begin(), context.end(), [](const char base) {
            return base_index_for_f1r2(base) < 0;
        })) return std::nullopt;
    return context;
}

fastgatk::somatic::F1R2ArchiveSample make_f1r2_archive_sample(
    const fastgatk::io::HeaderSummary& header,
    const fastgatk::io::ReadBatch& tumor_reads,
    const fastgatk::calling::Result& tumor,
    const std::vector<std::string>& references,
    const fastgatk::io::HtsInterval* output_core = nullptr) {
    fastgatk::somatic::F1R2ArchiveSample sample;
    sample.sample = header.samples.empty() ? "TUMOR" : header.samples.front();
    sample.max_depth = 200;
    for (auto& histogram : sample.ref_hist)
        histogram.assign(static_cast<std::size_t>(sample.max_depth + 1), 0);
    for (auto& histogram : sample.alt_hist)
        histogram.assign(static_cast<std::size_t>(sample.max_depth + 1), 0);
    const auto contexts = f1r2_contexts();
    const auto context_index = [&](const std::string& context) {
        const auto found = std::find(contexts.begin(), contexts.end(), context);
        return found == contexts.end() ? std::numeric_limits<std::size_t>::max()
                                       : static_cast<std::size_t>(found - contexts.begin());
    };
    for (const auto& block : tumor.reference_blocks) {
        if (output_core != nullptr && block.tid != output_core->tid) continue;
        const auto block_start = output_core == nullptr ? block.start : std::max(
            block.start, static_cast<std::int32_t>(output_core->start));
        const auto block_end = output_core == nullptr ? block.end : std::min(
            block.end, static_cast<std::int32_t>(output_core->end - 1));
        if (block_end < block_start) continue;
        // CollectF1R2Counts records one reference histogram observation per
        // covered reference locus.  Mutect2's result stores coalesced blocks,
        // so reconstruct the locus stream from the same filtered ReadBatch.
        for (std::int32_t position = block_start; position <= block_end; ++position) {
            fastgatk::calling::AssemblyCandidate representative;
            representative.tid = block.tid;
            representative.position = position;
            const auto context = f1r2_context(references, representative);
            if (!context.has_value()) continue;
            const auto index = context_index(*context);
            if (index == std::numeric_limits<std::size_t>::max()) continue;
            const auto raw_depth = depth_at_position(tumor_reads, block.tid, position);
            if (raw_depth == 0) continue;
            const auto depth = std::min<std::uint32_t>(raw_depth,
                static_cast<std::uint32_t>(sample.max_depth));
            ++sample.ref_hist[index][depth];
            if (position == std::numeric_limits<std::int32_t>::max()) break;
        }
    }
    for (const auto& call : tumor.calls) {
        if (output_core != nullptr &&
            (call.candidate.tid != output_core->tid ||
             call.candidate.position < output_core->start ||
             call.candidate.position >= output_core->end))
            continue;
        const auto reference = allele_ref(call.candidate);
        const auto alternate = allele_alt(call.candidate);
        if (reference.size() != 1 || alternate.size() != 1) continue;
        const auto context = f1r2_context(references, call.candidate);
        if (!context.has_value()) continue;
        const auto index = context_index(*context);
        if (index == std::numeric_limits<std::size_t>::max()) continue;
        const auto alt_index = base_index_for_f1r2(alternate.front());
        const auto ref_index = base_index_for_f1r2(reference.front());
        if (alt_index < 0 || ref_index < 0 || alt_index == ref_index) continue;
        const auto alt_orientation = orientation_counts_for_base(
            tumor_reads, call.candidate, alternate.front());
        const auto ref_orientation = orientation_counts_for_base(
            tumor_reads, call.candidate, reference.front());
        const auto alt_count = call.candidate.alternate_count;
        // A concrete ALT can survive joint genotyping with no informative
        // BestAllele owner (for example, when a sibling ALT owns every read).
        // GATK's standard F1R2 alt table omits such zero-support rows; writing
        // one would violate its depth/count contract and be rejected by
        // LearnReadOrientationModel.
        if (alt_count == 0) continue;
        const auto observed_depth = std::max<std::uint32_t>(call.candidate.depth,
            call.candidate.reference_count + alt_count);
        const auto ref_count = call.candidate.reference_count != 0
            ? call.candidate.reference_count
            : observed_depth > alt_count ? observed_depth - alt_count : 0;
        // The standard alt table is a compact locus table: depth must equal
        // the concrete REF+ALT counts (unlike Mutect2's wider DP annotation).
        const auto table_depth = ref_count + alt_count;
        const auto capped_depth = std::min<std::uint32_t>(table_depth,
            static_cast<std::uint32_t>(sample.max_depth));
        if (alt_count == 1) {
            int rank = 0;
            const auto middle = base_index_for_f1r2((*context)[1]);
            for (int candidate_alt = 0; candidate_alt < alt_index; ++candidate_alt)
                if (candidate_alt != middle) ++rank;
            const auto histogram_index = index * 6 + static_cast<std::size_t>(rank * 2 +
                (alt_orientation.f1r2 == 1 ? 0 : 1));
            if (histogram_index < sample.alt_hist.size())
                ++sample.alt_hist[histogram_index][capped_depth];
        } else {
            sample.alt_rows.push_back(fastgatk::somatic::F1R2ArchiveAltRow{
                *context, static_cast<int>(ref_count), static_cast<int>(alt_count),
                static_cast<int>(std::min<std::uint32_t>(ref_orientation.f1r2, ref_count)),
                static_cast<int>(std::min<std::uint32_t>(alt_orientation.f1r2, alt_count)),
                static_cast<int>(table_depth), alternate.front()});
        }
    }
    return sample;
}

std::optional<fastgatk::kernels::SomaticPosteriorResult> calculate_posterior(
    const Options& options,
    const fastgatk::calling::Result& tumor,
    const fastgatk::io::ReadBatch& tumor_reads,
    const fastgatk::calling::Result* normal) {
    // As for TLOD, posterior inference belongs to each complete
    // AssemblyRegion likelihood matrix.  Merge only the resulting
    // per-allele scalars into the flattened call set, preferring the active
    // core when padded regional halos contain the same event.
    if (!tumor.assembly_region_likelihood_results.empty()) {
        const bool debug_somatic_region_reduce =
            std::getenv("FASTGATK_DEBUG_REGION_SCHEDULING") != nullptr;
        fastgatk::kernels::SomaticPosteriorResult output;
        const auto count = tumor.candidates.size();
        output.somatic_probability.assign(count, 0.0);
        output.germline_probability.assign(count, 0.0);
        output.artifact_probability.assign(count, 1.0);
        output.contamination_adjusted_allele_fraction.assign(count, 0.0);
        output.orientation_bias_probability.assign(count, 0.0);
        output.somatic_log10_evidence.assign(count,
            -std::numeric_limits<double>::infinity());
        output.germline_log10_evidence.assign(count,
            -std::numeric_limits<double>::infinity());
        output.artifact_log10_evidence.assign(count,
            -std::numeric_limits<double>::infinity());
        output.informative_reads.assign(count, 0);
        std::vector<int> ownership(count, -1);
        // Region results are flattened by the exact EventMap allele key.  A
        // linear search here repeats the full candidate scan for every local
        // region candidate and is particularly expensive on a dense whole
        // chromosome HCC1143 traversal. `emplace` preserves the first equal
        // candidate, exactly matching the former `find_if` behavior when a
        // compatibility caller supplies duplicate rows.
        using MergedAlleleKey =
            std::tuple<std::int32_t, std::int32_t, std::string, std::string>;
        const auto allele_key = [](const auto& candidate) {
            return MergedAlleleKey{candidate.tid, candidate.position,
                                   allele_ref(candidate), allele_alt(candidate)};
        };
        std::map<MergedAlleleKey, std::size_t> merged_candidate_index;
        for (std::size_t index = 0; index < tumor.candidates.size(); ++index)
            merged_candidate_index.emplace(allele_key(tumor.candidates[index]), index);
        if (debug_somatic_region_reduce) {
            std::cerr << "[FASTGATK_POSTERIOR_REGION_REDUCE] stage=begin regions="
                      << tumor.assembly_region_likelihood_results.size()
                      << " candidates=" << tumor.candidates.size() << '\n';
        }
        for (std::size_t part_index = 0;
             part_index < tumor.assembly_region_likelihood_results.size(); ++part_index) {
            const auto& tumor_part_pointer = tumor.assembly_region_likelihood_results[part_index];
            if (!tumor_part_pointer) continue;
            const fastgatk::calling::Result* normal_part = nullptr;
            if (normal != nullptr &&
                part_index < normal->assembly_region_likelihood_results.size() &&
                normal->assembly_region_likelihood_results[part_index])
                normal_part = normal->assembly_region_likelihood_results[part_index].get();
            const auto local = calculate_posterior(
                options, *tumor_part_pointer, tumor_reads, normal_part);
            if (!local.has_value()) continue;
            output.prepare_seconds += local->prepare_seconds;
            output.seconds += local->seconds;
            if (output.execution_space.empty())
                output.execution_space = local->execution_space;
            const auto& part = *tumor_part_pointer;
            for (std::size_t local_index = 0; local_index < part.candidates.size(); ++local_index) {
                if (local_index >= local->somatic_probability.size()) continue;
                const auto merged = merged_candidate_index.find(
                    allele_key(part.candidates[local_index]));
                if (merged == merged_candidate_index.end()) continue;
                const auto merged_index = merged->second;
                int priority = 1;
                if (part.calling_regions.size() == 1) {
                    const auto& region = part.calling_regions.front();
                    if (part.candidates[local_index].tid == region.tid &&
                        part.candidates[local_index].position >= region.active_start &&
                        part.candidates[local_index].position <= region.active_end)
                        priority = 2;
                }
                if (priority < ownership[merged_index]) continue;
                ownership[merged_index] = priority;
                output.somatic_probability[merged_index] = local->somatic_probability[local_index];
                if (local_index < local->germline_probability.size())
                    output.germline_probability[merged_index] = local->germline_probability[local_index];
                if (local_index < local->artifact_probability.size())
                    output.artifact_probability[merged_index] = local->artifact_probability[local_index];
                if (local_index < local->contamination_adjusted_allele_fraction.size())
                    output.contamination_adjusted_allele_fraction[merged_index] =
                        local->contamination_adjusted_allele_fraction[local_index];
                if (local_index < local->orientation_bias_probability.size())
                    output.orientation_bias_probability[merged_index] =
                        local->orientation_bias_probability[local_index];
                if (local_index < local->somatic_log10_evidence.size())
                    output.somatic_log10_evidence[merged_index] =
                        local->somatic_log10_evidence[local_index];
                if (local_index < local->germline_log10_evidence.size())
                    output.germline_log10_evidence[merged_index] =
                        local->germline_log10_evidence[local_index];
                if (local_index < local->artifact_log10_evidence.size())
                    output.artifact_log10_evidence[merged_index] =
                        local->artifact_log10_evidence[local_index];
                if (local_index < local->informative_reads.size())
                    output.informative_reads[merged_index] = local->informative_reads[local_index];
            }
            if (debug_somatic_region_reduce &&
                ((part_index + 1U) % 128U == 0U ||
                 part_index + 1U == tumor.assembly_region_likelihood_results.size())) {
                std::cerr << "[FASTGATK_POSTERIOR_REGION_REDUCE] stage=progress regions="
                          << (part_index + 1U) << '/' << tumor.assembly_region_likelihood_results.size()
                          << '\n';
            }
        }
        if (debug_somatic_region_reduce) {
            std::cerr << "[FASTGATK_POSTERIOR_REGION_REDUCE] stage=end regions="
                      << tumor.assembly_region_likelihood_results.size()
                      << " candidates=" << tumor.candidates.size() << '\n';
        }
        return output;
    }
    if (tumor.candidates.empty() || tumor.reference_read_likelihoods.size() != tumor.candidates.size() ||
        tumor.allele_read_likelihoods.size() != tumor.candidates.size())
        return std::nullopt;
    const auto read_count = tumor.reference_read_likelihoods.front().size();
    if (read_count == 0) return std::nullopt;
    const auto candidate_count = tumor.candidates.size();
    std::vector<std::uint32_t> f1r2(candidate_count, 0);
    std::vector<std::uint32_t> r1f2(candidate_count, 0);
    for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
        const auto orientation = orientation_counts(tumor, tumor_reads, candidate);
        f1r2[candidate] = orientation.f1r2;
        r1f2[candidate] = orientation.r1f2;
    }

    std::vector<double> normal_reference;
    std::vector<double> normal_alternate;
    constexpr double missing = -1.0e300;
    if (normal != nullptr) {
        // Tumor and normal are separate AlleleLikelihoods collections.  A
        // normal commonly has a different number of retained fragments than
        // the tumor (for example, after sample selection or independent
        // coverage).  Preserve that native normal width so the Kokkos
        // posterior consumes it instead of dropping all normal evidence when
        // it does not happen to equal the tumor width.
        std::size_t normal_read_count = 0;
        for (const auto& row : normal->reference_read_likelihoods) {
            if (!row.empty()) {
                normal_read_count = row.size();
                break;
            }
        }
        if (normal_read_count == 0) {
            for (const auto& row : normal->allele_read_likelihoods) {
                if (!row.empty()) {
                    normal_read_count = row.size();
                    break;
                }
            }
        }
        if (normal_read_count > 0)
            normal_reference.reserve(candidate_count * normal_read_count);
        if (normal_read_count > 0)
            normal_alternate.reserve(candidate_count * normal_read_count);
        for (const auto& candidate : tumor.candidates) {
            const auto index = candidate_index(*normal, candidate);
            if (index == std::numeric_limits<std::size_t>::max() ||
                index >= normal->reference_read_likelihoods.size() ||
                index >= normal->allele_read_likelihoods.size() ||
                normal_read_count == 0 ||
                normal->reference_read_likelihoods[index].size() != normal_read_count ||
                normal->allele_read_likelihoods[index].size() != normal_read_count) {
                normal_reference.insert(normal_reference.end(), normal_read_count, missing);
                normal_alternate.insert(normal_alternate.end(), normal_read_count, missing);
                continue;
            }
            normal_reference.insert(normal_reference.end(),
                                    normal->reference_read_likelihoods[index].begin(),
                                    normal->reference_read_likelihoods[index].end());
            normal_alternate.insert(normal_alternate.end(),
                                    normal->allele_read_likelihoods[index].begin(),
                                    normal->allele_read_likelihoods[index].end());
        }
    }
    const auto flatten_with_fallback = [&](const auto& rows, const bool alternate) {
        std::vector<double> values;
        values.reserve(candidate_count * read_count);
        for (std::size_t candidate = 0; candidate < candidate_count; ++candidate) {
            const auto& row = rows[candidate];
            bool observed = false;
            for (const auto value : row) {
                if (std::isfinite(value) && value > -1.0e299) {
                    observed = true;
                    break;
                }
            }
            if (observed) {
                values.insert(values.end(), row.begin(), row.end());
                continue;
            }
            const auto fallback = candidate < tumor.likelihoods.size()
                ? tumor.likelihoods[candidate] : fastgatk::calling::Likelihoods{};
            values.push_back(alternate ? fallback.hom_alt : fallback.hom_ref);
            values.insert(values.end(), read_count > 1 ? read_count - 1 : 0, -1.0e300);
        }
        return values;
    };
    return fastgatk::kernels::calculate_somatic_posterior_kokkos(
        flatten_with_fallback(tumor.reference_read_likelihoods, false),
        flatten_with_fallback(tumor.allele_read_likelihoods, true),
        normal_reference, normal_alternate, f1r2, r1f2, candidate_count, read_count,
        options.contamination, options.somatic_prior, options.germline_prior,
        options.artifact_prior);
}

using Key = std::pair<std::int32_t, std::int32_t>;

// A normal sample keeps its own Host evidence/result for VCF serialization,
// while Mutect2's somatic model also receives one aggregate normal result
// (the same combined-likelihood boundary used by GATK's
// SomaticGenotypingEngine).  Pointers are invocation-local and never escape
// the writer call or a streaming tile.
struct NormalSampleView {
    std::string name;
    const fastgatk::io::ReadBatch* reads = nullptr;
    const fastgatk::calling::Result* result = nullptr;
    const fastgatk::kernels::SomaticLikelihoodResult* somatic = nullptr;
};

// GATK's TandemRepeat annotation is evaluated at the VCF boundary from
// normalized alleles plus forward reference context.  It deliberately does
// not participate in graph traversal or any Kokkos numerical kernel.
struct TandemRepeatAnnotation {
    std::vector<std::size_t> copies_per_allele;
    std::string repeat_unit;
};

std::size_t smallest_repeating_prefix_length(const std::string_view bases) {
    for (std::size_t length = 1; length <= bases.size(); ++length) {
        bool repeated = true;
        for (std::size_t start = length; start < bases.size(); start += length) {
            // Arrays.copyOfRange() in GATK pads a short final slice, so it
            // cannot equal the full candidate unit.  Require a whole unit
            // here rather than accepting a matching prefix of one.
            if (start + length > bases.size() ||
                bases.substr(0, length) != bases.substr(start, length)) {
                repeated = false;
                break;
            }
        }
        if (repeated) return length;
    }
    return bases.size();
}

std::size_t leading_repeat_count(const std::string_view unit,
                                 const std::string_view bases) {
    if (unit.empty()) return 0;
    std::size_t count = 0;
    for (std::size_t start = 0; start + unit.size() <= bases.size();
         start += unit.size()) {
        if (bases.substr(start, unit.size()) != unit) break;
        ++count;
    }
    return count;
}

std::optional<TandemRepeatAnnotation> tandem_repeat_annotation(
    const std::string& reference_allele,
    const std::vector<std::string>& alternate_alleles,
    const std::string_view reference_after_anchor) {
    if (reference_allele.empty() || alternate_alleles.empty()) return std::nullopt;
    if (std::all_of(alternate_alleles.begin(), alternate_alleles.end(),
                    [&](const auto& alternate) {
                        return alternate.size() == reference_allele.size();
                    })) return std::nullopt;
    const auto unpadded_reference = std::string_view(reference_allele).substr(1);
    TandemRepeatAnnotation annotation;
    for (const auto& alternate_allele : alternate_alleles) {
        if (alternate_allele.empty()) return std::nullopt;
        const auto unpadded_alternate = std::string_view(alternate_allele).substr(1);
        const auto longer = unpadded_alternate.size() > unpadded_reference.size()
            ? unpadded_alternate : unpadded_reference;
        if (longer.empty()) return std::nullopt;
        const auto repeat_unit = longer.substr(0, smallest_repeating_prefix_length(longer));
        const auto repeats_inside_reference = leading_repeat_count(repeat_unit, unpadded_reference);
        std::string reference_with_context(unpadded_reference);
        reference_with_context.append(reference_after_anchor);
        std::string alternate_with_context(unpadded_alternate);
        alternate_with_context.append(reference_after_anchor);
        const auto reference_repeats = leading_repeat_count(repeat_unit, reference_with_context);
        const auto alternate_repeats = leading_repeat_count(repeat_unit, alternate_with_context);
        if (reference_repeats <= repeats_inside_reference ||
            alternate_repeats <= repeats_inside_reference) return std::nullopt;
        if (annotation.copies_per_allele.empty())
            annotation.copies_per_allele.push_back(reference_repeats - repeats_inside_reference);
        annotation.copies_per_allele.push_back(alternate_repeats - repeats_inside_reference);
        annotation.repeat_unit.assign(repeat_unit);
    }
    return annotation;
}

std::string vcf_text(const Options& options, const fastgatk::io::HeaderSummary& header,
                     const fastgatk::io::ReadBatch& tumor_reads,
                     const fastgatk::calling::Result& tumor,
                     const fastgatk::calling::Result* normal,
                     const fastgatk::kernels::SomaticLikelihoodResult* tumor_somatic,
                     const fastgatk::kernels::SomaticLikelihoodResult* normal_somatic,
                     const fastgatk::kernels::SomaticPosteriorResult* posterior,
                     const std::vector<NormalSampleView>* normal_views = nullptr,
                     const std::vector<NormalSampleView>* tumor_views = nullptr,
                     const std::vector<std::string>* references = nullptr,
                     const MutectFeatureResources* feature_resources_input = nullptr) {
    using GroupKey = std::tuple<std::int32_t, std::int32_t, std::string>;
    using AlleleKey = std::tuple<std::int32_t, std::int32_t, std::string, std::string>;
    const MutectFeatureResources empty_feature_resources;
    const auto& feature_resources = feature_resources_input == nullptr
        ? empty_feature_resources : *feature_resources_input;
    std::vector<NormalSampleView> fallback_normal_views;
    if (normal_views != nullptr && !normal_views->empty()) {
        fallback_normal_views = *normal_views;
    } else if (normal != nullptr) {
        fallback_normal_views.push_back(NormalSampleView{
            options.normal_sample.empty() ? std::string("NORMAL") : options.normal_sample,
            nullptr, normal, normal_somatic});
    }
    const auto& output_normal_views = fallback_normal_views;
    const auto tumor_name = options.tumor_sample.empty()
        ? (header.samples.empty() ? std::string("TUMOR") : header.samples.front())
        : options.tumor_sample;
    std::vector<NormalSampleView> fallback_tumor_views;
    if (tumor_views != nullptr && !tumor_views->empty()) {
        fallback_tumor_views = *tumor_views;
    } else {
        fallback_tumor_views.push_back(NormalSampleView{tumor_name, &tumor_reads,
                                                         &tumor, tumor_somatic});
    }
    const auto& output_tumor_views = fallback_tumor_views;
    const bool multi_tumor_output = output_tumor_views.size() > 1;
    std::map<GroupKey, std::vector<const fastgatk::calling::GenotypeCall*>> tumor_groups;
    for (const auto& call : tumor.calls) {
        tumor_groups[GroupKey{call.candidate.tid, call.candidate.position,
                              allele_ref(call.candidate)}].push_back(&call);
    }
    // getVariantsFromActiveHaplotypes() discovers each EventMap event while
    // walking the already sorted haplotype list.  simpleMerge() preserves
    // that first-observed order in the final ALT vector.  Alphabetically
    // sorting the strings is tempting but changes every Number=A/R field at
    // a mixed SNP/indel locus (and, in particular, cannot reproduce a GATK
    // multiallelic tumor GT).  PairHMM retains the compact, region-local
    // membership of each EventMap allele precisely for this Host-side
    // serialization decision.
    const auto first_eventmap_haplotype = [](
        const fastgatk::calling::Result& result,
        const fastgatk::calling::AssemblyCandidate& candidate) {
        const auto missing = std::numeric_limits<std::uint32_t>::max();
        const auto matching_index = [&](const auto& source)
            -> std::optional<std::size_t> {
            for (std::size_t index = 0; index < source.candidates.size(); ++index) {
                const auto& item = source.candidates[index];
                if (item.tid == candidate.tid && item.position == candidate.position &&
                    allele_ref(item) == allele_ref(candidate) &&
                    allele_alt(item) == allele_alt(candidate))
                    return index;
            }
            return std::nullopt;
        };
        const auto rank_in = [&](const auto& source) {
            const auto index = matching_index(source);
            if (!index.has_value() || *index >= source.somatic_candidate_haplotype_indices.size() ||
                source.somatic_candidate_haplotype_indices[*index].empty())
                return missing;
            return *std::min_element(
                source.somatic_candidate_haplotype_indices[*index].begin(),
                source.somatic_candidate_haplotype_indices[*index].end());
        };
        auto rank = rank_in(result);
        if (rank != missing) return rank;
        // A flattened walker result intentionally owns no global haplotype
        // numbering.  Every call belongs to one AssemblyRegion, whose local
        // membership is the exact ordering source used by GATK.
        for (const auto& part : result.assembly_region_likelihood_results) {
            if (!part) continue;
            rank = rank_in(*part);
            if (rank != missing) return rank;
        }
        return missing;
    };
    for (auto& [key, calls] : tumor_groups) {
        std::sort(calls.begin(), calls.end(), [&](const auto* left, const auto* right) {
            const auto left_rank = first_eventmap_haplotype(tumor, left->candidate);
            const auto right_rank = first_eventmap_haplotype(tumor, right->candidate);
            if (left_rank != right_rank) return left_rank < right_rank;
            const auto left_alt = allele_alt(left->candidate);
            const auto right_alt = allele_alt(right->candidate);
            if (left_alt != right_alt) return left_alt < right_alt;
            return left->candidate.alternate_count > right->candidate.alternate_count;
        });
    }
    std::map<AlleleKey, const fastgatk::calling::GenotypeCall*> normal_allele_calls;
    if (normal != nullptr) {
        for (const auto& call : normal->calls) {
            const auto key = AlleleKey{call.candidate.tid, call.candidate.position,
                                      allele_ref(call.candidate), allele_alt(call.candidate)};
            normal_allele_calls.emplace(key, &call);
        }
    }
    std::vector<std::map<AlleleKey, const fastgatk::calling::GenotypeCall*>>
        output_normal_allele_calls;
    output_normal_allele_calls.reserve(output_normal_views.size());
    for (const auto& view : output_normal_views) {
        std::map<AlleleKey, const fastgatk::calling::GenotypeCall*> calls;
        if (view.result != nullptr) {
            for (const auto& call : view.result->calls) {
                const auto key = AlleleKey{call.candidate.tid, call.candidate.position,
                                          allele_ref(call.candidate), allele_alt(call.candidate)};
                calls.emplace(key, &call);
            }
        }
        output_normal_allele_calls.push_back(std::move(calls));
    }
    std::vector<std::map<AlleleKey, const fastgatk::calling::GenotypeCall*>>
        output_tumor_allele_calls;
    if (multi_tumor_output) {
        output_tumor_allele_calls.reserve(output_tumor_views.size());
        for (const auto& view : output_tumor_views) {
            std::map<AlleleKey, const fastgatk::calling::GenotypeCall*> calls;
            if (view.result != nullptr) {
                for (const auto& call : view.result->calls) {
                    const auto key = AlleleKey{call.candidate.tid, call.candidate.position,
                                              allele_ref(call.candidate), allele_alt(call.candidate)};
                    calls.emplace(key, &call);
                }
            }
            output_tumor_allele_calls.push_back(std::move(calls));
        }
    }
    // Mirror the stable GATK 4.6 Mutect2 header definitions even when a
    // definition is not populated in this particular callset.  GATK exposes
    // the complete writer schema to downstream tools, and omitting unused
    // fields changes the VCF interface despite identical data rows.
    std::ostringstream output;
    output << "##fileformat=VCFv4.2\n"
           << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
           << "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions of alternate alleles in the tumor\">\n"
           << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n"
           << "##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description=\"Count of reads in F1R2 pair orientation supporting each allele\">\n"
           << "##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description=\"Count of reads in F2R1 pair orientation supporting each allele\">\n"
           << "##FORMAT=<ID=FAD,Number=R,Type=Integer,Description=\"Count of fragments supporting each allele.\">\n"
           << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
           << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
           << "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles\">\n"
           << "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">\n"
           << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n"
           << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">\n"
           << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">\n";
    // This deliberately records the native executable rather than claiming
    // the Java command was run.  It has the same GATK structured-header
    // contract and is omitted, together with ##source, by the common switch.
    if (options.add_output_vcf_command_line)
        output << "##GATKCommandLine=<ID=Mutect2,Version=fastgatk-native,"
                  "CommandLine=\"fastgatk Mutect2\">\n";
    output << "##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.\">\n"
           << "##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description=\"Number of reads with unique start and mate end positions for each alt at a variant site\">\n"
           << "##INFO=<ID=CONTQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to contamination\">\n"
           << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n"
           << "##INFO=<ID=ECNT,Number=1,Type=Integer,Description=\"Number of potential somatic events in the assembly region\">\n"
           << "##INFO=<ID=ECNTH,Number=A,Type=Integer,Description=\"Number of somatic events in best supporting haplotype for each alt allele\">\n"
           << "##INFO=<ID=GERMQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not germline variants\">\n"
           << "##INFO=<ID=MBQ,Number=R,Type=Integer,Description=\"median base quality by allele\">\n"
           << "##INFO=<ID=MFRL,Number=R,Type=Integer,Description=\"median fragment length by allele\">\n"
           << "##INFO=<ID=MMQ,Number=R,Type=Integer,Description=\"median mapping quality by allele\">\n"
           << "##INFO=<ID=MPOS,Number=A,Type=Integer,Description=\"median distance from end of read\">\n"
           << "##INFO=<ID=NALOD,Number=A,Type=Float,Description=\"Log 10 odds of artifact in normal with same allele fraction as tumor\">\n"
           << "##INFO=<ID=NCount,Number=1,Type=Integer,Description=\"Count of N bases in the pileup\">\n"
           << "##INFO=<ID=NLOD,Number=A,Type=Float,Description=\"Normal log 10 likelihood ratio of diploid het or hom alt genotypes\">\n"
           << "##INFO=<ID=OCM,Number=1,Type=Integer,Description=\"Number of alt reads whose original alignment doesn't match the current contig.\">\n"
           << "##INFO=<ID=PON,Number=0,Type=Flag,Description=\"site found in panel of normals\">\n"
           << "##INFO=<ID=POPAF,Number=A,Type=Float,Description=\"negative log 10 population allele frequencies of alt alleles\">\n"
           << "##INFO=<ID=ROQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to read orientation artifact\">\n"
           << "##INFO=<ID=RPA,Number=R,Type=Integer,Description=\"Number of times tandem repeat unit is repeated, for each allele (including reference)\">\n"
           << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Tandem repeat unit (bases)\">\n"
           << "##INFO=<ID=SEQQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not sequencing errors\">\n"
           << "##INFO=<ID=STR,Number=0,Type=Flag,Description=\"Variant is a short tandem repeat\">\n"
           << "##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description=\"Phred-scaled quality of strand bias artifact\">\n"
           << "##INFO=<ID=STRQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors\">\n"
           << "##INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">\n"
           << "##MutectVersion=2.2\n";
    for (std::size_t index = 0; index < header.contigs.size(); ++index)
        output << "##contig=<ID=" << header.contigs[index] << ",length=" << header.contig_lengths[index]
               << (index < header.contig_assemblies.size() && !header.contig_assemblies[index].empty()
                       ? ",assembly=" + header.contig_assemblies[index] : std::string{})
               << ">\n";
    output << "##filtering_status=Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives.\n";
    // SomaticVariantOutputVCFWriter keeps the selected normal names as
    // ordinary metadata lines.  This is not redundant with #CHROM: tools
    // that consume a sites-only callset still need to know which role was
    // the matched normal.
    std::vector<std::string> normal_header_names;
    normal_header_names.reserve(output_normal_views.size());
    for (const auto& view : output_normal_views) {
        if (!view.name.empty()) normal_header_names.push_back(view.name);
    }
    std::sort(normal_header_names.begin(), normal_header_names.end());
    normal_header_names.erase(
        std::unique(normal_header_names.begin(), normal_header_names.end()),
        normal_header_names.end());
    for (const auto& normal_name : normal_header_names)
        output << "##normal_sample=" << normal_name << "\n";
    if (options.add_output_vcf_command_line)
        output << "##source=Mutect2\n";
    std::vector<std::string> tumor_header_names;
    tumor_header_names.reserve(output_tumor_views.size());
    for (const auto& view : output_tumor_views) {
        if (!view.name.empty()) tumor_header_names.push_back(view.name);
    }
    std::sort(tumor_header_names.begin(), tumor_header_names.end());
    tumor_header_names.erase(
        std::unique(tumor_header_names.begin(), tumor_header_names.end()),
        tumor_header_names.end());
    for (const auto& name : tumor_header_names)
        output << "##tumor_sample=" << name << "\n";
    struct OutputSampleOrder {
        std::string name;
        bool is_tumor = false;
        std::size_t format_index = std::numeric_limits<std::size_t>::max();
    };
    // Mutect2 constructs its SampleList from ReadUtils.getSamplesFromHeader(),
    // which is a TreeSet. VCF samples therefore use lexicographic SAM sample
    // order, not a fixed tumor-first order: HCC1143 sorts before HCC1143 BL,
    // while the DREAM normal sample sorts before its tumor. Keep the
    // evidence-role ownership separate from this writer order.
    std::vector<OutputSampleOrder> output_sample_order;
    output_sample_order.reserve(output_normal_views.size() + output_tumor_views.size());
    for (std::size_t index = 0; index < output_tumor_views.size(); ++index) {
        const auto& view = output_tumor_views[index];
        output_sample_order.push_back(OutputSampleOrder{
            view.name.empty() ? std::string("TUMOR") : view.name, true, index});
    }
    for (std::size_t index = 0; index < output_normal_views.size(); ++index) {
        const auto& view = output_normal_views[index];
        output_sample_order.push_back(OutputSampleOrder{
            view.name.empty() ? std::string("NORMAL") : view.name, false,
            (multi_tumor_output ? output_tumor_views.size() : 0) + index});
    }
    std::sort(output_sample_order.begin(), output_sample_order.end(),
              [](const auto& left, const auto& right) {
        return left.name < right.name;
    });
    output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if (!options.sites_only_vcf_output) {
        output << "\tFORMAT";
        for (const auto& sample : output_sample_order)
            output << '\t' << sample.name;
    }
    output << '\n';
    const auto find_depth_row = [](const fastgatk::calling::Result& result,
                                   const GroupKey& key,
                                   const std::vector<const fastgatk::calling::GenotypeCall*>& calls)
        -> const fastgatk::calling::MultiallelicDepth* {
        if (calls.empty()) return nullptr;
        for (const auto& row : result.multiallelic_depths) {
            if (row.tid != std::get<0>(key) || row.position != std::get<1>(key) ||
                row.reference != std::get<2>(key) || row.alternates.empty() ||
                row.depths.size() != row.alternates.size() + 1)
                continue;
            bool matched = true;
            for (const auto* call : calls) {
                if (std::find(row.alternates.begin(), row.alternates.end(),
                              allele_alt(call->candidate)) == row.alternates.end()) {
                    matched = false;
                    break;
                }
            }
            if (matched) return &row;
        }
        return nullptr;
    };
    const auto join_strings = [](const std::vector<std::string>& values) {
        std::ostringstream joined;
        for (std::size_t index = 0; index < values.size(); ++index) {
            if (index != 0) joined << ',';
            joined << values[index];
        }
        return joined.str();
    };
    // htsjdk VCFEncoder.formatVCFDouble(): INFO Double values below 0.01
    // use a three-decimal scientific representation (including negative
    // values), values below one use fixed three decimals, and all remaining
    // values use fixed two decimals.
    const auto vcf_double_text = [](const double value) {
        std::ostringstream text;
        if (value < 0.01) {
            if (std::abs(value) < 1.0e-20) return std::string("0.00");
            text << std::scientific << std::setprecision(3) << value;
        } else if (value < 1.0) {
            text << std::fixed << std::setprecision(3) << value;
        } else {
            text << std::fixed << std::setprecision(2) << value;
        }
        return text.str();
    };
    // htsjdk serializes the double[] carried by Mutect2's FORMAT/AF with
    // three decimal places at the VCF writer boundary.
    const auto format_af_text = [](const double value) {
        std::ostringstream text;
        text << std::fixed << std::setprecision(3) << value;
        return text.str();
    };
    // Mutect2 delegates physical phasing to
    // AssemblyBasedCallerUtils.phaseCalls(returnCalls, calledHaplotypes).
    // The phase relation is a relation between *sets of EventMap
    // haplotypes*, not a shared-read heuristic: a pair can be together on all
    // alternate haplotypes, apart on all alternate haplotypes, or genuinely
    // unphasable. PairHMM computes the read x haplotype scores in Kokkos; the
    // compact local EventMap memberships retained above let this Host writer
    // reproduce GATK's phaseCalls set state machine exactly.
    struct PhaseInfo {
        bool phased = false;
        std::string pgt;
        std::string pid;
        std::int32_t ps = 0;
    };
    std::map<GroupKey, PhaseInfo> phase_by_group;
    struct PhaseSite {
        GroupKey key;
        const fastgatk::calling::GenotypeCall* call = nullptr;
        const fastgatk::calling::Result* owner = nullptr;
        std::size_t candidate = std::numeric_limits<std::size_t>::max();
        std::set<std::uint32_t> haplotypes;
    };
    const auto phase_owner_for_candidate = [&](const fastgatk::calling::AssemblyCandidate& candidate) {
        const fastgatk::calling::Result* owner = &tumor;
        auto owner_index = candidate_index(tumor, candidate);
        int owner_priority = owner_index == std::numeric_limits<std::size_t>::max() ? -1 : 0;
        for (const auto& part_pointer : tumor.assembly_region_likelihood_results) {
            if (part_pointer == nullptr) continue;
            const auto& part = *part_pointer;
            const auto local_index = candidate_index(part, candidate);
            if (local_index == std::numeric_limits<std::size_t>::max()) continue;
            int priority = 1;
            for (const auto& region : part.calling_regions) {
                if (candidate.tid == region.tid && candidate.position >= region.active_start &&
                    candidate.position <= region.active_end) {
                    priority = 2;
                    break;
                }
            }
            if (priority > owner_priority) {
                owner = &part;
                owner_index = local_index;
                owner_priority = priority;
            }
        }
        return std::make_pair(owner, owner_index);
    };
    std::map<const fastgatk::calling::Result*, std::vector<PhaseSite>> phase_sites_by_owner;
    for (const auto& [key, calls] : tumor_groups) {
        if (calls.size() != 1 || calls.front() == nullptr) continue;
        const auto [owner, candidate] = phase_owner_for_candidate(calls.front()->candidate);
        if (owner == nullptr || candidate >= owner->somatic_candidate_haplotype_indices.size())
            continue;
        PhaseSite site;
        site.key = key;
        site.call = calls.front();
        site.owner = owner;
        site.candidate = candidate;
        const auto& memberships = owner->somatic_candidate_haplotype_indices[candidate];
        site.haplotypes.insert(memberships.begin(), memberships.end());
        if (!site.haplotypes.empty())
            phase_sites_by_owner[owner].push_back(std::move(site));
    }
    struct PhaseAssignment {
        int group = -1;
        bool phase_01 = true;
    };
    for (auto& [owner, sites] : phase_sites_by_owner) {
        (void)owner;
        if (sites.size() < 2) continue;
        std::sort(sites.begin(), sites.end(), [](const auto& left, const auto& right) {
            return left.key < right.key;
        });
        std::set<std::uint32_t> all_alternate_haplotypes;
        for (const auto& site : sites)
            all_alternate_haplotypes.insert(site.haplotypes.begin(), site.haplotypes.end());
        if (all_alternate_haplotypes.empty()) continue;
        std::vector<std::optional<PhaseAssignment>> assignments(sites.size());
        int unique_counter = 0;
        bool unphasable = false;
        for (std::size_t left = 0; left + 1 < sites.size() && !unphasable; ++left) {
            const auto& left_haplotypes = sites[left].haplotypes;
            bool left_is_on_all = left_haplotypes.size() == all_alternate_haplotypes.size();
            auto available_for_left = left_haplotypes;
            for (std::size_t right = left + 1; right < sites.size(); ++right) {
                const auto& right_haplotypes = sites[right].haplotypes;
                const bool right_is_on_all =
                    right_haplotypes.size() == all_alternate_haplotypes.size();
                const bool together =
                    (left_haplotypes.size() == right_haplotypes.size() &&
                     std::includes(left_haplotypes.begin(), left_haplotypes.end(),
                                   right_haplotypes.begin(), right_haplotypes.end())) ||
                    (left_is_on_all && std::includes(available_for_left.begin(),
                                                      available_for_left.end(),
                                                      right_haplotypes.begin(),
                                                      right_haplotypes.end())) ||
                    right_is_on_all;
                if (together) {
                    if (!assignments[left].has_value()) {
                        if (assignments[right].has_value()) {
                            unphasable = true;
                            break;
                        }
                        assignments[left] = PhaseAssignment{unique_counter, true};
                        assignments[right] = PhaseAssignment{unique_counter, true};
                        if (left_is_on_all && !right_is_on_all) {
                            std::set<std::uint32_t> intersection;
                            std::set_intersection(available_for_left.begin(), available_for_left.end(),
                                                  right_haplotypes.begin(), right_haplotypes.end(),
                                                  std::inserter(intersection, intersection.begin()));
                            available_for_left = std::move(intersection);
                        }
                        ++unique_counter;
                    } else if (!assignments[right].has_value()) {
                        assignments[right] = PhaseAssignment{
                            assignments[left]->group, assignments[left]->phase_01};
                    }
                    continue;
                }
                if (left_haplotypes.size() + right_haplotypes.size() !=
                    all_alternate_haplotypes.size())
                    continue;
                std::vector<std::uint32_t> intersection;
                std::set_intersection(left_haplotypes.begin(), left_haplotypes.end(),
                                      right_haplotypes.begin(), right_haplotypes.end(),
                                      std::back_inserter(intersection));
                if (!intersection.empty()) continue;
                if (!assignments[left].has_value()) {
                    if (assignments[right].has_value()) {
                        unphasable = true;
                        break;
                    }
                    assignments[left] = PhaseAssignment{unique_counter, true};
                    assignments[right] = PhaseAssignment{unique_counter, false};
                    ++unique_counter;
                } else if (!assignments[right].has_value()) {
                    assignments[right] = PhaseAssignment{
                        assignments[left]->group, !assignments[left]->phase_01};
                }
            }
        }
        if (unphasable) continue;
        std::map<int, std::vector<std::size_t>> members_by_group;
        for (std::size_t index = 0; index < assignments.size(); ++index) {
            if (assignments[index].has_value())
                members_by_group[assignments[index]->group].push_back(index);
        }
        for (const auto& [group, members] : members_by_group) {
            (void)group;
            if (members.size() < 2) continue;
            const auto first = members.front();
            const auto phase_set = std::get<1>(sites[first].key) + 1;
            const auto phase_id = std::to_string(phase_set) + "_" +
                allele_ref(sites[first].call->candidate) + "_" +
                allele_alt(sites[first].call->candidate);
            for (const auto member : members) {
                phase_by_group[sites[member].key] = PhaseInfo{
                    true, assignments[member]->phase_01 ? "0|1" : "1|0", phase_id, phase_set};
            }
        }
    }

    // SomaticGenotypingEngine computes ECNT from the EventMap events that
    // pass its tumor/normal emission predicates in this AssemblyRegion, not
    // from the number of ALT alleles in one serialized VCF row. ECNTH then
    // counts that same potential-event set on the called ALT's best
    // read-supported haplotype. Resolve the owning region before flattening so
    // overlapping padded halos cannot turn this into a chromosome-wide count.
    const auto event_map_owns_candidate = [](const fastgatk::calling::Result& source,
                                             const fastgatk::calling::AssemblyCandidate& candidate) {
        const auto matches = [&](const fastgatk::calling::ForcedAllele& event) {
            return event.tid == candidate.tid && event.position == candidate.position &&
                event.reference == allele_ref(candidate) && event.alternate == allele_alt(candidate);
        };
        if (!source.graph_haplotype_event_maps.empty()) {
            return std::any_of(source.graph_haplotype_event_maps.begin(),
                               source.graph_haplotype_event_maps.end(),
                               [&](const auto& event_map) {
                return std::any_of(event_map.events.begin(), event_map.events.end(), matches);
            });
        }
        // Older/no-graph compatibility inputs have no EventMap vector. Their
        // explicit graph-derived candidate marker is the only available
        // ownership proof; production Mutect2 region results take the path
        // above.
        return candidate.graph_derived;
    };
    const auto owner_for_candidate = [&](const fastgatk::calling::AssemblyCandidate& candidate) {
        const fastgatk::calling::Result* owner = &tumor;
        auto owner_index = candidate_index(tumor, candidate);
        int owner_priority = owner_index == std::numeric_limits<std::size_t>::max() ? -1 : 0;
        for (const auto& part_pointer : tumor.assembly_region_likelihood_results) {
            if (part_pointer == nullptr) continue;
            const auto& part = *part_pointer;
            const auto local_index = candidate_index(part, candidate);
            if (local_index == std::numeric_limits<std::size_t>::max()) continue;
            int priority = 1;
            for (const auto& region : part.calling_regions) {
                if (candidate.tid == region.tid && candidate.position >= region.active_start &&
                    candidate.position <= region.active_end) {
                    priority = 2;
                    break;
                }
            }
            if (priority > owner_priority) {
                owner = &part;
                owner_index = local_index;
                owner_priority = priority;
            }
        }
        return std::make_pair(owner, owner_index);
    };
    std::map<const fastgatk::calling::Result*, std::vector<std::uint8_t>>
        potential_somatic_events_by_owner;
    const auto potential_somatic_events = [&](const fastgatk::calling::Result& owner)
        -> const std::vector<std::uint8_t>& {
        const auto cached = potential_somatic_events_by_owner.find(&owner);
        if (cached != potential_somatic_events_by_owner.end()) return cached->second;
        std::vector<std::uint8_t> potential(owner.candidates.size(), 0U);
        for (std::size_t local = 0; local < owner.candidates.size(); ++local) {
            const auto& candidate = owner.candidates[local];
            if (!event_map_owns_candidate(owner, candidate)) continue;
            bool in_active_window = owner.calling_regions.empty();
            for (const auto& region : owner.calling_regions) {
                if (candidate.tid == region.tid && candidate.position >= region.active_start &&
                    candidate.position <= region.active_end) {
                    in_active_window = true;
                    break;
                }
            }
            if (!in_active_window) continue;
            const auto tumor_index = candidate_index(tumor, candidate);
            if (tumor_somatic == nullptr || tumor_index >= tumor_somatic->tlod.size() ||
                !std::isfinite(tumor_somatic->tlod[tumor_index]) ||
                tumor_somatic->tlod[tumor_index] <= effective_tumor_lod_emit(options))
                continue;
            if (normal != nullptr) {
                const auto normal_index = candidate_index(*normal, candidate);
                // A shared AssemblyResult normally forces every tumor EventMap
                // event into the normal likelihood collection. Preserve the
                // established missing-row behavior for legacy inputs, but use
                // the source NLOD predicate whenever the row exists.
                if (normal_index != std::numeric_limits<std::size_t>::max()) {
                    double nlod = normal_index < normal->likelihoods.size()
                        ? normal->likelihoods[normal_index].hom_ref -
                            normal->likelihoods[normal_index].het
                        : -std::numeric_limits<double>::infinity();
                    if (normal_somatic != nullptr &&
                        normal_index < normal_somatic->normal_log10_odds.size())
                        nlod = normal_somatic->normal_log10_odds[normal_index];
                    if (!std::isfinite(nlod) || nlod <= options.normal_log10_odds) continue;
                }
            }
            potential[local] = 1U;
        }
        return potential_somatic_events_by_owner.emplace(&owner, std::move(potential)).first->second;
    };
    const auto event_count_annotation =
        [&](const fastgatk::calling::GenotypeCall& call) -> std::pair<std::size_t, std::size_t> {
        const auto [owner, target] = owner_for_candidate(call.candidate);
        if (owner == nullptr || target >= owner->candidates.size()) return {1U, 1U};
        const auto& potential = potential_somatic_events(*owner);
        const auto concrete_event_count = static_cast<std::size_t>(
            std::count(potential.begin(), potential.end(), static_cast<std::uint8_t>(1U)));
        // Reference-confidence calls add <NON_REF> to each merged EventMap
        // VariantContext *before* SomaticGenotypingEngine records its
        // potential-event set.  With the negative-infinity gVCF emission
        // threshold both the concrete and symbolic alternate satisfy that
        // predicate, so ECNT includes one symbolic event per concrete one.
        // ECNTH remains the concrete assembled-haplotype count below.
        const auto event_count = concrete_event_count *
            (options.reference_confidence_mode == ReferenceConfidenceMode::None ? 1U : 2U);
        if (target >= owner->somatic_best_haplotype_candidate_indices.size())
            return {event_count, event_count};
        std::size_t haplotype_event_count = 0;
        for (const auto event : owner->somatic_best_haplotype_candidate_indices[target]) {
            if (event < potential.size() && potential[event] != 0) ++haplotype_event_count;
        }
        // An EventMap event can be injected under a legacy sparse-likelihood
        // route that has no materialized best-haplotype vector. Keep ECNTH
        // defined in that narrow case rather than serialize an invalid Number=A
        // zero; normal source-aligned runs take the explicit vector above.
        if (haplotype_event_count == 0 && event_count != 0)
            haplotype_event_count = event_count;
        return {event_count, haplotype_event_count};
    };
    for (const auto& [group_key, calls] : tumor_groups) {
        if (calls.empty()) continue;
        const auto tid_value = std::get<0>(group_key);
        const auto tid = static_cast<std::size_t>(std::max<std::int32_t>(0, tid_value));
        const auto chrom = tid < header.contigs.size() ? header.contigs[tid] : std::to_string(tid_value);
        const auto ref = std::get<2>(group_key);
        const auto position = std::get<1>(group_key);
        const bool in_panel_of_normals = resource_has_locus(
            feature_resources.panel_of_normals, chrom, position);
        // The germline-resource/PON early veto is applied by the Host before
        // Kokkos ActivityProfile evaluation, exactly where Mutect2Engine
        // applies it. Do not repeat it in the writer: --force-active must be
        // able to traverse an inactive resource locus and emit its normal
        // POPAF/PON annotations, just as GATK does.
        const auto tumor_row = find_depth_row(tumor, group_key, calls);
        std::vector<std::uint32_t> alt_depths(calls.size(), 0);
        std::uint32_t ref_depth = calls.front()->candidate.reference_count;
        std::uint32_t depth = calls.front()->candidate.depth;
        if (tumor_row != nullptr) {
            ref_depth = tumor_row->depths.front();
            depth = tumor_row->depth;
            for (std::size_t index = 0; index < calls.size(); ++index) {
                const auto found = std::find(tumor_row->alternates.begin(), tumor_row->alternates.end(),
                                             allele_alt(calls[index]->candidate));
                if (found != tumor_row->alternates.end())
                    alt_depths[index] = tumor_row->depths[static_cast<std::size_t>(found - tumor_row->alternates.begin()) + 1];
            }
        } else {
            for (std::size_t index = 0; index < calls.size(); ++index) {
                alt_depths[index] = calls[index]->candidate.alternate_count;
                ref_depth = std::max(ref_depth, calls[index]->candidate.reference_count);
                depth = std::max(depth, calls[index]->candidate.depth);
            }
            std::uint64_t alt_sum = 0;
            for (const auto value : alt_depths) alt_sum += value;
            if (ref_depth == 0 && depth > alt_sum) ref_depth = depth - static_cast<std::uint32_t>(alt_sum);
        }
        std::uint64_t represented_depth = ref_depth;
        // GATK's DepthPerAlleleBySample sees only the alleles retained in
        // this VariantContext.  Candidate-level scalars above are useful
        // fallbacks for incomplete/direct Results, but a complete PairHMM
        // matrix must be reduced over the final emitted allele list here.
        std::vector<fastgatk::calling::AssemblyCandidate> emitted_candidates;
        emitted_candidates.reserve(calls.size());
        for (const auto* call : calls) emitted_candidates.push_back(call->candidate);
        const auto emitted_tumor_annotation_counts = somatic_annotation_counts(
            tumor, tumor_reads, &emitted_candidates);
        const auto emitted_tumor_fragment_depth = emitted_somatic_fragment_depth(
            tumor, emitted_candidates);
        std::vector<std::optional<SomaticAnnotationCounts>> emitted_normal_annotation_counts;
        emitted_normal_annotation_counts.reserve(output_normal_views.size());
        for (const auto& view : output_normal_views) {
            emitted_normal_annotation_counts.push_back(view.result != nullptr && view.reads != nullptr
                ? somatic_annotation_counts(*view.result, *view.reads, &emitted_candidates)
                : std::optional<SomaticAnnotationCounts>{});
        }
        std::vector<std::optional<SomaticAnnotationCounts>> emitted_output_tumor_annotation_counts;
        if (multi_tumor_output) {
            emitted_output_tumor_annotation_counts.reserve(output_tumor_views.size());
            for (const auto& view : output_tumor_views) {
                emitted_output_tumor_annotation_counts.push_back(
                    view.result != nullptr && view.reads != nullptr
                    ? somatic_annotation_counts(*view.result, *view.reads, &emitted_candidates)
                    : std::optional<SomaticAnnotationCounts>{});
            }
        }
        if (const auto emitted_depth = emitted_somatic_read_depth(
                tumor, tumor_reads, emitted_candidates);
            emitted_depth.has_value() && emitted_depth->allele_depths.size() == calls.size() + 1) {
            ref_depth = emitted_depth->allele_depths.front();
            for (std::size_t index = 0; index < calls.size(); ++index)
                alt_depths[index] = emitted_depth->allele_depths[index + 1];
            depth = emitted_depth->depth;
        }
        represented_depth = ref_depth;
        for (const auto value : alt_depths) represented_depth += value;
        depth = std::max<std::uint32_t>(1, std::max(depth,
            static_cast<std::uint32_t>(std::min<std::uint64_t>(represented_depth,
                std::numeric_limits<std::uint32_t>::max()))));

        std::vector<std::string> alternates;
        std::vector<std::string> tlods, nalods, nlods, format_afs;
        std::vector<std::string> format_f1r2s, format_f2r1s;
        const auto first_tumor_candidate = candidate_index(tumor, calls.front()->candidate);
        std::uint64_t coverage_depth = coverage_depth_for_candidate(
            tumor, first_tumor_candidate);
        for (const auto& view : output_normal_views) {
            if (view.result == nullptr) continue;
            const auto candidate = candidate_index(*view.result, calls.front()->candidate);
            if (candidate < view.result->candidates.size())
                coverage_depth += coverage_depth_for_candidate(*view.result, candidate);
        }
        const auto info_depth = static_cast<std::uint32_t>(std::min<std::uint64_t>(
            coverage_depth, std::numeric_limits<std::uint32_t>::max()));
        if (emitted_tumor_annotation_counts.has_value() &&
            first_tumor_candidate < emitted_tumor_annotation_counts->fragment_f1r2_reference.size()) {
            format_f1r2s.push_back(std::to_string(
                emitted_tumor_annotation_counts->fragment_f1r2_reference[first_tumor_candidate]));
            format_f2r1s.push_back(std::to_string(
                emitted_tumor_annotation_counts->fragment_f2r1_reference[first_tumor_candidate]));
        } else {
            format_f1r2s.push_back(".");
            format_f2r1s.push_back(".");
        }
        std::vector<std::uint8_t> genotypes;
        double site_somatic_probability = 0.0;
        std::size_t best_alt = 0;
        for (std::size_t index = 0; index < calls.size(); ++index) {
            const auto& call = *calls[index];
            alternates.push_back(allele_alt(call.candidate));
            const auto tumor_candidate = candidate_index(tumor, call.candidate);
            const auto normal_candidate = normal == nullptr
                ? std::numeric_limits<std::size_t>::max()
                : candidate_index(*normal, call.candidate);
            const bool tumor_has_somatic = tumor_somatic != nullptr &&
                tumor_candidate < tumor_somatic->tlod.size();
            const double tlod = tumor_has_somatic ? tumor_somatic->tlod[tumor_candidate] :
                call.likelihoods.hom_alt - call.likelihoods.hom_ref;
            const bool normal_has_somatic = normal_somatic != nullptr &&
                normal_candidate < normal_somatic->normal_log10_odds.size();
            const double nlod = normal_candidate == std::numeric_limits<std::size_t>::max()
                ? 0.0 : (normal_has_somatic
                    ? normal_somatic->normal_log10_odds[normal_candidate] :
                    (normal_candidate < normal->likelihoods.size()
                        ? normal->likelihoods[normal_candidate].hom_ref -
                            normal->likelihoods[normal_candidate].het
                        : 0.0));
            // SomaticGenotypingEngine writes NALOD from somaticLogOdds() on
            // the normal fragment likelihood matrix, independently of the
            // fixed-diploid comparison that produces NLOD.
            const double nalod = normal_candidate == std::numeric_limits<std::size_t>::max()
                ? 0.0 : (normal_somatic != nullptr && normal_candidate < normal_somatic->tlod.size()
                    ? normal_somatic->tlod[normal_candidate] : 0.0);
            const double af = tumor_somatic != nullptr &&
                tumor_candidate < tumor_somatic->best_allele_fraction.size()
                ? tumor_somatic->best_allele_fraction[tumor_candidate]
                : static_cast<double>(alt_depths[index]) / static_cast<double>(depth);
            const bool has_posterior = posterior != nullptr &&
                tumor_candidate < posterior->somatic_probability.size();
            const auto somatic_probability = has_posterior
                ? posterior->somatic_probability[tumor_candidate] : 0.0;
            tlods.push_back(vcf_double_text(tlod));
            nalods.push_back(vcf_double_text(nalod));
            nlods.push_back(vcf_double_text(nlod));
            format_afs.push_back(format_af_text(af));
            if (emitted_tumor_annotation_counts.has_value() &&
                tumor_candidate < emitted_tumor_annotation_counts->fragment_f1r2_alternate.size()) {
                format_f1r2s.push_back(std::to_string(
                    emitted_tumor_annotation_counts->fragment_f1r2_alternate[tumor_candidate]));
                format_f2r1s.push_back(std::to_string(
                    emitted_tumor_annotation_counts->fragment_f2r1_alternate[tumor_candidate]));
            } else {
                format_f1r2s.push_back(".");
                format_f2r1s.push_back(".");
            }
            genotypes.push_back(call.genotype);
            if (somatic_probability > site_somatic_probability ||
                (somatic_probability == site_somatic_probability &&
                 alt_depths[index] > alt_depths[best_alt])) {
                site_somatic_probability = somatic_probability;
                best_alt = index;
            }
        }
        std::optional<TandemRepeatAnnotation> tandem_repeat;
        if (references != nullptr && tid_value >= 0 &&
            static_cast<std::size_t>(tid_value) < references->size()) {
            const auto& reference_sequence = (*references)[static_cast<std::size_t>(tid_value)];
            const auto position = std::get<1>(group_key);
            if (position >= 0 && static_cast<std::size_t>(position) < reference_sequence.size()) {
                const auto after_anchor = static_cast<std::size_t>(position) + 1;
                tandem_repeat = tandem_repeat_annotation(ref, alternates,
                    std::string_view(reference_sequence).substr(after_anchor));
            }
        }
        std::string tumor_sb = ".";
        if (emitted_tumor_annotation_counts.has_value()) {
            const auto candidate = candidate_index(tumor, calls.front()->candidate);
            if (candidate < emitted_tumor_annotation_counts->reference_forward.size()) {
                std::uint32_t alt_forward = 0;
                std::uint32_t alt_reverse = 0;
                for (const auto* call : calls) {
                    const auto index = candidate_index(tumor, call->candidate);
                    if (index < emitted_tumor_annotation_counts->alternate_forward.size()) {
                        alt_forward += emitted_tumor_annotation_counts->alternate_forward[index];
                        alt_reverse += emitted_tumor_annotation_counts->alternate_reverse[index];
                    }
                }
                tumor_sb = std::to_string(emitted_tumor_annotation_counts->reference_forward[candidate]) + ',' +
                    std::to_string(emitted_tumor_annotation_counts->reference_reverse[candidate]) + ',' +
                    std::to_string(alt_forward) + ',' + std::to_string(alt_reverse);
            }
        }
        std::ostringstream popaf;
        for (std::size_t index = 0; index < calls.size(); ++index) {
            if (index != 0) popaf << ',';
            const auto population_af = resource_allele_frequency(
                feature_resources.germline, chrom, position, ref,
                allele_alt(calls[index]->candidate)).value_or(
                    options.population_allele_frequency);
            // SomaticGenotypingEngine initializes every emitted ALT with the
            // configured default, then replaces only resource-matched alleles
            // and writes -log10(AF) through HTSJDK's VCF double formatter.
            popaf << vcf_double_text(-std::log10(population_af));
        }
        // Mutect2 itself does not assign a site QUAL.  GATK's
        // SomaticVariantOutputVCFWriter emits a missing QUAL ('.') and
        // leaves final confidence calculation to FilterMutectCalls.  A
        // posterior-derived numeric value here is tempting, but it changes
        // the GATK output contract and can make downstream filtering treat a
        // pre-filter score as a calibrated QUAL.  Keep the posterior in the
        // INFO fields and preserve the writer boundary exactly.
        // Mutect2 emits an unfiltered callset. Germline/orientation and all
        // other filtering decisions belong to the separate
        // FilterMutectCalls tool, so this column remains missing here.
        const auto filter = ".";
        std::ostringstream genotype_text;
        const auto phase_iter = phase_by_group.find(group_key);
        const auto phased = phase_iter != phase_by_group.end() && phase_iter->second.phased;
        // SomaticGenotypingEngine does not emit a diploid germline genotype
        // for the tumor.  For a single concrete ALT its tumor genotype list
        // is the complete biallelic allele list, which HTSlib serializes as
        // 0/1 even when the likelihood winner is hom-alt.  Keep the native
        // PairHMM winner for QUAL/filtering, but preserve this GATK-shaped
        // tumor GT boundary instead of leaking a germline 1/1 call.
        if (calls.size() == 1) {
            // AssemblyBasedCallerUtils.phaseVC() changes the concrete GT
            // allele order as well as adding PGT/PID/PS.  In particular,
            // the second member of a 1|0 phase group must serialize GT=1|0,
            // not the pre-phasing 0|1 order.  `pgt` is the same PhaseGroup
            // description GATK uses for a biallelic tumor genotype.
            genotype_text << (phased ? phase_iter->second.pgt : "0/1");
        } else {
            // SomaticGenotypingEngine.addGenotypes() does not emit a
            // likelihood-derived diploid germline GT for a tumor sample. It
            // writes the complete local AlleleList instead: REF followed by
            // every emitted tumor ALT in EventMap merge order.  Preserve
            // that variable-ploidy serialization rather than choosing the
            // best biallelic genotype, which loses low-AF siblings and makes
            // an otherwise correct multiallelic record differ from GATK.
            genotype_text << '0';
            for (std::size_t allele = 1; allele <= calls.size(); ++allele)
                genotype_text << '/' << allele;
        }
        std::vector<std::string> ad_text{std::to_string(ref_depth)};
        for (const auto value : alt_depths) ad_text.push_back(std::to_string(value));
        // FragmentDepthPerAlleleBySample independently annotates FAD from
        // grouped fragment BestAlleles.  It is not a duplicate spelling of
        // the read-level AD vector, even for biallelic sites.
        std::vector<std::string> fad_text = ad_text;
        if (emitted_tumor_fragment_depth.has_value() &&
            emitted_tumor_fragment_depth->allele_depths.size() == calls.size() + 1) {
            for (std::size_t allele = 0; allele < fad_text.size(); ++allele)
                fad_text[allele] = std::to_string(emitted_tumor_fragment_depth->allele_depths[allele]);
        } else {
            for (std::size_t index = 0; index < calls.size(); ++index) {
                const auto candidate = candidate_index(tumor, calls[index]->candidate);
                if (candidate >= tumor.somatic_fragment_reference_count.size() ||
                    candidate >= tumor.somatic_fragment_alternate_count.size())
                    continue;
                if (index == 0)
                    fad_text.front() = std::to_string(tumor.somatic_fragment_reference_count[candidate]);
                fad_text[index + 1] = std::to_string(tumor.somatic_fragment_alternate_count[candidate]);
            }
        }
        std::uint32_t as_reference_forward = 0;
        std::uint32_t as_reference_reverse = 0;
        std::vector<std::uint32_t> as_alternate_forward(calls.size(), 0);
        std::vector<std::uint32_t> as_alternate_reverse(calls.size(), 0);
        const auto add_allele_specific_strand = [&](const SomaticAnnotationCounts& annotation,
                                                     const std::vector<std::size_t>& candidates) {
            if (candidates.empty()) return;
            const auto reference_candidate = candidates.front();
            if (reference_candidate >= annotation.reference_forward.size() ||
                reference_candidate >= annotation.reference_reverse.size())
                return;
            std::uint64_t sample_total = annotation.reference_forward[reference_candidate] +
                annotation.reference_reverse[reference_candidate];
            for (std::size_t index = 0; index < candidates.size(); ++index) {
                const auto candidate = candidates[index];
                if (candidate >= annotation.alternate_forward.size() ||
                    candidate >= annotation.alternate_reverse.size())
                    continue;
                sample_total += annotation.alternate_forward[candidate] +
                    annotation.alternate_reverse[candidate];
            }
            // StrandBiasUtils excludes an entire sample table unless its
            // informative BestAllele population is strictly greater than 2.
            if (sample_total <= 2) return;
            as_reference_forward += annotation.reference_forward[reference_candidate];
            as_reference_reverse += annotation.reference_reverse[reference_candidate];
            for (std::size_t index = 0; index < candidates.size(); ++index) {
                const auto candidate = candidates[index];
                if (candidate >= annotation.alternate_forward.size() ||
                    candidate >= annotation.alternate_reverse.size())
                    continue;
                as_alternate_forward[index] += annotation.alternate_forward[candidate];
                as_alternate_reverse[index] += annotation.alternate_reverse[candidate];
            }
        };
        if (emitted_tumor_annotation_counts.has_value()) {
            std::vector<std::size_t> candidates;
            candidates.reserve(calls.size());
            for (const auto* call : calls)
                candidates.push_back(candidate_index(tumor, call->candidate));
            add_allele_specific_strand(*emitted_tumor_annotation_counts, candidates);
        }
        for (std::size_t normal_index = 0;
             normal_index < emitted_normal_annotation_counts.size(); ++normal_index) {
            const auto& annotation = emitted_normal_annotation_counts[normal_index];
            const auto* result = output_normal_views[normal_index].result;
            if (!annotation.has_value() || result == nullptr) continue;
            std::vector<std::size_t> candidates;
            candidates.reserve(calls.size());
            for (const auto* call : calls)
                candidates.push_back(candidate_index(*result, call->candidate));
            add_allele_specific_strand(*annotation, candidates);
        }
        std::ostringstream allele_specific_strand;
        allele_specific_strand << as_reference_forward << ',' << as_reference_reverse;
        for (std::size_t index = 0; index < calls.size(); ++index)
            allele_specific_strand << '|' << as_alternate_forward[index] << ','
                                  << as_alternate_reverse[index];
        PerAlleleAnnotationValues per_allele_reference;
        std::vector<PerAlleleAnnotationValues> per_allele_alternates(calls.size());
        const auto add_per_allele_annotation_values =
            [&](const SomaticAnnotationCounts& annotation,
                const std::vector<std::size_t>& candidates) {
                if (candidates.empty() ||
                    candidates.front() >= annotation.reference_annotation_values.size())
                    return;
                append_per_allele_annotation_values(per_allele_reference,
                    annotation.reference_annotation_values[candidates.front()]);
                for (std::size_t index = 0; index < candidates.size(); ++index) {
                    const auto candidate = candidates[index];
                    if (candidate < annotation.alternate_annotation_values.size())
                        append_per_allele_annotation_values(per_allele_alternates[index],
                            annotation.alternate_annotation_values[candidate]);
                }
            };
        if (emitted_tumor_annotation_counts.has_value()) {
            std::vector<std::size_t> candidates;
            candidates.reserve(calls.size());
            for (const auto* call : calls)
                candidates.push_back(candidate_index(tumor, call->candidate));
            add_per_allele_annotation_values(*emitted_tumor_annotation_counts, candidates);
        }
        for (std::size_t normal_index = 0;
             normal_index < emitted_normal_annotation_counts.size(); ++normal_index) {
            const auto& annotation = emitted_normal_annotation_counts[normal_index];
            const auto* result = output_normal_views[normal_index].result;
            if (!annotation.has_value() || result == nullptr) continue;
            std::vector<std::size_t> candidates;
            candidates.reserve(calls.size());
            for (const auto* call : calls)
                candidates.push_back(candidate_index(*result, call->candidate));
            add_per_allele_annotation_values(*annotation, candidates);
        }
        // ECNT is site-wide, while ECNTH is Number=A in the GATK header. Both
        // originate from the owner AssemblyRegion's EventMap rather than the
        // current VCF site's ALT cardinality.
        const auto first_event_count = event_count_annotation(*calls.front());
        std::vector<std::string> haplotype_event_counts;
        haplotype_event_counts.reserve(calls.size());
        for (const auto* call : calls)
            haplotype_event_counts.push_back(std::to_string(event_count_annotation(*call).second));
        output << chrom << '\t' << std::get<1>(group_key) + 1 << "\t.\t" << ref << '\t'
               << join_strings(alternates) << "\t.\t" << filter
               << "\tAS_SB_TABLE=" << allele_specific_strand.str()
               << ";DP=" << info_depth << ";ECNT=" << first_event_count.first << ";ECNTH="
               << join_strings(haplotype_event_counts)
               << ";MBQ=" << per_allele_annotation_values(
                   per_allele_reference, per_allele_alternates, 'b')
               << ";MFRL=" << per_allele_annotation_values(
                   per_allele_reference, per_allele_alternates, 'f')
               << ";MMQ=" << per_allele_annotation_values(
                   per_allele_reference, per_allele_alternates, 'm')
               << ";MPOS=" << alternate_per_allele_annotation_values(
                   per_allele_alternates, 'p');
        // SomaticGenotypingEngine only attaches normal-evidence annotations
        // when a matched normal is present.  Keep their header definitions
        // available, but do not manufacture zero-valued INFO fields for a
        // tumor-only callset.
        if (normal != nullptr)
            output << ";NALOD=" << join_strings(nalods) << ";NLOD=" << join_strings(nlods);
        if (in_panel_of_normals) output << ";PON";
        output << ";POPAF=" << popaf.str();
        if (tandem_repeat.has_value()) {
            std::vector<std::string> repeat_counts;
            repeat_counts.reserve(tandem_repeat->copies_per_allele.size());
            for (const auto count : tandem_repeat->copies_per_allele)
                repeat_counts.push_back(std::to_string(count));
            output << ";RPA=" << join_strings(repeat_counts)
                   << ";RU=" << tandem_repeat->repeat_unit << ";STR";
        }
        output << ";TLOD=" << join_strings(tlods);
        // The site-only mode is intentionally decided after all site INFO
        // values have been computed, so it cannot perturb TLOD/posterior or
        // F1R2 generation.  This matches GATK's writer option, which only
        // changes serialization of FORMAT/sample fields.
        if (options.sites_only_vcf_output) {
            output << '\n';
            continue;
        }
        const auto format = phased ? "GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB"
                                   : "GT:AD:AF:DP:F1R2:F2R1:FAD:SB";
        std::ostringstream tumor_sample_text;
        tumor_sample_text << genotype_text.str() << ':'
                          << join_strings(ad_text) << ':' << join_strings(format_afs) << ':' << depth << ':'
                          << join_strings(format_f1r2s) << ':' << join_strings(format_f2r1s) << ':'
                          << join_strings(fad_text);
        if (phased) {
            tumor_sample_text << ":" << phase_iter->second.pgt << ":" << phase_iter->second.pid
                              << ":" << phase_iter->second.ps;
        }
        tumor_sample_text << ':' << tumor_sb;
        struct FormatSampleView {
            const NormalSampleView* view = nullptr;
            bool is_tumor = false;
            std::size_t evidence_index = std::numeric_limits<std::size_t>::max();
        };
        std::vector<FormatSampleView> format_views;
        format_views.reserve(output_normal_views.size() +
                             (multi_tumor_output ? output_tumor_views.size() : 0));
        if (multi_tumor_output) {
            for (std::size_t index = 0; index < output_tumor_views.size(); ++index)
                format_views.push_back(FormatSampleView{&output_tumor_views[index], true, index});
        }
        for (std::size_t index = 0; index < output_normal_views.size(); ++index)
            format_views.push_back(FormatSampleView{&output_normal_views[index], false, index});
        std::vector<std::string> format_sample_texts;
        format_sample_texts.reserve(format_views.size());
        for (const auto& format_view : format_views) {
            const auto& view = *format_view.view;
            const auto* normal_view = view.result;
            const auto& normal_calls = format_view.is_tumor
                ? output_tumor_allele_calls[format_view.evidence_index]
                : output_normal_allele_calls[format_view.evidence_index];
            std::vector<const fastgatk::calling::GenotypeCall*> normal_group;
            normal_group.reserve(calls.size());
            std::vector<std::uint32_t> normal_alts(calls.size(), 0);
            std::vector<std::size_t> normal_candidate_indices(
                calls.size(), std::numeric_limits<std::size_t>::max());
            std::uint32_t normal_ref = 0, normal_depth = 0;
            for (std::size_t index = 0; index < calls.size(); ++index) {
                const auto iter = normal_calls.find(
                    AlleleKey{calls[index]->candidate.tid, calls[index]->candidate.position,
                              ref, allele_alt(calls[index]->candidate)});
                if (iter == normal_calls.end()) continue;
                const auto* normal_call = iter->second;
                normal_group.push_back(normal_call);
                normal_alts[index] = normal_call->candidate.alternate_count;
                normal_ref = std::max(normal_ref, normal_call->candidate.reference_count);
                normal_depth = std::max(normal_depth, normal_call->candidate.depth);
                if (normal_view != nullptr)
                    normal_candidate_indices[index] = candidate_index(
                        *normal_view, normal_call->candidate);
            }
            const auto normal_key = group_key;
            if (normal_view != nullptr) {
                if (const auto* row = find_depth_row(*normal_view, normal_key, normal_group);
                    row != nullptr && !normal_group.empty()) {
                    normal_ref = row->depths.front();
                    normal_depth = row->depth;
                    for (std::size_t index = 0; index < calls.size(); ++index) {
                        const auto found = std::find(row->alternates.begin(), row->alternates.end(),
                                                     allele_alt(calls[index]->candidate));
                        if (found != row->alternates.end())
                            normal_alts[index] = row->depths[static_cast<std::size_t>(found - row->alternates.begin()) + 1];
                    }
                }
            }
            const auto sample_emitted_read_depth = normal_view != nullptr && view.reads != nullptr
                ? emitted_somatic_read_depth(*normal_view, *view.reads, emitted_candidates)
                : std::optional<EmittedSomaticReadDepth>{};
            const auto sample_emitted_fragment_depth = normal_view != nullptr
                ? emitted_somatic_fragment_depth(*normal_view, emitted_candidates)
                : std::optional<EmittedSomaticFragmentDepth>{};
            if (sample_emitted_read_depth.has_value() &&
                sample_emitted_read_depth->allele_depths.size() == calls.size() + 1) {
                normal_ref = sample_emitted_read_depth->allele_depths.front();
                normal_depth = sample_emitted_read_depth->depth;
                for (std::size_t index = 0; index < calls.size(); ++index)
                    normal_alts[index] = sample_emitted_read_depth->allele_depths[index + 1];
            }
            std::vector<std::string> normal_ad{std::to_string(normal_ref)};
            std::vector<std::string> normal_af;
            for (std::size_t index = 0; index < normal_alts.size(); ++index) {
                const auto value = normal_alts[index];
                normal_ad.push_back(std::to_string(value));
                const auto normal_somatic_af = normal_view != nullptr && view.somatic != nullptr &&
                    normal_candidate_indices[index] < view.somatic->best_allele_fraction.size()
                    ? std::optional<double>(
                        view.somatic->best_allele_fraction[normal_candidate_indices[index]])
                    : std::optional<double>{};
                normal_af.push_back(normal_somatic_af.has_value()
                    ? format_af_text(*normal_somatic_af)
                    // SomaticGenotypingEngine evaluates even an empty
                    // matched-normal likelihood matrix.  Its flat allele
                    // pseudocount prior gives every allele equal AF rather
                    // than an undefined/zero value (0.500 for a biallelic
                    // no-read normal).
                    : (normal_depth == 0
                        ? format_af_text(1.0 / static_cast<double>(calls.size() + 1))
                        : format_af_text(
                        static_cast<double>(value) / static_cast<double>(normal_depth))));
            }
            std::vector<std::string> normal_fad = normal_ad;
            if (sample_emitted_fragment_depth.has_value() &&
                sample_emitted_fragment_depth->allele_depths.size() == calls.size() + 1) {
                for (std::size_t allele = 0; allele < normal_fad.size(); ++allele)
                    normal_fad[allele] = std::to_string(
                        sample_emitted_fragment_depth->allele_depths[allele]);
            } else if (normal_view != nullptr) {
                for (std::size_t index = 0; index < calls.size(); ++index) {
                    const auto iter = normal_calls.find(
                        AlleleKey{calls[index]->candidate.tid, calls[index]->candidate.position,
                                  ref, allele_alt(calls[index]->candidate)});
                    if (iter == normal_calls.end()) continue;
                    const auto candidate = candidate_index(*normal_view, iter->second->candidate);
                    if (candidate >= normal_view->somatic_fragment_reference_count.size() ||
                        candidate >= normal_view->somatic_fragment_alternate_count.size())
                        continue;
                    if (index == 0) {
                        normal_fad.front() = std::to_string(
                            normal_view->somatic_fragment_reference_count[candidate]);
                    }
                    normal_fad[index + 1] = std::to_string(
                        normal_view->somatic_fragment_alternate_count[candidate]);
                }
            }
            const auto& normal_annotation_counts = format_view.is_tumor
                ? emitted_output_tumor_annotation_counts[format_view.evidence_index]
                : emitted_normal_annotation_counts[format_view.evidence_index];
            // A selected sample with no reads remains a real output matrix
            // member.  Its count annotations are concrete zeroes, not
            // missing values; this applies equally to a second tumor.
            std::vector<std::string> normal_f1r2(calls.size() + 1, "0");
            std::vector<std::string> normal_f2r1(calls.size() + 1, "0");
            std::string normal_sb = "0,0,0,0";
            if (normal_annotation_counts.has_value()) {
                for (std::size_t index = 0; index < calls.size(); ++index) {
                    const auto candidate = normal_candidate_indices[index];
                    if (candidate >= normal_annotation_counts->fragment_f1r2_alternate.size())
                        continue;
                    normal_f1r2[index + 1] = std::to_string(
                        normal_annotation_counts->fragment_f1r2_alternate[candidate]);
                    normal_f2r1[index + 1] = std::to_string(
                        normal_annotation_counts->fragment_f2r1_alternate[candidate]);
                }
                const auto candidate = normal_candidate_indices.empty()
                    ? std::numeric_limits<std::size_t>::max() : normal_candidate_indices.front();
                if (candidate < normal_annotation_counts->reference_forward.size()) {
                    normal_f1r2.front() = std::to_string(
                        normal_annotation_counts->fragment_f1r2_reference[candidate]);
                    normal_f2r1.front() = std::to_string(
                        normal_annotation_counts->fragment_f2r1_reference[candidate]);
                    std::uint32_t alt_forward = 0;
                    std::uint32_t alt_reverse = 0;
                    for (const auto index : normal_candidate_indices) {
                        if (index < normal_annotation_counts->alternate_forward.size()) {
                            alt_forward += normal_annotation_counts->alternate_forward[index];
                            alt_reverse += normal_annotation_counts->alternate_reverse[index];
                        }
                    }
                    normal_sb = std::to_string(normal_annotation_counts->reference_forward[candidate]) + ',' +
                        std::to_string(normal_annotation_counts->reference_reverse[candidate]) + ',' +
                        std::to_string(alt_forward) + ',' + std::to_string(alt_reverse);
                }
            }
            // SomaticGenotypingEngine gives every non-normal sample the
            // common tumor allele genotype, while matched normals are
            // intentionally hom-ref.  phaseVC() subsequently phases both
            // roles, so their separators must follow the call-level phase.
            const auto normal_gt = format_view.is_tumor
                ? genotype_text.str()
                : std::string(phased ? "0|0" : "0/0");
            std::ostringstream normal_sample_text;
            normal_sample_text << normal_gt << ':' << join_strings(normal_ad) << ':'
                               << join_strings(normal_af) << ':' << normal_depth << ':'
                               << join_strings(normal_f1r2) << ':' << join_strings(normal_f2r1) << ':'
                               << join_strings(normal_fad);
            if (phased) normal_sample_text << ':' << phase_iter->second.pgt << ':'
                                           << phase_iter->second.pid << ':'
                                           << phase_iter->second.ps;
            normal_sample_text << ':' << normal_sb;
            format_sample_texts.push_back(normal_sample_text.str());
        }
        output << '\t' << format;
        for (const auto& sample : output_sample_order) {
            if (!multi_tumor_output && sample.is_tumor) {
                output << '\t' << tumor_sample_text.str();
            } else if (sample.format_index < format_sample_texts.size()) {
                output << '\t' << format_sample_texts[sample.format_index];
            } else {
                output << "\t.";
            }
        }
        output << '\n';
    }
    return output.str();
}

void write_text(const std::string& path, const std::string& text) {
    if (path.empty()) return;
    if (suffix(path, ".gz")) {
#if FASTGATK_HAS_HTSLIB
        BGZF* file = bgzf_open(path.c_str(), "w");
        if (!file) throw std::runtime_error("cannot open compressed output: " + path);
        const auto written = bgzf_write(file, text.data(), text.size());
        const auto close_status = bgzf_close(file);
        if (written != static_cast<ssize_t>(text.size()) || close_status != 0)
            throw std::runtime_error("cannot write compressed output: " + path);
#else
        throw std::runtime_error("compressed output requires HTSlib: " + path);
#endif
    } else {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot open output: " + path);
        output << text;
    }
}

void write_assembly_region_igv(
    const std::string& path,
    const fastgatk::io::HeaderSummary& header,
    const fastgatk::calling::Result& result) {
    if (path.empty()) return;
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write assembly-region output: " + path);
    output << "#track graphType=line\n"
           << "Chromosome\tStart\tEnd\tFeature\tAssemblyRegions\n";
    output << std::fixed << std::setprecision(5);
    for (const auto& region : result.activity_profile_regions) {
        if (region.tid < 0 || static_cast<std::size_t>(region.tid) >= header.contigs.size() ||
            region.start < 0 || region.end < region.start) continue;
        const auto& contig = header.contigs[static_cast<std::size_t>(region.tid)];
        const auto end = static_cast<std::int64_t>(region.end) + 1;
        output << contig << '\t' << region.start << '\t' << region.start + 1
               << "\tend-marker\t0.00000\n"
               << contig << '\t' << region.start << '\t' << end
               << "\tsize=" << (region.end - region.start + 1) << '\t'
               << (region.active ? 1.0 : -1.0) << '\n';
    }
    if (!output) throw std::runtime_error("cannot finalize assembly-region output: " + path);
}

// Mutect2Engine increments its callable counter after the standard Mutect2
// ReadFilters but before it partitions the pileup into tumor and normal
// samples.  Reuse the Kokkos ordinal filter here, then materialize the small
// Host CIGAR payload needed for the stats sweep.  This prevents a sidecar-only
// filter interpretation from drifting from the regular calling path.
fastgatk::io::ReadBatch filtered_mutect_callable_reads(
    const fastgatk::io::ReadBatch& reads,
    const fastgatk::calling::Options& calling) {
    fastgatk::kernels::ReadFilterInput input;
    input.offsets = reads.offsets;
    input.tids = reads.tids;
    input.positions = reads.positions;
    input.mapq = reads.mapq;
    input.flags = reads.flags;
    input.cigar_offsets = reads.cigar_offsets;
    input.cigar_ops = reads.cigar_ops;
    input.read_group_offsets = reads.read_group_offsets;
    input.read_groups = reads.read_groups;
    input.original_alignment_offsets = reads.original_alignment_offsets;
    input.original_alignments = reads.original_alignments;
    input.original_alignment_present = reads.original_alignment_present;
    input.mate_contig_offsets = reads.mate_contig_offsets;
    input.mate_contigs = reads.mate_contigs;
    input.mate_contig_present = reads.mate_contig_present;

    fastgatk::kernels::ReadFilterOptions options;
    options.min_mapq = calling.minimum_mapping_quality;
    options.exclude_mapping_quality_unavailable =
        calling.exclude_mapping_quality_unavailable;
    options.exclude_mapping_quality_zero = calling.exclude_mapping_quality_zero;
    options.exclude_unmapped = calling.exclude_unmapped;
    options.exclude_secondary = calling.exclude_secondary;
    options.exclude_supplementary = calling.exclude_supplementary;
    options.exclude_duplicates = calling.exclude_duplicates;
    options.exclude_qcfail = calling.exclude_qcfail;
    options.require_good_cigar = calling.require_good_cigar;
    options.require_nonzero_reference_span = calling.require_nonzero_reference_span;
    options.require_non_chimeric_original_alignment =
        calling.require_non_chimeric_original_alignment;
    options.require_no_n_cigar = calling.require_no_n_cigar;
    options.require_read_group = calling.require_read_group;
    options.require_read_length = calling.require_read_length;
    options.min_read_length = calling.min_read_length;
    options.max_read_length = calling.max_read_length;
    options.max_reads_per_locus = calling.max_reads_per_locus;
    options.seed = calling.downsampling_seed;
    options.mutect2_downsampling = calling.mutect2_downsampling;
    options.mutect2_downsampling_stride = calling.mutect2_downsampling_stride;
    options.mutect2_max_suspicious_reads_per_alignment_start =
        calling.mutect2_max_suspicious_reads_per_alignment_start;
    const auto filter = fastgatk::kernels::filter_reads_kokkos(input, options);

    fastgatk::io::ReadBatch filtered;
    filtered.offsets.push_back(0);
    for (std::size_t record = 0; record < reads.records(); ++record)
        if (filter.keep[record] != 0U) append_record(filtered, reads, record);
    return filtered;
}

// Mutect2Engine increments callableSites once for every traversal locus whose
// base pileup reaches --callable-depth.  The caller's ActivityProfile owns the
// exact traversal spans, while the flat CIGAR batch remains Host-owned.  Use a
// sweep-line rather than expanding a depth array so whole-contig invocations
// retain bounded O(read-CIGAR + profile-segment) Host state.
std::uint64_t count_mutect_callable_sites(
    const fastgatk::io::ReadBatch& reads,
    const fastgatk::calling::Result& result,
    const std::uint32_t minimum_depth) {
    if (minimum_depth == 0 || !reads.has_cigar() ||
        reads.cigar_offsets.size() != reads.records() + 1 ||
        reads.positions.size() != reads.records() || reads.tids.size() != reads.records())
        return 0;

    using Coordinate = std::int64_t;
    using Span = std::pair<Coordinate, Coordinate>;  // half-open
    std::map<std::int32_t, std::vector<Span>> traversal;
    for (const auto& region : result.activity_profile_regions) {
        if (region.tid < 0 || region.start < 0 || region.end < region.start) continue;
        traversal[region.tid].emplace_back(
            static_cast<Coordinate>(region.start), static_cast<Coordinate>(region.end) + 1);
    }
    for (auto& [_, spans] : traversal) {
        std::sort(spans.begin(), spans.end());
        std::vector<Span> merged;
        for (const auto& span : spans) {
            if (merged.empty() || span.first > merged.back().second)
                merged.push_back(span);
            else
                merged.back().second = std::max(merged.back().second, span.second);
        }
        spans = std::move(merged);
    }

    // delta[tid][position] is the coverage change immediately at position.
    // M/=/X and D participate in Java's base pileup; a reference skip (N)
    // advances the alignment but is not a pileup observation.
    std::map<std::int32_t, std::map<Coordinate, std::int64_t>> delta;
    for (std::size_t record = 0; record < reads.records(); ++record) {
        if (reads.tids[record] < 0 || reads.positions[record] < 0) continue;
        auto reference_position = static_cast<Coordinate>(reads.positions[record]);
        for (std::size_t cigar = reads.cigar_offsets[record];
             cigar < reads.cigar_offsets[record + 1]; ++cigar) {
            const auto operation = fastgatk::io::CigarOp::unpack(reads.cigar_ops[cigar]);
            if (!operation.valid()) break;
            const bool contributes_to_base_pileup =
                operation.code == fastgatk::io::CigarOpCode::Match ||
                operation.code == fastgatk::io::CigarOpCode::SequenceMatch ||
                operation.code == fastgatk::io::CigarOpCode::SequenceMismatch ||
                operation.code == fastgatk::io::CigarOpCode::Deletion;
            if (contributes_to_base_pileup) {
                const auto end = reference_position + static_cast<Coordinate>(operation.length);
                delta[reads.tids[record]][reference_position] += 1;
                delta[reads.tids[record]][end] -= 1;
            }
            if (operation.consumes_reference())
                reference_position += static_cast<Coordinate>(operation.length);
        }
    }

    // A no-profile result only occurs in a deliberately minimal/no-reference
    // caller path.  Its decoded read batch has already been interval-filtered,
    // so sweep the observed coordinate extent rather than silently writing an
    // unusable zero-callable table.
    if (traversal.empty()) {
        for (const auto& [tid, changes] : delta) {
            if (changes.empty()) continue;
            traversal[tid].emplace_back(changes.begin()->first, changes.rbegin()->first);
        }
    }

    // The source increments `callableSites` for every pileup at or above the
    // threshold, then increments once more only when its *raw* isActive()
    // decision accepted a below-threshold locus. ActivityRegion segments are
    // band-pass-expanded and therefore cannot stand in for that second,
    // sparse condition.
    std::map<std::int32_t, std::vector<Coordinate>> raw_active;
    for (const auto& [tid, position] : result.activity_active_loci) {
        if (tid >= 0 && position >= 0)
            raw_active[tid].push_back(static_cast<Coordinate>(position));
    }
    for (auto& [_, positions] : raw_active) {
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    }

    std::uint64_t callable = 0;
    for (const auto& [tid, spans] : traversal) {
        const auto deltas = delta.find(tid);
        if (deltas == delta.end()) continue;
        const auto& changes = deltas->second;
        auto cursor = changes.begin();
        std::int64_t coverage = 0;
        for (const auto& span : spans) {
            while (cursor != changes.end() && cursor->first <= span.first) {
                coverage += cursor->second;
                ++cursor;
            }
            Coordinate position = span.first;
            while (cursor != changes.end() && cursor->first < span.second) {
                if (coverage >= static_cast<std::int64_t>(minimum_depth))
                    callable += static_cast<std::uint64_t>(cursor->first - position);
                position = cursor->first;
                coverage += cursor->second;
                ++cursor;
            }
            if (coverage >= static_cast<std::int64_t>(minimum_depth))
                callable += static_cast<std::uint64_t>(span.second - position);
        }
    }

    for (const auto& [tid, positions] : raw_active) {
        const auto spans = traversal.find(tid);
        const auto deltas = delta.find(tid);
        if (spans == traversal.end() || deltas == delta.end()) continue;
        const auto& changes = deltas->second;
        auto cursor = changes.begin();
        std::int64_t coverage = 0;
        std::size_t span_index = 0;
        for (const auto position : positions) {
            while (cursor != changes.end() && cursor->first <= position) {
                coverage += cursor->second;
                ++cursor;
            }
            while (span_index < spans->second.size() &&
                   spans->second[span_index].second <= position)
                ++span_index;
            if (span_index == spans->second.size() ||
                position < spans->second[span_index].first ||
                coverage >= static_cast<std::int64_t>(minimum_depth))
                continue;
            ++callable;
        }
    }
    return callable;
}

void write_gatk_mutect_stats_table(const std::string& path, const std::uint64_t callable_sites) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write Mutect2 GATK stats: " + path);
    // MutectStats.writeToFile() produces this two-column GATK table.  Emit a
    // decimal value because the Java DataLine writer serializes the statistic
    // as a double (for example, `callable\t373.0`).
    output << "statistic\tvalue\ncallable\t" << std::fixed << std::setprecision(1)
           << static_cast<double>(callable_sites) << "\n";
    output.flush();
    if (!file_complete(path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Mutect2 GATK stats is empty");
}

class IncrementalVcfWriter {
public:
    explicit IncrementalVcfWriter(const std::string& path)
        : path_(path), compressed_(suffix(path, ".gz")) {
        if (path_ == "-" || path_.empty())
            throw std::invalid_argument("Mutect2 region streaming requires file output");
        if (compressed_) {
#if FASTGATK_HAS_HTSLIB
            bgzf_ = bgzf_open(path_.c_str(), "w");
            if (!bgzf_) throw std::runtime_error("cannot open BGZF output: " + path_);
#else
            throw std::runtime_error("compressed output requires HTSlib: " + path_);
#endif
        } else {
            plain_.open(path_);
            if (!plain_) throw std::runtime_error("cannot open output: " + path_);
        }
    }

    IncrementalVcfWriter(const IncrementalVcfWriter&) = delete;
    IncrementalVcfWriter& operator=(const IncrementalVcfWriter&) = delete;

    void write(const std::string& text) {
        if (text.empty()) return;
        if (!compressed_) {
            plain_ << text;
            if (!plain_) throw std::runtime_error("cannot write output: " + path_);
            return;
        }
#if FASTGATK_HAS_HTSLIB
        const auto written = bgzf_write(bgzf_, text.data(), text.size());
        if (written != static_cast<ssize_t>(text.size()))
            throw std::runtime_error("cannot write BGZF output: " + path_);
#endif
    }

    void close() {
        if (closed_) return;
        if (!compressed_) {
            plain_.close();
            if (!plain_) throw std::runtime_error("cannot finalize output: " + path_);
        } else {
#if FASTGATK_HAS_HTSLIB
            if (bgzf_close(bgzf_) != 0)
                throw std::runtime_error("cannot finalize BGZF output: " + path_);
            bgzf_ = nullptr;
#endif
        }
        closed_ = true;
    }

    ~IncrementalVcfWriter() {
        try { close(); } catch (...) {}
    }

private:
    std::string path_;
    bool compressed_ = false;
    bool closed_ = false;
    std::ofstream plain_;
#if FASTGATK_HAS_HTSLIB
    BGZF* bgzf_ = nullptr;
#endif
};

std::pair<std::string, std::string> split_vcf_header_body(const std::string& text) {
    const auto marker = text.find("#CHROM\t");
    if (marker == std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: VCF header is missing #CHROM");
    const auto line_end = text.find('\n', marker);
    if (line_end == std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: VCF header is unterminated");
    return {text.substr(0, line_end + 1), text.substr(line_end + 1)};
}

std::vector<std::string> split_text_field(const std::string& value, const char delimiter) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find(delimiter, begin);
        fields.push_back(value.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1U;
    }
    return fields;
}

std::string join_text_fields(const std::vector<std::string>& fields, const char delimiter) {
    std::ostringstream text;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0U) text << delimiter;
        text << fields[index];
    }
    return text.str();
}

std::string somatic_gvcf_double_text(const double value) {
    std::ostringstream text;
    if (value < 0.01) {
        if (std::abs(value) < 1.0e-20) return "0.00";
        text << std::scientific << std::setprecision(3) << value;
    } else if (value < 1.0) {
        text << std::fixed << std::setprecision(3) << value;
    } else {
        text << std::fixed << std::setprecision(2) << value;
    }
    return text.str();
}

int somatic_gvcf_partition_precision(const std::vector<double>& bands) {
    double smallest_delta = std::numeric_limits<double>::infinity();
    for (std::size_t index = 1; index < bands.size(); ++index)
        smallest_delta = std::min(smallest_delta, bands[index] - bands[index - 1U]);
    // GATK's default contains 0.5-wide bands and therefore uses one decimal
    // of fixed-point precision.  A one-bound custom list has no meaningful
    // pairwise delta; its only partition boundary is representable without a
    // fractional comparison, so use the minimum stable precision of zero.
    if (!std::isfinite(smallest_delta) || smallest_delta >= 1.0) return 0;
    return std::max(0, static_cast<int>(std::ceil(-std::log10(smallest_delta))));
}

std::int64_t somatic_gvcf_lod_bin(const double lod, const int precision) {
    const auto scale = std::pow(10.0, static_cast<double>(precision));
    if (!std::isfinite(lod)) return lod < 0.0
        ? std::numeric_limits<std::int32_t>::min()
        : std::numeric_limits<std::int32_t>::max();
    const auto scaled = std::floor(lod * scale);
    return std::clamp<std::int64_t>(static_cast<std::int64_t>(scaled),
                                    std::numeric_limits<std::int32_t>::min(),
                                    std::numeric_limits<std::int32_t>::max());
}

std::string somatic_reference_confidence_header(
    const Options& options, const std::string& ordinary_header) {
    if (options.reference_confidence_mode == ReferenceConfidenceMode::None)
        return ordinary_header;
    const auto column = ordinary_header.find("#CHROM\t");
    if (column == std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: VCF header is missing #CHROM");
    std::string metadata = ordinary_header.substr(0, column);
    const auto column_header = ordinary_header.substr(column);
    const auto insert_after = [&](const std::string_view anchor,
                                  const std::string_view addition) {
        const auto start = metadata.find(anchor);
        if (start == std::string::npos)
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: VCF header is missing required schema anchor " +
                std::string(anchor));
        const auto end = metadata.find('\n', start);
        if (end == std::string::npos)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: malformed VCF schema header");
        metadata.insert(end + 1U, addition);
    };
    if (metadata.find("##ALT=<ID=NON_REF,") == std::string::npos)
        insert_after("##fileformat=", "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele not already represented at this location by REF and ALT\">\n");
    if (metadata.find("##FORMAT=<ID=GQ,") == std::string::npos)
        insert_after("##FORMAT=<ID=FAD,", "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    if (metadata.find("##FORMAT=<ID=PL,") == std::string::npos)
        insert_after("##FORMAT=<ID=PID,", "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n");
    // SomaticGVCFBlock writes TLOD as a sample FORMAT field for both
    // reference blocks and <NON_REF> records.  The ordinary Mutect2 header
    // already declares the distinct INFO/TLOD annotation, but VCF scopes
    // INFO and FORMAT identifiers independently, so retain GATK's second
    // declaration rather than relying on an undeclared sample column.
    if (metadata.find("##FORMAT=<ID=TLOD,") == std::string::npos)
        insert_after("##FORMAT=<ID=SB,", "##FORMAT=<ID=TLOD,Number=A,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">\n");
    if (options.reference_confidence_mode == ReferenceConfidenceMode::Gvcf) {
        if (metadata.find("##FORMAT=<ID=MIN_DP,") == std::string::npos)
            insert_after("##FORMAT=<ID=GT,", "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">\n");
        if (metadata.find("##INFO=<ID=END,") == std::string::npos)
            insert_after("##INFO=<ID=ECNTH,", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n");
        const auto precision = somatic_gvcf_partition_precision(options.gvcf_lod_bands);
        std::vector<std::string> block_headers;
        std::int64_t lower = std::numeric_limits<std::int32_t>::min();
        for (const auto boundary : options.gvcf_lod_bands) {
            const auto upper = somatic_gvcf_lod_bin(boundary, precision);
            block_headers.push_back("##GVCFBlock" + std::to_string(lower) + "-" +
                std::to_string(upper) + "=minGQ=" + std::to_string(lower) +
                "(inclusive),maxGQ=" + std::to_string(upper) + "(exclusive)\n");
            lower = upper;
        }
        block_headers.push_back("##GVCFBlock" + std::to_string(lower) + "-" +
            std::to_string(std::numeric_limits<std::int32_t>::max()) + "=minGQ=" +
            std::to_string(lower) + "(inclusive),maxGQ=" +
            std::to_string(std::numeric_limits<std::int32_t>::max()) + "(exclusive)\n");
        std::sort(block_headers.begin(), block_headers.end());
        std::ostringstream block_text;
        for (const auto& block_header : block_headers) block_text << block_header;
        insert_after("##FORMAT=<ID=TLOD,", block_text.str());
    }
    return metadata + column_header;
}

struct SomaticGvcfReferenceMetrics {
    std::uint32_t depth = 0;
    std::uint32_t reference_depth = 0;
    std::uint32_t non_reference_depth = 0;
    double lod = 0.0;
};

std::string somatic_reference_confidence_vcf_text(
    const Options& options, const fastgatk::io::HeaderSummary& header,
    const fastgatk::io::ReadBatch& tumor_reads, const fastgatk::calling::Result& tumor,
    const fastgatk::calling::Result* normal,
    const fastgatk::kernels::SomaticLikelihoodResult* tumor_somatic,
    const fastgatk::kernels::SomaticLikelihoodResult* normal_somatic,
    const fastgatk::kernels::SomaticPosteriorResult* posterior,
    const std::vector<NormalSampleView>* normal_views,
    const std::vector<NormalSampleView>* tumor_views,
    const std::vector<std::string>* references,
    const MutectFeatureResources* feature_resources = nullptr) {
    const auto ordinary = vcf_text(options, header, tumor_reads, tumor, normal,
                                   tumor_somatic, normal_somatic, posterior,
                                   normal_views, tumor_views, references,
                                   feature_resources);
    if (options.reference_confidence_mode == ReferenceConfidenceMode::None)
        return ordinary;
    const auto ordinary_parts = split_vcf_header_body(ordinary);
    const auto gvcf_header = somatic_reference_confidence_header(options, ordinary_parts.first);

    std::map<std::string, std::int32_t> contig_ids;
    for (std::size_t index = 0; index < header.contigs.size(); ++index)
        contig_ids.emplace(header.contigs[index], static_cast<std::int32_t>(index));
    std::map<std::pair<std::int32_t, std::int32_t>, SomaticGvcfReferenceMetrics> metrics;
    for (const auto& locus : tumor.somatic_reference_confidence_loci) {
        if (locus.tid < 0 || locus.position < 0 || !std::isfinite(locus.log10_lod)) continue;
        metrics[std::make_pair(locus.tid, locus.position)] = SomaticGvcfReferenceMetrics{
            locus.depth, locus.reference_depth, locus.non_reference_depth, locus.log10_lod};
    }
    // A concrete gVCF record's symbolic <NON_REF> TLOD is not its pileup
    // reference-confidence score.  SomaticGenotypingEngine appends the
    // synthetic row to the local fragment likelihood matrix and evaluates it
    // with the same all-alleles-vs-without-allele model as each EventMap ALT.
    // Keep the regional active-core owner's Kokkos result when flattened
    // candidate metadata contains an overlapping halo duplicate.
    struct NonReferenceLikelihood {
        int priority = 0;
        double tlod = 0.0;
        double allele_fraction = 0.0;
    };
    std::map<std::pair<std::int32_t, std::int32_t>, NonReferenceLikelihood>
        variant_non_reference_likelihoods;
    if (tumor_somatic != nullptr) {
        const bool explicit_ownership =
            tumor.likelihood_candidate_event_map_owned.size() == tumor.candidates.size();
        for (std::size_t candidate = 0;
             candidate < tumor.candidates.size() &&
             candidate < tumor_somatic->non_reference_tlod.size(); ++candidate) {
            const auto lod = tumor_somatic->non_reference_tlod[candidate];
            if (!std::isfinite(lod)) continue;
            const auto& allele = tumor.candidates[candidate];
            int priority = allele.graph_derived ? 1 : 0;
            if (explicit_ownership && tumor.likelihood_candidate_event_map_owned[candidate] != 0U)
                priority = 2;
            const auto key = std::make_pair(allele.tid, allele.position);
            const auto existing = variant_non_reference_likelihoods.find(key);
            if (existing == variant_non_reference_likelihoods.end() ||
                priority > existing->second.priority) {
                variant_non_reference_likelihoods[key] = NonReferenceLikelihood{
                    priority, lod,
                    candidate < tumor_somatic->non_reference_allele_fraction.size()
                        ? tumor_somatic->non_reference_allele_fraction[candidate] : 0.0};
            }
        }
    }

    struct VariantRow {
        std::string text;
        std::int32_t end = -1;
    };
    std::map<std::pair<std::int32_t, std::int32_t>, VariantRow> variants;
    for (const auto& line : split_text_field(ordinary_parts.second, '\n')) {
        if (line.empty()) continue;
        auto fields = split_text_field(line, '\t');
        if (fields.size() < 8U) continue;
        const auto contig = contig_ids.find(fields[0]);
        if (contig == contig_ids.end()) continue;
        const auto position_one_based = std::stoll(fields[1]);
        if (position_one_based < 1 || position_one_based > std::numeric_limits<std::int32_t>::max())
            continue;
        const auto position = static_cast<std::int32_t>(position_one_based - 1);
        const auto metric = metrics.find(std::make_pair(contig->second, position));
        const auto variant_likelihood = variant_non_reference_likelihoods.find(
            std::make_pair(contig->second, position));
        const auto non_ref_lod = variant_likelihood != variant_non_reference_likelihoods.end()
            ? variant_likelihood->second.tlod
            : (metric == metrics.end() ? 0.0 : metric->second.lod);
        const auto non_ref_allele_fraction =
            variant_likelihood != variant_non_reference_likelihoods.end()
            ? variant_likelihood->second.allele_fraction : 0.0;
        const auto existing_alt_count = fields[4] == "." ? std::size_t{0} :
            split_text_field(fields[4], ',').size();
        fields[4] = fields[4] == "." ? "<NON_REF>" : fields[4] + ",<NON_REF>";
        auto info = split_text_field(fields[7], ';');
        for (auto& item : info) {
            const auto equals = item.find('=');
            const auto key = equals == std::string::npos ? item : item.substr(0, equals);
            if (equals == std::string::npos) continue;
            auto value = item.substr(equals + 1U);
            if (key == "AS_SB_TABLE") value += "|0,0";
            else if (key == "MBQ" || key == "MFRL") value += ",0";
            // PerAlleleAnnotation assigns benign defaults to an allele with
            // no BestAllele evidence: MappingQuality=60 and ReadPosition=50
            // (whereas base/fragment quality naturally remain zero).  The
            // synthetic <NON_REF> row has no directly supporting read, so
            // preserve those source defaults at the gVCF writer boundary.
            else if (key == "MMQ") value += ",60";
            else if (key == "MPOS") value += ",50";
            else if (key == "POPAF") value += "," + somatic_gvcf_double_text(
                -std::log10(options.population_allele_frequency));
            else if (key == "TLOD") value += "," + somatic_gvcf_double_text(non_ref_lod);
            item = key + "=" + value;
        }
        fields[7] = join_text_fields(info, ';');
        if (fields.size() > 9U && fields[8] != ".") {
            const auto keys = split_text_field(fields[8], ':');
            for (std::size_t sample = 9; sample < fields.size(); ++sample) {
                auto values = split_text_field(fields[sample], ':');
                if (values.size() < keys.size()) values.resize(keys.size(), ".");
                for (std::size_t key = 0; key < keys.size(); ++key) {
                    if (keys[key] == "GT" && values[key] != "." && values[key] != "./." &&
                        values[key] != "0/0" && values[key] != "0|0") {
                        values[key] += "/" + std::to_string(existing_alt_count + 1U);
                    } else if (keys[key] == "AF") {
                        values[key] += "," + somatic_gvcf_double_text(non_ref_allele_fraction);
                    } else if (keys[key] == "AD" ||
                               keys[key] == "F1R2" || keys[key] == "F2R1" || keys[key] == "FAD") {
                        values[key] += ",0";
                    }
                }
                fields[sample] = join_text_fields(values, ':');
            }
        }
        const auto reference_length = std::max<std::size_t>(1U, fields[3].size());
        const auto end = reference_length - 1U >
                static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max() - position)
            ? std::numeric_limits<std::int32_t>::max()
            : position + static_cast<std::int32_t>(reference_length - 1U);
        variants[std::make_pair(contig->second, position)] = VariantRow{
            join_text_fields(fields, '\t'), end};
    }

    std::vector<fastgatk::io::HtsInterval> domains = options.calling.locus_intervals.empty()
        ? options.calling.intervals : options.calling.locus_intervals;
    if (domains.empty() && references != nullptr) {
        for (std::size_t tid = 0; tid < references->size(); ++tid) {
            if (!(*references)[tid].empty() &&
                tid <= static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
                domains.push_back(fastgatk::io::HtsInterval{
                    static_cast<std::int32_t>(tid), 0,
                    static_cast<std::int64_t>((*references)[tid].size())});
        }
    }
    std::sort(domains.begin(), domains.end(), [](const auto& left, const auto& right) {
        return std::tie(left.tid, left.start, left.end) < std::tie(right.tid, right.start, right.end);
    });
    std::vector<fastgatk::io::HtsInterval> merged_domains;
    for (const auto& domain : domains) {
        if (domain.tid < 0 || domain.start < 0 || domain.end <= domain.start) continue;
        if (!merged_domains.empty() && merged_domains.back().tid == domain.tid &&
            domain.start <= merged_domains.back().end) {
            merged_domains.back().end = std::max(merged_domains.back().end, domain.end);
        } else {
            merged_domains.push_back(domain);
        }
    }

    struct OutputRow { std::int32_t tid; std::int32_t position; int order; std::string text; };
    std::vector<OutputRow> rows;
    const auto precision = somatic_gvcf_partition_precision(options.gvcf_lod_bands);
    const auto lod_band = [&](const double lod) {
        const auto value = somatic_gvcf_lod_bin(lod, precision);
        for (std::size_t index = 0; index < options.gvcf_lod_bands.size(); ++index)
            if (value < somatic_gvcf_lod_bin(options.gvcf_lod_bands[index], precision))
                return static_cast<int>(index);
        return static_cast<int>(options.gvcf_lod_bands.size());
    };
    const auto median_depth = [](std::vector<std::uint32_t> values) {
        if (values.empty()) return std::uint32_t{0};
        std::sort(values.begin(), values.end());
        const auto middle = values.size() / 2U;
        if ((values.size() & 1U) != 0U) return values[middle];
        return static_cast<std::uint32_t>((static_cast<std::uint64_t>(values[middle - 1U]) +
                                            values[middle] + 1U) / 2U);
    };
    for (const auto& domain : merged_domains) {
        if (references == nullptr || domain.tid < 0 ||
            static_cast<std::size_t>(domain.tid) >= references->size() ||
            static_cast<std::size_t>(domain.tid) >= header.contigs.size())
            continue;
        const auto& sequence = (*references)[static_cast<std::size_t>(domain.tid)];
        const auto start = std::max<std::int64_t>(0, domain.start);
        const auto end = std::min<std::int64_t>(domain.end, sequence.size());
        struct Block {
            bool valid = false; std::int32_t start = -1; std::int32_t end = -1; int band = -1;
            double min_lod = std::numeric_limits<double>::infinity();
            std::uint32_t min_depth = std::numeric_limits<std::uint32_t>::max();
            std::vector<std::uint32_t> depths;
        } block;
        const auto flush_block = [&] {
            if (!block.valid) return;
            const auto reference = sequence[static_cast<std::size_t>(block.start)];
            std::ostringstream text;
            text << header.contigs[static_cast<std::size_t>(domain.tid)] << '\t' << block.start + 1
                 << "\t.\t" << reference << "\t<NON_REF>\t.\t.\tEND=" << block.end + 1;
            if (!options.sites_only_vcf_output)
                text << "\tGT:DP:MIN_DP:TLOD\t0/0:" << median_depth(block.depths) << ':'
                     << (block.min_depth == std::numeric_limits<std::uint32_t>::max()
                         ? 0U : block.min_depth) << ':'
                     << somatic_gvcf_double_text(block.min_lod);
            rows.push_back(OutputRow{domain.tid, block.start, 0, text.str()});
            block = Block{};
        };
        for (std::int64_t position = start; position < end;) {
            const auto key = std::make_pair(domain.tid, static_cast<std::int32_t>(position));
            const auto variant = variants.find(key);
            if (variant != variants.end()) {
                flush_block();
                rows.push_back(OutputRow{domain.tid, static_cast<std::int32_t>(position), 1,
                                         variant->second.text});
                position = std::max<std::int64_t>(position + 1,
                    static_cast<std::int64_t>(variant->second.end) + 1);
                continue;
            }
            const auto found = metrics.find(key);
            const auto metric = found == metrics.end() ? SomaticGvcfReferenceMetrics{} : found->second;
            if (options.reference_confidence_mode == ReferenceConfidenceMode::BpResolution) {
                std::ostringstream text;
                text << header.contigs[static_cast<std::size_t>(domain.tid)] << '\t' << position + 1
                    << "\t.\t" << sequence[static_cast<std::size_t>(position)]
                    << "\t<NON_REF>\t.\t.\t.";
                if (!options.sites_only_vcf_output)
                    text << "\tGT:AD:DP:TLOD\t0/0:" << metric.reference_depth << ','
                         << metric.non_reference_depth << ':' << metric.depth << ':'
                         << somatic_gvcf_double_text(metric.lod);
                rows.push_back(OutputRow{domain.tid, static_cast<std::int32_t>(position), 0,
                                         text.str()});
            } else {
                const auto band = lod_band(metric.lod);
                if (!block.valid || block.band != band || block.end + 1 != position) {
                    flush_block();
                    block.valid = true;
                    block.start = static_cast<std::int32_t>(position);
                    block.end = static_cast<std::int32_t>(position);
                    block.band = band;
                } else {
                    block.end = static_cast<std::int32_t>(position);
                }
                block.min_lod = std::min(block.min_lod, metric.lod);
                block.min_depth = std::min(block.min_depth, metric.depth);
                block.depths.push_back(metric.depth);
            }
            ++position;
        }
        flush_block();
    }
    std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
        return std::tie(left.tid, left.position, left.order) <
               std::tie(right.tid, right.position, right.order);
    });
    std::ostringstream output;
    output << gvcf_header;
    for (const auto& row : rows) output << row.text << '\n';
    return output.str();
}

struct Mutect2RegionDecoded {
    fastgatk::io::ReadBatch tumor;
    fastgatk::io::ReadBatch normal;
    std::vector<std::string> tumor_sample_names;
    std::vector<fastgatk::io::ReadBatch> tumor_sample_reads;
    std::vector<std::string> normal_sample_names;
    std::vector<fastgatk::io::ReadBatch> normal_sample_reads;
    bool has_normal = false;
    std::int32_t tid = -1;
    std::int64_t core_start = 0;
    std::int64_t core_end = 0;
    // This is the tile's read/activity visibility, not its output core.
    // It is intersected with the original -L selection before decode so a
    // streamed invocation cannot see reads that the corresponding aggregate
    // invocation was never permitted to consume.
    std::vector<fastgatk::io::HtsInterval> locus_intervals;
};

struct Mutect2RegionComputed {
    std::string body;
    std::string f1r2_text;
    fastgatk::somatic::F1R2ArchiveSample f1r2_archive;
    std::uint64_t tumor_reads = 0;
    std::uint64_t normal_reads = 0;
    std::uint64_t candidate_sites = 0;
    std::uint64_t tumor_calls = 0;
    std::uint64_t normal_calls = 0;
    std::uint64_t pairhmm_pairs = 0;
    std::uint64_t itr_artifact_clipped_reads = 0;
    std::uint64_t downsampled_reads = 0;
    std::uint64_t overlapping_pairs = 0;
    std::uint64_t overlapping_bases = 0;
    std::uint64_t overlapping_conflicting_bases = 0;
    std::uint64_t overlapping_quality_caps = 0;
    std::uint64_t overlapping_indel_quality_caps = 0;
    std::uint64_t f1r2_alt_support = 0;
    std::uint64_t r1f2_alt_support = 0;
    std::uint64_t somatic_evidence_groups = 0;
    double somatic_likelihood_prepare_seconds = 0.0;
    double somatic_likelihood_seconds = 0.0;
    double somatic_posterior_prepare_seconds = 0.0;
    double somatic_posterior_seconds = 0.0;
    std::string somatic_likelihood_execution_space;
    std::string somatic_posterior_execution_space;
};

struct Mutect2StreamSummary {
    std::uint64_t regions = 0;
    std::uint64_t tumor_reads = 0;
    std::uint64_t normal_reads = 0;
    std::uint64_t candidate_sites = 0;
    std::uint64_t tumor_calls = 0;
    std::uint64_t normal_calls = 0;
    std::uint64_t pairhmm_pairs = 0;
    std::uint64_t itr_artifact_clipped_reads = 0;
    std::uint64_t downsampled_reads = 0;
    std::uint64_t overlapping_pairs = 0;
    std::uint64_t overlapping_bases = 0;
    std::uint64_t overlapping_conflicting_bases = 0;
    std::uint64_t overlapping_quality_caps = 0;
    std::uint64_t overlapping_indel_quality_caps = 0;
    std::uint64_t f1r2_alt_support = 0;
    std::uint64_t r1f2_alt_support = 0;
    std::uint64_t somatic_evidence_groups = 0;
    double somatic_likelihood_prepare_seconds = 0.0;
    double somatic_likelihood_seconds = 0.0;
    double somatic_posterior_prepare_seconds = 0.0;
    double somatic_posterior_seconds = 0.0;
    std::string somatic_likelihood_execution_space;
    std::string somatic_posterior_execution_space;
    fastgatk::somatic::F1R2ArchiveSample f1r2_archive;
};

std::uint64_t mutect2_region_decoded_bytes(const Mutect2RegionDecoded& item) {
    return std::max<std::uint64_t>(1, sizeof(Mutect2RegionDecoded) + item.tumor.bytes() +
                                      item.normal.bytes());
}

std::uint64_t mutect2_region_computed_bytes(const Mutect2RegionComputed& item) {
    std::uint64_t bytes = sizeof(Mutect2RegionComputed) + item.body.size() +
                          item.f1r2_text.size();
    for (const auto& histogram : item.f1r2_archive.ref_hist) bytes += histogram.size() * sizeof(std::uint64_t);
    for (const auto& histogram : item.f1r2_archive.alt_hist) bytes += histogram.size() * sizeof(std::uint64_t);
    for (const auto& row : item.f1r2_archive.alt_rows)
        bytes += sizeof(row) + row.context.size();
    return std::max<std::uint64_t>(1, bytes);
}

void merge_f1r2_archive(fastgatk::somatic::F1R2ArchiveSample& destination,
                        fastgatk::somatic::F1R2ArchiveSample source) {
    if (destination.sample.empty()) destination.sample = source.sample;
    destination.max_depth = std::max(destination.max_depth, source.max_depth);
    for (std::size_t index = 0; index < destination.ref_hist.size(); ++index) {
        if (destination.ref_hist[index].size() < source.ref_hist[index].size())
            destination.ref_hist[index].resize(source.ref_hist[index].size(), 0);
        for (std::size_t depth = 0; depth < source.ref_hist[index].size(); ++depth)
            destination.ref_hist[index][depth] += source.ref_hist[index][depth];
    }
    for (std::size_t index = 0; index < destination.alt_hist.size(); ++index) {
        if (destination.alt_hist[index].size() < source.alt_hist[index].size())
            destination.alt_hist[index].resize(source.alt_hist[index].size(), 0);
        for (std::size_t depth = 0; depth < source.alt_hist[index].size(); ++depth)
            destination.alt_hist[index][depth] += source.alt_hist[index][depth];
    }
    destination.alt_rows.insert(destination.alt_rows.end(),
                                std::make_move_iterator(source.alt_rows.begin()),
                                std::make_move_iterator(source.alt_rows.end()));
}

std::string mutect2_region_f1r2_text(const fastgatk::io::ReadBatch& reads,
                                     const fastgatk::calling::Result& result,
                                     std::uint64_t& f1r2_total,
                                     std::uint64_t& r1f2_total,
                                     const fastgatk::io::HtsInterval* output_core = nullptr) {
    std::ostringstream text;
    for (const auto& call : result.calls) {
        if (output_core != nullptr &&
            (call.candidate.tid != output_core->tid ||
             call.candidate.position < output_core->start ||
             call.candidate.position >= output_core->end))
            continue;
        const auto orientation = orientation_counts(reads, call.candidate);
        f1r2_total += orientation.f1r2;
        r1f2_total += orientation.r1f2;
        text << call.candidate.tid << '\t' << call.candidate.position + 1 << '\t'
             << allele_alt(call.candidate) << '\t' << orientation.f1r2 << '\t'
             << orientation.r1f2 << '\t' << (orientation.f1r2 + orientation.r1f2) << '\n';
    }
    return text.str();
}

std::string vcf_body_for_core(const std::string& body, const std::string& contig,
                              const std::int64_t core_start,
                              const std::int64_t core_end) {
    std::istringstream input(body);
    std::ostringstream output;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto first_tab = line.find('\t');
        if (first_tab == std::string::npos || line.substr(0, first_tab) != contig) continue;
        const auto second_tab = line.find('\t', first_tab + 1U);
        if (second_tab == std::string::npos) continue;
        std::int64_t position_one_based = 0;
        const auto position_begin = line.data() + first_tab + 1U;
        const auto position_end = line.data() + second_tab;
        const auto parsed = std::from_chars(position_begin, position_end, position_one_based);
        if (parsed.ec != std::errc{} || parsed.ptr != position_end) continue;
        const auto position_zero_based = position_one_based - 1;
        if (position_zero_based < core_start || position_zero_based >= core_end) continue;
        output << line << '\n';
    }
    return output.str();
}

Mutect2StreamSummary run_region_streaming(
    Options& options,
    const fastgatk::io::HtsReader& metadata_reader,
    const fastgatk::io::HeaderSummary& merged_header,
    const std::vector<std::string>& tumor_inputs,
    const std::vector<std::string>& references,
    const MutectFeatureResources& feature_resources,
    const fastgatk::runtime::ResourceSnapshot& resources) {
    if (options.output == "-" || options.output.empty() ||
        !(suffix(options.output, ".vcf") || suffix(options.output, ".vcf.gz")))
        throw std::invalid_argument("--stream-by-region requires VCF/VCF.GZ file output");
    if (!metadata_reader.has_index())
        throw std::runtime_error(
            "BACKEND_UNAVAILABLE: --stream-by-region requires an indexed BAM/CRAM input");
    configure_mutect2_itr_artifact_transform(options, references);

    struct CoreInterval {
        std::int32_t tid = -1;
        std::int64_t start = 0;
        std::int64_t end = 0;
    };
    std::vector<CoreInterval> cores;
    const auto append_cores = [&](const std::int32_t tid,
                                  const std::int64_t start,
                                  const std::int64_t end) {
        if (tid < 0 || start < 0 || end <= start ||
            static_cast<std::size_t>(tid) >= merged_header.contig_lengths.size()) return;
        const auto length = merged_header.contig_lengths[static_cast<std::size_t>(tid)];
        const auto bounded_end = std::min(end, length);
        const auto bounded_start = std::min(start, bounded_end);
        for (auto cursor = bounded_start; cursor < bounded_end;) {
            const auto next = std::min<std::int64_t>(
                bounded_end, cursor + static_cast<std::int64_t>(options.stream_region_size));
            cores.push_back(CoreInterval{tid, cursor, next});
            cursor = next;
        }
    };
    if (!metadata_reader.intervals().empty()) {
        for (const auto& interval : metadata_reader.intervals())
            append_cores(interval.tid, interval.start, interval.end);
    } else {
        for (std::size_t tid = 0; tid < merged_header.contig_lengths.size(); ++tid)
            append_cores(static_cast<std::int32_t>(tid), 0,
                         merged_header.contig_lengths[tid]);
    }
    if (cores.empty()) throw std::runtime_error("BAD_INPUT: no reference intervals to stream");

    const auto tumor_sample_names = selected_tumor_samples(options, merged_header);
    if (tumor_sample_names.empty())
        throw std::runtime_error(
            "BAD_INPUT: every tumor input sample was selected as a normal sample");
    if (options.reference_confidence_mode != ReferenceConfidenceMode::None &&
        tumor_sample_names.size() != 1)
        throw std::runtime_error(
            "BAD_INPUT: Mutect2 reference confidence requires exactly one sample; "
            "multi-tumor inputs are unsupported by GATK 4.6.2");

    IncrementalVcfWriter writer(options.output);
    fastgatk::io::ReadBatch empty_reads;
    fastgatk::calling::Result empty_result;
    fastgatk::calling::Result empty_normal;
    std::vector<NormalSampleView> header_normal_views;
    if (!options.normal.empty() || !options.normal_sample.empty()) {
        const auto names = selected_normal_samples(options);
        if (names.empty()) {
            header_normal_views.push_back(NormalSampleView{"NORMAL", nullptr, nullptr, nullptr});
        } else {
            for (const auto& name : names)
                header_normal_views.push_back(NormalSampleView{name, nullptr, nullptr, nullptr});
        }
    }
    std::vector<NormalSampleView> header_tumor_views;
    if (tumor_sample_names.size() > 1) {
        header_tumor_views.reserve(tumor_sample_names.size());
        for (const auto& name : tumor_sample_names)
            header_tumor_views.push_back(NormalSampleView{name, nullptr, nullptr, nullptr});
    }
    const auto header_parts = split_vcf_header_body(
        somatic_reference_confidence_vcf_text(
            options, merged_header, empty_reads, empty_result,
            (!options.normal.empty() || !options.normal_sample.empty())
                ? &empty_normal : nullptr,
            nullptr, nullptr, nullptr,
            header_normal_views.empty() ? nullptr : &header_normal_views,
            header_tumor_views.empty() ? nullptr : &header_tumor_views,
            &references));
    writer.write(header_parts.first);

    const auto stats_path = options.stats.empty() ? options.output + ".stats.json" : options.stats;
    const auto f1r2_path = options.f1r2.empty() ? options.output + ".f1r2.tsv" : options.f1r2;
    const bool standard_f1r2_archive = suffix(f1r2_path, ".tar.gz");
    std::ofstream f1r2_tsv;
    if (!standard_f1r2_archive) {
        f1r2_tsv.open(f1r2_path);
        if (!f1r2_tsv) throw std::runtime_error("cannot write Mutect2 F1R2 sidecar: " + f1r2_path);
        f1r2_tsv << "# FASTGATK-MUTECT2-F1R2 v2\n# contig\tposition\talt\tf1r2\tr1f2\ttotal\n";
    }

    const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    const auto stage_capacity = safe_memory_budget > 0
        ? std::max<std::uint64_t>(safe_memory_budget, 1)
        : 64ULL * 1024ULL * 1024ULL;
    const auto context = static_cast<std::int64_t>(std::max({
        static_cast<std::uint64_t>(options.calling.assembly_region_padding),
        static_cast<std::uint64_t>(options.calling.max_probability_propagation_distance) + 50U,
        static_cast<std::uint64_t>(options.calling.indel_padding_for_genotyping)}));
    // GATK decides ActivityRegion ownership across the complete selected
    // traversal, then invokes callRegion() on those regions.  An arbitrary
    // streaming tile is not such a boundary.  For an explicit -L selector,
    // replay that same selector-wide assembly input for every output core;
    // each decoded batch is still released before the next tile, while the
    // core filter below prevents duplicate or out-of-interval VCF records.
    // A no--L traversal retains the ordinary local halo so it can stream an
    // entire contig without staging it as one input window.
    const auto selected_input_halo = expanded_assembly_fetch_intervals(
        merged_header, metadata_reader.intervals(), options.calling.assembly_region_padding);
    const auto input_halo_intervals = [&](const CoreInterval& core) {
        if (!selected_input_halo.empty()) return selected_input_halo;
        const auto contig_length = merged_header.contig_lengths[
            static_cast<std::size_t>(core.tid)];
        const auto halo_start = std::max<std::int64_t>(0, core.start - context);
        const auto halo_end = std::min<std::int64_t>(contig_length, core.end + context);
        // AssemblyRegionWalker expands an indexed -L query by its assembly
        // context before building the graph.  The VCF core is constrained
        // separately below by tile_options.intervals; clipping this fetch
        // back to metadata_reader.intervals() drops read starts just outside
        // -L and changes the EventMap at the boundary.
        return std::vector<fastgatk::io::HtsInterval>{
            fastgatk::io::HtsInterval{core.tid, halo_start, halo_end}};
    };
    std::deque<CoreInterval> pending(cores.begin(), cores.end());
    std::atomic<std::uint64_t> peak_host_bytes{0};
    const auto update_peak = [&](const std::uint64_t bytes) {
        auto observed = peak_host_bytes.load(std::memory_order_relaxed);
        while (observed < bytes &&
               !peak_host_bytes.compare_exchange_weak(
                   observed, bytes, std::memory_order_relaxed, std::memory_order_relaxed)) {}
    };
    std::size_t decode_batch_records = options.batch_records;
    std::uint64_t adaptive_reductions = 0;
    Mutect2StreamSummary summary;

    using Pipeline = fastgatk::runtime::ThreeStagePipeline<
        Mutect2RegionDecoded, Mutect2RegionComputed, Mutect2RegionComputed>;
    auto decode = [&]() -> std::optional<Mutect2RegionDecoded> {
        while (!pending.empty()) {
            const auto core = pending.front();
            pending.pop_front();
            auto halo_intervals = input_halo_intervals(core);
            if (halo_intervals.empty()) continue;
            const auto& contig = merged_header.contigs[static_cast<std::size_t>(core.tid)];
            std::vector<std::string> tile_regions;
            tile_regions.reserve(halo_intervals.size());
            for (const auto& halo : halo_intervals) {
                tile_regions.push_back(contig + ":" + std::to_string(halo.start + 1) +
                                       "-" + std::to_string(halo.end));
            }
            std::size_t effective_tumor_batch = decode_batch_records;
            std::uint64_t tumor_reductions = 0;
            fastgatk::io::HeaderSummary tile_header;
            fastgatk::io::ReadBatch tumor;
            std::vector<fastgatk::io::ReadBatch> tumor_sample_reads;
            if (tumor_sample_names.size() == 1) {
                tumor = read_all_many(tumor_inputs, options.reference, tile_regions,
                                      decode_batch_records, tile_header, 0,
                                      options.tumor_sample, false, nullptr, nullptr, nullptr,
                                      &tumor_reductions, &effective_tumor_batch, nullptr, true,
                                      options.interval_set_rule);
            } else {
                tumor.offsets.push_back(0);
                tumor_sample_reads.reserve(tumor_sample_names.size());
                for (std::size_t sample_index = 0;
                     sample_index < tumor_sample_names.size(); ++sample_index) {
                    const auto& sample = tumor_sample_names[sample_index];
                    fastgatk::io::HeaderSummary sample_header;
                    std::uint64_t sample_reductions = 0;
                    std::size_t sample_effective_batch = decode_batch_records;
                    auto sample_reads = read_all_many(
                        tumor_inputs, options.reference, tile_regions,
                        decode_batch_records, sample_header, 0, sample, false,
                        nullptr, nullptr, nullptr, &sample_reductions,
                        &sample_effective_batch, nullptr, true,
                        options.interval_set_rule);
                    sort_batch_by_coordinate(sample_reads);
                    label_batch_sample(sample_reads,
                                       static_cast<std::int32_t>(sample_index));
                    append_batch(tumor, sample_reads);
                    tumor_sample_reads.push_back(std::move(sample_reads));
                    tumor_reductions += sample_reductions;
                    effective_tumor_batch = std::min(effective_tumor_batch,
                                                      sample_effective_batch);
                    tile_header = std::move(sample_header);
                }
            }
            sort_batch_by_coordinate(tumor);
            fastgatk::io::ReadBatch normal;
            std::vector<std::string> normal_sample_names;
            std::vector<fastgatk::io::ReadBatch> normal_sample_reads;
            const auto tumor_role_count = std::max<std::size_t>(1, tumor_sample_names.size());
            std::size_t effective_normal_batch = decode_batch_records;
            std::uint64_t normal_reductions = 0;
            // GATK accepts one multi-sample input with --normal-sample; the
            // normal path is then the same HTSlib stream, filtered by its
            // @RG/SM mapping.  --normal-input remains a compatibility alias
            // for callers that stage tumor and normal files separately.
            if (!options.normal.empty() || !options.normal_sample.empty()) {
                fastgatk::io::HeaderSummary normal_header;
                const auto normal_inputs = options.normal.empty()
                    ? tumor_inputs : std::vector<std::string>{options.normal};
                const auto requested_normal_samples = selected_normal_samples(options);
                if (requested_normal_samples.empty()) {
                    normal = read_all_many(normal_inputs, options.reference, tile_regions,
                                           decode_batch_records, normal_header, 0,
                                           {}, false, nullptr, nullptr, nullptr,
                                           &normal_reductions, &effective_normal_batch, nullptr, true,
                                           options.interval_set_rule);
                    sort_batch_by_coordinate(normal);
                    label_batch_sample(normal,
                                       static_cast<std::int32_t>(tumor_role_count));
                } else {
                    normal.offsets.push_back(0);
                    for (std::size_t normal_index = 0;
                         normal_index < requested_normal_samples.size(); ++normal_index) {
                        const auto& sample = requested_normal_samples[normal_index];
                        fastgatk::io::HeaderSummary sample_header;
                        std::uint64_t sample_reductions = 0;
                        std::size_t sample_effective_batch = decode_batch_records;
                        auto sample_reads = read_all_many(
                            normal_inputs, options.reference, tile_regions,
                            decode_batch_records, sample_header, 0, sample, false,
                            nullptr, nullptr, nullptr, &sample_reductions,
                            &sample_effective_batch, nullptr, true,
                            options.interval_set_rule);
                        sort_batch_by_coordinate(sample_reads);
                        label_batch_sample(sample_reads, static_cast<std::int32_t>(
                            tumor_role_count + normal_index));
                        append_batch(normal, sample_reads);
                        normal_sample_names.push_back(sample);
                        normal_sample_reads.push_back(std::move(sample_reads));
                        normal_reductions += sample_reductions;
                        effective_normal_batch = std::min(effective_normal_batch,
                                                          sample_effective_batch);
                    }
                }
            }
            adaptive_reductions += tumor_reductions + normal_reductions;
            decode_batch_records = std::max<std::size_t>(1, std::min(
                effective_tumor_batch, effective_normal_batch));
            const auto bytes = std::max<std::uint64_t>(
                1, sizeof(Mutect2RegionDecoded) + tumor.bytes() + normal.bytes());
            update_peak(bytes);
            if ((safe_memory_budget > 0 && bytes > safe_memory_budget) || bytes > stage_capacity) {
                const auto span = core.end - core.start;
                if (span <= 1)
                    throw std::runtime_error(
                        "RESOURCE_EXHAUSTED: single-base Mutect2 tile exceeds safe memory budget");
                const auto midpoint = core.start + span / 2;
                pending.push_front(CoreInterval{core.tid, midpoint, core.end});
                pending.push_front(CoreInterval{core.tid, core.start, midpoint});
                options.stream_region_size = std::max<std::size_t>(1, options.stream_region_size / 2);
                continue;
            }
            if (tumor.records() == 0) continue;
            return Mutect2RegionDecoded{std::move(tumor), std::move(normal),
                                        tumor_sample_names, std::move(tumor_sample_reads),
                                        std::move(normal_sample_names),
                                        std::move(normal_sample_reads),
                                        !options.normal.empty() || !options.normal_sample.empty(),
                                        core.tid, core.start, core.end,
                                        std::move(halo_intervals)};
        }
        return std::nullopt;
    };
    auto compute = [&](Mutect2RegionDecoded decoded) -> std::optional<Mutect2RegionComputed> {
        auto tile_options = options.calling;
        // Mutect2Engine calls one complete assembly/likelihood pipeline per
        // active AssemblyRegion even in tumor-only mode.  Keep the Kokkos
        // graph and PairHMM launches local to that unit; the outer tile only
        // owns decode/output aggregation.
        tile_options.partition_somatic_assembly_regions = true;
        // SomaticGenotypingEngine runs on each region-local PairHMM matrix.
        // Preserve those matrices through the tile merge so the final TLOD
        // gate consumes the Kokkos somatic likelihood result rather than the
        // legacy pileup genotype-summary fallback.
        tile_options.retain_assembly_region_likelihood_results = true;
        // Reproduce aggregate AssemblyRegion scheduling before applying the
        // streaming core at the writer boundary.  `CallingOptions::intervals`
        // is not merely a serialization filter: it controls ActivityProfile
        // construction and region ownership.  Replacing it with an arbitrary
        // tile would change an otherwise identical GATK graph.  The explicit
        // VCF crop below is therefore the only per-tile restriction.
        tile_options.intervals = metadata_reader.requested_intervals();
        tile_options.locus_intervals.clear();
        tile_options.interval_tid = -1;
        tile_options.interval_start = -1;
        tile_options.interval_end = -1;
        auto initial_tumor_options = tile_options;
        auto tumor_result = fastgatk::calling::run(
            decoded.tumor, references, initial_tumor_options);
        auto tumor_somatic = calculate_somatic(
            tumor_result, options.minimum_allele_fraction,
            options.reference_confidence_mode != ReferenceConfidenceMode::None);
        fastgatk::calling::Result empty_normal;
        std::optional<fastgatk::calling::Result> normal_result;
        std::vector<fastgatk::calling::Result> normal_sample_results;
        std::vector<std::optional<fastgatk::kernels::SomaticLikelihoodResult>>
            normal_sample_somatics;
        std::vector<NormalSampleView> normal_views;
        auto normal_options = tile_options;
        // GATK builds one AssemblyRegion graph from all samples and scores
        // both likelihood matrices against that shared haplotype population.
        // Preserve that AssemblyResult ownership here; each sample still has
        // its own read filtering, clipping and PairHMM likelihood pass.
        normal_options.restrict_to_forced_alleles = false;
        if (decoded.has_normal || decoded.tumor_sample_reads.size() > 1) {
            fastgatk::io::ReadBatch joint_reads;
            joint_reads.offsets.push_back(0);
            const auto tumor_role_count = std::max<std::size_t>(
                1, decoded.tumor_sample_names.size());
            if (decoded.tumor.sample_ids.size() != decoded.tumor.records())
                label_batch_sample(decoded.tumor, 0);
            if (decoded.normal.sample_ids.size() != decoded.normal.records())
                label_batch_sample(decoded.normal,
                                   static_cast<std::int32_t>(tumor_role_count));
            append_batch(joint_reads, decoded.tumor);
            append_batch(joint_reads, decoded.normal);
            // Read names may recur across the tumor and normal streams.
            // Assembly input is joint, while read-pair and fragment handling
            // stays sample-local exactly as in GATK's splitReadsBySample().
            auto joint_options = tile_options;
            // Mutect2Engine.isActive() evaluates the tumor pileup, then
            // suppresses a state if the normal supports the same dominant
            // allele.  Let the joint traversal form those regions itself;
            // pinning the tumor-only schedule here loses that source rule.
            joint_options.activity_tumor_sample_id = -1;
            joint_options.activity_normal_sample_id = -1;
            joint_options.activity_tumor_sample_ids.clear();
            joint_options.activity_normal_sample_ids.clear();
            joint_options.activity_tumor_sample_ids.reserve(tumor_role_count);
            for (std::size_t index = 0; index < tumor_role_count; ++index)
                joint_options.activity_tumor_sample_ids.push_back(
                    static_cast<std::int32_t>(index));
            if (decoded.has_normal) {
                const auto normal_role_count = std::max<std::size_t>(
                    1, decoded.normal_sample_names.empty()
                        ? 1 : decoded.normal_sample_names.size());
                joint_options.activity_normal_sample_ids.reserve(normal_role_count);
                for (std::size_t index = 0; index < normal_role_count; ++index)
                    joint_options.activity_normal_sample_ids.push_back(
                        static_cast<std::int32_t>(tumor_role_count + index));
            }
            // A low-quality recovery is an EventMap event owned by one
            // AssemblyRegion, never a reason to flatten every region in a
            // large tile into one graph.  Its linked insertion context is
            // carried by the region-local EventMap/PairHMM path below.  This
            // preserves Mutect2Engine.callRegion() ownership and keeps the
            // Kokkos graph kernel bounded to one GATK ActivityRegion.
            joint_options.partition_somatic_assembly_regions = true;
            joint_options.retain_assembly_region_results =
                joint_options.partition_somatic_assembly_regions;
            // Keep stream and aggregate traversal on the same source-shaped
            // ownership boundary: this joint pass produces AssemblyResults;
            // only the following per-sample passes own likelihood matrices.
            joint_options.assembly_only = true;
            auto joint_assembly = fastgatk::calling::run(
                joint_reads, references, joint_options);
            auto shared_options = tile_options;
            if (!joint_assembly.assembly_region_results.empty()) {
                // Mutect2Engine.callRegion preserves one AssemblyResult per
                // active region, then splitReadsBySample() evaluates tumor
                // and normal likelihoods on that exact graph/EventMap pair.
                // Do not concatenate these regional graph paths: doing so
                // invents cross-region competitors that do not exist in GATK.
                shared_options.shared_assembly_regions =
                    joint_assembly.assembly_region_results;
                retain_eventmap_projection_and_tumor_recovery(
                    shared_options.shared_assembly_regions, tumor_result,
                    tumor_somatic ? &*tumor_somatic : nullptr,
                    effective_tumor_lod_emit(options));
                shared_options.retain_assembly_region_likelihood_results = true;
                shared_options.forced_calling_regions.clear();
                shared_options.forced_calling_regions.reserve(
                    shared_options.shared_assembly_regions.size());
                for (const auto& assembly : shared_options.shared_assembly_regions)
                    shared_options.forced_calling_regions.push_back(assembly.region);
            } else {
                // Historical experimental fallback for one unpartitioned
                // graph.  Region AssemblyResult capture takes this branch out
                // of the GATK-aligned multi-region execution path.
                shared_options.forced_calling_regions = joint_assembly.calling_regions;
                shared_options.forced_alleles = event_map_forced_alleles(
                    joint_assembly, tumor_result,
                    tumor_somatic ? &*tumor_somatic : nullptr,
                    effective_tumor_lod_emit(options));
                shared_options.shared_assembly_graph = assembly_graph(joint_assembly);
                shared_options.restrict_to_forced_alleles = true;
            }
            // shared_options owns the selected AssemblyResult graph/EventMap
            // handoff. Release both superseded likelihood Result payloads
            // before the final tumor/normal PairHMM passes.
            joint_assembly = fastgatk::calling::Result{};
            tumor_result = fastgatk::calling::Result{};
            tumor_somatic.reset();
            tumor_result = fastgatk::calling::run(
                decoded.tumor, references, shared_options);
            tumor_somatic = calculate_somatic(
                tumor_result, options.minimum_allele_fraction,
                options.reference_confidence_mode != ReferenceConfidenceMode::None);
            normal_options = std::move(shared_options);
        }
        if (decoded.has_normal)
            normal_result = fastgatk::calling::run(decoded.normal, references, normal_options);
        const auto normal_somatic = normal_result
            ? calculate_somatic(*normal_result, options.minimum_allele_fraction)
            : std::optional<fastgatk::kernels::SomaticLikelihoodResult>{};
        materialize_low_quality_assembly_recovery_calls(
            tumor_result, tumor_somatic ? &*tumor_somatic : nullptr,
            effective_tumor_lod_emit(options));
        materialize_reference_confidence_eventmap_calls(
            tumor_result, tumor_somatic ? &*tumor_somatic : nullptr,
            options.reference_confidence_mode != ReferenceConfidenceMode::None);
        if (tumor_somatic) apply_grouped_somatic_depth(tumor_result);
        if (normal_somatic && normal_result)
            apply_grouped_somatic_depth(*normal_result);
        std::vector<fastgatk::calling::Result> tumor_sample_results;
        std::vector<std::optional<fastgatk::kernels::SomaticLikelihoodResult>>
            tumor_sample_somatics;
        std::vector<NormalSampleView> tumor_views;
        if (decoded.tumor_sample_reads.size() > 1 &&
            decoded.tumor_sample_names.size() == decoded.tumor_sample_reads.size()) {
            tumor_sample_results.reserve(decoded.tumor_sample_reads.size());
            tumor_sample_somatics.reserve(decoded.tumor_sample_reads.size());
            tumor_views.reserve(decoded.tumor_sample_reads.size());
            for (std::size_t index = 0; index < decoded.tumor_sample_reads.size(); ++index) {
                tumor_sample_results.push_back(fastgatk::calling::run(
                    decoded.tumor_sample_reads[index], references, normal_options));
                tumor_sample_somatics.push_back(calculate_somatic(
                    tumor_sample_results.back(), options.minimum_allele_fraction));
                materialize_low_quality_assembly_recovery_calls(
                    tumor_sample_results.back(),
                    tumor_sample_somatics.back() ? &*tumor_sample_somatics.back() : nullptr,
                    effective_tumor_lod_emit(options));
                materialize_reference_confidence_eventmap_calls(
                    tumor_sample_results.back(),
                    tumor_sample_somatics.back() ? &*tumor_sample_somatics.back() : nullptr,
                    options.reference_confidence_mode != ReferenceConfidenceMode::None);
                if (tumor_sample_somatics.back())
                    apply_grouped_somatic_depth(tumor_sample_results.back());
            }
            for (std::size_t index = 0; index < tumor_sample_results.size(); ++index) {
                tumor_views.push_back(NormalSampleView{
                    decoded.tumor_sample_names[index], &decoded.tumor_sample_reads[index],
                    &tumor_sample_results[index],
                    tumor_sample_somatics[index]
                        ? &*tumor_sample_somatics[index] : nullptr});
            }
        }
        if (decoded.has_normal && !decoded.normal_sample_reads.empty() &&
            decoded.normal_sample_names.size() == decoded.normal_sample_reads.size()) {
            normal_sample_results.reserve(decoded.normal_sample_reads.size());
            normal_sample_somatics.reserve(decoded.normal_sample_reads.size());
            normal_views.reserve(decoded.normal_sample_reads.size());
            for (std::size_t index = 0; index < decoded.normal_sample_reads.size(); ++index) {
                normal_sample_results.push_back(fastgatk::calling::run(
                    decoded.normal_sample_reads[index], references, normal_options));
                normal_sample_somatics.push_back(calculate_somatic(
                    normal_sample_results.back(), options.minimum_allele_fraction));
                if (normal_sample_somatics.back())
                    apply_grouped_somatic_depth(normal_sample_results.back());
            }
            for (std::size_t index = 0; index < normal_sample_results.size(); ++index) {
                normal_views.push_back(NormalSampleView{
                    decoded.normal_sample_names[index], &decoded.normal_sample_reads[index],
                    &normal_sample_results[index],
                    normal_sample_somatics[index]
                        ? &*normal_sample_somatics[index] : nullptr});
            }
        }
        const auto posterior = calculate_posterior(
            options, tumor_result, decoded.tumor,
            normal_result ? &*normal_result : &empty_normal);
        if (tumor_somatic.has_value()) {
            tumor_result.calls.erase(std::remove_if(tumor_result.calls.begin(), tumor_result.calls.end(),
                [&](const auto& call) {
                    const auto index = candidate_index(tumor_result, call.candidate);
                    const auto has_tlod = index < tumor_somatic->tlod.size() &&
                        std::isfinite(tumor_somatic->tlod[index]);
                    const auto passes_tlod = has_tlod &&
                        tumor_somatic->tlod[index] > effective_tumor_lod_emit(options);
                    const auto passes_normal = normal_lod_allows_call(
                        options.normal_log10_odds, call,
                        normal_result ? &*normal_result : nullptr,
                        normal_somatic ? &*normal_somatic : nullptr);
                    return !passes_tlod || !passes_normal;
                }), tumor_result.calls.end());
        }
        // Each decode tile replays the selected AssemblyRegion source
        // window, so its result owns the same complete EventMap/RCM output
        // as aggregate traversal.  Select emitted rows by record start below
        // (which preserves a GVCF block crossing a tile boundary exactly
        // once), but never downgrade the tile to the ordinary VCF writer.
        // That would silently discard <NON_REF> rows and TLOD block headers.
        const auto rendered = somatic_reference_confidence_vcf_text(
            options, merged_header, decoded.tumor, tumor_result,
            decoded.has_normal ? &*normal_result : nullptr,
            tumor_somatic ? &*tumor_somatic : nullptr,
            normal_somatic ? &*normal_somatic : nullptr,
            posterior ? &*posterior : nullptr,
            normal_views.empty() ? nullptr : &normal_views,
            tumor_views.empty() ? nullptr : &tumor_views,
            &references, &feature_resources);
        const auto parts = split_vcf_header_body(rendered);
        Mutect2RegionComputed result;
        const fastgatk::io::HtsInterval output_core{
            decoded.tid, decoded.core_start, decoded.core_end};
        result.body = vcf_body_for_core(
            parts.second, merged_header.contigs[static_cast<std::size_t>(decoded.tid)],
            decoded.core_start, decoded.core_end);
        result.tumor_reads = decoded.tumor.records();
        result.normal_reads = decoded.normal.records();
        const auto is_in_output_core = [&](const auto& candidate) {
            return candidate.tid == output_core.tid &&
                   candidate.position >= output_core.start &&
                   candidate.position < output_core.end;
        };
        result.candidate_sites = std::count_if(
            tumor_result.candidates.begin(), tumor_result.candidates.end(),
            [&](const auto& candidate) { return is_in_output_core(candidate); });
        result.tumor_calls = std::count_if(
            tumor_result.calls.begin(), tumor_result.calls.end(),
            [&](const auto& call) { return is_in_output_core(call.candidate); });
        result.normal_calls = normal_result ? std::count_if(
            normal_result->calls.begin(), normal_result->calls.end(),
            [&](const auto& call) { return is_in_output_core(call.candidate); }) : 0;
        result.pairhmm_pairs = tumor_result.pairhmm_pairs;
        result.itr_artifact_clipped_reads = tumor_result.post_read_filter_transformed_reads;
        result.downsampled_reads = tumor_result.downsampled_reads;
        result.overlapping_pairs = tumor_result.overlapping_pairs;
        result.overlapping_bases = tumor_result.overlapping_bases;
        result.overlapping_conflicting_bases = tumor_result.overlapping_conflicting_bases;
        result.overlapping_quality_caps = tumor_result.overlapping_quality_caps;
        result.overlapping_indel_quality_caps = tumor_result.overlapping_indel_quality_caps;
        if (tumor_somatic) {
            result.somatic_evidence_groups = tumor_somatic->evidence_groups;
            result.somatic_likelihood_prepare_seconds = tumor_somatic->prepare_seconds;
            result.somatic_likelihood_seconds = tumor_somatic->seconds;
            result.somatic_likelihood_execution_space = tumor_somatic->execution_space;
        }
        if (posterior) {
            result.somatic_posterior_prepare_seconds = posterior->prepare_seconds;
            result.somatic_posterior_seconds = posterior->seconds;
            result.somatic_posterior_execution_space = posterior->execution_space;
        }
        if (!standard_f1r2_archive)
            result.f1r2_text = mutect2_region_f1r2_text(
                decoded.tumor, tumor_result, result.f1r2_alt_support, result.r1f2_alt_support,
                &output_core);
        result.f1r2_archive = make_f1r2_archive_sample(
            merged_header, decoded.tumor, tumor_result, references, &output_core);
        if (!options.tumor_sample.empty()) result.f1r2_archive.sample = options.tumor_sample;
        return result;
    };
    auto encode = [](Mutect2RegionComputed computed) -> std::optional<Mutect2RegionComputed> {
        return computed;
    };
    auto sink = [&](Mutect2RegionComputed encoded) {
        writer.write(encoded.body);
        if (f1r2_tsv.is_open()) f1r2_tsv << encoded.f1r2_text;
        merge_f1r2_archive(summary.f1r2_archive, std::move(encoded.f1r2_archive));
        summary.tumor_reads += encoded.tumor_reads;
        summary.normal_reads += encoded.normal_reads;
        summary.candidate_sites += encoded.candidate_sites;
        summary.tumor_calls += encoded.tumor_calls;
        summary.normal_calls += encoded.normal_calls;
        summary.pairhmm_pairs += encoded.pairhmm_pairs;
        summary.itr_artifact_clipped_reads += encoded.itr_artifact_clipped_reads;
        summary.downsampled_reads += encoded.downsampled_reads;
        summary.overlapping_pairs += encoded.overlapping_pairs;
        summary.overlapping_bases += encoded.overlapping_bases;
        summary.overlapping_conflicting_bases += encoded.overlapping_conflicting_bases;
        summary.overlapping_quality_caps += encoded.overlapping_quality_caps;
        summary.overlapping_indel_quality_caps += encoded.overlapping_indel_quality_caps;
        summary.f1r2_alt_support += encoded.f1r2_alt_support;
        summary.r1f2_alt_support += encoded.r1f2_alt_support;
        summary.somatic_evidence_groups += encoded.somatic_evidence_groups;
        summary.somatic_likelihood_prepare_seconds += encoded.somatic_likelihood_prepare_seconds;
        summary.somatic_likelihood_seconds += encoded.somatic_likelihood_seconds;
        summary.somatic_posterior_prepare_seconds += encoded.somatic_posterior_prepare_seconds;
        summary.somatic_posterior_seconds += encoded.somatic_posterior_seconds;
        if (summary.somatic_likelihood_execution_space.empty())
            summary.somatic_likelihood_execution_space = encoded.somatic_likelihood_execution_space;
        if (summary.somatic_posterior_execution_space.empty())
            summary.somatic_posterior_execution_space = encoded.somatic_posterior_execution_space;
        ++summary.regions;
    };
    Pipeline pipeline(
        Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
        std::move(decode), std::move(compute), std::move(encode), std::move(sink),
        mutect2_region_decoded_bytes, mutect2_region_computed_bytes,
        mutect2_region_computed_bytes);
    const auto metrics = pipeline.run();
    writer.close();
    if (f1r2_tsv.is_open()) f1r2_tsv.close();
    if (standard_f1r2_archive) {
        if (summary.f1r2_archive.sample.empty())
            summary.f1r2_archive.sample = options.tumor_sample.empty()
                ? (merged_header.samples.empty() ? "TUMOR" : merged_header.samples.front())
                : options.tumor_sample;
        fastgatk::somatic::write_f1r2_archive(f1r2_path, {summary.f1r2_archive});
    }
    options.stream_indexed = true;
    options.streamed_regions = summary.regions;
    options.streamed_peak_host_bytes = peak_host_bytes.load(std::memory_order_relaxed);
    options.effective_batch_records = decode_batch_records;
    options.adaptive_batch_reductions = adaptive_reductions;
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
    return summary;
}

void write_mutect2_streaming_metadata(
    const Options& options,
    const fastgatk::runtime::ResourceSnapshot& resources,
    const fastgatk::io::HeaderSummary& header,
    const Mutect2StreamSummary& summary,
    const std::string& stats_path,
    const std::string& f1r2_path,
    const std::string& index_path) {
    const auto complete = [](const std::string& path) {
        return file_complete(path);
    };
    std::ofstream stats(stats_path);
    if (!stats) throw std::runtime_error("cannot write Mutect2 stats: " + stats_path);
    stats << std::setprecision(17)
          << "{\"schema_version\":1,\"tool\":\"Mutect2\",\"status\":\"prototype\"," 
          << "\"tumor_candidate_sites\":" << summary.candidate_sites
          << ",\"tumor_calls\":" << summary.tumor_calls
          << ",\"normal_calls\":" << summary.normal_calls
          << ",\"somatic_candidates\":" << summary.tumor_calls
          << ",\"tumor_reads\":" << summary.tumor_reads
          << ",\"normal_reads\":" << summary.normal_reads
          << ",\"tumor_sample\":\""
          << json_escape(options.tumor_sample.empty()
                             ? (header.samples.empty() ? std::string("TUMOR") : header.samples.front())
                             : options.tumor_sample)
          << "\",\"normal_sample\":\""
          << json_escape(options.normal_sample.empty() ? std::string("NORMAL") : options.normal_sample)
          << "\",\"normal_samples\":" << normal_sample_list_json(options)
          << ",\"sample_name_filtering\":"
          << ((!options.tumor_sample.empty() || !options.normal_sample.empty()) ? "true" : "false")
          << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
          << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
          << ",\"itr_artifact_clipping_enabled\":"
          << (options.ignore_itr_artifacts ? "false" : "true")
          << ",\"itr_artifact_clipped_reads\":" << summary.itr_artifact_clipped_reads
          << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
          << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
          << ",\"effective_tumor_lod_emit\":"
          << json_number_or_null(effective_tumor_lod_emit(options))
          << ",\"normal_lod\":" << options.normal_log10_odds
          << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
          << ",\"mitochondria_pruning_lod_threshold\":"
          << options.calling.graph_pruning_log_odds_threshold
          << ",\"f1r2_alt_support\":" << summary.f1r2_alt_support
          << ",\"r1f2_alt_support\":" << summary.r1f2_alt_support
          << ",\"multiallelic_locus_writer\":true"
          << ",\"population_allele_frequency\":" << options.population_allele_frequency
          << ",\"population_allele_frequency_log10\":"
          << -std::log10(options.population_allele_frequency)
          << ",\"f1r2_format\":\""
          << (suffix(f1r2_path, ".tar.gz") ? "collect-f1r2-counts-tar-gz" : "legacy-tsv") << "\""
          << ",\"somatic_likelihood_model\":\"per-read-dirichlet-variational-evidence-v1\""
          << ",\"somatic_evidence_groups\":" << summary.somatic_evidence_groups
          << ",\"somatic_evidence_grouping\":true"
          << ",\"independent_mates\":"
          << (options.independent_mates ? "true" : "false")
          << ",\"somatic_informative_read_overlap_margin\":2"
          << ",\"somatic_grouped_responsibility_depth\":true"
          << ",\"somatic_likelihood_execution_space\":\""
          << json_escape(summary.somatic_likelihood_execution_space) << "\""
          << ",\"somatic_likelihood_prepare_seconds\":" << summary.somatic_likelihood_prepare_seconds
          << ",\"somatic_likelihood_seconds\":" << summary.somatic_likelihood_seconds
          << ",\"somatic_posterior_model\":\"somatic-germline-artifact-orientation-v1\""
          << ",\"somatic_posterior_execution_space\":\""
          << json_escape(summary.somatic_posterior_execution_space) << "\""
          << ",\"somatic_posterior_prepare_seconds\":" << summary.somatic_posterior_prepare_seconds
          << ",\"somatic_posterior_seconds\":" << summary.somatic_posterior_seconds
          << ",\"requested_batch_records\":" << options.requested_batch_records
          << ",\"effective_batch_records\":" << options.effective_batch_records
          << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
          << ",\"stream_by_region\":true,\"streamed_regions\":" << options.streamed_regions
          << ",\"streamed_peak_host_bytes\":" << options.streamed_peak_host_bytes
          << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\""
          << ",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
          << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
          << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
          << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
          << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
          << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
          << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
          << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
          << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
          << ",\"pipeline_stage_capacity_bytes\":" << options.pipeline_stage_capacity_bytes
          << ",\"pairhmm_pairs\":" << summary.pairhmm_pairs
          << ",\"downsampled_reads\":" << summary.downsampled_reads
          << ",\"overlapping_pairs\":" << summary.overlapping_pairs
          << ",\"overlapping_bases\":" << summary.overlapping_bases
          << ",\"overlapping_conflicting_bases\":" << summary.overlapping_conflicting_bases
          << ",\"overlapping_quality_caps\":" << summary.overlapping_quality_caps
          << ",\"overlapping_indel_quality_caps\":" << summary.overlapping_indel_quality_caps
          << ",\"overlapping_pcr_snv_quality\":"
          << static_cast<unsigned>(options.calling.activity_pcr_snv_quality)
          << ",\"overlapping_pcr_indel_quality\":"
          << static_cast<unsigned>(options.calling.overlapping_pcr_indel_quality)
          << ",\"sites_only_vcf_output\":"
          << (options.sites_only_vcf_output ? "true" : "false")
          << ",\"create_output_variant_index\":"
          << (options.create_output_variant_index ? "true" : "false")
          << ",\"bit_identical_to_gatk\":false}\n";
    stats.flush();
    if (!complete(stats_path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Mutect2 stats is empty");

    const auto manifest_path = options.manifest.empty()
        ? options.output + ".manifest.json" : options.manifest;
    std::ofstream manifest(manifest_path);
    if (!manifest) throw std::runtime_error("cannot write Mutect2 manifest: " + manifest_path);
    manifest << "{\"schema_version\":1,\"tool\":\"Mutect2\","
             << "\"implementation\":\"fastgatk-mutect2\",\"status\":\"prototype\","
             << "\"primary_output\":\"" << json_escape(options.output)
             << "\",\"primary_output_kind\":\"vcf\",\"compatibility\":{"
             << "\"tumor_normal\":"
             << ((!options.normal.empty() || !options.normal_sample.empty()) ? "true" : "false")
             << ",\"normal_samples\":" << normal_sample_list_json(options)
             << ",\"somatic_filter\":true,\"f1r2\":true,\"multiallelic_locus_writer\":true"
             << ",\"somatic_read_likelihoods\":true,\"somatic_fragment_grouping\":true"
             << ",\"independent_mates\":"
             << (options.independent_mates ? "true" : "false")
             << ",\"somatic_informative_read_overlap\":true"
             << ",\"somatic_grouped_responsibility_depth\":true"
             << ",\"somatic_pileup_candidate_boundary\":true,\"somatic_posterior\":true"
             << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
             << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
             << ",\"itr_artifact_clipping_enabled\":"
             << (options.ignore_itr_artifacts ? "false" : "true")
             << ",\"mitochondria_defaults_applied\":" << (options.mitochondria_mode ? "true" : "false")
             << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
             << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
             << ",\"effective_tumor_lod_emit\":"
             << json_number_or_null(effective_tumor_lod_emit(options))
             << ",\"normal_lod\":" << options.normal_log10_odds
             << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
             << ",\"orientation_bias_posterior\":true,\"contamination_posterior\":true"
             << ",\"read_haplotype_request_deduplication\":true"
             << ",\"stream_by_region\":true,\"indexed_input\":true"
             << ",\"f1r2_standard_tar\":" << (suffix(f1r2_path, ".tar.gz") ? "true" : "false")
             << ",\"f1r2_legacy_tsv\":" << (suffix(f1r2_path, ".tar.gz") ? "false" : "true")
             << ",\"vcf_index\":" << (index_path.empty() ? "false" : "true")
             << ",\"create_output_variant_index\":"
             << (options.create_output_variant_index ? "true" : "false")
             << ",\"bit_identical_to_gatk\":false},\"outputs\":["
             << "{\"path\":\"" << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":"
             << (complete(options.output) ? "true" : "false") << "}"
             << ", {\"path\":\"" << json_escape(stats_path)
             << "\",\"kind\":\"mutect2-stats\",\"complete\":"
             << (complete(stats_path) ? "true" : "false") << "}"
             << ", {\"path\":\"" << json_escape(f1r2_path)
             << "\",\"kind\":\""
             << (suffix(f1r2_path, ".tar.gz") ? "f1r2-collect-tar-gz" : "f1r2-orientation")
             << "\",\"complete\":" << (complete(f1r2_path) ? "true" : "false") << "}"
             << (index_path.empty() ? "" : ", {\"path\":\"" + json_escape(index_path) +
                 "\",\"kind\":\"vcf-index\",\"complete\":" +
                 (complete(index_path) ? "true}" : "false}"))
             << "],\"telemetry\":{\"resources\":" << resources.to_json()
             << ",\"requested_pairhmm_threads\":" << options.threads
             << ",\"effective_pairhmm_threads\":"
             << resources.effective_threads(static_cast<std::size_t>(options.threads))
             << ",\"interval_list_inputs\":" << options.interval_file_inputs
             << ",\"interval_list_records\":" << options.interval_file_records
             << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
             << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
             << ",\"itr_artifact_clipping_enabled\":"
             << (options.ignore_itr_artifacts ? "false" : "true")
             << ",\"itr_artifact_clipped_reads\":" << summary.itr_artifact_clipped_reads
             << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
             << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
             << ",\"effective_tumor_lod_emit\":"
             << json_number_or_null(effective_tumor_lod_emit(options))
             << ",\"normal_lod\":" << options.normal_log10_odds
             << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
             << ",\"population_allele_frequency\":" << options.population_allele_frequency
             << ",\"mitochondria_pruning_lod_threshold\":"
             << options.calling.graph_pruning_log_odds_threshold
             << ",\"tumor_reads\":" << summary.tumor_reads
             << ",\"normal_reads\":" << summary.normal_reads
             << ",\"tumor_candidate_sites\":" << summary.candidate_sites
             << ",\"tumor_calls\":" << summary.tumor_calls
             << ",\"normal_calls\":" << summary.normal_calls
             << ",\"stream_by_region\":true,\"streamed_regions\":" << options.streamed_regions
             << ",\"stream_indexed\":true,\"streamed_peak_host_bytes\":"
             << options.streamed_peak_host_bytes
             << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\""
             << ",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
             << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
             << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
             << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
             << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
             << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
             << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
             << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
             << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
             << ",\"pipeline_stage_capacity_bytes\":" << options.pipeline_stage_capacity_bytes
             << ",\"somatic_likelihood_execution_space\":\""
             << json_escape(summary.somatic_likelihood_execution_space)
             << "\",\"somatic_posterior_execution_space\":\""
             << json_escape(summary.somatic_posterior_execution_space)
             << "\",\"somatic_evidence_groups\":" << summary.somatic_evidence_groups
             << ",\"somatic_likelihood_prepare_seconds\":" << summary.somatic_likelihood_prepare_seconds
             << ",\"somatic_likelihood_seconds\":" << summary.somatic_likelihood_seconds
             << ",\"somatic_posterior_prepare_seconds\":" << summary.somatic_posterior_prepare_seconds
             << ",\"somatic_posterior_seconds\":" << summary.somatic_posterior_seconds
             << ",\"sites_only_vcf_output\":"
             << (options.sites_only_vcf_output ? "true" : "false")
             << ",\"create_output_variant_index\":"
             << (options.create_output_variant_index ? "true" : "false")
             << ",\"pairhmm_pairs\":" << summary.pairhmm_pairs
             << ",\"f1r2_alt_support\":" << summary.f1r2_alt_support
             << ",\"r1f2_alt_support\":" << summary.r1f2_alt_support
             << ",\"requested_batch_records\":" << options.requested_batch_records
             << ",\"effective_batch_records\":" << options.effective_batch_records
             << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
             << ",\"bit_identical_to_gatk\":false}}\n";
    manifest.flush();
    if (!complete(manifest_path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Mutect2 manifest is empty");
}

}  // namespace

#if !FASTGATK_HAS_HTSLIB
int main() {
    std::cerr << "error: BACKEND_UNAVAILABLE: build with HTSlib for Mutect2\n";
    return 69;
}
#else

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        auto options = parse(argc, argv);
        options.requested_batch_records = options.batch_records;
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        if (options.stream_region_size > 0) {
            // The region streamer deliberately gives each bounded core its
            // own halo and ActivityProfile evaluation.  GATK's
            // --assembly-region-out is one continuous Walker traversal;
            // concatenating independently band-pass-filtered tile profiles
            // would manufacture boundaries that the source walker never
            // observed.  Reject this combination until the streamer owns a
            // cross-tile ActivityProfile state machine, rather than silently
            // succeeding with a missing or misleading IGV file.
            if (!options.assembly_region_out.empty())
                throw std::runtime_error(
                    "BACKEND_UNAVAILABLE: --assembly-region-out requires aggregate Mutect2 traversal");
            const auto tumor_inputs = options.tumor_inputs.empty()
                ? std::vector<std::string>{options.tumor} : options.tumor_inputs;
            // Probe every shard before constructing the metadata reader.  This
            // matches GATK's sequence-dictionary validation and guarantees
            // that an unindexed shard cannot be silently omitted from a tile.
            fastgatk::io::HeaderSummary merged_stream_header;
            std::size_t stream_interval_inputs = 0;
            std::size_t stream_interval_records = 0;
            for (const auto& path : tumor_inputs) {
                fastgatk::io::HtsReader probe(
                    path, options.reference, options.regions,
                    options.batch_records, {}, {}, 0, 0,
                    options.interval_set_rule);
                if (!probe.has_index())
                    throw std::runtime_error(
                        "BACKEND_UNAVAILABLE: --stream-by-region requires an indexed BAM/CRAM input");
                merge_input_headers(merged_stream_header, probe.header());
                stream_interval_inputs += probe.interval_file_inputs();
                stream_interval_records += probe.interval_file_records();
            }
            // Resolve Mutect2's default tumor sample against the merged header
            // rather than whichever shard happened to be argv[0].
            if (options.tumor_sample.empty()) {
                if (merged_stream_header.samples.empty())
                    throw std::runtime_error(
                        "BAD_INPUT: tumor input header has no @RG SM sample (GATK requires a sample)");
                options.tumor_sample = merged_stream_header.samples.front();
            }
            if (std::find(merged_stream_header.samples.begin(),
                          merged_stream_header.samples.end(), options.tumor_sample) ==
                merged_stream_header.samples.end())
                throw std::runtime_error(
                    "BAD_INPUT: requested tumor sample is absent from tumor input headers: " +
                    options.tumor_sample);
            const auto requested_normal_samples = selected_normal_samples(options);
            if (options.normal.empty()) {
                for (const auto& sample : requested_normal_samples) {
                    if (std::find(merged_stream_header.samples.begin(),
                                  merged_stream_header.samples.end(), sample) ==
                        merged_stream_header.samples.end())
                        throw std::runtime_error(
                            "BAD_INPUT: requested normal sample is absent from tumor input headers: " +
                            sample);
                }
            }
            if (!options.normal.empty()) {
                fastgatk::io::HtsReader normal_probe(
                    options.normal, options.reference, options.regions,
                    options.batch_records, {}, {}, 0, 0,
                    options.interval_set_rule);
                if (!normal_probe.has_index())
                    throw std::runtime_error(
                        "BACKEND_UNAVAILABLE: --stream-by-region requires an indexed BAM/CRAM input");
                auto dictionary_check = merged_stream_header;
                merge_input_headers(dictionary_check, normal_probe.header());
                for (const auto& sample : requested_normal_samples) {
                    if (std::find(normal_probe.header().samples.begin(),
                                  normal_probe.header().samples.end(), sample) ==
                        normal_probe.header().samples.end())
                        throw std::runtime_error(
                            "BAD_INPUT: requested normal sample is absent from normal input header: " +
                            sample);
                }
            }
            options.interval_file_inputs = stream_interval_inputs;
            options.interval_file_records = stream_interval_records;
            fastgatk::io::HtsReader metadata_reader(
                tumor_inputs.front(), options.reference, options.regions,
                options.batch_records, {}, {}, 0, 0,
                options.interval_set_rule);
            const auto references = load_reference(options.reference, merged_stream_header);
            const auto feature_resources = load_mutect_feature_resources(options);
            configure_mutect2_activity_feature_masks(
                options.calling, options, merged_stream_header, feature_resources,
                !options.normal.empty() || !options.normal_sample.empty());
            const auto stream_summary = run_region_streaming(
                options, metadata_reader, merged_stream_header, tumor_inputs,
                references, feature_resources, resources);
            std::string index;
            if (options.create_output_variant_index) {
                if (suffix(options.output, ".gz")) {
                    index = options.output + ".tbi";
                    if (tbx_index_build3(options.output.c_str(), index.c_str(), 0, 0,
                                         &tbx_conf_vcf) != 0)
                        throw std::runtime_error("cannot build Mutect2 VCF index: " + options.output);
                } else if (suffix(options.output, ".vcf")) {
                    index = options.output + ".idx";
                    fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index);
                }
            }
            const auto stats_path = options.stats.empty()
                ? options.output + ".stats.json" : options.stats;
            const auto f1r2_path = options.f1r2.empty()
                ? options.output + ".f1r2.tsv" : options.f1r2;
            write_mutect2_streaming_metadata(
                options, resources, metadata_reader.header(), stream_summary,
                stats_path, f1r2_path, index);
            std::cout << "{\"tool\":\"Mutect2\",\"status\":\"prototype\","
                      << "\"stream_by_region\":true,\"streamed_regions\":"
                      << stream_summary.regions << ",\"tumor_reads\":"
                      << stream_summary.tumor_reads << ",\"tumor_calls\":"
                      << stream_summary.tumor_calls << ",\"normal_calls\":"
                      << stream_summary.normal_calls << ",\"pipeline_lifecycle\":\""
                      << "Host decode->bounded queue->Kokkos compute->encode->sink\","
                      << "\"pipeline_decoded_items\":" << options.pipeline_decoded_items
                      << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
                      << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
                      << ",\"pipeline_stage_capacity_bytes\":"
                      << options.pipeline_stage_capacity_bytes << "}\n";
            Kokkos::finalize();
            initialized = false;
            return 0;
        }
        fastgatk::io::HeaderSummary tumor_header;
        std::vector<fastgatk::io::HtsInterval> requested_intervals;
        // GATK's AssemblyRegionWalker reads the assembly padding around an
        // explicit -L selector.  Discover the normalized core first, then
        // use an expanded Host fetch only for decoding.  CallingOptions keeps
        // requested_intervals below, so no candidate outside -L can be
        // serialized merely because its read contributed graph context.
        fastgatk::io::HtsReader tumor_interval_metadata(
            options.tumor_inputs.front(), options.reference, options.regions,
            options.batch_records, {}, {}, 0, 0, options.interval_set_rule);
        const auto assembly_read_regions = assembly_fetch_regions(
            tumor_interval_metadata.header(), tumor_interval_metadata.intervals(),
            options.calling.assembly_region_padding);
        const auto* assembly_read_region_selectors = assembly_read_regions.empty()
            ? nullptr : &assembly_read_regions;
        std::string tumor_sample_effective;
        auto tumor_reads = read_all_many(options.tumor_inputs, options.reference, options.regions,
                                         options.batch_records, tumor_header,
                                         resources.safe_memory_budget_bytes(), options.tumor_sample,
                                         true, &tumor_sample_effective,
                                         &options.interval_file_inputs, &options.interval_file_records,
                                         &options.adaptive_batch_reductions,
                                         &options.effective_batch_records,
                                         &requested_intervals, true, options.interval_set_rule,
                                         assembly_read_region_selectors);
        // ReadsPathDataSource uses a coordinate comparator when multiple
        // indexed inputs are supplied.  Re-establish that order before the
        // shared overlap/fragment stage; single-input HTSlib order is already
        // sorted, while repeated shards need a deterministic merge.
        sort_batch_by_coordinate(tumor_reads);
        if (!tumor_sample_effective.empty()) options.tumor_sample = tumor_sample_effective;
        const auto tumor_sample_names = selected_tumor_samples(options, tumor_header);
        if (tumor_sample_names.empty())
            throw std::runtime_error(
                "BAD_INPUT: every tumor input sample was selected as a normal sample");
        if (options.reference_confidence_mode != ReferenceConfidenceMode::None &&
            tumor_sample_names.size() != 1)
            throw std::runtime_error(
                "BAD_INPUT: Mutect2 reference confidence requires exactly one sample; "
                "multi-tumor inputs are unsupported by GATK 4.6.2");
        // A legacy single-tumor invocation keeps its established decode path
        // byte-for-byte.  With multiple non-normal samples, decode each one
        // against its source-local @RG namespace and then form the aggregate
        // tumor matrix explicitly.  Filtering a single --tumor-sample here
        // would silently discard valid tumor evidence before Kokkos PairHMM.
        std::vector<fastgatk::io::ReadBatch> tumor_sample_reads;
        if (tumor_sample_names.size() > 1) {
            tumor_reads = fastgatk::io::ReadBatch{};
            tumor_reads.offsets.push_back(0);
            tumor_sample_reads.reserve(tumor_sample_names.size());
            options.interval_file_inputs = 0;
            options.interval_file_records = 0;
            options.adaptive_batch_reductions = 0;
            options.effective_batch_records = 0;
            bool captured_intervals = false;
            for (std::size_t sample_index = 0;
                 sample_index < tumor_sample_names.size(); ++sample_index) {
                const auto& sample = tumor_sample_names[sample_index];
                fastgatk::io::HeaderSummary sample_header;
                std::size_t sample_interval_inputs = 0;
                std::size_t sample_interval_records = 0;
                std::uint64_t sample_reductions = 0;
                std::size_t sample_effective_batch = 0;
                std::vector<fastgatk::io::HtsInterval> sample_intervals;
                auto sample_reads = read_all_many(
                    options.tumor_inputs, options.reference, options.regions,
                    options.batch_records, sample_header,
                    resources.safe_memory_budget_bytes(), sample, false, nullptr,
                    &sample_interval_inputs, &sample_interval_records,
                    &sample_reductions, &sample_effective_batch,
                    &sample_intervals, true, options.interval_set_rule,
                    assembly_read_region_selectors);
                sort_batch_by_coordinate(sample_reads);
                label_batch_sample(sample_reads,
                                   static_cast<std::int32_t>(sample_index));
                append_batch(tumor_reads, sample_reads);
                tumor_sample_reads.push_back(std::move(sample_reads));
                options.interval_file_inputs += sample_interval_inputs;
                options.interval_file_records += sample_interval_records;
                options.adaptive_batch_reductions += sample_reductions;
                options.effective_batch_records = options.effective_batch_records == 0
                    ? sample_effective_batch
                    : std::min(options.effective_batch_records, sample_effective_batch);
                if (!captured_intervals) {
                    requested_intervals = std::move(sample_intervals);
                    captured_intervals = true;
                }
                tumor_header = std::move(sample_header);
            }
            sort_batch_by_coordinate(tumor_reads);
        }
        // PairHMM/assembly uses a halo, but Mutect2's output domain is the
        // requested core interval.  Feed the raw selectors to the shared
        // caller so graph candidates projected from the halo cannot leak past
        // -L (the previous path could emit a site at 70052 for -L 69000-70000).
        options.calling.intervals = requested_intervals;
        if (!options.tumor_sample.empty() &&
            std::find(tumor_header.samples.begin(), tumor_header.samples.end(), options.tumor_sample) ==
                tumor_header.samples.end())
            throw std::runtime_error("BAD_INPUT: requested tumor sample is absent from tumor input header: " +
                                     options.tumor_sample);
        const auto references = load_reference(options.reference, tumor_header);
        const auto feature_resources = load_mutect_feature_resources(options);
        configure_mutect2_activity_feature_masks(
            options.calling, options, tumor_header, feature_resources,
            !options.normal.empty() || !options.normal_sample.empty());
        configure_mutect2_itr_artifact_transform(options, references);
        // The non-streaming path has the same AssemblyRegion ownership rule
        // as the tiled path above.  Without this, a full -L interval feeds
        // several active islands to one compact graph and creates KBest
        // haplotypes that GATK never presents to PairHMM.
        options.calling.partition_somatic_assembly_regions = true;
        options.calling.retain_assembly_region_likelihood_results = true;
        auto initial_tumor_options = options.calling;
        auto tumor_result = fastgatk::calling::run(
            tumor_reads, references, initial_tumor_options);
        auto tumor_somatic = calculate_somatic(
            tumor_result, options.minimum_allele_fraction,
            options.reference_confidence_mode != ReferenceConfidenceMode::None);
        auto normal_options = options.calling;
        // Mutect2Engine shares one AssemblyResultSet across samples, then
        // evaluates each sample's likelihoods against its region-local
        // graph/EventMap. Keep that ownership in the native aggregate path.
        normal_options.restrict_to_forced_alleles = false;
        std::optional<fastgatk::calling::Result> normal_result;
        // Keep the aggregate normal batch outside the input branch.  A joint
        // AssemblyRegion appends it after the tumor batch, which gives the
        // post-filter source ordinals an unambiguous sample boundary.
        fastgatk::io::ReadBatch normal_reads;
        std::vector<std::string> normal_sample_names;
        std::vector<fastgatk::io::ReadBatch> normal_sample_reads;
        std::size_t normal_read_count = 0;
        const auto tumor_role_count = std::max<std::size_t>(1, tumor_sample_names.size());
        if (!options.normal.empty() || !options.normal_sample.empty()) {
            fastgatk::io::HeaderSummary normal_header;
            const auto normal_inputs = options.normal.empty()
                ? options.tumor_inputs : std::vector<std::string>{options.normal};
            const auto requested_normal_samples = selected_normal_samples(options);
            if (requested_normal_samples.empty()) {
                normal_reads = read_all_many(normal_inputs, options.reference, options.regions,
                                             options.batch_records, normal_header,
                                             resources.safe_memory_budget_bytes(), {}, false, nullptr,
                                             nullptr, nullptr, &options.adaptive_batch_reductions,
                                             &options.effective_batch_records, nullptr, false,
                                             options.interval_set_rule,
                                             assembly_read_region_selectors);
                sort_batch_by_coordinate(normal_reads);
                label_batch_sample(normal_reads,
                                   static_cast<std::int32_t>(tumor_role_count));
            } else {
                normal_reads.offsets.push_back(0);
                normal_sample_names.reserve(requested_normal_samples.size());
                normal_sample_reads.reserve(requested_normal_samples.size());
                for (std::size_t normal_index = 0;
                     normal_index < requested_normal_samples.size(); ++normal_index) {
                    const auto& sample = requested_normal_samples[normal_index];
                    fastgatk::io::ReadBatch sample_reads;
                    fastgatk::io::HeaderSummary sample_header;
                    std::uint64_t sample_reductions = 0;
                    std::size_t sample_effective_batch = options.effective_batch_records;
                    sample_reads = read_all_many(
                        normal_inputs, options.reference, options.regions,
                        options.batch_records, sample_header,
                        resources.safe_memory_budget_bytes(), sample, false, nullptr,
                        nullptr, nullptr, &sample_reductions, &sample_effective_batch,
                        // A selected matched-normal sample may simply have
                        // no reads in this interval.  GATK still carries that
                        // sample through the joint likelihood/output matrix;
                        // it is not a malformed input.  The aggregate normal
                        // batch remains non-empty when other selected normals
                        // overlap, while this sample-local batch is allowed
                        // to be empty just as in the streaming path.
                        nullptr, true, options.interval_set_rule,
                        assembly_read_region_selectors);
                    sort_batch_by_coordinate(sample_reads);
                    label_batch_sample(sample_reads, static_cast<std::int32_t>(
                        tumor_role_count + normal_index));
                    append_batch(normal_reads, sample_reads);
                    normal_sample_names.push_back(sample);
                    normal_sample_reads.push_back(std::move(sample_reads));
                    options.adaptive_batch_reductions += sample_reductions;
                    options.effective_batch_records = std::min(options.effective_batch_records,
                                                               sample_effective_batch);
                }
            }
        }
        normal_read_count = normal_reads.records();
        if (!options.normal.empty() || !options.normal_sample.empty() ||
            tumor_sample_reads.size() > 1) {
            {
                fastgatk::io::ReadBatch joint_reads;
                joint_reads.offsets.push_back(0);
                if (tumor_reads.sample_ids.size() != tumor_reads.records())
                    label_batch_sample(tumor_reads, 0);
                if (normal_reads.sample_ids.size() != normal_reads.records())
                    label_batch_sample(normal_reads,
                                       static_cast<std::int32_t>(tumor_role_count));
                append_batch(joint_reads, tumor_reads);
                append_batch(joint_reads, normal_reads);
                auto joint_options = options.calling;
                // Mutect2Engine.isActive() uses tumor evidence to choose its
                // leading allele and then applies its matched-normal gate.
                // Preserve that joint traversal rather than reusing the
                // tumor-only region schedule.
                joint_options.activity_tumor_sample_id = -1;
                joint_options.activity_normal_sample_id = -1;
                joint_options.activity_tumor_sample_ids.clear();
                joint_options.activity_normal_sample_ids.clear();
                joint_options.activity_tumor_sample_ids.reserve(tumor_role_count);
                for (std::size_t index = 0; index < tumor_role_count; ++index)
                    joint_options.activity_tumor_sample_ids.push_back(
                        static_cast<std::int32_t>(index));
                if (!options.normal.empty() || !options.normal_sample.empty()) {
                    const auto normal_role_count = std::max<std::size_t>(
                        1, normal_sample_names.empty() ? 1 : normal_sample_names.size());
                    joint_options.activity_normal_sample_ids.reserve(normal_role_count);
                    for (std::size_t index = 0; index < normal_role_count; ++index)
                        joint_options.activity_normal_sample_ids.push_back(
                            static_cast<std::int32_t>(tumor_role_count + index));
                }
            // Mutect2Engine.callRegion() owns one AssemblyResult per
            // ActivityRegion.  Low-BQ recovery remains an event on its
            // owning regional graph; it must not collapse an entire tile
            // into a cross-region graph just because one recovery exists.
            joint_options.partition_somatic_assembly_regions = true;
                joint_options.retain_assembly_region_results =
                    joint_options.partition_somatic_assembly_regions;
                // The joint traversal is Mutect2Engine's AssemblyResultSet
                // producer.  Its graph/EventMap is replayed below against
                // separate tumor and normal PairHMM matrices; no caller reads
                // the joint likelihoods, so skip that otherwise duplicated
                // Kokkos numerical pass without changing the Host ownership
                // or the subsequent per-sample scoring contract.
                joint_options.assembly_only = true;
                auto joint_assembly = fastgatk::calling::run(
                    joint_reads, references, joint_options);
                auto shared_options = options.calling;
                if (!joint_assembly.assembly_region_results.empty()) {
                    shared_options.shared_assembly_regions =
                        joint_assembly.assembly_region_results;
                    retain_eventmap_projection_and_tumor_recovery(
                        shared_options.shared_assembly_regions, tumor_result,
                        tumor_somatic ? &*tumor_somatic : nullptr,
                        effective_tumor_lod_emit(options));
                    shared_options.retain_assembly_region_likelihood_results = true;
                    shared_options.forced_calling_regions.clear();
                    shared_options.forced_calling_regions.reserve(
                        shared_options.shared_assembly_regions.size());
                    for (const auto& assembly : shared_options.shared_assembly_regions)
                        shared_options.forced_calling_regions.push_back(assembly.region);
                    // The normal sample receives immutable AssemblyResults
                    // from the joint traversal.  Its read batch therefore
                    // needs only AssemblyRegion windows (including their
                    // halo), not a full-BAM Host ActivityProfile prepass.
                    // This leaves C++ ownership of region/read mapping and
                    // invokes the same Kokkos PairHMM kernels per region.
                    normal_reads = select_reads_for_assembly_regions(
                        normal_reads, shared_options.shared_assembly_regions);
                    for (auto& sample_reads : tumor_sample_reads)
                        sample_reads = select_reads_for_assembly_regions(
                            sample_reads, shared_options.shared_assembly_regions);
                    for (auto& sample_reads : normal_sample_reads)
                        sample_reads = select_reads_for_assembly_regions(
                            sample_reads, shared_options.shared_assembly_regions);
                } else {
                    // Backwards-compatible single AssemblyResult fallback.
                    // The region-scoped branch above is the source-aligned
                    // Mutect2 path and deliberately does not flatten graphs.
                    shared_options.forced_calling_regions = joint_assembly.calling_regions;
                    shared_options.forced_alleles = event_map_forced_alleles(
                        joint_assembly, tumor_result,
                        tumor_somatic ? &*tumor_somatic : nullptr,
                        effective_tumor_lod_emit(options));
                    shared_options.shared_assembly_graph = assembly_graph(joint_assembly);
                    shared_options.restrict_to_forced_alleles = true;
                }
                // The graph/EventMap contract was copied into shared_options;
                // retaining the completed joint and preliminary-tumor
                // likelihood matrices through rescoring only doubles Host
                // residency and changes no GATK state.
                joint_assembly = fastgatk::calling::Result{};
                tumor_result = fastgatk::calling::Result{};
                tumor_somatic.reset();
                tumor_result = fastgatk::calling::run(tumor_reads, references, shared_options);
                tumor_somatic = calculate_somatic(
                    tumor_result, options.minimum_allele_fraction,
                    options.reference_confidence_mode != ReferenceConfidenceMode::None);
                normal_options = std::move(shared_options);
            }
            if (!options.normal.empty() || !options.normal_sample.empty())
                normal_result = fastgatk::calling::run(normal_reads, references, normal_options);
        }
        write_assembly_region_igv(options.assembly_region_out, tumor_header, tumor_result);
        const auto normal_somatic = normal_result
            ? calculate_somatic(*normal_result, options.minimum_allele_fraction)
            : std::optional<fastgatk::kernels::SomaticLikelihoodResult>{};
        materialize_low_quality_assembly_recovery_calls(
            tumor_result, tumor_somatic ? &*tumor_somatic : nullptr,
            effective_tumor_lod_emit(options));
        materialize_reference_confidence_eventmap_calls(
            tumor_result, tumor_somatic ? &*tumor_somatic : nullptr,
            options.reference_confidence_mode != ReferenceConfidenceMode::None);
        if (tumor_somatic) apply_grouped_somatic_depth(tumor_result);
        if (normal_somatic && normal_result)
            apply_grouped_somatic_depth(*normal_result);
        std::vector<fastgatk::calling::Result> tumor_sample_results;
        std::vector<std::optional<fastgatk::kernels::SomaticLikelihoodResult>>
            tumor_sample_somatics;
        std::vector<NormalSampleView> tumor_views;
        if (tumor_sample_reads.size() > 1 &&
            tumor_sample_names.size() == tumor_sample_reads.size()) {
            tumor_sample_results.reserve(tumor_sample_reads.size());
            tumor_sample_somatics.reserve(tumor_sample_reads.size());
            tumor_views.reserve(tumor_sample_reads.size());
            for (std::size_t index = 0; index < tumor_sample_reads.size(); ++index) {
                tumor_sample_results.push_back(fastgatk::calling::run(
                    tumor_sample_reads[index], references, normal_options));
                tumor_sample_somatics.push_back(calculate_somatic(
                    tumor_sample_results.back(), options.minimum_allele_fraction));
                materialize_low_quality_assembly_recovery_calls(
                    tumor_sample_results.back(),
                    tumor_sample_somatics.back() ? &*tumor_sample_somatics.back() : nullptr,
                    effective_tumor_lod_emit(options));
                materialize_reference_confidence_eventmap_calls(
                    tumor_sample_results.back(),
                    tumor_sample_somatics.back() ? &*tumor_sample_somatics.back() : nullptr,
                    options.reference_confidence_mode != ReferenceConfidenceMode::None);
                if (tumor_sample_somatics.back())
                    apply_grouped_somatic_depth(tumor_sample_results.back());
            }
            for (std::size_t index = 0; index < tumor_sample_results.size(); ++index) {
                tumor_views.push_back(NormalSampleView{
                    tumor_sample_names[index], &tumor_sample_reads[index],
                    &tumor_sample_results[index],
                    tumor_sample_somatics[index]
                        ? &*tumor_sample_somatics[index] : nullptr});
            }
        }
        std::vector<fastgatk::calling::Result> normal_sample_results;
        std::vector<std::optional<fastgatk::kernels::SomaticLikelihoodResult>>
            normal_sample_somatics;
        std::vector<NormalSampleView> normal_views;
        if (normal_result && !normal_sample_reads.empty() &&
            normal_sample_names.size() == normal_sample_reads.size()) {
            normal_sample_results.reserve(normal_sample_reads.size());
            normal_sample_somatics.reserve(normal_sample_reads.size());
            normal_views.reserve(normal_sample_reads.size());
            for (std::size_t index = 0; index < normal_sample_reads.size(); ++index) {
                normal_sample_results.push_back(fastgatk::calling::run(
                    normal_sample_reads[index], references, normal_options));
                normal_sample_somatics.push_back(calculate_somatic(
                    normal_sample_results.back(), options.minimum_allele_fraction));
                if (normal_sample_somatics.back())
                    apply_grouped_somatic_depth(normal_sample_results.back());
            }
            for (std::size_t index = 0; index < normal_sample_results.size(); ++index) {
                normal_views.push_back(NormalSampleView{
                    normal_sample_names[index], &normal_sample_reads[index],
                    &normal_sample_results[index],
                    normal_sample_somatics[index]
                        ? &*normal_sample_somatics[index] : nullptr});
            }
        }
        const auto posterior = calculate_posterior(
            options, tumor_result, tumor_reads, normal_result ? &*normal_result : nullptr);
        if (tumor_somatic.has_value()) {
            // Apply the native compatibility emission gate (the stricter of
            // GATK's initial activity and final emission thresholds) only after all
            // likelihood/posterior arrays have been materialized.  Keeping
            // the arrays intact lets stats/telemetry explain rejected
            // candidates while the VCF/F1R2/standard archive contain exactly
            // the records handed to FilterMutectCalls.
            tumor_result.calls.erase(std::remove_if(tumor_result.calls.begin(), tumor_result.calls.end(),
                [&](const auto& call) {
                    const auto index = candidate_index(tumor_result, call.candidate);
                    const auto has_tlod = index < tumor_somatic->tlod.size() &&
                        std::isfinite(tumor_somatic->tlod[index]);
                    const auto passes_tlod = has_tlod &&
                        tumor_somatic->tlod[index] > effective_tumor_lod_emit(options);
                    const auto passes_normal = normal_lod_allows_call(
                        options.normal_log10_odds, call,
                        normal_result ? &*normal_result : nullptr,
                        normal_somatic ? &*normal_somatic : nullptr);
                    return !passes_tlod || !passes_normal;
                }), tumor_result.calls.end());
        }
        const auto vcf = somatic_reference_confidence_vcf_text(
            options, tumor_header, tumor_reads, tumor_result,
            normal_result ? &*normal_result : nullptr,
            tumor_somatic ? &*tumor_somatic : nullptr,
            normal_somatic ? &*normal_somatic : nullptr,
            posterior ? &*posterior : nullptr,
            normal_views.empty() ? nullptr : &normal_views,
            tumor_views.empty() ? nullptr : &tumor_views,
            &references, &feature_resources);
        write_text(options.output, vcf);
        std::string index;
        if (options.create_output_variant_index) {
            if (suffix(options.output, ".gz")) {
                index = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("cannot build Mutect2 VCF index: " + options.output);
            } else if (suffix(options.output, ".vcf")) {
                index = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index);
            }
        }
        // Preserve the diagnostic JSON sidecar (including legacy explicit
        // --stats *.json callers) and always publish the table that GATK
        // FilterMutectCalls derives from the driving VCF path.
        const auto stats_path = options.stats.empty() ? options.output + ".stats.json" : options.stats;
        const auto gatk_stats_path = options.output + ".stats";
        // Mutect2Engine counts the full traversal pileup before separating
        // tumor and normal reads, so a matched-normal contributes to its
        // callable statistic even though its PairHMM likelihoods are scored
        // in a separate Kokkos invocation below.
        fastgatk::io::ReadBatch callable_reads;
        callable_reads.offsets.push_back(0);
        append_batch(callable_reads, tumor_reads);
        append_batch(callable_reads, normal_reads);
        auto filtered_callable_reads =
            filtered_mutect_callable_reads(callable_reads, options.calling);
        if (!options.ignore_itr_artifacts)
            (void)apply_mutect2_itr_artifact_clipping(filtered_callable_reads, references);
        const auto callable_sites = count_mutect_callable_sites(
            filtered_callable_reads, tumor_result, options.callable_depth);
        write_gatk_mutect_stats_table(gatk_stats_path, callable_sites);
        const auto f1r2_path = options.f1r2.empty() ? options.output + ".f1r2.tsv" : options.f1r2;
        const bool standard_f1r2_archive = suffix(f1r2_path, ".tar.gz");
        std::ofstream stats(stats_path);
        if (!stats) throw std::runtime_error("cannot write Mutect2 stats: " + stats_path);
        std::uint64_t f1r2_total = 0;
        std::uint64_t r1f2_total = 0;
        stats << std::setprecision(17)
              << "{\"schema_version\":1,\"tool\":\"Mutect2\",\"status\":\"prototype\"," 
              << "\"tumor_candidate_sites\":" << tumor_result.candidate_sites
              << ",\"tumor_calls\":" << tumor_result.calls.size() << ",\"normal_calls\":"
              << (normal_result ? normal_result->calls.size() : 0) << ",\"somatic_candidates\":"
              << tumor_result.calls.size()
              << ",\"tumor_reads\":" << tumor_reads.records()
              << ",\"normal_reads\":" << normal_read_count
              << ",\"callable_sites\":" << callable_sites
              << ",\"callable_depth\":" << options.callable_depth
              << ",\"tumor_sample\":\""
              << json_escape(options.tumor_sample.empty()
                                 ? (tumor_header.samples.empty() ? std::string("TUMOR")
                                                                  : tumor_header.samples.front())
                                 : options.tumor_sample)
              << "\",\"normal_sample\":\""
              << json_escape(options.normal_sample.empty() ? std::string("NORMAL")
                                                            : options.normal_sample)
              << "\",\"normal_samples\":" << normal_sample_list_json(options)
              << ",\"sample_name_filtering\":"
              << ((!options.tumor_sample.empty() || !options.normal_sample.empty()) ? "true" : "false")
              << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
              << ",\"itr_artifact_clipping_enabled\":"
              << (options.ignore_itr_artifacts ? "false" : "true")
              << ",\"itr_artifact_clipped_reads\":"
              << tumor_result.post_read_filter_transformed_reads
              << ",\"read_filter_filtered_reads\":" << tumor_result.filtered_reads
              << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
              << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
              << ",\"effective_tumor_lod_emit\":"
              << json_number_or_null(effective_tumor_lod_emit(options))
              << ",\"normal_lod\":" << options.normal_log10_odds
              << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
              << ",\"population_allele_frequency\":" << options.population_allele_frequency
              << ",\"mitochondria_pruning_lod_threshold\":"
              << options.calling.graph_pruning_log_odds_threshold
              << ",\"use_soft_clipped_bases\":"
              << (options.calling.use_soft_clipped_bases ? "true" : "false")
              << ",\"min_base_quality\":"
              << static_cast<unsigned>(options.calling.min_base_quality)
              << ",\"pcr_indel_model\":\""
              << (options.calling.pairhmm_conservative_indel_model
                      ? options.calling.pairhmm_pcr_indel_model : "NONE")
              << "\",\"pcr_error_rate_factor\":"
              << (options.calling.pairhmm_conservative_indel_model
                      ? options.calling.pairhmm_pcr_error_rate_factor : 0.0)
              << ",\"candidate_softclip_suppressed\":" << tumor_result.candidate_softclip_suppressed
              << ",\"candidate_fragment_suppressed\":" << tumor_result.candidate_fragment_suppressed
              << ",\"multiallelic_candidates_pruned\":" << tumor_result.multiallelic_candidates_pruned
              << ",\"somatic_pileup_fallback_candidates\":"
              << tumor_result.somatic_pileup_fallback_candidates
              << ",\"candidate_summary\":[";
        for (std::size_t candidate = 0; candidate < tumor_result.candidates.size(); ++candidate) {
            if (candidate != 0) stats << ',';
            const auto& item = tumor_result.candidates[candidate];
            const auto likelihood = candidate < tumor_result.likelihoods.size()
                ? tumor_result.likelihoods[candidate] : fastgatk::calling::Likelihoods{};
            stats << "{\"tid\":" << item.tid
                  << ",\"position\":" << item.position
                  << ",\"ref\":\"" << json_escape(allele_ref(item))
                  << "\",\"alt\":\"" << json_escape(allele_alt(item))
                  << "\",\"depth\":" << item.depth
                  << ",\"alt_count\":" << item.alternate_count
                  << ",\"hom_alt_minus_ref\":"
                  << json_number_or_null(likelihood.hom_alt - likelihood.hom_ref) << '}';
        }
        stats << ']';
        std::ofstream f1r2;
        if (!standard_f1r2_archive) {
            f1r2.open(f1r2_path);
            if (!f1r2) throw std::runtime_error("cannot write Mutect2 F1R2 sidecar: " + f1r2_path);
            f1r2 << "# FASTGATK-MUTECT2-F1R2 v2\n# contig\tposition\talt\tf1r2\tr1f2\ttotal\n";
        }
        for (const auto& call : tumor_result.calls) {
            const auto orientation = orientation_counts(tumor_reads, call.candidate);
            f1r2_total += orientation.f1r2;
            r1f2_total += orientation.r1f2;
            if (!standard_f1r2_archive) {
                f1r2 << call.candidate.tid << '\t' << call.candidate.position + 1 << '\t'
                     << allele_alt(call.candidate) << '\t' << orientation.f1r2 << '\t'
                     << orientation.r1f2 << '\t' << (orientation.f1r2 + orientation.r1f2) << '\n';
            }
        }
        if (standard_f1r2_archive) {
            const auto archive_sample = make_f1r2_archive_sample(
                tumor_header, tumor_reads, tumor_result, references);
            fastgatk::somatic::write_f1r2_archive(f1r2_path, {archive_sample});
        }
        stats << ",\"f1r2_alt_support\":" << f1r2_total
              << ",\"r1f2_alt_support\":" << r1f2_total
              << ",\"multiallelic_locus_writer\":true"
              << ",\"requested_batch_records\":" << options.requested_batch_records
              << ",\"effective_batch_records\":" << options.effective_batch_records
              << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
              << ",\"pairhmm_pairs\":" << tumor_result.pairhmm_pairs
              << ",\"pairhmm_error_model\":\"" << tumor_result.pairhmm_error_model << "\""
              << ",\"pairhmm_skip_reason\":\"" << json_escape(tumor_result.pairhmm_skip_reason) << "\""
              << ",\"pairhmm_default_indel_quality\":" << static_cast<unsigned>(tumor_result.pairhmm_default_indel_quality)
              << ",\"pairhmm_default_indel_terminal_quality\":" << static_cast<unsigned>(tumor_result.pairhmm_default_indel_terminal_quality)
              << ",\"pairhmm_pcr_indel_model\":\"" << tumor_result.pairhmm_pcr_indel_model << "\""
              << ",\"pairhmm_pcr_error_rate_factor\":" << tumor_result.pairhmm_pcr_error_rate_factor
              << ",\"phred_scaled_global_read_mismapping_rate\":"
              << options.calling.phred_scaled_global_read_mismapping_rate
              << ",\"pairhmm_pcr_adjusted_positions\":" << tumor_result.pairhmm_pcr_adjusted_positions
              << ",\"pairhmm_flow_reads\":" << tumor_result.pairhmm_flow_reads
              << ",\"pairhmm_flow_haplotypes\":" << tumor_result.pairhmm_flow_haplotypes
              << ",\"pairhmm_flow_haplotypes_clipped\":" << tumor_result.pairhmm_flow_haplotypes_clipped
              << ",\"pairhmm_flow_haplotypes_collapsed\":" << tumor_result.pairhmm_flow_haplotypes_collapsed
              << ",\"pairhmm_flow_haplotypes_uncollapsed\":" << tumor_result.pairhmm_flow_haplotypes_uncollapsed
              << ",\"pairhmm_flow_haplotype_remaps\":" << tumor_result.pairhmm_flow_haplotype_remaps
              << ",\"pairhmm_flow_identical_haplotype_groups\":" << tumor_result.pairhmm_flow_identical_haplotype_groups
              << ",\"pairhmm_flow_clipping_fallbacks\":" << tumor_result.pairhmm_flow_clipping_fallbacks
              << ",\"pairhmm_insertion_quality_reads\":" << tumor_result.pairhmm_insertion_quality_reads
              << ",\"pairhmm_deletion_quality_reads\":" << tumor_result.pairhmm_deletion_quality_reads
              << ",\"pairhmm_request_links\":" << tumor_result.pairhmm_request_links
              << ",\"pairhmm_haplotypes\":" << tumor_result.pairhmm_haplotypes
              << ",\"pairhmm_graph_haplotypes\":" << tumor_result.pairhmm_graph_haplotypes
              << ",\"pairhmm_graph_haplotypes_considered\":"
              << tumor_result.pairhmm_graph_haplotypes_considered
              << ",\"pairhmm_graph_haplotypes_pruned\":"
              << tumor_result.pairhmm_graph_haplotypes_pruned
              << ",\"pairhmm_haplotypes_pruned\":"
              << tumor_result.pairhmm_haplotypes_pruned
              << ",\"pairhmm_graph_snp_posterior_pairs\":"
              << tumor_result.pairhmm_graph_snp_posterior_pairs
              << ",\"pairhmm_reads_clipped\":" << tumor_result.pairhmm_reads_clipped
              << ",\"pairhmm_reads_dropped_after_clipping\":"
              << tumor_result.pairhmm_reads_dropped_after_clipping
              << ",\"pairhmm_minimum_read_length_after_trimming\":"
              << tumor_result.pairhmm_minimum_read_length_after_trimming
              << ",\"pairhmm_haplotype_combination_blocks\":"
              << tumor_result.pairhmm_haplotype_combination_blocks
              << ",\"pairhmm_haplotype_combination_chunked\":"
              << (tumor_result.pairhmm_haplotype_combination_chunked ? "true" : "false")
              << ",\"pairhmm_assembly_region_groups\":"
              << tumor_result.pairhmm_assembly_region_groups
              << ",\"pairhmm_unassigned_candidates\":"
              << tumor_result.pairhmm_unassigned_candidates
              << ",\"pairhmm_assembly_region_partitioned\":"
              << (tumor_result.pairhmm_assembly_region_partitioned ? "true" : "false")
              << ",\"downsampled_reads\":" << tumor_result.downsampled_reads
              << ",\"overlapping_quality_correction_used\":"
              << (tumor_result.overlapping_quality_correction_used ? "true" : "false")
              << ",\"overlapping_quality_correction_metadata_available\":"
              << (tumor_result.overlapping_quality_correction_metadata_available ? "true" : "false")
              << ",\"overlapping_pairs\":" << tumor_result.overlapping_pairs
              << ",\"overlapping_bases\":" << tumor_result.overlapping_bases
              << ",\"overlapping_conflicting_bases\":"
              << tumor_result.overlapping_conflicting_bases
              << ",\"overlapping_quality_caps\":" << tumor_result.overlapping_quality_caps
              << ",\"overlapping_indel_quality_caps\":"
              << tumor_result.overlapping_indel_quality_caps
              << ",\"overlapping_pcr_snv_quality\":"
              << static_cast<unsigned>(tumor_result.overlapping_pcr_snv_quality)
              << ",\"overlapping_pcr_indel_quality\":"
              << static_cast<unsigned>(tumor_result.overlapping_pcr_indel_quality)
              << ",\"read_filter_max_reads_per_locus\":"
              << options.calling.max_reads_per_locus
              << ",\"read_filter_downsampling_seed\":"
              << options.calling.downsampling_seed
              << ",\"assembly_unassigned_candidates\":"
              << tumor_result.assembly_unassigned_candidates
              << ",\"assembly_cross_region_candidates\":"
              << tumor_result.assembly_cross_region_candidates
              << ",\"assembly_region_partitioned\":"
              << (tumor_result.assembly_region_partitioned ? "true" : "false")
              << ",\"active_probability_threshold\":"
              << options.calling.active_probability_threshold
              << ",\"force_active\":"
              << (options.force_active ? "true" : "false")
              << ",\"activity_quality_aware_somatic\":"
              << (tumor_result.activity_quality_aware_somatic ? "true" : "false")
              << ",\"activity_initial_tumor_log10_odds\":"
              << json_number_or_null(options.calling.activity_initial_tumor_log10_odds)
              << ",\"activity_pcr_snv_quality\":"
              << static_cast<unsigned>(options.calling.activity_pcr_snv_quality)
              << ",\"overlapping_pcr_indel_quality\":"
              << static_cast<unsigned>(options.calling.overlapping_pcr_indel_quality)
              << ",\"activity_multiple_substitution_quality_correction\":"
              << static_cast<unsigned>(options.calling.activity_multiple_substitution_quality_correction)
              << ",\"assembly_region_padding\":"
              << options.calling.assembly_region_padding
              << ",\"max_assembly_region_size\":"
              << options.calling.max_assembly_region_size
              << ",\"max_probability_propagation_distance\":"
              << options.calling.max_probability_propagation_distance
              << ",\"orientation_bias_model\":\"strand-counts-only\""
              << ",\"f1r2_format\":\""
              << (standard_f1r2_archive ? "collect-f1r2-counts-tar-gz" : "legacy-tsv") << "\""
              << ",\"somatic_likelihood_model\":\""
              << (tumor_somatic ? "per-read-dirichlet-variational-evidence-v1" : "genotype-summary-fallback") << "\""
              << ",\"somatic_evidence_groups\":"
              << (tumor_somatic ? tumor_somatic->evidence_groups : 0)
              << ",\"somatic_evidence_grouping\":"
              << (tumor_somatic && tumor_somatic->evidence_groups <
                          tumor_result.likelihood_read_names.size() ? "true" : "false")
              << ",\"independent_mates\":"
              << (options.independent_mates ? "true" : "false")
              << ",\"somatic_informative_read_overlap_margin\":2"
              << ",\"somatic_grouped_responsibility_depth\":true"
              << ",\"somatic_likelihood_execution_space\":\""
              << (tumor_somatic ? tumor_somatic->execution_space : "") << "\""
              << ",\"somatic_likelihood_prepare_seconds\":"
              << (tumor_somatic ? tumor_somatic->prepare_seconds : 0.0)
              << ",\"somatic_likelihood_seconds\":"
              << (tumor_somatic ? tumor_somatic->seconds : 0.0)
              << ",\"somatic_posterior_model\":\""
              << (posterior ? "somatic-germline-artifact-orientation-v1" : "unavailable") << "\""
              << ",\"somatic_posterior_execution_space\":\""
              << (posterior ? posterior->execution_space : "") << "\""
              << ",\"somatic_posterior_prepare_seconds\":"
              << (posterior ? posterior->prepare_seconds : 0.0)
              << ",\"somatic_posterior_seconds\":"
              << (posterior ? posterior->seconds : 0.0)
              << ",\"contamination_fraction\":" << options.contamination
              << ",\"sites_only_vcf_output\":"
              << (options.sites_only_vcf_output ? "true" : "false")
              << ",\"create_output_variant_index\":"
              << (options.create_output_variant_index ? "true" : "false")
              << ",\"bit_identical_to_gatk\":false}\n";
        stats.flush();
        if (f1r2.is_open()) f1r2.flush();
        if (!file_complete(options.output) || !file_complete(stats_path) ||
            !file_complete(gatk_stats_path) ||
            !file_complete(f1r2_path) || (!index.empty() && !file_complete(index)) ||
            (!options.assembly_region_out.empty() &&
             !file_complete(options.assembly_region_out)))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Mutect2 output or sidecar is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write Mutect2 manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"Mutect2\",\"implementation\":\"fastgatk-mutect2\"," 
                 << "\"status\":\"prototype\",\"primary_output\":\"" << json_escape(options.output)
                 << "\",\"primary_output_kind\":\"vcf\",\"compatibility\":{\"tumor_normal\":"
                 << (normal_result ? "true" : "false") << ",\"somatic_filter\":true,\"f1r2\":true,"
                 << "\"multiallelic_locus_writer\":true,"
                 << "\"tumor_sample_selection\":\""
                 << json_escape(options.tumor_sample.empty()
                                    ? (tumor_header.samples.empty() ? std::string("TUMOR") : tumor_header.samples.front())
                                    : options.tumor_sample)
                 << "\",\"normal_sample_selection\":\""
                 << json_escape(options.normal_sample.empty() ? "NORMAL" : options.normal_sample)
                 << "\",\"normal_samples\":" << normal_sample_list_json(options) << ","
                 << "\"sample_name_filtering\":"
                 << ((!options.tumor_sample.empty() || !options.normal_sample.empty()) ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
                 << ",\"itr_artifact_clipping_enabled\":"
                 << (options.ignore_itr_artifacts ? "false" : "true")
                 << ",\"mitochondria_defaults_applied\":" << (options.mitochondria_mode ? "true" : "false")
                 << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
                 << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
                 << ",\"effective_tumor_lod_emit\":"
                 << json_number_or_null(effective_tumor_lod_emit(options))
                 << ",\"normal_lod\":" << options.normal_log10_odds
                 << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
                 << ",\"tumor_reads\":" << tumor_reads.records()
                 << ",\"normal_reads\":" << normal_read_count
                 << ",\"orientation_counts\":true,"
                 << "\"f1r2_standard_tar\":" << (standard_f1r2_archive ? "true" : "false") << ","
                 << "\"f1r2_legacy_tsv\":" << (standard_f1r2_archive ? "false" : "true") << ","
                 << "\"somatic_read_likelihoods\":" << (tumor_somatic ? "true" : "false") << ","
                 << "\"somatic_fragment_grouping\":"
                 << (tumor_somatic && tumor_somatic->evidence_groups <
                             tumor_result.likelihood_read_names.size() ? "true" : "false") << ","
                 << "\"somatic_informative_read_overlap\":true,"
                 << "\"somatic_grouped_responsibility_depth\":true,"
                 << "\"somatic_pileup_candidate_boundary\":true,"
                 << "\"somatic_posterior\":" << (posterior ? "true" : "false") << ","
                 << "\"orientation_bias_posterior\":" << (posterior ? "true" : "false") << ","
                 << "\"contamination_posterior\":" << (posterior ? "true" : "false") << ","
                 << "\"haplotype_combination_chunking\":"
                 << (tumor_result.pairhmm_haplotype_combination_chunked ? "true" : "false") << ","
                 << "\"assembly_region_pairhmm_partitioning\":"
                 << (tumor_result.pairhmm_assembly_region_partitioned ? "true" : "false") << ","
                 << "\"assembly_region_candidate_partitioning\":"
                 << (tumor_result.assembly_region_partitioned ? "true" : "false") << ","
                 << "\"assembly_region_out\":"
                 << (!options.assembly_region_out.empty() ? "true" : "false") << ","
                 << "\"force_active\":"
                 << (options.force_active ? "true" : "false") << ","
                 << "\"activity_profile_controls\":true,"
                 << "\"activity_bandpass\":true,"
                 << "\"active_probability_threshold\":"
                 << options.calling.active_probability_threshold << ","
                 << "\"force_active\":"
                 << (options.force_active ? "true" : "false") << ","
                 << "\"activity_quality_aware_somatic\":"
                 << (tumor_result.activity_quality_aware_somatic ? "true" : "false") << ","
                 << "\"activity_initial_tumor_log10_odds\":"
                 << json_number_or_null(options.calling.activity_initial_tumor_log10_odds) << ","
                 << "\"activity_pcr_snv_quality\":"
                 << static_cast<unsigned>(options.calling.activity_pcr_snv_quality) << ","
                 << "\"overlapping_pcr_indel_quality\":"
                 << static_cast<unsigned>(options.calling.overlapping_pcr_indel_quality) << ","
                 << "\"activity_multiple_substitution_quality_correction\":"
                 << static_cast<unsigned>(options.calling.activity_multiple_substitution_quality_correction) << ","
                 << "\"assembly_region_padding\":"
                 << options.calling.assembly_region_padding << ","
                 << "\"min_assembly_region_size\":"
                 << options.calling.min_assembly_region_size << ","
                 << "\"max_assembly_region_size\":"
                 << options.calling.max_assembly_region_size << ","
                 << "\"max_probability_propagation_distance\":"
                 << options.calling.max_probability_propagation_distance << ","
                 << "\"read_haplotype_request_deduplication\":true,"
                 << "\"overlapping_base_quality_correction\":"
                 << (tumor_result.overlapping_quality_correction_used ? "true" : "false") << ","
                 << "\"independent_mates\":"
                 << (options.independent_mates ? "true" : "false") << ","
                 << "\"read_filter_gatk_defaults\":"
                 << (options.calling.minimum_mapping_quality == 20 && options.calling.exclude_duplicates &&
                     options.calling.exclude_unmapped && options.calling.exclude_secondary &&
                     !options.calling.exclude_supplementary && options.calling.exclude_qcfail &&
                     options.calling.exclude_mapping_quality_unavailable &&
                     options.calling.exclude_mapping_quality_zero &&
                     options.calling.require_good_cigar && options.calling.require_nonzero_reference_span &&
                     options.calling.require_non_chimeric_original_alignment &&
                     options.calling.require_no_n_cigar && options.calling.require_read_group &&
                     options.calling.require_read_length && options.calling.min_read_length == 30
                         ? "true" : "false") << ","
                 << "\"use_soft_clipped_bases\":"
                 << (options.calling.use_soft_clipped_bases ? "true" : "false") << ","
                 << "\"min_base_quality\":"
                 << static_cast<unsigned>(options.calling.min_base_quality) << ","
                 << "\"pcr_indel_model\":\""
                 << (options.calling.pairhmm_conservative_indel_model
                         ? options.calling.pairhmm_pcr_indel_model : "NONE")
                 << "\",\"pcr_error_rate_factor\":"
                 << (options.calling.pairhmm_conservative_indel_model
                         ? options.calling.pairhmm_pcr_error_rate_factor : 0.0) << ","
                 << "\"phred_scaled_global_read_mismapping_rate\":"
                 << options.calling.phred_scaled_global_read_mismapping_rate << ","
                 << "\"disable_tool_default_read_filters\":"
                 << (options.calling.disable_tool_default_read_filters ? "true" : "false") << ","
                 << "\"read_filter_require_read_length\":"
                 << (options.calling.require_read_length ? "true" : "false") << ","
                 << "\"read_filter_min_read_length\":" << options.calling.min_read_length << ","
                 << "\"read_filter_max_read_length\":" << options.calling.max_read_length << ","
                 << "\"deterministic_downsampling\":"
                 << (options.calling.max_reads_per_locus > 0 ? "true" : "false") << ","
                 << "\"stats\":true,\"gatk_stats_table\":true,\"callable_sites\":"
                 << callable_sites << ",\"callable_depth\":" << options.callable_depth
                 << ",\"vcf_index\":" << (index.empty() ? "false" : "true")
                 << ",\"sites_only_vcf_output\":"
                 << (options.sites_only_vcf_output ? "true" : "false")
                 << ",\"create_output_variant_index\":"
                 << (options.create_output_variant_index ? "true" : "false")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":"
                 << (file_complete(options.output) ? "true" : "false") << "},{\"path\":\""
                 << json_escape(stats_path) << "\",\"kind\":\"mutect2-stats\",\"complete\":"
                 << (file_complete(stats_path) ? "true" : "false") << "},{\"path\":\""
                 << json_escape(gatk_stats_path) << "\",\"kind\":\"mutect2-gatk-stats\",\"complete\":"
                 << (file_complete(gatk_stats_path) ? "true" : "false") << "},{\"path\":\""
                 << json_escape(f1r2_path) << "\",\"kind\":\""
                 << (standard_f1r2_archive ? "f1r2-collect-tar-gz" : "f1r2-orientation")
                 << "\",\"complete\":"
                 << (file_complete(f1r2_path) ? "true" : "false") << "}"
                 << (options.assembly_region_out.empty() ? "" :
                     ",{\"path\":\"" + json_escape(options.assembly_region_out) +
                     "\",\"kind\":\"assembly-region-igv\",\"complete\":" +
                     (file_complete(options.assembly_region_out) ? "true}" : "false}"))
                 << (index.empty() ? "" : ",{\"path\":\"" + json_escape(index) + "\",\"kind\":\"vcf-index\",\"complete\":" +
                     (file_complete(index) ? "true}" : "false}"))
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"requested_pairhmm_threads\":" << options.threads
                 << ",\"effective_pairhmm_threads\":"
                 << resources.effective_threads(static_cast<std::size_t>(options.threads))
                 << ",\"interval_list_inputs\":" << options.interval_file_inputs
                 << ",\"interval_list_records\":" << options.interval_file_records
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"mitochondria_mode\":" << (options.mitochondria_mode ? "true" : "false")
                 << ",\"itr_artifact_clipping_enabled\":"
                 << (options.ignore_itr_artifacts ? "false" : "true")
                 << ",\"itr_artifact_clipped_reads\":"
                 << tumor_result.post_read_filter_transformed_reads
                 << ",\"initial_tumor_lod\":" << options.initial_tumor_lod
                 << ",\"tumor_lod_to_emit\":" << options.tumor_lod_to_emit
                 << ",\"effective_tumor_lod_emit\":"
                 << json_number_or_null(effective_tumor_lod_emit(options))
                 << ",\"normal_lod\":" << options.normal_log10_odds
                 << ",\"minimum_allele_fraction\":" << options.minimum_allele_fraction
                 << ",\"callable_sites\":" << callable_sites
                 << ",\"callable_depth\":" << options.callable_depth
                 << ",\"population_allele_frequency\":" << options.population_allele_frequency
                 << ",\"mitochondria_pruning_lod_threshold\":"
                 << options.calling.graph_pruning_log_odds_threshold
                 << ",\"requested_batch_records\":" << options.requested_batch_records
                 << ",\"effective_batch_records\":" << options.effective_batch_records
                 << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
                 << ",\"read_filter_min_mapping_quality\":"
                 << static_cast<unsigned>(options.calling.minimum_mapping_quality)
                 << ",\"read_filter_exclude_mapping_quality_unavailable\":"
                 << (options.calling.exclude_mapping_quality_unavailable ? "true" : "false")
                 << ",\"read_filter_exclude_mapping_quality_zero\":"
                 << (options.calling.exclude_mapping_quality_zero ? "true" : "false")
                 << ",\"read_filter_exclude_duplicates\":"
                 << (options.calling.exclude_duplicates ? "true" : "false")
                 << ",\"read_filter_exclude_supplementary\":"
                 << (options.calling.exclude_supplementary ? "true" : "false")
                 << ",\"use_soft_clipped_bases\":"
                 << (options.calling.use_soft_clipped_bases ? "true" : "false")
                 << ",\"min_base_quality\":"
                 << static_cast<unsigned>(options.calling.min_base_quality)
                 << ",\"pcr_indel_model\":\""
                 << (options.calling.pairhmm_conservative_indel_model
                         ? options.calling.pairhmm_pcr_indel_model : "NONE")
                 << "\",\"pcr_error_rate_factor\":"
                 << (options.calling.pairhmm_conservative_indel_model
                         ? options.calling.pairhmm_pcr_error_rate_factor : 0.0)
                 << ",\"phred_scaled_global_read_mismapping_rate\":"
                 << options.calling.phred_scaled_global_read_mismapping_rate
                 << ",\"read_filter_require_good_cigar\":"
                 << (options.calling.require_good_cigar ? "true" : "false")
                 << ",\"read_filter_require_nonzero_reference_span\":"
                 << (options.calling.require_nonzero_reference_span ? "true" : "false")
                 << ",\"read_filter_require_non_chimeric_original_alignment\":"
                 << (options.calling.require_non_chimeric_original_alignment ? "true" : "false")
                 << ",\"read_filter_require_no_n_cigar\":"
                 << (options.calling.require_no_n_cigar ? "true" : "false")
                 << ",\"read_filter_require_read_group\":"
                 << (options.calling.require_read_group ? "true" : "false")
                 << ",\"read_filter_require_read_length\":"
                 << (options.calling.require_read_length ? "true" : "false")
                 << ",\"read_filter_min_read_length\":"
                 << options.calling.min_read_length
                 << ",\"read_filter_max_read_length\":"
                 << options.calling.max_read_length
                 << ",\"read_filter_max_reads_per_locus\":"
                 << options.calling.max_reads_per_locus
                 << ",\"read_filter_downsampling_seed\":"
                 << options.calling.downsampling_seed
                 << ",\"downsampled_reads\":" << tumor_result.downsampled_reads
                 << ",\"overlapping_quality_correction_used\":"
                 << (tumor_result.overlapping_quality_correction_used ? "true" : "false")
                 << ",\"overlapping_quality_correction_metadata_available\":"
                 << (tumor_result.overlapping_quality_correction_metadata_available ? "true" : "false")
                 << ",\"overlapping_pairs\":" << tumor_result.overlapping_pairs
                 << ",\"overlapping_bases\":" << tumor_result.overlapping_bases
                 << ",\"overlapping_conflicting_bases\":"
                 << tumor_result.overlapping_conflicting_bases
                 << ",\"overlapping_quality_caps\":" << tumor_result.overlapping_quality_caps
                 << ",\"overlapping_indel_quality_caps\":"
                 << tumor_result.overlapping_indel_quality_caps
                 << ",\"overlapping_pcr_snv_quality\":"
                 << static_cast<unsigned>(tumor_result.overlapping_pcr_snv_quality)
                 << ",\"overlapping_pcr_indel_quality\":"
                 << static_cast<unsigned>(tumor_result.overlapping_pcr_indel_quality)
                 << ",\"tumor_candidate_sites\":" << tumor_result.candidate_sites
                 << ",\"somatic_pileup_fallback_candidates\":"
                 << tumor_result.somatic_pileup_fallback_candidates
                 << ",\"tumor_reads\":" << tumor_result.reads << ",\"tumor_calls\":"
                 << tumor_result.calls.size() << ",\"somatic_likelihood_execution_space\":\""
                 << (tumor_somatic ? tumor_somatic->execution_space : "") << "\",\"somatic_likelihood_prepare_seconds\":"
                 << (tumor_somatic ? tumor_somatic->prepare_seconds : 0.0)
                 << ",\"somatic_likelihood_seconds\":"
                 << (tumor_somatic ? tumor_somatic->seconds : 0.0)
                 << ",\"somatic_evidence_groups\":"
                 << (tumor_somatic ? tumor_somatic->evidence_groups : 0)
                 << ",\"somatic_evidence_grouping\":"
                 << (tumor_somatic && tumor_somatic->evidence_groups <
                             tumor_result.likelihood_read_names.size() ? "true" : "false")
                 << ",\"somatic_informative_read_overlap_margin\":2"
                 << ",\"somatic_grouped_responsibility_depth\":true"
                 << ",\"somatic_pileup_fallback_candidates\":"
                 << tumor_result.somatic_pileup_fallback_candidates
                 << ",\"somatic_posterior_execution_space\":\""
                 << (posterior ? posterior->execution_space : "") << "\",\"somatic_posterior_prepare_seconds\":"
                 << (posterior ? posterior->prepare_seconds : 0.0)
                 << ",\"somatic_posterior_seconds\":"
                 << (posterior ? posterior->seconds : 0.0)
                 << ",\"somatic_posterior_model\":\""
                 << (posterior ? "somatic-germline-artifact-orientation-v1" : "unavailable")
                 << "\",\"contamination_fraction\":" << options.contamination
                 << ",\"pairhmm_pairs\":" << tumor_result.pairhmm_pairs
                 << ",\"pairhmm_error_model\":\"" << tumor_result.pairhmm_error_model << "\""
                 << ",\"pairhmm_skip_reason\":\"" << json_escape(tumor_result.pairhmm_skip_reason) << "\""
                 << ",\"pairhmm_default_indel_quality\":" << static_cast<unsigned>(tumor_result.pairhmm_default_indel_quality)
                 << ",\"pairhmm_default_indel_terminal_quality\":" << static_cast<unsigned>(tumor_result.pairhmm_default_indel_terminal_quality)
                 << ",\"pairhmm_pcr_indel_model\":\"" << tumor_result.pairhmm_pcr_indel_model << "\""
                 << ",\"pairhmm_pcr_error_rate_factor\":" << tumor_result.pairhmm_pcr_error_rate_factor
                 << ",\"pairhmm_pcr_adjusted_positions\":" << tumor_result.pairhmm_pcr_adjusted_positions
                 << ",\"pairhmm_flow_reads\":" << tumor_result.pairhmm_flow_reads
                 << ",\"pairhmm_flow_haplotypes\":" << tumor_result.pairhmm_flow_haplotypes
                 << ",\"pairhmm_flow_haplotypes_clipped\":" << tumor_result.pairhmm_flow_haplotypes_clipped
                 << ",\"pairhmm_flow_haplotypes_collapsed\":" << tumor_result.pairhmm_flow_haplotypes_collapsed
                 << ",\"pairhmm_flow_haplotypes_uncollapsed\":" << tumor_result.pairhmm_flow_haplotypes_uncollapsed
                 << ",\"pairhmm_flow_haplotype_remaps\":" << tumor_result.pairhmm_flow_haplotype_remaps
                 << ",\"pairhmm_flow_identical_haplotype_groups\":" << tumor_result.pairhmm_flow_identical_haplotype_groups
                 << ",\"pairhmm_flow_clipping_fallbacks\":" << tumor_result.pairhmm_flow_clipping_fallbacks
                 << ",\"pairhmm_insertion_quality_reads\":" << tumor_result.pairhmm_insertion_quality_reads
                 << ",\"pairhmm_deletion_quality_reads\":" << tumor_result.pairhmm_deletion_quality_reads
                 << ",\"pairhmm_request_links\":" << tumor_result.pairhmm_request_links
                 << ",\"pairhmm_haplotypes\":" << tumor_result.pairhmm_haplotypes
                 << ",\"pairhmm_graph_haplotypes\":" << tumor_result.pairhmm_graph_haplotypes
                 << ",\"pairhmm_graph_haplotypes_considered\":"
                 << tumor_result.pairhmm_graph_haplotypes_considered
                 << ",\"pairhmm_graph_haplotypes_pruned\":"
                 << tumor_result.pairhmm_graph_haplotypes_pruned
                 << ",\"pairhmm_haplotypes_pruned\":"
                 << tumor_result.pairhmm_haplotypes_pruned
                 << ",\"pairhmm_graph_snp_posterior_pairs\":"
                 << tumor_result.pairhmm_graph_snp_posterior_pairs
                 << ",\"pairhmm_haplotype_combination_blocks\":"
                 << tumor_result.pairhmm_haplotype_combination_blocks
                 << ",\"pairhmm_haplotype_combination_chunked\":"
                 << (tumor_result.pairhmm_haplotype_combination_chunked ? "true" : "false")
                 << ",\"pairhmm_assembly_region_groups\":"
                 << tumor_result.pairhmm_assembly_region_groups
                 << ",\"pairhmm_unassigned_candidates\":"
                 << tumor_result.pairhmm_unassigned_candidates
                 << ",\"pairhmm_assembly_region_partitioned\":"
                 << (tumor_result.pairhmm_assembly_region_partitioned ? "true" : "false")
                 << ",\"assembly_unassigned_candidates\":"
                 << tumor_result.assembly_unassigned_candidates
                 << ",\"assembly_cross_region_candidates\":"
                 << tumor_result.assembly_cross_region_candidates
                 << ",\"assembly_region_partitioned\":"
                 << (tumor_result.assembly_region_partitioned ? "true" : "false")
                 << ",\"active_probability_threshold\":"
                 << options.calling.active_probability_threshold
                 << ",\"force_active\":"
                 << (options.force_active ? "true" : "false")
                 << ",\"activity_quality_aware_somatic\":"
                 << (tumor_result.activity_quality_aware_somatic ? "true" : "false")
                 << ",\"activity_initial_tumor_log10_odds\":"
                 << json_number_or_null(options.calling.activity_initial_tumor_log10_odds)
                 << ",\"activity_pcr_snv_quality\":"
                 << static_cast<unsigned>(options.calling.activity_pcr_snv_quality)
                 << ",\"overlapping_pcr_indel_quality\":"
                 << static_cast<unsigned>(options.calling.overlapping_pcr_indel_quality)
                 << ",\"activity_multiple_substitution_quality_correction\":"
                 << static_cast<unsigned>(options.calling.activity_multiple_substitution_quality_correction)
                 << ",\"assembly_region_padding\":"
                 << options.calling.assembly_region_padding
                 << ",\"min_assembly_region_size\":"
                 << options.calling.min_assembly_region_size
                 << ",\"max_assembly_region_size\":"
                 << options.calling.max_assembly_region_size
                 << ",\"max_probability_propagation_distance\":"
                 << options.calling.max_probability_propagation_distance
                 << ",\"graph_kmer_size\":" << options.calling.graph_kmer_size
                 << ",\"graph_kmer_sizes\":[";
        for (std::size_t index = 0; index < options.calling.graph_kmer_sizes.size(); ++index) {
            if (index != 0) manifest << ',';
            manifest << options.calling.graph_kmer_sizes[index];
        }
        manifest << ']'
                 << ",\"graph_min_kmer_count\":" << options.calling.graph_min_kmer_count
                 << ",\"graph_min_pruning\":" << options.calling.graph_min_pruning
                 << ",\"graph_num_pruning_samples\":" << options.calling.graph_num_pruning_samples
                 << ",\"graph_use_adaptive_pruning\":"
                 << (options.calling.graph_use_adaptive_pruning ? "true" : "false")
                 << ",\"graph_initial_error_rate_for_pruning\":"
                 << options.calling.graph_initial_error_rate_for_pruning
                 << ",\"graph_pruning_log_odds_threshold\":"
                 << options.calling.graph_pruning_log_odds_threshold
                 << ",\"graph_pruning_seeding_log_odds_threshold\":"
                 << options.calling.graph_pruning_seeding_log_odds_threshold
                 << ",\"graph_max_unpruned_variants\":"
                 << options.calling.graph_max_unpruned_variants
                 << ",\"graph_adaptive_pruned_nodes\":"
                 << tumor_result.graph_adaptive_pruned_nodes
                 << ",\"graph_linked_de_bruijn\":"
                 << (options.calling.graph_linked_de_bruijn ? "true" : "false")
                 << ",\"graph_disable_artificial_haplotype_recovery\":"
                 << (options.calling.graph_disable_artificial_haplotype_recovery ? "true" : "false")
                 << ",\"graph_enable_legacy_cycle_detection\":"
                 << (options.calling.graph_enable_legacy_cycle_detection ? "true" : "false")
                 << ",\"graph_min_dangling_branch_length\":"
                 << options.calling.graph_min_dangling_branch_length
                 << ",\"graph_min_dangling_matching_bases\":"
                 << options.calling.graph_min_dangling_matching_bases
                 << ",\"graph_allow_non_unique_kmers_in_ref\":"
                 << (options.calling.graph_allow_non_unique_kmers_in_ref ? "true" : "false")
                 << ",\"graph_recover_all_dangling_branches\":"
                 << (options.calling.graph_recover_all_dangling_branches ? "true" : "false")
                 << ",\"graph_dont_increase_kmer_sizes_for_cycles\":"
                 << (options.calling.graph_dont_increase_kmer_sizes_for_cycles ? "true" : "false")
                 << ",\"graph_kmer_size_selected\":"
                 << tumor_result.graph_kmer_size_selected
                 << ",\"graph_kmer_iterations\":"
                 << tumor_result.graph_kmer_iterations
                 << ",\"graph_has_non_reference_cycles\":"
                 << (tumor_result.graph_has_non_reference_cycles ? "true" : "false")
                 << ",\"error_correct_reads\":"
                 << (options.calling.error_correct_reads ? "true" : "false")
                 << ",\"error_correction_kmer_length\":"
                 << options.calling.error_correction_kmer_length
                 << ",\"error_correction_min_solid_observations\":"
                 << options.calling.error_correction_min_solid_observations
                 << ",\"error_correction_max_mismatches\":2"
                 << ",\"error_correction_max_sparse_observations\":1"
                 << ",\"error_correction_pileup_log_odds\":"
                 << (std::isfinite(options.calling.pileup_error_correction_log_odds)
                     ? std::to_string(options.calling.pileup_error_correction_log_odds) : "null")
                 << ",\"error_correction_corrected_reads\":"
                 << tumor_result.error_correction_corrected_reads
                 << ",\"error_correction_corrected_bases\":"
                 << tumor_result.error_correction_corrected_bases
                 << ",\"error_correction_solid_kmers\":"
                 << tumor_result.error_correction_solid_kmers
                 << ",\"error_correction_corrected_kmers\":"
                 << tumor_result.error_correction_corrected_kmers
                 << ",\"error_correction_uncorrectable_kmers\":"
                 << tumor_result.error_correction_uncorrectable_kmers
                 << ",\"error_correction_pileup_loci\":"
                 << tumor_result.error_correction_pileup_loci
                 << ",\"error_correction_pileup_skipped_indel_adjacent_bases\":"
                 << tumor_result.error_correction_pileup_skipped_indel_adjacent_bases
                 << ",\"error_correction_mode\":\""
                 << tumor_result.error_correction_mode << "\""
                 << ",\"error_correction_execution_space\":\""
                 << tumor_result.error_correction_execution_space << "\""
                 << ",\"error_correction_prepare_seconds\":"
                 << tumor_result.error_correction_prepare_seconds
                 << ",\"error_correction_seconds\":"
                 << tumor_result.error_correction_seconds
                 << ",\"graph_dangling_branch_paths\":"
                 << tumor_result.graph_dangling_branch_paths
                 << ",\"graph_dangling_branch_bases\":"
                 << tumor_result.graph_dangling_branch_bases
                 << ",\"graph_dangling_recovered_paths\":"
                 << tumor_result.graph_dangling_recovered_paths
                 << ",\"graph_dangling_recovered_bases\":"
                 << tumor_result.graph_dangling_recovered_bases
                 << ",\"graph_artificial_haplotype_recovery_paths\":"
                 << tumor_result.graph_artificial_haplotype_recovery_paths
                 << ",\"graph_artificial_haplotype_recovery_bases\":"
                 << tumor_result.graph_artificial_haplotype_recovery_bases
                 << ",\"graph_seqgraph_nodes\":"
                 << tumor_result.graph_seqgraph_nodes
                 << ",\"graph_seqgraph_edges\":"
                 << tumor_result.graph_seqgraph_edges
                 << ",\"graph_seqgraph_linear_chain_merges\":"
                 << tumor_result.graph_seqgraph_linear_chain_merges
                 << ",\"graph_seqgraph_diamond_merges\":"
                 << tumor_result.graph_seqgraph_diamond_merges
                 << ",\"graph_seqgraph_tail_merges\":"
                 << tumor_result.graph_seqgraph_tail_merges
                 << ",\"graph_seqgraph_suffix_splits\":"
                 << tumor_result.graph_seqgraph_suffix_splits
                 << ",\"graph_seqgraph_suffix_merges\":"
                 << tumor_result.graph_seqgraph_suffix_merges
                 << ",\"graph_non_unique_kmers\":"
                 << tumor_result.graph_non_unique_kmers
                 << ",\"graph_reference_non_unique_kmers\":"
                 << tumor_result.graph_reference_non_unique_kmers
                 << ",\"graph_reference_kmer_rejected\":"
                 << (tumor_result.graph_reference_kmer_rejected ? "true" : "false")
                 << ",\"graph_max_paths\":" << options.calling.graph_max_paths
                 << ",\"graph_max_depth\":" << options.calling.graph_max_depth
                 << ",\"flow_assembly_collapse_hmer_size\":"
                 << options.calling.flow_assembly_collapse_hmer_size
                 << ",\"flow_assembly_collapse_partial_mode\":"
                 << (options.calling.flow_assembly_collapse_partial_mode ? "true" : "false")
                 << ",\"max_haplotype_combination_alleles\":"
                 << options.calling.max_haplotype_combination_alleles
                 << ",\"max_mnp_distance\":" << options.calling.max_mnp_distance
                 << ",\"graph_mnp_candidates\":" << tumor_result.graph_mnp_candidates << "}}\n";
        std::cout << "{\"tool\":\"Mutect2\",\"status\":\"prototype\",\"tumor_reads\":"
                  << tumor_result.reads << ",\"tumor_calls\":" << tumor_result.calls.size()
                  << ",\"normal_reads\":" << normal_read_count
                  << ",\"normal_calls\":" << (normal_result ? normal_result->calls.size() : 0)
                  << ",\"tumor_sample\":\""
                  << json_escape(options.tumor_sample.empty()
                                     ? (tumor_header.samples.empty() ? std::string("TUMOR")
                                                                      : tumor_header.samples.front())
                                     : options.tumor_sample)
                  << "\",\"normal_sample\":\""
                  << json_escape(options.normal_sample.empty() ? std::string("NORMAL")
                                                                : options.normal_sample)
                  << "\",\"normal_samples\":" << normal_sample_list_json(options)
                  << ",\"sample_name_filtering\":"
                  << ((!options.tumor_sample.empty() || !options.normal_sample.empty()) ? "true" : "false")
                  << "}\n";
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        if (initialized) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}

#endif
