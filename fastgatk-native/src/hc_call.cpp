#include "fastgatk/calling/pipeline.hpp"
#include "fastgatk/kernels/genotype.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/output.hpp"
#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif
#if FASTGATK_HAS_HTSLIB
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    // GATK accepts repeatable -I/--input values (typically multiple shards for
    // one sample). Keep `input` as the first-input compatibility field used by
    // existing telemetry/rendering code, while `inputs` carries the complete
    // ordered set through the Host aggregation path.
    std::string input;
    std::vector<std::string> inputs;
    std::string reference;
    // HaplotypeCaller's GenotypeGivenAlleles feature input.  The VCF is
    // parsed on Host and supplies eligible concrete events to the existing
    // ActivityProfile/assembly/PairHMM pipeline.
    std::string alleles;
    bool force_call_filtered_alleles = false;
    std::string sample_name;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> exclusions;
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    std::string output;
    // AssemblyRegionWalker debug contract.  Like GATK's inherited
    // --assembly-region-out surface, this records the raw active/inactive
    // profile segments; assembly and PairHMM still use padded active windows.
    std::string assembly_region_out;
    std::string manifest;
    std::string telemetry;
    std::string index_path;
    bool create_output_variant_index = true;
    bool add_output_vcf_command_line = true;
    // GATK's writer-only sites-only mode still computes genotype/annotation
    // evidence before dropping FORMAT/sample columns at serialization.
    bool sites_only_vcf_output = false;
    // GATK's --floor-blocks is a reference-block writer mode: it retains
    // GT/DP/GQ but removes block-only MIN_DP/PL fields after genotyping.
    bool floor_blocks = false;
    bool native_pair_hmm_use_double_precision = false;
    bool native_pair_hmm_use_double_precision_explicit = false;
    bool disable_sequence_dictionary_validation = false;
    bool gvcf = false;
    std::string gvcf_mode;
    bool gvcf_mode_explicit = false;
    bool min_depth_explicit = false;
    bool min_alt_support_explicit = false;
    // GATK's --sample-ploidy controls Number=G genotype ordering and the
    // allele-count denominator.  Both the reference-backed VCF writer and
    // the reference-confidence block model consume the shared bounded
    // arbitrary-ploidy posterior path (the CLI keeps the production VCF
    // width bounded to eight copies to avoid unbounded Number=G rows).
    int sample_ploidy = 2;
    // GenotypeCalculationArgumentCollection.maxGenotypeCount. This is an
    // upstream HaplotypeCaller EventMap/haplotype selection limit, distinct
    // from the later GL-based --max-alternate-alleles writer subsetting.
    std::size_t max_genotype_count = 1024;
    // GenotypingEngine applies this limit after a reference-confidence
    // candidate has gained its symbolic <NON_REF> allele, but before the AF
    // posterior.  Keep it at the Host/output boundary: selecting an allele
    // subset and permuting ragged Number=G/Number=R fields is structural
    // work, while the resulting likelihood/AF calculations remain Kokkos.
    std::size_t max_alternate_alleles = 6;
    std::size_t batch_records = 4096;
    std::size_t requested_batch_records = 4096;
    std::size_t effective_batch_records = 0;
    std::uint64_t adaptive_batch_reductions = 0;
    // Opt-in coordinate-sorted contig streaming.  The default keeps the
    // historical aggregate HostBatch for strict compatibility; streaming
    // flushes one contig before decoding the next and never retains the full
    // input read payload.
    bool stream_by_contig = false;
    bool stream_by_region = false;
    std::size_t stream_region_size = 1U << 20;
    std::size_t streamed_contigs = 0;
    std::size_t streamed_regions = 0;
    std::uint64_t streamed_region_splits = 0;
    std::uint64_t streamed_peak_host_bytes = 0;
    bool stream_indexed = false;
    // Runtime ThreeStagePipeline telemetry for the bounded contig path.
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
    int threads = 1;
    fastgatk::calling::Options calling;
    // Keep the CLI's GATK default list (10,25) separate from the shared
    // scalar graph_kmer_size ABI field. The first explicit --kmer-size clears
    // this default and subsequent occurrences remain repeatable.
    bool graph_kmer_sizes_explicit = false;
    std::vector<std::string> read_filters;
    std::vector<std::string> disabled_read_filters;
    std::vector<std::string> compatibility_options;
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

// Keep the native filter surface deliberately explicit.  These are the
// built-in filters that map to fields in the shared Kokkos ordinal mask; a
// class we cannot represent must fail closed so the dispatcher can route the
// original argv to the Java fallback instead of silently ignoring it.
void enable_read_filter(fastgatk::calling::Options& calling,
                        const std::string& filter_name) {
    if (filter_name == "MappingQualityReadFilter") {
        // The filter uses the same --minimum-mapping-quality parameter as the
        // GATK tool-default instance, so the explicit class is already active.
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
    // Unlike the old implementation, disabling tool defaults does not disable
    // explicitly requested filters or the shared mask itself.
    options.calling.apply_read_filters = true;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto ch : text) {
        switch (ch) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default: escaped << ch; break;
        }
    }
    return escaped.str();
}

std::string candidate_reference(const fastgatk::calling::AssemblyCandidate& candidate) {
    return candidate.reference_allele.empty()
        ? std::string(1, static_cast<char>(candidate.reference)) : candidate.reference_allele;
}

std::string candidate_alternate(const fastgatk::calling::AssemblyCandidate& candidate) {
    return candidate.alternate_allele.empty()
        ? std::string(1, static_cast<char>(candidate.alternate)) : candidate.alternate_allele;
}

std::string output_sample_name(const fastgatk::io::HeaderSummary& header) {
    return header.samples.empty() ? "FASTGATK" : header.samples.front();
}

bool is_bp_resolution_mode(const std::string& mode) {
    std::string normalized;
    normalized.reserve(mode.size());
    for (const auto character : mode)
        normalized.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(character))));
    return normalized == "BP_RESOLUTION";
}

std::string normalize_ref_confidence_mode(const std::string& mode) {
    std::string normalized;
    normalized.reserve(mode.size());
    for (const auto character : mode)
        normalized.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(character))));
    return normalized;
}

std::array<int, 3> phred_likelihoods(const fastgatk::calling::Likelihoods& likelihoods) {
    const double maximum = std::max({likelihoods.hom_ref, likelihoods.het,
                                     likelihoods.hom_alt});
    const std::array<double, 3> values{
        likelihoods.hom_ref, likelihoods.het, likelihoods.hom_alt};
    std::array<int, 3> result{};
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double raw = std::isfinite(values[i]) && std::isfinite(maximum)
            ? -10.0 * (values[i] - maximum) : 999.0;
        result[i] = static_cast<int>(std::clamp<long long>(
            std::llround(std::max(0.0, raw)), 0LL,
            static_cast<long long>(std::numeric_limits<int>::max())));
    }
    return result;
}

std::size_t genotype_width(int allele_count, int ploidy) {
    if (allele_count < 1 || ploidy < 0) return 0;
    std::size_t width = 1;
    for (int copy = 1; copy <= ploidy; ++copy) {
        const auto numerator = static_cast<std::size_t>(allele_count + copy - 1);
        if (width > std::numeric_limits<std::size_t>::max() / numerator) return 0;
        width = (width * numerator) / static_cast<std::size_t>(copy);
    }
    return width;
}

// VCF Number=G ordering is the colexicographic order of non-decreasing allele
// vectors (for diploid this is 0/0, 0/1, 1/1, 0/2, ...): the rightmost copy
// is the major index, followed by the next copy.  Keeping the enumerator here,
// rather than relying on a vendor genotype helper, gives the same ordering to
// scalar/OpenMP/CUDA Kokkos builds and to the text writer.
std::vector<std::vector<int>> enumerate_genotypes(int allele_count, int ploidy) {
    std::vector<std::vector<int>> result;
    const auto width = genotype_width(allele_count, ploidy);
    if (width == 0 || width > 1'000'000) return result;
    result.reserve(width);
    std::vector<int> current(static_cast<std::size_t>(ploidy), 0);
    std::function<void(int, int)> visit = [&](const int position, const int maximum) {
        if (position < 0) {
            result.push_back(current);
            return;
        }
        for (int allele = 0; allele <= maximum; ++allele) {
            current[static_cast<std::size_t>(position)] = allele;
            visit(position - 1, allele);
        }
    };
    visit(ploidy - 1, allele_count - 1);
    return result;
}

// GenotypingEngine delegates --max-alternate-alleles selection to
// AlleleSubsettingUtils.calculateMostLikelyAlleles().  For every sample that
// has GLs, GATK finds the most likely *complete* genotype (REF and
// <NON_REF> are both eligible), adds its distance from hom-ref to each proper
// ALT present in that genotype, then retains the highest-scoring proper ALTs.
// Java's ordered stream sort is stable, so equal scores retain lower original
// allele indexes.  The caller is single-sample today, but this helper keeps
// the score representation additive so a future multi-sample Host can sum
// another sample's GL contribution without changing the rule.
std::vector<std::size_t> gatk_most_likely_alt_indices(
    const std::vector<int>& pl,
    int allele_count,
    int ploidy,
    std::size_t proper_alt_count,
    std::size_t max_alternate_alleles) {
    if (proper_alt_count == 0 || max_alternate_alleles >= proper_alt_count)
        return [&]() {
            std::vector<std::size_t> all(proper_alt_count);
            std::iota(all.begin(), all.end(), 0U);
            return all;
        }();
    const auto genotypes = enumerate_genotypes(allele_count, ploidy);
    if (genotypes.empty() || pl.size() != genotypes.size())
        throw std::runtime_error(
            "NUMERICAL_CONTRACT_FAILURE: incomplete GL matrix for --max-alternate-alleles");
    const auto best = static_cast<std::size_t>(std::distance(
        pl.begin(), std::min_element(pl.begin(), pl.end())));
    const auto reference = pl.front();
    const auto score = static_cast<std::uint64_t>(std::llabs(
        static_cast<long long>(pl[best]) - static_cast<long long>(reference)));
    std::vector<std::uint64_t> likelihood_sums(proper_alt_count, 0U);
    int previous_allele = -1;
    for (const auto allele : genotypes[best]) {
        // GenotypeAlleleCounts.containsAllele() contributes a supporting
        // genotype once per ALT, not once per chromosome copy.  A 0/1/2/2
        // call therefore ties ALT-1 and ALT-2 here; multiplying ALT-2 by
        // its copy count would incorrectly reverse GATK's stable tie-break.
        if (allele == previous_allele) continue;
        previous_allele = allele;
        if (allele <= 0 || static_cast<std::size_t>(allele) > proper_alt_count)
            continue;
        auto& sum = likelihood_sums[static_cast<std::size_t>(allele - 1)];
        sum = std::min<std::uint64_t>(
            std::numeric_limits<std::uint64_t>::max() - sum < score
                ? std::numeric_limits<std::uint64_t>::max() : sum + score,
            std::numeric_limits<std::uint64_t>::max());
    }
    std::vector<std::size_t> ranked(proper_alt_count);
    std::iota(ranked.begin(), ranked.end(), 0U);
    std::stable_sort(ranked.begin(), ranked.end(), [&](const auto left, const auto right) {
        return likelihood_sums[left] > likelihood_sums[right];
    });
    ranked.resize(max_alternate_alleles);
    std::sort(ranked.begin(), ranked.end());
    return ranked;
}

// Recode a Number=G PL array exactly as AlleleSubsettingUtils.subsetAlleles:
// copy the likelihood of the old genotype represented by each new genotype,
// then normalize the subset so its best likelihood is zero.  `new_to_old`
// includes REF at element zero and may retain symbolic `*`/`<NON_REF>` at a
// later element.
std::vector<int> remap_pl_for_allele_subset(const std::vector<int>& pl,
                                            int old_allele_count,
                                            int ploidy,
                                            const std::vector<std::size_t>& new_to_old) {
    const auto old_genotypes = enumerate_genotypes(old_allele_count, ploidy);
    const auto new_genotypes = enumerate_genotypes(
        static_cast<int>(new_to_old.size()), ploidy);
    if (old_genotypes.empty() || new_genotypes.empty() || pl.size() != old_genotypes.size())
        throw std::runtime_error(
            "NUMERICAL_CONTRACT_FAILURE: cannot remap incomplete Number=G likelihoods");
    std::map<std::vector<int>, std::size_t> old_indexes;
    for (std::size_t index = 0; index < old_genotypes.size(); ++index)
        old_indexes.emplace(old_genotypes[index], index);
    std::vector<int> remapped(new_genotypes.size(), 999);
    for (std::size_t index = 0; index < new_genotypes.size(); ++index) {
        std::vector<int> old_genotype;
        old_genotype.reserve(new_genotypes[index].size());
        for (const auto new_allele : new_genotypes[index]) {
            if (new_allele < 0 || static_cast<std::size_t>(new_allele) >= new_to_old.size())
                throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: invalid allele recoding");
            old_genotype.push_back(static_cast<int>(new_to_old[static_cast<std::size_t>(new_allele)]));
        }
        std::sort(old_genotype.begin(), old_genotype.end());
        const auto old = old_indexes.find(old_genotype);
        if (old == old_indexes.end())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: missing old genotype recoding");
        remapped[index] = pl[old->second];
    }
    const auto minimum = *std::min_element(remapped.begin(), remapped.end());
    for (auto& value : remapped)
        value = std::clamp(value - minimum, 0, std::numeric_limits<int>::max());
    return remapped;
}

template <typename Value>
std::vector<Value> remap_r_length_values(const std::vector<Value>& values,
                                          const std::vector<std::size_t>& new_to_old) {
    std::vector<Value> remapped;
    remapped.reserve(new_to_old.size());
    for (const auto old : new_to_old) {
        if (old >= values.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: incomplete Number=R field");
        remapped.push_back(values[old]);
    }
    return remapped;
}

std::size_t candidate_index_for_candidate(
    const fastgatk::calling::Result& result,
    const fastgatk::calling::AssemblyCandidate& candidate) {
    const auto ref = candidate_reference(candidate);
    const auto alt = candidate_alternate(candidate);
    for (std::size_t index = 0; index < result.candidates.size(); ++index) {
        const auto& item = result.candidates[index];
        if (item.tid == candidate.tid && item.position == candidate.position &&
            candidate_reference(item) == ref && candidate_alternate(item) == alt)
            return index;
    }
    return std::numeric_limits<std::size_t>::max();
}

// FORMAT/AD for an emitted allele set is an allele partition, not one
// independent biallelic count per ALT.  GATK's AlleleLikelihoods first
// chooses the best allele for each read and only keeps a read when the best
// versus second-best margin is informative.  The PairHMM matrix already
// contains the concrete REF/ALT likelihoods in read-major coordinates, so
// keep this reduction at the Host/output boundary while sharing the same
// likelihood source as PL/GT.  This prevents a read that supports ALT-1 from
// being counted as REF for ALT-2 (the old path inflated REF AD on tetra-allelic
// sites).
std::optional<std::pair<std::vector<int>, int>> derive_multiallelic_depths(
    const fastgatk::calling::Result& result,
    const std::vector<const fastgatk::calling::GenotypeCall*>& calls,
    bool include_spanning_deletion = false) {
    if (calls.empty() || result.allele_read_likelihoods.empty() ||
        result.reference_read_likelihoods.empty())
        return std::nullopt;
    std::vector<std::size_t> indices;
    indices.reserve(calls.size());
    std::size_t read_count = std::numeric_limits<std::size_t>::max();
    for (const auto* call : calls) {
        if (call == nullptr) return std::nullopt;
        const auto index = candidate_index_for_candidate(result, call->candidate);
        if (index == std::numeric_limits<std::size_t>::max() ||
            index >= result.allele_read_likelihoods.size() ||
            index >= result.reference_read_likelihoods.size())
            return std::nullopt;
        indices.push_back(index);
        read_count = std::min(read_count, result.allele_read_likelihoods[index].size());
        read_count = std::min(read_count, result.reference_read_likelihoods[index].size());
    }
    if (read_count == 0 || read_count == std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    constexpr double kInformativeMargin = 0.2;
    const auto spanning_index = calls.size() + 1;
    std::vector<int> depths(calls.size() + 1 + (include_spanning_deletion ? 1U : 0U), 0);
    int depth = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        if (read < result.annotation_read_qualified.size() &&
            result.annotation_read_qualified[read] == 0)
            continue;
        std::vector<double> scores(calls.size() + 1 + (include_spanning_deletion ? 1U : 0U),
                                   -std::numeric_limits<double>::infinity());
        for (const auto index : indices)
            scores[0] = std::max(scores[0], result.reference_read_likelihoods[index][read]);
        for (std::size_t alt = 0; alt < indices.size(); ++alt)
            scores[alt + 1] = result.allele_read_likelihoods[indices[alt]][read];
        if (include_spanning_deletion) {
            double spanning = -std::numeric_limits<double>::infinity();
            for (const auto index : indices) {
                if (index < result.spanning_deletion_read_likelihoods.size() &&
                    read < result.spanning_deletion_read_likelihoods[index].size()) {
                    spanning = std::max(spanning,
                        result.spanning_deletion_read_likelihoods[index][read]);
                }
            }
            scores[spanning_index] = spanning;
        }
        std::size_t best = 0;
        double best_score = scores[0];
        double second_score = -std::numeric_limits<double>::infinity();
        for (std::size_t allele = 1; allele < scores.size(); ++allele) {
            if (scores[allele] > best_score) {
                second_score = best_score;
                best_score = scores[allele];
                best = allele;
            } else {
                second_score = std::max(second_score, scores[allele]);
            }
        }
        if (!std::isfinite(best_score) || !std::isfinite(second_score) ||
            best_score - second_score <= kInformativeMargin)
            continue;
        ++depth;
        ++depths[best];
    }
    if (depth == 0) return std::nullopt;
    return std::make_pair(std::move(depths), depth);
}

// Build the concrete diploid PL matrix directly from PairHMM's per-read
// allele likelihoods.  This is the same likelihood model used by GATK's
// GenotypeLikelihoodCalculator for ploidy two; unlike the old writer path it
// evaluates ALT_i/ALT_j genotypes jointly rather than taking a max of two
// independently normalized homozygous penalties.
bool joint_candidate_pl(
    const fastgatk::calling::Result& result,
    const std::vector<const fastgatk::calling::AssemblyCandidate*>& candidates,
    int ploidy,
    std::vector<int>& output,
    int max_phred = 0,
    bool include_spanning_deletion = false) {
    if (candidates.empty() || result.allele_read_likelihoods.empty() ||
        result.reference_read_likelihoods.empty() || ploidy < 1) return false;
    const std::size_t allele_count = candidates.size() + 1U +
        (include_spanning_deletion ? 1U : 0U);
    const std::size_t spanning_index = candidates.size() + 1U;
    const auto genotypes = enumerate_genotypes(static_cast<int>(allele_count), ploidy);
    if (genotypes.empty()) return false;
    std::vector<std::size_t> candidate_indices;
    candidate_indices.reserve(candidates.size());
    std::size_t read_count = std::numeric_limits<std::size_t>::max();
    for (const auto* candidate : candidates) {
        if (candidate == nullptr) return false;
        const auto index = candidate_index_for_candidate(result, *candidate);
        if (index == std::numeric_limits<std::size_t>::max() ||
            index >= result.allele_read_likelihoods.size() ||
            index >= result.reference_read_likelihoods.size()) return false;
        candidate_indices.push_back(index);
        read_count = std::min(read_count, result.reference_read_likelihoods[index].size());
        read_count = std::min(read_count, result.allele_read_likelihoods[index].size());
    }
    if (read_count == 0 || read_count == std::numeric_limits<std::size_t>::max()) return false;
    std::vector<double> likelihoods(allele_count * read_count,
        -std::numeric_limits<double>::infinity());
    std::size_t informative_reads = 0;
    std::size_t spanning_reads = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        double ref_likelihood = -std::numeric_limits<double>::infinity();
        for (const auto index : candidate_indices) {
            if (std::isfinite(result.reference_read_likelihoods[index][read])) {
                ref_likelihood = result.reference_read_likelihoods[index][read];
                break;
            }
        }
        if (!std::isfinite(ref_likelihood)) continue;
        likelihoods[read] = ref_likelihood;
        bool has_alt = false;
        for (std::size_t alt = 0; alt < candidate_indices.size(); ++alt) {
            const auto value = result.allele_read_likelihoods[candidate_indices[alt]][read];
            likelihoods[(alt + 1) * read_count + read] = value;
            has_alt = has_alt || std::isfinite(value);
        }
        if (include_spanning_deletion) {
            double spanning = -std::numeric_limits<double>::infinity();
            for (const auto index : candidate_indices) {
                if (index < result.spanning_deletion_read_likelihoods.size() &&
                    read < result.spanning_deletion_read_likelihoods[index].size()) {
                    spanning = std::max(spanning,
                        result.spanning_deletion_read_likelihoods[index][read]);
                }
            }
            likelihoods[spanning_index * read_count + read] = spanning;
            if (std::isfinite(spanning)) ++spanning_reads;
        }
        if (has_alt) ++informative_reads;
    }
    if (informative_reads == 0 || (include_spanning_deletion && spanning_reads == 0))
        return false;
    const auto joint = fastgatk::kernels::calculate_joint_genotype_pl_kokkos(
        likelihoods, read_count, static_cast<int>(allele_count), ploidy, max_phred);
    output.assign(joint.pl.begin(), joint.pl.end());
    return true;
}

// Materialize GenotypePriorCalculator.assumingHW for an ALT group.
// `Result::candidate_prior_het/hom_alt` already contains the per-allele
// log10 values (including the SNP normalization constant).  When requested,
// the final symbolic <NON_REF> allele is appended as GATK's OTHER type, whose
// heterozygosity is max(snpHet, indelHet) and is not SNP-normalized.  The
// genotype allele-count rule is shared for all ploidies:
//   count=1 -> het, count>1 -> het + (count-1)*(hom-het).
// REF has the zero-point prior.  The resulting vector follows the same
// Number=G ordering as enumerate_genotypes and is consumed by the shared
// Kokkos posterior-assignment API.
std::vector<double> genotype_priors_for_group(
    fastgatk::calling::Result& result,
    const std::vector<const fastgatk::calling::AssemblyCandidate*>& candidates,
    int ploidy,
    bool include_non_ref,
    bool include_spanning_deletion = false) {
    const auto allele_count = static_cast<int>(
        candidates.size() + 1 + (include_spanning_deletion ? 1 : 0) +
        (include_non_ref ? 1 : 0));
    const auto genotypes = enumerate_genotypes(allele_count, ploidy);
    std::vector<double> priors(genotypes.size(), 0.0);
    if (!result.genotype_priors_used ||
        result.candidate_prior_het.size() != result.candidates.size() ||
        result.candidate_prior_hom_alt.size() != result.candidates.size())
        return priors;
    std::vector<std::pair<double, double>> allele_priors;
    allele_priors.reserve(candidates.size() + (include_spanning_deletion ? 1 : 0) +
                          (include_non_ref ? 1 : 0));
    for (const auto* candidate : candidates) {
        if (candidate == nullptr) return {};
        const auto index = candidate_index_for_candidate(result, *candidate);
        if (index >= result.candidate_prior_het.size()) return {};
        allele_priors.emplace_back(result.candidate_prior_het[index],
                                   result.candidate_prior_hom_alt[index]);
    }
    if (include_spanning_deletion) {
        const auto spanning_heterozygosity = result.genotype_indel_heterozygosity;
        if (!(spanning_heterozygosity > 0.0) || !std::isfinite(spanning_heterozygosity))
            return {};
        const auto log10_spanning_het = std::log10(std::max(spanning_heterozygosity, 1e-300));
        allele_priors.emplace_back(log10_spanning_het, 2.0 * log10_spanning_het);
    }
    if (include_non_ref) {
        const auto other_het = std::max(result.genotype_snp_heterozygosity,
                                        result.genotype_indel_heterozygosity);
        if (!(other_het > 0.0) || !std::isfinite(other_het))
            return {};
        const auto log10_other_het = std::log10(std::max(other_het, 1e-300));
        allele_priors.emplace_back(log10_other_het, 2.0 * log10_other_het);
    }
    std::vector<double> het(static_cast<std::size_t>(allele_count), 0.0);
    std::vector<double> hom(static_cast<std::size_t>(allele_count), 0.0);
    for (std::size_t allele = 1; allele < static_cast<std::size_t>(allele_count); ++allele) {
        het[allele] = allele_priors[allele - 1].first;
        hom[allele] = allele_priors[allele - 1].second;
    }
    const auto calculated = fastgatk::kernels::calculate_genotype_priors_kokkos(
        het, hom, ploidy);
    ++result.genotype_prior_kernel_calls;
    result.genotype_prior_prepare_seconds += calculated.prepare_seconds;
    result.genotype_prior_seconds += calculated.seconds;
    result.genotype_prior_execution_space = calculated.execution_space;
    return calculated.log10_priors;
}

// Extend the concrete diploid matrix with the GATK
// AlleleLikelihoods.updateNonRefAlleleLikelihoods result.  The symbolic
// <NON_REF> allele is not an ordinary alternate: its per-read value is the
// median of qualified concrete likelihoods, so it must be calculated before
// the genotype-combination sum rather than filled with a post-hoc envelope.
bool joint_candidate_pl_with_nonref(
    const fastgatk::calling::Result& result,
    const std::vector<const fastgatk::calling::AssemblyCandidate*>& candidates,
    int ploidy,
    std::vector<int>& output,
    int max_phred = 0,
    bool include_spanning_deletion = false) {
    if (candidates.empty() || result.allele_read_likelihoods.empty() ||
        result.reference_read_likelihoods.empty() ||
        result.non_ref_read_likelihoods.empty() || ploidy < 1) return false;
    const auto concrete_count = candidates.size() + 1;
    const auto spanning_index = concrete_count;
    const auto non_ref_index = concrete_count + (include_spanning_deletion ? 1U : 0U);
    const auto allele_count = non_ref_index + 1U;
    const auto genotypes = enumerate_genotypes(static_cast<int>(allele_count), ploidy);
    if (genotypes.empty()) return false;
    std::vector<std::size_t> candidate_indices;
    candidate_indices.reserve(candidates.size());
    std::size_t read_count = std::numeric_limits<std::size_t>::max();
    for (const auto* candidate : candidates) {
        if (candidate == nullptr) return false;
        const auto index = candidate_index_for_candidate(result, *candidate);
        if (index == std::numeric_limits<std::size_t>::max() ||
            index >= result.allele_read_likelihoods.size() ||
            index >= result.reference_read_likelihoods.size() ||
            index >= result.non_ref_read_likelihoods.size()) return false;
        candidate_indices.push_back(index);
        read_count = std::min(read_count, result.allele_read_likelihoods[index].size());
        read_count = std::min(read_count, result.reference_read_likelihoods[index].size());
        read_count = std::min(read_count, result.non_ref_read_likelihoods[index].size());
    }
    if (read_count == 0 || read_count == std::numeric_limits<std::size_t>::max()) return false;
    std::vector<std::vector<double>> likelihoods(allele_count,
        std::vector<double>(read_count, -std::numeric_limits<double>::infinity()));
    std::size_t non_ref_reads = 0;
    std::size_t spanning_reads = 0;
    for (std::size_t read = 0; read < read_count; ++read) {
        double reference = -std::numeric_limits<double>::infinity();
        for (const auto index : candidate_indices)
            reference = std::max(reference,
                                 result.reference_read_likelihoods[index][read]);
        likelihoods[0][read] = reference;
        for (std::size_t alt = 0; alt < candidate_indices.size(); ++alt)
            likelihoods[alt + 1][read] =
                result.allele_read_likelihoods[candidate_indices[alt]][read];
        if (include_spanning_deletion) {
            double spanning = -std::numeric_limits<double>::infinity();
            for (const auto index : candidate_indices) {
                if (index < result.spanning_deletion_read_likelihoods.size() &&
                    read < result.spanning_deletion_read_likelihoods[index].size()) {
                    spanning = std::max(spanning,
                        result.spanning_deletion_read_likelihoods[index][read]);
                }
            }
            likelihoods[spanning_index][read] = spanning;
            if (std::isfinite(spanning)) ++spanning_reads;
        }
        likelihoods[non_ref_index][read] =
            result.non_ref_read_likelihoods[candidate_indices.front()][read];
        if (std::isfinite(likelihoods[non_ref_index][read])) ++non_ref_reads;
    }
    if (non_ref_reads == 0 || (include_spanning_deletion && spanning_reads == 0)) return false;
    std::vector<double> flattened(allele_count * read_count,
        -std::numeric_limits<double>::infinity());
    for (std::size_t allele = 0; allele < allele_count; ++allele)
        std::copy(likelihoods[allele].begin(), likelihoods[allele].end(),
                  flattened.begin() + allele * read_count);
    const auto joint = fastgatk::kernels::calculate_joint_genotype_pl_kokkos(
        flattened, read_count, static_cast<int>(allele_count), ploidy, max_phred);
    output.assign(joint.pl.begin(), joint.pl.end());
    return true;
}

bool joint_concrete_pl(const fastgatk::calling::Result& result,
                       const std::vector<const fastgatk::calling::GenotypeCall*>& calls,
                       int ploidy,
                       std::vector<int>& output) {
    std::vector<const fastgatk::calling::AssemblyCandidate*> candidates;
    candidates.reserve(calls.size());
    for (const auto* call : calls) {
        if (call == nullptr) return false;
        candidates.push_back(&call->candidate);
    }
    return joint_candidate_pl(result, candidates, ploidy, output);
}

bool has_suffix(const std::string& path, const char* suffix) {
    const std::string ending(suffix);
    return path.size() >= ending.size() &&
           path.compare(path.size() - ending.size(), ending.size(), ending) == 0;
}

bool is_vcf_output(const std::string& path) {
    return has_suffix(path, ".vcf") || has_suffix(path, ".vcf.gz") ||
           has_suffix(path, ".g.vcf") || has_suffix(path, ".g.vcf.gz");
}

Options parse(int argc, char** argv) {
    Options options;
    // GATK ReadThreadingAssembler starts with k-mer 10 (and may retry with
    // larger values on cycles).  Keep the shared API's historical k=5 smoke
    // default, but make the production HC CLI start at the GATK value.
    options.calling.graph_kmer_size = 10;
    // HaplotypeCaller applies GoodCigarReadFilter and
    // NonZeroReferenceLengthAlignmentReadFilter as implicit standard
    // filters. Keep these enabled for the CLI while leaving the shared API
    // configurable for low-level tests and custom callers.
    options.calling.require_good_cigar = true;
    options.calling.require_nonzero_reference_span = true;
    options.calling.require_no_n_cigar = true;
    options.calling.require_read_group = true;
    options.calling.exclude_mapping_quality_unavailable = true;
    // GATK HC does not impose a fixed pileup-support gate before assembly;
    // leave low-support loci available to PairHMM/AF confidence.  The shared
    // calling API retains stricter smoke defaults for direct library users.
    options.calling.min_depth = 1;
    options.calling.min_alt_support = 1;
    options.calling.graph_kmer_sizes = {10, 25};
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-hc-call (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 input SAM/BAM/CRAM (repeatable; shards are combined)\n"
                         "  -R, --reference FILE             reference FASTA\n"
                         "      --alleles VCF                force-call concrete eligible alleles from VCF\n"
                         "      --force-call-filtered-alleles[=BOOL]  include filtered --alleles records\n"
                         "      --sample-name NAME           select one @RG SM sample from a multi-sample BAM\n"
                         "  -L, --intervals REGION           contig:start-end\n"
                         "      -isr, --interval-set-rule R  UNION (default) or INTERSECTION\n"
                         "  -XL, --exclude-intervals REGION  exclude loci after -L selection\n"
                         "      -ip, --interval-padding N    include interval padding\n"
                         "      -ixp, --interval-exclusion-padding N  exclusion padding\n"
                         "  -O, --output FILE                VCF or JSON summary\n"
                         "      --assembly-region-out FILE  IGV activity/profile regions\n"
                         "      --force-active[=BOOL]       force all activity segments active (default false)\n"
                         "      --create-output-variant-index [BOOL]  create .tbi for BGZF VCF/GVCF (default true)\n"
                         "      --sites-only-vcf-output [BOOL]  omit FORMAT/sample columns (default false)\n"
                         "      --add-output-vcf-command-line [BOOL]  accept GATK provenance switch\n"
                         "      --disable-sequence-dictionary-validation [BOOL]  skip BAM/reference dictionary check\n"
                         "      --output-manifest FILE       OutputManifest JSON\n"
                         "      --telemetry FILE             JSON telemetry summary\n"
                         "      --batch-records N             streaming batch size\n"
                         "      --stream-by-contig            flush coordinate-sorted contigs incrementally\n"
                         "      --stream-by-region N          indexed halo/core AssemblyRegion tiles (fixed Host bound)\n"
                         "      --threads N                  Kokkos worker count\n"
                         "      --min-depth N                candidate depth threshold\n"
                         "      --min-alt-support N          alternate support threshold\n"
                         "      --active-probability-threshold P  AssemblyRegion activity threshold\n"
                         "      --assembly-region-padding N  AssemblyRegion halo/padding\n"
                         "      --min-assembly-region-size N  minimum AssemblyRegion size (default 50)\n"
                         "      --max-assembly-region-size N  hard AssemblyRegion length cap\n"
                         "      --max-prob-propagation-distance N  active-locus grouping distance\n"
                         "      --minimum-mapping-quality Q  GATK read-filter MAPQ floor (default 20)\n"
                         "      --min-read-length N          ReadLengthReadFilter lower bound\n"
                         "      --max-read-length N          ReadLengthReadFilter upper bound\n"
                         "      --include-duplicates         keep duplicate reads (debug/compatibility)\n"
                         "      --dont-use-soft-clipped-bases  exclude soft-clipped bases from assembly evidence\n"
                         "      --soft-clip-low-quality-ends  preserve terminal low-quality bases as assembly soft clips\n"
                         "      --use-haplotype-realignment-for-rcm  project best PairHMM haplotypes into RCM\n"
                         "      --do-not-correct-overlapping-quality  disable GATK overlapping-mate quality correction\n"
                         "      --max-reads-per-locus N     positional reservoir downsampling cap (0=unbounded)\n"
                         "      --max-reads-per-alignment-start N  GATK positional downsampling cap (default 50)\n"
                         "      --downsampling-seed N      Java-compatible reservoir seed\n"
                         "      --min-base-quality-score Q, -mbq  base quality threshold\n"
                         "      --reference-model-deletion-quality Q  gVCF deletion evidence quality (default 30)\n"
                         "      --base-quality-score-threshold Q  PairHMM base-quality floor (default 18)\n"
                         "      --disable-cap-base-qualities-to-map-quality  disable PairHMM MAPQ cap\n"
                         "      --phred-scaled-global-read-mismapping-rate Q  per-read likelihood cap (default 45; negative disables)\n"
                         "      --pcr-indel-model M       NONE, HOSTILE, AGGRESSIVE, or CONSERVATIVE (default)\n"
                         "      --flow-assembly-collapse-hmer-size N  cap flow haplotype homopolymers (-1=auto, 0=off)\n"
                         "      --flow-assembly-collapse-partial-mode  stop HMER restoration at reference mismatch\n"
                         "      --max-candidates N           assembly candidate cap\n"
                         "      --max-mnp-distance N         merge phased graph SNPs into an MNP (default 0)\n"
                         "      --mnp-dist N                 short alias for --max-mnp-distance\n"
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
                         "      --min-dangling-matching-bases N  minimum exact bases for dangling recovery (default 3)\n"
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
                         "      --haplotype-pruning-log10 D graph-path likelihood pruning delta\n"
                         "      --sample-ploidy N            genotype ploidy (default 2; VCF 1..8)\n"
                         "      --max-genotype-count N       maximum genotypes per HC EventMap site (default 1024)\n"
                         "      --max-alternate-alleles N    maximum concrete alternate alleles to genotype (default 6)\n"
                         "      --heterozygosity P           SNP genotype prior (default 1e-3)\n"
                         "      --indel-heterozygosity P     indel genotype prior (default 1/8000)\n"
                         "      --heterozygosity-stdev P     AF prior spread (default 0.01)\n"
                         "      --standard-min-confidence-threshold-for-calling Q (default 30; GVCF=0)\n"
                         "      --stand-call-conf Q          short alias for the call threshold\n"
                         "      --genotype-assignment-method M  USE_PLS_TO_ASSIGN (default) or USE_POSTERIOR_PROBABILITIES\n"
                         "      --indel-size-to-eliminate-in-ref-model N  RCM indel bound (default 10)\n"
                         "      --allele-informative-reads-overlap-margin N  retain-evidence overlap margin (default 2)\n"
                         "      --no-genotype-priors        disable native genotype priors\n"
                         "      --gvcf / -ERC GVCF           emit GVCF reference blocks\n"
                         "      --floor-blocks [BOOL]        emit compact GT:DP:GQ reference blocks\n"
                         "      -ERC BP_RESOLUTION           emit one reference-confidence record per base\n"
                         "      -GQB, --GVCF-GQ-bands N       repeatable GVCF block upper bound (1..100)\n"
                         "      --emit-ref-confidence MODE  alias for -ERC/--gvcf\n"
                         "The output is a deterministic native prototype, not yet bit-identical\n"
                         "to GATK HaplotypeCaller.\n";
            std::exit(0);
        } else if (is_option(argument, "--input") || argument == "-I") {
            const auto input = require_value(i, argc, argv, argument, "--input", "-I");
            if (options.input.empty()) options.input = input;
            options.inputs.push_back(input);
        }
        else if (is_option(argument, "--reference") || argument == "-R")
            options.reference = require_value(i, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--alleles")) {
            options.alleles = require_value(i, argc, argv, argument, "--alleles");
            options.compatibility_options.push_back("--alleles=" + options.alleles);
        }
        else if (argument == "--force-call-filtered-alleles" ||
                 argument.rfind("--force-call-filtered-alleles=", 0) == 0 ||
                 argument == "--genotype-filtered-alleles" ||
                 argument.rfind("--genotype-filtered-alleles=", 0) == 0) {
            const char* name = argument.rfind("--genotype-filtered-alleles", 0) == 0
                ? "--genotype-filtered-alleles" : "--force-call-filtered-alleles";
            options.force_call_filtered_alleles =
                fastgatk::native::parse_optional_boolean(i, argc, argv, argument, name);
            options.compatibility_options.push_back(
                std::string(name) + "=" +
                (options.force_call_filtered_alleles ? "true" : "false"));
        }
        else if (is_option(argument, "--sample-name"))
            options.sample_name = require_value(i, argc, argv, argument, "--sample-name");
        else if (is_option(argument, "--region") || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || argument == "-L") {
            if (argument == "-L" || argument == "--intervals" || !option_value(argument, "--intervals").empty())
                options.regions.push_back(require_value(i, argc, argv, argument, "--intervals", "-L"));
            else if (argument == "--interval" || !option_value(argument, "--interval").empty())
                options.regions.push_back(require_value(i, argc, argv, argument, "--interval"));
            else options.regions.push_back(require_value(i, argc, argv, argument, "--region"));
        } else if (argument == "-XL" || is_option(argument, "--exclude-intervals")) {
            options.exclusions.push_back(require_value(i, argc, argv, argument,
                                                        "--exclude-intervals", "-XL"));
        } else if (argument == "-isr" || is_option(argument, "--interval-set-rule")) {
            auto rule = require_value(i, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char value) {
                return static_cast<char>(std::toupper(value));
            });
            if (rule == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (rule == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument(
                "BAD_INPUT: --interval-set-rule must be UNION or INTERSECTION");
        } else if (argument == "-ip" || is_option(argument, "--interval-padding")) {
            options.interval_padding = std::stoll(require_value(
                i, argc, argv, argument, "--interval-padding", "-ip"));
            if (options.interval_padding < 0)
                throw std::invalid_argument("--interval-padding must be non-negative");
        } else if (argument == "-ixp" || is_option(argument, "--interval-exclusion-padding")) {
            options.exclusion_padding = std::stoll(require_value(
                i, argc, argv, argument, "--interval-exclusion-padding", "-ixp"));
            if (options.exclusion_padding < 0)
                throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
        } else if (is_option(argument, "--output") || argument == "-O")
            options.output = require_value(i, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--assembly-region-out"))
            options.assembly_region_out = require_value(
                i, argc, argv, argument, "--assembly-region-out");
        else if (argument == "--force-active" ||
                 argument.rfind("--force-active=", 0) == 0) {
            options.calling.force_active = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--force-active");
            options.compatibility_options.push_back(
                std::string("--force-active=") +
                (options.calling.force_active ? "true" : "false"));
        } else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(i, argc, argv, argument,
                argument == "--manifest" ? "--manifest" : "--output-manifest");
        else if (is_option(argument, "--telemetry"))
            options.telemetry = require_value(i, argc, argv, argument, "--telemetry");
        else if (is_option(argument, "--batch-records"))
            options.batch_records = std::stoull(require_value(i, argc, argv, argument, "--batch-records"));
        else if (argument == "--stream-by-contig")
            options.stream_by_contig = true;
        else if (is_option(argument, "--stream-by-region")) {
            options.stream_by_region = true;
            options.stream_region_size = std::stoull(require_value(
                i, argc, argv, argument, "--stream-by-region"));
            if (options.stream_region_size == 0)
                throw std::invalid_argument("--stream-by-region tile size must be positive");
        }
        else if (is_option(argument, "--threads") || is_option(argument, "--native-pair-hmm-threads")) {
            const bool native_threads = argument == "--native-pair-hmm-threads" ||
                                        !option_value(argument, "--native-pair-hmm-threads").empty();
            options.threads = std::stoi(require_value(i, argc, argv, argument,
                native_threads ? "--native-pair-hmm-threads" : "--threads"));
        }
        else if (is_option(argument, "--min-depth")) {
            options.calling.min_depth = std::stoull(require_value(i, argc, argv, argument, "--min-depth"));
            options.min_depth_explicit = true;
        }
        else if (is_option(argument, "--min-alt-support")) {
            options.calling.min_alt_support = std::stoull(require_value(i, argc, argv, argument, "--min-alt-support"));
            options.min_alt_support_explicit = true;
        }
        else if (is_option(argument, "--active-probability-threshold")) {
            options.calling.active_probability_threshold = std::stod(require_value(
                i, argc, argv, argument, "--active-probability-threshold"));
            if (!std::isfinite(options.calling.active_probability_threshold) ||
                options.calling.active_probability_threshold < 0.0 ||
                options.calling.active_probability_threshold > 1.0)
                throw std::invalid_argument(
                    "--active-probability-threshold must be finite and in [0,1]");
        }
        else if (is_option(argument, "--assembly-region-padding"))
            options.calling.assembly_region_padding = static_cast<std::uint32_t>(std::stoul(
                require_value(i, argc, argv, argument, "--assembly-region-padding")));
        else if (is_option(argument, "--min-assembly-region-size"))
            options.calling.min_assembly_region_size = static_cast<std::uint32_t>(std::stoul(
                require_value(i, argc, argv, argument, "--min-assembly-region-size")));
        else if (is_option(argument, "--max-assembly-region-size"))
            options.calling.max_assembly_region_size = static_cast<std::uint32_t>(std::stoul(
                require_value(i, argc, argv, argument, "--max-assembly-region-size")));
        else if (is_option(argument, "--max-prob-propagation-distance"))
            options.calling.max_probability_propagation_distance = static_cast<std::uint32_t>(std::stoul(
                require_value(i, argc, argv, argument, "--max-prob-propagation-distance")));
        else if (is_option(argument, "--minimum-mapping-quality")) {
            const auto value = std::stoul(require_value(i, argc, argv, argument,
                                                        "--minimum-mapping-quality"));
            if (value > 255) throw std::invalid_argument("--minimum-mapping-quality must be <=255");
            options.calling.minimum_mapping_quality = static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--min-read-length")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--min-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-read-length is out of range");
            options.calling.min_read_length = static_cast<std::uint32_t>(value);
            options.calling.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        }
        else if (is_option(argument, "--max-read-length")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--max-read-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-read-length is out of range");
            options.calling.max_read_length = static_cast<std::uint32_t>(value);
            options.calling.require_read_length = true;
            options.read_filters.push_back("ReadLengthReadFilter");
        }
        else if (argument == "--include-duplicates")
            options.calling.exclude_duplicates = false;
        else if (argument == "--dont-use-soft-clipped-bases" ||
                 argument.rfind("--dont-use-soft-clipped-bases=", 0) == 0) {
            const auto disabled = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--dont-use-soft-clipped-bases");
            options.calling.use_soft_clipped_bases = !disabled;
            options.compatibility_options.push_back(
                std::string("--dont-use-soft-clipped-bases=") +
                (disabled ? "true" : "false"));
        }
        else if (argument == "--soft-clip-low-quality-ends" ||
                 argument.rfind("--soft-clip-low-quality-ends=", 0) == 0) {
            const auto enabled = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--soft-clip-low-quality-ends");
            options.calling.soft_clip_low_quality_ends = enabled;
            options.compatibility_options.push_back(
                std::string("--soft-clip-low-quality-ends=") +
                (enabled ? "true" : "false"));
        }
        else if (argument == "--use-haplotype-realignment-for-rcm") {
            options.calling.use_haplotype_realignment_for_rcm = true;
            options.compatibility_options.push_back(argument + "=true");
        }
        else if (argument == "--do-not-correct-overlapping-quality" ||
                 argument.rfind("--do-not-correct-overlapping-quality=", 0) == 0 ||
                 argument == "--do-not-correct-overlapping-base-qualities" ||
                 argument.rfind("--do-not-correct-overlapping-base-qualities=", 0) == 0) {
            const char* name = argument.rfind("--do-not-correct-overlapping-base-qualities", 0) == 0
                ? "--do-not-correct-overlapping-base-qualities"
                : "--do-not-correct-overlapping-quality";
            options.calling.do_not_correct_overlapping_base_qualities =
                fastgatk::native::parse_optional_boolean(i, argc, argv, argument, name);
        }
        else if (is_option(argument, "--max-reads-per-locus") ||
                 is_option(argument, "--max-reads-per-alignment-start")) {
            const auto name = argument.rfind("--max-reads-per-alignment-start", 0) == 0
                ? "--max-reads-per-alignment-start" : "--max-reads-per-locus";
            options.calling.max_reads_per_locus = std::stoull(
                require_value(i, argc, argv, argument, name));
        }
        else if (is_option(argument, "--downsampling-seed"))
            options.calling.downsampling_seed = std::stoull(
                require_value(i, argc, argv, argument, "--downsampling-seed"));
        else if (is_option(argument, "--min-base-quality") ||
                 is_option(argument, "--min-base-quality-score") || argument == "-mbq") {
            const char* name = argument == "-mbq" ? "--min-base-quality-score" :
                (argument.rfind("--min-base-quality-score", 0) == 0
                    ? "--min-base-quality-score" : "--min-base-quality");
            const auto value = std::stoul(require_value(
                i, argc, argv, argument, name, argument == "-mbq" ? "-mbq" : nullptr));
            if (value > 255)
                throw std::invalid_argument("--min-base-quality-score must be <=255");
            options.calling.min_base_quality = static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--reference-model-deletion-quality")) {
            const auto value = std::stoul(require_value(
                i, argc, argv, argument, "--reference-model-deletion-quality"));
            // GATK exposes this as a Java byte.  Keep unsupported negative
            // and out-of-range values fail-closed instead of wrapping them
            // into an unrelated unsigned phred score.
            if (value > 127)
                throw std::invalid_argument(
                    "--reference-model-deletion-quality must be in [0,127]");
            options.calling.reference_model_deletion_quality =
                static_cast<std::uint8_t>(value);
        }
        else if (is_option(argument, "--base-quality-score-threshold")) {
            const auto value = std::stoul(require_value(i, argc, argv, argument,
                                                         "--base-quality-score-threshold"));
            if (value < 6 || value > 255)
                throw std::invalid_argument("--base-quality-score-threshold must be in [6,255]");
            options.calling.pairhmm_base_quality_score_threshold =
                static_cast<std::uint8_t>(value);
        }
        else if (argument == "--disable-cap-base-qualities-to-map-quality")
            options.calling.pairhmm_disable_cap_base_qualities_to_mapq = true;
        else if (is_option(argument, "--phred-scaled-global-read-mismapping-rate")) {
            options.calling.phred_scaled_global_read_mismapping_rate = std::stod(
                require_value(i, argc, argv, argument,
                              "--phred-scaled-global-read-mismapping-rate"));
            if (!std::isfinite(options.calling.phred_scaled_global_read_mismapping_rate))
                throw std::invalid_argument(
                    "--phred-scaled-global-read-mismapping-rate must be finite");
        }
        else if (is_option(argument, "--pcr-indel-model")) {
            auto model = require_value(i, argc, argv, argument, "--pcr-indel-model");
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
            options.compatibility_options.push_back("--pcr-indel-model=" + model);
        }
        else if (is_option(argument, "--flow-assembly-collapse-hmer-size")) {
            const auto value = std::stoll(require_value(i, argc, argv, argument,
                                                         "--flow-assembly-collapse-hmer-size"));
            if (value < -1 || value > 255)
                throw std::invalid_argument(
                    "--flow-assembly-collapse-hmer-size must be in [-1,255]");
            options.calling.flow_assembly_collapse_hmer_size =
                static_cast<std::int32_t>(value);
            options.compatibility_options.push_back(
                "--flow-assembly-collapse-hmer-size=" + std::to_string(value));
        }
        else if (argument == "--flow-assembly-collapse-partial-mode") {
            options.calling.flow_assembly_collapse_partial_mode = true;
            options.compatibility_options.push_back(argument);
        }
        else if (is_option(argument, "--expected-mismatch-rate-for-read-disqualification")) {
            options.calling.pairhmm_expected_error_rate_per_base = std::stod(
                require_value(i, argc, argv, argument,
                              "--expected-mismatch-rate-for-read-disqualification"));
            if (!std::isfinite(options.calling.pairhmm_expected_error_rate_per_base) ||
                options.calling.pairhmm_expected_error_rate_per_base < 0.0)
                throw std::invalid_argument(
                    "--expected-mismatch-rate-for-read-disqualification must be finite and non-negative");
        }
        else if (is_option(argument, "--max-candidates"))
            options.calling.max_candidates = std::stoull(require_value(i, argc, argv, argument, "--max-candidates"));
        else if (is_option(argument, "--max-mnp-distance") || is_option(argument, "--mnp-dist")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument,
                                                          argument.rfind("--mnp-dist", 0) == 0
                                                              ? "--mnp-dist" : "--max-mnp-distance"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-mnp-distance must fit uint32");
            options.calling.max_mnp_distance = static_cast<std::uint32_t>(value);
            options.compatibility_options.push_back(argument +
                (argument.find('=') == std::string::npos ? "=" + std::to_string(value) : ""));
        }
        else if (is_option(argument, "--max-alternate-alleles")) {
            const auto value = std::stoull(require_value(
                i, argc, argv, argument, "--max-alternate-alleles"));
            if (value == 0)
                throw std::invalid_argument("--max-alternate-alleles must be positive");
            options.max_alternate_alleles = static_cast<std::size_t>(value);
            options.compatibility_options.push_back(argument +
                (argument.find('=') == std::string::npos ? "=" + std::to_string(value) : ""));
        }
        else if (is_option(argument, "--max-genotype-count")) {
            const auto value = std::stoull(require_value(
                i, argc, argv, argument, "--max-genotype-count"));
            if (value == 0 || value > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
                throw std::invalid_argument(
                    "--max-genotype-count must be a positive signed 32-bit integer");
            options.max_genotype_count = static_cast<std::size_t>(value);
            options.compatibility_options.push_back(argument +
                (argument.find('=') == std::string::npos ? "=" + std::to_string(value) : ""));
        }
        else if (is_option(argument, "--allele-informative-reads-overlap-margin")) {
            const auto value = std::stoull(require_value(
                i, argc, argv, argument, "--allele-informative-reads-overlap-margin"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument(
                    "--allele-informative-reads-overlap-margin must fit uint32");
            options.calling.informative_read_overlap_margin =
                static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--kmer-size")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--kmer-size"));
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
                    i, argc, argv, argument, "--dont-increase-kmer-sizes-for-cycles");
        else if (is_option(argument, "--min-kmer-count")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--min-kmer-count"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-kmer-count must be positive and fit uint32");
            options.calling.graph_min_kmer_count = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-pruning")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--min-pruning"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-pruning must fit uint32");
            options.calling.graph_min_pruning = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--num-pruning-samples")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--num-pruning-samples"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--num-pruning-samples must be positive and fit uint32");
            options.calling.graph_num_pruning_samples = static_cast<std::uint32_t>(value);
        }
        else if (argument == "--adaptive-pruning")
            options.calling.graph_use_adaptive_pruning = true;
        else if (argument == "--disable-adaptive-pruning")
            options.calling.graph_use_adaptive_pruning = false;
        else if (argument == "--linked-de-bruijn-graph")
            options.calling.graph_linked_de_bruijn = true;
        else if (argument == "--disable-artificial-haplotype-recovery")
            options.calling.graph_disable_artificial_haplotype_recovery = true;
        else if (argument == "--enable-legacy-graph-cycle-detection")
            options.calling.graph_enable_legacy_cycle_detection = true;
        else if (is_option(argument, "--adaptive-pruning-initial-error-rate"))
            options.calling.graph_initial_error_rate_for_pruning = std::stod(
                require_value(i, argc, argv, argument, "--adaptive-pruning-initial-error-rate"));
        else if (is_option(argument, "--pruning-lod-threshold"))
            options.calling.graph_pruning_log_odds_threshold = std::stod(
                require_value(i, argc, argv, argument, "--pruning-lod-threshold"));
        else if (is_option(argument, "--pruning-seeding-lod-threshold"))
            options.calling.graph_pruning_seeding_log_odds_threshold = std::stod(
                require_value(i, argc, argv, argument, "--pruning-seeding-lod-threshold"));
        else if (is_option(argument, "--max-unpruned-variants")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--max-unpruned-variants"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--max-unpruned-variants must be positive and fit uint32");
            options.calling.graph_max_unpruned_variants = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-dangling-branch-length")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument,
                                                          "--min-dangling-branch-length"));
            if (value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-dangling-branch-length must fit uint32");
            options.calling.graph_min_dangling_branch_length = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-dangling-matching-bases")) {
            const auto value = std::stoll(require_value(i, argc, argv, argument,
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
                    i, argc, argv, argument, "--recover-all-dangling-branches");
        else if (argument == "--error-correct-reads" ||
                 argument.rfind("--error-correct-reads=", 0) == 0)
            options.calling.error_correct_reads = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--error-correct-reads");
        else if (is_option(argument, "--error-correction-log-odds"))
            options.calling.pileup_error_correction_log_odds = std::stod(
                require_value(i, argc, argv, argument, "--error-correction-log-odds"));
        else if (is_option(argument, "--kmer-length-for-read-error-correction")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument,
                                                          "--kmer-length-for-read-error-correction"));
            if (value == 0 || value > 31)
                throw std::invalid_argument("--kmer-length-for-read-error-correction must be in [1,31]");
            options.calling.error_correction_kmer_length = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--min-observations-for-kmer-to-be-solid")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument,
                                                          "--min-observations-for-kmer-to-be-solid"));
            if (value == 0 || value > std::numeric_limits<std::uint32_t>::max())
                throw std::invalid_argument("--min-observations-for-kmer-to-be-solid must be positive and fit uint32");
            options.calling.error_correction_min_solid_observations = static_cast<std::uint32_t>(value);
        }
        else if (is_option(argument, "--max-haplotype-paths") ||
                 is_option(argument, "--max-num-haplotypes-in-population")) {
            const auto name = is_option(argument, "--max-num-haplotypes-in-population")
                ? "--max-num-haplotypes-in-population" : "--max-haplotype-paths";
            const auto value = std::stoull(require_value(i, argc, argv, argument, name));
            if (value == 0) throw std::invalid_argument("--max-haplotype-paths must be positive");
            options.calling.graph_max_paths = value;
        }
        else if (is_option(argument, "--max-haplotype-depth")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument, "--max-haplotype-depth"));
            if (value == 0) throw std::invalid_argument("--max-haplotype-depth must be positive");
            options.calling.graph_max_depth = value;
        }
        else if (is_option(argument, "--max-haplotype-combination-alleles")) {
            const auto value = std::stoull(require_value(i, argc, argv, argument,
                                                          "--max-haplotype-combination-alleles"));
            if (value == 0 || value > 16)
                throw std::invalid_argument("--max-haplotype-combination-alleles must be in [1,16]");
            options.calling.max_haplotype_combination_alleles = value;
        }
        else if (is_option(argument, "--haplotype-pruning-log10"))
            options.calling.haplotype_pruning_log10 = std::stod(
                require_value(i, argc, argv, argument, "--haplotype-pruning-log10"));
        else if (is_option(argument, "--sample-ploidy")) {
            const auto value = require_value(i, argc, argv, argument, "--sample-ploidy");
            try { options.sample_ploidy = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --sample-ploidy: " + value); }
            if (options.sample_ploidy < 1 || options.sample_ploidy > 8)
                throw std::invalid_argument("--sample-ploidy must be in [1,8]");
        }
        else if (is_option(argument, "--heterozygosity"))
            options.calling.heterozygosity = std::stod(require_value(i, argc, argv, argument, "--heterozygosity"));
        else if (is_option(argument, "--indel-heterozygosity"))
            options.calling.indel_heterozygosity = std::stod(require_value(i, argc, argv, argument, "--indel-heterozygosity"));
        else if (is_option(argument, "--heterozygosity-stdev"))
            options.calling.heterozygosity_stdev = std::stod(require_value(
                i, argc, argv, argument, "--heterozygosity-stdev"));
        else if (is_option(argument, "--standard-min-confidence-threshold-for-calling") ||
                 is_option(argument, "--stand-call-conf")) {
            const auto name = is_option(argument, "--stand-call-conf")
                ? "--stand-call-conf" : "--standard-min-confidence-threshold-for-calling";
            options.calling.standard_confidence_for_calling = std::stod(
                require_value(i, argc, argv, argument, name));
            if (!std::isfinite(options.calling.standard_confidence_for_calling) ||
                options.calling.standard_confidence_for_calling < 0.0)
                throw std::invalid_argument("call confidence threshold must be finite and non-negative");
        }
        else if (is_option(argument, "--genotype-assignment-method")) {
            const auto method = require_value(i, argc, argv, argument,
                                               "--genotype-assignment-method");
            if (method == "USE_POSTERIOR_PROBABILITIES") {
                options.calling.use_posterior_genotype_assignment = true;
            } else if (method == "USE_PLS_TO_ASSIGN" || method == "PREFER_PLS") {
                options.calling.use_posterior_genotype_assignment = false;
            } else {
                throw std::runtime_error(
                    "UNSUPPORTED_PARAMETER: --genotype-assignment-method " + method);
            }
            options.compatibility_options.push_back(
                "--genotype-assignment-method=" + method);
        }
        else if (is_option(argument, "--indel-size-to-eliminate-in-ref-model"))
            options.calling.indel_size_to_eliminate_in_ref_model = std::stoi(require_value(
                i, argc, argv, argument, "--indel-size-to-eliminate-in-ref-model"));
        else if (argument == "--no-genotype-priors")
            options.calling.use_genotype_priors = false;
        else if (argument == "-GQB" || is_option(argument, "--GVCF-GQ-bands") ||
                 is_option(argument, "--gvcf-gq-bands")) {
            const auto long_name = is_option(argument, "--gvcf-gq-bands")
                ? "--gvcf-gq-bands" : "--GVCF-GQ-bands";
            const auto text = require_value(i, argc, argv, argument, long_name, "-GQB");
            std::size_t consumed = 0;
            long long value = 0;
            try { value = std::stoll(text, &consumed); }
            catch (...) { throw std::invalid_argument("invalid GVCF GQ band: " + text); }
            if (consumed != text.size() || value <= 0 || value > 100)
                throw std::invalid_argument("GVCF GQ bands must be positive, <=100, and strictly increasing");
            if (!options.calling.gvcf_gq_bands.empty() &&
                value <= options.calling.gvcf_gq_bands.back())
                throw std::invalid_argument("GVCF GQ bands must be strictly increasing");
            options.calling.gvcf_gq_bands.push_back(static_cast<int>(value));
            options.compatibility_options.push_back(
                std::string("--GVCF-GQ-bands=") + std::to_string(value));
        }
        else if (argument == "--gvcf" || argument == "-ERC" ||
                 argument == "--emit-ref-confidence" ||
                 argument.rfind("--gvcf=", 0) == 0 || argument.rfind("-ERC=", 0) == 0 ||
                 argument.rfind("--emit-ref-confidence=", 0) == 0) {
            options.gvcf_mode_explicit = true;
            const auto equals = argument.find('=');
            if (equals != std::string::npos) options.gvcf_mode = argument.substr(equals + 1);
            else if (i + 1 < argc && argv[i + 1][0] != '-') options.gvcf_mode = argv[++i];
            else options.gvcf_mode = "GVCF";
            const auto normalized_mode = normalize_ref_confidence_mode(options.gvcf_mode);
            if (normalized_mode == "NONE") {
                options.gvcf = false;
                options.calling.emit_reference_confidence_bp_resolution = false;
            } else if (normalized_mode == "GVCF" || normalized_mode == "BP_RESOLUTION") {
                options.gvcf = true;
                options.calling.emit_reference_confidence_bp_resolution =
                    is_bp_resolution_mode(normalized_mode);
            } else {
                throw std::invalid_argument("unsupported reference-confidence mode: " +
                                            options.gvcf_mode);
            }
            options.compatibility_options.push_back(argument +
                (options.gvcf_mode.empty() ? "" : "=" + options.gvcf_mode));
        } else if (argument == "--floor-blocks" ||
                   argument.rfind("--floor-blocks=", 0) == 0) {
            options.floor_blocks = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--floor-blocks");
            options.compatibility_options.push_back(
                std::string("--floor-blocks=") + (options.floor_blocks ? "true" : "false"));
        } else if (is_option(argument, "--read-filter")) {
            const auto filter_name = require_value(i, argc, argv, argument, "--read-filter");
            enable_read_filter(options.calling, filter_name);
            options.read_filters.push_back(filter_name);
            options.compatibility_options.push_back(argument +
                (argument.find('=') == std::string::npos ? "=" + filter_name : ""));
        } else if (is_option(argument, "--disable-read-filter")) {
            // GATK takes a filter class as the value.  Require it here so a
            // misspelled/missing value cannot accidentally disable the whole
            // default mask.  Known classes are applied selectively.
            const auto filter_name = require_value(i, argc, argv, argument,
                                                    "--disable-read-filter");
            options.disabled_read_filters.push_back(filter_name);
            options.compatibility_options.push_back(argument +
                (argument.find('=') == std::string::npos ? "=" + filter_name : ""));
            disable_read_filter(options.calling, filter_name);
        } else if (argument == "--disable-tool-default-read-filters" ||
                   argument.rfind("--disable-tool-default-read-filters=", 0) == 0) {
            options.calling.disable_tool_default_read_filters =
                fastgatk::native::parse_optional_boolean(
                    i, argc, argv, argument, "--disable-tool-default-read-filters");
            options.compatibility_options.push_back(
                std::string("--disable-tool-default-read-filters=") +
                (options.calling.disable_tool_default_read_filters ? "true" : "false"));
        } else if (argument == "--create-output-variant-index" ||
                   argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_output_variant_index = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--create-output-variant-index");
            options.compatibility_options.push_back(
                std::string("--create-output-variant-index=") +
                (options.create_output_variant_index ? "true" : "false"));
        } else if (argument == "--add-output-vcf-command-line" ||
                   argument.rfind("--add-output-vcf-command-line=", 0) == 0) {
            options.add_output_vcf_command_line = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--add-output-vcf-command-line");
            options.compatibility_options.push_back(
                std::string("--add-output-vcf-command-line=") +
                (options.add_output_vcf_command_line ? "true" : "false"));
        } else if (argument == "--sites-only-vcf-output" ||
                   argument.rfind("--sites-only-vcf-output=", 0) == 0) {
            options.sites_only_vcf_output = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--sites-only-vcf-output");
            options.compatibility_options.push_back(
                std::string("--sites-only-vcf-output=") +
                (options.sites_only_vcf_output ? "true" : "false"));
        } else if (argument == "--native-pair-hmm-use-double-precision" ||
                   argument.rfind("--native-pair-hmm-use-double-precision=", 0) == 0) {
            options.native_pair_hmm_use_double_precision = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--native-pair-hmm-use-double-precision");
            options.native_pair_hmm_use_double_precision_explicit = true;
            options.compatibility_options.push_back(
                std::string("--native-pair-hmm-use-double-precision=") +
                (options.native_pair_hmm_use_double_precision ? "true" : "false"));
        } else if (argument == "--disable-sequence-dictionary-validation" ||
                   argument.rfind("--disable-sequence-dictionary-validation=", 0) == 0) {
            options.disable_sequence_dictionary_validation = fastgatk::native::parse_optional_boolean(
                i, argc, argv, argument, "--disable-sequence-dictionary-validation");
            options.compatibility_options.push_back(
                std::string("--disable-sequence-dictionary-validation=") +
                (options.disable_sequence_dictionary_validation ? "true" : "false"));
        } else if (argument == "--quiet") {
            // Launcher compatibility switch; keep it visible in the manifest.
            options.compatibility_options.push_back(argument + "=true");
        } else if (is_option(argument, "--seconds-between-progress-updates") ||
                   is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--tmp-dir")) {
            const auto equals = argument.find('=');
            const std::string name = argument.substr(0, equals);
            options.compatibility_options.push_back(argument);
            if (equals == std::string::npos) options.compatibility_options.back() += "=" + require_value(i, argc, argv, argument,
                argument == "--java-options" ? "--java-options" :
                argument == "--verbosity" ? "--verbosity" :
                argument == "--tmp-dir" ? "--tmp-dir" : name.c_str());
        } else {
            throw std::runtime_error("unknown option: " + argument);
        }
    }
    if (options.inputs.empty() || options.input.empty())
        throw std::invalid_argument("--input is required");
    if (options.output.empty()) options.output = "-";
    if (!options.gvcf_mode_explicit &&
        (has_suffix(options.output, ".g.vcf") || has_suffix(options.output, ".g.vcf.gz")))
        options.gvcf = true;
    // The no-reference path is an explicitly documented smoke fallback, not
    // the production HC caller. Preserve its historical bounded candidate
    // floor unless a caller intentionally supplied the native knobs.
    if (options.reference.empty()) {
        if (!options.min_depth_explicit) options.calling.min_depth = 4;
        if (!options.min_alt_support_explicit) options.calling.min_alt_support = 3;
    }
    // HaplotypeCaller deliberately emits low-confidence candidate envelopes in
    // GVCF/BP_RESOLUTION mode; GenotypeGVCFs owns the later cohort threshold.
    if (options.gvcf) options.calling.standard_confidence_for_calling = 0.0;
    if (!options.calling.gvcf_gq_bands.empty() &&
        options.calling.gvcf_gq_bands.back() != 100)
        options.calling.gvcf_gq_bands.push_back(100);
    // HTSJDK normalizes the repeatable ReadThreadingAssembler k-mer list by
    // sorting it before graph construction. Keep the final explicit value in
    // the legacy scalar field, but expose the normalized list to the shared
    // Kokkos graph API and manifest.
    std::sort(options.calling.graph_kmer_sizes.begin(),
              options.calling.graph_kmer_sizes.end());
    options.calling.graph_kmer_sizes.erase(
        std::unique(options.calling.graph_kmer_sizes.begin(),
                    options.calling.graph_kmer_sizes.end()),
        options.calling.graph_kmer_sizes.end());
    if (options.batch_records == 0 || options.threads < 1) throw std::invalid_argument("invalid batch/thread count");
    // The Float32 GKL recurrence is pairwise within the numeric replay bound,
    // but is not yet a full HC-output bitwise replacement on every assembled
    // haplotype set.  Keep the established exact-output HC default strict
    // Float64; an explicit `=false` selects Float32 plus the GKL low-
    // likelihood Float64 fallback.
    options.calling.pairhmm_use_double_precision =
        !options.native_pair_hmm_use_double_precision_explicit ||
        options.native_pair_hmm_use_double_precision;
    resolve_read_filter_options(options);
    return options;
}

void append_batch(fastgatk::io::ReadBatch& all, const fastgatk::io::ReadBatch& batch) {
    if (all.bases.size() > std::numeric_limits<std::uint32_t>::max() - batch.bases.size())
        throw std::runtime_error("RESOURCE_EXHAUSTED: calling batch offsets exceed uint32 range");
    if (!batch.cigar_offsets.empty() && !batch.cigar_layout_valid())
        throw std::runtime_error("BAD_INPUT: malformed CIGAR metadata in HTS batch");
    if (!batch.indel_quality_layout_valid())
        throw std::runtime_error("BAD_INPUT: malformed BI/BD indel quality metadata in HTS batch");
    const auto base_offset = static_cast<std::uint32_t>(all.bases.size());
    const auto record_offset = all.records();
    const auto record_count = batch.records();
    const bool has_cigar = !batch.cigar_offsets.empty();
    const bool has_names = batch.name_offsets.size() == record_count + 1;
    const bool has_read_groups = batch.read_group_offsets.size() == record_count + 1;
    const bool has_oa = batch.original_alignment_offsets.size() == record_count + 1;
    const bool has_xm = batch.mate_contig_offsets.size() == record_count + 1;
    const bool has_oa_presence = batch.original_alignment_present.size() == record_count;
    const bool has_xm_presence = batch.mate_contig_present.size() == record_count;
    if ((!batch.original_alignment_present.empty() && !has_oa_presence) ||
        (!batch.mate_contig_present.empty() && !has_xm_presence))
        throw std::runtime_error("BAD_INPUT: malformed OA/XM presence metadata in HTS batch");
    const bool has_insertion_qualities =
        batch.insertion_quality_offsets.size() == record_count + 1;
    const bool has_deletion_qualities =
        batch.deletion_quality_offsets.size() == record_count + 1;
    const bool has_flow_tp = batch.flow_tp_offsets.size() == record_count + 1;
    const bool has_flow_t0 = batch.flow_t0_offsets.size() == record_count + 1;
    const bool has_flow_order = batch.flow_order_offsets.size() == record_count + 1;
    const bool has_source_input_ids = batch.source_input_ids.size() == record_count;
    if (!batch.source_input_ids.empty() && !has_source_input_ids)
        throw std::runtime_error("BAD_INPUT: malformed input-source metadata in HTS batch");
    if (has_cigar && all.cigar_offsets.empty()) all.cigar_offsets.push_back(0);
    if (has_names && all.name_offsets.empty()) all.name_offsets.push_back(0);
    if (has_read_groups && all.read_group_offsets.empty()) all.read_group_offsets.push_back(0);
    if (has_oa && all.original_alignment_offsets.empty()) all.original_alignment_offsets.push_back(0);
    if (has_xm && all.mate_contig_offsets.empty()) all.mate_contig_offsets.push_back(0);
    if (has_insertion_qualities && all.insertion_quality_offsets.empty())
        all.insertion_quality_offsets.push_back(0);
    if (has_deletion_qualities && all.deletion_quality_offsets.empty())
        all.deletion_quality_offsets.push_back(0);
    if (has_flow_tp && all.flow_tp_offsets.empty()) all.flow_tp_offsets.push_back(0);
    if (has_flow_t0 && all.flow_t0_offsets.empty()) all.flow_t0_offsets.push_back(0);
    if (has_flow_order && all.flow_order_offsets.empty()) all.flow_order_offsets.push_back(0);
    if ((!has_cigar && !all.cigar_offsets.empty()) ||
        (!has_names && !all.name_offsets.empty()) ||
        (!has_read_groups && !all.read_group_offsets.empty()) ||
        (!has_oa && !all.original_alignment_offsets.empty()) ||
        (!has_xm && !all.mate_contig_offsets.empty()) ||
        (!has_oa_presence && !all.original_alignment_present.empty()) ||
        (!has_xm_presence && !all.mate_contig_present.empty()) ||
        (!has_insertion_qualities && !all.insertion_quality_offsets.empty()) ||
        (!has_deletion_qualities && !all.deletion_quality_offsets.empty()) ||
        (!has_flow_tp && !all.flow_tp_offsets.empty()) ||
        (!has_flow_t0 && !all.flow_t0_offsets.empty()) ||
        (!has_flow_order && !all.flow_order_offsets.empty()) ||
        (!has_source_input_ids && !all.source_input_ids.empty()) ||
        (has_source_input_ids && all.records() != 0 && all.source_input_ids.empty()))
        throw std::runtime_error("BAD_INPUT: inconsistent optional metadata across HTS batches");
    if ((has_oa_presence && !all.original_alignment_present.empty() &&
         all.original_alignment_present.size() != all.records()) ||
        (has_xm_presence && !all.mate_contig_present.empty() &&
         all.mate_contig_present.size() != all.records()) ||
        (has_oa_presence && all.records() != 0 && all.original_alignment_present.empty()) ||
        (has_xm_presence && all.records() != 0 && all.mate_contig_present.empty()))
        throw std::runtime_error("BAD_INPUT: inconsistent OA/XM presence metadata across HTS batches");

    auto append_string_span = [](std::vector<std::uint8_t>& payload,
                                 std::vector<std::uint32_t>& offsets,
                                 const std::vector<std::uint8_t>& source,
                                 const std::vector<std::uint32_t>& source_offsets,
                                 std::size_t index) {
        const auto begin = source_offsets[index];
        const auto end = source_offsets[index + 1];
        if (begin > end || end > source.size() ||
            payload.size() > std::numeric_limits<std::uint32_t>::max() - (end - begin))
            throw std::runtime_error("RESOURCE_EXHAUSTED: calling metadata offsets exceed uint32 range");
        payload.insert(payload.end(), source.begin() + begin, source.begin() + end);
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    };

    if (all.cigar_ops.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("RESOURCE_EXHAUSTED: calling CIGAR payload exceeds uint32 range");
    const auto cigar_base = static_cast<std::uint32_t>(all.cigar_ops.size());
    for (std::size_t i = 0; i < record_count; ++i) {
        const auto begin = batch.offsets[i];
        const auto end = batch.offsets[i + 1];
        all.bases.insert(all.bases.end(), batch.bases.begin() + begin, batch.bases.begin() + end);
        all.qualities.insert(all.qualities.end(), batch.qualities.begin() + begin, batch.qualities.begin() + end);
        all.offsets.push_back(base_offset + static_cast<std::uint32_t>(end));
        all.positions.push_back(batch.positions[i]);
        all.tids.push_back(batch.tids[i]);
        all.mapq.push_back(batch.mapq[i]);
        const auto output_record = record_offset + i;
        if (output_record > std::numeric_limits<std::uint32_t>::max())
            throw std::runtime_error("RESOURCE_EXHAUSTED: calling record ordinals exceed uint32 range");
        all.source_records.push_back(static_cast<std::uint32_t>(output_record));
        if (has_source_input_ids) all.source_input_ids.push_back(batch.source_input_ids[i]);
        if (has_cigar) {
            const auto cigar_begin = batch.cigar_offsets[i];
            const auto cigar_end = batch.cigar_offsets[i + 1];
            if (all.cigar_ops.size() > std::numeric_limits<std::uint32_t>::max() - (cigar_end - cigar_begin))
                throw std::runtime_error("RESOURCE_EXHAUSTED: calling CIGAR offsets exceed uint32 range");
            all.cigar_ops.insert(all.cigar_ops.end(), batch.cigar_ops.begin() + cigar_begin,
                                 batch.cigar_ops.begin() + cigar_end);
            all.cigar_offsets.push_back(cigar_base + static_cast<std::uint32_t>(all.cigar_ops.size() - cigar_base));
        }
        if (has_names)
            append_string_span(all.names, all.name_offsets, batch.names, batch.name_offsets, i);
        if (has_read_groups)
            append_string_span(all.read_groups, all.read_group_offsets,
                               batch.read_groups, batch.read_group_offsets, i);
        if (has_oa)
            append_string_span(all.original_alignments, all.original_alignment_offsets,
                               batch.original_alignments, batch.original_alignment_offsets, i);
        if (has_xm)
            append_string_span(all.mate_contigs, all.mate_contig_offsets,
                               batch.mate_contigs, batch.mate_contig_offsets, i);
        if (has_oa_presence) all.original_alignment_present.push_back(batch.original_alignment_present[i]);
        if (has_xm_presence) all.mate_contig_present.push_back(batch.mate_contig_present[i]);
        if (has_insertion_qualities)
            append_string_span(all.insertion_qualities, all.insertion_quality_offsets,
                               batch.insertion_qualities, batch.insertion_quality_offsets, i);
        if (has_deletion_qualities)
            append_string_span(all.deletion_qualities, all.deletion_quality_offsets,
                               batch.deletion_qualities, batch.deletion_quality_offsets, i);
        if (has_flow_tp) {
            const auto begin = batch.flow_tp_offsets[i];
            const auto end = batch.flow_tp_offsets[i + 1];
            all.flow_tp.insert(all.flow_tp.end(), batch.flow_tp.begin() + begin, batch.flow_tp.begin() + end);
            all.flow_tp_offsets.push_back(static_cast<std::uint32_t>(all.flow_tp.size()));
        }
        if (has_flow_t0) {
            const auto begin = batch.flow_t0_offsets[i];
            const auto end = batch.flow_t0_offsets[i + 1];
            all.flow_t0_phred.insert(all.flow_t0_phred.end(), batch.flow_t0_phred.begin() + begin, batch.flow_t0_phred.begin() + end);
            all.flow_t0_offsets.push_back(static_cast<std::uint32_t>(all.flow_t0_phred.size()));
        }
        if (has_flow_order) {
            const auto begin = batch.flow_order_offsets[i];
            const auto end = batch.flow_order_offsets[i + 1];
            all.flow_orders.insert(all.flow_orders.end(), batch.flow_orders.begin() + begin, batch.flow_orders.begin() + end);
            all.flow_order_offsets.push_back(static_cast<std::uint32_t>(all.flow_orders.size()));
        }
        all.flow_max_hmer.push_back(batch.flow_max_hmer.size() == record_count ? batch.flow_max_hmer[i] : 0);
        all.flags.push_back(batch.flags.size() == record_count ? batch.flags[i] : 0);
        all.mate_tids.push_back(batch.mate_tids.size() == record_count ? batch.mate_tids[i] : -1);
        all.mate_positions.push_back(batch.mate_positions.size() == record_count ? batch.mate_positions[i] : -1);
        all.template_lengths.push_back(batch.template_lengths.size() == record_count ? batch.template_lengths[i] : 0);
        if (batch.sample_ids.size() == record_count) all.sample_ids.push_back(batch.sample_ids[i]);
    }
}

// Materialize one record without exposing HTSlib's bam1_t lifetime.  This is
// used only by the opt-in contig streamer when a decode batch straddles a
// contig boundary; the normal path still appends whole batches.  Keeping the
// copy at the shared ReadBatch boundary means CIGAR, BI/BD and flow tags use
// exactly the same validation/layout rules in both modes.
void append_record(fastgatk::io::ReadBatch& destination,
                   const fastgatk::io::ReadBatch& source,
                   const std::size_t index,
                   const std::uint32_t source_input_id) {
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
    one.source_input_ids = {source_input_id};
    append_batch(destination, one);
}

std::uint64_t stream_signature_mix(std::uint64_t hash, std::uint64_t value) {
    hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    return hash;
}

void merge_stream_result(fastgatk::calling::Result& destination,
                         fastgatk::calling::Result&& source,
                         bool& initialized) {
    if (!initialized) {
        destination = std::move(source);
        // The per-contig matrices and graph payloads have already been
        // materialized into the output body.  Drop them before the next
        // contig is decoded; scalar counters below remain the manifest source
        // of truth for the aggregate stream result.
        destination.reference_blocks.clear();
        destination.graph_haplotype_sequences.clear();
        destination.graph_haplotype_path_tids.clear();
        destination.graph_haplotype_path_starts.clear();
        destination.graph_haplotype_path_ends.clear();
        destination.graph_haplotype_path_support.clear();
        destination.graph_haplotype_path_scores.clear();
        destination.graph_haplotype_path_has_non_reference_edge.clear();
        destination.graph_haplotype_path_alt_read_starts.clear();
        destination.graph_haplotype_path_alt_read_ends.clear();
        destination.graph_haplotype_path_cigars.clear();
        destination.graph_haplotype_path_alignment_scores.clear();
        destination.graph_haplotype_path_alignment_offsets.clear();
        destination.graph_haplotype_event_maps.clear();
        destination.allele_read_likelihoods.clear();
        destination.reference_read_likelihoods.clear();
        destination.spanning_deletion_read_likelihoods.clear();
        destination.multiallelic_depths.clear();
        destination.non_ref_read_likelihoods.clear();
        destination.candidate_prior_het.clear();
        destination.candidate_prior_hom_alt.clear();
        destination.candidates.clear();
        destination.likelihoods.clear();
        initialized = true;
        return;
    }

#define FASTGATK_SUM_SIZE(field) destination.field += source.field
#define FASTGATK_SUM_DOUBLE(field) destination.field += source.field
    FASTGATK_SUM_SIZE(reads);
    FASTGATK_SUM_SIZE(observations);
    FASTGATK_SUM_SIZE(projected_observations);
    FASTGATK_SUM_SIZE(loci);
    FASTGATK_SUM_SIZE(candidate_sites);
    FASTGATK_SUM_SIZE(candidate_softclip_suppressed);
    FASTGATK_SUM_SIZE(candidate_fragment_suppressed);
    FASTGATK_SUM_SIZE(candidate_low_support_suppressed);
    FASTGATK_SUM_SIZE(rcm_realigned_observations);
    FASTGATK_SUM_SIZE(rcm_realignment_fallback_observations);
    destination.rcm_haplotype_realignment_used =
        destination.rcm_haplotype_realignment_used || source.rcm_haplotype_realignment_used;
    FASTGATK_SUM_SIZE(multiallelic_candidates_pruned);
    FASTGATK_SUM_SIZE(softclip_candidates);
    FASTGATK_SUM_SIZE(compound_cigar_candidates);
    FASTGATK_SUM_SIZE(assembly_unassigned_candidates);
    FASTGATK_SUM_SIZE(assembly_cross_region_candidates);
    FASTGATK_SUM_SIZE(assembly_region_union_count);
    FASTGATK_SUM_SIZE(assembly_region_union_merges);
    FASTGATK_SUM_SIZE(assembly_cross_region_rescued);
    FASTGATK_SUM_SIZE(graph_variant_candidates);
    FASTGATK_SUM_SIZE(graph_snp_candidates);
    FASTGATK_SUM_SIZE(graph_mnp_candidates);
    FASTGATK_SUM_SIZE(graph_deletion_candidates);
    FASTGATK_SUM_SIZE(variant_calls);
    FASTGATK_SUM_DOUBLE(assembly_prepare_seconds);
    FASTGATK_SUM_DOUBLE(assembly_seconds);
    FASTGATK_SUM_DOUBLE(likelihood_prepare_seconds);
    FASTGATK_SUM_DOUBLE(likelihood_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_prepare_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_seconds);
    FASTGATK_SUM_SIZE(pairhmm_pairs);
    FASTGATK_SUM_SIZE(pairhmm_request_links);
    FASTGATK_SUM_SIZE(pairhmm_haplotypes);
    FASTGATK_SUM_SIZE(pairhmm_graph_haplotypes);
    FASTGATK_SUM_SIZE(pairhmm_graph_snp_posterior_pairs);
    FASTGATK_SUM_SIZE(pairhmm_graph_haplotypes_considered);
    FASTGATK_SUM_SIZE(pairhmm_haplotypes_pruned);
    FASTGATK_SUM_SIZE(pairhmm_graph_haplotypes_pruned);
    FASTGATK_SUM_SIZE(pairhmm_haplotype_combination_blocks);
    FASTGATK_SUM_SIZE(pairhmm_assembly_region_groups);
    FASTGATK_SUM_SIZE(pairhmm_unassigned_candidates);
    FASTGATK_SUM_SIZE(pairhmm_marginalized_candidates);
    FASTGATK_SUM_SIZE(pairhmm_cached_shapes);
    FASTGATK_SUM_SIZE(pairhmm_cache_hits);
    FASTGATK_SUM_DOUBLE(pairhmm_normalization_prepare_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_normalization_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_marginalization_prepare_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_marginalization_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_uncertainty_prepare_seconds);
    FASTGATK_SUM_DOUBLE(pairhmm_uncertainty_seconds);
    FASTGATK_SUM_SIZE(pairhmm_flow_reads);
    FASTGATK_SUM_SIZE(pairhmm_flow_haplotypes);
    FASTGATK_SUM_SIZE(pairhmm_flow_haplotypes_clipped);
    FASTGATK_SUM_SIZE(pairhmm_flow_haplotypes_collapsed);
    FASTGATK_SUM_SIZE(pairhmm_flow_haplotypes_uncollapsed);
    FASTGATK_SUM_SIZE(pairhmm_flow_haplotype_remaps);
    FASTGATK_SUM_SIZE(pairhmm_flow_identical_haplotype_groups);
    FASTGATK_SUM_SIZE(pairhmm_flow_reads_clipped);
    FASTGATK_SUM_SIZE(pairhmm_flow_clipping_fallbacks);
    FASTGATK_SUM_SIZE(pairhmm_insertion_quality_reads);
    FASTGATK_SUM_SIZE(pairhmm_deletion_quality_reads);
    FASTGATK_SUM_SIZE(pairhmm_pcr_adjusted_positions);
    FASTGATK_SUM_SIZE(pairhmm_reads_clipped);
    FASTGATK_SUM_SIZE(pairhmm_reads_dropped_after_clipping);
    FASTGATK_SUM_SIZE(pairhmm_reads_disqualified);
    // This is a configuration/telemetry scalar, not a per-contig work
    // counter.  Summing it would manufacture a value that no invocation
    // used once more than one contig is streamed; retain the strictest
    // threshold observed across chunks instead.
    destination.pairhmm_read_disqualification_threshold = std::max(
        destination.pairhmm_read_disqualification_threshold,
        source.pairhmm_read_disqualification_threshold);
    FASTGATK_SUM_DOUBLE(sw_prepare_seconds);
    FASTGATK_SUM_DOUBLE(sw_seconds);
    FASTGATK_SUM_SIZE(sw_pairs);
    FASTGATK_SUM_SIZE(sw_simd_groups);
    FASTGATK_SUM_SIZE(sw_traceback_pairs);
    FASTGATK_SUM_SIZE(sw_indel_alignments);
    FASTGATK_SUM_DOUBLE(sw_traceback_seconds);
    FASTGATK_SUM_SIZE(read_haplotype_cigar_pairs);
    FASTGATK_SUM_SIZE(read_haplotype_cigar_valid_pairs);
    FASTGATK_SUM_SIZE(read_haplotype_cigar_indel_pairs);
    FASTGATK_SUM_SIZE(read_haplotype_cigar_filtered_pairs);
    FASTGATK_SUM_SIZE(read_haplotype_softclip_filtered_pairs);
    FASTGATK_SUM_SIZE(read_haplotype_uncertain_reads);
    FASTGATK_SUM_SIZE(read_haplotype_informative_reads);
    FASTGATK_SUM_DOUBLE(read_haplotype_margin_sum);
    FASTGATK_SUM_SIZE(graph_nodes);
    FASTGATK_SUM_SIZE(graph_edges);
    FASTGATK_SUM_SIZE(graph_branching_nodes);
    FASTGATK_SUM_SIZE(graph_reference_nodes);
    FASTGATK_SUM_SIZE(graph_reference_edges);
    FASTGATK_SUM_SIZE(graph_reference_connected_nodes);
    FASTGATK_SUM_SIZE(graph_dangling_nodes);
    FASTGATK_SUM_SIZE(graph_pruned_nodes);
    FASTGATK_SUM_SIZE(graph_dangling_branch_paths);
    FASTGATK_SUM_SIZE(graph_dangling_branch_bases);
    FASTGATK_SUM_SIZE(graph_dangling_recovered_paths);
    FASTGATK_SUM_SIZE(graph_dangling_recovered_bases);
    FASTGATK_SUM_SIZE(graph_artificial_haplotype_recovery_paths);
    FASTGATK_SUM_SIZE(graph_artificial_haplotype_recovery_bases);
    FASTGATK_SUM_SIZE(graph_seqgraph_nodes);
    FASTGATK_SUM_SIZE(graph_seqgraph_edges);
    FASTGATK_SUM_SIZE(graph_seqgraph_linear_chain_merges);
    FASTGATK_SUM_SIZE(graph_seqgraph_diamond_merges);
    FASTGATK_SUM_SIZE(graph_seqgraph_tail_merges);
    FASTGATK_SUM_SIZE(graph_seqgraph_suffix_splits);
    FASTGATK_SUM_SIZE(graph_seqgraph_suffix_merges);
    FASTGATK_SUM_SIZE(graph_adaptive_pruned_nodes);
    FASTGATK_SUM_SIZE(graph_non_unique_kmers);
    FASTGATK_SUM_SIZE(graph_reference_non_unique_kmers);
    FASTGATK_SUM_SIZE(graph_reference_paths);
    FASTGATK_SUM_SIZE(graph_haplotype_paths);
    FASTGATK_SUM_SIZE(graph_haplotype_sequence_count);
    FASTGATK_SUM_SIZE(graph_haplotype_provenance_count);
    FASTGATK_SUM_SIZE(graph_haplotype_event_count);
    FASTGATK_SUM_SIZE(graph_haplotype_sw_simd_groups);
    FASTGATK_SUM_DOUBLE(graph_haplotype_sw_prepare_seconds);
    FASTGATK_SUM_DOUBLE(graph_haplotype_sw_seconds);
    FASTGATK_SUM_SIZE(activity_loci);
    FASTGATK_SUM_SIZE(active_loci);
    FASTGATK_SUM_SIZE(assembly_regions);
    // `--assembly-region-out` is a traversal trace, not a scalar counter.
    // Preserve each contig's raw profile segments in coordinate order so the
    // streaming path emits the same complete IGV surface as aggregate mode.
    destination.activity_profile_regions.insert(
        destination.activity_profile_regions.end(),
        source.activity_profile_regions.begin(), source.activity_profile_regions.end());
    destination.activity_active_loci.insert(
        destination.activity_active_loci.end(),
        source.activity_active_loci.begin(), source.activity_active_loci.end());
    destination.activity_filter_size = std::max(
        destination.activity_filter_size, source.activity_filter_size);
    destination.activity_effective_max_probability_propagation_distance = std::max(
        destination.activity_effective_max_probability_propagation_distance,
        source.activity_effective_max_probability_propagation_distance);
    FASTGATK_SUM_DOUBLE(activity_prepare_seconds);
    FASTGATK_SUM_DOUBLE(activity_seconds);
    FASTGATK_SUM_SIZE(filtered_reads);
    FASTGATK_SUM_SIZE(downsampled_reads);
    FASTGATK_SUM_DOUBLE(read_filter_prepare_seconds);
    FASTGATK_SUM_DOUBLE(read_filter_seconds);
    FASTGATK_SUM_SIZE(overlapping_pairs);
    FASTGATK_SUM_SIZE(overlapping_bases);
    FASTGATK_SUM_SIZE(overlapping_conflicting_bases);
    FASTGATK_SUM_SIZE(overlapping_quality_caps);
    FASTGATK_SUM_SIZE(error_correction_solid_kmers);
    FASTGATK_SUM_SIZE(error_correction_corrected_kmers);
    FASTGATK_SUM_SIZE(error_correction_uncorrectable_kmers);
    FASTGATK_SUM_SIZE(error_correction_pileup_loci);
    FASTGATK_SUM_SIZE(error_correction_pileup_skipped_indel_adjacent_bases);
    FASTGATK_SUM_SIZE(error_correction_corrected_reads);
    FASTGATK_SUM_SIZE(error_correction_corrected_bases);
    FASTGATK_SUM_DOUBLE(error_correction_prepare_seconds);
    FASTGATK_SUM_DOUBLE(error_correction_seconds);
    FASTGATK_SUM_DOUBLE(reference_confidence_prepare_seconds);
    FASTGATK_SUM_DOUBLE(reference_confidence_seconds);
    FASTGATK_SUM_SIZE(reference_block_count);
    FASTGATK_SUM_SIZE(genotype_prior_kernel_calls);
    FASTGATK_SUM_DOUBLE(genotype_prior_prepare_seconds);
    FASTGATK_SUM_DOUBLE(genotype_prior_seconds);
    FASTGATK_SUM_DOUBLE(genotyping_seconds);
#undef FASTGATK_SUM_SIZE
#undef FASTGATK_SUM_DOUBLE

    destination.signature = stream_signature_mix(destination.signature, source.signature);
    destination.sw_cigar_signature = stream_signature_mix(
        destination.sw_cigar_signature, source.sw_cigar_signature);
    destination.read_haplotype_uncertainty_signature = stream_signature_mix(
        destination.read_haplotype_uncertainty_signature,
        source.read_haplotype_uncertainty_signature);
    destination.graph_haplotype_cigar_signature = stream_signature_mix(
        destination.graph_haplotype_cigar_signature, source.graph_haplotype_cigar_signature);
    destination.sw_simd_width = std::max(destination.sw_simd_width, source.sw_simd_width);
    destination.graph_haplotype_sequence_max_length = std::max(
        destination.graph_haplotype_sequence_max_length,
        source.graph_haplotype_sequence_max_length);
    destination.graph_haplotype_sw_simd_width = std::max(
        destination.graph_haplotype_sw_simd_width,
        source.graph_haplotype_sw_simd_width);
    destination.graph_kmer_size_selected = destination.graph_kmer_size_selected == 0
        ? source.graph_kmer_size_selected
        : (destination.graph_kmer_size_selected == source.graph_kmer_size_selected
            ? destination.graph_kmer_size_selected : 0);
    destination.graph_kmer_iterations += source.graph_kmer_iterations;
    destination.pairhmm_haplotype_combination_chunked =
        destination.pairhmm_haplotype_combination_chunked || source.pairhmm_haplotype_combination_chunked;
    destination.pairhmm_used = destination.pairhmm_used || source.pairhmm_used;
    destination.sw_used = destination.sw_used || source.sw_used;
    destination.graph_haplotype_sw_used =
        destination.graph_haplotype_sw_used || source.graph_haplotype_sw_used;
    if (destination.graph_haplotype_sw_execution_space.empty())
        destination.graph_haplotype_sw_execution_space = source.graph_haplotype_sw_execution_space;
    else if (!source.graph_haplotype_sw_execution_space.empty() &&
             destination.graph_haplotype_sw_execution_space != source.graph_haplotype_sw_execution_space)
        destination.graph_haplotype_sw_execution_space = "mixed";
    destination.graph_used = destination.graph_used || source.graph_used;
    destination.activity_used = destination.activity_used || source.activity_used;
    destination.activity_reference_aware = destination.activity_reference_aware && source.activity_reference_aware;
    destination.read_filter_used = destination.read_filter_used || source.read_filter_used;
    destination.read_filter_applied = destination.read_filter_applied && source.read_filter_applied;
    destination.overlapping_quality_correction_used =
        destination.overlapping_quality_correction_used || source.overlapping_quality_correction_used;
    destination.overlapping_quality_correction_metadata_available =
        destination.overlapping_quality_correction_metadata_available &&
        source.overlapping_quality_correction_metadata_available;
    destination.error_correction_used = destination.error_correction_used || source.error_correction_used;
    destination.graph_reference_kmer_rejected =
        destination.graph_reference_kmer_rejected || source.graph_reference_kmer_rejected;
    destination.graph_has_non_reference_cycles =
        destination.graph_has_non_reference_cycles || source.graph_has_non_reference_cycles;
    destination.non_ref_likelihoods_present =
        destination.non_ref_likelihoods_present || source.non_ref_likelihoods_present;
    destination.reference_block_lookup_indexed =
        destination.reference_block_lookup_indexed && source.reference_block_lookup_indexed;

    const auto merge_text = [](std::string& target, const std::string& value,
                               const char* mixed) {
        if (target.empty()) target = value;
        else if (target != value) target = mixed;
    };
    merge_text(destination.pairhmm_error_model, source.pairhmm_error_model, "mixed");
    merge_text(destination.pairhmm_precision, source.pairhmm_precision, "mixed");
    merge_text(destination.pairhmm_pcr_indel_model, source.pairhmm_pcr_indel_model, "mixed");
    destination.pairhmm_pcr_error_rate_factor =
        destination.pairhmm_pcr_error_rate_factor == 0.0
            ? source.pairhmm_pcr_error_rate_factor
            : (source.pairhmm_pcr_error_rate_factor == 0.0 ||
               destination.pairhmm_pcr_error_rate_factor == source.pairhmm_pcr_error_rate_factor
                   ? destination.pairhmm_pcr_error_rate_factor : -1.0);
    if (destination.pairhmm_default_indel_quality == 0)
        destination.pairhmm_default_indel_quality = source.pairhmm_default_indel_quality;
    else if (source.pairhmm_default_indel_quality != 0 &&
             destination.pairhmm_default_indel_quality != source.pairhmm_default_indel_quality)
        destination.pairhmm_default_indel_quality = 0;
    if (destination.pairhmm_default_indel_terminal_quality == 0)
        destination.pairhmm_default_indel_terminal_quality = source.pairhmm_default_indel_terminal_quality;
    else if (source.pairhmm_default_indel_terminal_quality != 0 &&
             destination.pairhmm_default_indel_terminal_quality != source.pairhmm_default_indel_terminal_quality)
        destination.pairhmm_default_indel_terminal_quality = 0;
    merge_text(destination.pairhmm_skip_reason, source.pairhmm_skip_reason, "mixed");
    merge_text(destination.pairhmm_execution_space, source.pairhmm_execution_space, "mixed");
    merge_text(destination.pairhmm_normalization_execution_space,
               source.pairhmm_normalization_execution_space, "mixed");
    merge_text(destination.pairhmm_marginalization_execution_space,
               source.pairhmm_marginalization_execution_space, "mixed");
    merge_text(destination.pairhmm_uncertainty_execution_space,
               source.pairhmm_uncertainty_execution_space, "mixed");
    merge_text(destination.sw_execution_space, source.sw_execution_space, "mixed");
    merge_text(destination.graph_execution_space, source.graph_execution_space, "mixed");
    merge_text(destination.activity_execution_space, source.activity_execution_space, "mixed");
    merge_text(destination.read_filter_execution_space, source.read_filter_execution_space, "mixed");
    merge_text(destination.error_correction_execution_space, source.error_correction_execution_space, "mixed");
    merge_text(destination.error_correction_mode, source.error_correction_mode, "mixed");
    merge_text(destination.reference_confidence_execution_space,
               source.reference_confidence_execution_space, "mixed");
    merge_text(destination.genotype_prior_execution_space,
               source.genotype_prior_execution_space, "mixed");
    destination.calls.insert(destination.calls.end(),
                             std::make_move_iterator(source.calls.begin()),
                             std::make_move_iterator(source.calls.end()));
    destination.gvcf_hom_ref_calls.insert(
        destination.gvcf_hom_ref_calls.end(),
        std::make_move_iterator(source.gvcf_hom_ref_calls.begin()),
        std::make_move_iterator(source.gvcf_hom_ref_calls.end()));
    destination.multiallelic_depths.insert(
        destination.multiallelic_depths.end(),
        std::make_move_iterator(source.multiallelic_depths.begin()),
        std::make_move_iterator(source.multiallelic_depths.end()));
}

std::vector<std::string> load_reference(const std::string& path,
                                        const fastgatk::io::HeaderSummary& header) {
    std::vector<std::string> sequences(header.contigs.size());
    if (path.empty()) return sequences;
#if FASTGATK_HAS_HTSLIB
    // Match GATK's indexed ReferenceFileSource for both plain and bgzip
    // FASTA.  In particular, std::ifstream cannot parse a .fa.gz payload and
    // would silently route a reference-backed HaplotypeCaller run through the
    // no-reference PairHMM fallback.  Limit materialization to contigs in the
    // input header so targeted inputs do not require loading an entire genome.
    // Do not let faidx_load create an index as an incidental side effect for
    // short unindexed smoke FASTAs.  Besides mutating the input directory,
    // that would make a later invocation validate against a dictionary the
    // caller did not supply.  Production/bgzipped references provide the
    // FAI explicitly and take the random-access path below.
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
                // Dictionary validation is performed separately unless the
                // caller explicitly disabled it.  In the latter case GATK
                // can legally traverse a bounded interval on the subset of
                // contigs that *is* present (the tetra-ploid oracle is such
                // a fixture).  Leave unrelated header contigs empty rather
                // than rejecting the whole run before the selected interval
                // reaches the reference-backed pipeline.
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
    std::map<std::string, std::size_t> contig_ids;
    for (std::size_t i = 0; i < header.contigs.size(); ++i) contig_ids.emplace(header.contigs[i], i);
    std::string line;
    std::string current;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            std::istringstream name(line.substr(1));
            name >> current;
            continue;
        }
        if (current.empty()) continue;
        auto id = contig_ids.find(current);
        if (id != contig_ids.end()) sequences[id->second] += line;
    }
    return sequences;
}

void set_calling_intervals(Options& options, const fastgatk::io::HtsReader& reader) {
    options.calling.intervals = reader.intervals();
    options.calling.exclusion_intervals = reader.exclusion_intervals();
    if (options.calling.intervals.size() == 1) {
        const auto& interval = options.calling.intervals.front();
        if (interval.start < 0 || interval.end > std::numeric_limits<std::int32_t>::max())
            throw std::invalid_argument("BAD_INPUT: interval exceeds native int32 coordinate range");
        options.calling.interval_tid = interval.tid;
        options.calling.interval_start = static_cast<std::int32_t>(interval.start);
        options.calling.interval_end = static_cast<std::int32_t>(interval.end);
    }
}

#if FASTGATK_HAS_HTSLIB
bool feature_record_is_filtered(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && std::string(name) != "PASS" && std::string(name) != ".")
            return true;
    }
    return false;
}

bool concrete_feature_allele(const char* allele) {
    if (allele == nullptr || allele[0] == '\0') return false;
    const std::string value(allele);
    return value != "." && value != "*" && value.front() != '<' && value.back() != '>';
}

std::vector<fastgatk::calling::ForcedAllele> load_forced_alleles(
    const std::string& path, const fastgatk::io::HeaderSummary& read_header,
    const bool include_filtered) {
    std::vector<fastgatk::calling::ForcedAllele> result;
    if (path.empty()) return result;
    htsFile* input = bcf_open(path.c_str(), "r");
    if (input == nullptr)
        throw std::runtime_error("BAD_INPUT: cannot open --alleles VCF/BCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    bcf1_t* record = bcf_init();
    if (header == nullptr || record == nullptr) {
        if (record != nullptr) bcf_destroy(record);
        if (header != nullptr) bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read --alleles VCF/BCF header: " + path);
    }
    try {
        std::map<std::string, std::int32_t> contig_ids;
        for (std::size_t index = 0; index < read_header.contigs.size(); ++index) {
            if (index > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
                throw std::runtime_error("BAD_INPUT: input sequence dictionary exceeds native coordinate range");
            contig_ids.emplace(read_header.contigs[index], static_cast<std::int32_t>(index));
        }
        std::set<std::tuple<std::int32_t, std::int32_t, std::string, std::string>> seen;
        while (true) {
            const auto status = bcf_read(input, header, record);
            if (status < 0) {
                if (status != -1)
                    throw std::runtime_error("BAD_INPUT: cannot read --alleles VCF/BCF record: " + path);
                break;
            }
            bcf_unpack(record, BCF_UN_STR | BCF_UN_FLT);
            const auto* contig = record->rid < 0 ? nullptr : bcf_hdr_id2name(header, record->rid);
            if (contig == nullptr || record->pos < 0 ||
                record->pos > std::numeric_limits<std::int32_t>::max() ||
                record->d.allele == nullptr ||
                record->n_allele < 2) {
                bcf_clear(record);
                continue;
            }
            const auto target = contig_ids.find(contig);
            if (target == contig_ids.end())
                throw std::runtime_error(
                    "BAD_INPUT: --alleles VCF contig is absent from input reads: " + std::string(contig));
            if (!include_filtered && feature_record_is_filtered(header, record)) {
                bcf_clear(record);
                continue;
            }
            if (!concrete_feature_allele(record->d.allele[0])) {
                bcf_clear(record);
                continue;
            }
            const std::string reference(record->d.allele[0]);
            for (int alternate_index = 1; alternate_index < record->n_allele; ++alternate_index) {
                const auto* alternate = record->d.allele[alternate_index];
                if (!concrete_feature_allele(alternate)) continue;
                const std::string alternate_text(alternate);
                if (reference == alternate_text) continue;
                const auto key = std::make_tuple(target->second, record->pos,
                                                 reference, alternate_text);
                if (!seen.insert(key).second) continue;
                fastgatk::calling::ForcedAllele forced;
                forced.tid = target->second;
                forced.position = static_cast<std::int32_t>(record->pos);
                forced.reference = reference;
                forced.alternate = alternate_text;
                result.push_back(std::move(forced));
            }
            bcf_clear(record);
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(input);
        return result;
    } catch (...) {
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(input);
        throw;
    }
}
#else
std::vector<fastgatk::calling::ForcedAllele> load_forced_alleles(
    const std::string& path, const fastgatk::io::HeaderSummary&, const bool) {
    if (!path.empty())
        throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for HaplotypeCaller --alleles");
    return {};
}
#endif

std::string summary_json(const Options& options, const fastgatk::io::HtsReader& reader,
                 const fastgatk::calling::Result& result,
                 const fastgatk::runtime::ResourceSnapshot& resources) {
    std::ostringstream out;
    const auto assembly_tags = static_cast<std::size_t>(std::count_if(
        reader.header().contig_assemblies.begin(), reader.header().contig_assemblies.end(),
        [](const std::string& value) { return !value.empty(); }));
    const bool interval_subset = !options.calling.intervals.empty() ||
        (options.calling.interval_tid >= 0 && options.calling.interval_end > options.calling.interval_start);
    // Mirror AssemblyResultSet rather than the raw k-best graph traversal in
    // externally visible haplotype provenance.  GATK drops paths whose
    // SOFTCLIP CIGAR cannot span the complete reference haplotype after the
    // INDEL retry confirms a boundary-resolution failure.  The raw graph
    // remains in Result for candidate/PairHMM projection, but it must not be
    // serialized as a sequence-only, reference-connected haplotype.
    std::vector<std::size_t> retained_graph_paths;
    retained_graph_paths.reserve(result.graph_haplotype_sequences.size());
    for (std::size_t path = 0; path < result.graph_haplotype_sequences.size(); ++path) {
        if (path < result.graph_haplotype_path_cigars.size() &&
            !result.graph_haplotype_path_cigars[path].empty())
            retained_graph_paths.push_back(path);
    }
    std::size_t retained_graph_haplotype_max_length = 0;
    for (const auto path : retained_graph_paths)
        retained_graph_haplotype_max_length = std::max(
            retained_graph_haplotype_max_length,
            result.graph_haplotype_sequences[path].size());
    out << "{\"tool\":\"HaplotypeCallerNativeCall\",\"status\":\"prototype\""
        << ",\"pipeline\":\"assembly-likelihood-genotyping\""
        << ",\"reads\":" << result.reads
        << ",\"htslib_backend\":\"" << reader.backend_description() << "\""
        << ",\"sample_name\":\"" << json_escape(
            options.sample_name.empty() ? (reader.header().samples.empty()
                ? std::string("FASTGATK") : reader.header().samples.front())
                : options.sample_name) << "\""
        << ",\"sample_name_selected\":" << (!options.sample_name.empty() ? "true" : "false")
        << ",\"input_count\":" << options.inputs.size()
        << ",\"inputs\":[";
    for (std::size_t input_index = 0; input_index < options.inputs.size(); ++input_index) {
        if (input_index != 0) out << ',';
        out << "\"" << json_escape(options.inputs[input_index]) << "\"";
    }
    out << "]"
        << ",\"interval_set_rule\":\""
        << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
            ? "INTERSECTION" : "UNION") << "\""
        << ",\"create_output_variant_index\":"
        << (options.create_output_variant_index ? "true" : "false")
        << ",\"add_output_vcf_command_line\":"
        << (options.add_output_vcf_command_line ? "true" : "false")
        << ",\"sites_only_vcf_output\":"
        << (options.sites_only_vcf_output ? "true" : "false")
        << ",\"floor_blocks\":" << (options.floor_blocks ? "true" : "false")
        << ",\"alleles_feature\":" << (!options.alleles.empty() ? "true" : "false")
        << ",\"forced_alleles\":" << options.calling.forced_alleles.size()
        << ",\"force_call_filtered_alleles\":"
        << (options.force_call_filtered_alleles ? "true" : "false")
        << ",\"native_pair_hmm_use_double_precision\":"
        << (options.native_pair_hmm_use_double_precision ? "true" : "false")
        << ",\"pairhmm_effective_precision\":\""
        << json_escape(result.pairhmm_precision) << "\""
        << ",\"disable_sequence_dictionary_validation\":"
        << (options.disable_sequence_dictionary_validation ? "true" : "false")
        << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
        << ",\"resources\":" << resources.to_json()
        << ",\"requested_threads\":" << options.threads
        << ",\"effective_threads\":" << resources.effective_threads(static_cast<std::size_t>(options.threads))
        << ",\"requested_batch_records\":" << options.requested_batch_records
        << ",\"effective_batch_records\":" << (options.effective_batch_records == 0 ? options.batch_records : options.effective_batch_records)
        << ",\"adaptive_batch_reductions\":" << options.adaptive_batch_reductions
        << ",\"stream_by_contig\":" << (options.stream_by_contig ? "true" : "false")
        << ",\"stream_by_region\":" << (options.stream_by_region ? "true" : "false")
        << ",\"stream_region_size\":" << options.stream_region_size
        << ",\"streamed_contigs\":" << options.streamed_contigs
        << ",\"streamed_regions\":" << options.streamed_regions
        << ",\"streamed_region_splits\":" << options.streamed_region_splits
        << ",\"streamed_peak_host_bytes\":" << options.streamed_peak_host_bytes
        << ",\"stream_indexed\":" << (options.stream_indexed ? "true" : "false")
        << ",\"interval_exclusions\":" << options.calling.exclusion_intervals.size()
        << ",\"interval_padding\":" << options.interval_padding
        << ",\"interval_exclusion_padding\":" << options.exclusion_padding
        << ",\"pipeline_lifecycle\":\""
        << (options.pipeline_used
            ? "Host decode->bounded queue->Kokkos compute->encode->sink" : "not-used") << "\""
        << ",\"pipeline_decoded_items\":" << options.pipeline_decoded_items
        << ",\"pipeline_computed_items\":" << options.pipeline_computed_items
        << ",\"pipeline_encoded_items\":" << options.pipeline_encoded_items
        << ",\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes
        << ",\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes
        << ",\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes
        << ",\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes
        << ",\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes
        << ",\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes
        << ",\"contig_assembly_tags\":" << assembly_tags
        << ",\"pairhmm_used\":" << (result.pairhmm_used ? "true" : "false")
        << ",\"pairhmm_error_model\":\"" << json_escape(result.pairhmm_error_model) << "\""
        << ",\"pairhmm_default_indel_quality\":" << static_cast<unsigned>(result.pairhmm_default_indel_quality)
        << ",\"pairhmm_default_indel_terminal_quality\":" << static_cast<unsigned>(result.pairhmm_default_indel_terminal_quality)
        << ",\"pairhmm_pcr_indel_model\":\"" << json_escape(result.pairhmm_pcr_indel_model) << "\""
        << ",\"pairhmm_pcr_error_rate_factor\":" << result.pairhmm_pcr_error_rate_factor
        << ",\"phred_scaled_global_read_mismapping_rate\":"
        << options.calling.phred_scaled_global_read_mismapping_rate
        << ",\"pairhmm_pcr_adjusted_positions\":" << result.pairhmm_pcr_adjusted_positions
        << ",\"pairhmm_flow_reads\":" << result.pairhmm_flow_reads
        << ",\"pairhmm_flow_haplotypes\":" << result.pairhmm_flow_haplotypes
        << ",\"pairhmm_flow_haplotypes_clipped\":" << result.pairhmm_flow_haplotypes_clipped
        << ",\"pairhmm_flow_haplotypes_collapsed\":" << result.pairhmm_flow_haplotypes_collapsed
        << ",\"pairhmm_flow_haplotypes_uncollapsed\":" << result.pairhmm_flow_haplotypes_uncollapsed
        << ",\"pairhmm_flow_haplotype_remaps\":" << result.pairhmm_flow_haplotype_remaps
        << ",\"pairhmm_flow_identical_haplotype_groups\":" << result.pairhmm_flow_identical_haplotype_groups
        << ",\"pairhmm_flow_reads_clipped\":" << result.pairhmm_flow_reads_clipped
        << ",\"pairhmm_flow_clipping_fallbacks\":" << result.pairhmm_flow_clipping_fallbacks
        << ",\"pairhmm_insertion_quality_reads\":" << result.pairhmm_insertion_quality_reads
        << ",\"pairhmm_deletion_quality_reads\":" << result.pairhmm_deletion_quality_reads
        << ",\"pairhmm_reads_clipped\":" << result.pairhmm_reads_clipped
        << ",\"pairhmm_reads_dropped_after_clipping\":"
        << result.pairhmm_reads_dropped_after_clipping
        << ",\"pairhmm_reads_disqualified\":"
        << result.pairhmm_reads_disqualified
        << ",\"pairhmm_read_disqualification_threshold\":"
        << result.pairhmm_read_disqualification_threshold
        << ",\"pairhmm_skip_reason\":\"" << json_escape(result.pairhmm_skip_reason) << "\""
        << ",\"pairhmm_pruning_policy\":\"support-floor-with-low-support-siblings-metadata-only\""
        << ",\"pairhmm_combination_support_floor\":"
        << std::max<std::size_t>(3, options.calling.min_alt_support)
        << ",\"pairhmm_pairs\":" << result.pairhmm_pairs
        << ",\"pairhmm_request_links\":" << result.pairhmm_request_links
        << ",\"pairhmm_haplotypes\":" << result.pairhmm_haplotypes
        << ",\"pairhmm_graph_haplotypes\":" << result.pairhmm_graph_haplotypes
        << ",\"pairhmm_graph_snp_posterior_pairs\":"
        << result.pairhmm_graph_snp_posterior_pairs
        << ",\"pairhmm_graph_haplotypes_considered\":" << result.pairhmm_graph_haplotypes_considered
        << ",\"pairhmm_haplotypes_pruned\":" << result.pairhmm_haplotypes_pruned
        << ",\"pairhmm_graph_haplotypes_pruned\":" << result.pairhmm_graph_haplotypes_pruned
        << ",\"pairhmm_haplotype_combination_blocks\":"
        << result.pairhmm_haplotype_combination_blocks
        << ",\"pairhmm_haplotype_combination_chunked\":"
        << (result.pairhmm_haplotype_combination_chunked ? "true" : "false")
        << ",\"pairhmm_assembly_region_groups\":" << result.pairhmm_assembly_region_groups
        << ",\"pairhmm_unassigned_candidates\":" << result.pairhmm_unassigned_candidates
        << ",\"pairhmm_assembly_region_partitioned\":"
        << (result.pairhmm_assembly_region_partitioned ? "true" : "false")
        << ",\"pairhmm_marginalized_candidates\":" << result.pairhmm_marginalized_candidates
        << ",\"pairhmm_cached_shapes\":" << result.pairhmm_cached_shapes
        << ",\"pairhmm_cache_hits\":" << result.pairhmm_cache_hits
        << ",\"pairhmm_normalization_execution_space\":\""
        << json_escape(result.pairhmm_normalization_execution_space) << "\""
        << ",\"pairhmm_normalization_prepare_seconds\":"
        << result.pairhmm_normalization_prepare_seconds
        << ",\"pairhmm_normalization_seconds\":"
        << result.pairhmm_normalization_seconds
        << ",\"pairhmm_marginalization_execution_space\":\""
        << json_escape(result.pairhmm_marginalization_execution_space) << "\""
        << ",\"pairhmm_marginalization_prepare_seconds\":"
        << result.pairhmm_marginalization_prepare_seconds
        << ",\"pairhmm_marginalization_seconds\":"
        << result.pairhmm_marginalization_seconds
        << ",\"pairhmm_uncertainty_execution_space\":\""
        << json_escape(result.pairhmm_uncertainty_execution_space) << "\""
        << ",\"pairhmm_uncertainty_prepare_seconds\":"
        << result.pairhmm_uncertainty_prepare_seconds
        << ",\"pairhmm_uncertainty_seconds\":"
        << result.pairhmm_uncertainty_seconds
        << ",\"pairhmm_execution_space\":\"" << json_escape(result.pairhmm_execution_space) << "\""
        << ",\"pairhmm_prepare_seconds\":" << result.pairhmm_prepare_seconds
        << ",\"pairhmm_seconds\":" << result.pairhmm_seconds
        << ",\"sw_used\":" << (result.sw_used ? "true" : "false")
        << ",\"sw_pairs\":" << result.sw_pairs
        << ",\"sw_simd_width\":" << result.sw_simd_width
        << ",\"sw_simd_groups\":" << result.sw_simd_groups
        << ",\"sw_execution_space\":\"" << json_escape(result.sw_execution_space) << "\""
        << ",\"sw_prepare_seconds\":" << result.sw_prepare_seconds
        << ",\"sw_seconds\":" << result.sw_seconds
        << ",\"sw_traceback_pairs\":" << result.sw_traceback_pairs
        << ",\"sw_indel_alignments\":" << result.sw_indel_alignments
        << ",\"sw_cigar_signature\":" << result.sw_cigar_signature
        << ",\"sw_traceback_seconds\":" << result.sw_traceback_seconds
        << ",\"read_haplotype_cigar_pairs\":" << result.read_haplotype_cigar_pairs
        << ",\"read_haplotype_cigar_valid_pairs\":"
        << result.read_haplotype_cigar_valid_pairs
        << ",\"read_haplotype_cigar_indel_pairs\":"
        << result.read_haplotype_cigar_indel_pairs
        << ",\"read_haplotype_cigar_filtered_pairs\":"
        << result.read_haplotype_cigar_filtered_pairs
        << ",\"read_haplotype_softclip_filtered_pairs\":"
        << result.read_haplotype_softclip_filtered_pairs
        << ",\"read_haplotype_uncertain_reads\":"
        << result.read_haplotype_uncertain_reads
        << ",\"read_haplotype_informative_reads\":"
        << result.read_haplotype_informative_reads
        << ",\"read_haplotype_margin_sum\":"
        << result.read_haplotype_margin_sum
        << ",\"read_haplotype_uncertainty_signature\":"
        << result.read_haplotype_uncertainty_signature
        << ",\"graph_used\":" << (result.graph_used ? "true" : "false")
        << ",\"graph_nodes\":" << result.graph_nodes
        << ",\"graph_edges\":" << result.graph_edges
        << ",\"graph_branching_nodes\":" << result.graph_branching_nodes
        << ",\"graph_reference_nodes\":" << result.graph_reference_nodes
        << ",\"graph_reference_edges\":" << result.graph_reference_edges
        << ",\"graph_reference_connected_nodes\":" << result.graph_reference_connected_nodes
        << ",\"graph_dangling_nodes\":" << result.graph_dangling_nodes
        << ",\"graph_pruned_nodes\":" << result.graph_pruned_nodes
        << ",\"graph_dangling_branch_paths\":" << result.graph_dangling_branch_paths
        << ",\"graph_dangling_branch_bases\":" << result.graph_dangling_branch_bases
        << ",\"graph_dangling_recovered_paths\":" << result.graph_dangling_recovered_paths
        << ",\"graph_dangling_recovered_bases\":" << result.graph_dangling_recovered_bases
        << ",\"graph_artificial_haplotype_recovery_paths\":" << result.graph_artificial_haplotype_recovery_paths
        << ",\"graph_artificial_haplotype_recovery_bases\":" << result.graph_artificial_haplotype_recovery_bases
        << ",\"graph_seqgraph_nodes\":" << result.graph_seqgraph_nodes
        << ",\"graph_seqgraph_edges\":" << result.graph_seqgraph_edges
        << ",\"graph_seqgraph_linear_chain_merges\":" << result.graph_seqgraph_linear_chain_merges
        << ",\"graph_seqgraph_diamond_merges\":" << result.graph_seqgraph_diamond_merges
        << ",\"graph_seqgraph_tail_merges\":" << result.graph_seqgraph_tail_merges
        << ",\"graph_seqgraph_suffix_splits\":" << result.graph_seqgraph_suffix_splits
        << ",\"graph_seqgraph_suffix_merges\":" << result.graph_seqgraph_suffix_merges
        << ",\"graph_non_unique_kmers\":" << result.graph_non_unique_kmers
        << ",\"graph_reference_non_unique_kmers\":" << result.graph_reference_non_unique_kmers
        << ",\"graph_reference_kmer_rejected\":"
        << (result.graph_reference_kmer_rejected ? "true" : "false")
        << ",\"graph_reference_paths\":" << result.graph_reference_paths
        << ",\"graph_haplotype_paths\":" << result.graph_haplotype_paths
        << ",\"graph_haplotype_sequence_count\":"
        << retained_graph_paths.size()
        << ",\"graph_haplotype_sequence_max_length\":"
        << retained_graph_haplotype_max_length
        << ",\"graph_haplotype_path_details\":[";
    for (std::size_t detail = 0; detail < retained_graph_paths.size(); ++detail) {
        const auto path = retained_graph_paths[detail];
        if (detail != 0) out << ',';
        out << "{\"tid\":"
               << (path < result.graph_haplotype_path_tids.size()
                       ? result.graph_haplotype_path_tids[path] : -1)
               << ",\"start\":"
               << (path < result.graph_haplotype_path_starts.size()
                       ? result.graph_haplotype_path_starts[path] : -1)
               << ",\"end\":"
               << (path < result.graph_haplotype_path_ends.size()
                       ? result.graph_haplotype_path_ends[path] : -1)
               << ",\"support\":"
               << (path < result.graph_haplotype_path_support.size()
                       ? result.graph_haplotype_path_support[path] : 0)
               << ",\"non_reference_edge\":"
               << (path < result.graph_haplotype_path_has_non_reference_edge.size() &&
                   result.graph_haplotype_path_has_non_reference_edge[path] != 0 ? "true" : "false")
               << ",\"cigar\":\""
               << (path < result.graph_haplotype_path_cigars.size()
                       ? json_escape(result.graph_haplotype_path_cigars[path]) : std::string{})
               << "\",\"alignment_score\":"
               << (path < result.graph_haplotype_path_alignment_scores.size()
                       ? result.graph_haplotype_path_alignment_scores[path] : 0)
               << ",\"alignment_offset\":"
               << (path < result.graph_haplotype_path_alignment_offsets.size()
                       ? result.graph_haplotype_path_alignment_offsets[path] : 0)
               << ",\"sequence\":\"" << json_escape(result.graph_haplotype_sequences[path])
               << "\"}";
    }
    out << "]"
        << ",\"graph_haplotype_provenance_count\":"
        << retained_graph_paths.size()
        << ",\"graph_haplotype_event_map_count\":"
        << retained_graph_paths.size()
        << ",\"graph_haplotype_event_count\":"
        << result.graph_haplotype_event_count
        << ",\"graph_haplotype_cigar_signature\":"
        << result.graph_haplotype_cigar_signature
        << ",\"graph_haplotype_sw_used\":"
        << (result.graph_haplotype_sw_used ? "true" : "false")
        << ",\"graph_haplotype_sw_execution_space\":\""
        << json_escape(result.graph_haplotype_sw_execution_space) << "\""
        << ",\"graph_haplotype_sw_simd_width\":"
        << result.graph_haplotype_sw_simd_width
        << ",\"graph_haplotype_sw_simd_groups\":"
        << result.graph_haplotype_sw_simd_groups
        << ",\"graph_haplotype_sw_prepare_seconds\":"
        << result.graph_haplotype_sw_prepare_seconds
        << ",\"graph_haplotype_sw_seconds\":"
        << result.graph_haplotype_sw_seconds
        << ",\"graph_execution_space\":\"" << json_escape(result.graph_execution_space) << "\""
        << ",\"graph_prepare_seconds\":" << result.graph_prepare_seconds
        << ",\"graph_seconds\":" << result.graph_seconds
        << ",\"activity_used\":" << (result.activity_used ? "true" : "false")
        << ",\"activity_reference_aware\":"
        << (result.activity_reference_aware ? "true" : "false")
        << ",\"activity_loci\":" << result.activity_loci
        << ",\"active_loci\":" << result.active_loci
        << ",\"assembly_regions\":" << result.assembly_regions
        << ",\"activity_filter_size\":" << result.activity_filter_size
        << ",\"activity_effective_max_probability_propagation_distance\":"
        << result.activity_effective_max_probability_propagation_distance
        << ",\"active_probability_threshold\":"
        << options.calling.active_probability_threshold
        << ",\"force_active\":" << (options.calling.force_active ? "true" : "false")
        << ",\"assembly_region_padding\":"
        << options.calling.assembly_region_padding
        << ",\"min_assembly_region_size\":"
        << options.calling.min_assembly_region_size
        << ",\"max_assembly_region_size\":"
        << options.calling.max_assembly_region_size
        << ",\"max_probability_propagation_distance\":"
        << options.calling.max_probability_propagation_distance
        << ",\"assembly_unassigned_candidates\":" << result.assembly_unassigned_candidates
        << ",\"assembly_cross_region_candidates\":" << result.assembly_cross_region_candidates
        << ",\"assembly_region_union_count\":" << result.assembly_region_union_count
        << ",\"assembly_region_union_merges\":" << result.assembly_region_union_merges
        << ",\"assembly_cross_region_rescued\":" << result.assembly_cross_region_rescued
        << ",\"assembly_region_partitioned\":"
        << (result.assembly_region_partitioned ? "true" : "false")
        << ",\"activity_execution_space\":\"" << json_escape(result.activity_execution_space) << "\""
        << ",\"activity_prepare_seconds\":" << result.activity_prepare_seconds
        << ",\"activity_seconds\":" << result.activity_seconds
        << ",\"read_filter_used\":" << (result.read_filter_used ? "true" : "false")
        << ",\"read_filter_applied\":" << (result.read_filter_applied ? "true" : "false")
        << ",\"min_base_quality\":"
        << static_cast<unsigned>(options.calling.min_base_quality)
        << ",\"reference_model_deletion_quality\":"
        << static_cast<unsigned>(options.calling.reference_model_deletion_quality)
        << ",\"read_filter_min_mapping_quality\":"
        << static_cast<unsigned>(options.calling.minimum_mapping_quality)
        << ",\"read_filter_exclude_mapping_quality_unavailable\":"
        << (options.calling.exclude_mapping_quality_unavailable ? "true" : "false")
        << ",\"read_filter_exclude_mapping_quality_zero\":"
        << (options.calling.exclude_mapping_quality_zero ? "true" : "false")
        << ",\"read_filter_exclude_duplicates\":"
        << (options.calling.exclude_duplicates ? "true" : "false")
        << ",\"read_filter_exclude_unmapped\":"
        << (options.calling.exclude_unmapped ? "true" : "false")
        << ",\"read_filter_exclude_secondary\":"
        << (options.calling.exclude_secondary ? "true" : "false")
        << ",\"read_filter_exclude_supplementary\":"
        << (options.calling.exclude_supplementary ? "true" : "false")
        << ",\"read_filter_exclude_qcfail\":"
        << (options.calling.exclude_qcfail ? "true" : "false")
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
        << ",\"use_soft_clipped_bases\":"
        << (options.calling.use_soft_clipped_bases ? "true" : "false")
        << ",\"soft_clip_low_quality_ends\":"
        << (options.calling.soft_clip_low_quality_ends ? "true" : "false")
        << ",\"read_filter_max_reads_per_locus\":"
        << options.calling.max_reads_per_locus
        << ",\"read_filter_downsampling_seed\":"
        << options.calling.downsampling_seed
        << ",\"filtered_reads\":" << result.filtered_reads
        << ",\"downsampled_reads\":" << result.downsampled_reads
        << ",\"read_filter_execution_space\":\"" << json_escape(result.read_filter_execution_space) << "\""
        << ",\"read_filter_prepare_seconds\":" << result.read_filter_prepare_seconds
        << ",\"read_filter_seconds\":" << result.read_filter_seconds
        << ",\"overlapping_quality_correction_used\":"
        << (result.overlapping_quality_correction_used ? "true" : "false")
        << ",\"overlapping_quality_correction_metadata_available\":"
        << (result.overlapping_quality_correction_metadata_available ? "true" : "false")
        << ",\"overlapping_pairs\":" << result.overlapping_pairs
        << ",\"overlapping_bases\":" << result.overlapping_bases
        << ",\"overlapping_conflicting_bases\":"
        << result.overlapping_conflicting_bases
        << ",\"overlapping_quality_caps\":" << result.overlapping_quality_caps
        << ",\"do_not_correct_overlapping_base_qualities\":"
        << (options.calling.do_not_correct_overlapping_base_qualities ? "true" : "false")
        << ",\"error_correction_used\":"
        << (result.error_correction_used ? "true" : "false")
        << ",\"error_correction_solid_kmers\":" << result.error_correction_solid_kmers
        << ",\"error_correction_corrected_kmers\":" << result.error_correction_corrected_kmers
        << ",\"error_correction_uncorrectable_kmers\":" << result.error_correction_uncorrectable_kmers
        << ",\"error_correction_pileup_loci\":" << result.error_correction_pileup_loci
        << ",\"error_correction_pileup_skipped_indel_adjacent_bases\":"
        << result.error_correction_pileup_skipped_indel_adjacent_bases
        << ",\"error_correction_corrected_reads\":" << result.error_correction_corrected_reads
        << ",\"error_correction_corrected_bases\":" << result.error_correction_corrected_bases
        << ",\"error_correction_mode\":\""
        << json_escape(result.error_correction_mode) << "\""
        << ",\"error_correction_execution_space\":\""
        << json_escape(result.error_correction_execution_space) << "\""
        << ",\"error_correction_prepare_seconds\":" << result.error_correction_prepare_seconds
        << ",\"error_correction_seconds\":" << result.error_correction_seconds
        << ",\"reference_confidence_execution_space\":\""
        << json_escape(result.reference_confidence_execution_space) << "\""
        << ",\"reference_confidence_prepare_seconds\":"
        << result.reference_confidence_prepare_seconds
        << ",\"reference_confidence_seconds\":"
        << result.reference_confidence_seconds
        << ",\"reference_block_lookup_indexed\":"
        << (result.reference_block_lookup_indexed ? "true" : "false")
        << ",\"gvcf_reference_blocks\":"
        << (result.reference_blocks.empty() ? result.reference_block_count
                                             : result.reference_blocks.size())
        << ",\"gvcf_interval_reference_blocks\":"
        << (interval_subset ? "true" : "false")
        << ",\"interval_subset\":" << (interval_subset ? "true" : "false")
        << ",\"intervals\":" << options.calling.intervals.size()
        << ",\"interval_list_inputs\":" << reader.interval_file_inputs()
        << ",\"interval_list_records\":" << reader.interval_file_records()
        << ",\"gvcf_requested\":" << (options.gvcf ? "true" : "false")
        << ",\"gvcf_gq_bands\":[";
    for (std::size_t index = 0; index < options.calling.gvcf_gq_bands.size(); ++index) {
        if (index != 0) out << ',';
        out << options.calling.gvcf_gq_bands[index];
    }
    out << "]"
        << ",\"gvcf_semantics\":\""
        << (options.gvcf ? (options.calling.emit_reference_confidence_bp_resolution
            ? "base-pair-resolution+candidate-sites" : "reference-blocks+candidate-sites")
                         : "not-requested") << "\""
        << ",\"gvcf_standard_fields\":"
        << (options.gvcf && !options.calling.emit_reference_confidence_bp_resolution ? "true" : "false")
        << ",\"gvcf_candidate_standard_fields\":"
        << (options.gvcf ? "true" : "false")
        << ",\"gvcf_candidate_rebuild_after_genotyping\":"
        << (options.gvcf ? "true" : "false")
        << ",\"sample_ploidy\":" << options.sample_ploidy
        << ",\"genotype_priors\":" << (options.calling.use_genotype_priors ? "true" : "false")
        << ",\"genotype_assignment_method\":\""
        << (options.calling.use_posterior_genotype_assignment
            ? "USE_POSTERIOR_PROBABILITIES" : "USE_PLS_TO_ASSIGN") << "\""
        << ",\"joint_genotype_priors_used\":"
        << (result.genotype_priors_used ? "true" : "false")
        << ",\"genotype_prior_model\":\"GATK-assumingHW-v1\""
        << ",\"genotype_prior_snp_normalization\":\"log10(3)\""
        << ",\"genotype_prior_candidates\":[";
    for (std::size_t index = 0; index < result.candidates.size(); ++index) {
        if (index != 0) out << ',';
        out << "{\"ref\":\"" << json_escape(candidate_reference(result.candidates[index]))
            << "\",\"alt\":\"" << json_escape(candidate_alternate(result.candidates[index]))
            << "\",\"het\":"
            << std::setprecision(17)
            << (index < result.candidate_prior_het.size() ? result.candidate_prior_het[index] : 0.0)
            << ",\"hom_alt\":"
            << (index < result.candidate_prior_hom_alt.size() ? result.candidate_prior_hom_alt[index] : 0.0)
            << '}';
    }
    out << ']' << std::setprecision(6)
        << ",\"heterozygosity\":" << options.calling.heterozygosity
        << ",\"indel_heterozygosity\":" << options.calling.indel_heterozygosity
        << ",\"heterozygosity_stdev\":" << options.calling.heterozygosity_stdev
        << ",\"standard_confidence_for_calling\":"
        << options.calling.standard_confidence_for_calling
        << ",\"indel_size_to_eliminate_in_ref_model\":"
        << options.calling.indel_size_to_eliminate_in_ref_model
        << ",\"informative_read_overlap_margin\":"
        << options.calling.informative_read_overlap_margin
        << ",\"graph_kmer_size\":" << options.calling.graph_kmer_size
        << ",\"graph_kmer_sizes\":[";
    for (std::size_t index = 0; index < options.calling.graph_kmer_sizes.size(); ++index) {
        if (index != 0) out << ',';
        out << options.calling.graph_kmer_sizes[index];
    }
    out << ']'
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
        << ",\"graph_max_unpruned_variants\":" << options.calling.graph_max_unpruned_variants
        << ",\"graph_adaptive_pruned_nodes\":" << result.graph_adaptive_pruned_nodes
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
        << ",\"graph_kmer_size_selected\":" << result.graph_kmer_size_selected
        << ",\"graph_kmer_iterations\":" << result.graph_kmer_iterations
        << ",\"graph_has_non_reference_cycles\":"
        << (result.graph_has_non_reference_cycles ? "true" : "false")
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
        << ",\"graph_max_paths\":" << options.calling.graph_max_paths
        << ",\"graph_max_depth\":" << options.calling.graph_max_depth
        << ",\"max_haplotype_combination_alleles\":"
        << options.calling.max_haplotype_combination_alleles
        << ",\"flow_assembly_collapse_hmer_size\":"
        << options.calling.flow_assembly_collapse_hmer_size
        << ",\"flow_assembly_collapse_partial_mode\":"
        << (options.calling.flow_assembly_collapse_partial_mode ? "true" : "false")
        << ",\"max_mnp_distance\":" << options.calling.max_mnp_distance
        << ",\"max_alternate_alleles\":" << options.max_alternate_alleles
        << ",\"max_genotype_count\":" << options.max_genotype_count
        << ",\"haplotype_pruning_log10\":" << options.calling.haplotype_pruning_log10
        << ",\"input\":\"" << json_escape(options.input) << "\""
        << ",\"reads\":" << result.reads
        << ",\"observations\":" << result.observations
        << ",\"projected_observations\":" << result.projected_observations
        << ",\"loci\":" << result.loci
        << ",\"candidate_sites\":" << result.candidate_sites
        << ",\"candidate_softclip_suppressed\":" << result.candidate_softclip_suppressed
        << ",\"candidate_fragment_suppressed\":" << result.candidate_fragment_suppressed
        << ",\"candidate_low_support_suppressed\":" << result.candidate_low_support_suppressed
        << ",\"rcm_haplotype_realignment_used\":"
        << (result.rcm_haplotype_realignment_used ? "true" : "false")
        << ",\"rcm_realigned_observations\":" << result.rcm_realigned_observations
        << ",\"rcm_realignment_fallback_observations\":"
        << result.rcm_realignment_fallback_observations
        << ",\"multiallelic_candidates_pruned\":" << result.multiallelic_candidates_pruned
        << ",\"softclip_candidates\":" << result.softclip_candidates
        << ",\"compound_cigar_candidates\":" << result.compound_cigar_candidates
        << ",\"graph_variant_candidates\":" << result.graph_variant_candidates
        << ",\"graph_snp_candidates\":" << result.graph_snp_candidates
        << ",\"graph_mnp_candidates\":" << result.graph_mnp_candidates
        << ",\"graph_deletion_candidates\":" << result.graph_deletion_candidates
        << ",\"variant_calls\":" << result.variant_calls
        << ",\"signature\":" << result.signature
        << ",\"assembly_prepare_seconds\":" << std::setprecision(12) << result.assembly_prepare_seconds
        << ",\"assembly_seconds\":" << result.assembly_seconds
        << ",\"likelihood_prepare_seconds\":" << result.likelihood_prepare_seconds
        << ",\"likelihood_seconds\":" << result.likelihood_seconds
        << ",\"genotyping_seconds\":" << result.genotyping_seconds
        << ",\"genotype_prior_kernel_calls\":" << result.genotype_prior_kernel_calls
        << ",\"genotype_prior_prepare_seconds\":" << result.genotype_prior_prepare_seconds
        << ",\"genotype_prior_seconds\":" << result.genotype_prior_seconds
        << ",\"genotype_prior_execution_space\":\""
        << json_escape(result.genotype_prior_execution_space) << "\""
        << ",\"call_confidence_execution_space\":\""
        << json_escape(result.call_confidence_execution_space) << "\""
        << ",\"call_confidence_prepare_seconds\":"
        << result.call_confidence_prepare_seconds
        << ",\"call_confidence_seconds\":" << result.call_confidence_seconds
        << ",\"calls\":[";
    for (std::size_t i = 0; i < result.calls.size(); ++i) {
        const auto& call = result.calls[i];
        if (i != 0) out << ',';
        const char* genotype = call.genotype == 1 ? "0/1" : "1/1";
        const auto ref = candidate_reference(call.candidate);
        const auto alt = candidate_alternate(call.candidate);
        out << "{\"tid\":" << call.candidate.tid
            << ",\"position\":" << call.candidate.position
            << ",\"ref\":\"" << json_escape(ref)
            << "\",\"alt\":\"" << json_escape(alt)
            << "\",\"genotype\":\"" << genotype << "\",\"gq\":" << static_cast<int>(call.gq)
            << ",\"depth\":" << call.candidate.depth
            << ",\"variant_depth\":" << call.candidate.variant_depth
            << ",\"ref_count\":" << call.candidate.reference_count
            << ",\"alt_count\":" << call.candidate.alternate_count
            << ",\"qual\":" << call.qual
            << ",\"mq\":" << (std::isfinite(call.annotations.mq) ? call.annotations.mq : -1.0)
            << ",\"qd\":" << (std::isfinite(call.annotations.qd) ? call.annotations.qd : -1.0)
            << ",\"fs\":" << (std::isfinite(call.annotations.fs) ? call.annotations.fs : -1.0)
            << ",\"sor\":" << (std::isfinite(call.annotations.sor) ? call.annotations.sor : -1.0)
            << ",\"mq_rank_sum\":"
            << (std::isfinite(call.annotations.mq_rank_sum) ? call.annotations.mq_rank_sum : -1.0)
            << ",\"read_pos_rank_sum\":"
            << (std::isfinite(call.annotations.read_pos_rank_sum) ? call.annotations.read_pos_rank_sum : -1.0)
            << ",\"base_q_rank_sum\":"
            << (std::isfinite(call.annotations.base_q_rank_sum) ? call.annotations.base_q_rank_sum : -1.0)
            << ",\"log10_likelihoods\":[" << call.likelihoods.hom_ref << ','
            << call.likelihoods.het << ',' << call.likelihoods.hom_alt << "]}";
    }
    out << "]}\n";
    return out.str();
}

// Materialize the GATK AlleleFrequencyCalculator input and delegate the
// bounded EM reduction for INFO MLEAC/MLEAF to the shared Kokkos cohort
// primitive.  PLs remain prior-free; the AF calculation applies the
// assuming-HW Dirichlet pseudocounts only while estimating allele counts.
// Host still owns the small VCF/allele-classification boundary, but the
// numerically important EM/posterior loop is now the same API used by
// GenotypeGVCFs and the standalone cohort kernels on every backend.
std::vector<int> estimate_mle_allele_counts(
    const fastgatk::calling::Result& result,
    const std::vector<const fastgatk::calling::AssemblyCandidate*>& candidates,
    const std::vector<std::vector<int>>& genotype_vectors,
    const std::vector<int>& pl,
    int sample_ploidy,
    bool include_non_ref = false) {
    if (candidates.empty() || genotype_vectors.empty() || pl.size() != genotype_vectors.size() ||
        sample_ploidy <= 0 || !std::isfinite(result.genotype_heterozygosity_stdev) ||
        result.genotype_heterozygosity_stdev <= 0.0)
        return std::vector<int>(candidates.size(), 0);
    const auto allele_count = candidates.size() + 1 + (include_non_ref ? 1U : 0U);
    const auto ref_pseudocount = result.genotype_snp_heterozygosity /
        (result.genotype_heterozygosity_stdev * result.genotype_heterozygosity_stdev);
    if (!(ref_pseudocount > 0.0) || !std::isfinite(ref_pseudocount))
        return std::vector<int>(candidates.size(), 0);
    std::vector<double> prior(allele_count, result.genotype_snp_heterozygosity * ref_pseudocount);
    prior[0] = ref_pseudocount;
    for (std::size_t alt = 0; alt < candidates.size(); ++alt) {
        const auto& candidate = *candidates[alt];
        const bool indel = candidate_reference(candidate).size() !=
                           candidate_alternate(candidate).size();
        prior[alt + 1] = (indel ? result.genotype_indel_heterozygosity
                                : result.genotype_snp_heterozygosity) * ref_pseudocount;
    }
    if (include_non_ref)
        prior.back() = std::max(result.genotype_snp_heterozygosity,
                                result.genotype_indel_heterozygosity) * ref_pseudocount;
    // The public kernel takes a sample-major matrix.  The HC writer always
    // computes one sample at a time, so the existing Number=G vector is
    // already exactly the required one-row layout.
    const auto expected_width = genotype_width(static_cast<int>(allele_count), sample_ploidy);
    if (expected_width == 0 || genotype_vectors.size() != expected_width ||
        pl.size() != expected_width)
        return std::vector<int>(candidates.size(), 0);
    fastgatk::kernels::AlleleFrequencyResult af;
    try {
        af = fastgatk::kernels::calculate_allele_frequency_kokkos(
            pl, 1, static_cast<int>(allele_count), sample_ploidy, prior);
    } catch (const std::exception&) {
        // Preserve the writer's historical fail-closed behavior for malformed
        // or numerically unusable Number=G rows; a VCF record is still
        // structurally emitted, while the unsupported AF path is not allowed
        // to fabricate allele counts.
        return std::vector<int>(candidates.size(), 0);
    }
    std::vector<int> result_counts(candidates.size(), 0);
    for (std::size_t alt = 0; alt < candidates.size(); ++alt)
        if (alt + 1 < af.integer_allele_counts.size())
            result_counts[alt] = std::clamp(
                static_cast<int>(af.integer_allele_counts[alt + 1]), 0, sample_ploidy);
    return result_counts;
}

std::string format_allele_frequency(const double frequency) {
    std::ostringstream value;
    if (std::abs(frequency - 1.0) < 1.0e-12)
        value << std::fixed << std::setprecision(2) << frequency;
    else
        value << std::fixed << std::setprecision(3) << frequency;
    return value.str();
}

// QualByDepth intentionally de-jitters unusually high values so they do not
// form an artificial VQSR tail.  The implementation uses Utils' process-wide
// java.util.Random, seeded to 47382911, rather than a hash of the variant.
// Keep a small Host-local implementation here: annotations are serialized in
// coordinate order and the numerical PairHMM/Kokkos paths remain untouched.
class GatkJavaRandom {
public:
    explicit GatkJavaRandom(const std::uint64_t seed = 47382911ULL)
        : seed_((seed ^ kMultiplier) & kMask) {}

    double next_gaussian() {
        if (has_cached_gaussian_) {
            has_cached_gaussian_ = false;
            return cached_gaussian_;
        }
        double first = 0.0;
        double second = 0.0;
        double squared = 0.0;
        do {
            first = 2.0 * next_double() - 1.0;
            second = 2.0 * next_double() - 1.0;
            squared = first * first + second * second;
        } while (squared >= 1.0 || squared == 0.0);
        const auto scale = std::sqrt(-2.0 * std::log(squared) / squared);
        cached_gaussian_ = second * scale;
        has_cached_gaussian_ = true;
        return first * scale;
    }

private:
    static constexpr std::uint64_t kMultiplier = 0x5DEECE66DULL;
    static constexpr std::uint64_t kAddend = 0xBULL;
    static constexpr std::uint64_t kMask = (1ULL << 48U) - 1ULL;

    std::uint64_t next_bits(const unsigned int bits) {
        seed_ = (seed_ * kMultiplier + kAddend) & kMask;
        return seed_ >> (48U - bits);
    }

    double next_double() {
        const auto value = (next_bits(26U) << 27U) + next_bits(27U);
        return static_cast<double>(value) / static_cast<double>(1ULL << 53U);
    }

    std::uint64_t seed_;
    bool has_cached_gaussian_ = false;
    double cached_gaussian_ = 0.0;
};

double gatk_qual_by_depth(const double qd, GatkJavaRandom& random) {
    constexpr double kMaxQdBeforeFixing = 35.0;
    constexpr double kIdealHighQd = 30.0;
    constexpr double kJitterSigma = 3.0;
    return std::isfinite(qd) && qd >= kMaxQdBeforeFixing
        ? kIdealHighQd + random.next_gaussian() * kJitterSigma
        : qd;
}

std::string vcf_text(const fastgatk::io::HtsReader& reader,
                fastgatk::calling::Result& result,
                int sample_ploidy,
                bool add_output_vcf_command_line = true,
                bool sites_only_vcf_output = false,
                const fastgatk::io::ReadBatch* annotation_reads = nullptr,
                std::uint32_t informative_read_overlap_margin = 2,
                std::size_t max_alternate_alleles = 6) {
    std::ostringstream out;
    // GATK registers XC only when Flow HMER collapsing produced a collapsed
    // EventMap event. Ordinary HaplotypeCaller output must not gain this
    // flow-specific schema line.
    const bool has_flow_hmer_collapsed = std::any_of(
        result.calls.begin(), result.calls.end(), [](const auto& call) {
            return call.candidate.flow_hmer_collapsed;
        });
    out << "##fileformat=VCFv4.2\n"
        << "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
        << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n"
        << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n"
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
        << "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n"
        << "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">\n"
        << "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
        << "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">\n"
        << "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n"
        << "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
        << "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">\n"
        << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n"
        << "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">\n";
    if (has_flow_hmer_collapsed)
        out << "##INFO=<ID=XC,Number=1,Type=Integer,Description=\"Indicates longer hmer collapsing took place (this is a flow-based specific tag)\">\n";
    const auto& header = reader.header();
    for (std::size_t i = 0; i < header.contigs.size(); ++i)
        out << "##contig=<ID=" << header.contigs[i]
            << ",length=" << header.contig_lengths[i]
            << (i < header.contig_assemblies.size() && !header.contig_assemblies[i].empty()
                    ? ",assembly=" + header.contig_assemblies[i] : std::string{})
            << ">\n";
    // GATK adds this source/provenance pair together when
    // --add-output-vcf-command-line is enabled.  Native implementation
    // details remain in the manifest; this stable source line is GATK's own
    // ordinary-HC header contract.
    if (add_output_vcf_command_line)
        out << "##source=HaplotypeCaller\n";
    if (add_output_vcf_command_line)
        out << "##GATKCommandLine=<ID=HaplotypeCaller,Version=fastgatk-native,"
               "CommandLine=\"fastgatk HaplotypeCaller\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if (!sites_only_vcf_output)
        out << "\tFORMAT\t" << output_sample_name(header);
    out << "\n";
    using GroupKey = std::tuple<std::int32_t, std::int32_t, std::string>;
    std::map<GroupKey, std::vector<const fastgatk::calling::GenotypeCall*>> groups;
    for (const auto& call : result.calls)
        groups[GroupKey{call.candidate.tid, call.candidate.position,
                        candidate_reference(call.candidate)}].push_back(&call);
    const auto genotype_index = [](std::size_t left, std::size_t right) {
        if (left > right) std::swap(left, right);
        return right * (right + 1) / 2 + left;
    };
    const auto qual_text = [](double qual) {
        std::ostringstream value;
        value << std::fixed << std::setprecision(2)
              << std::max(0.0, qual);
        return value.str();
    };
    const auto annotation_text = [](double value, int precision) {
        if (!std::isfinite(value)) return std::string{};
        // HTSJDK serializes a rounded signed zero as 0.000.  Canonicalize at
        // the VCF boundary rather than changing the rank-sum calculation,
        // whose last bit is allowed to be negative for a tied distribution.
        if (std::abs(value) < 0.5 * std::pow(10.0, -precision)) value = 0.0;
        std::ostringstream output;
        // Rank-sum annotations are signed; only the statistical annotations
        // that are intrinsically non-negative are clamped at calculation
        // time.  Do not turn a valid negative ALT-vs-REF z score into zero at
        // VCF serialization.
        output << std::fixed << std::setprecision(precision) << value;
        return output.str();
    };
    GatkJavaRandom qd_random;
    for (const auto& [key, grouped_calls] : groups) {
        // GATK's first allele-subsetting pass may replace the group's
        // concrete ALT list before any AF posterior is calculated.  Keep a
        // local Host-owned vector so the immutable event map remains usable
        // for other records/streaming tiles.
        std::vector<const fastgatk::calling::GenotypeCall*> calls = grouped_calls;
        if (calls.empty()) continue;
        const auto tid_value = std::get<0>(key);
        const auto tid = static_cast<std::size_t>(std::max<std::int32_t>(0, tid_value));
        const std::string chrom = tid < header.contigs.size() ? header.contigs[tid]
                                                               : std::to_string(tid_value);
        const auto position = std::get<1>(key);
        const auto& ref = std::get<2>(key);
        // A full contig result keeps each AssemblyRegion's PairHMM matrix in
        // its owning Result to avoid retaining a dense global copy.  Prefer
        // the direct matrix for a non-streamed result, then resolve the
        // candidate tuple against the owning region exactly as the joint-PL
        // path below does.
        const auto likelihood_depths_from_owner =
            [&](const std::vector<const fastgatk::calling::GenotypeCall*>& selected_calls)
                -> std::optional<std::pair<std::vector<int>, int>> {
                if (const auto direct = derive_multiallelic_depths(result, selected_calls);
                    direct.has_value())
                    return direct;
                for (const auto& owner : result.assembly_region_likelihood_results) {
                    if (owner == nullptr) continue;
                    if (const auto owned = derive_multiallelic_depths(*owner, selected_calls);
                        owned.has_value())
                        return owned;
                }
                return std::nullopt;
            };
        std::vector<int> alt_counts;
        alt_counts.reserve(calls.size());
        std::size_t best = 0;
        for (std::size_t i = 0; i < calls.size(); ++i) {
            alt_counts.push_back(static_cast<int>(calls[i]->candidate.alternate_count));
            if (calls[i]->gq > calls[best]->gq ||
                (calls[i]->gq == calls[best]->gq &&
                 calls[i]->candidate.alternate_count > calls[best]->candidate.alternate_count))
                best = i;
        }
        std::size_t allele_count = calls.size() + 1;
        const auto genotype_vectors = enumerate_genotypes(static_cast<int>(allele_count), sample_ploidy);
        if (genotype_vectors.empty())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: genotype PL width exceeds bounded sample-ploidy limit");
        std::vector<int> pl(genotype_vectors.size(), 999);
        std::vector<int> ad;
        ad.push_back(static_cast<int>(calls[best]->candidate.reference_count));
        int depth = static_cast<int>(calls[best]->candidate.depth);
        int variant_depth = static_cast<int>(calls[best]->candidate.variant_depth);
        if (variant_depth == 0 && depth != 0) variant_depth = depth;
        bool multiallelic_ad = false;
        for (const auto& row : result.multiallelic_depths) {
            if (row.tid != std::get<0>(key) || row.position != position ||
                row.reference != ref || row.depths.size() != calls.size() + 1)
                continue;
            std::vector<int> projected(calls.size() + 1, 0);
            if (row.depths[0] > static_cast<std::uint32_t>(std::numeric_limits<int>::max()))
                continue;
            projected[0] = static_cast<int>(row.depths[0]);
            bool matched = true;
            for (std::size_t call_index = 0; call_index < calls.size(); ++call_index) {
                const auto alt = candidate_alternate(calls[call_index]->candidate);
                const auto alt_iter = std::find(row.alternates.begin(), row.alternates.end(), alt);
                if (alt_iter == row.alternates.end()) {
                    matched = false;
                    break;
                }
                const auto row_index = static_cast<std::size_t>(alt_iter - row.alternates.begin()) + 1;
                if (row_index >= row.depths.size() ||
                    row.depths[row_index] > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
                    matched = false;
                    break;
                }
                projected[call_index + 1] = static_cast<int>(row.depths[row_index]);
            }
            if (matched) {
                ad = std::move(projected);
                depth = static_cast<int>(row.depth);
                multiallelic_ad = true;
            }
            break;
        }
        if (!multiallelic_ad) {
            if (const auto likelihood_depths = likelihood_depths_from_owner(calls);
                likelihood_depths.has_value()) {
                ad = likelihood_depths->first;
                depth = likelihood_depths->second;
                multiallelic_ad = true;
            }
        }
        const int an = sample_ploidy;
        double site_qual = 0.0;
        for (const auto* call : calls)
            if (call != nullptr) site_qual = std::max(site_qual, call->qual);
        auto annotations = calls[best] != nullptr
            ? calls[best]->annotations : fastgatk::calling::GenotypeCall::Annotations{};
        if (!multiallelic_ad)
            for (const auto* call : calls)
                ad.push_back(static_cast<int>(call->candidate.alternate_count));
        // The partition reducer deliberately flattens calls, while a
        // reference-backed PairHMM matrix remains owned by the AssemblyRegion
        // that produced them.  Probe that owner before falling back: this is
        // a Host ownership lookup, and keeps normal VCF PLs on the Kokkos
        // joint-genotyping result just as the gVCF <NON_REF> writer does.
        const auto joint_from_owner = [&](std::vector<int>& destination,
                                          bool include_spanning_deletion = false) {
            if (sample_ploidy != 2) {
                std::vector<const fastgatk::calling::AssemblyCandidate*> candidates;
                candidates.reserve(calls.size());
                for (const auto* call : calls) {
                    if (call == nullptr) return false;
                    candidates.push_back(&call->candidate);
                }
                if (joint_candidate_pl(result, candidates, sample_ploidy, destination, 0,
                                       include_spanning_deletion)) return true;
                for (const auto& owner : result.assembly_region_likelihood_results) {
                    if (owner != nullptr &&
                        joint_candidate_pl(*owner, candidates, sample_ploidy, destination, 0,
                                           include_spanning_deletion))
                        return true;
                }
                return false;
            }
            // The established diploid writer has its own spanning-deletion
            // pre-subsetting path below.  Keep this ordinary-VCF ownership
            // probe concrete-only for the ploidy-two contract.
            if (include_spanning_deletion) return false;
            if (joint_concrete_pl(result, calls, sample_ploidy, destination)) return true;
            for (const auto& owner : result.assembly_region_likelihood_results) {
                if (owner != nullptr &&
                    joint_concrete_pl(*owner, calls, sample_ploidy, destination))
                    return true;
            }
            return false;
        };
        std::vector<int> spanning_pl;
        const bool has_spanning_deletion = sample_ploidy != 2 &&
            joint_from_owner(spanning_pl, true);
        bool joint = false;
        if (has_spanning_deletion) {
            // GenotypingEngine calculates GL/AF on REF/concrete ALT(s)/*,
            // then drops a spurious `*` from the ordinary VCF output.  Its
            // allele-subsetting utility selects the concrete Number=G rows
            // and renormalizes them; genotype enumeration guarantees those
            // rows precede every genotype containing the final `*` allele.
            const auto output_width = genotype_width(static_cast<int>(allele_count), sample_ploidy);
            if (output_width == 0 || spanning_pl.size() < output_width)
                throw std::runtime_error(
                    "NUMERICAL_CONTRACT_FAILURE: spanning-deletion PL prefix is incomplete");
            pl.assign(spanning_pl.begin(), spanning_pl.begin() + output_width);
            const auto minimum = *std::min_element(pl.begin(), pl.end());
            for (auto& value : pl)
                value = std::clamp(value - minimum, 0,
                                   std::numeric_limits<int>::max());
            joint = true;
        } else {
            joint = joint_from_owner(pl);
        }
        if (!joint && sample_ploidy != 2)
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: --sample-ploidy requires reference-backed PairHMM likelihoods");
        if (!joint && sample_ploidy == 2) {
            // Inputs without a reference-backed PairHMM matrix retain the
            // deterministic biallelic fallback.  Reference-backed calls use
            // the joint path above, so cross-ALT PLs are not fabricated from
            // independently normalized penalties.
            std::vector<std::array<int, 3>> biallelic_pl;
            biallelic_pl.reserve(calls.size());
            for (const auto* call : calls) biallelic_pl.push_back(phred_likelihoods(call->likelihoods));
            pl[genotype_index(0, 0)] = 999;
            for (std::size_t i = 0; i < calls.size(); ++i) {
                pl[genotype_index(0, 0)] = std::min(pl[genotype_index(0, 0)], biallelic_pl[i][0]);
                pl[genotype_index(0, i + 1)] = biallelic_pl[i][1];
                pl[genotype_index(i + 1, i + 1)] = biallelic_pl[i][2];
            }
            for (std::size_t left = 0; left < calls.size(); ++left)
                for (std::size_t right = left + 1; right < calls.size(); ++right)
                    pl[genotype_index(left + 1, right + 1)] =
                        std::max(biallelic_pl[left][2], biallelic_pl[right][2]);
            const auto minimum_pl = *std::min_element(pl.begin(), pl.end());
            for (auto& value : pl) value = std::clamp(value - minimum_pl, 0,
                std::numeric_limits<int>::max());
        }
        // GenotypingEngine reduces a high-ALT VariantContext *before* it
        // invokes AFCalculator.  Its selection uses the complete GL matrix,
        // then AlleleSubsettingUtils recodes Number=G and Number=R fields.
        // Do that structural work on the Host and leave the subsequent AF
        // posterior to calculate_allele_frequency_kokkos.
        bool max_allele_subset_changed = false;
        if (calls.size() > max_alternate_alleles) {
            const auto retained_alt_indices = gatk_most_likely_alt_indices(
                pl, static_cast<int>(allele_count), sample_ploidy, calls.size(),
                max_alternate_alleles);
            std::vector<std::size_t> new_to_old{0U};
            new_to_old.reserve(retained_alt_indices.size() + 1U);
            for (const auto alt : retained_alt_indices) new_to_old.push_back(alt + 1U);
            pl = remap_pl_for_allele_subset(
                pl, static_cast<int>(allele_count), sample_ploidy, new_to_old);
            if (has_spanning_deletion) {
                auto spanning_new_to_old = new_to_old;
                // In ordinary VCF mode `*` is held only for the AF
                // calculation.  Preserve it while remapping that internal
                // likelihood matrix, exactly as GATK does before suppressing
                // it from the public ALT list.
                spanning_new_to_old.push_back(allele_count);
                spanning_pl = remap_pl_for_allele_subset(
                    spanning_pl, static_cast<int>(allele_count + 1U), sample_ploidy,
                    spanning_new_to_old);
            }
            if (ad.size() == calls.size() + 1U)
                ad = remap_r_length_values(ad, new_to_old);
            std::vector<const fastgatk::calling::GenotypeCall*> retained_calls;
            retained_calls.reserve(retained_alt_indices.size());
            for (const auto alt : retained_alt_indices) retained_calls.push_back(calls[alt]);
            calls = std::move(retained_calls);
            allele_count = calls.size() + 1U;
            max_allele_subset_changed = true;
        }
        const auto frequency_allele_count = allele_count + (has_spanning_deletion ? 1U : 0U);
        const auto& frequency_pl = has_spanning_deletion ? spanning_pl : pl;
        std::optional<fastgatk::kernels::AlleleFrequencyResult> joint_frequency;
        std::optional<double> joint_qd;
        if ((calls.size() > 1 || sample_ploidy != 2) && !frequency_pl.empty()) {
            // GATK's site QUAL for a multi-ALT record, and for every
            // arbitrary-ploidy record, is the AFCalculator posterior that
            // the sample is not hom-ref rather than the maximum of
            // independently normalized biallelic ALT qualities.  The former
            // branch condition accidentally skipped a single-ALT triploid
            // call, leaving its QUAL/MLEAC/QD at a diploid provisional value.
            // Reuse the shared Kokkos EM/posterior kernel for both shapes.
            // Reuse the shared Kokkos EM/posterior kernel so the result is
            // backend-independent and follows the same 0.1 convergence rule
            // as GenotypeGVCFs.
            std::vector<double> prior_pseudocounts(frequency_allele_count, 0.0);
            const auto ref_pseudocount = result.genotype_snp_heterozygosity /
                (result.genotype_heterozygosity_stdev * result.genotype_heterozygosity_stdev);
            prior_pseudocounts[0] = ref_pseudocount;
            for (std::size_t alt = 0; alt < calls.size(); ++alt) {
                const auto& candidate = calls[alt]->candidate;
                const auto indel = candidate_reference(candidate).size() !=
                                   candidate_alternate(candidate).size();
                const auto heterozygosity = indel
                    ? result.genotype_indel_heterozygosity
                    : result.genotype_snp_heterozygosity;
                prior_pseudocounts[alt + 1] = heterozygosity * ref_pseudocount;
            }
            if (has_spanning_deletion)
                prior_pseudocounts.back() = result.genotype_indel_heterozygosity *
                    ref_pseudocount;
            joint_frequency = fastgatk::kernels::calculate_allele_frequency_kokkos(
                frequency_pl, 1, static_cast<int>(frequency_allele_count), sample_ploidy,
                prior_pseudocounts, {}, {}, has_spanning_deletion
                    ? static_cast<int>(calls.size() + 1U) : -1);
            if (joint_frequency->samples_with_likelihoods != 0 &&
                std::isfinite(joint_frequency->qual))
                site_qual = joint_frequency->qual;
            // GenotypingEngine retains a given ALT even when no ALT passes
            // the AF plausibility threshold. In that monomorphic forced-site
            // case its log10PError is P(variant present), the complement of
            // the joint Kokkos AF no-variant posterior; ordinary and
            // plausibly variant multi-ALT sites keep the standard QUAL.
            const bool has_forced_feature_allele = std::any_of(
                calls.begin(), calls.end(), [](const auto* call) {
                    return call != nullptr && call->candidate.forced_by_alleles_feature;
                });
            if (has_forced_feature_allele &&
                joint_frequency->log10_p_allele_absent.size() == frequency_allele_count) {
                constexpr double kGatkAfThresholdEpsilon = 1.0e-10;
                const auto threshold = -0.1 * result.genotype_standard_confidence_for_calling;
                bool any_plausible = false;
                for (std::size_t alt = 0; alt < calls.size(); ++alt) {
                    const auto absent = joint_frequency->log10_p_allele_absent[alt + 1U];
                    any_plausible = any_plausible || (std::isfinite(absent) &&
                        absent + kGatkAfThresholdEpsilon < threshold);
                }
                const auto log10_no_variant = joint_frequency->log10_p_no_variant;
                if (!any_plausible && std::isfinite(log10_no_variant) &&
                    log10_no_variant < 0.0) {
                    const auto p_variant_present = 1.0 - std::pow(10.0, log10_no_variant);
                    if (p_variant_present > 0.0)
                        site_qual = -10.0 * std::log10(p_variant_present);
                }
            }
            int variant_sample_depth = 0;
            for (const auto value : ad)
                if (value > 0 && variant_sample_depth <=
                    std::numeric_limits<int>::max() - value)
                    variant_sample_depth += value;
            if (variant_sample_depth > 0) {
                annotations.qd = site_qual / static_cast<double>(variant_sample_depth);
                joint_qd = annotations.qd;
            }
        }
        // GenotypingEngine's multi-ALT boundary is not "emit every
        // independently biallelic call".  It asks AFCalculationResult if
        // each concrete ALT is plausible in the *joint* genotype matrix and
        // subsets the VC before it materializes GT/AD/PL/AC fields.  In
        // particular, a weak A sibling beside a confident C/C call must be
        // removed even if A passed an earlier biallelic candidate gate.  The
        // posterior and EM work stay in calculate_allele_frequency_kokkos;
        // this Host code only performs the structural VCF remap.
        std::vector<const fastgatk::calling::GenotypeCall*> output_calls = calls;
        std::vector<std::size_t> output_alt_indices(calls.size(), 0);
        std::iota(output_alt_indices.begin(), output_alt_indices.end(), 0U);
        bool joint_af_allele_subset_changed = false;
        if (joint_frequency.has_value() &&
            joint_frequency->log10_p_allele_absent.size() == frequency_allele_count) {
            constexpr double kGatkAfThresholdEpsilon = 1.0e-10;
            const auto threshold = -0.1 * result.genotype_standard_confidence_for_calling;
            std::vector<const fastgatk::calling::GenotypeCall*> retained_calls;
            std::vector<std::size_t> retained_alt_indices;
            retained_calls.reserve(calls.size());
            retained_alt_indices.reserve(calls.size());
            for (std::size_t alt = 0; alt < calls.size(); ++alt) {
                const auto absent = joint_frequency->log10_p_allele_absent[alt + 1];
                const bool plausible = std::isfinite(absent) &&
                    absent + kGatkAfThresholdEpsilon < threshold;
                const bool forced_feature = calls[alt] != nullptr &&
                    calls[alt]->candidate.forced_by_alleles_feature;
                if (plausible || forced_feature) {
                    retained_calls.push_back(calls[alt]);
                    retained_alt_indices.push_back(alt);
                }
            }
            // HaplotypeCaller in EMIT_VARIANTS_ONLY mode does not serialize a
            // now-monomorphic context.
            if (retained_calls.empty()) continue;
            if (retained_calls.size() != calls.size()) {
                const auto remapped_genotypes = enumerate_genotypes(
                    static_cast<int>(retained_calls.size() + 1), sample_ploidy);
                std::vector<int> remapped_pl(remapped_genotypes.size(), 999);
                for (std::size_t genotype = 0; genotype < remapped_genotypes.size(); ++genotype) {
                    const auto old_left = remapped_genotypes[genotype][0] == 0 ? 0U :
                        retained_alt_indices[static_cast<std::size_t>(
                            remapped_genotypes[genotype][0] - 1)] + 1U;
                    const auto old_right = remapped_genotypes[genotype][1] == 0 ? 0U :
                        retained_alt_indices[static_cast<std::size_t>(
                            remapped_genotypes[genotype][1] - 1)] + 1U;
                    const auto old_index = genotype_index(old_left, old_right);
                    if (old_index < pl.size()) remapped_pl[genotype] = pl[old_index];
                }
                pl = std::move(remapped_pl);
                if (ad.size() == calls.size() + 1) {
                    std::vector<int> remapped_ad;
                    remapped_ad.reserve(retained_calls.size() + 1);
                    remapped_ad.push_back(ad.front());
                    for (const auto alt : retained_alt_indices)
                        remapped_ad.push_back(ad[alt + 1]);
                    ad = std::move(remapped_ad);
                }
                output_calls = std::move(retained_calls);
                output_alt_indices = std::move(retained_alt_indices);
                joint_af_allele_subset_changed = true;
            }
        }
        // GATK annotates the final VariantContext's allele set, which can be
        // narrower than both the candidate list and `calls`: a weak sibling
        // can fail its earlier biallelic gate before joint AF subsetting.  In
        // either case, re-run BestAllele over exactly the published alleles
        // rather than retaining candidate-local AD columns.  A read that
        // favored removed A in a T,A,C locus can still be informative for
        // the emitted T/C pair and must be counted as C-supporting.
        if (const auto output_likelihood_depths =
                likelihood_depths_from_owner(output_calls);
            output_likelihood_depths.has_value()) {
            ad = output_likelihood_depths->first;
            depth = output_likelihood_depths->second;
        }
        const auto output_allele_count = output_calls.size() + 1;
        const auto output_genotype_vectors = enumerate_genotypes(
            static_cast<int>(output_allele_count), sample_ploidy);
        if (output_genotype_vectors.empty())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: output genotype PL width exceeds bounded sample-ploidy limit");
        // The original winner may be an ALT that joint-AF subsetting removed.
        // Re-select the annotation owner from published alleles.
        best = 0;
        for (std::size_t i = 1; i < output_calls.size(); ++i)
            if (output_calls[i]->gq > output_calls[best]->gq ||
                (output_calls[i]->gq == output_calls[best]->gq &&
                 output_calls[i]->candidate.alternate_count >
                    output_calls[best]->candidate.alternate_count))
                best = i;
        annotations = output_calls[best]->annotations;
        // HaplotypeCaller annotates the final VariantContext after joint AF
        // has selected its concrete ALT set.  Keep the original multi-ALT
        // likelihood matrix for BestAllele, but allow the Host annotation
        // boundary to reject a read whose winning ALT was removed.  The
        // owner lookup mirrors joint PL/AD because an aggregate streamed
        // Result deliberately leaves each ragged PairHMM collection in its
        // AssemblyRegion owner.  This never moves PairHMM or AF work out of
        // the existing Kokkos kernels.
        if (annotation_reads != nullptr &&
            (max_allele_subset_changed || joint_af_allele_subset_changed)) {
            std::vector<const fastgatk::calling::AssemblyCandidate*> output_candidates;
            output_candidates.reserve(output_calls.size());
            for (const auto* call : output_calls) {
                if (call == nullptr) {
                    output_candidates.clear();
                    break;
                }
                output_candidates.push_back(&call->candidate);
            }
            const auto annotations_from_owner = [&]()
                -> std::optional<fastgatk::calling::GenotypeCall::Annotations> {
                if (output_candidates.empty()) return std::nullopt;
                // A partitioned result has flattened candidate metadata but
                // each compact likelihood-row ordinal remains owned by its
                // AssemblyRegion Result.  Prefer that owner so its stable
                // source-record map resolves against `annotation_reads`.
                for (const auto& owner : result.assembly_region_likelihood_results) {
                    if (owner == nullptr) continue;
                    if (const auto owned = fastgatk::calling::calculate_output_variant_annotations(
                            *annotation_reads, *owner, output_candidates, site_qual,
                            informative_read_overlap_margin, false,
                            max_allele_subset_changed && !joint_af_allele_subset_changed);
                        owned.has_value())
                        return owned;
                }
                if (const auto direct = fastgatk::calling::calculate_output_variant_annotations(
                        *annotation_reads, result, output_candidates, site_qual,
                        informative_read_overlap_margin, false,
                        max_allele_subset_changed && !joint_af_allele_subset_changed);
                    direct.has_value())
                    return direct;
                return std::nullopt;
            };
            if (const auto final_annotations = annotations_from_owner();
                final_annotations.has_value())
                annotations = *final_annotations;
        }
        if (joint_qd.has_value()) annotations.qd = *joint_qd;
        std::vector<int> ac(output_calls.size(), 0);
        std::vector<const fastgatk::calling::AssemblyCandidate*> mle_candidates;
        mle_candidates.reserve(output_calls.size());
        for (const auto* call : output_calls) {
            if (call == nullptr)
                throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: missing call in MLE candidate group");
            mle_candidates.push_back(&call->candidate);
        }
        auto mle_ac = estimate_mle_allele_counts(
            result, mle_candidates, output_genotype_vectors, pl, sample_ploidy);
        // GenotypingEngine retains MLE counts from the original joint
        // AFCalculationResult, then projects only the published ALT entries.
        if (joint_frequency.has_value() &&
            joint_frequency->integer_allele_counts.size() == frequency_allele_count) {
            mle_ac.clear();
            mle_ac.reserve(output_alt_indices.size());
            for (const auto alt : output_alt_indices)
                mle_ac.push_back(joint_frequency->integer_allele_counts[alt + 1]);
        }
        std::size_t best_genotype = 0;
        int best_gq = 0;
        // Keep genotype selection on the shared Kokkos API for both diploid
        // and arbitrary-ploidy VCF output.  The Host owns only the ragged
        // per-read likelihood-to-PL materialization above.
        std::vector<const fastgatk::calling::AssemblyCandidate*> prior_candidates;
        prior_candidates.reserve(output_calls.size());
        for (const auto* call : output_calls) {
            if (call == nullptr)
                throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: missing call in VCF candidate group");
            prior_candidates.push_back(&call->candidate);
        }
        const auto genotype_priors = genotype_priors_for_group(
            result, prior_candidates, sample_ploidy, false);
        if (genotype_priors.size() != output_genotype_vectors.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: genotype prior width mismatch");
        const auto derived = result.genotype_priors_used
            ? fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                pl, 1, static_cast<int>(output_allele_count), sample_ploidy, genotype_priors)
            : fastgatk::kernels::derive_genotype_gt_gq_kokkos(
                pl, 1, static_cast<int>(output_allele_count), sample_ploidy);
        if (derived.alleles.size() != static_cast<std::size_t>(sample_ploidy))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Kokkos GT is incomplete");
        for (std::size_t genotype = 0; genotype < output_genotype_vectors.size(); ++genotype) {
            if (output_genotype_vectors[genotype] ==
                std::vector<int>(derived.alleles.begin(), derived.alleles.end())) {
                best_genotype = genotype;
                break;
            }
        }
        best_gq = derived.gq.empty() ? 0 : std::clamp(derived.gq.front(), 0, 99);
        annotations.qd = gatk_qual_by_depth(annotations.qd, qd_random);
        std::ostringstream genotype_text;
        for (std::size_t copy = 0;
             copy < output_genotype_vectors[best_genotype].size(); ++copy) {
            if (copy != 0) genotype_text << '/';
            const auto allele = output_genotype_vectors[best_genotype][copy];
            genotype_text << allele;
            if (allele > 0) ++ac[static_cast<std::size_t>(allele - 1)];
        }
        // GenotypeGivenAlleles keeps a feature ALT even when every sample is
        // hom-ref. In that one source-defined output shape GATK serializes
        // zero AF/MLEAF with two places and omits QD (there is no called ALT
        // depth to normalize). Ordinary HC rows retain their established
        // Number=A formatter and annotation policy.
        const bool forced_hom_ref_feature = std::all_of(
            output_calls.begin(), output_calls.end(), [](const auto* call) {
                return call != nullptr && call->candidate.forced_by_alleles_feature;
            }) && std::all_of(ac.begin(), ac.end(), [](const auto count) { return count == 0; });
        const auto append_annotation = [&](const char* name, double value, int precision) {
            const auto text = annotation_text(value, precision);
            if (!text.empty()) out << ';' << name << '=' << text;
        };
        out << chrom << '\t' << (position + 1) << "\t.\t" << ref << '\t';
        for (std::size_t i = 0; i < output_calls.size(); ++i) {
            if (i != 0) out << ',';
            out << candidate_alternate(output_calls[i]->candidate);
        }
        // Match the standard HaplotypeCaller record contract: AD is a FORMAT
        // field (not an INFO field), unfiltered calls use '.', and INFO
        // values follow GATK's registration order.  This makes the textual
        // boundary stable for consumers that stream without reordering keys.
        out << '\t' << qual_text(site_qual) << "\t.\tAC=";
        for (std::size_t i = 0; i < ac.size(); ++i) {
            if (i != 0) out << ',';
            out << ac[i];
        }
        out << ";AF=";
        for (std::size_t i = 0; i < ac.size(); ++i) {
            if (i != 0) out << ',';
            const auto frequency = static_cast<double>(ac[i]) / static_cast<double>(an);
            out << (forced_hom_ref_feature && ac[i] == 0
                ? std::string("0.00") : format_allele_frequency(frequency));
        }
        out << ";AN=" << an;
        append_annotation("BaseQRankSum", annotations.base_q_rank_sum, 3);
        out << ";DP=" << variant_depth;
        // GATK's ExcessHet annotation is a diploid/cohort exact-test field.
        // It is not populated for arbitrary-ploidy HC records (for example
        // the tetra-ploid multi-ALT oracle), so do not fabricate a zero
        // value that would make an otherwise compatible record diverge.
        if (sample_ploidy == 2)
            out << ";ExcessHet=0.0000";
        append_annotation("FS", annotations.fs, 3);
        out << ";MLEAC=";
        for (std::size_t i = 0; i < ac.size(); ++i) {
            if (i != 0) out << ',';
            out << (i < mle_ac.size() ? mle_ac[i] : 0);
        }
        out << ";MLEAF=";
        for (std::size_t i = 0; i < ac.size(); ++i) {
            if (i != 0) out << ',';
            const auto count = i < mle_ac.size() ? mle_ac[i] : 0;
            const auto frequency = static_cast<double>(count) / static_cast<double>(an);
            out << (forced_hom_ref_feature && count == 0
                ? std::string("0.00") : format_allele_frequency(frequency));
        }
        append_annotation("MQ", annotations.mq, 2);
        append_annotation("MQRankSum", annotations.mq_rank_sum, 3);
        if (!forced_hom_ref_feature)
            append_annotation("QD", annotations.qd, 2);
        append_annotation("ReadPosRankSum", annotations.read_pos_rank_sum, 3);
        append_annotation("SOR", annotations.sor, 3);
        if (std::any_of(output_calls.begin(), output_calls.end(), [](const auto* call) {
                return call != nullptr && call->candidate.flow_hmer_collapsed;
            }))
            out << ";XC=1";
        if (sites_only_vcf_output) {
            out << '\n';
            continue;
        }
        out << "\tGT:AD:DP:GQ:PL\t" << genotype_text.str() << ':';
        for (std::size_t i = 0; i < ad.size(); ++i) {
            if (i != 0) out << ',';
            out << ad[i];
        }
        out << ':' << depth << ':' << best_gq << ':';
        for (std::size_t i = 0; i < pl.size(); ++i) {
            if (i != 0) out << ',';
            out << pl[i];
        }
        out << '\n';
    }
    return out.str();
}

// GATK annotates a record from the AssemblyRegion whose retained PairHMM read
// set produced it.  A partitioned native Result can list the same allele in
// more than one AssemblyRegion owner, and an owner that never requested
// PairHMM for that allele keeps the UINT32_MAX context-ordinal sentinel.
// Annotating from such an owner empties `context_mapping_evidence` (its own
// likelihood rows are all -inf) and drops the MQ gate onto the raw-overlap
// predicate, so the evidence - and therefore the published annotations -
// depend on the `-L` window start.  Resolve every published allele to its
// identical-allele twin inside the owner and require that twin to carry a real
// PairHMM context ordinal, so both the likelihood-row lookups and the
// context-ordinal lookup are served by the twin.
bool owner_has_pairhmm_context(
    const fastgatk::calling::Result& owner,
    const std::vector<const fastgatk::calling::AssemblyCandidate*>& alleles) {
    if (alleles.empty()) return false;
    for (const auto* allele : alleles) {
        if (allele == nullptr) return false;
        std::size_t twin = std::numeric_limits<std::size_t>::max();
        for (std::size_t index = 0; index < owner.candidates.size(); ++index) {
            const auto& candidate = owner.candidates[index];
            if (candidate.tid == allele->tid && candidate.position == allele->position &&
                candidate.reference == allele->reference &&
                candidate.alternate == allele->alternate &&
                candidate.reference_allele == allele->reference_allele &&
                candidate.alternate_allele == allele->alternate_allele) {
                twin = index;
                break;
            }
        }
        if (twin == std::numeric_limits<std::size_t>::max()) return false;
        if (twin >= owner.likelihood_candidate_read_context_ordinals.size()) return false;
        if (owner.likelihood_candidate_read_context_ordinals[twin] ==
            std::numeric_limits<std::uint32_t>::max())
            return false;
    }
    return true;
}

std::string gvcf(const fastgatk::io::HtsReader& reader,
                 fastgatk::calling::Result& result,
                 int sample_ploidy,
                 bool bp_resolution,
                 std::size_t min_alt_support = 1,
                 bool add_output_vcf_command_line = true,
                 const std::vector<int>& requested_gq_bands = {},
                 bool sites_only_vcf_output = false,
                 bool floor_blocks = false,
                 const fastgatk::io::ReadBatch* annotation_reads = nullptr,
                 std::uint32_t informative_read_overlap_margin = 2,
                 std::size_t max_alternate_alleles = 6) {
    std::ostringstream out;
    // Keep the gVCF schema compatible with HTSJDK's 4.6.2.0 header.  The
    // native-only status lines remain explicit provenance, while unused
    // optional fields (PGT/PID/PS/SB and InbreedingCoeff/RAW_MQandDP) are
    // declared exactly as GATK does so downstream Combine/GenotypeGVCFs can
    // consume either producer without schema surprises.
    out << "##fileformat=VCFv4.2\n"
        << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele not already represented at this location by REF and ALT\">\n"
        << "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
        << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n"
        << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    // BP_RESOLUTION emits single-site records, so GATK does not declare the
    // reference-block-only END/MIN_DP header fields in that mode.
    if (!bp_resolution)
        out << "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">\n";
    out
        << "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles\">\n"
        << "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">\n"
        << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n"
        << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">\n"
        << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">\n";
    if (!bp_resolution) {
        // HTSJDK sorts these symbolic GVCFBlock keys lexicographically, so
        // emit the same stable order instead of numeric lower-bound order.
        // GATK's -GQB values are upper bounds; the implicit lower bound is
        // zero and the terminal 100 bound is added when omitted.
        std::vector<int> gq_bands = requested_gq_bands;
        if (gq_bands.empty()) {
            for (int upper = 1; upper <= 60; ++upper) gq_bands.push_back(upper);
            gq_bands.insert(gq_bands.end(), {70, 80, 90, 99, 100});
        } else if (gq_bands.back() != 100) {
            gq_bands.push_back(100);
        }
        std::vector<std::string> blocks;
        int lower = 0;
        for (const auto upper : gq_bands) {
            if (upper <= lower || upper < 1 || upper > 100) continue;
            blocks.push_back("##GVCFBlock" + std::to_string(lower) + '-' +
                std::to_string(upper) + "=minGQ=" + std::to_string(lower) +
                "(inclusive),maxGQ=" + std::to_string(upper) + "(exclusive)");
            lower = upper;
            if (upper == 100) break;
        }
        std::sort(blocks.begin(), blocks.end());
        for (const auto& block : blocks) out << block << '\n';
    }
    out << "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n";
    if (!bp_resolution)
        out << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n";
    out
        << "##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Phred-scaled p-value for exact test of excess heterozygosity\">\n"
        << "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">\n"
        << "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">\n"
        << "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
        << "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=\"Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.\">\n"
        << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n";
    const auto& header = reader.header();
    for (std::size_t i = 0; i < header.contigs.size(); ++i)
        out << "##contig=<ID=" << header.contigs[i]
            << ",length=" << header.contig_lengths[i]
            << (i < header.contig_assemblies.size() && !header.contig_assemblies[i].empty()
                    ? ",assembly=" + header.contig_assemblies[i] : std::string{})
            << ">\n";
    // GATK adds the source/provenance pair together for a gVCF as well.  Keep
    // any native implementation status in the manifest rather than adding
    // non-GATK header records.
    if (add_output_vcf_command_line)
        out << "##source=HaplotypeCaller\n";
    if (add_output_vcf_command_line)
        out << "##GATKCommandLine=<ID=HaplotypeCaller,Version=fastgatk-native,"
               "CommandLine=\"fastgatk HaplotypeCaller\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    if (!sites_only_vcf_output)
        out << "\tFORMAT\t" << output_sample_name(header);
    out << "\n";
    const auto chrom_for = [&](std::int32_t tid) {
        const auto index = static_cast<std::size_t>(std::max<std::int32_t>(0, tid));
        return index < header.contigs.size() ? header.contigs[index] : std::to_string(tid);
    };
    const auto annotation_text = [](double value, int precision) {
        if (!std::isfinite(value)) return std::string{};
        if (std::abs(value) < 0.5 * std::pow(10.0, -precision)) value = 0.0;
        std::ostringstream output;
        output << std::fixed << std::setprecision(precision) << value;
        return output.str();
    };
    const auto qual_text = [](double value) {
        if (!std::isfinite(value)) return std::string(".");
        std::ostringstream output;
        output << std::fixed << std::setprecision(2)
               << std::clamp(value, 0.0, 9999999.0);
        return output.str();
    };
    const auto floored_block_gq = [&](int gq) {
        if (!floor_blocks || bp_resolution) return gq;
        // --floor-blocks uses the lower edge of the GQ band, rather than
        // merely dropping MIN_DP/PL.  Build the same implicit partition used
        // by the GVCFBlock header above so custom -GQB values stay aligned.
        std::vector<int> gq_bands = requested_gq_bands;
        if (gq_bands.empty()) {
            for (int upper = 1; upper <= 60; ++upper) gq_bands.push_back(upper);
            gq_bands.insert(gq_bands.end(), {70, 80, 90, 99, 100});
        } else if (gq_bands.back() != 100) {
            gq_bands.push_back(100);
        }
        int lower = 0;
        for (const auto upper : gq_bands) {
            if (upper <= lower || upper < 1 || upper > 100) continue;
            if (gq < upper) return lower;
            lower = upper;
        }
        return lower;
    };
    const auto write_block = [&](const fastgatk::calling::ReferenceBlock& block) {
        const auto position = block.start + 1;
        const auto block_pl = block.genotype_pl.empty()
            ? std::vector<std::int32_t>(block.pl.begin(), block.pl.end())
            : block.genotype_pl;
        const auto block_ploidy = block_pl.empty()
            ? 2 : static_cast<int>(block_pl.size() - 1);
        out << chrom_for(block.tid) << '\t' << position << "\t.\t"
            << static_cast<char>(block.reference) << "\t<NON_REF>\t.\t.\t";
        if (!bp_resolution)
            out << "END=" << (block.end + 1);
        else
            out << ".";
        if (sites_only_vcf_output) {
            out << '\n';
            return;
        }
        if (!bp_resolution)
            out << (floor_blocks ? "\tGT:DP:GQ\t" : "\tGT:DP:GQ:MIN_DP:PL\t");
        else
            out << "\tGT:AD:DP:GQ:PL\t";
        for (int copy = 0; copy < block_ploidy; ++copy) {
            if (copy != 0) out << '/';
            out << "0";
        }
        if (bp_resolution) {
            out << ':' << block.reference_count << ',' << block.non_ref_count
                << ':' << block.depth << ':' << static_cast<int>(block.gq) << ':';
        } else {
            out << ':' << block.depth << ':' << floored_block_gq(static_cast<int>(block.gq));
            if (!floor_blocks) out << ':' << block.min_depth << ':';
        }
        if (!floor_blocks || bp_resolution) {
            for (std::size_t index = 0; index < block_pl.size(); ++index) {
                if (index != 0) out << ',';
                out << block_pl[index];
            }
        }
        out << '\n';
    };
    struct CandidateGroup {
        std::int32_t tid = -1;
        std::int32_t position = -1;
        std::string reference;
        std::vector<const fastgatk::calling::AssemblyCandidate*> candidates;
        // GATK's --max-alternate-alleles pass is applied after <NON_REF> is
        // added, so its surviving Number=G likelihoods are a subset of the
        // original symbolic matrix, not a fresh median recomputation over
        // only the kept ALTs.  Cache that Host structural projection here;
        // PairHMM itself remains the source of the original Kokkos matrix.
        bool max_alt_subset = false;
        // `*` is a proper alternate allele for GenotypingEngine's cap (only
        // <NON_REF> is automatically retained). This records whether the
        // GL-based selection kept the upstream spanning-deletion row.
        bool max_alt_retains_spanning_deletion = false;
        std::vector<int> max_alt_symbolic_pl;
    };
    std::vector<CandidateGroup> candidate_groups;
    for (const auto& candidate : result.candidates) {
        // An ordinary hom-ref/no-call candidate is represented by a
        // reference-confidence block. A retained EventMap candidate recorded
        // by the calling pipeline is a concrete GATK gVCF row even when its
        // final GT is 0/0.  The pipeline has already rebuilt reference blocks
        // using GenotypingEngine's full REF/ALT/<NON_REF> AF decision, so the
        // rendered group is filtered against that authoritative block span
        // below rather than guessing from this provisional biallelic call.
        const auto called = std::find_if(result.calls.begin(), result.calls.end(),
            [&](const auto& call) {
                return call.candidate.tid == candidate.tid &&
                       call.candidate.position == candidate.position &&
                       candidate_alternate(call.candidate) == candidate_alternate(candidate);
            });
        const auto gvcf_hom_ref = std::find_if(
            result.gvcf_hom_ref_calls.begin(), result.gvcf_hom_ref_calls.end(),
            [&](const auto& call) {
                return call.candidate.tid == candidate.tid &&
                       call.candidate.position == candidate.position &&
                       candidate_alternate(call.candidate) == candidate_alternate(candidate);
            });
        const auto indel = candidate_reference(candidate).size() !=
                           candidate_alternate(candidate).size();
        const bool retain_gvcf_hom_ref = gvcf_hom_ref != result.gvcf_hom_ref_calls.end();
        const bool retain_called = called != result.calls.end();
        if (!retain_called && !retain_gvcf_hom_ref &&
            (!indel || candidate.alternate_count < min_alt_support)) continue;
        const auto reference = candidate_reference(candidate);
        if (candidate_groups.empty() || candidate_groups.back().tid != candidate.tid ||
            candidate_groups.back().position != candidate.position ||
            candidate_groups.back().reference != reference) {
            candidate_groups.push_back(CandidateGroup{candidate.tid, candidate.position,
                                                      reference, {}, false, false, {}});
        }
        candidate_groups.back().candidates.push_back(&candidate);
    }
    // ReferenceConfidenceModel emits a block at every locus that did not
    // survive GenotypingEngine as a VariantContext.  In particular, a
    // provisional biallelic EventMap call can become monomorphic after its
    // <NON_REF> likelihood is incorporated.  The calling pipeline rebuilds
    // blocks after that Kokkos PL/AF decision; use their span as the single
    // Host rendering authority so a discarded candidate cannot reappear as a
    // concrete gVCF record.
    candidate_groups.erase(std::remove_if(candidate_groups.begin(), candidate_groups.end(),
        [&](const CandidateGroup& group) {
            return std::any_of(result.reference_blocks.begin(), result.reference_blocks.end(),
                [&](const auto& block) {
                    return block.tid == group.tid && block.start <= group.position &&
                           group.position <= block.end;
                });
        }), candidate_groups.end());
    // Determine `*` ownership before concrete ALT subsetting as well as for
    // final rendering. GATK has already inserted a spanning-deletion allele
    // into the EventMap allele mapper when its subsequent GenotypingEngine
    // --max-alternate-alleles pass examines the full symbolic GL matrix. The
    // interval depends on the group's REF span, not which concrete ALT later
    // survives, so this pre-pass is stable across the structural subset.
    std::vector<std::uint8_t> gvcf_group_covered_before_alt_subsetting(
        candidate_groups.size(), 0U);
    std::vector<std::int32_t> pre_subset_deletion_ends;
    std::int32_t pre_subset_deletion_tid = -1;
    for (std::size_t group_index = 0; group_index < candidate_groups.size(); ++group_index) {
        const auto& group = candidate_groups[group_index];
        if (group.tid != pre_subset_deletion_tid) {
            pre_subset_deletion_ends.clear();
            pre_subset_deletion_tid = group.tid;
        }
        pre_subset_deletion_ends.erase(
            std::remove_if(pre_subset_deletion_ends.begin(), pre_subset_deletion_ends.end(),
                [&](const auto end) { return end < group.position; }),
            pre_subset_deletion_ends.end());
        gvcf_group_covered_before_alt_subsetting[group_index] = std::any_of(
            pre_subset_deletion_ends.begin(), pre_subset_deletion_ends.end(),
            [&](const auto end) { return end >= group.position; }) ? 1U : 0U;
        if (group.reference.size() > 1U) {
            const auto end = static_cast<std::int64_t>(group.position) +
                static_cast<std::int64_t>(group.reference.size()) - 1;
            if (end <= std::numeric_limits<std::int32_t>::max())
                pre_subset_deletion_ends.push_back(static_cast<std::int32_t>(end));
        }
    }
    // Reference-confidence genotyping adds <NON_REF> before it reaches
    // GenotypingEngine.calculateGenotypes().  Resolve the full symbolic
    // PairHMM matrix now, use the source GL selection rule on proper ALTs,
    // and retain the subsetted matrix for serialization/AF below.  This is
    // intentionally before the gVCF deletion-state and phase writers: both
    // observe the same final concrete ALT list as GATK.
    for (std::size_t group_index = 0; group_index < candidate_groups.size(); ++group_index) {
        auto& group = candidate_groups[group_index];
        const bool has_spanning_deletion =
            group_index < gvcf_group_covered_before_alt_subsetting.size() &&
            gvcf_group_covered_before_alt_subsetting[group_index] != 0;
        const auto proper_alt_count = group.candidates.size() +
            (has_spanning_deletion ? 1U : 0U);
        if (proper_alt_count <= max_alternate_alleles) continue;
        std::vector<int> full_pl;
        const bool include_spanning_deletion = has_spanning_deletion;
        const auto expected_width = genotype_width(
            static_cast<int>(group.candidates.size() + 2U +
                             (include_spanning_deletion ? 1U : 0U)), sample_ploidy);
        const auto load_symbolic_pl = [&](const fastgatk::calling::Result& owner) {
            return joint_candidate_pl_with_nonref(
                owner, group.candidates, sample_ploidy, full_pl, 0,
                include_spanning_deletion) &&
                full_pl.size() == expected_width;
        };
        bool loaded = load_symbolic_pl(result);
        if (!loaded) {
            for (const auto& owner : result.assembly_region_likelihood_results) {
                if (owner != nullptr && load_symbolic_pl(*owner)) {
                    loaded = true;
                    break;
                }
            }
        }
        if (!loaded)
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: --max-alternate-alleles gVCF requires reference-backed PairHMM likelihoods");
        const auto retained = gatk_most_likely_alt_indices(
            full_pl, static_cast<int>(group.candidates.size() + 2U +
                                      (include_spanning_deletion ? 1U : 0U)), sample_ploidy,
            proper_alt_count, max_alternate_alleles);
        std::vector<std::size_t> new_to_old{0U};
        new_to_old.reserve(retained.size() + 2U +
                           (include_spanning_deletion ? 1U : 0U));
        const auto original_concrete_count = group.candidates.size();
        std::vector<const fastgatk::calling::AssemblyCandidate*> retained_candidates;
        retained_candidates.reserve(retained.size());
        for (const auto alt : retained) {
            if (alt < original_concrete_count) {
                new_to_old.push_back(alt + 1U);
                retained_candidates.push_back(group.candidates[alt]);
            } else if (include_spanning_deletion && alt == original_concrete_count) {
                // REF, concrete ALTs, `*`, <NON_REF>.
                new_to_old.push_back(original_concrete_count + 1U);
                group.max_alt_retains_spanning_deletion = true;
            }
        }
        // The old symbolic <NON_REF> allele is always final.
        new_to_old.push_back(original_concrete_count + 1U +
                             (include_spanning_deletion ? 1U : 0U));
        group.max_alt_symbolic_pl = remap_pl_for_allele_subset(
            full_pl, static_cast<int>(original_concrete_count + 2U +
                                      (include_spanning_deletion ? 1U : 0U)), sample_ploidy,
            new_to_old);
        group.candidates = std::move(retained_candidates);
        group.max_alt_subset = true;
    }
    // GenotypingEngine records previously emitted deletion intervals before it
    // visits the next event.  In reference-confidence mode every concrete
    // candidate carries `<NON_REF>`; for a multi-base reference that symbolic
    // one-base allele itself records a deletion interval.  Consequently an
    // EventMap `*` at the following covered locus is retained in gVCF, while
    // the same `*` is spurious and omitted from an ordinary VCF.  Compute this
    // small ordered Host state machine once from the already selected output
    // groups; PairHMM/PL/AF work remains in the Kokkos kernels below.
    std::vector<std::uint8_t> gvcf_group_covered_by_upstream_deletion(
        candidate_groups.size(), 0U);
    std::vector<std::int32_t> upstream_deletion_ends;
    std::int32_t upstream_deletion_tid = -1;
    for (std::size_t group_index = 0; group_index < candidate_groups.size(); ++group_index) {
        const auto& group = candidate_groups[group_index];
        if (group.tid != upstream_deletion_tid) {
            upstream_deletion_ends.clear();
            upstream_deletion_tid = group.tid;
        }
        upstream_deletion_ends.erase(
            std::remove_if(upstream_deletion_ends.begin(), upstream_deletion_ends.end(),
                [&](const auto end) { return end < group.position; }),
            upstream_deletion_ends.end());
        gvcf_group_covered_by_upstream_deletion[group_index] = std::any_of(
            upstream_deletion_ends.begin(), upstream_deletion_ends.end(),
            [&](const auto end) { return end >= group.position; }) ? 1U : 0U;
        // `<NON_REF>` has one reference base in HTSJDK's allele length model.
        // It is part of every emitted gVCF candidate record, so any longer REF
        // records an interval through its final reference base.
        if (group.reference.size() > 1U) {
            const auto span = static_cast<std::int64_t>(group.reference.size()) - 1;
            const auto end = static_cast<std::int64_t>(group.position) + span;
            if (end <= std::numeric_limits<std::int32_t>::max())
                upstream_deletion_ends.push_back(static_cast<std::int32_t>(end));
        }
    }
    // GATK's AssemblyBasedCallerUtils.phaseCalls() phases only biallelic,
    // concrete calls and derives relationships from EventMap haplotype
    // membership.  Keep that graph bookkeeping on the Host: the PairHMM
    // likelihoods and GT/GQ calculation remain the existing Kokkos work.
    struct PhysicalPhase {
        bool enabled = false;
        bool phase_01 = true;
        std::string id;
        int phase_set = 0;
    };
    std::vector<PhysicalPhase> physical_phases(candidate_groups.size());
    if (std::getenv("FASTGATK_DEBUG_PHASE") != nullptr) {
        std::cerr << "[FASTGATK_PHASE_INPUT] candidates=" << result.candidates.size()
                  << " candidate_groups=" << candidate_groups.size()
                  << " haplotype_rows=" << result.somatic_candidate_haplotype_indices.size()
                  << " calls=" << result.calls.size() << '\n';
    }
    if (sample_ploidy == 2) {
        std::vector<const fastgatk::calling::Result*> phase_contexts;
        phase_contexts.push_back(&result);
        for (const auto& owner : result.assembly_region_likelihood_results)
            if (owner != nullptr) phase_contexts.push_back(owner.get());
        for (const auto* phase_context : phase_contexts) {
            if (phase_context == nullptr ||
                phase_context->somatic_candidate_haplotype_indices.size() !=
                    phase_context->candidates.size())
                continue;
        std::vector<std::set<std::uint32_t>> haplotype_map(candidate_groups.size());
        for (std::size_t group_index = 0; group_index < candidate_groups.size(); ++group_index) {
            const auto& group = candidate_groups[group_index];
            if (group.candidates.size() != 1) continue;
            const auto* candidate = group.candidates.front();
            const auto called = std::find_if(result.calls.begin(), result.calls.end(),
                [&](const auto& call) {
                    return call.candidate.tid == candidate->tid &&
                           call.candidate.position == candidate->position &&
                           candidate_alternate(call.candidate) == candidate_alternate(*candidate);
            });
            const auto gvcf_hom_ref_called = std::find_if(
                result.gvcf_hom_ref_calls.begin(), result.gvcf_hom_ref_calls.end(),
                [&](const auto& call) {
                    return call.candidate.tid == candidate->tid &&
                           call.candidate.position == candidate->position &&
                           candidate_alternate(call.candidate) == candidate_alternate(*candidate);
                });
            // In reference-confidence mode GenotypingEngine returns concrete
            // EventMap calls with GT=0/0 as well. They still contribute their
            // haplotypes to `calledHaplotypes` before phaseCalls(), so leaving
            // them out here can make an otherwise unphasable set appear to
            // phase cleanly.
            if (called == result.calls.end() &&
                gvcf_hom_ref_called == result.gvcf_hom_ref_calls.end()) continue;
            const auto owner_candidate = std::find_if(
                phase_context->candidates.begin(), phase_context->candidates.end(),
                [&](const auto& item) {
                    return item.tid == candidate->tid && item.position == candidate->position &&
                           candidate_reference(item) == candidate_reference(*candidate) &&
                           candidate_alternate(item) == candidate_alternate(*candidate);
                });
            if (owner_candidate == phase_context->candidates.end()) continue;
            const auto owner_call = std::find_if(
                phase_context->calls.begin(), phase_context->calls.end(),
                [&](const auto& item) {
                    return item.candidate.tid == candidate->tid &&
                           item.candidate.position == candidate->position &&
                           candidate_alternate(item.candidate) == candidate_alternate(*candidate);
                });
            const auto owner_gvcf_hom_ref_call = std::find_if(
                phase_context->gvcf_hom_ref_calls.begin(),
                phase_context->gvcf_hom_ref_calls.end(),
                [&](const auto& item) {
                    return item.candidate.tid == candidate->tid &&
                           item.candidate.position == candidate->position &&
                           candidate_alternate(item.candidate) == candidate_alternate(*candidate);
                });
            if (owner_call == phase_context->calls.end() &&
                owner_gvcf_hom_ref_call == phase_context->gvcf_hom_ref_calls.end()) continue;
            const auto candidate_index = static_cast<std::size_t>(
                std::distance(phase_context->candidates.begin(), owner_candidate));
            const auto& members =
                phase_context->somatic_candidate_haplotype_indices[candidate_index];
            haplotype_map[group_index].insert(members.begin(), members.end());
            if (std::getenv("FASTGATK_DEBUG_PHASE") != nullptr) {
                std::cerr << "[FASTGATK_PHASE_EVENT] pos=" << (candidate->position + 1)
                          << " ref=" << candidate_reference(*candidate)
                          << " alt=" << candidate_alternate(*candidate) << " haps=";
                for (const auto haplotype : haplotype_map[group_index])
                    std::cerr << haplotype << ',';
                std::cerr << '\n';
            }
        }

        std::set<std::uint32_t> all_called_haplotypes;
        for (const auto& members : haplotype_map)
            all_called_haplotypes.insert(members.begin(), members.end());
        const auto total_available_haplotypes = all_called_haplotypes.size();
        std::vector<int> phase_group(candidate_groups.size(), -1);
        std::vector<bool> phase_01(candidate_groups.size(), true);
        int unique_counter = 0;
        bool unphasable = false;
        for (std::size_t i = 0; i + 1 < candidate_groups.size() && !unphasable; ++i) {
            const auto& call_haplotypes = haplotype_map[i];
            if (call_haplotypes.empty()) continue;
            const bool call_is_on_all_alt_haplotypes =
                call_haplotypes.size() == total_available_haplotypes;
            auto call_haplotypes_available_for_phasing = call_haplotypes;
            for (std::size_t j = i + 1; j < candidate_groups.size(); ++j) {
                const auto& comp_haplotypes = haplotype_map[j];
                if (comp_haplotypes.empty()) continue;
                const bool comp_is_on_all_alt_haplotypes =
                    comp_haplotypes.size() == total_available_haplotypes;
                const bool always_together = call_haplotypes == comp_haplotypes ||
                    (call_is_on_all_alt_haplotypes &&
                     std::includes(call_haplotypes_available_for_phasing.begin(),
                                   call_haplotypes_available_for_phasing.end(),
                                   comp_haplotypes.begin(), comp_haplotypes.end())) ||
                    comp_is_on_all_alt_haplotypes;
                if (always_together) {
                    if (phase_group[i] == -1) {
                        if (phase_group[j] != -1) {
                            unphasable = true;
                            break;
                        }
                        phase_group[i] = unique_counter;
                        phase_01[i] = true;
                        phase_group[j] = unique_counter;
                        phase_01[j] = true;
                        call_haplotypes_available_for_phasing.clear();
                        std::set_intersection(call_haplotypes.begin(), call_haplotypes.end(),
                                              comp_haplotypes.begin(), comp_haplotypes.end(),
                                              std::inserter(call_haplotypes_available_for_phasing,
                                                            call_haplotypes_available_for_phasing.end()));
                        ++unique_counter;
                    } else if (phase_group[j] == -1) {
                        phase_group[j] = phase_group[i];
                        phase_01[j] = phase_01[i];
                    }
                    continue;
                }
                const bool always_apart = call_haplotypes.size() + comp_haplotypes.size() ==
                    total_available_haplotypes &&
                    std::none_of(call_haplotypes.begin(), call_haplotypes.end(),
                        [&](const auto haplotype) {
                            return comp_haplotypes.count(haplotype) != 0;
                        });
                if (!always_apart) continue;
                if (phase_group[i] == -1) {
                    if (phase_group[j] != -1) {
                        unphasable = true;
                        break;
                    }
                    phase_group[i] = unique_counter;
                    phase_01[i] = true;
                    phase_group[j] = unique_counter;
                    phase_01[j] = false;
                    ++unique_counter;
                } else if (phase_group[j] == -1) {
                    phase_group[j] = phase_group[i];
                    phase_01[j] = !phase_01[i];
                }
            }
        }
        if (!unphasable) {
            for (int group_id = 0; group_id < unique_counter; ++group_id) {
                std::vector<std::size_t> members;
                for (std::size_t index = 0; index < phase_group.size(); ++index)
                    if (phase_group[index] == group_id) members.push_back(index);
                if (members.size() < 2) continue;
                const auto first = members.front();
                const auto& first_candidate = *candidate_groups[first].candidates.front();
                const auto id = std::to_string(candidate_groups[first].position + 1) + "_" +
                    candidate_groups[first].reference + "_" + candidate_alternate(first_candidate);
                const auto phase_set = candidate_groups[first].position + 1;
                for (const auto member : members) {
                    physical_phases[member].enabled = true;
                    physical_phases[member].phase_01 = phase_01[member];
                    physical_phases[member].id = id;
                    physical_phases[member].phase_set = phase_set;
                    if (std::getenv("FASTGATK_DEBUG_PHASE") != nullptr) {
                        std::cerr << "[FASTGATK_PHASE_OUTPUT] pos="
                                  << (candidate_groups[member].position + 1)
                                  << " phase=" << (phase_01[member] ? "0|1" : "1|0")
                                  << " pid=" << id << '\n';
                    }
                }
            }
        }
        }
    }
    const auto write_candidate = [&](const CandidateGroup& group,
                                     const PhysicalPhase& physical_phase,
                                     bool covered_by_upstream_deletion) {
        if (group.candidates.empty()) return;
        const bool has_indel = std::any_of(group.candidates.begin(), group.candidates.end(),
            [&](const auto* candidate) {
                return candidate_reference(*candidate).size() !=
                       candidate_alternate(*candidate).size();
            });
        const auto tid = static_cast<std::size_t>(std::max<std::int32_t>(0, group.tid));
        const std::string chrom = tid < header.contigs.size() ? header.contigs[tid]
                                                               : std::to_string(group.tid);
        const auto position = group.position + 1;
        const auto find_call = [&](const auto* candidate)
            -> const fastgatk::calling::GenotypeCall* {
            const auto ordinary = std::find_if(result.calls.begin(), result.calls.end(),
                [&](const fastgatk::calling::GenotypeCall& item) {
                    return item.candidate.tid == candidate->tid &&
                           item.candidate.position == candidate->position &&
                           candidate_alternate(item.candidate) == candidate_alternate(*candidate);
                });
            if (ordinary != result.calls.end()) return &*ordinary;
            const auto gvcf_hom_ref = std::find_if(
                result.gvcf_hom_ref_calls.begin(), result.gvcf_hom_ref_calls.end(),
                [&](const fastgatk::calling::GenotypeCall& item) {
                    return item.candidate.tid == candidate->tid &&
                           item.candidate.position == candidate->position &&
                           candidate_alternate(item.candidate) == candidate_alternate(*candidate);
                });
            return gvcf_hom_ref == result.gvcf_hom_ref_calls.end()
                ? nullptr : &*gvcf_hom_ref;
        };
        std::vector<const fastgatk::calling::GenotypeCall*> calls;
        calls.reserve(group.candidates.size());
        std::size_t best = 0;
        for (const auto* candidate : group.candidates) {
            const auto call = find_call(candidate);
            calls.push_back(call);
            if (calls.back() != nullptr &&
                (calls[best] == nullptr || calls.back()->gq > calls[best]->gq))
                best = calls.size() - 1;
        }
        if (sample_ploidy != 2) {
            std::vector<int> pl;
            // Partitioned callers intentionally flatten concrete calls, while
            // the symbolic <NON_REF> (and EventMap spanning-deletion) rows
            // remain in their AssemblyRegion PairHMM owner.  The old
            // arbitrary-ploidy branch only tried the flattened Result, so a
            // perfectly valid candidate could fail before any Kokkos PL work
            // ran.  Use the same owner lookup as the diploid writer; it is
            // Host-side ownership resolution and retains the common Kokkos
            // joint-genotype kernel for the numerical calculation.
            const auto joint_non_ref_from_owner = [&](bool include_spanning_deletion,
                                                      std::vector<int>& destination) {
                const auto allele_count = group.candidates.size() + 2U +
                    (include_spanning_deletion ? 1U : 0U);
                const auto expected_width = genotype_width(
                    static_cast<int>(allele_count), sample_ploidy);
                if (joint_candidate_pl_with_nonref(
                        result, group.candidates, sample_ploidy, destination, 0,
                        include_spanning_deletion) && destination.size() == expected_width)
                    return true;
                for (const auto& owner : result.assembly_region_likelihood_results) {
                    if (owner != nullptr && joint_candidate_pl_with_nonref(
                            *owner, group.candidates, sample_ploidy, destination, 0,
                            include_spanning_deletion) && destination.size() == expected_width)
                        return true;
                }
                return false;
            };
            // A spanning-deletion state may exist in the original gVCF
            // matrix yet be removed by --max-alternate-alleles.  Only render
            // (or size a Number=G matrix for) `*` when it survived that
            // source-compatible ALT selection.
            const bool span_is_selected = covered_by_upstream_deletion &&
                (!group.max_alt_subset || group.max_alt_retains_spanning_deletion);
            const bool include_spanning_deletion = span_is_selected &&
                joint_non_ref_from_owner(true, pl);
            // A single-ploidy GATK gVCF may fold an upstream EventMap allele
            // into a reference block, yet retain that allele as an internal
            // spanning-deletion state for the AF/QUAL calculation of the
            // overlapping concrete site.  It is not serialized as `*`, but
            // omitting it from AF changes chr20:10020680 from GATK's 89.07
            // to 140.04.  Keep Host ownership/output selection separate from
            // the Kokkos likelihood and AF work.
            const bool has_hidden_spanning_deletion = sample_ploidy == 1 &&
                !include_spanning_deletion && std::any_of(
                    result.candidates.begin(), result.candidates.end(),
                    [&](const auto& candidate) {
                        const auto ref_length = candidate_reference(candidate).size();
                        return candidate.tid == group.tid &&
                               candidate.position < group.position && ref_length > 1U &&
                               candidate.position + static_cast<std::int32_t>(ref_length) - 1 >=
                                   group.position;
                    });
            std::vector<int> hidden_spanning_pl;
            const bool hidden_spanning_deletion = has_hidden_spanning_deletion &&
                joint_non_ref_from_owner(true, hidden_spanning_pl);
            const auto concrete_allele_count = group.candidates.size() + 1U +
                (include_spanning_deletion ? 1U : 0U); // REF + ALTs + `*`
            const auto allele_count = concrete_allele_count + 1U; // plus <NON_REF>
            const auto frequency_allele_count = sample_ploidy == 1
                ? concrete_allele_count + (hidden_spanning_deletion ? 1U : 0U)
                : allele_count;
            const auto genotype_vectors = enumerate_genotypes(
                static_cast<int>(allele_count), sample_ploidy);
            const auto concrete_genotype_vectors = enumerate_genotypes(
                static_cast<int>(concrete_allele_count), sample_ploidy);
            bool joint_non_ref = include_spanning_deletion ||
                joint_non_ref_from_owner(false, pl);
            if (group.max_alt_subset) {
                pl = group.max_alt_symbolic_pl;
                joint_non_ref = true;
            }
            if (genotype_vectors.empty() || !joint_non_ref ||
                pl.size() != genotype_vectors.size())
                throw std::runtime_error(
                    "BACKEND_UNAVAILABLE: polyploid gVCF candidate requires reference-backed PairHMM likelihoods");
            const auto genotype_priors = genotype_priors_for_group(
                result, group.candidates, sample_ploidy, false, include_spanning_deletion);
            if (genotype_priors.size() != concrete_genotype_vectors.size())
                throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: gVCF genotype prior width mismatch");
            // `<NON_REF>` is emitted as an evidence envelope, never as the
            // called allele at a concrete candidate site.  Assign GT/GQ from
            // the concrete PL prefix; retain the full Number=G matrix in the
            // record for downstream GenotypeGVCFs.
            const auto concrete_width = concrete_genotype_vectors.size();
            if (pl.size() < concrete_width)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: polyploid concrete PL prefix is incomplete");
            const std::vector<int> concrete_pl(pl.begin(), pl.begin() + concrete_width);
            const auto derived = result.genotype_priors_used
                ? fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                    concrete_pl, 1, static_cast<int>(concrete_allele_count), sample_ploidy, genotype_priors)
                : fastgatk::kernels::derive_genotype_gt_gq_kokkos(
                    concrete_pl, 1, static_cast<int>(concrete_allele_count), sample_ploidy);
            if (derived.alleles.size() != static_cast<std::size_t>(sample_ploidy))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: polyploid gVCF GT is incomplete");
            std::size_t best_genotype = 0;
            for (std::size_t genotype = 0; genotype < genotype_vectors.size(); ++genotype) {
                if (genotype_vectors[genotype] == std::vector<int>(derived.alleles.begin(), derived.alleles.end())) {
                    best_genotype = genotype;
                    break;
                }
            }
            const int best_gq = derived.gq.empty() ? 0 : std::clamp(derived.gq.front(), 0, 99);
            // A gVCF candidate's QUAL and MLEAC are both defined on the full
            // symbolic matrix.  In particular, when `*` is retained in a
            // triploid record it is a real indel allele for the prior/AF
            // calculation, rather than a display-only ALT.  Keep the bounded
            // posterior in the existing Kokkos AF kernel; this Host code only
            // materializes the GATK allele ordering REF, concrete ALT(s), *,
            // <NON_REF>.
            std::optional<fastgatk::kernels::AlleleFrequencyResult> gvcf_frequency;
            if (result.genotype_heterozygosity_stdev > 0.0 &&
                std::isfinite(result.genotype_heterozygosity_stdev)) {
                const auto ref_pseudocount = result.genotype_snp_heterozygosity /
                    (result.genotype_heterozygosity_stdev *
                     result.genotype_heterozygosity_stdev);
                if (ref_pseudocount > 0.0 && std::isfinite(ref_pseudocount)) {
                    // GATK's haploid GenotypingEngine obtains site QUAL and
                    // MLE from its concrete EventMap AF calculation.  It adds
                    // <NON_REF> only when serializing the gVCF likelihood
                    // envelope.  Including that symbolic row in the haploid
                    // AF matrix changes chr20:10020680 from 89.07 to 140.04
                    // despite an otherwise identical GT/PL.  Keep the full
                    // Number=G matrix for output, but hand the existing
                    // concrete PL prefix to the same Kokkos AF kernel for
                    // this source-defined boundary.
                    std::vector<int> hidden_spanning_frequency_pl;
                    if (hidden_spanning_deletion) {
                        const auto hidden_width = genotype_width(
                            static_cast<int>(frequency_allele_count), sample_ploidy);
                        if (hidden_spanning_pl.size() < hidden_width)
                            throw std::runtime_error(
                                "OUTPUT_CONTRACT_FAILURE: hidden spanning-deletion PL prefix is incomplete");
                        hidden_spanning_frequency_pl.assign(
                            hidden_spanning_pl.begin(),
                            hidden_spanning_pl.begin() + static_cast<std::ptrdiff_t>(hidden_width));
                    }
                    const auto& frequency_pl = hidden_spanning_deletion
                        ? hidden_spanning_frequency_pl
                        : (sample_ploidy == 1 ? concrete_pl : pl);
                    std::vector<double> prior_pseudocounts(
                        frequency_allele_count,
                        result.genotype_snp_heterozygosity * ref_pseudocount);
                    prior_pseudocounts[0] = ref_pseudocount;
                    for (std::size_t alt = 0; alt < group.candidates.size(); ++alt) {
                        const auto& candidate = *group.candidates[alt];
                        const bool indel = candidate_reference(candidate).size() !=
                                           candidate_alternate(candidate).size();
                        prior_pseudocounts[alt + 1] =
                            (indel ? result.genotype_indel_heterozygosity
                                   : result.genotype_snp_heterozygosity) * ref_pseudocount;
                    }
                    if (include_spanning_deletion || hidden_spanning_deletion) {
                        const auto spanning = group.candidates.size() + 1U;
                        prior_pseudocounts[spanning] =
                            result.genotype_indel_heterozygosity * ref_pseudocount;
                    }
                    if (sample_ploidy != 1) {
                        // HTSJDK reports Allele.NON_REF_ALLELE.length() as
                        // zero because it is symbolic.  GATK's
                        // AlleleFrequencyCalculator therefore assigns it the
                        // indel pseudocount, not the SNP one.
                        prior_pseudocounts.back() =
                            result.genotype_indel_heterozygosity * ref_pseudocount;
                    }
                    try {
                        const auto calculated = fastgatk::kernels::calculate_allele_frequency_kokkos(
                            frequency_pl, 1, static_cast<int>(frequency_allele_count), sample_ploidy,
                            prior_pseudocounts, {}, {}, include_spanning_deletion
                                || hidden_spanning_deletion
                                ? static_cast<int>(group.candidates.size() + 1U) : -1);
                        if (calculated.samples_with_likelihoods != 0 &&
                            std::isfinite(calculated.qual))
                            gvcf_frequency = calculated;
                    } catch (const std::exception&) {
                        // A malformed symbolic row must not be converted into
                        // a synthetic confidence value.  The established
                        // candidate-call fallback remains available below.
                    }
                }
            }
            std::vector<int> output_mle_ac(allele_count - 1U, 0);
            if (gvcf_frequency.has_value() &&
                gvcf_frequency->integer_allele_counts.size() ==
                    (sample_ploidy == 1 ? frequency_allele_count : allele_count)) {
                const auto mle_allele_count = sample_ploidy == 1
                    ? concrete_allele_count : allele_count;
                for (std::size_t allele = 1; allele < mle_allele_count; ++allele)
                    output_mle_ac[allele - 1U] = std::clamp(
                        static_cast<int>(gvcf_frequency->integer_allele_counts[allele]),
                        0, sample_ploidy);
            } else if (!include_spanning_deletion) {
                const auto mle_ac = estimate_mle_allele_counts(
                    result, group.candidates, genotype_vectors, pl, sample_ploidy, true);
                for (std::size_t alt = 0; alt < mle_ac.size() && alt < output_mle_ac.size(); ++alt)
                    output_mle_ac[alt] = mle_ac[alt];
            }
            std::vector<int> ad{static_cast<int>(group.candidates.front()->reference_count)};
            for (const auto* candidate : group.candidates)
                ad.push_back(static_cast<int>(candidate->alternate_count));
            int depth = static_cast<int>(group.candidates.front()->depth);
            int variant_depth = static_cast<int>(group.candidates.front()->variant_depth);
            if (variant_depth == 0 && depth != 0) variant_depth = depth;
            const auto likelihood_depths_from_owner = [&]()
                -> std::optional<std::pair<std::vector<int>, int>> {
                if (const auto direct = derive_multiallelic_depths(
                        result, calls, include_spanning_deletion); direct.has_value())
                    return direct;
                for (const auto& owner : result.assembly_region_likelihood_results) {
                    if (owner == nullptr) continue;
                    if (const auto owned = derive_multiallelic_depths(
                            *owner, calls, include_spanning_deletion); owned.has_value())
                        return owned;
                }
                return std::nullopt;
            };
            if (const auto likelihood_depths = likelihood_depths_from_owner();
                likelihood_depths.has_value()) {
                ad = likelihood_depths->first;
                depth = likelihood_depths->second;
            }
            auto annotations = calls[best] != nullptr
                ? calls[best]->annotations : fastgatk::calling::GenotypeCall::Annotations{};
            if ((include_spanning_deletion || group.max_alt_subset) &&
                annotation_reads != nullptr) {
                const auto annotation_qual = calls[best] == nullptr ? 0.0 : calls[best]->qual;
                const auto annotations_from_owner = [&]()
                    -> std::optional<fastgatk::calling::GenotypeCall::Annotations> {
                    for (const auto& owner : result.assembly_region_likelihood_results) {
                        if (owner == nullptr) continue;
                        if (const auto owned = fastgatk::calling::calculate_output_variant_annotations(
                                *annotation_reads, *owner, group.candidates, annotation_qual,
                                informative_read_overlap_margin, include_spanning_deletion,
                                group.max_alt_subset); owned.has_value())
                            return owned;
                    }
                    return fastgatk::calling::calculate_output_variant_annotations(
                        *annotation_reads, result, group.candidates, annotation_qual,
                        informative_read_overlap_margin, include_spanning_deletion,
                        group.max_alt_subset);
                };
                if (const auto recalculated = annotations_from_owner(); recalculated.has_value())
                    annotations = *recalculated;
            }
            out << chrom << '\t' << position << "\t.\t" << group.reference << '\t';
            for (std::size_t i = 0; i < group.candidates.size(); ++i) {
                if (i != 0) out << ',';
                out << candidate_alternate(*group.candidates[i]);
            }
            if (include_spanning_deletion) out << ",*";
            const auto reference_end = position + static_cast<int>(group.candidates.front()->reference_allele.empty()
                ? 1 : group.candidates.front()->reference_allele.size()) - 1;
            const auto candidate_qual = calls[best] == nullptr ? std::string(".") :
                (gvcf_frequency.has_value() ? qual_text(gvcf_frequency->qual)
                                             : qual_text(calls[best]->qual));
            out << ",<NON_REF>\t" << candidate_qual << "\t.\t";
            bool info_started = false;
            const auto append_annotation = [&](const char* name, double value, int precision) {
                const auto text = annotation_text(value, precision);
                if (!text.empty()) {
                    if (info_started) out << ';';
                    out << name << '=' << text;
                    info_started = true;
                }
            };
            if (!bp_resolution && has_indel) {
                out << "END=" << reference_end;
                info_started = true;
            }
            append_annotation("BaseQRankSum", annotations.base_q_rank_sum, 3);
            if (info_started) out << ';';
            out << "DP=" << variant_depth << ";MLEAC=";
            for (std::size_t i = 0; i < output_mle_ac.size(); ++i) {
                if (i != 0) out << ',';
                out << output_mle_ac[i];
            }
            out << ";MLEAF=";
            for (std::size_t i = 0; i < output_mle_ac.size(); ++i) {
                if (i != 0) out << ',';
                const auto count = output_mle_ac[i];
                out << (count == 0 ? std::string("0.00") : format_allele_frequency(
                    static_cast<double>(count) / static_cast<double>(sample_ploidy)));
            }
            append_annotation("MQRankSum", annotations.mq_rank_sum, 3);
            if (annotations.raw_mq_depth != 0)
                out << ";RAW_MQandDP=" << annotations.raw_mq_sum_square
                    << ',' << annotations.raw_mq_depth;
            append_annotation("ReadPosRankSum", annotations.read_pos_rank_sum, 3);
            if (sites_only_vcf_output) {
                out << '\n';
                return;
            }
            out << "\tGT:AD:DP:GQ:PL:SB\t";
            for (std::size_t copy = 0; copy < genotype_vectors[best_genotype].size(); ++copy) {
                if (copy != 0) out << '/';
                out << genotype_vectors[best_genotype][copy];
            }
            out << ':';
            for (std::size_t i = 0; i < ad.size(); ++i) {
                if (i != 0) out << ',';
                out << ad[i];
            }
            out << ",0:" << depth << ':' << best_gq << ':';
            for (std::size_t i = 0; i < pl.size(); ++i) {
                if (i != 0) out << ',';
                out << pl[i];
            }
            out << ':' << annotations.ref_forward << ',' << annotations.ref_reverse
                << ',' << annotations.alt_forward << ',' << annotations.alt_reverse << '\n';
            return;
        }
        const auto genotype_index = [](std::size_t left, std::size_t right) {
            if (left > right) std::swap(left, right);
            return right * (right + 1) / 2 + left;
        };
        std::vector<int> ad{static_cast<int>(group.candidates.front()->reference_count)};
        int depth = static_cast<int>(group.candidates.front()->depth);
        int variant_depth = static_cast<int>(group.candidates.front()->variant_depth);
        if (variant_depth == 0 && depth != 0) variant_depth = depth;
        const int an = 2;
        std::vector<std::array<int, 3>> biallelic_pl(group.candidates.size(),
                                                       std::array<int, 3>{0, 999, 999});
        std::vector<int> full_pl;
        // A partitioned gVCF keeps concrete calls flattened for normal VCF
        // semantics, while each candidate's symbolic `<NON_REF>` and, when
        // applicable, EventMap `*` likelihoods remain owned by its
        // AssemblyRegion-local PairHMM matrix.  Probe that owning result when
        // the flattened matrix is intentionally absent; this is a Host
        // ownership lookup, not a numerical fallback.
        const auto joint_non_ref_from_owner = [&](bool include_spanning_deletion) {
            const auto allele_count = group.candidates.size() + 2U +
                (include_spanning_deletion ? 1U : 0U);
            const auto expected_width = genotype_width(static_cast<int>(allele_count), 2);
            if (joint_candidate_pl_with_nonref(result, group.candidates, 2, full_pl, 0,
                                               include_spanning_deletion) &&
                full_pl.size() == expected_width)
                return true;
            for (const auto& owner : result.assembly_region_likelihood_results) {
                if (owner != nullptr &&
                    joint_candidate_pl_with_nonref(*owner, group.candidates, 2, full_pl, 0,
                                                    include_spanning_deletion) &&
                    full_pl.size() == expected_width)
                    return true;
            }
            return false;
        };
        // As above, coverage establishes the original symbolic state, while
        // max-ALT selection establishes whether it remains in this record.
        const bool span_is_selected = covered_by_upstream_deletion &&
            (!group.max_alt_subset || group.max_alt_retains_spanning_deletion);
        const bool include_spanning_deletion = span_is_selected &&
            joint_non_ref_from_owner(true);
        const auto concrete_alleles = group.candidates.size() + 1U +
            (include_spanning_deletion ? 1U : 0U); // REF + concrete ALTs + `*`
        const auto allele_count = concrete_alleles + 1U; // plus `<NON_REF>`
        const auto concrete_end = concrete_alleles * (concrete_alleles + 1U) / 2U;
        std::vector<int> pl(allele_count * (allele_count + 1U) / 2U, 999);
        std::vector<int> concrete_pl;
        bool joint_non_ref = include_spanning_deletion || joint_non_ref_from_owner(false);
        if (group.max_alt_subset) {
            full_pl = group.max_alt_symbolic_pl;
            joint_non_ref = true;
        }
        bool joint = false;
        if (joint_non_ref) {
            pl = full_pl;
            concrete_pl.assign(full_pl.begin(), full_pl.begin() + concrete_end);
            joint = true;
        } else {
            joint = joint_candidate_pl(result, group.candidates, 2, concrete_pl, 0) &&
                    concrete_pl.size() == concrete_end;
            if (joint) std::copy(concrete_pl.begin(), concrete_pl.end(), pl.begin());
        }
        const auto likelihood_depths_from_owner = [&]()
            -> std::optional<std::pair<std::vector<int>, int>> {
            if (const auto direct = derive_multiallelic_depths(
                    result, calls, include_spanning_deletion); direct.has_value())
                return direct;
            for (const auto& owner : result.assembly_region_likelihood_results) {
                if (owner == nullptr) continue;
                if (const auto owned = derive_multiallelic_depths(
                        *owner, calls, include_spanning_deletion); owned.has_value())
                    return owned;
            }
            return std::nullopt;
        };
        bool multiallelic_ad = false;
        if (const auto likelihood_depths = likelihood_depths_from_owner();
            likelihood_depths.has_value() &&
            likelihood_depths->first.size() == concrete_alleles) {
            ad = likelihood_depths->first;
            depth = likelihood_depths->second;
            multiallelic_ad = true;
        }
        pl[genotype_index(0, 0)] = joint ? pl[genotype_index(0, 0)] : 999;
        for (std::size_t i = 0; i < group.candidates.size(); ++i) {
            if (!multiallelic_ad)
                ad.push_back(static_cast<int>(group.candidates[i]->alternate_count));
            if (!joint) {
                if (calls[i] == nullptr) continue;
                biallelic_pl[i] = phred_likelihoods(calls[i]->likelihoods);
                pl[genotype_index(0, 0)] = std::min(pl[genotype_index(0, 0)], biallelic_pl[i][0]);
                pl[genotype_index(0, i + 1)] = biallelic_pl[i][1];
                pl[genotype_index(i + 1, i + 1)] = biallelic_pl[i][2];
            }
        }
        if (!multiallelic_ad && include_spanning_deletion)
            ad.push_back(0);
        if (!joint) {
            for (std::size_t left = 0; left < group.candidates.size(); ++left)
                for (std::size_t right = left + 1; right < group.candidates.size(); ++right)
                    pl[genotype_index(left + 1, right + 1)] =
                        std::max(biallelic_pl[left][2], biallelic_pl[right][2]);
            const auto minimum_pl = *std::min_element(pl.begin(), pl.begin() + concrete_end);
            for (std::size_t index = 0; index < concrete_end; ++index)
                pl[index] = std::clamp(pl[index] - minimum_pl, 0, 999);
            // The independent biallelic fallback constructs the concrete
            // matrix directly in `pl`; retain that prefix for the shared
            // Kokkos GT/GQ derivation below just like the joint PairHMM path.
            concrete_pl.assign(pl.begin(), pl.begin() + concrete_end);
        }
        // If the per-read non-reference matrix is unavailable (for example,
        // the no-reference pileup fallback), retain a deterministic finite
        // envelope so the output remains structurally valid.  The
        // reference-backed PairHMM path above uses the GATK median-based
        // <NON_REF> likelihood instead.
        const auto non_ref = allele_count - 1U;
        if (joint && !joint_non_ref) {
            auto min_concrete = [&](const std::vector<std::size_t>& indices) {
                int value = 999;
                for (const auto index : indices)
                    if (index < concrete_pl.size()) value = std::min(value, concrete_pl[index]);
                return value;
            };
            std::vector<std::size_t> ref_non_ref;
            for (std::size_t alt = 1; alt < concrete_alleles; ++alt)
                ref_non_ref.push_back(genotype_index(0, alt));
            pl[genotype_index(0, non_ref)] = min_concrete(ref_non_ref);
            for (std::size_t alt = 1; alt < concrete_alleles; ++alt) {
                std::vector<std::size_t> alt_non_ref;
                for (std::size_t other = 1; other < concrete_alleles; ++other)
                    alt_non_ref.push_back(genotype_index(alt, other));
                pl[genotype_index(alt, non_ref)] = min_concrete(alt_non_ref);
            }
            std::vector<std::size_t> any_non_ref;
            for (std::size_t index = 1; index < concrete_alleles; ++index)
                for (std::size_t other = index; other < concrete_alleles; ++other)
                    any_non_ref.push_back(genotype_index(index, other));
            pl[genotype_index(non_ref, non_ref)] = min_concrete(any_non_ref);
        }
        const auto full_genotype_vectors = enumerate_genotypes(
            static_cast<int>(allele_count), 2);
        // Candidate-site QUAL in a gVCF is calculated over the complete
        // REF/concrete-ALT/<NON_REF> likelihood matrix.  Reusing the ordinary
        // biallelic call QUAL loses the symbolic-allele posterior (chr17:69067
        // is 60.31 in GATK rather than 60.32), even though GT and the concrete
        // PL prefix are unchanged.  The Host only supplies the GATK-shaped
        // allele prior vector; the AF posterior itself remains the shared
        // Kokkos kernel calculation.
        std::optional<fastgatk::kernels::AlleleFrequencyResult> gvcf_frequency;
        if (joint_non_ref && result.genotype_heterozygosity_stdev > 0.0 &&
            std::isfinite(result.genotype_heterozygosity_stdev)) {
            const auto ref_pseudocount = result.genotype_snp_heterozygosity /
                (result.genotype_heterozygosity_stdev *
                 result.genotype_heterozygosity_stdev);
            if (ref_pseudocount > 0.0 && std::isfinite(ref_pseudocount)) {
                std::vector<double> prior_pseudocounts(
                    allele_count, result.genotype_snp_heterozygosity * ref_pseudocount);
                prior_pseudocounts[0] = ref_pseudocount;
                for (std::size_t alt = 0; alt < group.candidates.size(); ++alt) {
                    const auto& candidate = *group.candidates[alt];
                    const bool indel = candidate_reference(candidate).size() !=
                                       candidate_alternate(candidate).size();
                    prior_pseudocounts[alt + 1] =
                        (indel ? result.genotype_indel_heterozygosity
                               : result.genotype_snp_heterozygosity) * ref_pseudocount;
                }
                if (include_spanning_deletion) {
                    const auto spanning = group.candidates.size() + 1U;
                    prior_pseudocounts[spanning] =
                        result.genotype_indel_heterozygosity * ref_pseudocount;
                }
                // See the arbitrary-ploidy writer: symbolic <NON_REF> has
                // zero HTSJDK allele length and therefore uses GATK's indel
                // prior in the AF calculation.
                prior_pseudocounts[non_ref] =
                    result.genotype_indel_heterozygosity * ref_pseudocount;
                try {
                    const auto calculated = fastgatk::kernels::calculate_allele_frequency_kokkos(
                        pl, 1, static_cast<int>(allele_count), 2, prior_pseudocounts,
                        {}, {}, include_spanning_deletion
                            ? static_cast<int>(group.candidates.size() + 1U) : -1);
                    if (calculated.samples_with_likelihoods != 0 &&
                        std::isfinite(calculated.qual))
                        gvcf_frequency = calculated;
                } catch (const std::exception&) {
                    // The structural gVCF writer already has a fail-closed
                    // Number=G fallback below.  Do not manufacture a QUAL
                    // when its full symbolic AF matrix is unusable.
                }
            }
        }
        std::vector<int> output_mle_ac(allele_count - 1U, 0);
        if (gvcf_frequency.has_value() &&
            gvcf_frequency->integer_allele_counts.size() == allele_count) {
            for (std::size_t allele = 1; allele < allele_count; ++allele)
                output_mle_ac[allele - 1U] = std::clamp(
                    static_cast<int>(gvcf_frequency->integer_allele_counts[allele]), 0, an);
        } else if (!include_spanning_deletion) {
            const auto mle_ac = estimate_mle_allele_counts(
                result, group.candidates, full_genotype_vectors, pl, 2, true);
            for (std::size_t alt = 0; alt < mle_ac.size() && alt < output_mle_ac.size(); ++alt)
                output_mle_ac[alt] = mle_ac[alt];
        }
        // Candidate-site GT/GQ must use the same Kokkos primitive as ordinary
        // VCF output.  `<NON_REF>` is an evidence envelope, not a concrete
        // called allele, so assign from the concrete PL prefix and keep the
        // symbolic rows only in the emitted Number=G matrix.
        const auto concrete_genotype_vectors = enumerate_genotypes(
            static_cast<int>(concrete_alleles), 2);
        if (concrete_pl.size() < concrete_genotype_vectors.size())
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: diploid concrete PL prefix is incomplete");
        const std::vector<int> concrete_call_pl(
            concrete_pl.begin(), concrete_pl.begin() + concrete_genotype_vectors.size());
        const auto genotype_priors = genotype_priors_for_group(
            result, group.candidates, 2, false, include_spanning_deletion);
        if (genotype_priors.size() != concrete_genotype_vectors.size())
            throw std::runtime_error("NUMERICAL_CONTRACT_FAILURE: gVCF genotype prior width mismatch");
        const auto derived = result.genotype_priors_used
            ? fastgatk::kernels::derive_genotype_gt_gq_from_log10_priors_kokkos(
                concrete_call_pl, 1, static_cast<int>(concrete_alleles), 2, genotype_priors)
            : fastgatk::kernels::derive_genotype_gt_gq_kokkos(
                concrete_call_pl, 1, static_cast<int>(concrete_alleles), 2);
        if (derived.alleles.size() != 2)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Kokkos gVCF GT is incomplete");
        const std::size_t best_left = static_cast<std::size_t>(std::max(0, derived.alleles[0]));
        const std::size_t best_right = static_cast<std::size_t>(std::max(0, derived.alleles[1]));
        const int best_gq = derived.gq.empty() ? 0 : std::clamp(derived.gq.front(), 0, 99);
        auto annotations = calls[best] != nullptr
            ? calls[best]->annotations : fastgatk::calling::GenotypeCall::Annotations{};
        if ((include_spanning_deletion || group.max_alt_subset) &&
            annotation_reads != nullptr) {
            const auto annotation_qual = calls[best] == nullptr ? 0.0 : calls[best]->qual;
            const auto annotations_from_owner = [&]()
                -> std::optional<fastgatk::calling::GenotypeCall::Annotations> {
                // Owners that hold the identical-allele twin with a real
                // PairHMM context come first; owners that merely carry the
                // allele with the UINT32_MAX context sentinel keep their
                // historical first-match order as a fallback.  Selecting the
                // twin owner is what makes RAW_MQandDP/SB independent of the
                // `-L` window start.
                std::vector<const fastgatk::calling::Result*> ordered_owners;
                ordered_owners.reserve(result.assembly_region_likelihood_results.size());
                for (const auto& owner : result.assembly_region_likelihood_results)
                    if (owner != nullptr && owner_has_pairhmm_context(*owner, group.candidates))
                        ordered_owners.push_back(owner.get());
                for (const auto& owner : result.assembly_region_likelihood_results)
                    if (owner != nullptr && !owner_has_pairhmm_context(*owner, group.candidates))
                        ordered_owners.push_back(owner.get());
                for (const auto* owner : ordered_owners) {
                    if (const auto owned = fastgatk::calling::calculate_output_variant_annotations(
                            *annotation_reads, *owner, group.candidates, annotation_qual,
                            informative_read_overlap_margin, include_spanning_deletion,
                            group.max_alt_subset); owned.has_value())
                        return owned;
                }
                return fastgatk::calling::calculate_output_variant_annotations(
                    *annotation_reads, result, group.candidates, annotation_qual,
                    informative_read_overlap_margin, include_spanning_deletion,
                    group.max_alt_subset);
            };
            if (const auto recalculated = annotations_from_owner(); recalculated.has_value())
                annotations = *recalculated;
        }
        out << chrom << '\t' << position << "\t.\t" << group.reference << '\t';
        for (std::size_t i = 0; i < group.candidates.size(); ++i) {
            if (i != 0) out << ',';
            out << candidate_alternate(*group.candidates[i]);
        }
        if (include_spanning_deletion) out << ",*";
        const auto reference_end = position + static_cast<int>(group.candidates.front()->reference_allele.empty()
            ? 1 : group.candidates.front()->reference_allele.size()) - 1;
        // GVCF candidate sites are ordinary VariantContext records in GATK:
        // retain the site QUAL, use an unfiltered '.', and serialize only the
        // standard candidate INFO/FORMAT fields.  The native-only AC/AF/MQ/
        // QD/SOR/GQ INFO envelope made these rows difficult for downstream
        // CombineGVCFs to distinguish from a normal VCF call.
        // GenotypingEngine retains the AF log10PError on a concrete gVCF
        // record even when the selected concrete genotype is 0/0.  A forced
        // EventMap allele at chr20:10020493, for example, has QUAL 0.01;
        // collapsing every hom-ref concrete row to zero loses that posterior.
        // QUAL is printed to two decimals for concrete candidates. Preserve
        // GATK's literal `0` for a hom-ref EventMap posterior that rounds to
        // zero; only a representable nonzero full-posterior value (for
        // example 0.012 -> 0.01) survives as a gVCF candidate QUAL.
        constexpr double kGvcfPrintedNonzeroQual = 0.005;
        const auto candidate_qual = calls[best] == nullptr ? std::string(".") :
            (gvcf_frequency.has_value() &&
             gvcf_frequency->qual >= kGvcfPrintedNonzeroQual
                ? qual_text(gvcf_frequency->qual)
                : (calls[best]->genotype == 0
                    ? std::string("0") : qual_text(calls[best]->qual)));
        out << ",<NON_REF>\t" << candidate_qual << "\t.\t";
        bool info_started = false;
        if (!bp_resolution && has_indel) {
            out << "END=" << reference_end;
            info_started = true;
        }
        const auto append_annotation = [&](const char* name, double value, int precision) {
            const auto text = annotation_text(value, precision);
            if (!text.empty()) {
                if (info_started) out << ';';
                out << name << '=' << text;
                info_started = true;
            }
        };
        append_annotation("BaseQRankSum", annotations.base_q_rank_sum, 3);
        if (info_started) out << ';';
        out << "DP=" << variant_depth << ";ExcessHet=0.0000;MLEAC=";
        info_started = true;
        for (std::size_t i = 0; i < output_mle_ac.size(); ++i) {
            if (i != 0) out << ',';
            out << output_mle_ac[i];
        }
        out << ";MLEAF=";
        for (std::size_t i = 0; i < output_mle_ac.size(); ++i) {
            if (i != 0) out << ',';
            const auto count = output_mle_ac[i];
            // HTSJDK keeps the MLEAF zero representation at 0.00 but uses
            // three fractional places for non-zero fractional MLEAF values.
            // Keep this gVCF-only text rule local; ordinary VCF keeps its
            // established Number=A serialization contract.
            out << (count == 0 ? std::string("0.00") : format_allele_frequency(
                static_cast<double>(count) / static_cast<double>(an)));
        }
        append_annotation("MQRankSum", annotations.mq_rank_sum, 3);
        if (annotations.raw_mq_depth != 0)
            out << ";RAW_MQandDP=" << annotations.raw_mq_sum_square
                << ',' << annotations.raw_mq_depth;
        append_annotation("ReadPosRankSum", annotations.read_pos_rank_sum, 3);
        if (sites_only_vcf_output) {
            out << '\n';
            return;
        }
        const bool reverse_het_gt = physical_phase.enabled && !physical_phase.phase_01 &&
            best_left != best_right && !include_spanning_deletion;
        out << "\tGT:AD:DP:GQ";
        if (physical_phase.enabled)
            out << ":PGT:PID:PL:PS:SB\t";
        else
            out << ":PL:SB\t";
        out << (reverse_het_gt ? best_right : best_left)
            << (physical_phase.enabled ? '|' : '/')
            << (reverse_het_gt ? best_left : best_right) << ':';
        for (std::size_t i = 0; i < ad.size(); ++i) {
            if (i != 0) out << ',';
            out << ad[i];
        }
        out << ",0:" << depth << ':' << best_gq << ':';
        if (physical_phase.enabled)
            out << (physical_phase.phase_01 ? "0|1" : "1|0") << ':'
                << physical_phase.id << ':';
        for (std::size_t i = 0; i < pl.size(); ++i) {
            if (i != 0) out << ',';
            out << pl[i];
        }
        if (physical_phase.enabled) out << ':' << physical_phase.phase_set;
        out << ':' << annotations.ref_forward << ',' << annotations.ref_reverse
            << ',' << annotations.alt_forward << ',' << annotations.alt_reverse << '\n';
    };
    std::size_t block_index = 0;
    std::size_t candidate_index = 0;
    while (block_index < result.reference_blocks.size() ||
           candidate_index < candidate_groups.size()) {
        const bool candidate_first = candidate_index < candidate_groups.size() &&
            (block_index == result.reference_blocks.size() ||
             candidate_groups[candidate_index].tid < result.reference_blocks[block_index].tid ||
             (candidate_groups[candidate_index].tid == result.reference_blocks[block_index].tid &&
              candidate_groups[candidate_index].position <= result.reference_blocks[block_index].start));
        if (candidate_first) {
            write_candidate(candidate_groups[candidate_index], physical_phases[candidate_index],
                            gvcf_group_covered_by_upstream_deletion[candidate_index] != 0U);
            ++candidate_index;
        }
        else write_block(result.reference_blocks[block_index++]);
    }
    return out.str();
}

bool file_complete(const std::string& path) {
    if (path == "-") return true;
    std::error_code error;
    const auto regular = std::filesystem::is_regular_file(path, error);
    if (error || !regular) return false;
    const auto size = std::filesystem::file_size(path, error);
    return !error && size > 0;
}

bool output_contract_valid(const Options& options) {
    if (!file_complete(options.output)) return false;
    if (!options.index_path.empty() && !file_complete(options.index_path)) return false;
    return options.assembly_region_out.empty() || file_complete(options.assembly_region_out);
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

void validate_published_manifest(const Options& options, const std::string& manifest_path) {
    // Stdout is intentionally a stream-only mode and has no publishable
    // output bundle. Every file-boundary path with a manifest goes through
    // the shared runtime validator before success is reported.
    if (manifest_path.empty() || options.output == "-") return;
    fastgatk::runtime::OutputBundle bundle{
        std::filesystem::path(options.output),
        std::filesystem::path(options.index_path),
        std::filesystem::path(manifest_path),
        !options.index_path.empty(),
        true};
    fastgatk::runtime::require_complete_output(bundle);
}

void validate_input_headers(const fastgatk::io::HtsReader& first,
                            const fastgatk::io::HtsReader& candidate,
                            const bool explicit_sample) {
    const auto& expected = first.header();
    const auto& actual = candidate.header();
    if (expected.contigs != actual.contigs || expected.contig_lengths != actual.contig_lengths)
        throw std::runtime_error(
            "BAD_INPUT: all HaplotypeCaller -I inputs must share the same sequence dictionary");
    // HaplotypeCaller combines shards for one sample. An explicit
    // --sample-name lets each input contain additional samples, because the
    // Host reader has already filtered them. Without it, silently combining
    // different sample sets would produce a misleading single-sample VCF.
    if (!explicit_sample) {
        auto expected_samples = expected.samples;
        auto actual_samples = actual.samples;
        std::sort(expected_samples.begin(), expected_samples.end());
        std::sort(actual_samples.begin(), actual_samples.end());
        if (expected_samples != actual_samples)
            throw std::runtime_error(
                "BAD_INPUT: all HaplotypeCaller -I inputs must expose the same sample set; "
                "use --sample-name to select one sample explicitly");
    }
}

// GATK validates the read sequence dictionary against the reference before
// entering the walker, but accepts a reference dictionary that is a proper
// subset of a BAM dictionary (for example its bundled chrM-only fixture).
// Keep that name/length compatibility check at the Host boundary so every
// backend sees the same contig mapping. The explicit disable escape hatch is
// still not an implicit fallback for a mismatched or malformed reference.
void validate_reference_dictionary(const std::string& reference,
                                   const fastgatk::io::HeaderSummary& observed) {
#if FASTGATK_HAS_HTSLIB
    // The native smoke path deliberately accepts a plain, unindexed FASTA
    // (its sequence loader is a small Host parser).  Check for the sidecar
    // before calling fai_load: HTSlib may create a missing FAI on demand, and
    // that must not turn an unindexed smoke input into an implicit dictionary.
    std::error_code error;
    const auto fai_path = std::filesystem::path(reference).string() + ".fai";
    const bool has_fai = std::filesystem::exists(fai_path, error);
    if (error)
        throw std::runtime_error("BACKEND_UNAVAILABLE: cannot inspect reference FAI: " + fai_path);
    if (!has_fai) return;
    faidx_t* fai = fai_load(reference.c_str());
    if (!fai)
        throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA has an unreadable .fai: " + reference);
    try {
        std::map<std::string, std::int64_t> input_contigs;
        for (std::size_t index = 0; index < observed.contigs.size(); ++index)
            input_contigs.emplace(observed.contigs[index], observed.contig_lengths[index]);
        const auto expected_count = faidx_nseq(fai);
        for (int index = 0; index < expected_count; ++index) {
            const char* name = faidx_iseq(fai, index);
            const auto length = name == nullptr ? -1 : static_cast<std::int64_t>(faidx_seq_len64(fai, name));
            const auto input = name == nullptr ? input_contigs.end() : input_contigs.find(name);
            if (input == input_contigs.end() || input->second != length) {
                throw std::runtime_error(
                    "BAD_INPUT: reference sequence dictionary is incompatible with input BAM (reference=" +
                    (name == nullptr ? std::string("<missing>") : std::string(name)) + ", input=" +
                    (input == input_contigs.end() ? std::string("<missing>") : input->first) + ")");
            }
        }
        fai_destroy(fai);
    } catch (...) {
        fai_destroy(fai);
        throw;
    }
#else
    (void)reference;
    (void)observed;
    throw std::runtime_error(
        "BACKEND_UNAVAILABLE: sequence dictionary validation requires HTSlib FAIDX");
#endif
}

std::string make_manifest(const Options& options, const fastgatk::io::HtsReader& reader,
                          const fastgatk::calling::Result& result,
                          const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto telemetry = summary_json(options, reader, result, resources);
    const bool interval_subset = !options.calling.intervals.empty() ||
        (options.calling.interval_tid >= 0 && options.calling.interval_end > options.calling.interval_start);
    const bool complete = output_contract_valid(options);
    bool non_ref_likelihoods = result.non_ref_likelihoods_present;
    const auto assembly_tags = static_cast<std::size_t>(std::count_if(
        reader.header().contig_assemblies.begin(), reader.header().contig_assemblies.end(),
        [](const std::string& value) { return !value.empty(); }));
    for (const auto& candidate : result.non_ref_read_likelihoods)
        for (const auto likelihood : candidate)
            if (std::isfinite(likelihood)) { non_ref_likelihoods = true; break; }
        
    std::ostringstream out;
    out << "{\"schema_version\":1,\"tool\":\"HaplotypeCaller\","
        << "\"implementation\":\"fastgatk-hc-call\",\"status\":\"prototype\","
        << "\"reads\":" << result.reads << ","
        << "\"input_count\":" << options.inputs.size() << ","
        << "\"inputs\":[";
    for (std::size_t input_index = 0; input_index < options.inputs.size(); ++input_index) {
        if (input_index != 0) out << ',';
        out << "\"" << json_escape(options.inputs[input_index]) << "\"";
    }
    out << "],"
        << "\"interval_set_rule\":\""
        << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
            ? "INTERSECTION" : "UNION") << "\","
        << "\"sample_name\":\"" << json_escape(
            options.sample_name.empty() ? (reader.header().samples.empty()
                ? std::string("FASTGATK") : reader.header().samples.front())
                : options.sample_name) << "\","
        << "\"primary_output\":\"" << json_escape(options.output) << "\","
        << "\"primary_output_kind\":\"" << (options.gvcf ? "gvcf" : (is_vcf_output(options.output) ? "vcf" : "json")) << "\","
        << "\"compatibility\":{\"gatk_parameter_aliases\":true,"
        << "\"create_output_variant_index\":"
        << (options.create_output_variant_index ? "true" : "false") << ","
        << "\"add_output_vcf_command_line\":"
        << (options.add_output_vcf_command_line ? "true" : "false") << ","
        << "\"sites_only_vcf_output\":"
        << (options.sites_only_vcf_output ? "true" : "false") << ","
        << "\"floor_blocks\":" << (options.floor_blocks ? "true" : "false") << ","
        << "\"alleles_feature\":" << (!options.alleles.empty() ? "true" : "false") << ","
        << "\"forced_alleles\":" << options.calling.forced_alleles.size() << ","
        << "\"force_call_filtered_alleles\":"
        << (options.force_call_filtered_alleles ? "true" : "false") << ","
        << "\"native_pair_hmm_use_double_precision\":"
        << (options.native_pair_hmm_use_double_precision ? "true" : "false") << ","
        << "\"pairhmm_effective_precision\":\""
        << json_escape(result.pairhmm_precision) << "\","
        << "\"informative_read_overlap_margin\":"
        << options.calling.informative_read_overlap_margin << ","
        << "\"disable_sequence_dictionary_validation\":"
        << (options.disable_sequence_dictionary_validation ? "true" : "false") << ","
        << "\"interval_set_rule\":\""
        << (options.interval_set_rule == fastgatk::io::HtsIntervalSetRule::Intersection
            ? "INTERSECTION" : "UNION") << "\","
        << "\"sample_name_selection\":" << (!options.sample_name.empty() ? "true" : "false") << ","
        << "\"cigar_aware_projection\":true,"
        << "\"use_soft_clipped_bases\":"
        << (options.calling.use_soft_clipped_bases ? "true" : "false") << ","
        << "\"soft_clip_low_quality_ends\":"
        << (options.calling.soft_clip_low_quality_ends ? "true" : "false") << ","
        << "\"use_haplotype_realignment_for_rcm\":"
        << (options.calling.use_haplotype_realignment_for_rcm ? "true" : "false") << ","
        << "\"pcr_indel_model\":\"" << json_escape(options.calling.pairhmm_conservative_indel_model
            ? options.calling.pairhmm_pcr_indel_model : "NONE") << "\","
        << "\"pcr_error_rate_factor\":" << (options.calling.pairhmm_conservative_indel_model
            ? options.calling.pairhmm_pcr_error_rate_factor : 0.0) << ","
        << "\"phred_scaled_global_read_mismapping_rate\":"
        << options.calling.phred_scaled_global_read_mismapping_rate << ","
        << "\"contig_assembly_metadata\":" << (assembly_tags > 0 ? "true" : "false") << ","
        << "\"assembly_region_out\":"
        << (!options.assembly_region_out.empty() ? "true" : "false") << ","
        << "\"force_active\":" << (options.calling.force_active ? "true" : "false") << ","
        << "\"stream_by_contig\":" << (options.stream_by_contig ? "true" : "false") << ","
        << "\"stream_by_region\":" << (options.stream_by_region ? "true" : "false") << ","
        << "\"stream_region_size\":" << options.stream_region_size << ","
        << "\"streamed_contigs\":" << options.streamed_contigs << ","
        << "\"streamed_regions\":" << options.streamed_regions << ","
        << "\"streamed_region_splits\":" << options.streamed_region_splits << ","
        << "\"streamed_peak_host_bytes\":" << options.streamed_peak_host_bytes << ","
        << "\"stream_indexed\":" << (options.stream_indexed ? "true" : "false") << ","
        << "\"pipeline_lifecycle\":\""
        << (options.pipeline_used
            ? "Host decode->bounded queue->Kokkos compute->encode->sink" : "not-used") << "\","
        << "\"pipeline_decoded_items\":" << options.pipeline_decoded_items << ","
        << "\"pipeline_computed_items\":" << options.pipeline_computed_items << ","
        << "\"pipeline_encoded_items\":" << options.pipeline_encoded_items << ","
        << "\"pipeline_decoded_bytes\":" << options.pipeline_decoded_bytes << ","
        << "\"pipeline_computed_bytes\":" << options.pipeline_computed_bytes << ","
        << "\"pipeline_encoded_bytes\":" << options.pipeline_encoded_bytes << ","
        << "\"pipeline_peak_decoded_bytes\":" << options.pipeline_peak_decoded_bytes << ","
        << "\"pipeline_peak_computed_bytes\":" << options.pipeline_peak_computed_bytes << ","
        << "\"pipeline_peak_encoded_bytes\":" << options.pipeline_peak_encoded_bytes << ","
        << "\"pairhmm_haplotype_likelihoods\":" << (result.pairhmm_used ? "true" : "false") << ","
        << "\"pairhmm_likelihood_normalization\":"
        << (!result.pairhmm_normalization_execution_space.empty() ? "true" : "false") << ","
        << "\"pairhmm_allele_marginalization\":"
        << (!result.pairhmm_marginalization_execution_space.empty() ? "true" : "false") << ","
        << "\"pairhmm_read_allele_uncertainty\":"
        << (!result.pairhmm_uncertainty_execution_space.empty() ? "true" : "false") << ","
        << "\"pairhmm_error_model\":\"" << json_escape(result.pairhmm_error_model) << "\"," 
        << "\"pairhmm_default_indel_quality\":" << static_cast<unsigned>(result.pairhmm_default_indel_quality) << ","
        << "\"pairhmm_default_indel_terminal_quality\":" << static_cast<unsigned>(result.pairhmm_default_indel_terminal_quality) << ","
        << "\"pairhmm_pcr_indel_model\":\"" << json_escape(result.pairhmm_pcr_indel_model) << "\","
        << "\"pairhmm_pcr_error_rate_factor\":" << result.pairhmm_pcr_error_rate_factor << ","
        << "\"pairhmm_pcr_adjusted_positions\":" << result.pairhmm_pcr_adjusted_positions << ","
        << "\"pairhmm_flow_reads\":" << result.pairhmm_flow_reads << ","
        << "\"pairhmm_flow_haplotypes\":" << result.pairhmm_flow_haplotypes << ","
        << "\"pairhmm_flow_haplotypes_clipped\":" << result.pairhmm_flow_haplotypes_clipped << ","
        << "\"pairhmm_flow_haplotypes_collapsed\":" << result.pairhmm_flow_haplotypes_collapsed << ","
        << "\"pairhmm_flow_haplotypes_uncollapsed\":" << result.pairhmm_flow_haplotypes_uncollapsed << ","
        << "\"pairhmm_flow_haplotype_remaps\":" << result.pairhmm_flow_haplotype_remaps << ","
        << "\"pairhmm_flow_identical_haplotype_groups\":" << result.pairhmm_flow_identical_haplotype_groups << ","
        << "\"pairhmm_flow_reads_clipped\":" << result.pairhmm_flow_reads_clipped << ","
        << "\"pairhmm_flow_clipping_fallbacks\":" << result.pairhmm_flow_clipping_fallbacks << ","
        << "\"pairhmm_insertion_quality_reads\":" << result.pairhmm_insertion_quality_reads << ","
        << "\"pairhmm_deletion_quality_reads\":" << result.pairhmm_deletion_quality_reads << ","
        << "\"pairhmm_reads_clipped\":" << result.pairhmm_reads_clipped << ","
        << "\"pairhmm_reads_dropped_after_clipping\":"
        << result.pairhmm_reads_dropped_after_clipping << ","
        << "\"pairhmm_reads_disqualified\":"
        << result.pairhmm_reads_disqualified << ","
        << "\"pairhmm_read_disqualification_threshold\":"
        << result.pairhmm_read_disqualification_threshold << ","
        << "\"pairhmm_fallback_reason\":\"" << json_escape(result.pairhmm_skip_reason) << "\","
        << "\"pairhmm_haplotype_set_marginalization\":"
        << (result.pairhmm_marginalized_candidates > 0 ? "true" : "false") << ","
        << "\"likelihood_based_haplotype_pruning\":"
        << (result.pairhmm_graph_haplotypes_considered > 0 ? "true" : "false") << ","
        << "\"graph_haplotype_likelihoods\":"
        << (result.pairhmm_graph_haplotypes > 0 ? "true" : "false") << ","
        << "\"graph_snp_posterior_ownership\":"
        << (result.pairhmm_graph_snp_posterior_pairs > 0 ? "true" : "false") << ","
        << "\"haplotype_combination_chunking\":"
        << (result.pairhmm_haplotype_combination_chunked ? "true" : "false") << ","
        << "\"assembly_region_pairhmm_partitioning\":"
        << (result.pairhmm_assembly_region_partitioned ? "true" : "false") << ","
        << "\"assembly_region_candidate_partitioning\":"
        << (result.assembly_region_partitioned ? "true" : "false") << ","
        << "\"multiallelic_locus_grouping\":true,\"shared_observation_ranges\":true,"
        << "\"arbitrary_ploidy_vcf\":" << (options.sample_ploidy != 2 ? "true" : "false") << ","
        << "\"site_annotations\":true,"
        << "\"smith_waterman_scores\":" << (result.sw_used ? "true" : "false") << ","
        << "\"smith_waterman_simd\":" << (result.sw_simd_groups > 0 ? "true" : "false") << ","
        << "\"smith_waterman_traceback\":" << (result.sw_traceback_pairs > 0 ? "true" : "false") << ","
        << "\"read_haplotype_cigar_provenance\":"
        << (result.read_haplotype_cigar_valid_pairs > 0 ? "true" : "false") << ","
        << "\"read_haplotype_uncertainty_margin_telemetry\":"
        << (result.read_haplotype_informative_reads > 0 ? "true" : "false") << ","
        << "\"kmer_graph\":" << (result.graph_used ? "true" : "false") << ","
        << "\"reference_connected_graph\":" << (result.graph_reference_nodes > 0 ? "true" : "false") << ","
        << "\"haplotype_traversal\":" << (result.graph_haplotype_paths > 0 ? "true" : "false") << ","
        << "\"haplotype_sequence_materialization\":"
        << ((!result.graph_haplotype_sequences.empty() || result.graph_haplotype_sequence_count > 0)
            ? "true" : "false") << ","
        << "\"terminal_softclip_candidates\":" << (result.softclip_candidates > 0 ? "true" : "false") << ","
        << "\"compound_cigar_candidates\":" << (result.compound_cigar_candidates > 0 ? "true" : "false") << ","
        << "\"graph_variant_candidate_extraction\":" << (result.graph_variant_candidates > 0 ? "true" : "false") << ","
        << "\"graph_snp_candidate_extraction\":" << (result.graph_snp_candidates > 0 ? "true" : "false") << ","
        << "\"graph_mnp_candidate_extraction\":" << (result.graph_mnp_candidates > 0 ? "true" : "false") << ","
        << "\"graph_deletion_candidate_extraction\":" << (result.graph_deletion_candidates > 0 ? "true" : "false") << ","
        << "\"activity_profile\":" << (result.activity_used ? "true" : "false") << ","
        << "\"activity_reference_projection\":"
        << (result.activity_reference_aware ? "true" : "false") << ","
        << "\"activity_profile_controls\":true,"
        << "\"activity_bandpass\":true,"
        << "\"min_assembly_region_size\":"
        << options.calling.min_assembly_region_size << ","
        << "\"activity_filter_size\":" << result.activity_filter_size << ","
        << "\"activity_effective_max_probability_propagation_distance\":"
        << result.activity_effective_max_probability_propagation_distance << ","
        << "\"activity_probability_threshold\":"
        << options.calling.active_probability_threshold << ","
        << "\"assembly_region_padding\":"
        << options.calling.assembly_region_padding << ","
        << "\"max_assembly_region_size\":"
        << options.calling.max_assembly_region_size << ","
        << "\"max_probability_propagation_distance\":"
        << options.calling.max_probability_propagation_distance << ","
        << "\"reference_confidence_kokkos\":"
        << (!result.reference_confidence_execution_space.empty() ? "true" : "false") << ","
        << "\"rcm_haplotype_realignment\":"
        << (result.rcm_haplotype_realignment_used ? "true" : "false") << ","
        << "\"read_filter_order\":" << (result.read_filter_used ? "true" : "false") << ","
        << "\"read_filter_applied\":" << (result.read_filter_applied ? "true" : "false") << ","
        << "\"overlapping_base_quality_correction\":"
        << (result.overlapping_quality_correction_used ? "true" : "false") << ","
        << "\"overlapping_quality_correction_metadata_available\":"
        << (result.overlapping_quality_correction_metadata_available ? "true" : "false") << ","
        << "\"overlapping_pairs\":" << result.overlapping_pairs << ","
        << "\"overlapping_bases\":" << result.overlapping_bases << ","
        << "\"overlapping_conflicting_bases\":" << result.overlapping_conflicting_bases << ","
        << "\"overlapping_quality_caps\":" << result.overlapping_quality_caps << ","
        << "\"do_not_correct_overlapping_base_qualities\":"
        << (options.calling.do_not_correct_overlapping_base_qualities ? "true" : "false") << ","
        << "\"read_filter_gatk_defaults\":"
        << (options.calling.minimum_mapping_quality == 20 && options.calling.exclude_duplicates &&
            options.calling.exclude_unmapped && options.calling.exclude_secondary &&
            !options.calling.exclude_supplementary && options.calling.exclude_qcfail &&
            options.calling.exclude_mapping_quality_unavailable &&
            options.calling.require_good_cigar && options.calling.require_nonzero_reference_span &&
            options.calling.require_no_n_cigar && options.calling.require_read_group
                ? "true" : "false") << ","
        << "\"disable_tool_default_read_filters\":"
        << (options.calling.disable_tool_default_read_filters ? "true" : "false") << ","
        << "\"min_base_quality\":"
        << static_cast<unsigned>(options.calling.min_base_quality) << ","
        << "\"reference_model_deletion_quality\":"
        << static_cast<unsigned>(options.calling.reference_model_deletion_quality) << ","
        << "\"read_filter_require_read_length\":"
        << (options.calling.require_read_length ? "true" : "false") << ","
        << "\"read_filter_min_read_length\":" << options.calling.min_read_length << ","
        << "\"read_filter_max_read_length\":" << options.calling.max_read_length << ","
        << "\"deterministic_downsampling\":"
        << (options.calling.max_reads_per_locus > 0 ? "true" : "false") << ","
        << "\"vcf_index\":" << (!options.index_path.empty() ? "true" : "false") << ","
        << "\"gvcf\":" << (options.gvcf ? "true" : "false") << ","
        << "\"sample_ploidy\":" << options.sample_ploidy << ","
        << "\"gvcf_gq_bands\":[";
    for (std::size_t index = 0; index < options.calling.gvcf_gq_bands.size(); ++index) {
        if (index != 0) out << ',';
        out << options.calling.gvcf_gq_bands[index];
    }
    out << "],"
        << "\"gvcf_semantics\":\""
        << (options.gvcf ? (options.calling.emit_reference_confidence_bp_resolution
            ? "base-pair-resolution+candidate-sites" : "reference-blocks+candidate-sites")
                         : "not-requested") << "\"," 
        << "\"gvcf_standard_fields\":"
        << (options.gvcf && !options.calling.emit_reference_confidence_bp_resolution ? "true" : "false") << ","
        << "\"gvcf_floor_blocks\":" << (options.floor_blocks ? "true" : "false") << ","
        << "\"gvcf_candidate_standard_fields\":"
        << (options.gvcf ? "true" : "false") << ","
        << "\"gvcf_candidate_rebuild_after_genotyping\":"
        << (options.gvcf ? "true" : "false") << ","
        << "\"gvcf_reference_confidence_likelihoods\":" << (options.gvcf ? "true" : "false") << ","
        << "\"reference_block_lookup_indexed\":"
        << (result.reference_block_lookup_indexed ? "true" : "false") << ","
        << "\"gvcf_gq_band_state_machine\":" << (options.gvcf ? "true" : "false") << ","
        << "\"gvcf_non_ref_likelihood_model\":\""
        << (non_ref_likelihoods ? "median-qualified-concrete-alleles" : "concrete-envelope")
        << "\"," 
        << "\"gvcf_interval_reference_blocks\":"
        << (interval_subset ? "true" : "false") << ","
        << "\"interval_exclusions\":" << options.calling.exclusion_intervals.size() << ","
        << "\"interval_padding\":" << options.interval_padding << ","
        << "\"interval_exclusion_padding\":" << options.exclusion_padding << ","
        << "\"interval_list_inputs\":" << reader.interval_file_inputs() << ","
        << "\"interval_list_records\":" << reader.interval_file_records() << ","
        << "\"genotype_priors\":" << (options.calling.use_genotype_priors ? "true" : "false") << ","
        << "\"genotype_assignment_method\":\""
        << (options.calling.use_posterior_genotype_assignment
            ? "USE_POSTERIOR_PROBABILITIES" : "USE_PLS_TO_ASSIGN") << "\","
        << "\"joint_genotype_priors_used\":"
        << (result.genotype_priors_used ? "true" : "false") << ","
        << "\"genotype_prior_model\":\"GATK-assumingHW-v1\"," 
        << "\"genotype_prior_kernel_calls\":" << result.genotype_prior_kernel_calls << ","
        << "\"genotype_prior_prepare_seconds\":" << result.genotype_prior_prepare_seconds << ","
        << "\"genotype_prior_seconds\":" << result.genotype_prior_seconds << ","
        << "\"genotype_prior_execution_space\":\""
        << json_escape(result.genotype_prior_execution_space) << "\","
        << "\"heterozygosity\":" << options.calling.heterozygosity << ","
        << "\"indel_heterozygosity\":" << options.calling.indel_heterozygosity << ","
        << "\"indel_size_to_eliminate_in_ref_model\":"
        << options.calling.indel_size_to_eliminate_in_ref_model << ","
        << "\"graph_kmer_size\":" << options.calling.graph_kmer_size
        << ",\"graph_kmer_sizes\":[";
    for (std::size_t index = 0; index < options.calling.graph_kmer_sizes.size(); ++index) {
        if (index != 0) out << ',';
        out << options.calling.graph_kmer_sizes[index];
    }
    out << "],"
        << "\"graph_min_kmer_count\":" << options.calling.graph_min_kmer_count << ","
        << "\"graph_min_pruning\":" << options.calling.graph_min_pruning << ","
        << "\"graph_num_pruning_samples\":" << options.calling.graph_num_pruning_samples << ","
        << "\"graph_use_adaptive_pruning\":"
        << (options.calling.graph_use_adaptive_pruning ? "true" : "false") << ","
        << "\"graph_initial_error_rate_for_pruning\":"
        << options.calling.graph_initial_error_rate_for_pruning << ","
        << "\"graph_pruning_log_odds_threshold\":"
        << options.calling.graph_pruning_log_odds_threshold << ","
        << "\"graph_pruning_seeding_log_odds_threshold\":"
        << options.calling.graph_pruning_seeding_log_odds_threshold << ","
        << "\"graph_max_unpruned_variants\":"
        << options.calling.graph_max_unpruned_variants << ","
        << "\"graph_adaptive_pruned_nodes\":" << result.graph_adaptive_pruned_nodes << ","
        << "\"graph_linked_de_bruijn\":"
        << (options.calling.graph_linked_de_bruijn ? "true" : "false") << ","
        << "\"graph_disable_artificial_haplotype_recovery\":"
        << (options.calling.graph_disable_artificial_haplotype_recovery ? "true" : "false") << ","
        << "\"graph_enable_legacy_cycle_detection\":"
        << (options.calling.graph_enable_legacy_cycle_detection ? "true" : "false") << ","
        << "\"graph_min_dangling_branch_length\":"
        << options.calling.graph_min_dangling_branch_length << ","
        << "\"graph_min_dangling_matching_bases\":"
        << options.calling.graph_min_dangling_matching_bases << ","
        << "\"graph_allow_non_unique_kmers_in_ref\":"
        << (options.calling.graph_allow_non_unique_kmers_in_ref ? "true" : "false") << ","
        << "\"graph_recover_all_dangling_branches\":"
        << (options.calling.graph_recover_all_dangling_branches ? "true" : "false") << ","
        << "\"graph_dont_increase_kmer_sizes_for_cycles\":"
        << (options.calling.graph_dont_increase_kmer_sizes_for_cycles ? "true" : "false") << ","
        << "\"graph_kmer_size_selected\":" << result.graph_kmer_size_selected << ","
        << "\"graph_kmer_iterations\":" << result.graph_kmer_iterations << ","
        << "\"graph_has_non_reference_cycles\":"
        << (result.graph_has_non_reference_cycles ? "true" : "false") << ","
        << "\"graph_reference_non_unique_kmers\":"
        << result.graph_reference_non_unique_kmers << ","
        << "\"graph_reference_kmer_rejected\":"
        << (result.graph_reference_kmer_rejected ? "true" : "false") << ","
        << "\"error_correct_reads\":"
        << (options.calling.error_correct_reads ? "true" : "false") << ","
        << "\"error_correction_kmer_length\":"
        << options.calling.error_correction_kmer_length << ","
        << "\"error_correction_min_solid_observations\":"
        << options.calling.error_correction_min_solid_observations << ","
        << "\"error_correction_max_mismatches\":2,"
        << "\"error_correction_max_sparse_observations\":1,"
        << "\"error_correction_pileup_log_odds\":"
        << (std::isfinite(options.calling.pileup_error_correction_log_odds)
            ? std::to_string(options.calling.pileup_error_correction_log_odds) : "null") << ","
        << "\"graph_max_paths\":" << options.calling.graph_max_paths << ","
        << "\"graph_max_depth\":" << options.calling.graph_max_depth << ","
        << "\"max_haplotype_combination_alleles\":"
        << options.calling.max_haplotype_combination_alleles << ","
        << "\"haplotype_combination_blocks\":"
        << result.pairhmm_haplotype_combination_blocks << ","
        << "\"haplotype_pruning_log10\":" << options.calling.haplotype_pruning_log10 << ","
        << "\"biological_variant_calling\":\"prototype\",\"bit_identical_to_gatk\":false},"
        << "\"outputs\":[{\"path\":\"" << json_escape(options.output) << "\",\"kind\":\""
        << (options.gvcf ? "gvcf" : (is_vcf_output(options.output) ? "vcf" : "json"))
        << "\",\"complete\":" << (complete ? "true" : "false") << "}"
        << (options.index_path.empty() ? "" : ",{\"path\":\"" + json_escape(options.index_path) +
             "\",\"kind\":\"vcf-index\",\"complete\":" + (complete ? "true" : "false") + "}")
        << (options.assembly_region_out.empty() ? "" : ",{\"path\":\"" +
             json_escape(options.assembly_region_out) +
             "\",\"kind\":\"assembly-region-igv\",\"complete\":" +
             (file_complete(options.assembly_region_out) ? "true}" : "false}"))
        << "],"
        << "\"telemetry\":" << telemetry << "}\n";
    return out.str();
}

void write_text(const std::string& path, const std::string& text) {
    if (path == "-") return;
    if (has_suffix(path, ".gz")) {
#if FASTGATK_HAS_HTSLIB
        BGZF* file = bgzf_open(path.c_str(), "w");
        if (!file) throw std::runtime_error("cannot open BGZF output: " + path);
        const auto written = bgzf_write(file, text.data(), text.size());
        const auto close_status = bgzf_close(file);
        if (written != static_cast<ssize_t>(text.size()) || close_status != 0)
            throw std::runtime_error("cannot write BGZF output: " + path);
#elif FASTGATK_HAS_ZLIB
        gzFile file = gzopen(path.c_str(), "wb");
        if (!file) throw std::runtime_error("cannot open compressed output: " + path);
        const int written = gzwrite(file, text.data(), static_cast<unsigned int>(text.size()));
        const int close_status = gzclose(file);
        if (written != static_cast<int>(text.size()) || close_status != Z_OK)
            throw std::runtime_error("cannot write compressed output: " + path);
#else
        throw std::runtime_error("compressed output requires zlib: " + path);
#endif
    } else {
        std::ofstream file(path);
        if (!file) throw std::runtime_error("cannot open output: " + path);
        file << text;
        if (!file) throw std::runtime_error("cannot write output: " + path);
    }
}

class IncrementalVcfWriter {
public:
    explicit IncrementalVcfWriter(const std::string& path) : path_(path), compressed_(has_suffix(path, ".gz")) {
        if (path_ == "-") throw std::invalid_argument("--stream-by-contig requires a file output");
        if (compressed_) {
#if FASTGATK_HAS_HTSLIB
            bgzf_ = bgzf_open(path_.c_str(), "w");
            if (!bgzf_) throw std::runtime_error("cannot open BGZF output: " + path_);
#elif FASTGATK_HAS_ZLIB
            gz_ = gzopen(path_.c_str(), "wb");
            if (!gz_) throw std::runtime_error("cannot open compressed output: " + path_);
#else
            throw std::runtime_error("compressed output requires zlib: " + path_);
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
#elif FASTGATK_HAS_ZLIB
        const auto written = gzwrite(gz_, text.data(), static_cast<unsigned int>(text.size()));
        if (written != static_cast<int>(text.size()))
            throw std::runtime_error("cannot write compressed output: " + path_);
#endif
    }

    void close() {
        if (closed_) return;
        if (!compressed_) {
            plain_.close();
            if (!plain_) throw std::runtime_error("cannot finalize output: " + path_);
        } else {
#if FASTGATK_HAS_HTSLIB
            if (bgzf_close(bgzf_) != 0) throw std::runtime_error("cannot finalize BGZF output: " + path_);
            bgzf_ = nullptr;
#elif FASTGATK_HAS_ZLIB
            if (gzclose(gz_) != Z_OK) throw std::runtime_error("cannot finalize compressed output: " + path_);
            gz_ = nullptr;
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
#if FASTGATK_HAS_ZLIB
    gzFile gz_ = nullptr;
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

// Region streaming renders one independent GVCF tile at a time.  A tile
// boundary is not a biological boundary, however: GATK's HomRefBlockCombiner
// is allowed to merge adjacent simple REF/<NON_REF> records when their GQ
// bands agree.  Keep one pending block across tiles and apply the same
// element-wise PL/min-DP/rounded-median-DP reductions used by the in-memory
// caller.  Candidate records (and bp-resolution records without END) flush
// the pending block and are passed through byte-for-byte.
std::vector<std::string> split_tab_fields(const std::string& line) {
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

std::vector<std::string> split_colon_fields(const std::string& value) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find(':', begin);
        fields.push_back(value.substr(begin, end == std::string::npos
            ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

std::vector<std::string> split_comma_fields(const std::string& value) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find(',', begin);
        fields.push_back(value.substr(begin, end == std::string::npos
            ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

bool parse_stream_int(const std::string& text, std::int64_t& value) {
    if (text.empty() || text == ".") return false;
    try {
        std::size_t consumed = 0;
        const auto parsed = std::stoll(text, &consumed);
        if (consumed != text.size()) return false;
        value = parsed;
        return true;
    } catch (...) {
        return false;
    }
}

std::uint8_t stream_reference_confidence_band(const std::uint8_t gq,
                                              const std::vector<int>& gq_bands) {
    // Keep the streaming reducer's partition identity exactly aligned with
    // calling_pipeline.cpp.  GATK's -GQB values are exclusive upper bounds;
    // the stored identity is the inclusive lower bound.  In particular, GQ
    // 60 starts the 60-70 block and must not merge with GQ 59's 59-60 block.
    if (!gq_bands.empty()) {
        int lower = 0;
        for (const auto upper : gq_bands) {
            if (upper <= lower) continue;
            if (static_cast<int>(gq) < upper)
                return static_cast<std::uint8_t>(std::clamp(lower, 0, 99));
            lower = upper;
            if (upper >= 100) break;
        }
        return static_cast<std::uint8_t>(std::clamp(lower, 0, 99));
    }
    if (gq < 60) return gq;
    if (gq < 70) return 60;
    if (gq < 80) return 70;
    if (gq < 90) return 80;
    if (gq < 99) return 90;
    return 99;
}

std::uint32_t stream_rounded_median(
    const std::map<std::uint32_t, std::uint64_t>& histogram,
    const std::uint64_t count) {
    if (histogram.empty() || count == 0) return 0;
    const auto kth = [&](const std::uint64_t index) {
        std::uint64_t seen = 0;
        for (const auto& [depth, frequency] : histogram) {
            if (index < seen + frequency) return depth;
            seen += frequency;
        }
        return histogram.rbegin()->first;
    };
    const auto middle = count / 2;
    if ((count & 1U) != 0U) return kth(middle);
    const auto sum = static_cast<std::uint64_t>(kth(middle - 1)) +
                     static_cast<std::uint64_t>(kth(middle));
    return static_cast<std::uint32_t>((sum + 1U) / 2U);
}

std::string join_tab_fields(const std::vector<std::string>& fields) {
    std::ostringstream out;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) out << '\t';
        out << fields[index];
    }
    return out.str();
}

std::string join_colon_fields(const std::vector<std::string>& fields) {
    std::ostringstream out;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) out << ':';
        out << fields[index];
    }
    return out.str();
}

std::string replace_info_end(const std::string& info, const std::int64_t end) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    bool replaced = false;
    while (begin <= info.size()) {
        const auto stop = info.find(';', begin);
        auto item = info.substr(begin, stop == std::string::npos
            ? std::string::npos : stop - begin);
        if (item.rfind("END=", 0) == 0) {
            item = "END=" + std::to_string(end);
            replaced = true;
        }
        fields.push_back(std::move(item));
        if (stop == std::string::npos) break;
        begin = stop + 1;
    }
    if (!replaced) fields.push_back("END=" + std::to_string(end));
    std::ostringstream out;
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) out << ';';
        out << fields[index];
    }
    return out.str();
}

struct PendingGvcfBlock {
    std::vector<std::string> columns;
    std::vector<std::string> format_keys;
    std::vector<std::string> sample_values;
    std::vector<int> pl;
    std::map<std::uint32_t, std::uint64_t> depth_histogram;
    std::uint64_t depth_count = 0;
    std::string chrom;
    std::string alt;
    std::string format;
    std::string gt;
    std::int64_t start = -1;
    std::int64_t end = -1;
    std::uint8_t gq_band = 0;
    std::uint32_t min_depth = 0;
};

class RegionGvcfStitcher {
public:
    explicit RegionGvcfStitcher(IncrementalVcfWriter& writer,
                                const std::vector<int>& gq_bands)
        : writer_(writer), gq_bands_(gq_bands) {}

    void consume(const std::string& body,
                 const fastgatk::calling::Result& tile,
                 const fastgatk::io::HeaderSummary& header) {
        std::istringstream lines(body);
        std::string line;
        while (std::getline(lines, line)) {
            if (line.empty()) continue;
            const auto columns = split_tab_fields(line);
            auto block = parse_block(columns, tile, header, gq_bands_);
            if (!block.has_value()) {
                flush();
                writer_.write(line + "\n");
                continue;
            }
            if (!pending_.has_value()) {
                pending_ = std::move(*block);
            } else if (can_merge(*pending_, *block)) {
                merge(*pending_, *block);
            } else {
                flush();
                pending_ = std::move(*block);
            }
        }
    }

    void finish() { flush(); }

private:
    static std::optional<PendingGvcfBlock> parse_block(
        const std::vector<std::string>& columns,
        const fastgatk::calling::Result& tile,
        const fastgatk::io::HeaderSummary& header,
        const std::vector<int>& gq_bands) {
        if (columns.size() < 10 || columns[4] != "<NON_REF>") return std::nullopt;
        std::int64_t position = -1;
        if (!parse_stream_int(columns[1], position) || position < 1) return std::nullopt;
        std::int64_t end = -1;
        std::size_t info_begin = 0;
        while (info_begin <= columns[7].size()) {
            const auto info_end = columns[7].find(';', info_begin);
            const auto item = columns[7].substr(info_begin, info_end == std::string::npos
                ? std::string::npos : info_end - info_begin);
            if (item.rfind("END=", 0) == 0) {
                parse_stream_int(item.substr(4), end);
                break;
            }
            if (info_end == std::string::npos) break;
            info_begin = info_end + 1;
        }
        if (end < position) return std::nullopt; // bp-resolution/candidate record
        const auto format_keys = split_colon_fields(columns[8]);
        const auto sample_values = split_colon_fields(columns[9]);
        auto key_index = [&](const char* key) {
            return std::find(format_keys.begin(), format_keys.end(), key);
        };
        const auto gt_iter = key_index("GT");
        const auto dp_iter = key_index("DP");
        const auto gq_iter = key_index("GQ");
        const auto min_dp_iter = key_index("MIN_DP");
        const auto pl_iter = key_index("PL");
        if (gt_iter == format_keys.end() || dp_iter == format_keys.end() ||
            gq_iter == format_keys.end() || min_dp_iter == format_keys.end() ||
            pl_iter == format_keys.end()) return std::nullopt;
        const auto index_of = [&](const auto iter) {
            return static_cast<std::size_t>(iter - format_keys.begin());
        };
        const auto value_at = [&](const auto iter) -> std::string {
            const auto index = index_of(iter);
            return index < sample_values.size() ? sample_values[index] : std::string{};
        };
        std::int64_t gq = 0;
        std::int64_t min_dp = 0;
        if (!parse_stream_int(value_at(gq_iter), gq) || !parse_stream_int(value_at(min_dp_iter), min_dp))
            return std::nullopt;
        std::vector<int> pl;
        for (const auto& value : split_comma_fields(value_at(pl_iter))) {
            std::int64_t parsed = 0;
            if (!parse_stream_int(value, parsed)) return std::nullopt;
            pl.push_back(static_cast<int>(parsed));
        }
        if (pl.empty()) return std::nullopt;

        PendingGvcfBlock result;
        result.columns = columns;
        result.format_keys = format_keys;
        result.sample_values = sample_values;
        result.pl = std::move(pl);
        result.chrom = columns[0];
        result.alt = columns[4];
        result.format = columns[8];
        result.gt = value_at(gt_iter);
        result.start = position;
        result.end = end;
        result.gq_band = stream_reference_confidence_band(
            static_cast<std::uint8_t>(std::clamp<std::int64_t>(gq, 0, 99)), gq_bands);
        result.min_depth = static_cast<std::uint32_t>(std::max<std::int64_t>(0, min_dp));

        std::int32_t tid = -1;
        for (std::size_t index = 0; index < header.contigs.size(); ++index)
            if (header.contigs[index] == result.chrom) {
                tid = static_cast<std::int32_t>(index);
                break;
            }
        const auto block = std::find_if(tile.reference_blocks.begin(), tile.reference_blocks.end(),
            [&](const auto& candidate) {
                return candidate.tid == tid && candidate.start + 1 == position &&
                       candidate.end + 1 == end;
            });
        std::vector<std::uint32_t> depths;
        if (block != tile.reference_blocks.end() && !block->depth_values.empty()) {
            depths = block->depth_values;
        } else {
            // This fallback is only for unusual renderers that omit the
            // internal block metadata.  It preserves a valid merged record;
            // normal HaplotypeCaller tiles always take the exact path above.
            std::int64_t depth = 0;
            if (parse_stream_int(value_at(dp_iter), depth))
                depths.assign(static_cast<std::size_t>(end - position + 1),
                                     static_cast<std::uint32_t>(std::max<std::int64_t>(0, depth)));
        }
        if (depths.empty()) return std::nullopt;
        for (const auto depth : depths) ++result.depth_histogram[depth];
        result.depth_count = depths.size();
        result.min_depth = result.depth_histogram.begin()->first;
        return result;
    }

    static bool can_merge(const PendingGvcfBlock& left, const PendingGvcfBlock& right) {
        return left.chrom == right.chrom && left.end + 1 == right.start &&
               left.gq_band == right.gq_band && left.alt == right.alt &&
               left.format == right.format && left.gt == right.gt &&
               left.pl.size() == right.pl.size();
    }

    static void set_sample_value(PendingGvcfBlock& block, const char* key,
                                 const std::string& value) {
        const auto iter = std::find(block.format_keys.begin(), block.format_keys.end(), key);
        if (iter == block.format_keys.end()) return;
        const auto index = static_cast<std::size_t>(iter - block.format_keys.begin());
        if (block.sample_values.size() <= index) block.sample_values.resize(index + 1, ".");
        block.sample_values[index] = value;
    }

    static void merge(PendingGvcfBlock& left, const PendingGvcfBlock& right) {
        left.end = right.end;
        for (const auto& [depth, frequency] : right.depth_histogram)
            left.depth_histogram[depth] += frequency;
        left.depth_count += right.depth_count;
        left.min_depth = left.depth_histogram.begin()->first;
        for (std::size_t index = 0; index < left.pl.size(); ++index)
            left.pl[index] = std::min(left.pl[index], right.pl[index]);
        const auto dp = stream_rounded_median(left.depth_histogram, left.depth_count);
        int gq = 99;
        for (std::size_t index = 1; index < left.pl.size(); ++index)
            gq = std::min(gq, left.pl[index]);
        gq = std::clamp(gq, 0, 99);
        left.columns[7] = replace_info_end(left.columns[7], left.end);
        set_sample_value(left, "DP", std::to_string(dp));
        set_sample_value(left, "GQ", std::to_string(gq));
        set_sample_value(left, "MIN_DP", std::to_string(left.min_depth));
        std::ostringstream pl;
        for (std::size_t index = 0; index < left.pl.size(); ++index) {
            if (index != 0) pl << ',';
            pl << left.pl[index];
        }
        set_sample_value(left, "PL", pl.str());
        left.columns[9] = join_colon_fields(left.sample_values);
    }

    void flush() {
        if (!pending_.has_value()) return;
        writer_.write(join_tab_fields(pending_->columns) + "\n");
        pending_.reset();
    }

    IncrementalVcfWriter& writer_;
    const std::vector<int>& gq_bands_;
    std::optional<PendingGvcfBlock> pending_;
};

// A contig is the safe serialization boundary for the opt-in streaming HC
// path: assembly, PairHMM and reference-confidence state may span any reads
// within one contig, but no emitted record can depend on a later contig.  The
// pipeline keeps the expensive caller work off the decode thread while the
// final sink still writes contigs in input/header order.
struct HcContigDecodedBatch {
    fastgatk::io::ReadBatch reads;
    std::int32_t tid = -1;
};

struct HcContigComputedBatch {
    fastgatk::calling::Result result;
    // Final VCF annotation needs the original Host read/CIGAR view after
    // joint allele subsetting.  Retain it only through this bounded stream
    // pipeline; Kokkos still owns all numerical likelihood work.
    fastgatk::io::ReadBatch reads;
    std::int32_t tid = -1;
};

struct HcContigEncodedBatch {
    fastgatk::calling::Result result;
    fastgatk::io::ReadBatch reads;
    std::int32_t tid = -1;
};

struct HcRegionDecodedBatch {
    fastgatk::io::ReadBatch reads;
    std::int32_t tid = -1;
    std::int64_t core_start = 0;
    std::int64_t core_end = 0;
    // Indexed tile visibility is clipped to the original input interval set.
    // This prevents a halo from consuming reads that aggregate traversal
    // would never see past a user-supplied -L boundary.
    std::vector<fastgatk::io::HtsInterval> locus_intervals;
};

struct HcRegionComputedBatch {
    fastgatk::calling::Result result;
    std::string body;
    std::int32_t tid = -1;
    std::int64_t core_start = 0;
    std::int64_t core_end = 0;
};

struct HcRegionEncodedBatch {
    fastgatk::calling::Result result;
    std::string body;
    std::int32_t tid = -1;
    std::int64_t core_start = 0;
    std::int64_t core_end = 0;
};

std::uint64_t calling_result_bytes(const fastgatk::calling::Result& result) {
    std::uint64_t bytes = sizeof(fastgatk::calling::Result);
    bytes += result.candidates.size() * sizeof(fastgatk::calling::AssemblyCandidate);
    bytes += result.likelihoods.size() * sizeof(fastgatk::calling::Likelihoods);
    bytes += result.calls.size() * sizeof(fastgatk::calling::GenotypeCall);
    bytes += result.reference_blocks.size() * sizeof(fastgatk::calling::ReferenceBlock);
    bytes += result.multiallelic_depths.size() * sizeof(fastgatk::calling::MultiallelicDepth);
    bytes += result.graph_haplotype_path_tids.size() * sizeof(std::int32_t);
    bytes += result.graph_haplotype_path_starts.size() * sizeof(std::int32_t);
    bytes += result.graph_haplotype_path_ends.size() * sizeof(std::int32_t);
    bytes += result.graph_haplotype_path_support.size() * sizeof(std::uint32_t);
    bytes += result.graph_haplotype_path_scores.size() * sizeof(double);
    bytes += result.graph_haplotype_path_has_non_reference_edge.size();
    bytes += result.graph_haplotype_path_alt_read_starts.size() * sizeof(std::uint32_t);
    bytes += result.graph_haplotype_path_alt_read_ends.size() * sizeof(std::uint32_t);
    bytes += result.graph_haplotype_path_alignment_scores.size() * sizeof(std::int32_t);
    bytes += result.graph_haplotype_path_alignment_offsets.size() * sizeof(std::int32_t);
    bytes += result.candidate_prior_het.size() * sizeof(double);
    bytes += result.candidate_prior_hom_alt.size() * sizeof(double);
    bytes += result.annotation_read_qualified.size();
    bytes += result.likelihood_read_names.size() * sizeof(std::string);
    bytes += result.allele_read_likelihoods.size() * sizeof(std::vector<double>);
    bytes += result.reference_read_likelihoods.size() * sizeof(std::vector<double>);
    bytes += result.spanning_deletion_read_likelihoods.size() * sizeof(std::vector<double>);
    bytes += result.non_ref_read_likelihoods.size() * sizeof(std::vector<double>);
    for (const auto& value : result.graph_haplotype_sequences) bytes += value.size();
    for (const auto& value : result.graph_haplotype_path_cigars) bytes += value.size();
    for (const auto& value : result.likelihood_read_names) bytes += value.size();
    for (const auto& values : result.allele_read_likelihoods) bytes += values.size() * sizeof(double);
    for (const auto& values : result.reference_read_likelihoods) bytes += values.size() * sizeof(double);
    for (const auto& values : result.spanning_deletion_read_likelihoods) bytes += values.size() * sizeof(double);
    for (const auto& values : result.non_ref_read_likelihoods) bytes += values.size() * sizeof(double);
    for (const auto& block : result.reference_blocks) {
        bytes += block.depth_values.size() * sizeof(std::uint32_t);
        bytes += block.genotype_pl.size() * sizeof(std::int32_t);
    }
    for (const auto& depth : result.multiallelic_depths) {
        bytes += depth.reference.size();
        for (const auto& alternate : depth.alternates) bytes += alternate.size();
        bytes += depth.depths.size() * sizeof(std::uint32_t);
    }
    return std::max<std::uint64_t>(1, bytes);
}

std::uint64_t hc_contig_decoded_bytes(const HcContigDecodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcContigDecodedBatch) + batch.reads.bytes());
}

std::uint64_t hc_contig_computed_bytes(const HcContigComputedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcContigComputedBatch) +
                                   batch.reads.bytes() + calling_result_bytes(batch.result));
}

std::uint64_t hc_contig_encoded_bytes(const HcContigEncodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcContigEncodedBatch) +
                                   batch.reads.bytes() + calling_result_bytes(batch.result));
}

std::uint64_t hc_region_decoded_bytes(const HcRegionDecodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcRegionDecodedBatch) + batch.reads.bytes());
}

std::uint64_t hc_region_computed_bytes(const HcRegionComputedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcRegionComputedBatch) + batch.body.size() +
                                   calling_result_bytes(batch.result));
}

std::uint64_t hc_region_encoded_bytes(const HcRegionEncodedBatch& batch) {
    return std::max<std::uint64_t>(1, sizeof(HcRegionEncodedBatch) + batch.body.size() +
                                   calling_result_bytes(batch.result));
}

void trim_stream_result_payload(fastgatk::calling::Result& result,
                                const bool keep_reference_blocks) {
    if (!keep_reference_blocks) result.reference_blocks.clear();
    result.graph_haplotype_sequences.clear();
    result.graph_haplotype_path_tids.clear();
    result.graph_haplotype_path_starts.clear();
    result.graph_haplotype_path_ends.clear();
    result.graph_haplotype_path_support.clear();
    result.graph_haplotype_path_scores.clear();
    result.graph_haplotype_path_has_non_reference_edge.clear();
    result.graph_haplotype_path_alt_read_starts.clear();
    result.graph_haplotype_path_alt_read_ends.clear();
    result.graph_haplotype_path_cigars.clear();
    result.graph_haplotype_path_alignment_scores.clear();
    result.graph_haplotype_path_alignment_offsets.clear();
    result.graph_haplotype_event_maps.clear();
    result.allele_read_likelihoods.clear();
    result.reference_read_likelihoods.clear();
    result.spanning_deletion_read_likelihoods.clear();
    result.multiallelic_depths.clear();
    result.non_ref_read_likelihoods.clear();
    result.candidate_prior_het.clear();
    result.candidate_prior_hom_alt.clear();
    result.candidates.clear();
    result.likelihoods.clear();
}

void restrict_stream_calls_to_core(fastgatk::calling::Result& result,
                                   const std::int32_t tid,
                                   const std::int64_t start,
                                   const std::int64_t end) {
    result.calls.erase(std::remove_if(result.calls.begin(), result.calls.end(),
        [&](const fastgatk::calling::GenotypeCall& call) {
            return call.candidate.tid != tid || call.candidate.position < start ||
                   call.candidate.position >= end;
        }), result.calls.end());
    result.variant_calls = result.calls.size();
}

// A region-streamed GVCF evaluates a halo so that assembly and PairHMM retain
// the same evidence at a core edge as an aggregate traversal.  The resulting
// reference-confidence blocks can therefore extend outside the serialization
// core.  Crop that Host-side writer state before rendering: otherwise every
// adjacent halo writes its overlapping HomRef blocks again and the final VCF
// is no longer coordinate sorted.  The numerical evidence remains untouched;
// only the VCF ownership domain changes here.
void restrict_stream_gvcf_to_core(fastgatk::calling::Result& result,
                                  const std::int32_t tid,
                                  const std::int64_t start,
                                  const std::int64_t end) {
    const auto owns_position = [&](const std::int32_t candidate_tid,
                                   const std::int32_t position) {
        return candidate_tid == tid && position >= start && position < end;
    };
    result.candidates.erase(std::remove_if(result.candidates.begin(), result.candidates.end(),
        [&](const fastgatk::calling::AssemblyCandidate& candidate) {
            return !owns_position(candidate.tid, candidate.position);
        }), result.candidates.end());
    result.calls.erase(std::remove_if(result.calls.begin(), result.calls.end(),
        [&](const fastgatk::calling::GenotypeCall& call) {
            return !owns_position(call.candidate.tid, call.candidate.position);
        }), result.calls.end());
    result.variant_calls = result.calls.size();

    std::vector<fastgatk::calling::ReferenceBlock> core_blocks;
    core_blocks.reserve(result.reference_blocks.size());
    for (auto block : result.reference_blocks) {
        if (block.tid != tid || block.end < start || block.start >= end) continue;
        const auto source_start = block.start;
        const auto cropped_start = static_cast<std::int32_t>(std::max<std::int64_t>(block.start, start));
        const auto cropped_end = static_cast<std::int32_t>(std::min<std::int64_t>(block.end, end - 1));
        if (cropped_start > cropped_end) continue;
        if (!block.depth_values.empty()) {
            const auto first = static_cast<std::size_t>(cropped_start - source_start);
            const auto count = static_cast<std::size_t>(cropped_end - cropped_start + 1);
            if (first < block.depth_values.size()) {
                const auto last = std::min(block.depth_values.size(), first + count);
                block.depth_values = std::vector<std::uint32_t>(
                    block.depth_values.begin() + static_cast<std::ptrdiff_t>(first),
                    block.depth_values.begin() + static_cast<std::ptrdiff_t>(last));
                if (!block.depth_values.empty()) {
                    std::vector<std::uint32_t> ordered = block.depth_values;
                    std::sort(ordered.begin(), ordered.end());
                    const auto middle = ordered.size() / 2;
                    block.depth = (ordered.size() & 1U) != 0U
                        ? ordered[middle]
                        : static_cast<std::uint32_t>((static_cast<std::uint64_t>(ordered[middle - 1]) +
                                                       static_cast<std::uint64_t>(ordered[middle]) + 1U) / 2U);
                    block.min_depth = ordered.front();
                }
            }
        }
        block.start = cropped_start;
        block.end = cropped_end;
        core_blocks.push_back(std::move(block));
    }
    result.reference_blocks.swap(core_blocks);
    result.reference_block_count = result.reference_blocks.size();
}

fastgatk::calling::Result run_contig_streaming(
    Options& options,
    const std::vector<fastgatk::io::HtsReader*>& readers,
    const std::vector<std::string>& references,
    const fastgatk::runtime::ResourceSnapshot& resources,
    const fastgatk::runtime::AdaptiveController& controller,
    fastgatk::runtime::BatchLimits& batch_limits) {
    if (readers.empty()) throw std::invalid_argument("--stream-by-contig requires an input");
    auto& reader = *readers.front();
    if (!is_vcf_output(options.output) || options.output == "-")
        throw std::invalid_argument("--stream-by-contig requires VCF/GVCF file output");

    IncrementalVcfWriter writer(options.output);
    fastgatk::calling::Result empty;
    const auto header_and_body = options.gvcf
        ? split_vcf_header_body(gvcf(reader, empty, options.sample_ploidy,
                                     options.calling.emit_reference_confidence_bp_resolution,
                                     options.calling.min_alt_support,
                                     options.add_output_vcf_command_line,
                                     options.calling.gvcf_gq_bands,
                                     options.sites_only_vcf_output, options.floor_blocks))
        : split_vcf_header_body(vcf_text(reader, empty, options.sample_ploidy,
                                         options.add_output_vcf_command_line,
                                         options.sites_only_vcf_output));
    writer.write(header_and_body.first);

    fastgatk::calling::Result aggregate;
    bool aggregate_initialized = false;
    // Keep a private reader for serialization.  The decode stage advances the
    // caller's reader concurrently, while VCF/GVCF rendering only needs the
    // immutable header and interval metadata.  Separate ownership avoids
    // making HTSlib reader state part of the pipeline's synchronization model.
    fastgatk::io::HtsReader render_reader(
        options.input, options.reference, options.regions, 1, options.sample_name,
        options.exclusions, options.interval_padding, options.exclusion_padding,
        options.interval_set_rule);
    // Coordinate-sorted shards are merged with one bounded decoded batch per
    // input. Ties use input order, giving deterministic output without
    // concatenating complete files and re-entering an already flushed contig.
    struct Cursor {
        std::optional<fastgatk::io::ReadBatch> batch;
        std::size_t index = 0;
        bool eof = false;
    };
    std::vector<Cursor> cursors(readers.size());
    fastgatk::io::ReadBatch current;
    current.offsets.push_back(0);
    std::set<std::int32_t> flushed_tids;
    std::int32_t current_tid = -1;
    std::int32_t last_position = -1;
    const auto safe_memory_budget = resources.safe_memory_budget_bytes();

    const auto refill = [&](const std::size_t input_index) {
        auto& cursor = cursors[input_index];
        if (cursor.eof) return false;
        if (cursor.batch.has_value() && cursor.index < cursor.batch->records()) return true;
        fastgatk::io::ReadBatch next;
        if (!readers[input_index]->next(next)) {
            cursor.batch.reset();
            cursor.index = 0;
            cursor.eof = true;
            return false;
        }
        cursor.batch = std::move(next);
        cursor.index = 0;
        return cursor.batch->records() != 0;
    };

    using Pipeline = fastgatk::runtime::ThreeStagePipeline<
        HcContigDecodedBatch, HcContigComputedBatch, HcContigEncodedBatch>;
    auto decode = [&]() -> std::optional<HcContigDecodedBatch> {
        while (true) {
            std::size_t selected_input = readers.size();
            for (std::size_t input_index = 0; input_index < readers.size(); ++input_index) {
                if (!refill(input_index)) continue;
                auto& candidate = *cursors[input_index].batch;
                if (candidate.tids.size() != candidate.records() ||
                    candidate.positions.size() != candidate.records())
                    throw std::runtime_error("BAD_INPUT: contig streaming requires coordinate metadata");
                const auto candidate_tid = candidate.tids[cursors[input_index].index];
                const auto candidate_position = candidate.positions[cursors[input_index].index];
                if (candidate_tid < 0 || candidate_position < 0)
                    throw std::runtime_error(
                        "BAD_INPUT: --stream-by-contig requires mapped coordinate-sorted reads");
                if (selected_input == readers.size()) {
                    selected_input = input_index;
                    continue;
                }
                auto& selected = *cursors[selected_input].batch;
                const auto selected_tid = selected.tids[cursors[selected_input].index];
                const auto selected_position = selected.positions[cursors[selected_input].index];
                if (std::tie(candidate_tid, candidate_position, input_index) <
                    std::tie(selected_tid, selected_position, selected_input))
                    selected_input = input_index;
            }
            if (selected_input == readers.size()) {
                if (current.records() == 0) return std::nullopt;
                HcContigDecodedBatch result{std::move(current), current_tid};
                current = fastgatk::io::ReadBatch{};
                current.offsets.push_back(0);
                current_tid = -1;
                last_position = -1;
                return result;
            }
            auto& cursor = cursors[selected_input];
            auto& batch = *cursor.batch;
            const auto pending_index = cursor.index;
            if (batch.tids.size() != batch.records() || batch.positions.size() != batch.records())
                throw std::runtime_error("BAD_INPUT: contig streaming requires coordinate metadata");
            const auto tid = batch.tids[pending_index];
            const auto position = batch.positions[pending_index];
            if (tid < 0 || position < 0)
                throw std::runtime_error("BAD_INPUT: --stream-by-contig requires mapped coordinate-sorted reads");
            if (current_tid < 0) current_tid = tid;
            if (tid != current_tid) {
                if (tid < current_tid || flushed_tids.count(tid) != 0)
                    throw std::runtime_error(
                        "BAD_INPUT: --stream-by-contig requires coordinate-sorted input without contig re-entry");
                flushed_tids.insert(current_tid);
                HcContigDecodedBatch result{std::move(current), current_tid};
                current = fastgatk::io::ReadBatch{};
                current.offsets.push_back(0);
                current_tid = -1;
                last_position = -1;
                return result;
            }
            if (last_position > position)
                throw std::runtime_error(
                    "BAD_INPUT: --stream-by-contig requires coordinate-sorted records within each contig");
            append_record(current, batch, pending_index,
                          static_cast<std::uint32_t>(selected_input));
            ++cursor.index;
            last_position = position;
            std::uint64_t pending_bytes = 0;
            for (const auto& pending : cursors) {
                if (!pending.batch.has_value() || pending.index >= pending.batch->records()) continue;
                const auto bytes = pending.batch->bytes();
                pending_bytes = bytes > std::numeric_limits<std::uint64_t>::max() - pending_bytes
                    ? std::numeric_limits<std::uint64_t>::max() : pending_bytes + bytes;
            }
            const auto current_bytes = current.bytes();
            const auto estimated_bytes = current_bytes >
                    std::numeric_limits<std::uint64_t>::max() - pending_bytes
                ? std::numeric_limits<std::uint64_t>::max() : current_bytes + pending_bytes;
            options.streamed_peak_host_bytes = std::max(options.streamed_peak_host_bytes,
                                                        estimated_bytes);
            if (safe_memory_budget != 0 && estimated_bytes > safe_memory_budget)
                throw std::runtime_error("RESOURCE_EXHAUSTED: streamed contig staging exceeds safe memory budget");
            if (cursor.index == cursor.batch->records()) {
                const auto batch_bytes = cursor.batch->bytes();
                fastgatk::runtime::RuntimeTelemetry runtime;
                runtime.host_bytes = batch_bytes;
                runtime.inflight_bytes = batch_bytes;
                runtime.compute_queue_empty = true;
                const auto pressure = safe_memory_budget != 0 && batch_bytes > safe_memory_budget
                    ? fastgatk::runtime::Pressure::HostMemory
                    : fastgatk::runtime::Pressure::Normal;
                const auto next_limits = controller.next(runtime, batch_limits, pressure);
                if (next_limits.max_reads < readers[selected_input]->batch_records()) {
                    readers[selected_input]->set_batch_records(next_limits.max_reads);
                    ++options.adaptive_batch_reductions;
                }
                batch_limits = next_limits;
                options.effective_batch_records = readers[selected_input]->batch_records();
                cursor.batch.reset();
                cursor.index = 0;
            }
        }
    };
    auto compute = [&](HcContigDecodedBatch decoded) -> std::optional<HcContigComputedBatch> {
        return HcContigComputedBatch{
            fastgatk::calling::run(decoded.reads, references, options.calling),
            std::move(decoded.reads), decoded.tid};
    };
    // Rendering is intentionally kept at the deterministic sink below.  The
    // encoded stage still owns the value transition, while the sink is the
    // sole writer and the only stage allowed to merge Result metadata.
    auto encode = [](HcContigComputedBatch computed) -> std::optional<HcContigEncodedBatch> {
        return HcContigEncodedBatch{
            std::move(computed.result), std::move(computed.reads), computed.tid};
    };
    auto sink = [&](HcContigEncodedBatch encoded) {
        auto rendered = options.gvcf
            ? gvcf(render_reader, encoded.result, options.sample_ploidy,
                   options.calling.emit_reference_confidence_bp_resolution,
                   options.calling.min_alt_support,
                   options.add_output_vcf_command_line,
                   options.calling.gvcf_gq_bands,
                   options.sites_only_vcf_output, options.floor_blocks, &encoded.reads,
                   options.calling.informative_read_overlap_margin,
                   options.max_alternate_alleles)
            : vcf_text(render_reader, encoded.result, options.sample_ploidy,
                       options.add_output_vcf_command_line,
                       options.sites_only_vcf_output, &encoded.reads,
                       options.calling.informative_read_overlap_margin,
                       options.max_alternate_alleles);
        const auto parts = split_vcf_header_body(rendered);
        writer.write(parts.second);
        merge_stream_result(aggregate, std::move(encoded.result), aggregate_initialized);
        ++options.streamed_contigs;
    };
    // A decoded tile is bounded by the scheduler/cgroup budget, while a
    // computed Result may legitimately be larger than its input because it
    // carries ragged read×haplotype likelihoods and graph provenance. Keep
    // separate queue limits instead of rejecting a valid bounded tile merely
    // because its derived result has a different representation size.
    const auto decoded_capacity = safe_memory_budget > 0
        ? std::max<std::uint64_t>(safe_memory_budget, 1)
        : 64ULL * 1024ULL * 1024ULL;
    const auto result_capacity = std::max<std::uint64_t>(
        64ULL * 1024ULL * 1024ULL, safe_memory_budget);
    Pipeline pipeline(
        Pipeline::Limits{decoded_capacity, result_capacity, result_capacity},
        std::move(decode), std::move(compute), std::move(encode), std::move(sink),
        hc_contig_decoded_bytes, hc_contig_computed_bytes, hc_contig_encoded_bytes);
    const auto pipeline_metrics = pipeline.run();
    options.streamed_peak_host_bytes = std::max(
        options.streamed_peak_host_bytes, pipeline_metrics.peak_decoded_bytes);
    writer.close();
    if (!aggregate_initialized || aggregate.reads == 0)
        throw std::runtime_error("BAD_INPUT: no reads selected");
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
    return aggregate;
}

// Region streaming is deliberately driven by the HTS interval/index layer,
// not by a second read decoder.  A bounded core is expanded by the same
// activity/genotyping context used by the caller, decoded through an indexed
// HtsReader, then evaluated with the core interval as the output domain.  The
// halo supplies reads that cross a core boundary; the core interval prevents
// duplicate candidate/reference-confidence records in adjacent tiles.
fastgatk::calling::Result run_region_streaming(
    Options& options,
    const fastgatk::io::HtsReader& metadata_reader,
    const std::vector<std::string>& references,
    const fastgatk::runtime::ResourceSnapshot& resources,
    const fastgatk::runtime::AdaptiveController& controller,
    fastgatk::runtime::BatchLimits& batch_limits) {
    if (!is_vcf_output(options.output) || options.output == "-")
        throw std::invalid_argument("--stream-by-region requires VCF/GVCF file output");
    IncrementalVcfWriter writer(options.output);
    fastgatk::calling::Result empty;
    const auto header_and_body = options.gvcf
        ? split_vcf_header_body(gvcf(metadata_reader, empty, options.sample_ploidy,
                                     options.calling.emit_reference_confidence_bp_resolution,
                                     options.calling.min_alt_support,
                                     options.add_output_vcf_command_line,
                                     options.calling.gvcf_gq_bands,
                                     options.sites_only_vcf_output, options.floor_blocks))
        : split_vcf_header_body(vcf_text(metadata_reader, empty, options.sample_ploidy,
                                         options.add_output_vcf_command_line,
                                         options.sites_only_vcf_output));
    writer.write(header_and_body.first);
    std::optional<RegionGvcfStitcher> gvcf_stitcher;
    if (options.gvcf) gvcf_stitcher.emplace(writer, options.calling.gvcf_gq_bands);

    struct CoreInterval {
        std::int32_t tid = -1;
        std::int64_t start = 0;
        std::int64_t end = 0;
    };
    std::vector<CoreInterval> cores;
    const auto append_cores = [&](std::int32_t tid, std::int64_t start,
                                  std::int64_t end) {
        if (tid < 0 || start < 0 || end <= start ||
            static_cast<std::size_t>(tid) >= metadata_reader.header().contig_lengths.size()) return;
        const auto length = metadata_reader.header().contig_lengths[static_cast<std::size_t>(tid)];
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
        for (std::size_t tid = 0; tid < metadata_reader.header().contig_lengths.size(); ++tid)
            append_cores(static_cast<std::int32_t>(tid), 0,
                         metadata_reader.header().contig_lengths[tid]);
    }
    if (!metadata_reader.exclusion_intervals().empty()) {
        std::vector<CoreInterval> remaining;
        for (const auto core : cores) {
            std::vector<CoreInterval> segments{core};
            for (const auto exclusion : metadata_reader.exclusion_intervals()) {
                if (exclusion.tid != core.tid || exclusion.end <= core.start ||
                    exclusion.start >= core.end)
                    continue;
                std::vector<CoreInterval> next;
                for (const auto segment : segments) {
                    if (exclusion.end <= segment.start || exclusion.start >= segment.end) {
                        next.push_back(segment);
                        continue;
                    }
                    if (segment.start < exclusion.start)
                        next.push_back(CoreInterval{segment.tid, segment.start,
                            std::min(segment.end, exclusion.start)});
                    if (exclusion.end < segment.end)
                        next.push_back(CoreInterval{segment.tid,
                            std::max(segment.start, exclusion.end), segment.end});
                }
                segments.swap(next);
                if (segments.empty()) break;
            }
            remaining.insert(remaining.end(), segments.begin(), segments.end());
        }
        cores.swap(remaining);
    }
    if (cores.empty()) throw std::runtime_error("BAD_INPUT: no reference intervals to stream");
    const auto safe_memory_budget = resources.safe_memory_budget_bytes();
    const auto hard_memory_limit = resources.effective_memory_limit_bytes();
    // Leave a small fixed allowance for the decoded-batch object and HTS
    // metadata that are not part of ReadBatch::bytes().  Without this margin
    // a synthetic 155K allocation repeatedly split the same halo down to
    // single-base cores even though the hard allocation could hold the
    // indivisible tile; that changed AssemblyRegion context and call output.
    const auto split_memory_budget = safe_memory_budget > 0
        ? safe_memory_budget + std::min<std::uint64_t>(4096, std::numeric_limits<std::uint64_t>::max() - safe_memory_budget)
        : 0;
    const auto effective_activity_propagation =
        static_cast<std::uint64_t>(options.calling.max_probability_propagation_distance) + 50U;
    const auto context = static_cast<std::int64_t>(std::max({
        static_cast<std::uint64_t>(options.calling.assembly_region_padding),
        effective_activity_propagation,
        static_cast<std::uint64_t>(options.calling.indel_padding_for_genotyping)}));
    const auto input_halo_intervals = [&](const CoreInterval& core) {
        const auto contig_length = metadata_reader.header().contig_lengths[
            static_cast<std::size_t>(core.tid)];
        const auto halo_start = std::max<std::int64_t>(0, core.start - context);
        const auto halo_end = std::min<std::int64_t>(contig_length, core.end + context);
        std::vector<fastgatk::io::HtsInterval> result;
        // Aggregate traversal applies -L to input selection, not only to
        // VCF serialization.  Preserve that input domain for every stream
        // tile: reads outside it can create an extra active region and alter
        // a core graph/EventMap through pruning.
        if (metadata_reader.intervals().empty()) {
            result.push_back(fastgatk::io::HtsInterval{core.tid, halo_start, halo_end});
            return result;
        }
        for (const auto& input : metadata_reader.intervals()) {
            if (input.tid != core.tid) continue;
            const auto start = std::max(halo_start, input.start);
            const auto end = std::min(halo_end, input.end);
            if (end > start)
                result.push_back(fastgatk::io::HtsInterval{core.tid, start, end});
        }
        return result;
    };
    fastgatk::calling::Result aggregate;
    bool aggregate_initialized = false;
    // Keep the interval queue ordered. If a tile exceeds the allocation
    // budget, its core is split into left/right children and requeued in that
    // order. The halo is recomputed for each child, so reads crossing the
    // split remain available to both children while emitted records stay in
    // core coordinates and are therefore not duplicated.
    std::deque<CoreInterval> pending_cores(cores.begin(), cores.end());
    using Pipeline = fastgatk::runtime::ThreeStagePipeline<
        HcRegionDecodedBatch, HcRegionComputedBatch, HcRegionEncodedBatch>;
    auto decode = [&]() -> std::optional<HcRegionDecodedBatch> {
        while (!pending_cores.empty()) {
            const auto core = pending_cores.front();
            pending_cores.pop_front();
            auto halo_intervals = input_halo_intervals(core);
            if (halo_intervals.empty()) continue;
            const auto& contig = metadata_reader.header().contigs[static_cast<std::size_t>(core.tid)];
            std::vector<std::string> tile_regions;
            tile_regions.reserve(halo_intervals.size());
            for (const auto& halo : halo_intervals) {
                tile_regions.push_back(contig + ":" + std::to_string(halo.start + 1) +
                                       "-" + std::to_string(halo.end));
            }
            fastgatk::io::ReadBatch all;
            all.offsets.push_back(0);
            bool split_for_memory = false;
            // Merge the indexed tile heads by coordinate. This preserves the
            // same deterministic read order as aggregate/contig streaming;
            // appending one complete shard after another can change overlap
            // quality correction and therefore the emitted annotations.
            std::vector<std::unique_ptr<fastgatk::io::HtsReader>> tile_readers;
            tile_readers.reserve(options.inputs.size());
            for (const auto& input : options.inputs) {
                tile_readers.push_back(std::make_unique<fastgatk::io::HtsReader>(
                    input, options.reference, tile_regions,
                    options.batch_records, options.sample_name,
                    options.exclusions, 0, options.exclusion_padding));
                if (!tile_readers.back()->indexed())
                    throw std::runtime_error(
                        "BACKEND_UNAVAILABLE: --stream-by-region requires an HTSlib BAM/CRAM index");
                options.stream_indexed = true;
            }
            struct TileCursor {
                std::optional<fastgatk::io::ReadBatch> batch;
                std::size_t index = 0;
                bool eof = false;
            };
            std::vector<TileCursor> cursors(tile_readers.size());
            const auto refill = [&](const std::size_t input_index) {
                auto& cursor = cursors[input_index];
                if (cursor.eof) return false;
                if (cursor.batch.has_value() && cursor.index < cursor.batch->records()) return true;
                fastgatk::io::ReadBatch next;
                if (!tile_readers[input_index]->next(next)) {
                    cursor.batch.reset();
                    cursor.index = 0;
                    cursor.eof = true;
                    return false;
                }
                cursor.batch = std::move(next);
                cursor.index = 0;
                return cursor.batch->records() != 0;
            };
            while (!split_for_memory) {
                std::size_t selected_input = tile_readers.size();
                for (std::size_t input_index = 0; input_index < tile_readers.size(); ++input_index) {
                    if (!refill(input_index)) continue;
                    auto& candidate = *cursors[input_index].batch;
                    if (candidate.tids.size() != candidate.records() ||
                        candidate.positions.size() != candidate.records())
                        throw std::runtime_error("BAD_INPUT: region streaming requires coordinate metadata");
                    const auto candidate_tid = candidate.tids[cursors[input_index].index];
                    const auto candidate_position = candidate.positions[cursors[input_index].index];
                    if (candidate_tid < 0 || candidate_position < 0)
                        throw std::runtime_error(
                            "BAD_INPUT: --stream-by-region requires mapped coordinate reads");
                    if (selected_input == tile_readers.size()) {
                        selected_input = input_index;
                        continue;
                    }
                    auto& selected = *cursors[selected_input].batch;
                    const auto selected_tid = selected.tids[cursors[selected_input].index];
                    const auto selected_position = selected.positions[cursors[selected_input].index];
                    if (std::tie(candidate_tid, candidate_position, input_index) <
                        std::tie(selected_tid, selected_position, selected_input))
                        selected_input = input_index;
                }
                if (selected_input == tile_readers.size()) break;
                std::uint64_t resident_pending = 0;
                for (const auto& cursor : cursors) {
                    if (!cursor.batch.has_value() || cursor.index >= cursor.batch->records()) continue;
                    const auto bytes = cursor.batch->bytes();
                    resident_pending = bytes > std::numeric_limits<std::uint64_t>::max() - resident_pending
                        ? std::numeric_limits<std::uint64_t>::max() : resident_pending + bytes;
                }
                const auto accumulated_bytes = all.bytes();
                const auto estimated_bytes = accumulated_bytes >
                        std::numeric_limits<std::uint64_t>::max() - resident_pending
                    ? std::numeric_limits<std::uint64_t>::max() : accumulated_bytes + resident_pending;
                options.streamed_peak_host_bytes = std::max(options.streamed_peak_host_bytes,
                                                            estimated_bytes);
                if (split_memory_budget != 0 && estimated_bytes > split_memory_budget) {
                    const auto span = core.end - core.start;
                    if (span <= 1) {
                        // A single-base core plus its required halo is
                        // indivisible.  Permit it when it still fits inside
                        // the hard scheduler/cgroup allocation; only fail
                        // closed when even that allocation cannot hold the
                        // decoded tile.
                        if (hard_memory_limit == 0 || estimated_bytes > hard_memory_limit)
                            throw std::runtime_error(
                                "RESOURCE_EXHAUSTED: single-base AssemblyRegion tile exceeds memory allocation");
                        return HcRegionDecodedBatch{std::move(all), core.tid,
                                                    core.start, core.end,
                                                    std::move(halo_intervals)};
                    }
                    const auto midpoint = core.start + span / 2;
                    if (midpoint <= core.start || midpoint >= core.end)
                        throw std::runtime_error(
                            "RESOURCE_EXHAUSTED: AssemblyRegion tile cannot be split further");
                    pending_cores.push_front(CoreInterval{core.tid, midpoint, core.end});
                    pending_cores.push_front(CoreInterval{core.tid, core.start, midpoint});
                    ++options.streamed_region_splits;
                    options.batch_records = std::max<std::size_t>(1, options.batch_records / 2);
                    split_for_memory = true;
                    break;
                }
                auto& cursor = cursors[selected_input];
                auto& batch = *cursor.batch;
                const auto record_index = cursor.index;
                append_record(all, batch, record_index,
                              static_cast<std::uint32_t>(selected_input));
                ++cursor.index;
                if (cursor.index == cursor.batch->records()) {
                    const auto batch_bytes = cursor.batch->bytes();
                    fastgatk::runtime::RuntimeTelemetry runtime;
                    runtime.host_bytes = batch_bytes;
                    runtime.inflight_bytes = batch_bytes;
                    runtime.compute_queue_empty = true;
                    const auto pressure = safe_memory_budget != 0 && batch_bytes > safe_memory_budget
                        ? fastgatk::runtime::Pressure::HostMemory
                        : fastgatk::runtime::Pressure::Normal;
                    const auto next_limits = controller.next(runtime, batch_limits, pressure);
                    if (next_limits.max_reads < tile_readers[selected_input]->batch_records()) {
                        tile_readers[selected_input]->set_batch_records(next_limits.max_reads);
                        ++options.adaptive_batch_reductions;
                    }
                    batch_limits = next_limits;
                    options.effective_batch_records = tile_readers[selected_input]->batch_records();
                    cursor.batch.reset();
                    cursor.index = 0;
                }
            }
            if (split_for_memory) {
                all.clear();
                continue;
            }
            return HcRegionDecodedBatch{std::move(all), core.tid, core.start, core.end,
                                        std::move(halo_intervals)};
        }
        return std::nullopt;
    };
    auto compute = [&](HcRegionDecodedBatch decoded) -> std::optional<HcRegionComputedBatch> {
        auto tile_options = options.calling;
        // Evaluate a VCF tile over its complete read/activity halo, then
        // retain only calls owned by the serialization core.  Constraining
        // `intervals` to the core before run() changes its parent span and
        // therefore its AssemblyRegion partitioning, which made an arbitrary
        // stream boundary alter the graph/EventMap for an interior allele.
        // Reference-confidence calls likewise run over the complete halo;
        // restrict_stream_gvcf_to_core() then gives each serialized locus a
        // single owner while RegionGvcfStitcher reconstructs GATK blocks.
        tile_options.locus_intervals = decoded.locus_intervals;
        tile_options.intervals = decoded.locus_intervals;
        if (tile_options.intervals.empty()) return std::nullopt;
        tile_options.interval_tid = decoded.tid;
        tile_options.interval_start = static_cast<std::int32_t>(
            tile_options.intervals.front().start);
        tile_options.interval_end = static_cast<std::int32_t>(
            tile_options.intervals.back().end);
        auto result = fastgatk::calling::run(decoded.reads, references, tile_options);
        if (options.gvcf)
            restrict_stream_gvcf_to_core(result, decoded.tid,
                                         decoded.core_start, decoded.core_end);
        else
            restrict_stream_calls_to_core(result, decoded.tid,
                                          decoded.core_start, decoded.core_end);
        const auto rendered = options.gvcf
            ? gvcf(metadata_reader, result, options.sample_ploidy,
                   options.calling.emit_reference_confidence_bp_resolution,
                   options.calling.min_alt_support,
                   options.add_output_vcf_command_line,
                   options.calling.gvcf_gq_bands,
                   options.sites_only_vcf_output, options.floor_blocks, &decoded.reads,
                   tile_options.informative_read_overlap_margin,
                   options.max_alternate_alleles)
            : vcf_text(metadata_reader, result, options.sample_ploidy,
                       options.add_output_vcf_command_line,
                       options.sites_only_vcf_output, &decoded.reads,
                       tile_options.informative_read_overlap_margin,
                       options.max_alternate_alleles);
        const auto parts = split_vcf_header_body(rendered);
        // Serialization consumes the ragged read×haplotype matrices.  Retain
        // only the reference-block metadata required by GVCF stitching before
        // the Result enters the next byte-bounded queue.
        trim_stream_result_payload(result, options.gvcf);
        return HcRegionComputedBatch{
            std::move(result), parts.second, decoded.tid, decoded.core_start, decoded.core_end};
    };
    auto encode = [](HcRegionComputedBatch computed) -> std::optional<HcRegionEncodedBatch> {
        return HcRegionEncodedBatch{
            std::move(computed.result), std::move(computed.body), computed.tid,
            computed.core_start, computed.core_end};
    };
    auto sink = [&](HcRegionEncodedBatch encoded) {
        if (gvcf_stitcher.has_value())
            gvcf_stitcher->consume(encoded.body, encoded.result, metadata_reader.header());
        else
            writer.write(encoded.body);
        merge_stream_result(aggregate, std::move(encoded.result), aggregate_initialized);
        ++options.streamed_regions;
    };
    // The scheduler's 80% safe budget (plus the fixed decoded-object margin)
    // drives tile splitting, but an
    // indivisible one-base halo can still carry a few bytes of fixed
    // ReadBatch/metadata overhead above that conservative target.  Give the
    // decoded queue the hard allocation as its lower bound so such a tile is
    // not rejected after the splitter has reached its minimum granularity;
    // the telemetry/split guard remains tied to the safe budget and therefore
    // never permits the stream to exceed the scheduler allocation.
    const auto stage_limit = hard_memory_limit > safe_memory_budget
        ? hard_memory_limit : safe_memory_budget;
    const auto stage_capacity = stage_limit > 0
        ? std::max<std::uint64_t>(stage_limit, 1)
        : 64ULL * 1024ULL * 1024ULL;
    Pipeline pipeline(
        Pipeline::Limits{stage_capacity, stage_capacity, stage_capacity},
        std::move(decode), std::move(compute), std::move(encode), std::move(sink),
        hc_region_decoded_bytes, hc_region_computed_bytes, hc_region_encoded_bytes);
    const auto pipeline_metrics = pipeline.run();
    if (gvcf_stitcher.has_value()) gvcf_stitcher->finish();
    writer.close();
    if (!aggregate_initialized || (aggregate.reads == 0 && !options.gvcf))
        throw std::runtime_error("BAD_INPUT: no reads selected");
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
    return aggregate;
}

bool write_vcf_index(const std::string& path, std::string& index_path) {
    if (!has_suffix(path, ".gz")) {
        index_path = path + ".idx";
        fastgatk::io::write_uncompressed_vcf_tribble_index(path, index_path);
        return true;
    }
#if FASTGATK_HAS_HTSLIB
    index_path = path + ".tbi";
    const auto status = tbx_index_build3(path.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf);
    if (status != 0) {
        index_path.clear();
        throw std::runtime_error("cannot build VCF tabix index: " + path);
    }
    return true;
#else
    (void)index_path;
    return false;
#endif
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        Options options = parse(argc, argv);
        options.requested_batch_records = options.batch_records;
        // Probe the scheduler/cgroup envelope before creating Kokkos.  This
        // keeps a process from oversubscribing a SLURM allocation when a
        // caller passes a larger local thread count.
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        // Select a safe initial decode batch before HTSlib allocates records.
        // The estimate is deliberately conservative; actual staging is still
        // checked after every decoded batch and returns RESOURCE_EXHAUSTED if
        // an unusually long read exceeds the budget.
        const auto budget = resources.budget();
        const auto controller = fastgatk::runtime::AdaptiveController{};
        const auto initial_limits = controller.initial(
            budget, fastgatk::runtime::WorkEstimate{
                static_cast<std::uint32_t>(std::min<std::size_t>(
                    options.batch_records, std::numeric_limits<std::uint32_t>::max())),
                1, 1024, 256, 2048});
        options.batch_records = std::min<std::size_t>(
            options.batch_records, initial_limits.max_reads);
        options.effective_batch_records = options.batch_records;
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        if (options.inputs.empty()) options.inputs.push_back(options.input);
        std::vector<std::unique_ptr<fastgatk::io::HtsReader>> input_readers;
        input_readers.reserve(options.inputs.size());
        for (const auto& input : options.inputs) {
            input_readers.push_back(std::make_unique<fastgatk::io::HtsReader>(
                input, options.reference, options.regions, options.batch_records,
                options.sample_name, options.exclusions, options.interval_padding,
                options.exclusion_padding, options.interval_set_rule));
        }
        auto& reader = *input_readers.front();
        if (!options.reference.empty() && !options.disable_sequence_dictionary_validation)
            validate_reference_dictionary(options.reference, reader.header());
        for (std::size_t input_index = 1; input_index < input_readers.size(); ++input_index)
            validate_input_headers(reader, *input_readers[input_index],
                                   !options.sample_name.empty());
        std::vector<fastgatk::io::HtsReader*> reader_ptrs;
        reader_ptrs.reserve(input_readers.size());
        for (auto& input_reader : input_readers) reader_ptrs.push_back(input_reader.get());
        auto batch_limits = initial_limits;
        batch_limits.max_reads = static_cast<std::uint32_t>(std::min<std::size_t>(
            reader.batch_records(), std::numeric_limits<std::uint32_t>::max()));
        const auto references = load_reference(options.reference, reader.header());
        set_calling_intervals(options, reader);
        options.calling.forced_alleles = load_forced_alleles(
            options.alleles, reader.header(), options.force_call_filtered_alleles);
        options.calling.sample_ploidy = options.sample_ploidy;
        options.calling.max_genotype_count = options.max_genotype_count;
        options.calling.emit_reference_confidence = options.gvcf;
        // An aggregate normal VCF is written only after independent
        // AssemblyRegion results are reduced.  Preserve the region-owned
        // PairHMM matrices through that reduction so vcf_text() can recover
        // the exact joint PLs.  The Results are moved, not copied, and are
        // released immediately after this invocation has been serialized.
        options.calling.retain_assembly_region_likelihood_results =
            options.calling.retain_assembly_region_likelihood_results ||
            is_vcf_output(options.output);
        if (options.stream_by_contig && options.stream_by_region)
            throw std::invalid_argument("--stream-by-contig and --stream-by-region are mutually exclusive");
        // Region-tiled streaming deliberately re-evaluates a halo on both
        // sides of every tile.  It therefore cannot emit one authoritative
        // AssemblyRegionWalker trace without a cross-tile profile state
        // machine.  Reject rather than write a misleading GATK-shaped file.
        if (options.stream_by_region && !options.assembly_region_out.empty())
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: --assembly-region-out requires aggregate or --stream-by-contig traversal");
        if (options.stream_by_region) {
            const auto result = run_region_streaming(options, reader, references, resources,
                                                     controller, batch_limits);
            const auto summary = summary_json(options, reader, result, resources);
            if (is_vcf_output(options.output) && options.create_output_variant_index)
                write_vcf_index(options.output, options.index_path);
            if (!output_contract_valid(options))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: primary output or index is incomplete");
            const auto manifest_path = options.manifest.empty() && is_vcf_output(options.output)
                ? options.output + ".manifest.json" : options.manifest;
            if (!manifest_path.empty()) write_text(manifest_path,
                                                   make_manifest(options, reader, result, resources));
            validate_published_manifest(options, manifest_path);
            if (!options.telemetry.empty()) write_text(options.telemetry, summary);
            std::cout << summary;
            Kokkos::finalize();
            return 0;
        }
        if (options.stream_by_contig) {
            const auto result = run_contig_streaming(options, reader_ptrs, references, resources,
                                                     controller, batch_limits);
            write_assembly_region_igv(options.assembly_region_out, reader.header(), result);
            const auto summary = summary_json(options, reader, result, resources);
            if (is_vcf_output(options.output) && options.create_output_variant_index)
                write_vcf_index(options.output, options.index_path);
            if (!output_contract_valid(options))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: primary output or index is incomplete");
            const auto manifest_path = options.manifest.empty() && is_vcf_output(options.output)
                ? options.output + ".manifest.json" : options.manifest;
            if (!manifest_path.empty()) write_text(manifest_path,
                                                   make_manifest(options, reader, result, resources));
            validate_published_manifest(options, manifest_path);
            if (!options.telemetry.empty()) write_text(options.telemetry, summary);
            std::cout << summary;
            Kokkos::finalize();
            return 0;
        }
        fastgatk::io::ReadBatch all;
        all.offsets.push_back(0);
        const auto safe_memory_budget = resources.safe_memory_budget_bytes();
        // The aggregate path also performs a coordinate k-way merge. A plain
        // file-by-file append would make multi-shard input order-dependent and
        // could feed an unsorted stream into AssemblyRegion traversal. One
        // decoded batch per shard keeps the extra resident memory bounded.
        struct AggregateCursor {
            std::optional<fastgatk::io::ReadBatch> batch;
            std::size_t index = 0;
            bool eof = false;
        };
        std::vector<AggregateCursor> cursors(input_readers.size());
        const auto refill = [&](const std::size_t input_index) {
            auto& cursor = cursors[input_index];
            if (cursor.eof) return false;
            if (cursor.batch.has_value() && cursor.index < cursor.batch->records()) return true;
            fastgatk::io::ReadBatch next;
            if (!input_readers[input_index]->next(next)) {
                cursor.batch.reset();
                cursor.index = 0;
                cursor.eof = true;
                return false;
            }
            cursor.batch = std::move(next);
            cursor.index = 0;
            return cursor.batch->records() != 0;
        };
        const auto pending_bytes = [&]() {
            std::uint64_t total = 0;
            for (const auto& cursor : cursors) {
                if (!cursor.batch.has_value() || cursor.index >= cursor.batch->records()) continue;
                const auto bytes = cursor.batch->bytes();
                total = bytes > std::numeric_limits<std::uint64_t>::max() - total
                    ? std::numeric_limits<std::uint64_t>::max() : total + bytes;
            }
            return total;
        };
        while (true) {
            std::size_t selected_input = input_readers.size();
            for (std::size_t input_index = 0; input_index < input_readers.size(); ++input_index) {
                if (!refill(input_index)) continue;
                auto& candidate = *cursors[input_index].batch;
                if (candidate.tids.size() != candidate.records() ||
                    candidate.positions.size() != candidate.records())
                    throw std::runtime_error("BAD_INPUT: aggregate input requires coordinate metadata");
                const auto candidate_tid = candidate.tids[cursors[input_index].index];
                const auto candidate_position = candidate.positions[cursors[input_index].index];
                if (candidate_tid < 0 || candidate_position < 0)
                    throw std::runtime_error("BAD_INPUT: aggregate input requires mapped coordinate reads");
                if (selected_input == input_readers.size()) {
                    selected_input = input_index;
                    continue;
                }
                auto& selected = *cursors[selected_input].batch;
                const auto selected_tid = selected.tids[cursors[selected_input].index];
                const auto selected_position = selected.positions[cursors[selected_input].index];
                if (std::tie(candidate_tid, candidate_position, input_index) <
                    std::tie(selected_tid, selected_position, selected_input))
                    selected_input = input_index;
            }
            if (selected_input == input_readers.size()) break;
            const auto resident_pending = pending_bytes();
            const auto accumulated_bytes = all.bytes();
            const auto estimated_bytes = accumulated_bytes >
                    std::numeric_limits<std::uint64_t>::max() - resident_pending
                ? std::numeric_limits<std::uint64_t>::max()
                : accumulated_bytes + resident_pending;
            if (safe_memory_budget != 0 && estimated_bytes > safe_memory_budget)
                throw std::runtime_error("RESOURCE_EXHAUSTED: input staging exceeds safe memory budget");
            auto& cursor = cursors[selected_input];
            auto& batch = *cursor.batch;
            const auto record_index = cursor.index;
            append_record(all, batch, record_index,
                          static_cast<std::uint32_t>(selected_input));
            ++cursor.index;
            if (cursor.index == cursor.batch->records()) {
                const auto batch_bytes = cursor.batch->bytes();
                fastgatk::runtime::RuntimeTelemetry runtime;
                runtime.host_bytes = batch_bytes;
                runtime.inflight_bytes = batch_bytes;
                runtime.compute_queue_empty = true;
                const auto pressure = safe_memory_budget != 0 && batch_bytes > safe_memory_budget
                    ? fastgatk::runtime::Pressure::HostMemory
                    : fastgatk::runtime::Pressure::Normal;
                const auto next_limits = controller.next(runtime, batch_limits, pressure);
                if (next_limits.max_reads < input_readers[selected_input]->batch_records()) {
                    input_readers[selected_input]->set_batch_records(next_limits.max_reads);
                    ++options.adaptive_batch_reductions;
                }
                batch_limits = next_limits;
                options.effective_batch_records = input_readers[selected_input]->batch_records();
                cursor.batch.reset();
                cursor.index = 0;
            }
        }
        // Empty selected intervals are valid AssemblyRegionWalker input.
        // GATK's IntervalAlignmentContextIterator supplies zero pileups for
        // those coordinates, emits inactive max-sized profile chunks, and
        // --force-active may then send the chunks through callRegion().  Let
        // the shared Host/Kokkos calling path represent this instead of
        // rejecting the traversal before ActivityProfile is constructed.
        auto result = fastgatk::calling::run(all, references, options.calling);
        write_assembly_region_igv(options.assembly_region_out, reader.header(), result);
        const auto summary = summary_json(options, reader, result, resources);
        if (options.output == "-" || !is_vcf_output(options.output)) write_text(options.output, summary);
        else write_text(options.output, options.gvcf ? gvcf(reader, result, options.sample_ploidy,
                                                           options.calling.emit_reference_confidence_bp_resolution,
                                                           options.calling.min_alt_support,
                                                           options.add_output_vcf_command_line,
                                                           options.calling.gvcf_gq_bands,
                                                           options.sites_only_vcf_output, options.floor_blocks,
                                                           &all,
                                                           options.calling.informative_read_overlap_margin,
                                                           options.max_alternate_alleles) :
                        vcf_text(reader, result, options.sample_ploidy,
                                 options.add_output_vcf_command_line,
                                 options.sites_only_vcf_output, &all,
                                 options.calling.informative_read_overlap_margin,
                                 options.max_alternate_alleles));
        if (is_vcf_output(options.output) && options.create_output_variant_index)
            write_vcf_index(options.output, options.index_path);
        if (!output_contract_valid(options))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: primary output or index is incomplete");
        const auto manifest_path = options.manifest.empty() && is_vcf_output(options.output)
            ? options.output + ".manifest.json" : options.manifest;
        if (!manifest_path.empty()) write_text(manifest_path, make_manifest(options, reader, result, resources));
        validate_published_manifest(options, manifest_path);
        if (!options.telemetry.empty()) write_text(options.telemetry, summary);
        std::cout << summary;
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        return 2;
    }
}
