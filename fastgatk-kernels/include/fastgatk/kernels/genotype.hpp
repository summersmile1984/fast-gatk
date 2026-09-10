#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// Deterministic PL materialization. PL is row-major by sample and uses the
// VCF genotype-combination ordering (nondecreasing allele tuple, grouped by
// the highest allele). The diploid case is the familiar triangular ordering
// index(a,b) = max(a,b) * (max(a,b) + 1) / 2 + min(a,b).
// The kernel deliberately returns -1 for a missing call/GQ rather than
// depending on HTSlib sentinel constants; Host code owns VCF encoding.
struct GenotypePlResult {
    // Generic sample-major output. For a diploid call this contains the same
    // values as first_allele/second_allele.
    std::vector<std::int32_t> alleles;
    int ploidy = 0;
    std::vector<std::int32_t> first_allele;
    std::vector<std::int32_t> second_allele;
    std::vector<std::int32_t> gq;
    // For posterior assignment, optional sample-major Number=G phred-scaled
    // normalized posteriors and the corresponding prior vector are returned
    // in addition to GT/GQ.  They remain empty for the ordinary PL path.
    std::vector<double> posterior_phred;
    std::vector<double> prior_phred;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Score candidate alleles for GATK's GenotypingEngine maximum-ALT boundary.
// For each sample the kernel finds the lowest-PL genotype (ties keep the
// lowest VCF genotype index), then adds the PL distance from hom-ref to every
// ALT present in that genotype.  This is the likelihood-sum definition used
// by AlleleSubsettingUtils.calculateMostLikelyAlleles; reduction to Host is
// deliberately sample-ordered so OpenMP/CUDA/HIP remain deterministic.
struct AlleleLikelihoodScoreResult {
    std::vector<double> scores;
    std::size_t sample_count = 0;
    int allele_count = 0;
    int ploidy = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

AlleleLikelihoodScoreResult calculate_allele_likelihood_scores_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy);

// Materialize a joint Number=G PL row from a reference-backed allele/read
// likelihood matrix.  `log10_likelihoods` is allele-major and contains the
// reference at allele 0 followed by every concrete (or symbolic) ALT; a
// non-finite cell means that the read cannot support that allele.  The
// genotype mixture and normalization run in the selected Kokkos execution
// space, so HC multi-ALT and polyploid writers do not reimplement a Host-only
// likelihood reduction.  The returned PL uses the same VCF colexicographic
// ordering as the GT/GQ APIs above.
struct JointGenotypePlResult {
    std::vector<std::int32_t> pl;
    int ploidy = 0;
    int allele_count = 0;
    std::size_t read_count = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

JointGenotypePlResult calculate_joint_genotype_pl_kokkos(
    const std::vector<double>& log10_likelihoods,
    std::size_t read_count,
    int allele_count,
    int ploidy,
    // GenotypeLikelihoods.getAsPLs preserves finite phred likelihoods above
    // 999 (for example, a high-depth HaplotypeCaller site can emit PL=1356).
    // Callers may supply a positive compatibility bound for a format that
    // explicitly requires one; zero preserves the normal GATK VCF behavior.
    int max_phred = 0);

// Materialize GATK GenotypePriorCalculator.assumingHW for every unphased
// genotype in VCF Number=G order.  The Host supplies one log10 heterozygous
// and homozygous prior per allele (reference at index 0); Kokkos performs the
// arbitrary-ploidy allele-count enumeration and returns the deterministic
// prior vector plus execution telemetry.  Keeping this enumeration beside the
// joint PL kernel avoids separate Host-only genotype-order implementations in
// HC and GenotypeGVCFs and makes the same path available to CUDA/HIP.
struct GenotypePriorResult {
    std::vector<double> log10_priors;
    int ploidy = 0;
    int allele_count = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

GenotypePriorResult calculate_genotype_priors_kokkos(
    const std::vector<double>& log10_heterozygous_priors,
    const std::vector<double>& log10_homozygous_priors,
    int ploidy);

GenotypePlResult derive_genotype_gt_gq_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy);

// Select the maximum posterior genotype from a sample-major PL matrix.  The
// caller supplies one log10 prior/posterior-frequency term per genotype row;
// the kernel combines it with the phred-scaled likelihood (log10L = -0.1 PL)
// and returns the same deterministic GT/GQ representation as the PL path.
// This is the explicit USE_POSTERIOR_PROBABILITIES compatibility path.  The
// default GenotypeGVCFs assignment remains PL/PREFER_PLS and therefore does
// not pay for this kernel unless requested.
GenotypePlResult derive_genotype_gt_gq_from_log10_priors_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& log10_genotype_priors);

GenotypePlResult derive_diploid_gt_gq_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count);

// Remap a sample-major VCF Number=G PL vector after an allele projection.
// `kept_alleles[new_index]` contains the source allele index.  A -1 mapping
// marks a target allele that is absent from this source shard; every target
// genotype containing it is emitted as the missing PL sentinel.  Mapped
// source tuples are sorted in the kernel before ranking, so the API also
// handles shards whose ALT ordering differs from the joint union.  The
// genotype ordering and ploidy are handled by the same Kokkos path used for
// GT/GQ derivation, so callers do not need a diploid-only triangular
// shortcut. Negative source PL values are preserved as missing sentinels.
struct RemapGenotypePlResult {
    std::vector<std::int32_t> pl;
    int ploidy = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

RemapGenotypePlResult remap_genotype_pl_kokkos(
    const std::vector<std::int32_t>& source_pl,
    std::size_t sample_count,
    int source_allele_count,
    int target_allele_count,
    int ploidy,
    const std::vector<std::int32_t>& kept_alleles);

// Remap a sample-major Number=R allele field (for example FORMAT/AD) with
// the same target-to-source allele map.  A -1 mapping emits the missing
// sentinel.  This keeps integer allele-field projection on the same Kokkos
// execution-space/telemetry path as Number=G PL projection.
struct RemapAlleleFieldResult {
    std::vector<std::int32_t> values;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

RemapAlleleFieldResult remap_allele_field_kokkos(
    const std::vector<std::int32_t>& source_values,
    std::size_t sample_count,
    int source_allele_count,
    int target_allele_count,
    const std::vector<std::int32_t>& target_to_source);

// Count decoded, sample-major genotype alleles.  Host code owns HTSlib GT
// sentinel encoding and passes -1 for missing/vector-end alleles.  The
// kernel returns deterministic integer counts and AN; no floating-point
// reduction or device-side VCF/HTSlib dependency is involved.
struct AlleleCountResult {
    std::vector<std::int32_t> counts;
    std::int32_t an = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

AlleleCountResult count_alleles_kokkos(
    const std::vector<std::int32_t>& allele_indices,
    std::size_t sample_count,
    int ploidy,
    int allele_count);

// Compute the posterior probability that every observed sample is hom-ref
// under a per-sample Hardy--Weinberg genotype prior.  The input PL values are
// phred-scaled likelihoods (lower is better) in the same sample-major VCF
// ordering used by derive_genotype_gt_gq_kokkos.  `log10_priors` contains one
// log10 prior per genotype row and is deliberately supplied by Host code so
// allele classification (SNP/INDEL/OTHER) remains independent of HTSlib and
// can be shared by CPU/GPU backends.  Missing rows are excluded from the
// product.  This is the reusable Kokkos posterior primitive used by the
// GenotypeGVCFs `--use-posteriors-to-calculate-qual` compatibility path; it is
// intentionally separate from GT/GQ selection, which remains PL-based like
// GATK's USE_PLS_TO_ASSIGN default.
struct GenotypePosteriorResult {
    double log10_no_variant = 0.0;
    double qual = 0.0;
    std::size_t samples_with_likelihoods = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

// Compute one-sample biallelic AF posteriors for a batch of candidate sites.
// `pl` is candidate-major (three diploid PL values per candidate) and
// `prior_pseudocounts` is candidate-major (REF, ALT).  The kernel mirrors
// GATK AlleleFrequencyCalculator: flat allele frequencies seed an EM loop,
// Dirichlet pseudocounts update the frequencies, and QUAL is the posterior
// probability that the candidate is not hom-ref. `log10_p_alt_absent` is
// the corresponding posterior of the all-reference state and lets the Host
// apply GATK's AF monomorphic gate without reconstructing that posterior.
// Its complement is retained for GenotypeGivenAlleles, whose forced
// monomorphic records report P(variant present) as their site QUAL.
// A negative PL row is treated as unavailable and returns QUAL=0 and an
// absent-allele posterior of zero. Keeping the loop in Kokkos
// avoids a Host-only confidence gate after PairHMM and works unchanged on
// Serial, OpenMP, CUDA and HIP execution spaces.
struct BiallelicCallConfidenceResult {
    std::vector<double> qual;
    std::vector<double> log10_p_alt_absent;
    std::vector<double> log10_p_variant_present;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

BiallelicCallConfidenceResult calculate_biallelic_call_confidence_kokkos(
    const std::vector<std::int32_t>& pl,
    const std::vector<double>& prior_pseudocounts);

GenotypePosteriorResult calculate_site_posterior_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& log10_priors,
    int spanning_deletion_index = -1);

// Cross-sample reference-confidence posterior using one shared genotype-prior
// vector for the entire cohort.  The kernel returns both the per-sample
// hom-reference log10 posterior and the deterministic joint probability that
// every informative sample is hom-reference.  This is deliberately separate
// from the optional site-QUAL primitive: GenotypeGVCFs can expose the joint
// reference-confidence value without changing PL-based GT/GQ assignment.
struct CrossSampleReferenceConfidenceResult {
    std::vector<double> sample_log10_p_reference;
    double joint_log10_p_reference = 0.0;
    double joint_qual = 0.0;
    std::size_t samples_with_likelihoods = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

CrossSampleReferenceConfidenceResult
calculate_cross_sample_reference_confidence_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& shared_log10_priors);

// GATK AlleleFrequencyCalculator-style cohort likelihood model.  The Host
// supplies the Dirichlet prior pseudocount for each allele (reference first),
// while Kokkos evaluates normalized genotype posteriors and effective allele
// counts for every sample.  EM starts from a flat allele-frequency prior and
// iterates until the largest effective-count change is <= 0.1, matching the
// GATK convergence rule.  PL/GT/GQ assignment is intentionally separate; the
// result is used for cohort QUAL and MLE AC/AF annotations.
struct AlleleFrequencyResult {
    std::vector<std::int32_t> integer_allele_counts;
    std::vector<double> effective_allele_counts;
    // Final EM log10 allele frequencies (reference first).  Host callers may
    // use these to reproduce the same cohort posterior genotype assignment
    // without moving the EM state or HTSlib objects onto the device.
    std::vector<double> log10_allele_frequencies;
    std::vector<double> log10_p_allele_absent;
    double log10_p_no_variant = 0.0;
    double qual = 0.0;
    std::size_t samples_with_likelihoods = 0;
    // Diploid hom-ref genotypes with GQ but no PL can be used by GATK's
    // AlleleFrequencyCalculator through GenotypeUtils' approximate PL path.
    // This reports how many samples took that compatibility path.
    std::size_t approximate_gq_samples = 0;
    int iterations = 0;
    bool converged = false;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

AlleleFrequencyResult calculate_allele_frequency_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& prior_pseudocounts);

// Variant of the cohort calculator that also accepts decoded sample GT
// alleles and FORMAT/GQ.  When a diploid sample has an entirely missing PL
// row but is hom-ref and has GQ, the kernel synthesizes the same approximate
// likelihoods as GATK's GenotypeUtils.makeApproximateDiploidLog10LikelihoodsFromGQ:
// hom-ref=0, every ref-containing heterozygote=GQ, and every hom-var row=
// round(30 / (-10*log10(0.5))) * GQ (=10*GQ).  Samples with non-hom-ref GT or
// unsupported ploidy remain excluded, matching
// genotypeIsUsableForAFCalculation.
AlleleFrequencyResult calculate_allele_frequency_kokkos(
    const std::vector<std::int32_t>& pl,
    std::size_t sample_count,
    int allele_count,
    int ploidy,
    const std::vector<double>& prior_pseudocounts,
    const std::vector<std::int32_t>& sample_gq,
    const std::vector<std::int32_t>& sample_alleles,
    int spanning_deletion_index = -1);

}  // namespace fastgatk::kernels
