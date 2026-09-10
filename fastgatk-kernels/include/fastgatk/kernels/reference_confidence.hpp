#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::kernels {

// Per-observation reference-vs-any likelihood primitive used by GATK's
// ReferenceConfidenceModel.  Host code supplies only already projected,
// high-quality pileup observations; HTSlib/CIGAR and indel-informativeness
// policy remain Host responsibilities.  The output arrays are reduced on Host
// in input order so Strict mode keeps a deterministic summation order across
// Kokkos backends.
struct ReferenceConfidenceResult {
    // Per-observation likelihoods remain in the input order so Host can apply
    // the same deterministic genotype-combination rule for arbitrary sample
    // ploidy without moving CIGAR/HTSlib state into the kernel.
    std::vector<double> reference_likelihoods;
    std::vector<double> non_ref_likelihoods;
    std::vector<double> hom_ref;
    std::vector<double> het;
    std::vector<double> hom_alt;
    std::vector<std::uint32_t> depth;
    std::vector<std::uint32_t> reference_count;
    std::vector<std::uint32_t> non_ref_count;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

ReferenceConfidenceResult calculate_reference_confidence_kokkos(
    const std::vector<std::uint32_t>& locus_indices,
    const std::vector<std::uint8_t>& qualities,
    const std::vector<std::uint8_t>& is_alt,
    std::size_t locus_count);

// Materialize the biallelic REF-vs-any Number=G PL row for every locus from
// the projected per-observation likelihoods returned above.  The genotype
// index is the non-reference copy count (0..ploidy), which is the complete
// VCF ordering for two alleles at arbitrary fixed ploidy.  The mixture,
// hom-ref cap, PL rounding and GQ reduction execute in Kokkos so polyploid
// reference-confidence blocks do not fall back to a Host-only hot loop.
struct ReferenceConfidenceGenotypeResult {
    std::vector<std::int32_t> genotype_pl;  // locus-major, width ploidy + 1
    std::vector<std::uint8_t> gq;
    // Unrounded, uncapped hom-ref GQ in phred space.  GATK compares this
    // value against its indel cache before it materializes integer PLs.
    std::vector<double> hom_ref_gq_phred;
    int ploidy = 0;
    std::size_t locus_count = 0;
    double prepare_seconds = 0.0;
    double seconds = 0.0;
    std::string execution_space;
};

ReferenceConfidenceGenotypeResult calculate_reference_confidence_genotypes_kokkos(
    const std::vector<std::uint32_t>& locus_indices,
    const std::vector<double>& reference_likelihoods,
    const std::vector<double>& non_ref_likelihoods,
    std::size_t locus_count,
    int ploidy);

}  // namespace fastgatk::kernels
