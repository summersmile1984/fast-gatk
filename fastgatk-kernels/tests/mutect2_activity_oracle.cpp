#include "fastgatk/kernels/activity_profile.hpp"

#include <Kokkos_Core.hpp>

#include <stdexcept>

namespace {

fastgatk::kernels::ActivityProfileInput reference_only_pileup() {
    fastgatk::kernels::ActivityProfileInput input;
    input.tids = {0};
    input.positions = {302};
    input.counts = {78, 0, 0, 0};
    input.reference_bases = {0};
    input.indel_counts = {0};
    input.quality_offsets = {0, 78};
    input.quality_bases.assign(78, 0);
    input.quality_values.assign(78, 40);
    input.quality_alt_flags.assign(78, 0);
    input.indel_quality_offsets = {0, 0};
    return input;
}

fastgatk::kernels::ActivityProfileOptions mutect2_options(const double initial_lod) {
    fastgatk::kernels::ActivityProfileOptions options;
    options.min_depth = 1;
    options.halo = 100;
    options.min_region_size = 50;
    options.max_region_size = 300;
    options.max_prob_propagation_distance = 50;
    options.quality_aware_somatic = true;
    options.initial_tumor_log10_odds = initial_lod;
    options.traversal_spans = {{0, 0, 600}};
    return options;
}

}  // namespace

int main() {
    Kokkos::initialize();
    try {
        // Mutect2Engine.isActive() always invokes logLikelihoodRatio.  With
        // nAlt=0 its flat-Beta entropy is -log(nRef + 1), rather than zero;
        // therefore a reference-only pileup remains inactive even at
        // --initial-tumor-lod 0.
        const auto options = mutect2_options(0.0);
        const auto reference_only = fastgatk::kernels::compute_activity_profile_kokkos(
            reference_only_pileup(), options);
        if (reference_only.active_loci != 0 || reference_only.active.size() != 1 ||
            reference_only.active[0] != 0 || reference_only.activity[0] != 0.0 ||
            reference_only.profile_regions.size() != 2 ||
            reference_only.profile_regions[0].start != 0 ||
            reference_only.profile_regions[0].end != 299 ||
            reference_only.profile_regions[1].start != 300 ||
            reference_only.profile_regions[1].end != 599 ||
            reference_only.profile_regions[0].active || reference_only.profile_regions[1].active)
            throw std::runtime_error("reference-only zero-LOD ActivityProfile differs from Mutect2");

        // A separate one-base indel bucket must not be treated as an
        // unconditional seed.  The GATK PileupQualBuffer likelihood for one
        // Q30 indel among 78 pileup elements remains negative at zero LOD.
        auto sparse_indel_input = reference_only_pileup();
        sparse_indel_input.indel_counts = {1};
        sparse_indel_input.indel_quality_offsets = {0, 1};
        sparse_indel_input.indel_quality_values = {30};
        const auto sparse_indel = fastgatk::kernels::compute_activity_profile_kokkos(
            sparse_indel_input, options);
        if (sparse_indel.active_loci != 0 || sparse_indel.active.size() != 1 ||
            sparse_indel.active[0] != 0 || sparse_indel.activity[0] != 0.0 ||
            sparse_indel.profile_regions.size() != 2)
            throw std::runtime_error("sparse indel zero-LOD ActivityProfile differs from Mutect2");

        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
