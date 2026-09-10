#include "fastgatk/calling/pipeline.hpp"

#include <Kokkos_Core.hpp>

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
    Kokkos::initialize();
    try {
        // AssemblyRegionWalker creates its unbounded traversal from
        // IntervalUtils.getAllIntervalsForReference().  The four reads below
        // begin at zero-based 8 on a ten-base reference and each projects to
        // position 15.  GATK evaluates only positions 8 and 9 in its
        // ActivityProfile; the overhanging read bases still belong to the
        // final padded AssemblyRegion, but cannot create pileup loci past the
        // reference sequence dictionary.
        fastgatk::io::ReadBatch reads;
        reads.offsets.push_back(0);
        reads.cigar_offsets.push_back(0);
        for (std::size_t record = 0; record < 4; ++record) {
            reads.bases.insert(reads.bases.end(), 8, 'A');
            reads.qualities.insert(reads.qualities.end(), 8, 30);
            reads.offsets.push_back(static_cast<std::uint32_t>(reads.bases.size()));
            reads.positions.push_back(8);
            reads.tids.push_back(0);
            reads.mapq.push_back(60);
            reads.flags.push_back(0);
            reads.cigar_ops.push_back(fastgatk::io::CigarOp{
                8, fastgatk::io::CigarOpCode::Match}.pack());
            reads.cigar_offsets.push_back(static_cast<std::uint32_t>(reads.cigar_ops.size()));
        }

        fastgatk::calling::Options options;
        options.min_depth = 1;
        options.min_alt_support = 2;
        options.min_base_quality = 10;
        options.graph_allow_non_unique_kmers_in_ref = true;
        const auto result = fastgatk::calling::run(reads, {std::string(10, 'A')}, options);
        if (result.activity_loci != 2)
            throw std::runtime_error(
                "CIGAR bases beyond the reference traversal end reached ActivityProfile");
        if (result.activity_profile_regions.empty())
            throw std::runtime_error("reference-bounded traversal produced no activity profile regions");
        for (const auto& region : result.activity_profile_regions) {
            if (region.tid != 0 || region.start < 0 || region.end >= 10)
                throw std::runtime_error("activity profile region escaped reference boundary");
        }

        std::cout << "{\"status\":\"pass\",\"activity_loci\":"
                  << result.activity_loci << "}\n";
        Kokkos::finalize();
        return 0;
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
}
