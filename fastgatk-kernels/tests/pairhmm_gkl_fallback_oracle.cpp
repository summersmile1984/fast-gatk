#include "fastgatk/kernels/pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

constexpr std::size_t kLength = 120;
struct GklCase {
    std::size_t length;
    std::uint8_t base_quality;
    std::uint8_t insertion_quality;
    std::uint8_t deletion_quality;
    std::uint8_t continuation_quality;
    std::uint8_t read_base;
    std::uint8_t haplotype_base;
    double expected;
};

// Captured from GATK 4.6.2.0's bundled GKL 0.8.11 AVX-512 PairHMM.  All but
// the last row fall below GKL's Float32 scaled-sum threshold (1e-28f) and are
// therefore recomputed in Context<double>; the last row is deliberately just
// above the boundary and proves the ordinary Float32 path remains selected.
constexpr GklCase kCases[] = {
    {65, 30, 40, 40, 10, 'A', 'C', -70.520086228479270},
    {80, 30, 40, 40, 10, 'A', 'C', -85.520058668593600},
    {100, 30, 40, 40, 10, 'A', 'C', -105.52002563283520},
    {120, 30, 40, 40, 10, 'A', 'C', -125.51999489415991},
    {120, 20, 30, 35, 15, 'A', 'C', -182.21663105488668},
    {160, 25, 20, 50, 20, 'T', 'A', -314.17001793720806},
    {90, 35, 45, 20, 6, 'G', 'T', -61.402435302734375},
};

}  // namespace

int main() {
    bool initialized = false;
    try {
        std::vector<fastgatk::pairhmm::PairHmmRead> reads;
        std::vector<fastgatk::pairhmm::PairHmmHaplotype> haplotypes;
        std::vector<fastgatk::pairhmm::PairHmmRequest> requests;
        reads.reserve(std::size(kCases));
        haplotypes.reserve(std::size(kCases));
        requests.reserve(std::size(kCases));
        for (std::size_t index = 0; index < std::size(kCases); ++index) {
            const auto& item = kCases[index];
            fastgatk::pairhmm::PairHmmRead read;
            read.bases.assign(item.length, item.read_base);
            read.qualities.assign(item.length, item.base_quality);
            read.insertion_gop.assign(item.length, item.insertion_quality);
            read.deletion_gop.assign(item.length, item.deletion_quality);
            read.gap_continuation.assign(item.length, item.continuation_quality);
            reads.push_back(std::move(read));
            fastgatk::pairhmm::PairHmmHaplotype haplotype;
            haplotype.bases.assign(item.length, item.haplotype_base);
            haplotypes.push_back(std::move(haplotype));
            requests.push_back({static_cast<std::uint32_t>(index),
                                static_cast<std::uint32_t>(index)});
        }

        Kokkos::initialize();
        initialized = true;
        const auto result = fastgatk::pairhmm::compute_kokkos_bucketed(
            reads, haplotypes, requests, 1,
            fastgatk::pairhmm::PairHmmPrecision::Float32);
        if (result.likelihoods.size() != std::size(kCases))
            throw std::runtime_error("GKL PairHMM fallback result count changed");
        constexpr double kPrintedScientificTolerance = 5.0e-6;
        double maximum_delta = 0.0;
        for (std::size_t index = 0; index < std::size(kCases); ++index) {
            const auto observed = result.likelihoods[index];
            if (!std::isfinite(observed))
                throw std::runtime_error("GKL low-likelihood Float32 result did not take double fallback");
            const auto delta = std::abs(observed - kCases[index].expected);
            maximum_delta = std::max(maximum_delta, delta);
            if (delta > kPrintedScientificTolerance)
                throw std::runtime_error("GKL PairHMM fallback differs beyond GATK debug precision");
        }
        Kokkos::finalize();
        initialized = false;
        std::cout << "{\"status\":\"pass\",\"cases\":" << std::size(kCases)
                  << ",\"max_delta\":" << maximum_delta << "}\n";
        return 0;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
