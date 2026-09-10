#include "fastgatk/kernels/pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>

#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

int main() {
    bool initialized = false;
    try {
        Kokkos::initialize();
        initialized = true;
        using Request = fastgatk::pairhmm::FragmentHaplotypeAggregationRequest;
        const double negative_infinity = -std::numeric_limits<double>::infinity();
        // Cells are (fragment A/B, haplotype REF0/REF1/ALT), flattened only
        // for this oracle. Interleaving requests also verifies stable
        // per-cell encounter order rather than sorted-value reassociation.
        const std::vector<Request> requests{
            {0, -1.0}, {3, -2.0}, {1, -2.0}, {4, -4.0}, {2, -0.25},
            {5, negative_infinity}, {0, negative_infinity}, {3, -3.0},
            {1, -3.0}, {4, -5.0}, {2, -0.75}, {5, -0.5},
            {6, -0.1}, {6, -0.2}, {6, -0.3}};
        const auto result =
            fastgatk::pairhmm::aggregate_fragment_haplotype_likelihoods_kokkos(
                requests, 7);
        if (result.sums_by_cell.size() != 7 || result.execution_space.empty())
            throw std::runtime_error("fragment aggregation oracle shape changed");
        const auto alleles =
            fastgatk::pairhmm::marginalize_read_allele_likelihoods_kokkos(
                {{0, 0, result.sums_by_cell[0]}, {0, 0, result.sums_by_cell[1]},
                 {0, 1, result.sums_by_cell[2]}, {1, 0, result.sums_by_cell[3]},
                 {1, 0, result.sums_by_cell[4]}, {1, 1, result.sums_by_cell[5]}},
                2);
        std::cout << "{\"status\":\"pass\",\"execution_space\":\""
                  << result.execution_space << "\",\"bits\":[";
        for (std::size_t cell = 0; cell < result.sums_by_cell.size(); ++cell) {
            if (cell != 0) std::cout << ',';
            std::cout << '"' << std::hex << std::setfill('0') << std::setw(16)
                      << std::bit_cast<std::uint64_t>(result.sums_by_cell[cell])
                      << std::dec << '"';
        }
        std::cout << "],\"allele_bits\":[";
        for (std::size_t cell = 0; cell < alleles.best_by_row_allele.size(); ++cell) {
            if (cell != 0) std::cout << ',';
            std::cout << '\"' << std::hex << std::setfill('0') << std::setw(16)
                      << std::bit_cast<std::uint64_t>(alleles.best_by_row_allele[cell])
                      << std::dec << '\"';
        }
        std::cout << "]}\n";
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
