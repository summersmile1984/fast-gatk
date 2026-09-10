#include "fastgatk/contamination_kernel_segmenter.hpp"

#include <Kokkos_Core.hpp>

#include <iostream>
#include <vector>

int main() {
    Kokkos::initialize();
    int status = 0;
    {
        std::vector<double> two_steps;
        for (int i = 0; i < 90; ++i) two_steps.push_back(i < 30 ? 0.1 : (i < 60 ? 0.5 : 0.2));
        const auto steps = fastgatk::contamination::find_changepoints_with_windows(
            two_steps, 10, 6, {8, 16}, 1.0, 1.0);
        const std::vector<std::size_t> expected_steps{29, 59};
        if (steps != expected_steps) {
            std::cerr << "two-steps changepoints mismatch\n";
            status = 1;
        }
        std::vector<double> flat(90, 0.3);
        if (!fastgatk::contamination::find_changepoints_with_windows(
                flat, 10, 6, {8, 16}, 1.0, 1.0).empty()) {
            std::cerr << "flat series unexpectedly segmented\n";
            status = 1;
        }
        std::vector<double> ramp;
        for (int i = 0; i < 90; ++i) ramp.push_back(static_cast<double>(i) / 90.0);
        const auto ramp_points = fastgatk::contamination::find_changepoints_with_windows(
            ramp, 10, 6, {8, 16}, 1.0, 1.0);
        // This API uses ContaminationSegmenter's folded minor-allele kernel;
        // the result is therefore not the generic (unfolded) Gaussian golden
        // used by the standalone Java KernelSegmenter tests.
        const std::vector<std::size_t> expected_ramp{18, 29, 60, 71};
        if (ramp_points != expected_ramp) {
            std::cerr << "ramp changepoints mismatch:";
            for (const auto point : ramp_points) std::cerr << ' ' << point;
            std::cerr << '\n';
            status = 1;
        }
    }
    Kokkos::finalize();
    return status;
}
