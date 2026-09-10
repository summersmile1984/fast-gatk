#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace fastgatk::pairhmm {

struct PairInput {
    std::vector<std::uint8_t> haplotype;
    std::vector<std::uint8_t> read;
    std::vector<std::uint8_t> read_qual;
    std::vector<std::uint8_t> insertion_gop;
    std::vector<std::uint8_t> deletion_gop;
    std::vector<std::uint8_t> gap_continuation;
};

struct BatchResult {
    std::vector<double> likelihoods;
    double seconds = 0.0;
    double pairs_per_second = 0.0;
    double checksum = 0.0;
};

enum class SimdBackend { Auto, Avx2, Avx512 };

// Return the SIMD ABI compiled into this binary by Kokkos.  Runtime CPU
// feature probes alone are insufficient: a binary built with an AVX2 (or
// AVX-512) Kokkos SIMD ABI cannot safely change width after it has been
// loaded.  Callers that expose an explicit backend selector should compare
// against this value and fail closed on an ABI mismatch.
SimdBackend compiled_simd_backend();
std::size_t compiled_simd_width();
std::string compiled_simd_backend_name();

double log10_likelihood_scalar(const PairInput& input, bool tristate_correction = true);
// Compatibility batch entry points.  In the native target these forward to
// the Kokkos implementation; the selected Kokkos SIMD ABI is the portability
// boundary.  The standalone pairhmm-demo Makefile retains the historical
// scalar/raw-SIMD oracle separately.
BatchResult compute_simd(const std::vector<PairInput>& inputs, int threads,
                         int iterations = 1,
                         SimdBackend backend = SimdBackend::Auto);
BatchResult compute_scalar(const std::vector<PairInput>& inputs, int threads,
                           int iterations = 1);
std::string cpu_simd_description();
bool cpu_supports_avx2();
bool cpu_supports_avx512();

}  // namespace fastgatk::pairhmm
