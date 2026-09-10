#include "fastgatk/kernels/pairhmm.hpp"
#include "fastgatk/kernels/pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace fastgatk::pairhmm {
namespace {

class KokkosScope {
public:
    explicit KokkosScope(int threads) {
        if (Kokkos::is_initialized()) return;
        if (threads < 1) throw std::invalid_argument("threads must be positive");
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(threads);
        Kokkos::initialize(settings);
        owned_ = true;
    }

    KokkosScope(const KokkosScope&) = delete;
    KokkosScope& operator=(const KokkosScope&) = delete;

    ~KokkosScope() {
        if (owned_ && Kokkos::is_initialized()) Kokkos::finalize();
    }

private:
    bool owned_ = false;
};

BatchResult compute_compat(const std::vector<PairInput>& inputs, int threads, int iterations,
                           bool tristate_correction = true) {
    if (inputs.empty()) throw std::invalid_argument("PairHMM input batch must not be empty");
    KokkosScope scope(threads);
    PairIndexBatch pairs;
    pairs.read_ids.resize(inputs.size());
    pairs.haplotype_ids.resize(inputs.size());
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        pairs.read_ids[index] = static_cast<std::uint32_t>(index);
        pairs.haplotype_ids[index] = static_cast<std::uint32_t>(index);
    }
    const auto result = compute_kokkos(inputs, pairs, iterations, tristate_correction);
    BatchResult compatible;
    compatible.likelihoods = result.likelihoods;
    compatible.seconds = result.seconds;
    compatible.pairs_per_second = result.pairs_per_second;
    compatible.checksum = result.checksum;
    return compatible;
}

bool cpu_supports_avx2_impl() {
#if defined(__x86_64__) || defined(__i386__)
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
#else
    return false;
#endif
}

bool cpu_supports_avx512_impl() {
#if defined(__x86_64__) || defined(__i386__)
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f");
#else
    return false;
#endif
#else
    return false;
#endif
}

SimdBackend compiled_simd_backend_impl() {
#if defined(KOKKOS_ARCH_AVX512XEON)
    return SimdBackend::Avx512;
#elif defined(KOKKOS_ARCH_AVX2)
    return SimdBackend::Avx2;
#else
    return SimdBackend::Auto;
#endif
}

const char* compiled_simd_backend_name_impl() {
    switch (compiled_simd_backend_impl()) {
    case SimdBackend::Avx512: return "avx512";
    case SimdBackend::Avx2: return "avx2";
    case SimdBackend::Auto: return "scalar";
    }
    return "scalar";
}

}  // namespace

SimdBackend compiled_simd_backend() { return compiled_simd_backend_impl(); }

std::size_t compiled_simd_width() {
    return Kokkos::Experimental::simd<double>::size();
}

std::string compiled_simd_backend_name() {
    return compiled_simd_backend_name_impl();
}

double log10_likelihood_scalar(const PairInput& input, bool tristate_correction) {
    const auto result = compute_compat({input}, 1, 1, tristate_correction);
    return result.likelihoods.front();
}

BatchResult compute_scalar(const std::vector<PairInput>& inputs, int threads, int iterations) {
    return compute_compat(inputs, threads, iterations);
}

BatchResult compute_simd(const std::vector<PairInput>& inputs, int threads, int iterations,
                         SimdBackend backend) {
    const auto compiled = compiled_simd_backend_impl();
    if (backend != SimdBackend::Auto && backend != compiled) {
        throw std::runtime_error(
            "requested SIMD backend " +
            std::string(backend == SimdBackend::Avx512 ? "avx512" : "avx2") +
            " does not match the Kokkos ABI compiled into this binary (" +
            compiled_simd_backend_name_impl() + "); rebuild the matching target");
    }
    if (backend == SimdBackend::Avx2 && !cpu_supports_avx2_impl())
        throw std::runtime_error("AVX2 backend requested but avx2 is unavailable");
    if (backend == SimdBackend::Avx512 && !cpu_supports_avx512_impl())
        throw std::runtime_error("AVX-512 backend requested but avx512f is unavailable");
    // Kokkos owns the compile-time ABI selection.  The backend argument is
    // retained for source compatibility and runtime feature validation; it
    // never introduces an ISA-specific intrinsic into the production target.
    return compute_compat(inputs, threads, iterations);
}

std::string cpu_simd_description() {
    return std::string("Kokkos SIMD ABI=") + compiled_simd_backend_name_impl() +
           ", width=" + std::to_string(compiled_simd_width()) +
           "; CPU avx512=" + (cpu_supports_avx512_impl() ? "true" : "false") +
           ", avx2=" + (cpu_supports_avx2_impl() ? "true" : "false");
}

bool cpu_supports_avx2() { return cpu_supports_avx2_impl(); }
bool cpu_supports_avx512() { return cpu_supports_avx512_impl(); }

}  // namespace fastgatk::pairhmm
