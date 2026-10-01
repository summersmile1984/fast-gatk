#pragma once

// Cross-backend (GPU) portability policy for fastgatk-kernels.
//
// fastgatk is a Kokkos-native implementation.  Code under this include
// tree must compile and run correctly on every Kokkos 5.2 backend
// (OpenMP / Threads / Serial / HPX / CUDA / HIP / SYCL) without source
// edits.
//
// Allowed:
//   - Kokkos::* APIs (Kokkos::View, Kokkos::parallel_for, Kokkos::resize,
//     Kokkos::view_alloc, Kokkos::Experimental::simd)
//   - Kokkos::DefaultExecutionSpace / Kokkos::DefaultHostMirrorSpace
//   - glibc-only host calls gated by `#if defined(__GLIBC__)` (e.g.
//     malloc_trim).  Non-glibc targets compile to no-op.
//
// Disallowed (each breaks at least one Kokkos backend):
//   - x86 SIMD intrinsics: <immintrin.h>, _mm_*, _mm256_*, _mm512_*,
//     __builtin_ia32_*, inline asm volatile.  Lock the build to x86 host
//     and bypass Kokkos SIMD dispatch.
//   - Host-only third-party libraries: WFA2-lib, jemalloc, GKL JNI,
//     Parasail, IntelSmithWaterman, libsched, libnuma.  Bind the build
//     to specific host implementations and (often) drop GPU support.
//   - Process-global state mutation: _mm_setcsr (FTZ/DAZ), std::cout
//     inside kernels, std::printf inside kernels.  Pollute sibling
//     streams or non-PairHMM float code.
//   - Plain new/delete/malloc/mmap outside Kokkos::* allocation hooks.
//     Kokkos OpenMP / CUDA / HIP / SYCL backends route through their
//     own allocator (cudaMalloc / hipMalloc / sycl::malloc_device); raw
//     host malloc bypasses it.
//   - OpenMP-only pragmas inside Kokkos kernels.  Use Kokkos
//     parallel_for with the right execution policy instead.
//
// The static assertions below fire at compile time when one of the
// disallowed patterns is included unguarded.  They are advisory —
// they don't catch every host-only construct, but they catch the
// categories that broke previous revisions.

#if defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) || \
    defined(__AVX__) || defined(__AVX2__) || defined(__AVX512F__) || \
    defined(__MMX__) || defined(__3dNOW__)
// Use _Pragma + message for C++20 portability.  The Kokkos build is
// configured with -std=c++20 -Wpedantic, and `#warning` is a C++23
// feature.
#pragma message("fastgatk-kernels targets Kokkos cross-backend.  "\
               "Compile-time SSE/AVX macros are set by the compiler, "\
               "but any actual _mm_* intrinsic or <immintrin.h> include "\
               "will break CUDA / HIP / SYCL portability.  See "\
               "fastgatk/kernels/gpu_safety.hpp.")
#endif