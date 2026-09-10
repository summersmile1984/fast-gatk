#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "$0")/../.." && pwd)
cmake_bin="$project_root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake"
source_dir="$project_root/fastgatk-native"
build_dir="$source_dir/build"

bash "$source_dir/scripts/build_htslib.sh"

kokkos_arch_args=()
if [[ -n "${FASTGATK_KOKKOS_ARCH:-}" ]]; then
    if [[ ! "${FASTGATK_KOKKOS_ARCH}" =~ ^[A-Za-z0-9_]+$ ]]; then
        echo "FASTGATK_KOKKOS_ARCH must contain only letters, digits, and underscores" >&2
        exit 2
    fi
    # Kokkos owns the ISA selection. Examples are ZEN3/ZEN4 for AMD Zen
    # hosts, or SKX/ICL for Intel AVX-512 hosts. The same source/API is
    # rebuilt with the selected backend; no module-specific intrinsics are
    # introduced here.
    kokkos_arch_args+=("-DKokkos_ARCH_${FASTGATK_KOKKOS_ARCH}=ON")
fi

genomicsdb_bridge_args=()
if [[ -n "${FASTGATK_ENABLE_GENOMICSDB_BRIDGE:-}" ]]; then
    genomicsdb_bridge_args+=("-DFASTGATK_ENABLE_GENOMICSDB_BRIDGE=${FASTGATK_ENABLE_GENOMICSDB_BRIDGE}")
fi
if [[ -n "${FASTGATK_GENOMICSDB_PACKAGE:-}" ]]; then
    genomicsdb_bridge_args+=("-DFASTGATK_GENOMICSDB_PACKAGE=${FASTGATK_GENOMICSDB_PACKAGE}")
fi
if [[ -n "${FASTGATK_GENOMICSDB_LIBRARY:-}" ]]; then
    genomicsdb_bridge_args+=("-DFASTGATK_GENOMICSDB_LIBRARY=${FASTGATK_GENOMICSDB_LIBRARY}")
fi
if [[ -n "${FASTGATK_GENOMICSDB_JVM_LIBRARY:-}" ]]; then
    genomicsdb_bridge_args+=("-DFASTGATK_GENOMICSDB_JVM_LIBRARY=${FASTGATK_GENOMICSDB_JVM_LIBRARY}")
fi

# Build one artifact per Kokkos execution-space family.  Keeping this choice
# in the launcher makes the source/API identical across CPU and accelerator
# builds while avoiding accidental host-only defaults for a requested device
# backend.  Backend initialization/runtime fallback remains explicit in the
# dispatcher; this script only configures the artifact.
kokkos_backend="${FASTGATK_KOKKOS_BACKEND:-cpu}"
# Explicitly set every backend on each configure.  CMake caches options in the
# shared build directory; without the OFF values, switching from a CPU build
# to `serial` (or to a device backend) would silently retain the old backend.
kokkos_backend_args=(
    "-DKokkos_ENABLE_SERIAL=ON"
    "-DKokkos_ENABLE_OPENMP=OFF"
    "-DKokkos_ENABLE_CUDA=OFF"
    "-DKokkos_ENABLE_HIP=OFF"
    "-DKokkos_ENABLE_SYCL=OFF"
)
case "$kokkos_backend" in
    cpu|openmp)
        kokkos_backend_args+=("-DKokkos_ENABLE_OPENMP=ON")
        ;;
    serial)
        ;;
    cuda)
        kokkos_backend_args+=("-DKokkos_ENABLE_CUDA=ON")
        ;;
    hip)
        kokkos_backend_args+=("-DKokkos_ENABLE_HIP=ON")
        ;;
    sycl)
        kokkos_backend_args+=("-DKokkos_ENABLE_SYCL=ON")
        ;;
    cuda_openmp)
        kokkos_backend_args+=("-DKokkos_ENABLE_CUDA=ON" "-DKokkos_ENABLE_OPENMP=ON")
        ;;
    hip_openmp)
        kokkos_backend_args+=("-DKokkos_ENABLE_HIP=ON" "-DKokkos_ENABLE_OPENMP=ON")
        ;;
    sycl_openmp)
        kokkos_backend_args+=("-DKokkos_ENABLE_SYCL=ON" "-DKokkos_ENABLE_OPENMP=ON")
        ;;
    *)
        echo "FASTGATK_KOKKOS_BACKEND must be cpu|openmp|serial|cuda|hip|sycl|cuda_openmp|hip_openmp|sycl_openmp" >&2
        exit 2
        ;;
esac

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    "${kokkos_backend_args[@]}" \
    -DKokkos_ENABLE_TESTS=OFF \
    -DFAST_GATK_HTSLIB_ROOT="$project_root/third_party/htslib-build/htslib-src" \
    -DFAST_GATK_ZLIB_ROOT="$project_root/third_party/htslib-build/install" \
    "${genomicsdb_bridge_args[@]}" \
    "${kokkos_arch_args[@]}"
"$cmake_bin" --build "$build_dir" --parallel
