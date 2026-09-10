#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "$0")/../.." && pwd)
cmake_bin="$project_root/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake"
source_dir="$project_root/pairhmm-demo"

build_variant() {
    local name=$1
    shift
    local build_dir="$source_dir/build-kokkos-$name"
    "$cmake_bin" -S "$source_dir" -B "$build_dir" \
        -DCMAKE_BUILD_TYPE=Release \
        -DKokkos_ENABLE_OPENMP=ON \
        -DKokkos_ENABLE_SERIAL=ON \
        -DKokkos_ENABLE_TESTS=OFF \
        "$@"
    "$cmake_bin" --build "$build_dir" --parallel
}

build_variant scalar
build_variant avx2 -DKokkos_ARCH_ZEN3=ON
build_variant avx512 -DKokkos_ARCH_ZEN4=ON
