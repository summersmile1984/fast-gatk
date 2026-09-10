cmake_minimum_required(VERSION 3.22)
project(fast_gatk_pairhmm_kokkos LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FAST_GATK_KOKKOS_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/kokkos-src/kokkos-5.2.0"
    CACHE PATH "Pinned Kokkos source tree")
add_subdirectory("${FAST_GATK_KOKKOS_SOURCE}" "${CMAKE_CURRENT_BINARY_DIR}/kokkos")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core"
                 "${CMAKE_CURRENT_BINARY_DIR}/fastgatk-core")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels"
                 "${CMAKE_CURRENT_BINARY_DIR}/fastgatk-kernels")

add_executable(pairhmm-kokkos
    src/main_kokkos.cpp)
target_include_directories(pairhmm-kokkos PRIVATE include)
target_link_libraries(pairhmm-kokkos PRIVATE fastgatk-kernels)
target_compile_options(pairhmm-kokkos PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
