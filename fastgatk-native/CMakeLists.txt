cmake_minimum_required(VERSION 3.22)
project(fast_gatk_native LANGUAGES CXX)

include(GNUInstallDirs)

enable_testing()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(FAST_GATK_KOKKOS_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/kokkos-src/kokkos-5.2.0"
    CACHE PATH "Pinned Kokkos source tree")

# One source tree is compiled against the selected Kokkos execution space.
# AUTO preserves any Kokkos_* cache values supplied by a toolchain; explicit
# backends make reproducible Serial/OpenMP/CUDA/HIP/SYCL builds possible from
# the project interface without introducing backend-specific code paths.
set(FASTGATK_KOKKOS_BACKEND "AUTO" CACHE STRING
    "Kokkos backend: AUTO, SERIAL, OPENMP, CUDA, HIP, SYCL, or a *_OPENMP combination")
set_property(CACHE FASTGATK_KOKKOS_BACKEND PROPERTY STRINGS
    AUTO SERIAL OPENMP CUDA HIP SYCL CUDA_OPENMP HIP_OPENMP SYCL_OPENMP)
set(FASTGATK_KOKKOS_ARCH "AUTO" CACHE STRING
    "Kokkos architecture name (for example ZEN4, SKX, AMPERE80), or AUTO")
option(FASTGATK_INSTALL_NATIVE_TOOLS
    "Install native tool binaries and the relocatable gatk dispatcher" ON)
string(TOUPPER "${FASTGATK_KOKKOS_BACKEND}" FASTGATK_KOKKOS_BACKEND_UPPER)
string(TOUPPER "${FASTGATK_KOKKOS_ARCH}" FASTGATK_KOKKOS_ARCH_UPPER)

set(_fastgatk_kokkos_backends SERIAL OPENMP CUDA HIP SYCL)
set(_fastgatk_kokkos_backend_combinations CUDA_OPENMP HIP_OPENMP SYCL_OPENMP)
if (NOT FASTGATK_KOKKOS_ARCH_UPPER STREQUAL "AUTO")
    if (NOT FASTGATK_KOKKOS_ARCH_UPPER MATCHES "^[A-Z0-9_]+$")
        message(FATAL_ERROR
            "FASTGATK_KOKKOS_ARCH must contain only letters, digits, and underscores, or AUTO")
    endif()
    # Kokkos owns the architecture catalogue and compiler flags.  This
    # project-level alias keeps direct CMake invocations consistent with the
    # build_native.sh launcher without adding an ISA-specific code path.
    set(Kokkos_ARCH_${FASTGATK_KOKKOS_ARCH_UPPER} ON CACHE BOOL
        "Enable Kokkos architecture ${FASTGATK_KOKKOS_ARCH_UPPER}" FORCE)
endif()
if (NOT FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "AUTO")
    list(FIND _fastgatk_kokkos_backends "${FASTGATK_KOKKOS_BACKEND_UPPER}" _fastgatk_backend_index)
    list(FIND _fastgatk_kokkos_backend_combinations "${FASTGATK_KOKKOS_BACKEND_UPPER}" _fastgatk_combination_index)
    if (_fastgatk_backend_index EQUAL -1 AND _fastgatk_combination_index EQUAL -1)
        message(FATAL_ERROR
            "FASTGATK_KOKKOS_BACKEND must be AUTO, SERIAL, OPENMP, CUDA, HIP, SYCL, CUDA_OPENMP, HIP_OPENMP, or SYCL_OPENMP")
    endif()
    foreach (_backend IN LISTS _fastgatk_kokkos_backends)
        set(Kokkos_ENABLE_${_backend} OFF CACHE BOOL "Enable Kokkos ${_backend} backend" FORCE)
    endforeach()
    if (FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "CUDA" OR
        FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "CUDA_OPENMP")
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable Kokkos CUDA backend" FORCE)
    elseif (FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "HIP" OR
            FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "HIP_OPENMP")
        set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable Kokkos HIP backend" FORCE)
    elseif (FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "SYCL" OR
            FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "SYCL_OPENMP")
        set(Kokkos_ENABLE_SYCL ON CACHE BOOL "Enable Kokkos SYCL backend" FORCE)
    endif()
    if (FASTGATK_KOKKOS_BACKEND_UPPER MATCHES "_OPENMP$" OR
        FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "OPENMP")
        set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable Kokkos OpenMP backend" FORCE)
    endif()
    if (FASTGATK_KOKKOS_BACKEND_UPPER STREQUAL "SERIAL" OR
        FASTGATK_KOKKOS_BACKEND_UPPER MATCHES "^(CUDA|HIP|SYCL)$" OR
        FASTGATK_KOKKOS_BACKEND_UPPER MATCHES "_OPENMP$")
        # Device-only Kokkos builds still need a Host execution space for
        # decode/collect and deterministic output.  Serial is the portable
        # host fallback when no OpenMP host backend was requested.
        set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable Kokkos Serial backend" FORCE)
    endif()
endif()
add_subdirectory("${FAST_GATK_KOKKOS_SOURCE}" "${CMAKE_CURRENT_BINARY_DIR}/kokkos")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core"
                 "${CMAKE_CURRENT_BINARY_DIR}/fastgatk-core")

# Kokkos selects a device execution space before a host space when both are
# enabled.  Derive this from the configured targets so CTest and telemetry do
# not report a CUDA/HIP/SYCL build as Serial merely because Serial is also
# required for the Host decode/collect boundary.
set(FASTGATK_EXPECTED_EXECUTION_SPACE "Host")
if (Kokkos_ENABLE_CUDA)
    set(FASTGATK_EXPECTED_EXECUTION_SPACE "Cuda")
elseif (Kokkos_ENABLE_HIP)
    set(FASTGATK_EXPECTED_EXECUTION_SPACE "HIP")
elseif (Kokkos_ENABLE_SYCL)
    set(FASTGATK_EXPECTED_EXECUTION_SPACE "SYCL")
elseif (Kokkos_ENABLE_OPENMP)
    set(FASTGATK_EXPECTED_EXECUTION_SPACE "OpenMP")
elseif (Kokkos_ENABLE_SERIAL)
    set(FASTGATK_EXPECTED_EXECUTION_SPACE "Serial")
endif()
if (NOT DEFINED FASTGATK_RUNTIME_BUILD_TESTS)
    set(FASTGATK_RUNTIME_BUILD_TESTS ON)
endif()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-runtime"
                 "${CMAKE_CURRENT_BINARY_DIR}/fastgatk-runtime")
if (NOT DEFINED FASTGATK_KERNELS_BUILD_TESTS)
    set(FASTGATK_KERNELS_BUILD_TESTS ON)
endif()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels"
                 "${CMAKE_CURRENT_BINARY_DIR}/fastgatk-kernels")

set(FAST_GATK_HTSLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/htslib-build/htslib-src"
    CACHE PATH "HTSlib source/build prefix")
set(FAST_GATK_ZLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/htslib-build/install"
    CACHE PATH "zlib prefix used to build HTSlib")
set(FAST_GATK_HDF5_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/hdf5/usr"
    CACHE PATH "Optional pinned HDF5 prefix for CNV count collections")

find_path(FAST_GATK_HTS_INCLUDE_DIR htslib/hts.h
    PATHS "${FAST_GATK_HTSLIB_ROOT}" "${FAST_GATK_HTSLIB_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(FAST_GATK_HTS_LIBRARY NAMES hts libhts.a
    PATHS "${FAST_GATK_HTSLIB_ROOT}" "${FAST_GATK_HTSLIB_ROOT}/lib"
    NO_DEFAULT_PATH)
find_library(FAST_GATK_Z_LIBRARY NAMES z libz.a
    PATHS "${FAST_GATK_ZLIB_ROOT}/lib" NO_DEFAULT_PATH)
find_path(FAST_GATK_ZLIB_INCLUDE_DIR zlib.h
    PATHS "${FAST_GATK_ZLIB_ROOT}/include" "${FAST_GATK_ZLIB_ROOT}" NO_DEFAULT_PATH)
find_path(FAST_GATK_HDF5_INCLUDE_DIR hdf5.h
    PATHS "${FAST_GATK_HDF5_ROOT}/include/hdf5/serial" "${FAST_GATK_HDF5_ROOT}/include"
    NO_DEFAULT_PATH)
find_library(FAST_GATK_HDF5_LIBRARY NAMES hdf5_serial
    PATHS "${FAST_GATK_HDF5_ROOT}/lib/x86_64-linux-gnu" "${FAST_GATK_HDF5_ROOT}/lib"
    NO_DEFAULT_PATH)
find_library(FAST_GATK_SZIP_LIBRARY NAMES sz
    PATHS "${FAST_GATK_HDF5_ROOT}/lib/x86_64-linux-gnu" "${FAST_GATK_HDF5_ROOT}/lib"
    NO_DEFAULT_PATH)
if (NOT FAST_GATK_SZIP_LIBRARY)
    find_file(FAST_GATK_SZIP_LIBRARY NAMES libsz.so.2
        PATHS "${FAST_GATK_HDF5_ROOT}/lib/x86_64-linux-gnu" "${FAST_GATK_HDF5_ROOT}/lib"
        NO_DEFAULT_PATH)
endif()

if (FAST_GATK_HTS_INCLUDE_DIR AND FAST_GATK_HTS_LIBRARY AND FAST_GATK_Z_LIBRARY)
    set(FAST_GATK_HAS_HTSLIB ON)
    message(STATUS "HTSlib enabled: ${FAST_GATK_HTS_LIBRARY}")
else()
    set(FAST_GATK_HAS_HTSLIB OFF)
    message(WARNING "HTSlib not found; BAM/CRAM smoke will report BACKEND_UNAVAILABLE")
endif()

if (FAST_GATK_ZLIB_INCLUDE_DIR AND FAST_GATK_Z_LIBRARY)
    set(FAST_GATK_HAS_ZLIB ON)
else()
    set(FAST_GATK_HAS_ZLIB OFF)
endif()

if (FAST_GATK_HDF5_INCLUDE_DIR AND FAST_GATK_HDF5_LIBRARY AND FAST_GATK_SZIP_LIBRARY)
    set(FAST_GATK_HAS_HDF5 ON)
    get_filename_component(FAST_GATK_HDF5_LIBRARY_DIR "${FAST_GATK_HDF5_LIBRARY}" DIRECTORY)
    message(STATUS "HDF5 enabled: ${FAST_GATK_HDF5_LIBRARY}")
else()
    set(FAST_GATK_HAS_HDF5 OFF)
    message(STATUS "HDF5 not found; CNV HDF5 paths will report BACKEND_UNAVAILABLE")
endif()

# The normal fastgatk binaries remain entirely C++20/Kokkos/HTSlib based.
# When a pinned GATK package is available, build a deliberately isolated
# helper for reading *true* TileDB/GenomicsDB workspaces. GATK 4.6 ships its
# C++ library with the pre-C++11 libstdc++ ABI and with a libjvm dependency;
# keeping this helper in a separate process prevents either ABI from leaking
# into the Kokkos application. Users with another GATK release can point the
# cache variables at a matching libtiledbgenomicsdb.so and libjvm.so.
set(_fastgatk_gatk_package_default
    "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar")
set(FASTGATK_GENOMICSDB_PACKAGE "${_fastgatk_gatk_package_default}" CACHE FILEPATH
    "GATK package containing libtiledbgenomicsdb.so (optional)")
set(FASTGATK_GENOMICSDB_LIBRARY "" CACHE FILEPATH
    "Matching extracted libtiledbgenomicsdb.so (optional; overrides package extraction)")
set(FASTGATK_GENOMICSDB_JVM_LIBRARY "" CACHE FILEPATH
    "Matching libjvm.so required by libtiledbgenomicsdb.so (optional)")
set(_fastgatk_bridge_default OFF)
if (EXISTS "${FASTGATK_GENOMICSDB_LIBRARY}" OR EXISTS "${FASTGATK_GENOMICSDB_PACKAGE}")
    set(_fastgatk_bridge_default ON)
endif()
option(FASTGATK_ENABLE_GENOMICSDB_BRIDGE
    "Build the optional isolated true GenomicsDB workspace export helper"
    ${_fastgatk_bridge_default})

set(FASTGATK_HAS_GENOMICSDB_BRIDGE OFF)
if (FASTGATK_ENABLE_GENOMICSDB_BRIDGE)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if (NOT FASTGATK_GENOMICSDB_LIBRARY AND EXISTS "${FASTGATK_GENOMICSDB_PACKAGE}" AND Python3_Interpreter_FOUND)
        set(_fastgatk_bridge_dir "${CMAKE_CURRENT_BINARY_DIR}/genomicsdb-bridge")
        set(FASTGATK_GENOMICSDB_LIBRARY "${_fastgatk_bridge_dir}/libtiledbgenomicsdb.so")
        add_custom_command(
            OUTPUT "${FASTGATK_GENOMICSDB_LIBRARY}"
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/scripts/extract_genomicsdb_library.py"
                    --package "${FASTGATK_GENOMICSDB_PACKAGE}"
                    --output "${FASTGATK_GENOMICSDB_LIBRARY}"
            DEPENDS "${FASTGATK_GENOMICSDB_PACKAGE}"
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/extract_genomicsdb_library.py"
            COMMENT "Extracting libtiledbgenomicsdb.so from the GATK package"
            VERBATIM)
        add_custom_target(fastgatk-genomicsdb-library
            DEPENDS "${FASTGATK_GENOMICSDB_LIBRARY}")
    endif()
    if (NOT FASTGATK_GENOMICSDB_JVM_LIBRARY AND EXISTS
        "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/jdk17/lib/server/libjvm.so")
        set(FASTGATK_GENOMICSDB_JVM_LIBRARY
            "${CMAKE_CURRENT_SOURCE_DIR}/../third_party/jdk17/lib/server/libjvm.so")
    endif()
    if (FASTGATK_GENOMICSDB_LIBRARY AND FASTGATK_GENOMICSDB_JVM_LIBRARY AND
        (EXISTS "${FASTGATK_GENOMICSDB_LIBRARY}" OR TARGET fastgatk-genomicsdb-library))
        add_executable(fastgatk-genomicsdb-export src/genomicsdb_export_bridge.cpp)
        target_compile_features(fastgatk-genomicsdb-export PRIVATE cxx_std_20)
        target_compile_definitions(fastgatk-genomicsdb-export PRIVATE _GLIBCXX_USE_CXX11_ABI=0)
        target_compile_options(fastgatk-genomicsdb-export PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-O2;-Wall;-Wextra;-Wpedantic>)
        target_link_options(fastgatk-genomicsdb-export PRIVATE "LINKER:--no-as-needed")
        target_link_libraries(fastgatk-genomicsdb-export PRIVATE
            "${FASTGATK_GENOMICSDB_LIBRARY}" "${FASTGATK_GENOMICSDB_JVM_LIBRARY}" dl pthread)
        get_filename_component(_fastgatk_bridge_library_dir
            "${FASTGATK_GENOMICSDB_LIBRARY}" DIRECTORY)
        get_filename_component(_fastgatk_bridge_jvm_dir
            "${FASTGATK_GENOMICSDB_JVM_LIBRARY}" DIRECTORY)
        set_target_properties(fastgatk-genomicsdb-export PROPERTIES
            BUILD_RPATH "${_fastgatk_bridge_library_dir};${_fastgatk_bridge_jvm_dir}")
        if (TARGET fastgatk-genomicsdb-library)
            add_dependencies(fastgatk-genomicsdb-export fastgatk-genomicsdb-library)
        endif()
        set(FASTGATK_HAS_GENOMICSDB_BRIDGE ON)
        message(STATUS "True GenomicsDB workspace bridge enabled")
    else()
        message(STATUS "True GenomicsDB workspace bridge unavailable; native sparse index remains the fallback")
    endif()
endif()

add_library(fastgatk-hdf5 STATIC src/hdf5_count_collection.cpp)
target_include_directories(fastgatk-hdf5 PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include")
if (FAST_GATK_HAS_HDF5)
    target_compile_definitions(fastgatk-hdf5 PUBLIC FASTGATK_HAS_HDF5=1)
    target_include_directories(fastgatk-hdf5 PRIVATE "${FAST_GATK_HDF5_INCLUDE_DIR}")
    target_link_libraries(fastgatk-hdf5 PUBLIC "${FAST_GATK_HDF5_LIBRARY}" "${FAST_GATK_SZIP_LIBRARY}" m dl)
    target_link_options(fastgatk-hdf5 INTERFACE "LINKER:--no-as-needed")
    set_target_properties(fastgatk-hdf5 PROPERTIES BUILD_RPATH "${FAST_GATK_HDF5_LIBRARY_DIR}")
else()
    target_compile_definitions(fastgatk-hdf5 PUBLIC FASTGATK_HAS_HDF5=0)
endif()
target_compile_features(fastgatk-hdf5 PUBLIC cxx_std_20)
target_compile_options(fastgatk-hdf5 PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

add_library(fastgatk-hts STATIC ../fastgatk-core/src/hts_reader.cpp)
target_include_directories(fastgatk-hts PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
if (FAST_GATK_HAS_HTSLIB)
    target_compile_definitions(fastgatk-hts PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-hts PRIVATE "${FAST_GATK_HTS_INCLUDE_DIR}")
    target_link_libraries(fastgatk-hts PUBLIC "${FAST_GATK_HTS_LIBRARY}" "${FAST_GATK_Z_LIBRARY}" m pthread)
else()
    target_compile_definitions(fastgatk-hts PRIVATE FASTGATK_HAS_HTSLIB=0)
endif()
target_compile_features(fastgatk-hts PUBLIC cxx_std_20)
target_compile_options(fastgatk-hts PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

# Shared HTSJDK-compatible F1R2 archive serialization.  Mutect2 and the
# standalone CollectF1R2Counts tool use the same Host tar/metrics contract;
# numerical pileup remains in their respective Kokkos paths.
if (FAST_GATK_HAS_ZLIB)
    add_library(fastgatk-f1r2 STATIC src/f1r2_archive.cpp)
    target_include_directories(fastgatk-f1r2 PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    target_link_libraries(fastgatk-f1r2 PUBLIC "${FAST_GATK_Z_LIBRARY}")
    target_compile_features(fastgatk-f1r2 PUBLIC cxx_std_20)
    target_compile_options(fastgatk-f1r2 PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-hts-reader-test tests/hts_reader_projection_test.cpp)
    target_include_directories(fastgatk-hts-reader-test PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-hts-reader-test PRIVATE fastgatk-hts fastgatk-core)
    target_compile_features(fastgatk-hts-reader-test PRIVATE cxx_std_20)
    target_compile_options(fastgatk-hts-reader-test PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
    add_test(NAME fastgatk-hts-reader-indel-quality
        COMMAND fastgatk-hts-reader-test
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/flow_reader_fixture.sam")
    set(FASTGATK_HCC1143_INDEXED_BAM
        "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam")
    if (EXISTS "${FASTGATK_HCC1143_INDEXED_BAM}" AND
        EXISTS "${FASTGATK_HCC1143_INDEXED_BAM}.bai")
        # Streaming first opens an unbounded metadata reader.  It must still
        # detect the sidecar index before per-tile readers receive -L spans.
        add_test(NAME fastgatk-hts-reader-unbounded-index-capability
            COMMAND fastgatk-hts-reader-test "${FASTGATK_HCC1143_INDEXED_BAM}"
                    --require-index)
    endif()
    set(FASTGATK_FLOW_ORACLE_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/../gatk-source/src/test/resources/org/broadinstitute/hellbender/utils/read/flow/reads/outputs")
    set(FASTGATK_FLOW_ORACLE_BAM
        "${CMAKE_CURRENT_SOURCE_DIR}/../gatk-source/src/test/resources/org/broadinstitute/hellbender/utils/read/flow/reads/input/sample.bam")
    set(FASTGATK_FLOW_ORACLE_T0_BAM
        "${CMAKE_CURRENT_SOURCE_DIR}/../gatk-source/src/test/resources/org/broadinstitute/hellbender/utils/read/flow/reads/input/sample.t0.bam")
    if (EXISTS "${FASTGATK_FLOW_ORACLE_DIR}/sample.0.key.txt" AND
        EXISTS "${FASTGATK_FLOW_ORACLE_BAM}" AND
        EXISTS "${FASTGATK_FLOW_ORACLE_T0_BAM}")
        add_test(NAME fastgatk-flow-codec-gatk-oracle
            COMMAND fastgatk-hts-reader-test "${FASTGATK_FLOW_ORACLE_BAM}"
                    --flow-oracle "${FASTGATK_FLOW_ORACLE_DIR}" flow-order=GTAC)
        add_test(NAME fastgatk-flow-codec-gatk-t0-oracle
            COMMAND fastgatk-hts-reader-test "${FASTGATK_FLOW_ORACLE_T0_BAM}"
                    --flow-oracle "${FASTGATK_FLOW_ORACLE_DIR}" t0 flow-order=TGCA)
    endif()
endif()

add_executable(fastgatk-hc-smoke src/hc_smoke.cpp)
target_include_directories(fastgatk-hc-smoke PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-hc-smoke PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-hc-smoke PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
else()
    target_compile_definitions(fastgatk-hc-smoke PRIVATE FASTGATK_HAS_ZLIB=0)
endif()
target_link_libraries(fastgatk-hc-smoke PRIVATE fastgatk-kernels fastgatk-runtime fastgatk-hts)
target_compile_options(fastgatk-hc-smoke PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_library(fastgatk-calling STATIC
    src/calling_pipeline.cpp
    src/flow_hmer.cpp)
target_include_directories(fastgatk-calling PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-calling PUBLIC Kokkos::kokkos fastgatk-core fastgatk-kernels)
target_compile_features(fastgatk-calling PUBLIC cxx_std_20)
target_compile_options(fastgatk-calling PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-calling-reference-boundary-smoke
    tests/calling_reference_boundary_smoke.cpp)
target_include_directories(fastgatk-calling-reference-boundary-smoke PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-calling-reference-boundary-smoke PRIVATE
    fastgatk-calling fastgatk-kernels fastgatk-runtime fastgatk-hts)
target_compile_options(fastgatk-calling-reference-boundary-smoke PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
add_test(NAME fastgatk-calling-reference-boundary-smoke
    COMMAND fastgatk-calling-reference-boundary-smoke)

add_executable(fastgatk-haplotype-cigar-trim-smoke
    tests/haplotype_cigar_trim_smoke.cpp)
target_include_directories(fastgatk-haplotype-cigar-trim-smoke PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-haplotype-cigar-trim-smoke PRIVATE
    fastgatk-calling fastgatk-kernels fastgatk-runtime fastgatk-hts)
target_compile_options(fastgatk-haplotype-cigar-trim-smoke PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
add_test(NAME fastgatk-haplotype-cigar-trim-smoke
    COMMAND fastgatk-haplotype-cigar-trim-smoke)

add_executable(fastgatk-hc-call src/hc_call.cpp)
target_include_directories(fastgatk-hc-call PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-hc-call PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-hc-call PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
else()
    target_compile_definitions(fastgatk-hc-call PRIVATE FASTGATK_HAS_ZLIB=0)
endif()
if (FAST_GATK_HAS_HTSLIB)
    target_compile_definitions(fastgatk-hc-call PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-hc-call PRIVATE "${FAST_GATK_HTS_INCLUDE_DIR}")
else()
    target_compile_definitions(fastgatk-hc-call PRIVATE FASTGATK_HAS_HTSLIB=0)
endif()
target_link_libraries(fastgatk-hc-call PRIVATE fastgatk-kernels fastgatk-calling fastgatk-runtime fastgatk-hts)
target_compile_options(fastgatk-hc-call PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-calling-indel-smoke tests/calling_indel_smoke.cpp)
    target_include_directories(fastgatk-calling-indel-smoke PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-calling-indel-smoke PRIVATE
        fastgatk-calling fastgatk-kernels fastgatk-runtime fastgatk-hts)
    target_compile_options(fastgatk-calling-indel-smoke PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-calling-indel-smoke COMMAND fastgatk-calling-indel-smoke)
    add_test(NAME fastgatk-hc-reverse-strand-indel-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_reverse_strand_indel.py")
    set_tests_properties(fastgatk-hc-reverse-strand-indel-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CALLING_INDEL_SMOKE_BINARY=$<TARGET_FILE:fastgatk-calling-indel-smoke>")
endif()

foreach(bqsr_target IN ITEMS fastgatk-bqsr fastgatk-apply-bqsr)
    add_executable(${bqsr_target} src/bqsr_tool.cpp)
    target_include_directories(${bqsr_target} PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(${bqsr_target} PRIVATE fastgatk-kernels fastgatk-runtime fastgatk-hts)
    target_compile_features(${bqsr_target} PRIVATE cxx_std_20)
    target_compile_options(${bqsr_target} PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    if (FAST_GATK_HAS_HTSLIB)
        target_compile_definitions(${bqsr_target} PRIVATE FASTGATK_HAS_HTSLIB=1)
        target_include_directories(${bqsr_target} PRIVATE "${FAST_GATK_HTS_INCLUDE_DIR}")
    else()
        target_compile_definitions(${bqsr_target} PRIVATE FASTGATK_HAS_HTSLIB=0)
    endif()
endforeach()
target_compile_definitions(fastgatk-apply-bqsr PRIVATE FASTGATK_BQSR_APPLY=1)

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-genotype-gvcf src/genotype_gvcf_tool.cpp)
    target_compile_definitions(fastgatk-genotype-gvcf PRIVATE FASTGATK_HAS_HTSLIB=1)
    if (TARGET fastgatk-genomicsdb-export)
        # The path is embedded only as a default. FASTGATK_GENOMICSDB_BRIDGE
        # remains an override for deployments that ship a matching GATK
        # library outside this build tree.
        target_compile_definitions(fastgatk-genotype-gvcf PRIVATE
            FASTGATK_GENOMICSDB_BRIDGE_DEFAULT="$<TARGET_FILE:fastgatk-genomicsdb-export>")
        add_dependencies(fastgatk-genotype-gvcf fastgatk-genomicsdb-export)
    endif()
    target_include_directories(fastgatk-genotype-gvcf PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-genotype-gvcf PRIVATE fastgatk-runtime fastgatk-hts fastgatk-kernels)
    target_compile_features(fastgatk-genotype-gvcf PRIVATE cxx_std_20)
    target_compile_options(fastgatk-genotype-gvcf PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-reblock-gvcf src/reblock_gvcf_tool.cpp)
    target_compile_definitions(fastgatk-reblock-gvcf PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-reblock-gvcf PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-reblock-gvcf PRIVATE fastgatk-runtime fastgatk-hts fastgatk-kernels fastgatk-core)
    target_compile_features(fastgatk-reblock-gvcf PRIVATE cxx_std_20)
    target_compile_options(fastgatk-reblock-gvcf PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-select-variants src/select_variants_tool.cpp)
    target_compile_definitions(fastgatk-select-variants PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-select-variants PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include")
    target_link_libraries(fastgatk-select-variants PRIVATE fastgatk-runtime fastgatk-hts fastgatk-kernels)
    target_compile_features(fastgatk-select-variants PRIVATE cxx_std_20)
    target_compile_options(fastgatk-select-variants PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-gather-vcfs src/gather_vcfs_tool.cpp)
    target_compile_definitions(fastgatk-gather-vcfs PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-gather-vcfs PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-gather-vcfs PRIVATE fastgatk-runtime fastgatk-hts fastgatk-core)
    target_compile_features(fastgatk-gather-vcfs PRIVATE cxx_std_20)
    target_compile_options(fastgatk-gather-vcfs PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-left-align-trim src/left_align_tool.cpp)
    target_compile_definitions(fastgatk-left-align-trim PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-left-align-trim PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include")
    target_link_libraries(fastgatk-left-align-trim PRIVATE fastgatk-runtime fastgatk-hts fastgatk-kernels)
    target_compile_features(fastgatk-left-align-trim PRIVATE cxx_std_20)
    target_compile_options(fastgatk-left-align-trim PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-variant-filtration src/variant_filtration_tool.cpp)
    target_compile_definitions(fastgatk-variant-filtration PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-variant-filtration PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-variant-filtration PRIVATE fastgatk-runtime fastgatk-hts fastgatk-core)
    target_compile_features(fastgatk-variant-filtration PRIVATE cxx_std_20)
    target_compile_options(fastgatk-variant-filtration PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-apply-vqsr src/apply_vqsr_tool.cpp)
    target_compile_definitions(fastgatk-apply-vqsr PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-apply-vqsr PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-apply-vqsr PRIVATE fastgatk-runtime fastgatk-hts fastgatk-core)
    target_compile_features(fastgatk-apply-vqsr PRIVATE cxx_std_20)
    target_compile_options(fastgatk-apply-vqsr PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-variant-recalibrator src/variant_recalibrator_tool.cpp)
    target_compile_definitions(fastgatk-variant-recalibrator PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-variant-recalibrator PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    if (FAST_GATK_HAS_ZLIB)
        target_compile_definitions(fastgatk-variant-recalibrator PRIVATE FASTGATK_HAS_ZLIB=1)
        target_include_directories(fastgatk-variant-recalibrator PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    else()
        target_compile_definitions(fastgatk-variant-recalibrator PRIVATE FASTGATK_HAS_ZLIB=0)
    endif()
    target_link_libraries(fastgatk-variant-recalibrator PRIVATE fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    if (FAST_GATK_HAS_ZLIB)
        target_link_libraries(fastgatk-variant-recalibrator PRIVATE "${FAST_GATK_Z_LIBRARY}")
    endif()
    target_compile_features(fastgatk-variant-recalibrator PRIVATE cxx_std_20)
    target_compile_options(fastgatk-variant-recalibrator PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-sort-sam src/sort_sam_tool.cpp)
    target_compile_definitions(fastgatk-sort-sam PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-sort-sam PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}")
    target_link_libraries(fastgatk-sort-sam PRIVATE fastgatk-runtime fastgatk-hts)
    target_compile_features(fastgatk-sort-sam PRIVATE cxx_std_20)
    target_compile_options(fastgatk-sort-sam PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-mark-duplicates src/mark_duplicates_tool.cpp)
    target_compile_definitions(fastgatk-mark-duplicates PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-mark-duplicates PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}")
    target_link_libraries(fastgatk-mark-duplicates PRIVATE fastgatk-runtime fastgatk-hts)
    target_compile_features(fastgatk-mark-duplicates PRIVATE cxx_std_20)
    target_compile_options(fastgatk-mark-duplicates PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-mutect2 src/mutect2_tool.cpp)
    target_compile_definitions(fastgatk-mutect2 PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-mutect2 PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
        "${FAST_GATK_HTS_INCLUDE_DIR}")
    target_link_libraries(fastgatk-mutect2 PRIVATE
        fastgatk-calling fastgatk-kernels fastgatk-runtime fastgatk-hts
        $<$<BOOL:${FAST_GATK_HAS_ZLIB}>:fastgatk-f1r2>)
    target_compile_features(fastgatk-mutect2 PRIVATE cxx_std_20)
    target_compile_options(fastgatk-mutect2 PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-combine-gvcfs src/combine_gvcf_tool.cpp)
    target_compile_definitions(fastgatk-combine-gvcfs PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-combine-gvcfs PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-combine-gvcfs PRIVATE fastgatk-kernels fastgatk-runtime fastgatk-hts)
    target_compile_features(fastgatk-combine-gvcfs PRIVATE cxx_std_20)
    target_compile_options(fastgatk-combine-gvcfs PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-filter-mutect-calls src/filter_mutect_tool.cpp)
    target_compile_definitions(fastgatk-filter-mutect-calls PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-filter-mutect-calls PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include"
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    target_link_libraries(fastgatk-filter-mutect-calls PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-kernels fastgatk-core
        Kokkos::kokkos "${FAST_GATK_Z_LIBRARY}" m pthread)
    target_compile_features(fastgatk-filter-mutect-calls PRIVATE cxx_std_20)
    target_compile_options(fastgatk-filter-mutect-calls PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

# Full Java-parity learning-trajectory oracle (no GATK jar needed): replays
# the three-pass SomaticClusteringModel learning over the contamination-joint
# ten-ALT fixture via fmc_java_model.hpp + fastgatk-kernels and pins the final
# model state and threshold (R17: 5 EM(false) iterations per accepted split).
add_executable(fastgatk-fmc-trajectory-oracle tests/fmc_trajectory_oracle.cpp)
target_include_directories(fastgatk-fmc-trajectory-oracle PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-kernels/include")
target_link_libraries(fastgatk-fmc-trajectory-oracle PRIVATE
    fastgatk-kernels fastgatk-core Kokkos::kokkos)
target_compile_features(fastgatk-fmc-trajectory-oracle PRIVATE cxx_std_20)
target_compile_options(fastgatk-fmc-trajectory-oracle PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
add_test(NAME fastgatk-fmc-trajectory-oracle
    COMMAND fastgatk-fmc-trajectory-oracle)
# Determinism pins: the trajectory must be bit-identical regardless of the
# OpenMP worker count (single work item / independent RangePolicy datums).
foreach(tcount 1 4 16)
    add_test(NAME fastgatk-fmc-trajectory-oracle-omp${tcount}
        COMMAND fastgatk-fmc-trajectory-oracle)
    set_tests_properties(fastgatk-fmc-trajectory-oracle-omp${tcount} PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=${tcount}")
endforeach()

add_executable(fastgatk-genomicsdb-import src/genomicsdb_import_tool.cpp)
if (FAST_GATK_HAS_HTSLIB)
    target_compile_definitions(fastgatk-genomicsdb-import PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-genomicsdb-import PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include"
        "${FAST_GATK_HTS_INCLUDE_DIR}")
    target_link_libraries(fastgatk-genomicsdb-import PRIVATE fastgatk-runtime fastgatk-hts
        "${FAST_GATK_HTS_LIBRARY}" "${FAST_GATK_Z_LIBRARY}" m pthread)
else()
    target_compile_definitions(fastgatk-genomicsdb-import PRIVATE FASTGATK_HAS_HTSLIB=0)
    target_link_libraries(fastgatk-genomicsdb-import PRIVATE fastgatk-runtime)
endif()
target_compile_features(fastgatk-genomicsdb-import PRIVATE cxx_std_20)
target_compile_options(fastgatk-genomicsdb-import PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-gather-bqsr-reports src/gather_bqsr_tool.cpp)
target_link_libraries(fastgatk-gather-bqsr-reports PRIVATE fastgatk-runtime)
target_compile_features(fastgatk-gather-bqsr-reports PRIVATE cxx_std_20)
target_compile_options(fastgatk-gather-bqsr-reports PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-analyze-covariates src/analyze_covariates_tool.cpp)
if (FASTGATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-analyze-covariates PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-analyze-covariates PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
else()
    target_compile_definitions(fastgatk-analyze-covariates PRIVATE FASTGATK_HAS_ZLIB=0)
endif()
target_include_directories(fastgatk-analyze-covariates PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-analyze-covariates PRIVATE fastgatk-runtime fastgatk-core Kokkos::kokkos
    $<$<BOOL:${FAST_GATK_HAS_ZLIB}>:${FAST_GATK_Z_LIBRARY}>)
target_compile_features(fastgatk-analyze-covariates PRIVATE cxx_std_20)
target_compile_options(fastgatk-analyze-covariates PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-variants-to-table src/variants_to_table_tool.cpp)
    target_compile_definitions(fastgatk-variants-to-table PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-variants-to-table PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-variants-to-table PRIVATE fastgatk-runtime fastgatk-hts)
    target_compile_features(fastgatk-variants-to-table PRIVATE cxx_std_20)
    target_compile_options(fastgatk-variants-to-table PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-variant-eval src/variant_eval_tool.cpp)
    target_compile_definitions(fastgatk-variant-eval PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-variant-eval PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-variant-eval PRIVATE fastgatk-runtime fastgatk-hts Kokkos::kokkos)
    target_compile_features(fastgatk-variant-eval PRIVATE cxx_std_20)
    target_compile_options(fastgatk-variant-eval PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-validate-variants src/validate_variants_tool.cpp)
    target_compile_definitions(fastgatk-validate-variants PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-validate-variants PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(fastgatk-validate-variants PRIVATE fastgatk-runtime fastgatk-hts)
    target_compile_features(fastgatk-validate-variants PRIVATE cxx_std_20)
    target_compile_options(fastgatk-validate-variants PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-get-pileup-summaries src/get_pileup_summaries_tool.cpp)
    target_compile_definitions(fastgatk-get-pileup-summaries PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-get-pileup-summaries PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-get-pileup-summaries PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-kernels Kokkos::kokkos)
    target_compile_features(fastgatk-get-pileup-summaries PRIVATE cxx_std_20)
    target_compile_options(fastgatk-get-pileup-summaries PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

add_executable(fastgatk-calculate-contamination
    src/calculate_contamination_tool.cpp
    src/contamination_kernel_segmenter.cpp)
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-calculate-contamination PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-calculate-contamination PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
else()
    target_compile_definitions(fastgatk-calculate-contamination PRIVATE FASTGATK_HAS_ZLIB=0)
endif()
target_include_directories(fastgatk-calculate-contamination PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-calculate-contamination PRIVATE fastgatk-runtime fastgatk-core Kokkos::kokkos)
if (FAST_GATK_HAS_ZLIB)
    target_link_libraries(fastgatk-calculate-contamination PRIVATE "${FAST_GATK_Z_LIBRARY}")
endif()
target_compile_features(fastgatk-calculate-contamination PRIVATE cxx_std_20)
target_compile_options(fastgatk-calculate-contamination PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-contamination-kernel-segmenter-contract
    tests/contamination_kernel_segmenter_contract.cpp
    src/contamination_kernel_segmenter.cpp)
target_include_directories(fastgatk-contamination-kernel-segmenter-contract PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-contamination-kernel-segmenter-contract PRIVATE Kokkos::kokkos)
target_compile_features(fastgatk-contamination-kernel-segmenter-contract PRIVATE cxx_std_20)
target_compile_options(fastgatk-contamination-kernel-segmenter-contract PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
add_test(NAME fastgatk-contamination-kernel-segmenter-contract
    COMMAND fastgatk-contamination-kernel-segmenter-contract)

add_executable(fastgatk-gather-pileup-summaries src/gather_pileup_summaries_tool.cpp)
target_include_directories(fastgatk-gather-pileup-summaries PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-gather-pileup-summaries PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-gather-pileup-summaries PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    target_link_libraries(fastgatk-gather-pileup-summaries PRIVATE fastgatk-runtime "${FAST_GATK_Z_LIBRARY}")
else()
    target_compile_definitions(fastgatk-gather-pileup-summaries PRIVATE FASTGATK_HAS_ZLIB=0)
    target_link_libraries(fastgatk-gather-pileup-summaries PRIVATE fastgatk-runtime)
endif()
target_compile_features(fastgatk-gather-pileup-summaries PRIVATE cxx_std_20)
target_compile_options(fastgatk-gather-pileup-summaries PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-learn-read-orientation-model src/learn_read_orientation_model_tool.cpp)
target_include_directories(fastgatk-learn-read-orientation-model PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-learn-read-orientation-model PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-learn-read-orientation-model PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    target_link_libraries(fastgatk-learn-read-orientation-model PRIVATE fastgatk-runtime fastgatk-core Kokkos::kokkos "${FAST_GATK_Z_LIBRARY}")
else()
    target_compile_definitions(fastgatk-learn-read-orientation-model PRIVATE FASTGATK_HAS_ZLIB=0)
    target_link_libraries(fastgatk-learn-read-orientation-model PRIVATE fastgatk-runtime fastgatk-core Kokkos::kokkos)
endif()
target_compile_features(fastgatk-learn-read-orientation-model PRIVATE cxx_std_20)
target_compile_options(fastgatk-learn-read-orientation-model PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-denoise-read-counts src/denoise_read_counts_tool.cpp)
target_include_directories(fastgatk-denoise-read-counts PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-denoise-read-counts PRIVATE Kokkos::kokkos fastgatk-core fastgatk-hdf5)
target_compile_features(fastgatk-denoise-read-counts PRIVATE cxx_std_20)
target_compile_options(fastgatk-denoise-read-counts PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-create-read-count-panel-of-normals
    src/create_read_count_panel_of_normals_tool.cpp)
target_include_directories(fastgatk-create-read-count-panel-of-normals PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-create-read-count-panel-of-normals PRIVATE
    Kokkos::kokkos fastgatk-core fastgatk-hdf5)
target_compile_features(fastgatk-create-read-count-panel-of-normals PRIVATE cxx_std_20)
target_compile_options(fastgatk-create-read-count-panel-of-normals PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-call-copy-ratio-segments src/call_copy_ratio_segments_tool.cpp)
target_include_directories(fastgatk-call-copy-ratio-segments PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-call-copy-ratio-segments PRIVATE Kokkos::kokkos fastgatk-core)
target_compile_features(fastgatk-call-copy-ratio-segments PRIVATE cxx_std_20)
target_compile_options(fastgatk-call-copy-ratio-segments PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-model-segments src/model_segments_tool.cpp)
target_include_directories(fastgatk-model-segments PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
target_link_libraries(fastgatk-model-segments PRIVATE Kokkos::kokkos fastgatk-core)
target_compile_features(fastgatk-model-segments PRIVATE cxx_std_20)
target_compile_options(fastgatk-model-segments PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

add_executable(fastgatk-gather-tranches src/gather_tranches_tool.cpp)
if (FAST_GATK_HAS_ZLIB)
    target_compile_definitions(fastgatk-gather-tranches PRIVATE FASTGATK_HAS_ZLIB=1)
    target_include_directories(fastgatk-gather-tranches PRIVATE "${FAST_GATK_ZLIB_INCLUDE_DIR}")
    target_link_libraries(fastgatk-gather-tranches PRIVATE "${FAST_GATK_Z_LIBRARY}")
else()
    target_compile_definitions(fastgatk-gather-tranches PRIVATE FASTGATK_HAS_ZLIB=0)
endif()
target_compile_features(fastgatk-gather-tranches PRIVATE cxx_std_20)
target_compile_options(fastgatk-gather-tranches PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
target_link_libraries(fastgatk-gather-tranches PRIVATE fastgatk-runtime)

if (FAST_GATK_HAS_HTSLIB)
    add_executable(fastgatk-annotate-intervals src/annotate_intervals_tool.cpp)
    target_compile_definitions(fastgatk-annotate-intervals PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-annotate-intervals PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-annotate-intervals PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-annotate-intervals PRIVATE cxx_std_20)
    target_compile_options(fastgatk-annotate-intervals PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-count-bases-in-reference src/count_bases_in_reference_tool.cpp)
    target_compile_definitions(fastgatk-count-bases-in-reference PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-count-bases-in-reference PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-count-bases-in-reference PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-count-bases-in-reference PRIVATE cxx_std_20)
    target_compile_options(fastgatk-count-bases-in-reference PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-compare-references src/compare_references_tool.cpp)
    target_compile_definitions(fastgatk-compare-references PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-compare-references PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-compare-references PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-compare-references PRIVATE cxx_std_20)
    target_compile_options(fastgatk-compare-references PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-check-reference-compatibility src/check_reference_compatibility_tool.cpp)
    target_compile_definitions(fastgatk-check-reference-compatibility PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-check-reference-compatibility PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-check-reference-compatibility PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core)
    target_compile_features(fastgatk-check-reference-compatibility PRIVATE cxx_std_20)
    target_compile_options(fastgatk-check-reference-compatibility PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-fasta-reference-maker src/fasta_reference_tool.cpp)
    target_compile_definitions(fastgatk-fasta-reference-maker PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-fasta-reference-maker PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-fasta-reference-maker PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-fasta-reference-maker PRIVATE cxx_std_20)
    target_compile_options(fastgatk-fasta-reference-maker PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-fasta-alternate-reference-maker src/fasta_reference_tool.cpp)
    target_compile_definitions(fastgatk-fasta-alternate-reference-maker PRIVATE
        FAST_GATK_HAS_HTSLIB=1 FASTGATK_FASTA_ALTERNATE=1)
    target_include_directories(fastgatk-fasta-alternate-reference-maker PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-fasta-alternate-reference-maker PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-fasta-alternate-reference-maker PRIVATE cxx_std_20)
    target_compile_options(fastgatk-fasta-alternate-reference-maker PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-shift-fasta src/shift_fasta_tool.cpp)
    target_compile_definitions(fastgatk-shift-fasta PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-shift-fasta PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-shift-fasta PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-shift-fasta PRIVATE cxx_std_20)
    target_compile_options(fastgatk-shift-fasta PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-index-feature-file src/index_feature_file_tool.cpp)
    target_compile_definitions(fastgatk-index-feature-file PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-index-feature-file PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-index-feature-file PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-index-feature-file PRIVATE cxx_std_20)
    target_compile_options(fastgatk-index-feature-file PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-count-reads src/read_metrics_tool.cpp)
    target_compile_definitions(fastgatk-count-reads PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-count-reads PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-count-reads PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-count-reads PRIVATE cxx_std_20)
    target_compile_options(fastgatk-count-reads PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-flag-stat src/read_metrics_tool.cpp)
    target_compile_definitions(fastgatk-flag-stat PRIVATE FASTGATK_HAS_HTSLIB=1 FASTGATK_FLAGSTAT=1)
    target_include_directories(fastgatk-flag-stat PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-flag-stat PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-flag-stat PRIVATE cxx_std_20)
    target_compile_options(fastgatk-flag-stat PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-split-intervals src/split_intervals_tool.cpp)
    target_compile_definitions(fastgatk-split-intervals PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-split-intervals PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-split-intervals PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-split-intervals PRIVATE cxx_std_20)
    target_compile_options(fastgatk-split-intervals PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-filter-intervals src/filter_intervals_tool.cpp)
    target_compile_definitions(fastgatk-filter-intervals PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-filter-intervals PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-filter-intervals PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-hdf5 fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-filter-intervals PRIVATE cxx_std_20)
    target_compile_options(fastgatk-filter-intervals PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-preprocess-intervals src/preprocess_intervals_tool.cpp)
    target_compile_definitions(fastgatk-preprocess-intervals PRIVATE FAST_GATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-preprocess-intervals PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-preprocess-intervals PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-preprocess-intervals PRIVATE cxx_std_20)
    target_compile_options(fastgatk-preprocess-intervals PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-collect-read-counts src/collect_read_counts_tool.cpp)
    target_compile_definitions(fastgatk-collect-read-counts PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-collect-read-counts PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-collect-read-counts PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core fastgatk-hdf5 Kokkos::kokkos)
    target_compile_features(fastgatk-collect-read-counts PRIVATE cxx_std_20)
    target_compile_options(fastgatk-collect-read-counts PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-collect-f1r2-counts src/collect_f1r2_counts_tool.cpp)
    target_compile_definitions(fastgatk-collect-f1r2-counts PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-collect-f1r2-counts PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${FAST_GATK_ZLIB_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-collect-f1r2-counts PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos "${FAST_GATK_HTS_LIBRARY}"
        "${FAST_GATK_Z_LIBRARY}" m pthread)
    target_compile_features(fastgatk-collect-f1r2-counts PRIVATE cxx_std_20)
    target_compile_options(fastgatk-collect-f1r2-counts PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-collect-allelic-counts src/collect_allelic_counts_tool.cpp)
    target_compile_definitions(fastgatk-collect-allelic-counts PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-collect-allelic-counts PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-collect-allelic-counts PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-collect-allelic-counts PRIVATE cxx_std_20)
    target_compile_options(fastgatk-collect-allelic-counts PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

    add_executable(fastgatk-depth-of-coverage src/depth_of_coverage_tool.cpp)
    target_compile_definitions(fastgatk-depth-of-coverage PRIVATE FASTGATK_HAS_HTSLIB=1)
    target_include_directories(fastgatk-depth-of-coverage PRIVATE
        "${FAST_GATK_HTS_INCLUDE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-core/include")
    target_link_libraries(fastgatk-depth-of-coverage PRIVATE
        fastgatk-runtime fastgatk-hts fastgatk-core Kokkos::kokkos)
    target_compile_features(fastgatk-depth-of-coverage PRIVATE cxx_std_20)
    target_compile_options(fastgatk-depth-of-coverage PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
    add_test(NAME fastgatk-launcher-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_launcher_contract.py")
    add_test(NAME fastgatk-progress-score-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/compute_progress_score.py")
    add_test(NAME fastgatk-fixture-digest-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_fixture_digests.py")
    add_test(NAME fastgatk-kokkos-api-boundary
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_kokkos_api_boundary.py")
    # Configuration-only gate: validates the selected execution space,
    # architecture flags, generated Kokkos package metadata and compiler
    # provenance.  It deliberately does not initialize a device, so CPU-only
    # CI can catch backend-selection/fallback mistakes before GPU runners are
    # available.
    add_test(NAME fastgatk-kokkos-backend-config
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_kokkos_backend_config.py"
                --cache "${CMAKE_BINARY_DIR}/CMakeCache.txt"
                --kokkos-config "${CMAKE_CURRENT_BINARY_DIR}/kokkos/KokkosConfigCommon.cmake")
endif()
if (FAST_GATK_HAS_HTSLIB AND Python3_Interpreter_FOUND)
    add_test(NAME fastgatk-annotate-intervals-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_annotate_intervals.py")
    add_test(NAME fastgatk-count-bases-in-reference-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_count_bases_in_reference.py")
    add_test(NAME fastgatk-compare-references-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_compare_references.py")
    add_test(NAME fastgatk-check-reference-compatibility-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_check_reference_compatibility.py")
    add_test(NAME fastgatk-fasta-reference-maker-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_fasta_reference_tools.py")
    add_test(NAME fastgatk-fasta-alternate-iupac-hom-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_fasta_alternate_iupac_hom_gatk_oracle.py")
    add_test(NAME fastgatk-shift-fasta-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_shift_fasta.py")
    add_test(NAME fastgatk-index-feature-file-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_index_feature_file.py")
    set_tests_properties(fastgatk-index-feature-file-contract PROPERTIES
        ENVIRONMENT "FASTGATK_INDEX_FEATURE_FILE_BINARY=$<TARGET_FILE:fastgatk-index-feature-file>")
    add_test(NAME fastgatk-count-reads-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_count_reads.py")
    set_tests_properties(fastgatk-count-reads-contract PROPERTIES
        ENVIRONMENT "FASTGATK_COUNT_READS_BINARY=$<TARGET_FILE:fastgatk-count-reads>")
    add_test(NAME fastgatk-flag-stat-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_flag_stat.py")
    set_tests_properties(fastgatk-flag-stat-contract PROPERTIES
        ENVIRONMENT "FASTGATK_FLAG_STAT_BINARY=$<TARGET_FILE:fastgatk-flag-stat>")
    add_test(NAME fastgatk-split-intervals-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_split_intervals.py")
    set_tests_properties(fastgatk-split-intervals-contract PROPERTIES
        ENVIRONMENT "FASTGATK_SPLIT_INTERVALS_BINARY=$<TARGET_FILE:fastgatk-split-intervals>")
    add_test(NAME fastgatk-filter-intervals-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_intervals.py")
    add_test(NAME fastgatk-preprocess-intervals-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_preprocess_intervals.py")
    add_test(NAME fastgatk-remote-staging-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_remote_staging.py")
    add_test(NAME fastgatk-resource-limits-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_resource_limits.py")
    set_tests_properties(fastgatk-resource-limits-contract PROPERTIES
        ENVIRONMENT "FASTGATK_RESOURCE_LIMITS_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-slurm-resource-wrapper-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_slurm_wrapper.py")
    add_test(NAME fastgatk-scheduler-retry-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_scheduler_retry_contract.py")
    # End-to-end local fallback for the production Nextflow/SLURM boundary.
    # The script uses bundled Nextflow/JDK when available and otherwise still
    # exercises the same dispatcher, wrapper, scatter/gather, and manifest
    # contracts locally; it never requires a real cluster in CI.
    add_test(NAME fastgatk-workflow-local-contract
        COMMAND bash
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_pipeline_local.sh")
    set_tests_properties(fastgatk-workflow-local-contract PROPERTIES
        ENVIRONMENT "FASTGATK_WORKFLOW_BUILD_DIR=${CMAKE_BINARY_DIR}")
    add_test(NAME fastgatk-bqsr-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr.py")
    # The broad native BQSR contract is complemented by a focused pinned Java
    # oracle.  It compares all four GATKReport tables and Java-decoded
    # ApplyBQSR records, including --use-original-qualities/OQ, against both
    # configured Kokkos host backends.
    add_test(NAME fastgatk-bqsr-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>")
    # ContextCovariate width is owned by the recalibration report, not an
    # ApplyBQSR CLI default.  This high-coverage fixture makes a size-3 model
    # affect QUAL and compares both producer tables and consumer records.
    add_test(NAME fastgatk-bqsr-context-size-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_context_size_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-context-size-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-apply-bqsr-alias-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_bqsr_alias_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-bqsr-alias-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # BI/BD are optional SAM Z tags carrying per-base insertion/deletion
    # qualities.  Keep this focused oracle separate from the ordinary
    # Q45-fallback path so a future change cannot hide an event-table drift.
    add_test(NAME fastgatk-bqsr-indel-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_indel_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-indel-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Long reads exercise a real CycleCovariate memory/domain boundary.  The
    # oracle checks GATK's default rejection above cycle 500, then compares
    # the complete report after selecting a larger maximum on both Java and
    # native implementations (including inserted read bases).
    add_test(NAME fastgatk-bqsr-long-read-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_long_read_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-long-read-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # BaseRecalibrator's generic read-filter plugin is a high-impact
    # compatibility boundary: Java resolves defaults, -DF removal and -RF
    # additions after parsing.  The pinned oracle uses a synthetic SAM with
    # each MAPQ/flag case and compares all GATKReport tables.
    add_test(NAME fastgatk-bqsr-read-filter-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_read_filter_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-read-filter-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # ApplyBQSR's preserve threshold is an output semantic: low-quality bases
    # bypass both the hierarchical model and every quantizer.  Keep it as a
    # separate pinned gate so a transformed quality cannot hide this branch.
    add_test(NAME fastgatk-bqsr-preserve-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_preserve_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-preserve-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # The producer-side report gate makes Java ApplyBQSR consume the native
    # GATKReport directly, exercising fixed-width table serialization in
    # addition to the numerical table oracle above.
    add_test(NAME fastgatk-bqsr-report-roundtrip-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_report_roundtrip_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-report-roundtrip-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # CRAM/reference/index handling is a separate pinned oracle: unlike the
    # BAM-focused BQSR oracle, it consumes the shipped CRAM fixture, publishes
    # indexed CRAM, and reopens both outputs through HTSJDK.
    add_test(NAME fastgatk-bqsr-cram-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bqsr_cram_gatk_oracle.py")
    set_tests_properties(fastgatk-bqsr-cram-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_BQSR_BINARY=$<TARGET_FILE:fastgatk-bqsr>;FASTGATK_APPLY_BQSR_BINARY=$<TARGET_FILE:fastgatk-apply-bqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-gather-bqsr-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_bqsr_gatk_oracle.py")
    set_tests_properties(fastgatk-gather-bqsr-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_GATHER_BQSR_BINARY=$<TARGET_FILE:fastgatk-gather-bqsr-reports>")
    add_test(NAME fastgatk-analyze-covariates-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_analyze_covariates.py")
    add_test(NAME fastgatk-analyze-covariates-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_analyze_covariates.py")
    set_tests_properties(fastgatk-analyze-covariates-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_ANALYZE_COVARIATES_BINARY=$<TARGET_FILE:fastgatk-analyze-covariates>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-analyze-covariates-bqsr-alias-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_analyze_covariates_bqsr_alias_gatk_oracle.py")
    set_tests_properties(fastgatk-analyze-covariates-bqsr-alias-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_ANALYZE_COVARIATES_BINARY=$<TARGET_FILE:fastgatk-analyze-covariates>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    set_tests_properties(fastgatk-analyze-covariates-contract PROPERTIES
        ENVIRONMENT "FASTGATK_ANALYZE_COVARIATES_BINARY=$<TARGET_FILE:fastgatk-analyze-covariates>")
    add_test(NAME fastgatk-indel-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_indel.py")
    add_test(NAME fastgatk-hc-ploidy-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_ploidy.py")
    add_test(NAME fastgatk-hc-likelihood-filter-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_likelihood_filter.py")
    set_tests_properties(fastgatk-hc-likelihood-filter-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-variant-annotations-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_variant_annotations.py")
    set_tests_properties(fastgatk-hc-variant-annotations-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-genotype-priors-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_genotype_priors.py")
    set_tests_properties(fastgatk-hc-genotype-priors-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-broad-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_broad_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-broad-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # Pinned writer-boundary oracle: HC computes full genotype/reference-
    # confidence evidence, then optionally publishes an 8-column sites-only
    # VCF/GVCF.  Bare and separated Boolean forms are checked against GATK.
    add_test(NAME fastgatk-hc-sites-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_sites_only_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-sites-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # A small max-size window makes GATK's AssemblyRegion defaults observable:
    # pin the 100 bp halo/50 bp propagation boundary and the resulting
    # short-edge read ownership before the PairHMM oracle is run.
    add_test(NAME fastgatk-hc-assembly-region-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_assembly_region_boundary_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-assembly-region-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # CIGAR I/D events are explicit GATK AssemblyRegion activity seeds even
    # when aligned bases are reference-matching.  This bounded oracle catches
    # the ownership/fallback failure before graph/PairHMM numeric comparisons.
    add_test(NAME fastgatk-hc-cigar-indel-activity-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_cigar_indel_activity_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-cigar-indel-activity-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Given-alleles features force a single active locus and inject concrete
    # events into HC's local haplotype set; keep that Host/Kokkos boundary
    # pinned with a reference-only pileup, including filtered-feature policy.
    add_test(NAME fastgatk-hc-alleles-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_alleles_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-alleles-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-hc-min-pruning-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_min_pruning_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-min-pruning-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-min-base-quality-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_min_base_quality_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-min-base-quality-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Q0/Q1 is a separate source boundary: HC finalization passes a signed
    # minTailQuality and rebuilds FragmentCollection after soft-clip
    # reversion.  Pin whole-gVCF DP/block behavior, not only the Q20 graph
    # segmentation checked above.
    add_test(NAME fastgatk-hc-min-base-quality-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_min_base_quality_boundary_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-min-base-quality-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-hc-kmer-list-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_kmer_list_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-kmer-list-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # ``retainEvidence`` expands VariantContexts by this GATK-configurable
    # Host-side margin after PairHMM marginalization.  It must not be folded
    # into graph construction or a device-kernel window.
    add_test(NAME fastgatk-hc-informative-overlap-margin-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_informative_overlap_margin_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-informative-overlap-margin-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK issue3466: reverting a leading soft clip must first hard-clip the
    # prefix that would extend before coordinate one.  Keep this Host CIGAR
    # boundary separate from the device PairHMM numerical oracle.
    add_test(NAME fastgatk-hc-softclip-contig-start-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_softclip_contig_start_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-softclip-contig-start-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK issue3845 accepts a chrM-only reference dictionary with a complete
    # BAM dictionary.  It also pins default alignment-start reservoir sampling
    # and the candidate-local retained read set consumed by RMSMappingQuality.
    add_test(NAME fastgatk-hc-issue3845-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_issue3845_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-issue3845-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-hc-softclip-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_softclip_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-softclip-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-softclip-low-quality-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_softclip_low_quality_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-softclip-low-quality-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-hc-chr17-independent-regions
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr17_independent_regions.py")
    set_tests_properties(fastgatk-hc-chr17-independent-regions PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-gq-bands-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_gq_bands_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-gq-bands-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # This bounded all-sites oracle exercises ReferenceConfidenceModel after
    # AssemblyBasedCallerUtils.finalizeRegion() has handled clipped low-Q
    # mate overlap.  Keep it in CTest rather than only verify_all.sh: GVCF
    # RCM regressions otherwise evade the focused HC oracle suite.
    add_test(NAME fastgatk-hc-bp-resolution-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_bp_resolution_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-bp-resolution-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-rcm-pl-range-gatk-oracle
             COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_rcm_pl_range_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-rcm-pl-range-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-indel-zero-gatk-oracle
             COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_indel_zero_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-indel-zero-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-hc-rcm-realignment-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_rcm_realignment.py")
    set_tests_properties(fastgatk-hc-rcm-realignment-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-flow-hmer-option-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_flow_hmer_option.py")
    set_tests_properties(fastgatk-flow-hmer-option-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-polyploid-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_polyploid_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-polyploid-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-complex-multiallelic-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_complex_multiallelic_oracle.py")
    set_tests_properties(fastgatk-hc-complex-multiallelic-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-streaming-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_streaming.py")
    set_tests_properties(fastgatk-hc-streaming-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-hc-region-streaming-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_region_streaming.py")
    set_tests_properties(fastgatk-hc-region-streaming-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # Diagnostic oracle for the OPEN divergence documented in
    # fastgatk-native/evidence/gvcf-stream-by-region-divergence-20260910.md:
    # reference-confidence block granularity, span-deletion annotations
    # (RAW_MQandDP/SB), and spurious phasing differ between --stream-by-region
    # tiles, between -L windows, and from pinned GATK.  It deliberately reports
    # a JSON status of "diagnostic" and exits 0 so a known, unfixed divergence
    # does not turn the suite red; its output stays in the test log as the
    # standing reminder.  When the fix lands, convert its comparisons to hard
    # assertions (streamed == non-streamed == GATK) and add
    # FASTGATK_REQUIRE_GATK_ORACLE=1 to this property list.
    add_test(NAME fastgatk-hc-gvcf-stream-overlap-diagnostic
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gvcf_stream_overlapping_indels_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-gvcf-stream-overlap-diagnostic PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # STRICT gate for the fixed -L-window-invariance defect (D2): a published
    # allele must be annotated from the AssemblyRegion owner that actually
    # carries its PairHMM context, so RAW_MQandDP/SB at 20:10020680 must equal
    # pinned GATK's per-window value in every window.  This is a hard assertion
    # (the script exits non-zero on mismatch) and it is green as of the
    # twin-owner fix; it runs pinned GATK across 12 windows (~90 s).
    add_test(NAME fastgatk-hc-window-invariance-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_window_invariance_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-window-invariance-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Same defect family as the gate above, but in the arbitrary-ploidy gVCF
    # branch, where the identical owner-selection bug was found live after the
    # diploid fix landed.  Without this, that second instance would have no
    # regression guard at all.  Runs pinned GATK at --sample-ploidy 3 across 7
    # windows (~35 s).
    add_test(NAME fastgatk-hc-ploidy-window-invariance-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_ploidy_window_invariance_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-ploidy-window-invariance-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Gate for the --alleles overlapping forced-event emission defect: a forced
    # feature allele overlapping an assembled event used to be suppressed
    # because the spanning-deletion override clobbered its correctly computed
    # compensated QUAL (see the fix in calling_pipeline.cpp).  Strict: it exits
    # non-zero when a GATK row is missing from native output, and it also runs
    # the --drop-alleles controls so --alleles attribution is proven.
    add_test(NAME fastgatk-hc-alleles-overlap-gate-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_alleles_overlap_gate_oracle.py")
    set_tests_properties(fastgatk-hc-alleles-overlap-gate-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Gate for the spanning-deletion allele-frequency PRIOR: GATK's
    # AlleleFrequencyCalculator picks each allele's Dirichlet pseudocount from
    # its length, so the symbolic '*' deletion (1 bp) takes the SNP prior on a
    # 1 bp REF record.  Using the indel prior unconditionally changed QUAL/QD
    # for a supported SNP on a pileup that also carries an assembled deletion.
    add_test(NAME fastgatk-hc-span-del-qual-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_span_del_qual_oracle.py")
    set_tests_properties(fastgatk-hc-span-del-qual-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Same allele-length-prior class as the gate above, but on the gVCF
    # reference-confidence path: a gVCF record carrying the symbolic spanning
    # deletion '*' must use the SNP pseudocount on a 1 bp REF (htsjdk
    # Allele.SPAN_DEL.length() == 1), not the indel one.
    add_test(NAME fastgatk-hc-gvcf-symbolic-prior-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_gvcf_symbolic_prior_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-gvcf-symbolic-prior-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Same length-based prior rule as the gate above; arbitrary-ploidy ordinary-VCF '*' AF pseudocount (1 bp REF must take the SNP prior).
    add_test(NAME fastgatk-hc-arbitrary-ploidy-span-del-prior-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py")
    set_tests_properties(fastgatk-hc-arbitrary-ploidy-span-del-prior-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Same length-based prior rule as the gate above; polyploid / hidden-spanning gVCF '*' AF pseudocount.
    add_test(NAME fastgatk-hc-polyploid-gvcf-span-del-prior-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py")
    set_tests_properties(fastgatk-hc-polyploid-gvcf-span-del-prior-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Same length-based prior rule as the gate above; '*' genotype prior on a 1 bp REF (needs --genotype-assignment-method USE_POSTERIOR_PROBABILITIES).
    add_test(NAME fastgatk-hc-spanning-prior-genotype-gq-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py")
    set_tests_properties(fastgatk-hc-spanning-prior-genotype-gq-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # AF/MLEAF zero-value formatting: GATK's htsjdk VCFEncoder.formatVCFDouble
    # picks precision per value by magnitude, so an exact zero is "0.00" even
    # inside a list whose other entries are "0.500".  Strict gate.
    add_test(NAME fastgatk-hc-af-zero-format-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_af_zero_format_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-af-zero-format-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # ordinary-VCF multi-ALT owner-annotation guard (reached, no divergence observed).  Added in the registration-triage round; verified green on both backends.
    add_test(NAME fastgatk-hc-multialt-owner-annotation-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_multialt_owner_annotation_fixture_oracle.py")
    set_tests_properties(fastgatk-hc-multialt-owner-annotation-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Mutect2 tumor/normal replay parity on the DREAM chr20 fixture.  Added in the registration-triage round; verified green on both backends.
    add_test(NAME fastgatk-mutect2-recheck-normal-replay-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_recheck_normal_replay.py")
    set_tests_properties(fastgatk-mutect2-recheck-normal-replay-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Mutect2 joint AssemblyResultSet / matched-normal parity.  Added in the registration-triage round; verified green on both backends.
    add_test(NAME fastgatk-mutect2-recheck-joint-assembly-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_recheck_assembly_resultset_joint.py")
    set_tests_properties(fastgatk-mutect2-recheck-joint-assembly-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # A concrete gVCF variant record must NOT carry INFO/END (GATK hands the
    # variant VC straight to the block combiner; only reference blocks set END).
    add_test(NAME fastgatk-hc-gvcf-indel-end-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_gvcf_indel_end_gatk_oracle.py")
    set_tests_properties(fastgatk-hc-gvcf-indel-end-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # VariantRecalibrator INFO/culprit must be the NAME of the worst-performing
    # annotation (GATK VariantRecalibratorEngine.calculateWorstPerformingAnnotation),
    # not a model provenance string.  Strict gate.
    add_test(NAME fastgatk-variant-recalibrator-culprit-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_culprit_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-culprit-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # SelectVariants keeps an all-hom-ref record (ALT='.') under
    # --remove-unused-alternates and drops it only under --exclude-non-variants.
    add_test(NAME fastgatk-select-variants-refonly-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants_refonly_gatk_oracle.py")
    set_tests_properties(fastgatk-select-variants-refonly-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # VariantFiltration --apply-allele-specific-filters must use GATK's
    # AS_FilterStatus encoding (SITE placeholder + '|' separator).
    add_test(NAME fastgatk-variant-filtration-asfilterstatus-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration_asfilterstatus_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-filtration-asfilterstatus-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # With --apply-allele-specific-filters GATK routes EVERY rule (and mask/
    # cluster) through a per-ALT split context, so both FILTER and
    # AS_FilterStatus must match GATK in the flag path.
    add_test(NAME fastgatk-variant-filtration-flag-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration_flag_only_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-filtration-flag-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # ReblockGVCF --drop-low-quals: re-genotyping concrete variants must precede
    # the --rgq-threshold PL[0] conversion (ReblockGVCF.java:415-424 before :427).
    add_test(NAME fastgatk-reblock-gvcf-droplowqual-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_reblock_gvcf_droplowqual_gatk_oracle.py")
    set_tests_properties(fastgatk-reblock-gvcf-droplowqual-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GenotypeGVCFs dense mode must materialize GATK's monomorphic REF-only
    # no-call call; native used to abort when the output-allele subset
    # collapsed to REF.
    add_test(NAME fastgatk-genotype-gvcf-spandel-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_spandel_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-spandel-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # A malformed FORMAT record must fail loudly, not silently truncate the
    # output: the linear readers cannot distinguish htslib EOF (-1) from a
    # parse failure (-2), so an unchecked `bcf_read(...) == 0` loop ends the
    # traversal early while still exiting 0.
    add_test(NAME fastgatk-genotype-gvcf-malformed-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_malformed_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-malformed-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GenotypeGVCFs reverse-trims the alleles it emits
    # (GenotypeGVCFsEngine.java:167 -> GATKVariantContextUtils.java:1443-1521):
    # native published the merged REF/ALT untrimmed, so a locus whose alleles
    # share trailing bases (AAAA/AACA) came out two bases longer than GATK's
    # (AAA/AAC).  The gate pins the one-base-ALT guard, the clip-one-fewer
    # emptiness corner, the ordering against the output-allele subset, and the
    # untouched '*'.
    add_test(NAME fastgatk-genotype-gvcf-reverse-trim-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_reverse_trim_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-reverse-trim-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GenotypeGVCFs regenotypes a locus only when the merged INFO/DP is positive
    # (GenotypeGVCFsEngine.java:160 and again at :181, with the depth computed by
    # ReferenceConfidenceVariantContextMerger.calculateVCDepth()); native used to
    # publish a fully called variant row for a locus GATK never regenotypes.
    add_test(NAME fastgatk-genotype-gvcf-depth-gate-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_depth_gate_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-depth-gate-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GenotypingEngine refuses a locus whose only surviving ALT is the spanning
    # deletion unless --include-non-variant-sites is set (:172-175, before
    # recordDeletions at :178-179), and QualByDepth divides the sign-preserving
    # -10*log10PError, so a zero-confidence locus publishes QD=-0.00 while its
    # QUAL column reads 0.  Native published an extra '*' row in default mode
    # and a positive-zero QD in dense mode.
    add_test(NAME fastgatk-genotype-gvcf-star-only-locus-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_star_only_locus_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-star-only-locus-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # STRICT gate for the spanning-locus materialization gap: GATK's traversal
    # visits every reference coordinate a record spans (VariantLocusWalker.java:
    # 150-176) and dense mode publishes those loci, while native's locus list was
    # built from record starts only, so 139 positions of GATK's own chr20 gVCF
    # were missing.  Two cases inside it are REPORTED ONLY (a ~1e-16
    # sign-of-zero in QD of a regenotyped locus) and do not affect the exit
    # status; see the script docstring.
    add_test(NAME fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_dense_materialize_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-dense-spanning-loci-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # STRICT gate for the sample column of a materialized spanning locus: GATK
    # rewrites the projected call as hom-ref using the SOURCE record's GQ (99 at
    # 10008964/10008965, 41 at 10076991 -- the values differ on purpose, so a
    # constant cannot pass), while native re-read its already PL-derived GQ (0)
    # and, being guarded on gq > 0, emitted "./.:34:0" instead of "0/0:34:99".
    add_test(NAME fastgatk-genotype-gvcf-spanning-source-gq-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_spanning_source_gq_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-spanning-source-gq-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # STRICT gate for the allele-frequency posterior underflow: a confident
    # hom-alt sample's P(hom-ref) fell below the smallest positive double, so the
    # linear normalisation produced log10(0) = -inf, which collided with the
    # "no likelihoods" sentinel, pruned the ALT and dropped a real variant call
    # (18 of GATK's own 252 chr20 loci).
    add_test(NAME fastgatk-genotype-gvcf-confident-hom-alt-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_confident_hom_alt_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-confident-hom-alt-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Systemic guard: an htslib read failure (-2) must not be mistaken for EOF
    # (-1), which silently truncated output across many tools while exiting 0.
    add_test(NAME fastgatk-malformed-input-fail-loud-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_malformed_input_fail_loud_oracle.py")
    set_tests_properties(fastgatk-malformed-input-fail-loud-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # The GATK-compat header must be ordered the way htsjdk orders it (full-line
    # text sort, ##contig by input index, ##fileformat first, #CHROM last).
    add_test(NAME fastgatk-genotype-gvcf-header-order-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_header_order_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-header-order-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK applies its own site-level FILTER=LowQual to a dense row whose
    # recomputed (pre-posterior) call confidence is below the cutoff.
    add_test(NAME fastgatk-genotype-gvcf-lowqual-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_lowqual_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-lowqual-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # NOTE: verify_reblock_gvcf_triploid_gatk_oracle.py is deliberately NOT registered:
    # its divergence is confirmed real (GATK reblocks the triploid record to 0/0/0
    # END=30 with 4 PLs) but the fix is blocked by an unrelated, GATK-untrue block in
    # verify_reblock_gvcf.py; see .diag/round-reblock-gvcf.md.  Register it once that
    # fixture and its ~25 derived values are corrected.
    # A6 WGS benchmark oracle: runs the existing benchmark_end_to_end.py on the
    # canonical NA12878 chr17 69-70k fixture and verifies the JSON artifact
    # (schema_version=4, status=pass, 5 native cases + GATK baseline).
    add_test(NAME fastgatk-wgs-benchmark-a6-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_wgs_benchmark_a6_oracle.py")
    set_tests_properties(fastgatk-wgs-benchmark-a6-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR}/..")
    add_test(NAME fastgatk-hc-multi-input-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_multi_input.py")
    set_tests_properties(fastgatk-hc-multi-input-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    # Large real-data no-call oracle: this 100 kb CEUTrio window has 90k
    # reads and must remain byte-identical to GATK's empty HC VCF.  It guards
    # default read filtering and activity/assembly traversal at a scale that
    # the candidate-rich micro-fixtures cannot cover.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/chr20/reference/GRCh37.chr20.fa"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/ceutrio/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam")
        add_test(NAME fastgatk-hc-chr20-100k-nocall-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_100k_nocall_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-100k-nocall-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    add_test(NAME fastgatk-activity-region-controls-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_activity_region_controls.py")
    set_tests_properties(fastgatk-activity-region-controls-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-overlapping-quality-correction-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_overlapping_quality_correction.py")
    set_tests_properties(fastgatk-overlapping-quality-correction-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
    add_test(NAME fastgatk-real-assembly-graph-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_real_assembly_graph.py")
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../fixtures/chr20/ref20mnp.fasta"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../fixtures/chr20/mnp.bam")
        add_test(NAME fastgatk-hc-chr20-real-contract
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_real_contract.py")
        set_tests_properties(fastgatk-hc-chr20-real-contract PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>")
        add_test(NAME fastgatk-hc-chr20-min-pruning-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-min-pruning-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        # `--max-mnp-distance 1` creates an EventMap overlap whose gVCF
        # serialization retains `*` only after the preceding `<NON_REF>` row
        # enters GenotypingEngine's upstream-deletion state machine.  Pin the
        # complete REF/ALT/*/<NON_REF> matrix and phase boundary byte-for-byte.
        add_test(NAME fastgatk-hc-chr20-max-mnp-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-max-mnp-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        # The same EventMap overlap is also an arbitrary-ploidy boundary:
        # partitioned AssemblyRegion owners must supply the complete
        # REF/ALT/*/<NON_REF> matrix before the shared Kokkos PL/AF kernels
        # run.  Compare all data rows with a pinned triploid GATK oracle.
        add_test(NAME fastgatk-hc-chr20-max-mnp-polyploid-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-max-mnp-polyploid-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        # Haploid gVCF folds the weak predecessor to a reference block, but
        # retains it as an internal spanning-deletion state for the target
        # MNP's Kokkos AF/QUAL calculation.  The serialized target has no
        # `*`; GATK still reports QUAL 89.07 from that hidden state.
        add_test(NAME fastgatk-hc-chr20-max-mnp-haploid-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-max-mnp-haploid-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_HC_SAMPLE_PLOIDY=1")
        add_test(NAME fastgatk-hc-chr20-max-mnp-polyploid4-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-max-mnp-polyploid4-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_HC_SAMPLE_PLOIDY=4")
        # p8 catches two rounding-sensitive source boundaries: symbolic
        # <NON_REF> uses the HTSJDK indel prior, and reference blocks choose
        # their SNP/indel model from an unrounded hom-ref GQ.
        add_test(NAME fastgatk-hc-chr20-max-mnp-polyploid8-gvcf-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr20-max-mnp-polyploid8-gvcf-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_HC_SAMPLE_PLOIDY=8")
        # Dense default HC gVCF followed by GenotypeGVCFs.  This has phased
        # candidate islands separated by RCM blocks, so it catches both the
        # Host CIGAR projection and the phase-preserving joint-genotype
        # writer boundary that independent micro-fixtures cannot observe.
        add_test(NAME fastgatk-hc-dense-gvcf-genotype-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_dense_gvcf_genotype_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-dense-gvcf-genotype-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        add_test(NAME fastgatk-vqsr-scatter-joint-contract
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_vqsr_scatter_joint.py")
    endif()
    # B2 4-shard extension of the 2-shard scatter contract above: same
    # scatter→gather→apply pipeline with 4 shards × 5 variants, plus a
    # determinism check on the gather merge and on ApplyVQSR emission.
    add_test(NAME fastgatk-vqsr-scatter-joint-4shard-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_vqsr_scatter_joint_4shard_oracle.py")
    # A3 ≥5 real fixture pin (5th window): the canonical NA12878 chr17 69-70k
    # integration BAM emits 3 well-known GATK HaplotypeCaller calls
    # (17:69067 T/G 1/1, 17:69368 G/C 0/1, 17:69631 C/T 0/1) and fast-gatk
    # must match all three REF/ALT tuples exactly. Pairs with the 4 mnp.bam
    # windows in fastgatk-hc-chr20-real-contract to deliver the documented
    # "≥5 real fixture" oracle for HaplotypeCaller parity.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta")
        add_test(NAME fastgatk-hc-chr17-69k-70k-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_chr17_69k_70k_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-chr17-69k-70k-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        add_test(NAME fastgatk-hc-segmented-real-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_segmented_real_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-segmented-real-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    # A2 flow PairHMM real-data structural oracle: verifies the kernel-level
    # flow PairHMM GATK-raw-bit oracle passes (4 cases), the fastgatk-hc-call
    # CLI exposes --flow-assembly-collapse-hmer-size + --flow-assembly-collapse-
    # partial-mode, and the FlowBasedHaplotype chr9 BAM has flow-tag annotations.
    add_test(NAME fastgatk-flow-pairhmm-real-data-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_flow_pairhmm_real_data_oracle.py")
    set_tests_properties(fastgatk-flow-pairhmm-real-data-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-genotype-gvcf-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf.py")
    set_tests_properties(fastgatk-genotype-gvcf-contract PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>")
    add_test(NAME fastgatk-genotype-gvcf-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_genotype_gvcf.py")
    # The oracle script defaults to the OpenMP build directory for ad-hoc
    # invocation, but CTest must exercise the binary belonging to the current
    # configured backend as well (notably the Serial portability build).
    set_tests_properties(fastgatk-genotype-gvcf-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>")
    add_test(NAME fastgatk-genotype-gvcf-legacy-qual-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_genotype_gvcf_legacy_qual.py")
    set_tests_properties(fastgatk-genotype-gvcf-legacy-qual-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>")
    add_test(NAME fastgatk-genotype-gvcf-multisample-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_genotype_gvcf_multisample.py")
    set_tests_properties(fastgatk-genotype-gvcf-multisample-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>")
    add_test(NAME fastgatk-genotype-gvcf-multiallelic-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_genotype_gvcf_multiallelic.py")
    set_tests_properties(fastgatk-genotype-gvcf-multiallelic-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_COMBINE_BINARY=$<TARGET_FILE:fastgatk-combine-gvcfs>")
    add_test(NAME fastgatk-genotype-gvcf-inbreeding-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gatk_genotype_gvcf_inbreeding.py")
    set_tests_properties(fastgatk-genotype-gvcf-inbreeding-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_COMBINE_BINARY=$<TARGET_FILE:fastgatk-combine-gvcfs>")
    add_test(NAME fastgatk-genotype-gvcf-include-non-variant-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_include_non_variant_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-include-non-variant-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for GenotypeGVCFs' final output-allele subset.  The
    # AF posterior must remove an unsupported sibling ALT before PL/AD and
    # Number=A annotations are published; both accepted assignment spellings
    # are checked against GATK's PREFER_PLS engine policy.
    add_test(NAME fastgatk-genotype-gvcf-assignment-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_assignment_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-assignment-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for multi-sample reference-confidence blocks.  The
    # INFO/DP value must be recomputed from the merged sample-major FORMAT/DP
    # vector (not copied from the first source block), and aggregate and
    # --stream-by-locus must retain identical dense rows.
    add_test(NAME fastgatk-genotype-gvcf-multisample-reference-confidence-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_multisample_reference_confidence_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-multisample-reference-confidence-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for GenotypingEngine's pre-AF maximum-ALT boundary.
    # Likelihood-sum ranking, no-call/PL normalization for dropped alleles,
    # and aggregate versus --stream-by-locus output are compared exactly.
    add_test(NAME fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-max-alternate-alleles-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for upstream FORMAT/GP/PG vectors.  GenotypeGVCFs
    # forces PREFER_PLS and removes posterior fields during <NON_REF> allele
    # subsetting; aggregate and --stream-by-locus must match that final
    # cohort schema even when --gp-qual is explicitly supplied.
    add_test(NAME fastgatk-genotype-gvcf-gp-input-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_gp_input_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-gp-input-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for GATK's include-then-exclude interval operation.
    # The fixture places a concrete variant inside -XL and compares both
    # aggregate and bounded stream traversal, including manifest telemetry.
    add_test(NAME fastgatk-genotype-gvcf-exclude-intervals-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_exclude_intervals_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-exclude-intervals-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the deprecated GenotypeGVCFs STARTS_IN writer
    # mode. A long deletion starting before -L remains available to joint
    # traversal, but only finalized records whose POS is inside -L are emitted.
    add_test(NAME fastgatk-genotype-gvcf-starts-in-intervals-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-starts-in-intervals-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the bounded spanning-deletion merge boundary.
    # GATK removes an orphan '*' only after joint likelihood/QUAL calculation,
    # while retaining it when a concrete deletion record covers the locus.
    add_test(NAME fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genotype_gvcf_spanning_deletion_gatk_oracle.py")
    set_tests_properties(fastgatk-genotype-gvcf-spanning-deletion-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-reblock-gvcf-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_reblock_gvcf.py")
    # The synthetic contract covers option/error/index boundaries; keep the
    # real pinned GATK oracle in CTest as a separate regression so changes to
    # block merge, allele compaction, or annotation cleanup cannot silently
    # pass only the synthetic fixture.  The verifier skips when the bundled
    # GATK/reference corpus is unavailable, matching the other optional Java
    # oracle tests; CI can force it with FASTGATK_REQUIRE_GATK_ORACLE=1.
    add_test(NAME fastgatk-reblock-gvcf-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_reblock_gatk_oracle.py")
    set_tests_properties(fastgatk-reblock-gvcf-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REBLOCK_BINARY=$<TARGET_FILE:fastgatk-reblock-gvcf>")
    # Focused pinned oracle for the Java writer's buffered overlap behavior:
    # a concrete call must trim/split an input reference block before GQ-band
    # merging, including when the block and call arrive at the same locus.
    add_test(NAME fastgatk-reblock-gvcf-overlap-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_reblock_overlap_gatk_oracle.py")
    set_tests_properties(fastgatk-reblock-gvcf-overlap-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REBLOCK_BINARY=$<TARGET_FILE:fastgatk-reblock-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK's supported multi-input mode is non-overlapping shards from one
    # sample (multi-sample VCFs are rejected by the Java tool).  Keep that
    # boundary under a separate real oracle rather than conflating it with the
    # local multi-sample kernel contract.
    add_test(NAME fastgatk-reblock-gvcf-shards-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_reblock_gatk_shards.py")
    set_tests_properties(fastgatk-reblock-gvcf-shards-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REBLOCK_BINARY=$<TARGET_FILE:fastgatk-reblock-gvcf>")
    # Keep the Java/native file-boundary benchmark command surface aligned in
    # CTest without paying for JVM startup in the ordinary regression matrix.
    # The measured mode is reproducible on demand via the same script; its
    # output refuses a speedup claim for undersized workloads.
    add_test(NAME fastgatk-reblock-gvcf-java-benchmark-command-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/benchmark_reblock_gvcf_java.py"
                --dry-run)
    set_tests_properties(fastgatk-reblock-gvcf-java-benchmark-command-contract PROPERTIES
        ENVIRONMENT "FASTGATK_REBLOCK_BINARY=$<TARGET_FILE:fastgatk-reblock-gvcf>")
    add_test(NAME fastgatk-select-variants-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants.py")
    set_tests_properties(fastgatk-select-variants-contract PROPERTIES
        ENVIRONMENT "FASTGATK_SELECT_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-select-variants>;FASTGATK_EXPECTED_EXECUTION_SPACE=${FASTGATK_EXPECTED_EXECUTION_SPACE}")
    # Focused Java oracle for the common type/non-variant selection boundary,
    # including Barclay optional-boolean spellings of --exclude-non-variants.
    add_test(NAME fastgatk-select-variants-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants_gatk_oracle.py")
    set_tests_properties(fastgatk-select-variants-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SELECT_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-select-variants>")
    # Pinned Java oracle for the final sites-only writer boundary.  Keep
    # selection/subsetting on a full sample-bearing record, then strip FORMAT
    # and sample columns only when the output record is committed.
    add_test(NAME fastgatk-select-variants-sites-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants_sites_only_gatk_oracle.py")
    set_tests_properties(fastgatk-select-variants-sites-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SELECT_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-select-variants>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for genotype-level FT -> no-call conversion.  The
    # conversion must precede unused-ALT compaction so a filtered sole ALT
    # carrier is removed and AC/AN/AF are refreshed like GATK.
    add_test(NAME fastgatk-select-variants-filtered-nocall-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants_filtered_nocall_oracle.py")
    set_tests_properties(fastgatk-select-variants-filtered-nocall-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SELECT_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-select-variants>")
    # Pinned Java oracle for optional Boolean FILTER exclusion.  PASS and '.'
    # must remain selectable while any non-PASS FILTER is removed in true mode.
    add_test(NAME fastgatk-select-variants-filtered-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_select_variants_filtered_oracle.py")
    set_tests_properties(fastgatk-select-variants-filtered-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SELECT_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-select-variants>")
    add_test(NAME fastgatk-gather-vcfs-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_vcfs.py")
    # Pinned Picard/GATK oracle for uppercase INPUT/OUTPUT/QUIET/VERBOSITY
    # aliases together with comment, compression-level, and CREATE_INDEX
    # writer controls.
    add_test(NAME fastgatk-gather-vcfs-cli-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_vcfs_cli_boundary_gatk_oracle.py")
    set_tests_properties(fastgatk-gather-vcfs-cli-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-gather-vcfs>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-left-align-trim-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_left_align.py")
    set_tests_properties(fastgatk-left-align-trim-contract PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-left-align-trim>")
    # Release-pinned record oracle for split/right-trim, bounded left-shift,
    # genotype/annotation projection, and interval/END selection.
    add_test(NAME fastgatk-left-align-trim-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_left_align_gatk_oracle.py")
    set_tests_properties(fastgatk-left-align-trim-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-left-align-trim>")
    # Pinned Java oracle for the optional sites-only writer boundary.  The
    # normalization/splitting path still consumes full FORMAT/sample payload;
    # only the final VCF writer is reduced to eight site columns.
    add_test(NAME fastgatk-left-align-trim-sites-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_left_align_sites_only_gatk_oracle.py")
    set_tests_properties(fastgatk-left-align-trim-sites-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-left-align-trim>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned CLI/VariantWalker boundary oracle for include/exclude interval
    # padding, interval aliases/rules, Boolean-valued normalization switches,
    # and the -OVI index alias.
    add_test(NAME fastgatk-left-align-cli-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_left_align_cli_boundary_gatk_oracle.py")
    set_tests_properties(fastgatk-left-align-cli-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-left-align-trim>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-variant-filtration-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration.py")
    set_tests_properties(fastgatk-variant-filtration-contract PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_FILTRATION_BINARY=$<TARGET_FILE:fastgatk-variant-filtration>")
    # Keep a small release-pinned Java oracle separate from the broad native
    # contract.  It compares QUAL/INFO/FORMAT/GT/FT/FILTER semantics and the
    # inverted site-filter path on both configured Kokkos host backends.
    add_test(NAME fastgatk-variant-filtration-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-filtration-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_FILTRATION_BINARY=$<TARGET_FILE:fastgatk-variant-filtration>")
    # Pinned Java oracle for the Boolean-valued missing-annotation policy.
    # GATK accepts bare or separated true/false forms and rejects an embedded
    # '=' spelling; null checks remain independent of numeric missing policy.
    add_test(NAME fastgatk-variant-filtration-missing-boolean-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration_missing_boolean_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-filtration-missing-boolean-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_FILTRATION_BINARY=$<TARGET_FILE:fastgatk-variant-filtration>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned writer/semantic oracle for --set-filtered-genotype-to-no-call:
    # filtered called genotypes become unphased no-calls and AC/AN/AF are
    # recomputed only when at least one genotype changes.
    add_test(NAME fastgatk-variant-filtration-set-nocall-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_filtration_set_nocall_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-filtration-set-nocall-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_FILTRATION_BINARY=$<TARGET_FILE:fastgatk-variant-filtration>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-apply-vqsr-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr.py")
    set_tests_properties(fastgatk-apply-vqsr-contract PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-apply-vqsr>")
    # Release-pinned AS oracle: GATK's one-ALT recalibration records are
    # projected into exact Number=A vectors and scalar per-ALT culprits.
    add_test(NAME fastgatk-apply-vqsr-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-vqsr-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_APPLY_VQSR_BINARY=$<TARGET_FILE:fastgatk-apply-vqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned threshold oracle: a supplied tranches file does not activate
    # tranche walking unless an explicit truth-sensitivity target is given;
    # GATK falls back to VQSLOD >= 0 and LOW_VQSLOD otherwise.
    add_test(NAME fastgatk-apply-vqsr-default-cutoff-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr_default_cutoff_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-vqsr-default-cutoff-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_APPLY_VQSR_BINARY=$<TARGET_FILE:fastgatk-apply-vqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned writer-boundary oracle: GATK's --sites-only-vcf-output strips
    # FORMAT/sample columns only at publication, preserving VQSLOD/FILTER.
    add_test(NAME fastgatk-apply-vqsr-sites-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr_sites_only_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-vqsr-sites-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_APPLY_VQSR_BINARY=$<TARGET_FILE:fastgatk-apply-vqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned traversal oracle: ApplyVQSR applies -XL after the include set,
    # including -ixp's saturating interval padding, before score/filter work.
    add_test(NAME fastgatk-apply-vqsr-exclude-intervals-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr_exclude_intervals_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-vqsr-exclude-intervals-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_APPLY_VQSR_BINARY=$<TARGET_FILE:fastgatk-apply-vqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # ApplyVQSR's two filter-control switches are optional Booleans.  Their
    # false values are observable writer/score boundaries, not merely parser
    # aliases, so keep bare/true/false pinned to the Java command.
    add_test(NAME fastgatk-apply-vqsr-filter-booleans-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_apply_vqsr_filter_booleans_gatk_oracle.py")
    set_tests_properties(fastgatk-apply-vqsr-filter-booleans-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_APPLY_VQSR_BINARY=$<TARGET_FILE:fastgatk-apply-vqsr>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-variant-recalibrator-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator.py")
    # Release-pinned GATK model artifact oracle: native must load Java's
    # GATKReport tables, preserve annotation dimensions, and score the input
    # without claiming Java VBEM/raw-bit equivalence.
    add_test(NAME fastgatk-variant-recalibrator-gatk-model-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_gatk_model_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-gatk-model-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the actual VBEM update path.  The fixture uses a
    # single Gaussian to isolate JavaRandom initialization, Normal-Wishart
    # regularization, standard-deviation exclusion, and final model scoring.
    add_test(NAME fastgatk-variant-recalibrator-vbem-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_vbem_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-vbem-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for VariantDataManager's post-normalization
    # dimension reorder.  This order controls random covariance initialization
    # and must be preserved before Kokkos VBEM training.
    add_test(NAME fastgatk-variant-recalibrator-annotation-order-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_annotation_order_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-annotation-order-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the hidden sampling mode.  GATK samples at the
    # zero-based traversal boundary (first, then every Nth record) and emits
    # only those sampled records in the recalibration table; this is distinct
    # from merely omitting scores on unsampled records.
    add_test(NAME fastgatk-variant-recalibrator-sample-every-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_sample_every_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-sample-every-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned normalization-boundary oracle: GATK rejects a training
    # annotation with standard deviation < 1e-5; native must fail closed at
    # the same model-building boundary instead of substituting unit scale.
    add_test(NAME fastgatk-variant-recalibrator-zero-variance-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_zero_variance_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-zero-variance-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned CLI/runtime oracle: GATK keeps the model-build retry budget
    # (--max-attempts, default 1) separate from the VBEM iteration bound
    # (--max-iterations, default 150).  Native telemetry must preserve both
    # controls instead of treating attempts as iterations.
    add_test(NAME fastgatk-variant-recalibrator-attempt-iteration-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_recalibrator_attempt_iteration_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-recalibrator-attempt-iteration-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_RECALIBRATOR_BINARY=$<TARGET_FILE:fastgatk-variant-recalibrator>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-sort-sam-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_sort_sam.py")
    # Release-pinned Picard/GATK oracle for HTSJDK ordering and the
    # CREATE_INDEX default/explicit side-effect contract.
    add_test(NAME fastgatk-sort-sam-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_sort_sam_gatk_oracle.py")
    set_tests_properties(fastgatk-sort-sam-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SORT_SAM_BINARY=$<TARGET_FILE:fastgatk-sort-sam>")
    # Pinned HTSJDK SAMRecordDuplicateComparator oracle: library, unclipped
    # read/mate coordinates, orientation, mapped-end, score and name ordering.
    add_test(NAME fastgatk-sort-sam-duplicate-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_sort_sam_duplicate_gatk_oracle.py")
    set_tests_properties(fastgatk-sort-sam-duplicate-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SORT_SAM_BINARY=$<TARGET_FILE:fastgatk-sort-sam>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Picard/GATK oracle for the direct-replacement CLI boundary:
    # uppercase Picard aliases, bounded spill controls, writer compression
    # and validation settings, and optional Boolean switches.
    add_test(NAME fastgatk-sort-sam-cli-boundary-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_sort_sam_cli_boundary_gatk_oracle.py")
    set_tests_properties(fastgatk-sort-sam-cli-boundary-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SORT_SAM_BINARY=$<TARGET_FILE:fastgatk-sort-sam>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-mark-duplicates-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mark_duplicates.py")
    # Keep the release-pinned Picard comparison independent from the broad
    # native contract.  It checks duplicate flags, DT/PG tags, metrics and
    # duplicate-set histogram semantics, including stale-input cleanup and
    # CREATE_INDEX=false by default versus an explicit true sidecar.
    add_test(NAME fastgatk-mark-duplicates-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mark_duplicates_gatk_oracle.py")
    set_tests_properties(fastgatk-mark-duplicates-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MARK_DUPLICATES_BINARY=$<TARGET_FILE:fastgatk-mark-duplicates>")
    # Picard's OpticalOnly policy is an output-semantic boundary: only
    # optical members receive DT:Z:SQ while ordinary library duplicates remain
    # untagged.  Keep its pinned oracle separate from the broad flag/metrics
    # comparison so this policy cannot regress unnoticed.
    add_test(NAME fastgatk-mark-duplicates-optical-only-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mark_duplicates_tagging_policy_gatk_oracle.py")
    set_tests_properties(fastgatk-mark-duplicates-optical-only-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MARK_DUPLICATES_BINARY=$<TARGET_FILE:fastgatk-mark-duplicates>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Picard oracle for paired ReadEnds construction when MC is absent
    # and terminal clipping changes the unclipped five-prime coordinates.
    add_test(NAME fastgatk-mark-duplicates-pair-key-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mark_duplicates_pair_key_gatk_oracle.py")
    set_tests_properties(fastgatk-mark-duplicates-pair-key-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MARK_DUPLICATES_BINARY=$<TARGET_FILE:fastgatk-mark-duplicates>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-mutect2-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2.py")
    set_tests_properties(fastgatk-mutect2-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # Mutect2's default post-read transformer hard-clips terminal inverted
    # repeat artifacts.  This Host-boundary contract transcribes GATK's
    # PalindromeArtifactClipReadTransformer unit fixture and also pins the
    # public --ignore-itr-artifacts opt-out.
    add_test(NAME fastgatk-mutect2-itr-artifact-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_itr_artifact_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-itr-artifact-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # Mutect2Engine calls cleanOverlappingReadPairs with independent halves
    # of --pcr-snv-qual and --pcr-indel-qual.  Keep both Host-side caps pinned
    # before the Kokkos graph/PairHMM stages so the public indel option cannot
    # regress to the historical fixed Q20 behavior.
    add_test(NAME fastgatk-mutect2-pcr-overlap-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_pcr_overlap_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-pcr-overlap-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # SomaticGenotypingEngine changes its evidence key from read name to the
    # read itself under --independent-mates.  The fixture pins one paired
    # fragment becoming two Host groups without changing Kokkos PairHMM work.
    add_test(NAME fastgatk-mutect2-independent-mates-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_independent_mates_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-independent-mates-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # SomaticReferenceConfidenceModel produces <NON_REF> reference contexts
    # and SomaticGVCFBlockCombiner groups them by TLOD (not HC's GQ) bands.
    # Pin both default and explicit -LODB behavior to the bundled Java oracle.
    add_test(NAME fastgatk-mutect2-gvcf-reference-blocks-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_gvcf_reference_blocks_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-gvcf-reference-blocks-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # BP resolution shares the somatic reference model but serializes one
    # complete <NON_REF> record per locus instead of GVCF TLOD bands.
    add_test(NAME fastgatk-mutect2-bp-resolution-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_bp_resolution_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-bp-resolution-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # AssemblyResultSet.trimTo() assigns padded EventMap events to one active
    # core. Keep the concrete gVCF allele set pinned: halo events must remain
    # reference-confidence blocks, while low-TLOD core events are emitted.
    add_test(NAME fastgatk-mutect2-gvcf-eventmap-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_gvcf_eventmap_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-gvcf-eventmap-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Full EventMap traversal and a one-transaction direct interval must
    # preserve reference confidence on both sides of an assembled event: BP
    # writes every locus while GVCF bands the same Kokkos TLOD rows.
    add_test(NAME fastgatk-mutect2-reference-confidence-eventmap-full-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-reference-confidence-eventmap-full-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK 4.6.2's reference-confidence model accepts exactly one sample;
    # tumor+normal -ERC is an input error, not a two-sample gVCF mode.
    add_test(NAME fastgatk-mutect2-gvcf-matched-normal-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_gvcf_matched_normal_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-gvcf-matched-normal-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-mutect2-gvcf-multitumor-gatk-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_gvcf_multitumor_gatk_contract.py")
    set_tests_properties(fastgatk-mutect2-gvcf-multitumor-gatk-contract PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Keep the pinned Java integration oracle distinct from the synthetic CLI
    # contract.  It requires the real GATK 4.6.2.0 call set, full VCF/gVCF
    # data records and stable schema metadata (apart from execution
    # provenance), while separately exposing the PairHMM/assembly boundary.
    add_test(NAME fastgatk-mutect2-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # Mutect2 has a distinct Host downsampler from HaplotypeCaller: at an
    # overloaded alignment-start pool it reservoirs only reads with MAPQ > 50.
    # This chrM fixture pins the selected read set through an exact VCF oracle.
    add_test(NAME fastgatk-mutect2-issue3845-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_issue3845_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-issue3845-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # AssemblyRegionWalker expands an explicit -L selector by its 100bp Host
    # assembly context.  This narrow real mitochondrial fixture locks that
    # read-fetch boundary independently from the ongoing dense-graph VCF
    # alignment work, while ensuring the writer still crops back to -L.
    add_test(NAME fastgatk-mutect2-mito-interval-halo-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_mito_interval_halo_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-mito-interval-halo-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # In this high-depth chrM repeat GATK's 0.2-log10 best-haplotype tie
    # threshold selects an indel-bearing path and its composed CIGAR moves a
    # first mate before chrM:8877.  Pin the final orientation/fragment fields
    # so a log-domain conversion or a direct base-map realignment cannot add
    # the historical extra F2R1 alternate observation.
    add_test(NAME fastgatk-mutect2-mito-realign-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_mito_realign_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-mito-realign-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK's zero global-mismapping-rate value is invalid: negative values map
    # to the Double.NEGATIVE_INFINITY disabled-cap sentinel and positive values
    # enable the cap. Keep this release-pinned option boundary separate because
    # accepting zero as an uncapped mode would silently flatten all per-read
    # likelihoods and change TLOD; the verifier also probes the pinned Java
    # AlleleLikelihoods API to ensure -Infinity cells remain unchanged.
    add_test(NAME fastgatk-mutect2-mismapping-rate-boundary
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_mismapping_rate_boundary.py")
    set_tests_properties(fastgatk-mutect2-mismapping-rate-boundary PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    # Pinned GATK 4.6.2.0 oracle for --force-active: ActivityProfile
    # boundaries remain stable while every profile segment becomes active;
    # the native candidate site set must remain stable on the fixture.
    # A4 engine-parity contract (no GATK jar required): the native
    # fragment-matrix dump is replayed through a faithful Python transcription
    # of the Java SomaticLikelihoodsEngine so every emitted single-ALT TLOD
    # must be reproduced from the engine's own inputs to <1e-4.  This guards
    # the posterior/evidence math independently of upstream assembly parity.
    add_test(NAME fastgatk-mutect2-tlod-formula
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_tlod_formula.py")
    set_tests_properties(fastgatk-mutect2-tlod-formula PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>")
    add_test(NAME fastgatk-mutect2-force-active-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_force_active_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-force-active-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-mutect2-zero-lod-indel-activity-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_zero_lod_indel_activity_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-zero-lod-indel-activity-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-mutect2-feature-resource-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_feature_resource_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-feature-resource-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned GATK 4.6.2.0 oracle for --normal-lod: a matched-normal ALT
    # evidence threshold controls candidate emission without changing the
    # native caller's likelihood computation.
    add_test(NAME fastgatk-mutect2-normal-lod-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_normal_lod_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-normal-lod-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK's missing BI/BD fallback is raw Q45; the Conservative PCR model
    # adjusts only non-terminal positions, leaving terminal insertion/deletion
    # qualities at Q45.  Keep this PairHMM input-boundary oracle separate from
    # the broader Mutect2 posterior/assembly tolerance oracle.
    add_test(NAME fastgatk-mutect2-pairhmm-default-indel-quality-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_pairhmm_default_indel_quality_gatk_oracle.py")
    set_tests_properties(fastgatk-mutect2-pairhmm-default-indel-quality-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # A4 real-data oracle on the GATK 4.6.2.0 DREAM synthetic tumor/normal
    # fixture (b37, chr20, ~30k reads each, fully downloaded from
    # gatk-source/.../large/mutect/dream_synthetic_bams/). The oracle pins
    # end-to-end Mutect2 run + manifest schema + sample names + activity-range
    # coverage and full GATK call-position recall. Native-only calls are
    # reported by the oracle while the remaining precision gap is addressed.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/dream_synthetic/chr20/tumor.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/dream_synthetic/chr20/normal.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/downloads/reference/hs37d5.fa.gz")
        # Fast real-data guard for the low-BQ compound haplotype that GATK's
        # ReadThreadingAssembler retains at chr20:10004241 before normal
        # replay.  Keep this separate from the expensive full-chr20 oracle
        # so source-aligned recovery regressions fail quickly.
        add_test(NAME fastgatk-mutect2-dream-low-bq-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_dream_low_bq_gatk_oracle.py")
        set_tests_properties(fastgatk-mutect2-dream-low-bq-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        # Multiple -I files may reuse an @RG ID with distinct sample names.
        # GATK resolves that ID in the source file's header, and also retains
        # a selected matched-normal sample that has no reads in this window.
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/dream_synthetic/chr20/tumor_1.bam")
            add_test(NAME fastgatk-mutect2-multinormal-rg-namespace-gatk-oracle
                COMMAND ${Python3_EXECUTABLE}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_multinormal_rg_namespace_gatk_oracle.py")
            set_tests_properties(fastgatk-mutect2-multinormal-rg-namespace-gatk-oracle PROPERTIES
                ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
            # Mutect2's role model is not one --tumor-sample equals one
            # sample: all non-normal @RG/SM samples contribute to the joint
            # tumor likelihood matrix and retain individual FORMAT columns.
            add_test(NAME fastgatk-mutect2-multitumor-format-gatk-oracle
                COMMAND ${Python3_EXECUTABLE}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_multitumor_format_gatk_oracle.py")
            set_tests_properties(fastgatk-mutect2-multitumor-format-gatk-oracle PROPERTIES
                ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        endif()
        # Mutect2 inherits ReadThreadingAssembler's repeatable k-mer list.
        # Pin both its [10,25] default and an explicitly reordered request to
        # the Java caller so a single-k graph retry cannot regress silently.
        add_test(NAME fastgatk-mutect2-kmer-list-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_kmer_list_gatk_oracle.py")
        set_tests_properties(fastgatk-mutect2-kmer-list-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        add_test(NAME fastgatk-mutect2-dream-synthetic-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_dream_synthetic_oracle.py")
        set_tests_properties(fastgatk-mutect2-dream-synthetic-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_STRICT_DREAM_ORACLE=1")
    endif()
    # A4 real-data strict oracle on the GATK 4.6.2.0 HCC1143 tumor/normal
    # fixture (b37, chr20, downsampled to 1Mb synthetic Ns reference).
    # It pins end-to-end traversal plus the complete ordered VCF data-row
    # contract, rather than treating position overlap as sufficient.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/chr20/HCC1143_normal.bam")
        # HCC1143 includes CIGARs that overhang the deliberately truncated
        # 1Mb chr20 reference.  The walker must retain those reads for
        # padded-region ownership while restricting ActivityProfile traversal
        # to the reference dictionary end, exactly as Java does.
        add_test(NAME fastgatk-mutect2-hcc1143-reference-boundary-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_hcc1143_reference_boundary_oracle.py")
        set_tests_properties(fastgatk-mutect2-hcc1143-reference-boundary-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        # Fast, strict real tumor/normal probe.  Unlike the broader chr20
        # resource traversal below, this locks the complete GATK VCF row and
        # tumor-before-normal sample order at a known HCC1143 somatic site.
        add_test(NAME fastgatk-mutect2-hcc1143-spot-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_hcc1143_spot_gatk_oracle.py")
        set_tests_properties(fastgatk-mutect2-hcc1143-spot-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
        add_test(NAME fastgatk-mutect2-hcc1143-chr20-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_mutect2_hcc1143_chr20_oracle.py")
        set_tests_properties(fastgatk-mutect2-hcc1143-chr20-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    add_test(NAME fastgatk-combine-gvcfs-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_combine_gvcfs.py")
    # Pinned GATK 4.6.2.0 oracle for the supported single-sample-shard and
    # explicit --call-genotypes paths.  Keep it separate from the synthetic
    # contract so sample-union and FORMAT payload regressions cannot hide
    # behind self-authored expectations.
    add_test(NAME fastgatk-combine-gvcfs-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_combine_gvcfs_gatk_oracle.py")
    set_tests_properties(fastgatk-combine-gvcfs-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_COMBINE_GVCFS_BINARY=$<TARGET_FILE:fastgatk-combine-gvcfs>")
    # Pinned GATK 4.6.2.0 oracle for the interval/reference-confidence block
    # boundary: Java drops REF/<NON_REF> blocks whose END lies beyond an
    # interval end, while retaining the original full block when END is
    # contained.  Aggregate and --stream-merge must share this predicate.
    add_test(NAME fastgatk-combine-gvcfs-interval-refblock-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_combine_gvcfs_interval_refblock_gatk_oracle.py")
    set_tests_properties(fastgatk-combine-gvcfs-interval-refblock-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_COMBINE_GVCFS_BINARY=$<TARGET_FILE:fastgatk-combine-gvcfs>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned GATK 4.6.2.0 oracle for PL-less reference-confidence records.
    # CombineGVCFs must not synthesize FORMAT/PL values when every input call
    # lacks likelihoods, and pure reference blocks must not promote FORMAT/DP
    # to INFO/DP.  This is separate from the PL-bearing payload oracle.
    add_test(NAME fastgatk-combine-gvcfs-plless-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_combine_gvcfs_plless_oracle.py")
    set_tests_properties(fastgatk-combine-gvcfs-plless-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_COMBINE_GVCFS_BINARY=$<TARGET_FILE:fastgatk-combine-gvcfs>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-filter-mutect-calls-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_calls.py")
    set_tests_properties(fastgatk-filter-mutect-calls-contract PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>")
    add_test(NAME fastgatk-optional-boolean-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_optional_boolean_contract.py")
    set_tests_properties(fastgatk-optional-boolean-contract PROPERTIES
        ENVIRONMENT "FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-left-align-trim>")
    add_test(NAME fastgatk-filter-mutect-calls-contamination-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_calls_contamination_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-contamination-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>")
    # Pinned GATK 4.6.2.0 oracle for the NormalArtifactFilter posterior
    # threshold boundary.  A positive NALOD posterior below the effective
    # FilterMutectCalls error threshold must not become normal_artifact merely
    # because it is non-zero.
    add_test(NAME fastgatk-filter-mutect-calls-normal-artifact-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_normal_artifact_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-normal-artifact-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned GATK 4.6.2.0 oracle for the joint ErrorProbabilities learning
    # boundary.  ContaminationFilter is a NON_SOMATIC contributor and must be
    # present during OPTIMAL_F_SCORE threshold learning, not only final
    # filtering; the fixture checks the AD/POPAF crossing and learned range.
    add_test(NAME fastgatk-filter-mutect-calls-contamination-joint-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_contamination_joint_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-contamination-joint-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # ReadOrientationFilter is a variant-level ARTIFACT contributor to Java's
    # ErrorProbabilities and must participate in empirical threshold learning.
    # This oracle compares the orientation-enabled native threshold/pattern
    # with GATK 4.6.2.0 and a native no-prior baseline.
    add_test(NAME fastgatk-filter-mutect-calls-orientation-joint-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_orientation_joint_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-orientation-joint-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned GATK 4.6.2.0 oracle for the bounded GermlineFilter posterior
    # surface.  Keep this separate from the hard-filter/contamination oracle:
    # it checks exact population-AF corner handling, fail-closed missing POPAF,
    # and materialized Number=A GERMQ phred output on both backends.
    add_test(NAME fastgatk-filter-mutect-calls-germline-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_calls_germline_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-germline-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Release-pinned orientation writer oracle.  QualityUtils.errorProbToQual
    # bounds ROQ to [1,93]; keep this independent from the release-specific
    # full FilterMutectCalls posterior/threshold model.
    add_test(NAME fastgatk-filter-mutect-calls-orientation-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_orientation_gatk_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-orientation-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK's Barclay surface exposes the variant-index switch as both the
    # long name and -OVI.  Keep a release-pinned oracle for the separated
    # Boolean form used by existing wrappers, including the true/false
    # sidecar boundary.
    add_test(NAME fastgatk-filter-mutect-calls-variant-index-alias-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_variant_index_alias_gatk_oracle.py")
    set_tests_properties(fastgatk-filter-mutect-calls-variant-index-alias-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Java-parity (FASTGATK_FMC_JAVA_PASSES=1) acceptance surface: the env
    # SomaticClusteringModel learning path, validated against the same pinned
    # GATK oracles (contamination joint / contamination / germline / normal
    # artifact / orientation joint / orientation writer / variant-index).
    # The native-default mega-contract (verify_filter_mutect_calls.py) pins
    # default INFO-driven approximate-learning telemetry and is intentionally
    # NOT part of the Java-path surface.
    foreach(java_oracle
        contamination_joint_oracle
        contamination_oracle
        normal_artifact_oracle
        germline_oracle
        orientation_joint_oracle
        orientation_gatk_oracle
        variant_index_alias_gatk_oracle)
        if (java_oracle STREQUAL "contamination_oracle")
            set(java_oracle_script verify_filter_mutect_calls_contamination_oracle.py)
        elseif (java_oracle STREQUAL "germline_oracle")
            set(java_oracle_script verify_filter_mutect_calls_germline_oracle.py)
        else()
            set(java_oracle_script verify_filter_mutect_${java_oracle}.py)
        endif()
        add_test(NAME fastgatk-filter-mutect-calls-java-${java_oracle}
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/${java_oracle_script}")
        set_tests_properties(fastgatk-filter-mutect-calls-java-${java_oracle} PROPERTIES
            ENVIRONMENT "FASTGATK_FILTER_MUTECT_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_REQUIRE_GATK_ORACLE=1;FASTGATK_FMC_JAVA_PASSES=1")
    endforeach()
    # A5 real-data oracle: FilterMutectCalls joint empirical model + manifest
    # + index on the GATK 4.6.2.0 DREAM synthetic Mutect2 output. Self-contained:
    # regenerates the Mutect2 input from the DREAM BAMs + hs37d5.fa.gz.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/dream_synthetic/chr20/tumor.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/dream_synthetic/chr20/normal.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/downloads/reference/hs37d5.fa.gz")
        add_test(NAME fastgatk-filter-mutect-calls-dream-joint-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_dream_synthetic_joint_oracle.py")
        set_tests_properties(fastgatk-filter-mutect-calls-dream-joint-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_FILTER_MUTECT_CALLS_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    # A5 real-data oracle on the GATK 4.6.2.0 HCC1143 chr20 truncated fixture
    # (44 Mutect2 calls; exercises the joint empirical model on a denser set).
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/chr20/HCC1143_tumor.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/chr20/HCC1143_normal.bam")
        add_test(NAME fastgatk-filter-mutect-calls-hcc1143-joint-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_filter_mutect_hcc1143_joint_oracle.py")
        set_tests_properties(fastgatk-filter-mutect-calls-hcc1143-joint-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_FILTER_MUTECT_CALLS_BINARY=$<TARGET_FILE:fastgatk-filter-mutect-calls>;FASTGATK_MUTECT2_BINARY=$<TARGET_FILE:fastgatk-mutect2>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    add_test(NAME fastgatk-genomicsdb-import-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_import.py")
    # The adapter contract above uses a fake backend to exercise resource and
    # retry boundaries.  Keep a pinned Java oracle separate so publication is
    # also proven by reopening the resulting real TileDB/GenomicsDB workspace
    # with GATK GenotypeGVCFs.
    add_test(NAME fastgatk-genomicsdb-import-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_import_gatk_oracle.py")
    set_tests_properties(fastgatk-genomicsdb-import-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOMICSDB_IMPORT_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-import>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned sample-name-map oracle: GATK exposes the map's sample name in
    # callset.json and in downstream GenotypeGVCFs columns.  The native sparse
    # workspace must materialize a private header with the same name rather
    # than silently retaining the source VCF sample.
    add_test(NAME fastgatk-genomicsdb-import-sample-map-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_import_sample_map_gatk_oracle.py")
    set_tests_properties(fastgatk-genomicsdb-import-sample-map-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOMICSDB_IMPORT_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-import>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pin GATK's incremental-workspace boundary: update imports require an
    # existing workspace, reuse the initial interval set, and append callsets.
    # Native sparse workspaces use a deterministic atomic input-index rebuild;
    # this is useful for the native GenotypeGVCFs path but is not TileDB append.
    add_test(NAME fastgatk-genomicsdb-import-update-workspace-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_import_update_workspace_gatk_oracle.py")
    set_tests_properties(fastgatk-genomicsdb-import-update-workspace-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GENOMICSDB_IMPORT_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-import>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned interval-span oracle for the no-Java sparse workspace.  GATK
    # imports a complete record when its reference span overlaps -L; native
    # materialization must not silently retain the whole input or clip a
    # reference-confidence block at the writer boundary.  The oracle uses
    # the isolated true-workspace exporter, so only register it when that
    # optional ABI-matched helper was built.
    if (TARGET fastgatk-genomicsdb-export)
        add_test(NAME fastgatk-genomicsdb-import-native-interval-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_import_native_interval_gatk_oracle.py")
        set_tests_properties(fastgatk-genomicsdb-import-native-interval-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_GENOMICSDB_IMPORT_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-import>;FASTGATK_GENOMICSDB_BRIDGE_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-export>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    # Exercise the isolated legacy-ABI bridge in CTest whenever the optional
    # true GenomicsDB helper was built.  Native sparse workspaces are tested by
    # verify_genomicsdb_import.py and must not be mislabeled as TileDB.
    if (TARGET fastgatk-genomicsdb-export)
        add_test(NAME fastgatk-genomicsdb-bridge-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_bridge.py")
        set_tests_properties(fastgatk-genomicsdb-bridge-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_GENOMICSDB_BRIDGE_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-export>;FASTGATK_GENOTYPE_BINARY=$<TARGET_FILE:fastgatk-genotype-gvcf>;FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-genomicsdb-export>;FASTGATK_REQUIRE_GENOMICSDB_BRIDGE=1")
        # Keep the storage boundary executable: portable sparse workspaces are
        # consumed by the native reader, but neither the legacy GenomicsDB
        # query API nor GATK may mistake them for TileDB arrays.
        add_test(NAME fastgatk-genomicsdb-native-storage-boundary
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_genomicsdb_native_storage_boundary.py")
        set_tests_properties(fastgatk-genomicsdb-native-storage-boundary PROPERTIES
            ENVIRONMENT "FASTGATK_GENOMICSDB_IMPORT_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-import>;FASTGATK_GENOMICSDB_BRIDGE_BINARY=$<TARGET_FILE:fastgatk-genomicsdb-export>;FASTGATK_NATIVE_BUILD=$<TARGET_FILE_DIR:fastgatk-genomicsdb-export>;FASTGATK_REQUIRE_GENOMICSDB_BRIDGE=1")
    endif()
    add_test(NAME fastgatk-genomicsdb-import-benchmark
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/benchmark_genomicsdb_import.py")
    add_test(NAME fastgatk-variants-to-table-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variants_to_table.py")
    # Pin the molten writer's observed Java row-consumption behavior, AS
    # Number=R splitting, and -EMD presence semantics separately from the
    # broad native field contract above.
    add_test(NAME fastgatk-variants-to-table-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variants_to_table_gatk_oracle.py")
    set_tests_properties(fastgatk-variants-to-table-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANTS_TO_TABLE_BINARY=$<TARGET_FILE:fastgatk-variants-to-table>")
    add_test(NAME fastgatk-variant-eval-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_eval.py")
    # Release-pinned oracle for the bounded VariantEval report surface:
    # CountVariants aggregate/Sample rows are byte-identical to GATK 4.6.2.0,
    # and malformed VCF input must fail closed.
    add_test(NAME fastgatk-variant-eval-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_eval_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-eval-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_EVAL_BINARY=$<TARGET_FILE:fastgatk-variant-eval>")
    # Pinned Java oracle for genotype-aware --keep-ac0 promotion of
    # hom-ref/all-no-call SNP and INDEL records through CountVariants.
    add_test(NAME fastgatk-variant-eval-keep-ac0-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_eval_keep_ac0_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-eval-keep-ac0-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_EVAL_BINARY=$<TARGET_FILE:fastgatk-variant-eval>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for ValidationReport's comparison-driven site-status
    # confusion matrix, INFO/AC classification, and evaluation-sample subset.
    add_test(NAME fastgatk-variant-eval-validation-report-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_variant_eval_validation_report_gatk_oracle.py")
    set_tests_properties(fastgatk-variant-eval-validation-report-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VARIANT_EVAL_BINARY=$<TARGET_FILE:fastgatk-variant-eval>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-validate-variants-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_validate_variants.py")
    add_test(NAME fastgatk-validate-variants-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_validate_variants_gatk_oracle.py")
    set_tests_properties(fastgatk-validate-variants-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VALIDATE_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-validate-variants>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Release-pinned semantic guard: VariantContext.validateAlternateAlleles
    # excludes symbolic ALT alleles from the called-genotype usage check while
    # retaining the strict check for concrete alternates.
    add_test(NAME fastgatk-validate-variants-symbolic-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_validate_variants_symbolic_gatk_oracle.py")
    set_tests_properties(fastgatk-validate-variants-symbolic-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_VALIDATE_VARIANTS_BINARY=$<TARGET_FILE:fastgatk-validate-variants>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-get-pileup-summaries-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_get_pileup_summaries.py")
    set_tests_properties(fastgatk-get-pileup-summaries-contract PROPERTIES
        ENVIRONMENT "FASTGATK_GET_PILEUP_BINARY=$<TARGET_FILE:fastgatk-get-pileup-summaries>")
    add_test(NAME fastgatk-get-pileup-summaries-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_get_pileup_gatk_oracle.py")
    set_tests_properties(fastgatk-get-pileup-summaries-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GET_PILEUP_BINARY=$<TARGET_FILE:fastgatk-get-pileup-summaries>")
    add_test(NAME fastgatk-calculate-contamination-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_calculate_contamination.py")
    add_test(NAME fastgatk-calculate-contamination-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_calculate_contamination_gatk_oracle.py")
    add_test(NAME fastgatk-gather-pileup-summaries-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_pileup_summaries.py")
    add_test(NAME fastgatk-gather-pileup-summaries-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_pileup_gatk_oracle.py")
    set_tests_properties(fastgatk-gather-pileup-summaries-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_GATHER_PILEUP_BINARY=$<TARGET_FILE:fastgatk-gather-pileup-summaries>")
    add_test(NAME fastgatk-learn-read-orientation-model-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_learn_read_orientation_model.py")
    add_test(NAME fastgatk-learn-read-orientation-model-multisample-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_learn_read_orientation_model_multisample_gatk_oracle.py")
    set_tests_properties(fastgatk-learn-read-orientation-model-multisample-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_LEARN_ORIENTATION_BINARY=$<TARGET_FILE:fastgatk-learn-read-orientation-model>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-collect-read-counts-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_collect_read_counts.py")
    # Release-pinned GATK oracle for read-start counting and the default
    # mapped/non-duplicate/MAPQ>=30 filter.  Keep it separate from the broad
    # native contract so a regression in interval traversal cannot be hidden
    # by a self-authored expectation.
    add_test(NAME fastgatk-collect-read-counts-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_collect_read_counts_gatk_oracle.py")
    set_tests_properties(fastgatk-collect-read-counts-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_COLLECT_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-collect-read-counts>")
    # Focused HDF5SimpleCountCollection gate: compare sample, complete
    # sequence-dictionary metadata (including SQ tags), intervals and counts
    # through a Java 4.6.2.0 reader.  This is intentionally separate from the
    # read-start/filter oracle and does not promote the broader CNV/PoN scope.
    add_test(NAME fastgatk-hdf5-simple-count-collection-metadata
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hdf5_simple_count_collection.py")
    set_tests_properties(fastgatk-hdf5-simple-count-collection-metadata PROPERTIES
        ENVIRONMENT "FASTGATK_COLLECT_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-collect-read-counts>")
    add_test(NAME fastgatk-denoise-read-counts-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_denoise_read_counts.py")
    # Release-pinned GATK oracle for the HDF5 SVD PoN case-input boundary:
    # GATK requires the case interval list to be identical to the original
    # PoN interval list before filtering, rather than merely containing all
    # post-filter panel intervals.
    add_test(NAME fastgatk-denoise-read-counts-interval-identity-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_denoise_read_counts_interval_identity_gatk_oracle.py")
    set_tests_properties(fastgatk-denoise-read-counts-interval-identity-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_DENOISE_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-denoise-read-counts>;FASTGATK_CREATE_PON_BINARY=$<TARGET_FILE:fastgatk-create-read-count-panel-of-normals>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Release-pinned oracle for SimpleCountCollection's integer COUNT codec.
    # Decimal/scientific or signed-int-overflow tokens must be rejected before
    # fractional-coverage normalization, matching GATK DataLine.getInt().
    add_test(NAME fastgatk-denoise-read-counts-integer-input-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_denoise_read_counts_integer_input_gatk_oracle.py")
    set_tests_properties(fastgatk-denoise-read-counts-integer-input-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_DENOISE_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-denoise-read-counts>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Release-pinned oracle for DenoiseReadCounts HDF5 input metadata.  The
    # SimpleCountCollection sequence dictionary is stored outside the count
    # matrix; both CopyRatioCollection outputs must preserve the complete
    # @SQ metadata exactly as GATK does.
    add_test(NAME fastgatk-denoise-read-counts-hdf5-metadata-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py")
    set_tests_properties(fastgatk-denoise-read-counts-hdf5-metadata-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_DENOISE_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-denoise-read-counts>;FASTGATK_COLLECT_READ_COUNTS_BINARY=$<TARGET_FILE:fastgatk-collect-read-counts>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-create-read-count-panel-of-normals-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_create_read_count_panel_of_normals.py")
    # Release-pinned GATK oracle for the failure boundary where a
    # multi-sample, SVD-requested panel has no non-zero singular values.
    # Both engines must reject before writing an apparently valid PoN.
    add_test(NAME fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py")
    set_tests_properties(fastgatk-create-read-count-panel-of-normals-degenerate-svd-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CREATE_PON_BINARY=$<TARGET_FILE:fastgatk-create-read-count-panel-of-normals>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Release-pinned oracle for v7 PoN sample metadata: GATK stores absolute
    # input count-file paths in both original and filtered panel arrays, not
    # the @RG SM identifier.  This is required for reusable-panel row
    # provenance/alignment in downstream DenoiseReadCounts workflows.
     add_test(NAME fastgatk-create-read-count-panel-of-normals-sample-metadata-gatk-oracle
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py")
     set_tests_properties(fastgatk-create-read-count-panel-of-normals-sample-metadata-gatk-oracle PROPERTIES
         ENVIRONMENT "FASTGATK_CREATE_PON_BINARY=$<TARGET_FILE:fastgatk-create-read-count-panel-of-normals>;FASTGATK_REQUIRE_GATK_ORACLE=1")
     # B4 byte-level structure parity oracle: the GATK 4.6.2.0 reference ships
     # wes-no-gc.pon.hdf5 (the SM-74NEG/SM-74P4M b37 chr20-truncated SimpleCount
     # Collection HDF5 panel-of-normals, 2 samples × 90 panel intervals). This
     # oracle pins that the v7 layout (groups, dataset shapes, sample_filenames,
     # singular_values, sequence_dictionary) is byte-stable and matches the
     # truncated reference .dict's @SQ line. Acts as the upper bound of any
     # PoN consumer's structure expectations.
     if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/wes-no-gc.pon.hdf5")
         add_test(NAME fastgatk-pon-gatk-bundled-byte-structure-oracle
             COMMAND ${Python3_EXECUTABLE}
                     "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_pon_gatk_bundled_byte_structure_oracle.py")
         set_tests_properties(fastgatk-pon-gatk-bundled-byte-structure-oracle PROPERTIES
             ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
     endif()
    add_test(NAME fastgatk-call-copy-ratio-segments-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_call_copy_ratio_segments.py")
    add_test(NAME fastgatk-call-copy-ratio-segments-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_call_copy_ratio_segments_gatk_oracle.py")
    # The oracle also exercises GATK's full z-score option names in addition
    # to the native short aliases, while retaining sidecar/degenerate checks.
    set_tests_properties(fastgatk-call-copy-ratio-segments-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CALL_COPY_RATIO_BINARY=$<TARGET_FILE:fastgatk-call-copy-ratio-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-call-copy-ratio-segments-interval-validation-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py")
    set_tests_properties(fastgatk-call-copy-ratio-segments-interval-validation-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CALL_COPY_RATIO_BINARY=$<TARGET_FILE:fastgatk-call-copy-ratio-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # GATK's TSV decoder accepts NaN and +/-Infinity in MEAN_LOG2_COPY_RATIO;
    # SimpleCopyRatioCaller preserves those values through IEEE-754 statistics
    # and emits neutral calls.  Keep this non-finite pass-through pinned to the
    # release Java output on every backend.
    add_test(NAME fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py")
    set_tests_properties(fastgatk-call-copy-ratio-segments-nonfinite-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CALL_COPY_RATIO_BINARY=$<TARGET_FILE:fastgatk-call-copy-ratio-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # DoubleStream.sum returns its high-order accumulator minus the negated
    # compensation word.  This pinned threshold fixture changes CALL if the
    # low-order word is discarded, so compare both emitted tables byte-exact.
    add_test(NAME fastgatk-call-copy-ratio-segments-compensated-sum-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py")
    set_tests_properties(fastgatk-call-copy-ratio-segments-compensated-sum-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_CALL_COPY_RATIO_BINARY=$<TARGET_FILE:fastgatk-call-copy-ratio-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # B5 CNV somatic e2e pipeline oracle: runs CollectReadCounts →
    # DenoiseReadCounts → ModelSegments → CallCopyRatioSegments end-to-end on
    # the HCC1143 chr20-truncated fixture using the GATK 4.6.2.0 reference
    # PoN artifact (wes-no-gc.pon.hdf5) as input. Proves the native pipeline
    # plumbing works against a real GATK-bundled artifact.
    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/HCC1143-t1-chr20-downsampled.deduplicated.bam"
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../testdata/real/cnv_somatic/wes-no-gc.pon.hdf5")
        add_test(NAME fastgatk-cnv-somatic-e2e-hcc1143-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_cnv_somatic_e2e_hcc1143_oracle.py")
        set_tests_properties(fastgatk-cnv-somatic-e2e-hcc1143-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    endif()
    add_test(NAME fastgatk-model-segments-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments.py")
    # Release-pinned oracle for the advanced --segments input: GATK skips
    # kernel segmentation and models the supplied Picard interval partition.
    add_test(NAME fastgatk-model-segments-input-segments-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_input_segments_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-input-segments-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the copy-ratio conditional-model slice.  It locks
    # the midpoint partition and GATK 4.6.2.0 global-parameter report shape;
    # native MCMC draws remain deterministic and are not claimed bit-identical.
    add_test(NAME fastgatk-model-segments-copy-ratio-conditionals-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_copy_ratio_conditionals_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-copy-ratio-conditionals-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # The same pinned fixture includes real allelic-count evidence and checks
    # finite MAF posterior deciles plus the GATK-shaped AF global report.  Keep
    # an explicit AF test name so release gates can select this model slice
    # independently from copy-ratio conditional changes.
    add_test(NAME fastgatk-model-segments-allele-fraction-conditionals-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_copy_ratio_conditionals_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-allele-fraction-conditionals-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the real AlleleFractionInitializer starting
    # state: integrated alt-minor/ref-minor beta responsibility plus GATK's
    # flat-prior pseudocount, implemented with a Kokkos RangePolicy.
    add_test(NAME fastgatk-model-segments-allele-fraction-initialization-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_allele_fraction_initialization_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-allele-fraction-initialization-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Bounded oracle for the count-backed AlleleFractionLikelihoods Laplace
    # conditional.  It compares the native conditional mean-bias state with
    # Java's modelBegin report while intentionally leaving stochastic posterior
    # draws non-bit-identical.
    add_test(NAME fastgatk-model-segments-allele-fraction-likelihood-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_allele_fraction_likelihood_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-allele-fraction-likelihood-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for GATK MultidimensionalModeller's
    # credible-interval adjacent-segment smoothing.  The fixture checks both
    # merge and strict non-merge decisions on the same supplied partition.
    add_test(NAME fastgatk-model-segments-smoothing-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_smoothing_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-smoothing-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for GATK's default segmentation dispatch.  Without
    # --segments, ModelSegments must select KernelSegmenter; the fixture also
    # exercises the exact rank-one linear copy-ratio feature used by the
    # Kokkos path.
    add_test(NAME fastgatk-model-segments-default-kernel-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_default_kernel_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-default-kernel-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Combined CR+AF segmentation consumes the first coordinate-ordered
    # oriented ALT fraction per copy-ratio interval.  This differs from both
    # folding to MAF and averaging all overlapping allelic sites.
    add_test(NAME fastgatk-model-segments-first-alt-fraction-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_first_alt_fraction_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-first-alt-fraction-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Multisample KernelSegmenter consumes the complete per-sample vector.  A
    # mean-reduced implementation loses anti-correlated copy-ratio events;
    # this pinned fixture covers concordant and anti-correlated boundaries.
    add_test(NAME fastgatk-model-segments-multisample-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_multisample_gatk_oracle.py")
    set_tests_properties(fastgatk-model-segments-multisample-gatk-oracle PROPERTIES
         ENVIRONMENT "FASTGATK_MODEL_SEGMENTS_BINARY=$<TARGET_FILE:fastgatk-model-segments>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    # B3 ModelSegments MCMC determinism + convergence contract: 2 identical
    # runs are byte-equal; 25 vs 100 sample chains converge on the small
    # HCC1143 chr20 fixture; 4-run byte stability verifies the random
    # seed is independent of the output prefix.
    add_test(NAME fastgatk-model-segments-mcmc-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_model_segments_mcmc_oracle.py")
    set_tests_properties(fastgatk-model-segments-mcmc-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-gather-tranches-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_gather_tranches.py")
    add_test(NAME fastgatk-collect-allelic-counts-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_collect_allelic_counts.py")
    add_test(NAME fastgatk-collect-allelic-counts-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_collect_allelic_counts_gatk_oracle.py")
    set_tests_properties(fastgatk-collect-allelic-counts-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_COLLECT_ALLELIC_COUNTS_BINARY=$<TARGET_FILE:fastgatk-collect-allelic-counts>;FASTGATK_REQUIRE_GATK_ORACLE=1")
    add_test(NAME fastgatk-collect-f1r2-counts-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_collect_f1r2_counts.py")
    set_tests_properties(fastgatk-collect-f1r2-counts-contract PROPERTIES
        ENVIRONMENT "FASTGATK_COLLECT_F1R2_COUNTS_BINARY=$<TARGET_FILE:fastgatk-collect-f1r2-counts>")
    add_test(NAME fastgatk-depth-of-coverage-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_depth_of_coverage.py")
    add_test(NAME fastgatk-depth-of-coverage-multisample-contract
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_depth_of_coverage_multisample.py")
    # Pinned Java oracle for the DepthOfCoverage default filter boundary:
    # QC-fail reads are retained, while WellformedReadFilter rejects a whole
    # read containing the CIGAR reference-skip (N) operator.
    add_test(NAME fastgatk-depth-of-coverage-read-filter-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_depth_of_coverage_read_filter_gatk_oracle.py")
    set_tests_properties(fastgatk-depth-of-coverage-read-filter-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")
    # Pinned Java oracle for the include-deletions/ignore-deletion-sites
    # interaction.  GATK parses both switches independently but the latter
    # disables deletion pileups from effective depth when both are requested.
    add_test(NAME fastgatk-depth-of-coverage-ignore-deletion-sites-gatk-oracle
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py")
    set_tests_properties(fastgatk-depth-of-coverage-ignore-deletion-sites-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_REQUIRE_GATK_ORACLE=1")

    # Every contract script that resolves binaries through
    # FASTGATK_NATIVE_BUILD must exercise the binary belonging to this CMake
    # tree.  This is especially important for the Serial portability build:
    # without the shared environment, scripts silently fall back to the
    # sibling OpenMP build and a green Serial CTest is not evidence of Serial
    # execution.  Preserve per-test target overrides (HC/Genotype/etc.) while
    # appending the generic build root.
    # Contract scripts that assert a concrete Kokkos execution-space label
    # (currently SelectVariants) need the same label that the configured
    # backend emits.  Keep this in CMake so a Serial CTest cannot accidentally
    # inherit the OpenMP default used for ad-hoc invocations.
    if (Kokkos_ENABLE_CUDA)
        set(_fastgatk_expected_execution_space "Cuda")
    elseif (Kokkos_ENABLE_HIP)
        set(_fastgatk_expected_execution_space "HIP")
    elseif (Kokkos_ENABLE_SYCL)
        set(_fastgatk_expected_execution_space "SYCL")
    elseif (Kokkos_ENABLE_OPENMP)
        set(_fastgatk_expected_execution_space "OpenMP")
    elseif (Kokkos_ENABLE_SERIAL)
        set(_fastgatk_expected_execution_space "Serial")
    else()
        set(_fastgatk_expected_execution_space "Host")
    endif()
    get_property(_fastgatk_all_tests DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY TESTS)
    foreach(_fastgatk_test IN LISTS _fastgatk_all_tests)
        get_test_property(${_fastgatk_test} ENVIRONMENT _fastgatk_test_environment)
        if(_fastgatk_test_environment STREQUAL "NOTFOUND")
            set(_fastgatk_test_environment "")
        elseif(NOT _fastgatk_test_environment STREQUAL "")
            set(_fastgatk_test_environment "${_fastgatk_test_environment};")
        endif()
        set_tests_properties(${_fastgatk_test} PROPERTIES
            ENVIRONMENT "${_fastgatk_test_environment}FASTGATK_NATIVE_BUILD=${CMAKE_BINARY_DIR};FASTGATK_EXPECTED_EXECUTION_SPACE=${_fastgatk_expected_execution_space}")
    endforeach()
endif()

if (FASTGATK_INSTALL_NATIVE_TOOLS)
    # Install every production executable declared in this directory while
    # excluding test-only helpers.  The target query keeps this list in sync
    # as tools are added and remains valid when HTSlib-dependent targets are
    # disabled at configure time.
    get_property(_fastgatk_native_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
    set(_fastgatk_install_targets)
    foreach (_fastgatk_target IN LISTS _fastgatk_native_targets)
        get_target_property(_fastgatk_target_type "${_fastgatk_target}" TYPE)
        if (_fastgatk_target_type STREQUAL "EXECUTABLE" AND
            NOT _fastgatk_target MATCHES "(^|-)smoke$|(^|-)contract$|(^|-)test$" AND
            NOT _fastgatk_target STREQUAL "fastgatk-genomicsdb-export")
            list(APPEND _fastgatk_install_targets "${_fastgatk_target}")
        endif()
    endforeach()
    if (_fastgatk_install_targets)
        install(TARGETS ${_fastgatk_install_targets}
                RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    endif()

    # Keep the dispatcher and registry together under share/, while the
    # generated wrappers derive the prefix from their own location.  This
    # makes a tarball/Conda/module install relocatable and avoids embedding
    # the build directory into the installed registry.
    install(FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/dispatcher/fastgatk.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/dispatcher/tool_registry.json"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/fastgatk/dispatcher")
    install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/dispatcher/fastgatk-installed.in"
            DESTINATION "${CMAKE_INSTALL_BINDIR}" RENAME fastgatk)
    install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/dispatcher/fastgatk-installed.in"
            DESTINATION "${CMAKE_INSTALL_BINDIR}" RENAME gatk)
endif()
