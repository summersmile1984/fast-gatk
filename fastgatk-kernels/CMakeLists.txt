cmake_minimum_required(VERSION 3.22)
project(fastgatk_kernels LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FASTGATK_KERNELS_BUILD_TESTS "Build the reusable kernel API smoke test" OFF)

# PairHMMModel's Java initialization uses a finite Jacobian/transition table.
# Keep the pinned table beside the source tree as the deterministic default so
# a production build does not silently fall back to platform libm rounding.
# FAST_GATK_PAIRHMM_TABLES remains an explicit runtime override for packaged or
# regenerated tables.
set(FASTGATK_PAIRHMM_TABLES_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/../pairhmm-demo/results/gatk-tables.hex"
    CACHE FILEPATH "Java-generated GATK PairHMM transition table")
option(FASTGATK_USE_PINNED_PAIRHMM_TABLES
       "Use the pinned Java-generated PairHMM table by default" ON)
option(FASTGATK_INSTALL_PAIRHMM_TABLES
       "Install the pinned Java-generated PairHMM table with the kernel artifacts" ON)
set(FASTGATK_PAIRHMM_TABLES_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/fastgatk"
    CACHE PATH "Install directory for the pinned PairHMM table")

if (NOT TARGET Kokkos::kokkos)
    message(FATAL_ERROR "fastgatk-kernels requires a parent Kokkos::kokkos target")
endif()
if (NOT TARGET fastgatk-core)
    message(FATAL_ERROR "fastgatk-kernels requires a parent fastgatk-core target")
endif()

add_library(fastgatk-kernels STATIC
    src/pairhmm_compat.cpp
    src/pairhmm_kokkos.cpp
    src/smith_waterman.cpp
    src/kmer_graph.cpp
    src/seqgraph_pruning.cpp
    src/activity_profile.cpp
    src/read_filter.cpp
    src/read_error_correction.cpp
    src/pileup_error_correction.cpp
    src/bqsr.cpp
    src/smith_waterman_kokkos.cpp
    src/genotype.cpp
    src/reference_confidence.cpp
    src/somatic.cpp
    src/filter_model_kokkos.cpp)
target_include_directories(fastgatk-kernels PUBLIC include)
target_link_libraries(fastgatk-kernels PUBLIC Kokkos::kokkos fastgatk-core)
target_compile_features(fastgatk-kernels PUBLIC cxx_std_20)
if (FASTGATK_USE_PINNED_PAIRHMM_TABLES AND EXISTS "${FASTGATK_PAIRHMM_TABLES_FILE}")
    # The absolute source path is intentional for source-tree builds: the
    # table is an auditable artifact generated from the pinned GATK jar.  An
    # install/package build can set FASTGATK_PAIRHMM_TABLES_FILE to its staged
    # data path, or use FAST_GATK_PAIRHMM_TABLES at runtime.
    if (IS_ABSOLUTE "${FASTGATK_PAIRHMM_TABLES_INSTALL_DIR}")
        set(_fastgatk_pairhmm_installed_table
            "${FASTGATK_PAIRHMM_TABLES_INSTALL_DIR}/gatk-tables.hex")
    else()
        set(_fastgatk_pairhmm_installed_table
            "${CMAKE_INSTALL_PREFIX}/${FASTGATK_PAIRHMM_TABLES_INSTALL_DIR}/gatk-tables.hex")
    endif()
    # Keep a prefix-independent candidate as well.  ``cmake --install
    # --prefix`` may relocate a package after configure time, but the
    # executable and data directory retain this relative layout.
    file(RELATIVE_PATH _fastgatk_pairhmm_installed_relative
         "${CMAKE_INSTALL_FULL_BINDIR}" "${_fastgatk_pairhmm_installed_table}")
    target_compile_definitions(fastgatk-kernels PRIVATE
        FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH="${FASTGATK_PAIRHMM_TABLES_FILE}"
        FASTGATK_INSTALLED_PAIRHMM_TABLES_PATH="${_fastgatk_pairhmm_installed_table}"
        FASTGATK_INSTALLED_PAIRHMM_TABLES_RELATIVE_PATH="${_fastgatk_pairhmm_installed_relative}")
    if (FASTGATK_INSTALL_PAIRHMM_TABLES)
        # The table is data, not generated code: ship the exact Java artifact
        # beside installed binaries so strict mode survives leaving the
        # source tree.  CMAKE_INSTALL_FULL_DATADIR is also compiled into the
        # loader as a fallback for package layouts that do not retain sources.
        install(FILES "${FASTGATK_PAIRHMM_TABLES_FILE}"
                DESTINATION "${FASTGATK_PAIRHMM_TABLES_INSTALL_DIR}"
                RENAME gatk-tables.hex)
    endif()
endif()
target_compile_options(fastgatk-kernels PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

if (FASTGATK_KERNELS_BUILD_TESTS)
    enable_testing()
    add_executable(fastgatk-kernels-api-smoke tests/api_smoke.cpp)
    target_link_libraries(fastgatk-kernels-api-smoke PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-kernels-api-smoke PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-kernels-api-smoke COMMAND fastgatk-kernels-api-smoke)
    add_executable(fastgatk-mutect2-activity-kernel-oracle
        tests/mutect2_activity_oracle.cpp)
    target_link_libraries(fastgatk-mutect2-activity-kernel-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-mutect2-activity-kernel-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-mutect2-activity-kernel-oracle
        COMMAND fastgatk-mutect2-activity-kernel-oracle)
    add_executable(fastgatk-flow-pairhmm-gatk-oracle tests/flow_pairhmm_gatk_oracle.cpp)
    target_link_libraries(fastgatk-flow-pairhmm-gatk-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-flow-pairhmm-gatk-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-flow-pairhmm-gatk-oracle COMMAND fastgatk-flow-pairhmm-gatk-oracle)
    # Replay the exact variable-length read/haplotype/quality rows emitted by
    # GATK's --pair-hmm-results-file through the same bucketed Kokkos kernel.
    # GATK writes the expected value with `%e` (six decimals); the executable
    # enforces only that documented serialization bound, keeping assembly
    # request-set differences separate from PairHMM arithmetic differences.
    add_executable(fastgatk-pairhmm-results-oracle tests/pairhmm_results_oracle.cpp)
    target_link_libraries(fastgatk-pairhmm-results-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-pairhmm-results-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    # GKL runs Float32 first, then recomputes only pairs whose scaled sum is
    # below 1e-28 in Context<double>.  Keep that native numerical boundary
    # pinned separately from ordinary GATK debug rows.
    add_executable(fastgatk-pairhmm-gkl-fallback-oracle
        tests/pairhmm_gkl_fallback_oracle.cpp)
    target_link_libraries(fastgatk-pairhmm-gkl-fallback-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-pairhmm-gkl-fallback-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-pairhmm-gkl-fallback-oracle
        COMMAND fastgatk-pairhmm-gkl-fallback-oracle)
    add_executable(fastgatk-fmc-model-oracle tests/fmc_model_oracle.cpp)
    target_link_libraries(fastgatk-fmc-model-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-fmc-model-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O2;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-fmc-model-oracle COMMAND fastgatk-fmc-model-oracle)
    add_test(NAME fastgatk-pairhmm-results-gatk-oracle
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-native/scripts/verify_pairhmm_results_oracle.py)
    set_tests_properties(fastgatk-pairhmm-results-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_PAIRHMM_RESULTS_ORACLE_BINARY=$<TARGET_FILE:fastgatk-pairhmm-results-oracle>")
    add_executable(fastgatk-fragment-aggregation-gatk-oracle
        tests/fragment_aggregation_gatk_oracle.cpp)
    target_link_libraries(fastgatk-fragment-aggregation-gatk-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-fragment-aggregation-gatk-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-fragment-aggregation-gatk-oracle
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-native/scripts/verify_fragment_aggregation_gatk_oracle.py)
    set_tests_properties(fastgatk-fragment-aggregation-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_FRAGMENT_AGGREGATION_ORACLE_BINARY=$<TARGET_FILE:fastgatk-fragment-aggregation-gatk-oracle>")
    add_executable(fastgatk-kernels-benchmark tests/benchmark.cpp)
    target_link_libraries(fastgatk-kernels-benchmark PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-kernels-benchmark PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-O3;-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    # Keep a small machine-readable matrix oracle separate from the broad API
    # smoke test.  The verifier computes the same GATK
    # SomaticLikelihoodsEngine equations independently in Python and compares
    # the Kokkos result, so a release/library-math drift cannot be hidden by a
    # hard-coded C++ assertion alone.
    add_executable(fastgatk-somatic-likelihood-oracle tests/somatic_likelihood_oracle.cpp)
    target_link_libraries(fastgatk-somatic-likelihood-oracle PRIVATE fastgatk-kernels)
    target_compile_options(fastgatk-somatic-likelihood-oracle PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
    add_test(NAME fastgatk-somatic-likelihood-gatk-oracle
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-native/scripts/verify_somatic_likelihood_oracle.py)
    set_tests_properties(fastgatk-somatic-likelihood-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SOMATIC_ORACLE_BINARY=$<TARGET_FILE:fastgatk-somatic-likelihood-oracle>")
    # NLOD is the fixed diploid hom-ref-vs-het likelihood ratio from
    # SomaticGenotypingEngine, not the variational evidence used for TLOD.
    # Keep a dedicated pinned equation oracle so emission and VCF annotation
    # cannot regress together while preserving the same sign.
    add_test(NAME fastgatk-somatic-normal-lod-gatk-oracle
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-native/scripts/verify_somatic_normal_lod_oracle.py)
    set_tests_properties(fastgatk-somatic-normal-lod-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SOMATIC_ORACLE_BINARY=$<TARGET_FILE:fastgatk-somatic-likelihood-oracle>")
    add_test(NAME fastgatk-somatic-posterior-normal-count-gatk-oracle
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../fastgatk-native/scripts/verify_somatic_posterior_normal_oracle.py)
    set_tests_properties(fastgatk-somatic-posterior-normal-count-gatk-oracle PROPERTIES
        ENVIRONMENT "FASTGATK_SOMATIC_ORACLE_BINARY=$<TARGET_FILE:fastgatk-somatic-likelihood-oracle>")
endif()
