add_library(fastgatk-runtime STATIC src/resource.cpp src/output.cpp)
target_include_directories(fastgatk-runtime PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(fastgatk-runtime PUBLIC cxx_std_20)
target_compile_options(fastgatk-runtime PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)

# Keep the runtime independent of Kokkos headers while recording the same
# configured execution-space provenance as the kernels. A launcher allocation
# is usable only in explicit CUDA/HIP/SYCL builds; host binaries retain it as
# rejected provenance.
set(_fastgatk_runtime_execution_space "Host")
if (DEFINED FASTGATK_EXPECTED_EXECUTION_SPACE)
    set(_fastgatk_runtime_execution_space "${FASTGATK_EXPECTED_EXECUTION_SPACE}")
endif()
set(_fastgatk_runtime_device_backend 0)
if (_fastgatk_runtime_execution_space STREQUAL "Cuda" OR
    _fastgatk_runtime_execution_space STREQUAL "HIP" OR
    _fastgatk_runtime_execution_space STREQUAL "SYCL")
    set(_fastgatk_runtime_device_backend 1)
endif()
target_compile_definitions(fastgatk-runtime PRIVATE
    FASTGATK_RUNTIME_EXECUTION_SPACE=\"${_fastgatk_runtime_execution_space}\"
    FASTGATK_RUNTIME_DEVICE_BACKEND=${_fastgatk_runtime_device_backend})

option(FASTGATK_RUNTIME_BUILD_TESTS "Build runtime resource probe and pipeline smoke tests" OFF)
if (FASTGATK_RUNTIME_BUILD_TESTS)
    add_executable(fastgatk-runtime-smoke tests/resource_smoke.cpp)
    target_link_libraries(fastgatk-runtime-smoke PRIVATE fastgatk-runtime)
    target_compile_features(fastgatk-runtime-smoke PRIVATE cxx_std_20)
    target_compile_options(fastgatk-runtime-smoke PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
    add_test(NAME fastgatk-runtime-smoke COMMAND fastgatk-runtime-smoke)
    add_executable(fastgatk-runtime-pipeline-smoke tests/pipeline_smoke.cpp)
    target_link_libraries(fastgatk-runtime-pipeline-smoke PRIVATE fastgatk-runtime)
    target_compile_features(fastgatk-runtime-pipeline-smoke PRIVATE cxx_std_20)
    target_compile_options(fastgatk-runtime-pipeline-smoke PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
    add_test(NAME fastgatk-runtime-pipeline-smoke COMMAND fastgatk-runtime-pipeline-smoke)
    add_executable(fastgatk-runtime-output-smoke tests/output_smoke.cpp)
    target_link_libraries(fastgatk-runtime-output-smoke PRIVATE fastgatk-runtime)
    target_compile_features(fastgatk-runtime-output-smoke PRIVATE cxx_std_20)
    target_compile_options(fastgatk-runtime-output-smoke PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
    add_test(NAME fastgatk-runtime-output-smoke COMMAND fastgatk-runtime-output-smoke)
endif()
