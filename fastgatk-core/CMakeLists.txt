add_library(fastgatk-core STATIC
    src/plan.cpp
    src/flow_codec.cpp
    src/tribble_index.cpp)
target_include_directories(fastgatk-core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(fastgatk-core PUBLIC Kokkos::kokkos)
target_compile_features(fastgatk-core PUBLIC cxx_std_20)
target_compile_options(fastgatk-core PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-ffp-contract=off;-fno-fast-math>)
