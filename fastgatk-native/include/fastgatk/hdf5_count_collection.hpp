#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace fastgatk::hdf5 {

struct CountInterval {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

struct CountCollection {
    std::string sample = "UNKNOWN";
    std::string sequence_dictionary;
    std::vector<CountInterval> intervals;
    std::vector<double> counts;
};

// The GATK CNV panel-of-normals is a separate HDF5 container from a
// SimpleCountCollection.  Keep the representation flat at this boundary so
// the native Host can validate it before handing dense matrices to Kokkos.
struct PanelOfNormals {
    double version = 7.0;
    std::string command_line;
    std::string sequence_dictionary;
    std::vector<std::vector<double>> original_read_counts;  // samples x intervals
    std::vector<std::string> original_sample_filenames;
    std::vector<CountInterval> original_intervals;
    std::vector<double> original_interval_gc_content;
    std::vector<std::string> panel_sample_filenames;
    std::vector<CountInterval> panel_intervals;
    std::vector<double> panel_interval_fractional_medians;
    std::vector<double> singular_values;
    std::vector<std::vector<double>> eigensample_vectors;  // intervals x eigensamples
    // Maximum number of matrix cells per persisted HDF5 chunk.  This is an
    // I/O/layout control only; it never changes the numerical SVD result.
    // A value of zero selects the default single-chunk layout.
    std::size_t maximum_chunk_size = 0;
};

bool available();
CountCollection read_simple_count_collection(const std::string& path);
void write_simple_count_collection(const std::string& path, const CountCollection& collection);
PanelOfNormals read_panel_of_normals(const std::string& path);
void write_panel_of_normals(const std::string& path, const PanelOfNormals& panel);

}  // namespace fastgatk::hdf5
