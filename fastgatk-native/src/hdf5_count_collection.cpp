#include "fastgatk/hdf5_count_collection.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>

#if FASTGATK_HAS_HDF5
#include <hdf5.h>
#endif

namespace fastgatk::hdf5 {

#if FASTGATK_HAS_HDF5
namespace {

std::set<std::string> sequence_dictionary_contigs(const std::string& text) {
    // HDF5SimpleCountCollection delegates this field to SAMTextHeaderCodec.
    // Keep the native boundary deliberately narrow: only the @HD/@SQ subset
    // is a sequence dictionary.  In particular, accepting a whole BAM header
    // here would make an otherwise plausible file unreadable by Java GATK.
    std::set<std::string> contigs;
    bool saw_header = false;
    bool saw_sequence = false;
    std::size_t begin = 0;
    while (begin < text.size()) {
        const auto newline = text.find('\n', begin);
        const auto limit = newline == std::string::npos ? text.size() : newline;
        std::string line = text.substr(begin, limit - begin);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) {
            if (line.rfind("@HD\t", 0) == 0) {
                saw_header = true;
            } else if (line.rfind("@SQ\t", 0) == 0) {
                std::string contig;
                std::string length;
                std::size_t field_begin = 4;
                while (field_begin <= line.size()) {
                    const auto field_end = line.find('\t', field_begin);
                    const auto field_limit = field_end == std::string::npos ? line.size() : field_end;
                    const auto field = line.substr(field_begin, field_limit - field_begin);
                    if (field.rfind("SN:", 0) == 0) contig = field.substr(3);
                    if (field.rfind("LN:", 0) == 0) length = field.substr(3);
                    if (field_end == std::string::npos) break;
                    field_begin = field_end + 1;
                }
                if (contig.empty() || length.empty() ||
                    length.find_first_not_of("0123456789") != std::string::npos ||
                    length == "0" || !contigs.insert(contig).second)
                    throw std::runtime_error("BAD_INPUT: invalid HDF5 sequence dictionary @SQ record");
                saw_sequence = true;
            } else {
                throw std::runtime_error("BAD_INPUT: HDF5 sequence dictionary contains a non-SAM dictionary record");
            }
        }
        if (newline == std::string::npos) break;
        begin = newline + 1;
    }
    if (!saw_header || !saw_sequence)
        throw std::runtime_error("BAD_INPUT: HDF5 sequence dictionary is missing @HD or @SQ metadata");
    return contigs;
}

void check(herr_t status, const char* message) {
    if (status < 0) throw std::runtime_error(std::string("BAD_INPUT: HDF5 ") + message);
}

void ensure_group(hid_t file, const char* path) {
    hid_t group = -1;
    H5E_BEGIN_TRY {
        group = H5Gopen2(file, path, H5P_DEFAULT);
    } H5E_END_TRY;
    if (group >= 0) {
        H5Gclose(group);
        return;
    }
    group = H5Gcreate2(file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 group: ") + path);
    H5Gclose(group);
}

std::string read_string_array(hid_t file, const char* path) {
    const hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error(std::string("BAD_INPUT: missing HDF5 string dataset: ") + path);
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[1] = {0};
    check(H5Sget_simple_extent_dims(space, dims, nullptr), "reading string dimensions");
    if (dims[0] == 0) {
        H5Sclose(space); H5Dclose(dataset);
        throw std::runtime_error(std::string("BAD_INPUT: empty HDF5 string dataset: ") + path);
    }
    const hid_t type = H5Dget_type(dataset);
    char** values = new char*[dims[0]]{};
    check(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values), "reading string dataset");
    const std::string result = values[0] == nullptr ? std::string{} : std::string(values[0]);
    H5Dvlen_reclaim(type, space, H5P_DEFAULT, values);
    delete[] values;
    H5Tclose(type); H5Sclose(space); H5Dclose(dataset);
    return result;
}

std::vector<std::string> read_string_values(hid_t file, const char* path) {
    const hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error(std::string("BAD_INPUT: missing HDF5 string dataset: ") + path);
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[1] = {0};
    check(H5Sget_simple_extent_dims(space, dims, nullptr), "reading string dimensions");
    const hid_t type = H5Dget_type(dataset);
    char** values = new char*[dims[0]]{};
    check(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values), "reading string dataset");
    std::vector<std::string> result;
    result.reserve(dims[0]);
    for (hsize_t index = 0; index < dims[0]; ++index) result.emplace_back(values[index] == nullptr ? "" : values[index]);
    H5Dvlen_reclaim(type, space, H5P_DEFAULT, values);
    delete[] values;
    H5Tclose(type); H5Sclose(space); H5Dclose(dataset);
    return result;
}

std::vector<double> read_matrix(hid_t file, const char* path, hsize_t& rows, hsize_t& columns) {
    const hid_t dataset = H5Dopen2(file, path, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error(std::string("BAD_INPUT: missing HDF5 matrix: ") + path);
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[2] = {0, 0};
    if (H5Sget_simple_extent_ndims(space) != 2 || H5Sget_simple_extent_dims(space, dims, nullptr) < 0)
        throw std::runtime_error(std::string("BAD_INPUT: HDF5 matrix is not rank-2: ") + path);
    rows = dims[0]; columns = dims[1];
    std::vector<double> values(static_cast<std::size_t>(rows * columns));
    check(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), "reading matrix");
    H5Sclose(space); H5Dclose(dataset);
    return values;
}

bool link_exists(hid_t file, const std::string& path) {
    htri_t result = -1;
    H5E_BEGIN_TRY {
        result = H5Lexists(file, path.c_str(), H5P_DEFAULT);
    } H5E_END_TRY;
    return result > 0;
}

double read_scalar(hid_t file, const std::string& path) {
    const hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error("BAD_INPUT: missing HDF5 scalar: " + path);
    double value = 0.0;
    check(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value),
          "reading scalar");
    H5Dclose(dataset);
    return value;
}

std::vector<double> read_vector(hid_t file, const std::string& path) {
    const hid_t dataset = H5Dopen2(file, path.c_str(), H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error("BAD_INPUT: missing HDF5 vector: " + path);
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[1] = {0};
    if (H5Sget_simple_extent_ndims(space) != 1 || H5Sget_simple_extent_dims(space, dims, nullptr) < 0)
        throw std::runtime_error("BAD_INPUT: HDF5 vector is not rank-1: " + path);
    std::vector<double> values(static_cast<std::size_t>(dims[0]));
    if (!values.empty()) check(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
                                "reading vector");
    H5Sclose(space); H5Dclose(dataset);
    return values;
}

std::vector<double> read_chunked_matrix(hid_t file, const std::string& path,
                                        hsize_t& rows, hsize_t& columns) {
    if (!link_exists(file, path + "/num_rows")) {
        return read_matrix(file, path.c_str(), rows, columns);
    }
    rows = static_cast<hsize_t>(read_scalar(file, path + "/num_rows"));
    columns = static_cast<hsize_t>(read_scalar(file, path + "/num_columns"));
    const auto chunks = static_cast<hsize_t>(read_scalar(file, path + "/num_chunks"));
    if (rows == 0 || columns == 0 || chunks == 0)
        throw std::runtime_error("BAD_INPUT: HDF5 chunked matrix has invalid dimensions: " + path);
    std::vector<double> values(static_cast<std::size_t>(rows * columns));
    hsize_t row_offset = 0;
    for (hsize_t chunk = 0; chunk < chunks; ++chunk) {
        hsize_t chunk_rows = 0, chunk_columns = 0;
        const auto chunk_path = path + "/chunk_" + std::to_string(chunk);
        const auto chunk_values = read_matrix(file, chunk_path.c_str(),
                                              chunk_rows, chunk_columns);
        if (chunk_columns != columns || chunk_rows == 0 || row_offset + chunk_rows > rows)
            throw std::runtime_error("BAD_INPUT: invalid HDF5 matrix chunk: " + path);
        std::copy(chunk_values.begin(), chunk_values.end(), values.begin() +
                  static_cast<std::ptrdiff_t>(row_offset * columns));
        row_offset += chunk_rows;
    }
    if (row_offset != rows)
        throw std::runtime_error("BAD_INPUT: HDF5 matrix chunks do not cover matrix: " + path);
    return values;
}

void write_string_array(hid_t file, const char* path, const std::vector<std::string>& values) {
    hsize_t dims[1] = {values.size()};
    const hid_t space = H5Screate_simple(1, dims, nullptr);
    const hid_t type = H5Tcopy(H5T_C_S1);
    check(H5Tset_size(type, H5T_VARIABLE), "creating variable-length string type");
    const hid_t dataset = H5Dcreate2(file, path, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 string dataset: ") + path);
    std::vector<const char*> pointers;
    pointers.reserve(values.size());
    for (const auto& value : values) pointers.push_back(value.c_str());
    check(H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, pointers.data()), "writing string dataset");
    H5Dclose(dataset); H5Tclose(type); H5Sclose(space);
}

void write_matrix(hid_t file, const char* path, hsize_t rows, hsize_t columns,
                  const std::vector<double>& values) {
    hsize_t dims[2] = {rows, columns};
    const hid_t space = H5Screate_simple(2, dims, nullptr);
    const hid_t dataset = H5Dcreate2(file, path, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error(std::string("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 matrix: ") + path);
    check(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), "writing matrix");
    H5Dclose(dataset); H5Sclose(space);
}

void write_scalar(hid_t file, const std::string& path, double value) {
    // HDF5File.makeDouble (the Java binding used by GATK) stores a one-element
    // rank-1 array, not an HDF5 scalar dataspace.  Keep that distinction here;
    // Java HDF5File.readDouble rejects a true 0-D scalar.
    hsize_t dims[1] = {1};
    const hid_t space = H5Screate_simple(1, dims, nullptr);
    const hid_t dataset = H5Dcreate2(file, path.c_str(), H5T_IEEE_F64LE, space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 scalar: " + path);
    check(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value),
          "writing scalar");
    H5Dclose(dataset); H5Sclose(space);
}

void write_vector(hid_t file, const std::string& path, const std::vector<double>& values) {
    hsize_t dims[1] = {values.size()};
    const hid_t space = H5Screate_simple(1, dims, nullptr);
    const hid_t dataset = H5Dcreate2(file, path.c_str(), H5T_IEEE_F64LE, space,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 vector: " + path);
    if (!values.empty()) check(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                        H5P_DEFAULT, values.data()), "writing vector");
    H5Dclose(dataset); H5Sclose(space);
}

void write_intervals(hid_t file, const std::string& base,
                     const std::vector<CountInterval>& intervals) {
    if (intervals.empty()) throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 interval list is empty");
    ensure_group(file, base.c_str());
    std::vector<std::string> contigs;
    std::vector<double> matrix(3 * intervals.size());
    for (std::size_t index = 0; index < intervals.size(); ++index) {
        const auto& interval = intervals[index];
        if (interval.contig.empty() || interval.start < 1 || interval.end < interval.start)
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: invalid HDF5 interval");
        const auto found = std::find(contigs.begin(), contigs.end(), interval.contig);
        const double contig_index = found == contigs.end()
            ? static_cast<double>(contigs.size())
            : static_cast<double>(found - contigs.begin());
        if (found == contigs.end()) contigs.push_back(interval.contig);
        matrix[index] = contig_index;
        matrix[intervals.size() + index] = static_cast<double>(interval.start);
        matrix[2 * intervals.size() + index] = static_cast<double>(interval.end);
    }
    write_string_array(file, (base + "/indexed_contig_names").c_str(), contigs);
    write_matrix(file, (base + "/transposed_index_start_end").c_str(), 3, intervals.size(), matrix);
}

std::vector<CountInterval> read_intervals(hid_t file, const std::string& base) {
    const auto contigs = read_string_values(file, (base + "/indexed_contig_names").c_str());
    hsize_t rows = 0, columns = 0;
    const auto matrix = read_matrix(file, (base + "/transposed_index_start_end").c_str(), rows, columns);
    if (rows != 3 || columns == 0) throw std::runtime_error("BAD_INPUT: HDF5 interval matrix must be 3 x N");
    std::vector<CountInterval> intervals;
    intervals.reserve(columns);
    for (hsize_t index = 0; index < columns; ++index) {
        const auto contig_index = static_cast<std::size_t>(matrix[index]);
        const auto start = static_cast<std::int64_t>(matrix[columns + index]);
        const auto end = static_cast<std::int64_t>(matrix[2 * columns + index]);
        if (contig_index >= contigs.size() || start < 1 || end < start)
            throw std::runtime_error("BAD_INPUT: invalid HDF5 interval");
        intervals.push_back({contigs[contig_index], start, end});
    }
    return intervals;
}

void write_chunked_matrix(hid_t file, const std::string& path,
                          hsize_t rows, hsize_t columns,
                          const std::vector<double>& values,
                          std::size_t maximum_chunk_size = 0) {
    if (rows == 0 || columns == 0 || values.size() != rows * columns)
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: invalid HDF5 chunked matrix");
    ensure_group(file, path.c_str());
    write_scalar(file, path + "/num_rows", static_cast<double>(rows));
    write_scalar(file, path + "/num_columns", static_cast<double>(columns));
    const std::size_t cell_limit = maximum_chunk_size == 0
        ? static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns)
        : maximum_chunk_size;
    const std::size_t rows_per_chunk = std::max<std::size_t>(
        1, std::min<std::size_t>(rows, cell_limit / static_cast<std::size_t>(columns)));
    const std::size_t chunk_count = (static_cast<std::size_t>(rows) + rows_per_chunk - 1) / rows_per_chunk;
    write_scalar(file, path + "/num_chunks", static_cast<double>(chunk_count));
    for (std::size_t chunk = 0; chunk < chunk_count; ++chunk) {
        const std::size_t row_begin = chunk * rows_per_chunk;
        const std::size_t chunk_rows = std::min<std::size_t>(rows, row_begin + rows_per_chunk) - row_begin;
        std::vector<double> chunk_values(chunk_rows * static_cast<std::size_t>(columns));
        std::copy(values.begin() + static_cast<std::ptrdiff_t>(row_begin * static_cast<std::size_t>(columns)),
                  values.begin() + static_cast<std::ptrdiff_t>((row_begin + chunk_rows) * static_cast<std::size_t>(columns)),
                  chunk_values.begin());
        write_matrix(file, (path + "/chunk_" + std::to_string(chunk)).c_str(),
                     static_cast<hsize_t>(chunk_rows), columns, chunk_values);
    }
}

}  // namespace
#endif

bool available() {
#if FASTGATK_HAS_HDF5
    return true;
#else
    return false;
#endif
}

CountCollection read_simple_count_collection(const std::string& path) {
#if !FASTGATK_HAS_HDF5
    throw std::runtime_error("BACKEND_UNAVAILABLE: HDF5 support was not built");
#else
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("BAD_INPUT: cannot open HDF5 count collection: " + path);
    try {
        CountCollection collection;
        collection.sample = read_string_array(file, "/sample_metadata/sample_name");
        collection.sequence_dictionary = read_string_array(file, "/locatable_metadata/sequence_dictionary");
        if (collection.sample.empty())
            throw std::runtime_error("BAD_INPUT: HDF5 sample metadata is empty");
        const auto dictionary_contigs = sequence_dictionary_contigs(collection.sequence_dictionary);
        const auto contigs = read_string_values(file, "/intervals/indexed_contig_names");
        hsize_t interval_rows = 0, interval_columns = 0;
        const auto interval_matrix = read_matrix(file, "/intervals/transposed_index_start_end", interval_rows, interval_columns);
        if (interval_rows != 3 || interval_columns == 0)
            throw std::runtime_error("BAD_INPUT: HDF5 intervals matrix must be 3 x N");
        hsize_t count_rows = 0, count_columns = 0;
        collection.counts = read_matrix(file, "/counts/values", count_rows, count_columns);
        if (count_rows != 1 || count_columns != interval_columns)
            throw std::runtime_error("BAD_INPUT: HDF5 counts matrix does not match intervals");
        collection.intervals.reserve(interval_columns);
        for (hsize_t index = 0; index < interval_columns; ++index) {
            const auto raw_contig_index = interval_matrix[index];
            const auto raw_start = interval_matrix[interval_columns + index];
            const auto raw_end = interval_matrix[2 * interval_columns + index];
            // Validate before converting floating-point HDF5 cells to integer
            // coordinates; NaN/Inf-to-integer conversion is undefined in C++.
            const auto integer_cell = [](double value) {
                return std::isfinite(value) && value >= 0.0 &&
                       std::floor(value) == value &&
                       value <= static_cast<double>(std::numeric_limits<std::int64_t>::max());
            };
            if (!integer_cell(raw_contig_index) || !integer_cell(raw_start) ||
                !integer_cell(raw_end))
                throw std::runtime_error("BAD_INPUT: invalid HDF5 interval coordinate");
            const auto contig_index = static_cast<std::size_t>(raw_contig_index);
            const auto start = static_cast<std::int64_t>(raw_start);
            const auto end = static_cast<std::int64_t>(raw_end);
            if (contig_index >= contigs.size() || start < 1 || end < start ||
                dictionary_contigs.find(contigs[contig_index]) == dictionary_contigs.end() ||
                !std::isfinite(collection.counts[index]) || collection.counts[index] < 0.0)
                throw std::runtime_error("BAD_INPUT: invalid HDF5 count collection interval/value");
            collection.intervals.push_back({contigs[contig_index], start, end});
        }
        H5Fclose(file);
        return collection;
    } catch (...) {
        H5Fclose(file);
        throw;
    }
#endif
}

void write_simple_count_collection(const std::string& path, const CountCollection& collection) {
#if !FASTGATK_HAS_HDF5
    (void)collection;
    throw std::runtime_error("BACKEND_UNAVAILABLE: HDF5 support was not built");
#else
    if (collection.sample.empty() || collection.sequence_dictionary.empty() ||
        collection.intervals.empty() || collection.intervals.size() != collection.counts.size())
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 count collection dimensions are invalid");
    const auto dictionary_contigs = sequence_dictionary_contigs(collection.sequence_dictionary);
    std::set<std::tuple<std::string, std::int64_t, std::int64_t>> seen_intervals;
    for (std::size_t index = 0; index < collection.intervals.size(); ++index) {
        const auto& interval = collection.intervals[index];
        if (interval.contig.empty() || interval.start < 1 || interval.end < interval.start ||
            dictionary_contigs.find(interval.contig) == dictionary_contigs.end() ||
            !seen_intervals.emplace(interval.contig, interval.start, interval.end).second ||
            !std::isfinite(collection.counts[index]) || collection.counts[index] < 0.0)
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: invalid HDF5 count collection metadata");
    }
    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 count collection: " + path);
    try {
        ensure_group(file, "/sample_metadata");
        ensure_group(file, "/locatable_metadata");
        ensure_group(file, "/intervals");
        ensure_group(file, "/counts");
        write_string_array(file, "/sample_metadata/sample_name", {collection.sample});
        write_string_array(file, "/locatable_metadata/sequence_dictionary", {collection.sequence_dictionary});
        std::vector<std::string> contigs;
        std::vector<double> matrix(3 * collection.intervals.size());
        for (std::size_t index = 0; index < collection.intervals.size(); ++index) {
            const auto& interval = collection.intervals[index];
            const auto found = std::find(contigs.begin(), contigs.end(), interval.contig);
            const double contig_index = found == contigs.end()
                ? static_cast<double>(contigs.size()) : static_cast<double>(found - contigs.begin());
            if (found == contigs.end()) contigs.push_back(interval.contig);
            matrix[index] = contig_index;
            matrix[collection.intervals.size() + index] = static_cast<double>(interval.start);
            matrix[2 * collection.intervals.size() + index] = static_cast<double>(interval.end);
        }
        write_string_array(file, "/intervals/indexed_contig_names", contigs);
        write_matrix(file, "/intervals/transposed_index_start_end", 3, collection.intervals.size(), matrix);
        write_matrix(file, "/counts/values", 1, collection.counts.size(), collection.counts);
        check(H5Fclose(file), "closing HDF5 count collection");
    } catch (...) {
        H5Fclose(file);
        throw;
    }
    if (!std::filesystem::is_regular_file(path) || std::filesystem::file_size(path) == 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete HDF5 count collection: " + path);
#endif
}

PanelOfNormals read_panel_of_normals(const std::string& path) {
#if !FASTGATK_HAS_HDF5
    throw std::runtime_error("BACKEND_UNAVAILABLE: HDF5 support was not built");
#else
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("BAD_INPUT: cannot open HDF5 panel of normals: " + path);
    try {
        PanelOfNormals panel;
        if (!link_exists(file, "/version/value"))
            throw std::runtime_error("BAD_INPUT: HDF5 file is not a GATK SVD panel of normals");
        panel.version = read_scalar(file, "/version/value");
        if (!std::isfinite(panel.version) || panel.version < 7.0)
            throw std::runtime_error("BAD_INPUT: HDF5 panel of normals version is older than 7.0");
        panel.command_line = read_string_array(file, "/command_line/value");
        panel.sequence_dictionary = read_string_array(file, "/sequence_dictionary/value");
        hsize_t original_rows = 0, original_columns = 0;
        const auto original = read_chunked_matrix(file,
            "/original_data/read_counts_samples_by_intervals", original_rows, original_columns);
        panel.original_read_counts.resize(original_rows, std::vector<double>(original_columns));
        for (hsize_t row = 0; row < original_rows; ++row) {
            std::copy(original.begin() + static_cast<std::ptrdiff_t>(row * original_columns),
                      original.begin() + static_cast<std::ptrdiff_t>((row + 1) * original_columns),
                      panel.original_read_counts[row].begin());
        }
        panel.original_sample_filenames = read_string_values(file, "/original_data/sample_filenames");
        panel.original_intervals = read_intervals(file, "/original_data/intervals");
        if (panel.original_sample_filenames.size() != original_rows ||
            panel.original_intervals.size() != original_columns)
            throw std::runtime_error("BAD_INPUT: HDF5 panel original matrix metadata is misaligned");
        if (link_exists(file, "/original_data/interval_gc_content"))
            panel.original_interval_gc_content = read_vector(file, "/original_data/interval_gc_content");
        panel.panel_sample_filenames = read_string_values(file, "/panel/sample_filenames");
        panel.panel_intervals = read_intervals(file, "/panel/intervals");
        panel.panel_interval_fractional_medians = read_vector(file, "/panel/interval_fractional_medians");
        if (panel.panel_intervals.size() != panel.panel_interval_fractional_medians.size())
            throw std::runtime_error("BAD_INPUT: HDF5 panel interval medians are misaligned");
        if (link_exists(file, "/panel/singular_values"))
            panel.singular_values = read_vector(file, "/panel/singular_values");
        hsize_t eigen_rows = 0, eigen_columns = 0;
        if (link_exists(file, "/panel/transposed_eigensamples_samples_by_intervals/num_rows")) {
            const auto transposed = read_chunked_matrix(file,
                "/panel/transposed_eigensamples_samples_by_intervals", eigen_rows, eigen_columns);
            if (eigen_columns != panel.panel_intervals.size())
                throw std::runtime_error("BAD_INPUT: HDF5 panel eigensample dimensions are misaligned");
            panel.eigensample_vectors.assign(eigen_columns, std::vector<double>(eigen_rows));
            for (hsize_t eig = 0; eig < eigen_rows; ++eig)
                for (hsize_t interval = 0; interval < eigen_columns; ++interval)
                    panel.eigensample_vectors[interval][eig] = transposed[eig * eigen_columns + interval];
        }
        H5Fclose(file);
        return panel;
    } catch (...) {
        H5Fclose(file);
        throw;
    }
#endif
}

void write_panel_of_normals(const std::string& path, const PanelOfNormals& panel) {
#if !FASTGATK_HAS_HDF5
    (void)panel;
    throw std::runtime_error("BACKEND_UNAVAILABLE: HDF5 support was not built");
#else
    if (panel.original_read_counts.empty() || panel.original_sample_filenames.empty() ||
        panel.original_intervals.empty() || panel.panel_sample_filenames.empty() ||
        panel.panel_intervals.empty())
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel metadata is empty");
    const std::size_t original_rows = panel.original_read_counts.size();
    const std::size_t original_columns = panel.original_read_counts.front().size();
    if (panel.original_sample_filenames.size() != original_rows ||
        panel.original_intervals.size() != original_columns)
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel original dimensions are invalid");
    if (!panel.original_interval_gc_content.empty() &&
        panel.original_interval_gc_content.size() != original_columns)
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel GC annotation dimension is invalid");
    for (const auto gc : panel.original_interval_gc_content)
        if (!(std::isnan(gc) || (std::isfinite(gc) && gc >= 0.0 && gc <= 1.0)))
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel GC annotation is invalid");
    for (const auto& row : panel.original_read_counts)
        if (row.size() != original_columns)
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel matrix is ragged");
    if (panel.panel_intervals.size() != panel.panel_interval_fractional_medians.size())
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 panel interval medians are invalid");
    if (!panel.eigensample_vectors.empty()) {
        if (panel.eigensample_vectors.size() != panel.panel_intervals.size())
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 eigensample interval dimension is invalid");
        const auto eigensamples = panel.eigensample_vectors.front().size();
        for (const auto& row : panel.eigensample_vectors)
            if (row.size() != eigensamples)
                throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 eigensample matrix is ragged");
        if (!panel.singular_values.empty() && panel.singular_values.size() != eigensamples)
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: HDF5 singular-value dimension is invalid");
    }
    const hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create HDF5 panel of normals: " + path);
    try {
        ensure_group(file, "/version");
        ensure_group(file, "/command_line");
        ensure_group(file, "/sequence_dictionary");
        ensure_group(file, "/original_data");
        ensure_group(file, "/panel");
        write_scalar(file, "/version/value", panel.version);
        write_string_array(file, "/command_line/value", {panel.command_line});
        write_string_array(file, "/sequence_dictionary/value", {panel.sequence_dictionary});
        std::vector<double> original_values(original_rows * original_columns);
        for (std::size_t row = 0; row < original_rows; ++row)
            std::copy(panel.original_read_counts[row].begin(), panel.original_read_counts[row].end(),
                      original_values.begin() + static_cast<std::ptrdiff_t>(row * original_columns));
        write_chunked_matrix(file, "/original_data/read_counts_samples_by_intervals",
                             original_rows, original_columns, original_values,
                             panel.maximum_chunk_size);
        write_string_array(file, "/original_data/sample_filenames", panel.original_sample_filenames);
        write_intervals(file, "/original_data/intervals", panel.original_intervals);
        if (!panel.original_interval_gc_content.empty())
            write_vector(file, "/original_data/interval_gc_content", panel.original_interval_gc_content);
        write_string_array(file, "/panel/sample_filenames", panel.panel_sample_filenames);
        write_intervals(file, "/panel/intervals", panel.panel_intervals);
        write_vector(file, "/panel/interval_fractional_medians", panel.panel_interval_fractional_medians);
        if (!panel.singular_values.empty())
            write_vector(file, "/panel/singular_values", panel.singular_values);
        if (!panel.eigensample_vectors.empty()) {
            const std::size_t eigen_rows = panel.eigensample_vectors.front().size();
            const std::size_t eigen_columns = panel.eigensample_vectors.size();
            std::vector<double> transposed(eigen_rows * eigen_columns);
            for (std::size_t eig = 0; eig < eigen_rows; ++eig)
                for (std::size_t interval = 0; interval < eigen_columns; ++interval)
                    transposed[eig * eigen_columns + interval] = panel.eigensample_vectors[interval][eig];
            write_chunked_matrix(file, "/panel/transposed_eigensamples_samples_by_intervals",
                                 eigen_rows, eigen_columns, transposed,
                                 panel.maximum_chunk_size);
        }
        check(H5Fclose(file), "closing HDF5 panel of normals");
    } catch (...) {
        H5Fclose(file);
        throw;
    }
    if (!std::filesystem::is_regular_file(path) || std::filesystem::file_size(path) == 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete HDF5 panel of normals: " + path);
#endif
}

}  // namespace fastgatk::hdf5
