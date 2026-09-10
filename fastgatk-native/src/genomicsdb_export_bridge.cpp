// Optional bridge to the native GenomicsDB query API shipped in a GATK
// distribution.  This translation unit deliberately uses the legacy GNU C++
// ABI because GATK 4.6 packages libtiledbgenomicsdb.so with
// _GLIBCXX_USE_CXX11_ABI=0.  It is an isolated process boundary: the main
// Kokkos/HTSlib binaries never link to this ABI or load the library in-process.

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI != 0
#error "GenomicsDB bridge must be compiled with _GLIBCXX_USE_CXX11_ABI=0"
#endif
#ifndef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
#endif

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/types.h>
#include <unistd.h>
#endif

// Keep this declaration intentionally small, but preserve the ABI-visible
// constructor and member signatures.  The implementation's private state is
// owned by the shared library; the opaque storage is large enough for the
// current GenomicsDB query object and avoids vendoring the entire upstream
// C++ header into the public fastgatk include tree.
class GenomicsDB {
 public:
    GenomicsDB(const std::string& workspace,
               const std::string& callset_mapping_file,
               const std::string& vid_mapping_file,
               const std::vector<std::string> attributes = {},
               std::uint64_t segment_size = 10u * 1024u * 1024u);
    ~GenomicsDB();

    void generate_vcf(const std::string& array,
                      std::vector<std::pair<std::int64_t, std::int64_t>> column_ranges,
                      std::vector<std::pair<std::int64_t, std::int64_t>> row_ranges,
                      const std::string& reference_genome,
                      const std::string& vcf_header,
                      const std::string& output,
                      const std::string& output_format = "",
                      bool overwrite = false);

 private:
    // The current query object is only a few hundred bytes, but leave ample
    // room for a matching GATK release to add private bookkeeping without
    // letting the constructor write past the placement buffer. The bridge is
    // version-coupled and still validates the generated VCF contract.
    alignas(std::max_align_t) unsigned char opaque_state_[65536]{};
};

namespace {

std::string read_text(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("BAD_INPUT: cannot open GenomicsDB mapping file: " + path.string());
    std::ostringstream text;
    text << stream.rdbuf();
    if (!stream.good() && !stream.eof())
        throw std::runtime_error("BAD_INPUT: cannot read GenomicsDB mapping file: " + path.string());
    return text.str();
}

void write_text(const std::filesystem::path& path, const std::string& text) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write GenomicsDB bridge file: " + path.string());
    stream << text;
    if (!stream)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize GenomicsDB bridge file: " + path.string());
}

// GATK's JSON writer emits several numeric mapping values as JSON strings,
// while the standalone GenomicsDB C++ query API expects JSON numbers.  Only
// normalize the four schema keys known to be numeric; descriptors such as
// variable_length_descriptor="R" and contig names remain untouched.
std::string normalize_numeric_key(const std::string& input, const std::string& key) {
    const std::string marker = "\"" + key + "\"";
    std::string output = input;
    std::size_t search = 0;
    while ((search = output.find(marker, search)) != std::string::npos) {
        std::size_t cursor = search + marker.size();
        while (cursor < output.size() && std::isspace(static_cast<unsigned char>(output[cursor]))) ++cursor;
        if (cursor >= output.size() || output[cursor] != ':') {
            search += marker.size();
            continue;
        }
        ++cursor;
        while (cursor < output.size() && std::isspace(static_cast<unsigned char>(output[cursor]))) ++cursor;
        if (cursor >= output.size() || output[cursor] != '"') {
            search = cursor;
            continue;
        }
        const std::size_t quote = cursor++;
        const std::size_t digits = cursor;
        while (cursor < output.size() && std::isdigit(static_cast<unsigned char>(output[cursor]))) ++cursor;
        if (cursor == digits || cursor >= output.size() || output[cursor] != '"') {
            search = cursor;
            continue;
        }
        const auto numeric = output.substr(digits, cursor - digits);
        output.erase(quote, (cursor + 1) - quote);
        output.insert(quote, numeric);
        search = quote + numeric.size();
    }
    return output;
}

std::string normalize_callset(const std::string& input) {
    auto output = normalize_numeric_key(input, "row_idx");
    return normalize_numeric_key(output, "idx_in_file");
}

std::string normalize_vidmap(const std::string& input) {
    auto output = normalize_numeric_key(input, "length");
    return normalize_numeric_key(output, "tiledb_column_offset");
}

std::string object_string_field(const std::string& object, const std::string& key) {
    const std::string marker = "\"" + key + "\"";
    const auto begin = object.find(marker);
    if (begin == std::string::npos) return {};
    auto cursor = object.find(':', begin + marker.size());
    if (cursor == std::string::npos) return {};
    ++cursor;
    while (cursor < object.size() && std::isspace(static_cast<unsigned char>(object[cursor]))) ++cursor;
    if (cursor >= object.size()) return {};
    if (object[cursor] == '"') {
        ++cursor;
        const auto end = object.find('"', cursor);
        return end == std::string::npos ? std::string{} : object.substr(cursor, end - cursor);
    }
    const auto end = object.find_first_of(",}\n\r\t ", cursor);
    return object.substr(cursor, end == std::string::npos ? std::string::npos : end - cursor);
}

std::uint64_t object_uint_field(const std::string& object, const std::string& key) {
    const auto value = object_string_field(object, key);
    if (value.empty()) throw std::runtime_error("BAD_INPUT: GenomicsDB VID mapping is missing " + key);
    std::size_t consumed = 0;
    const auto parsed = std::stoull(value, &consumed);
    if (consumed != value.size())
        throw std::runtime_error("BAD_INPUT: GenomicsDB VID mapping has invalid " + key);
    return parsed;
}

struct ContigColumnInfo {
    std::uint64_t offset = 0;
    std::uint64_t length = 0;
};

std::map<std::string, ContigColumnInfo> read_contig_columns(const std::string& vidmap) {
    std::map<std::string, ContigColumnInfo> result;
    const auto contigs = vidmap.find("\"contigs\"");
    if (contigs == std::string::npos) return result;
    auto cursor = vidmap.find('{', contigs);
    while (cursor != std::string::npos) {
        const auto end = vidmap.find('}', cursor);
        if (end == std::string::npos) break;
        const auto object = vidmap.substr(cursor, end - cursor + 1);
        const auto name = object_string_field(object, "name");
        if (!name.empty()) {
            result.emplace(name, ContigColumnInfo{object_uint_field(object, "tiledb_column_offset"),
                                                  object_uint_field(object, "length")});
        }
        const auto next = vidmap.find('{', end + 1);
        const auto array_end = vidmap.find(']', end + 1);
        if (next == std::string::npos || (array_end != std::string::npos && next > array_end)) break;
        cursor = next;
    }
    return result;
}

struct RequestedInterval {
    std::string contig;
    std::uint64_t begin = 0;
    std::uint64_t end = 0;
};

std::optional<RequestedInterval> parse_requested_interval(const std::string& text) {
    const auto colon = text.find(':');
    if (colon == std::string::npos || colon == 0) return std::nullopt;
    const auto dash = text.find('-', colon + 1);
    if (dash == std::string::npos || dash == colon + 1 || dash + 1 >= text.size()) return std::nullopt;
    try {
        std::size_t consumed = 0;
        const auto begin = std::stoull(text.substr(colon + 1, dash - colon - 1), &consumed);
        if (consumed != dash - colon - 1 || begin == 0) return std::nullopt;
        consumed = 0;
        const auto end = std::stoull(text.substr(dash + 1), &consumed);
        if (consumed != text.size() - dash - 1 || end < begin) return std::nullopt;
        return RequestedInterval{text.substr(0, colon), begin, end};
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::vector<std::pair<std::int64_t, std::int64_t>> requested_column_ranges(
    const std::string& vidmap, const std::vector<std::string>& intervals) {
    const auto columns = read_contig_columns(vidmap);
    std::vector<std::pair<std::int64_t, std::int64_t>> result;
    for (const auto& text : intervals) {
        const auto parsed = parse_requested_interval(text);
        if (!parsed.has_value()) continue;
        const auto found = columns.find(parsed->contig);
        if (found == columns.end()) continue;
        const auto end = std::min(parsed->end, found->second.length);
        if (end < parsed->begin) continue;
        const auto first = found->second.offset + parsed->begin - 1;
        const auto last = found->second.offset + end - 1;
        if (first > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) ||
            last > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()))
            throw std::runtime_error("BAD_INPUT: GenomicsDB interval column overflows int64");
        result.emplace_back(static_cast<std::int64_t>(first), static_cast<std::int64_t>(last));
    }
    if (result.empty())
        result.emplace_back(0, std::numeric_limits<std::int64_t>::max() - 1);
    std::sort(result.begin(), result.end());
    std::vector<std::pair<std::int64_t, std::int64_t>> merged;
    for (const auto range : result) {
        if (!merged.empty() && range.first <= merged.back().second + 1)
            merged.back().second = std::max(merged.back().second, range.second);
        else
            merged.push_back(range);
    }
    return merged;
}

std::string process_id() {
#if defined(__unix__) || defined(__APPLE__)
    return std::to_string(static_cast<unsigned long long>(::getpid()));
#else
    return "0";
#endif
}

std::string unique_suffix() {
    const auto ticks = std::chrono::steady_clock::now().time_since_epoch().count();
    return process_id() + "-" + std::to_string(static_cast<long long>(ticks));
}

std::vector<std::filesystem::path> discover_arrays(const std::filesystem::path& workspace) {
    std::vector<std::filesystem::path> arrays;
    std::error_code error;
    for (const auto& entry : std::filesystem::directory_iterator(workspace, error)) {
        if (error) break;
        if (!entry.is_directory(error) || error) {
            error.clear();
            continue;
        }
        std::error_code schema_error;
        if (std::filesystem::is_regular_file(entry.path() / "__array_schema.tdb", schema_error) &&
            !schema_error)
            arrays.push_back(entry.path());
    }
    if (error)
        throw std::runtime_error("BAD_INPUT: cannot inspect GenomicsDB workspace: " + workspace.string());
    std::sort(arrays.begin(), arrays.end());
    return arrays;
}

std::filesystem::path find_required(const std::filesystem::path& workspace,
                                    std::initializer_list<const char*> names,
                                    const char* description) {
    for (const auto* name : names) {
        const auto path = workspace / name;
        std::error_code error;
        if (std::filesystem::is_regular_file(path, error) && !error)
            return path;
    }
    throw std::runtime_error(std::string("BAD_INPUT: GenomicsDB workspace is missing ") + description + ": " +
                             workspace.string());
}

void remove_if_exists(const std::filesystem::path& path) {
    std::error_code error;
    std::filesystem::remove(path, error);
}

struct TempDir {
    std::filesystem::path path;
    ~TempDir() {
        if (path.empty()) return;
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

void write_manifest(const std::filesystem::path& manifest,
                    const std::vector<std::pair<std::string, std::filesystem::path>>& outputs) {
    std::ofstream stream(manifest, std::ios::binary | std::ios::trunc);
    if (!stream)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write GenomicsDB export manifest: " +
                                 manifest.string());
    stream << "# fastgatk GenomicsDB native export manifest v1\n";
    for (const auto& [array, output] : outputs)
        stream << array << '\t' << output.string() << '\n';
    if (!stream)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize GenomicsDB export manifest: " +
                                 manifest.string());
}

int run(const std::filesystem::path& workspace, const std::filesystem::path& manifest,
        const std::vector<std::string>& intervals) {
    std::error_code error;
    if (!std::filesystem::is_directory(workspace, error) || error)
        throw std::runtime_error("BAD_INPUT: GenomicsDB workspace is not a directory: " + workspace.string());
    const auto arrays = discover_arrays(workspace);
    if (arrays.empty())
        throw std::runtime_error("BAD_INPUT: GenomicsDB workspace contains no queryable arrays: " + workspace.string());
    const auto callset = find_required(workspace, {"callset.json", "callsets.json"}, "callset mapping");
    const auto vidmap = find_required(workspace, {"vidmap.json", "vid_mapping.json"}, "VID mapping");
    const auto header = find_required(workspace, {"vcfheader.vcf", "vcf_header.vcf"}, "VCF header");

    const auto temp_root = std::filesystem::temp_directory_path() /
                           ("fastgatk-genomicsdb-bridge-" + unique_suffix());
    std::filesystem::create_directories(temp_root);
    TempDir cleanup{temp_root};
    const auto callset_normalized = temp_root / "callset.json";
    const auto vidmap_normalized = temp_root / "vidmap.json";
    write_text(callset_normalized, normalize_callset(read_text(callset)));
    const auto vidmap_text = normalize_vidmap(read_text(vidmap));
    write_text(vidmap_normalized, vidmap_text);
    const auto column_ranges = requested_column_ranges(vidmap_text, intervals);

    const auto manifest_directory = manifest.has_parent_path()
        ? manifest.parent_path() : std::filesystem::current_path();
    std::filesystem::create_directories(manifest_directory);
    std::vector<std::pair<std::string, std::filesystem::path>> outputs;
    outputs.reserve(arrays.size());
    try {
        for (std::size_t index = 0; index < arrays.size(); ++index) {
            const auto output = manifest_directory /
                (manifest.filename().string() + ".array-" + std::to_string(index) + ".vcf.gz");
            std::error_code existing_error;
            if (std::filesystem::exists(output, existing_error) && !existing_error)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: refusing to overwrite GenomicsDB bridge output: " +
                                         output.string());
            const std::string array_name = arrays[index].filename().string();
            std::vector<std::pair<std::int64_t, std::int64_t>> rows;
            void* storage = ::operator new(sizeof(GenomicsDB), std::align_val_t(alignof(GenomicsDB)));
            GenomicsDB* database = nullptr;
            try {
                database = new (storage) GenomicsDB(workspace.string(), callset_normalized.string(),
                                                    vidmap_normalized.string(), {}, 10u * 1024u * 1024u);
                // The standalone API uses "z" for BGZF output. This gives the
                // parent HTSlib reader a normal .tbi indexing boundary and
                // avoids decoding an entire workspace for a narrow -L query.
                database->generate_vcf(array_name, column_ranges, rows, "", header.string(), output.string(), "z", false);
                database->~GenomicsDB();
                ::operator delete(storage, std::align_val_t(alignof(GenomicsDB)));
                database = nullptr;
            } catch (...) {
                if (database != nullptr) database->~GenomicsDB();
                ::operator delete(storage, std::align_val_t(alignof(GenomicsDB)));
                throw;
            }
            std::error_code output_error;
            if (!std::filesystem::is_regular_file(output, output_error) || output_error ||
                std::filesystem::file_size(output, output_error) == 0 || output_error) {
                remove_if_exists(output);
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GenomicsDB export produced an empty VCF: " +
                                         output.string());
            }
            outputs.emplace_back(array_name, output);
        }
        write_manifest(manifest, outputs);
    } catch (...) {
        for (const auto& [unused, generated] : outputs) remove_if_exists(generated);
        throw;
    }
    std::cout << "{\"status\":\"pass\",\"arrays\":" << arrays.size()
              << ",\"manifest\":\"" << manifest.string() << "\"}\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
            std::cout << "fastgatk-genomicsdb-export WORKSPACE MANIFEST [--interval CONTIG:START-END ...]\n"
                         "  Export true GATK GenomicsDB arrays to temporary indexed VCF shards.\n";
            return 0;
        }
        if (argc < 3) {
            std::cerr << "usage: fastgatk-genomicsdb-export WORKSPACE MANIFEST [--interval CONTIG:START-END ...]\n";
            return 2;
        }
        std::vector<std::string> intervals;
        for (int index = 3; index < argc; ++index) {
            if (std::string(argv[index]) != "--interval" || index + 1 >= argc)
                throw std::invalid_argument("usage: --interval CONTIG:START-END");
            intervals.emplace_back(argv[++index]);
        }
        return run(std::filesystem::path(argv[1]), std::filesystem::path(argv[2]), intervals);
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-genomicsdb-export: " << error.what() << '\n';
        return 2;
    }
}
