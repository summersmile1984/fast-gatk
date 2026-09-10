#include "fastgatk/io/tribble_index.hpp"

#include "fastgatk/io/java_numeric.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <chrono>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <sys/stat.h>
#include <vector>

namespace fastgatk::io {
namespace {

struct Block { std::uint64_t start = 0; std::uint64_t end = 0; };
struct Contig {
    std::string name;
    std::vector<Block> blocks;
    int bin_width = 2000;
    int longest_feature = 0;
    int n_features = 0;
};

class TemporaryPath {
public:
    explicit TemporaryPath(std::string path) : path_(std::move(path)) {}
    TemporaryPath(const TemporaryPath&) = delete;
    TemporaryPath& operator=(const TemporaryPath&) = delete;
    ~TemporaryPath() {
        if (!committed_) {
            std::error_code error;
            std::filesystem::remove(path_, error);
        }
    }
    const std::string& path() const { return path_; }
    void commit() noexcept { committed_ = true; }

private:
    std::string path_;
    bool committed_ = false;
};

void write_i32(std::ofstream& out, std::int32_t value) {
    const auto bits = static_cast<std::uint32_t>(value);
    for (int shift = 0; shift < 32; shift += 8)
        out.put(static_cast<char>((bits >> shift) & 0xffU));
}

void write_i64(std::ofstream& out, std::int64_t value) {
    const auto bits = static_cast<std::uint64_t>(value);
    for (int shift = 0; shift < 64; shift += 8)
        out.put(static_cast<char>((bits >> shift) & 0xffU));
}

void write_string(std::ofstream& out, const std::string& value) {
    if (value.find('\0') != std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Tribble index strings cannot contain NUL");
    out.write(value.data(), static_cast<std::streamsize>(value.size()));
    out.put('\0');
}

std::string file_uri(const std::string& path) {
    const auto absolute = std::filesystem::absolute(path).generic_string();
    constexpr char hex[] = "0123456789ABCDEF";
    std::string result = "file://";
    for (const unsigned char c : absolute) {
        const bool unreserved = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.' || c == '~' || c == '/';
        if (unreserved) result.push_back(static_cast<char>(c));
        else {
            result.push_back('%');
            result.push_back(hex[(c >> 4) & 0xf]);
            result.push_back(hex[c & 0xf]);
        }
    }
    return result;
}

std::int64_t file_mtime_millis(const std::string& path) {
    struct stat status {};
    if (::stat(path.c_str(), &status) != 0)
        throw std::runtime_error("BAD_INPUT: cannot stat VCF for Tribble index: " + path);
    return static_cast<std::int64_t>(status.st_mtim.tv_sec) * 1000 +
        static_cast<std::int64_t>(status.st_mtim.tv_nsec / 1000000);
}

std::vector<std::string_view> split_tabs(std::string_view line) {
    std::vector<std::string_view> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        if (end == std::string_view::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
    return fields;
}

int positive_coordinate(std::string_view value, const char* field) {
    int result = 0;
    if (value.empty())
        throw std::runtime_error(std::string("BAD_INPUT: missing VCF ") + field);
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() || result < 1)
        throw std::runtime_error(std::string("BAD_INPUT: invalid VCF ") + field);
    return result;
}

int record_end(std::string_view info, int start, int ref_length) {
    int end = start + std::max(1, ref_length) - 1;
    std::size_t begin = 0;
    while (begin <= info.size()) {
        const auto next = info.find(';', begin);
        const auto token = info.substr(begin,
            next == std::string_view::npos ? info.size() - begin : next - begin);
        if (token.rfind("END=", 0) == 0) {
            end = positive_coordinate(token.substr(4), "END");
            break;
        }
        if (next == std::string_view::npos) break;
        begin = next + 1;
    }
    if (end < start) throw std::runtime_error("BAD_INPUT: VCF END precedes POS");
    return end;
}

void optimize(Contig& contig) {
    constexpr double max_features_per_bin = 100.0;
    constexpr int max_bin_width = 1000000000;
    constexpr int max_occupied_bin_width = 1024000;
    while (contig.blocks.size() > 1) {
        std::uint64_t total_size = 0;
        for (const auto& block : contig.blocks) total_size += block.end - block.start;
        const double average = static_cast<double>(total_size) /
            static_cast<double>(std::max(1, contig.n_features));
        if (average <= 0.0) break;
        double densest = -1.0;
        for (const auto& block : contig.blocks)
            densest = std::max(densest, static_cast<double>(block.end - block.start) / average);
        const bool bad_width = contig.bin_width > max_bin_width || contig.bin_width < 0 ||
            (contig.n_features > 1 && contig.bin_width >= max_occupied_bin_width);
        if (densest > max_features_per_bin || bad_width) break;
        std::vector<Block> merged;
        merged.reserve((contig.blocks.size() + 1) / 2);
        for (std::size_t i = 0; i < contig.blocks.size(); i += 2) {
            auto block = contig.blocks[i];
            if (i + 1 < contig.blocks.size())
                block.end += contig.blocks[i + 1].end - contig.blocks[i + 1].start;
            merged.push_back(block);
        }
        if (merged.size() == 1) break;
        contig.blocks.swap(merged);
        contig.bin_width = contig.bin_width > max_bin_width / 2
            ? max_bin_width + 1 : contig.bin_width * 2;
    }
}

}  // namespace

void write_uncompressed_vcf_tribble_index(const std::string& input_path,
                                           const std::string& index_path) {
    std::error_code error;
    const auto input_size = std::filesystem::file_size(input_path, error);
    if (error) throw std::runtime_error("BAD_INPUT: cannot determine VCF size: " + input_path);
    std::ifstream input(input_path, std::ios::binary);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot read VCF: " + input_path);

    struct Property { std::string key; std::string value; };
    std::vector<Property> dictionary;
    std::vector<Contig> contigs;
    std::map<std::string, std::size_t> seen_contigs;
    std::string current_contig;
    int last_start = 0;
    std::size_t records = 0;
    double mean = 0.0;
    double m2 = 0.0;
    std::string line;
    while (true) {
        const auto offset = input.tellg();
        if (!std::getline(input, line)) break;
        if (offset < 0) throw std::runtime_error("BAD_INPUT: invalid VCF byte offset");
        if (!line.empty() && line.back() == '\r') line.pop_back();
        constexpr std::string_view prefix = "##contig=<ID=";
        if (line.rfind(prefix, 0) == 0) {
            const auto id_begin = prefix.size();
            const auto id_end = line.find(',', id_begin);
            const auto length_tag = line.find("length=", id_end == std::string::npos ? id_begin : id_end);
            const auto length_end = line.find('>', length_tag == std::string::npos ? id_begin : length_tag);
            if (id_end != std::string::npos && length_tag != std::string::npos && length_end != std::string::npos)
                dictionary.push_back({"DICT:" + line.substr(id_begin, id_end - id_begin),
                                      line.substr(length_tag + 7, length_end - length_tag - 7)});
        }
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 5)
            throw std::runtime_error("BAD_INPUT: malformed VCF record while building Tribble index");
        const std::string contig_name(fields[0]);
        const int start = positive_coordinate(fields[1], "POS");
        const int end = record_end(fields.size() > 7 ? fields[7] : std::string_view{},
                                   start, static_cast<int>(fields[3].size()));
        if (contig_name == current_contig && start < last_start)
            throw std::runtime_error("BAD_INPUT: VCF records are not sorted by coordinate");
        if (contig_name != current_contig) {
            if (seen_contigs.find(contig_name) != seen_contigs.end())
                throw std::runtime_error("BAD_INPUT: VCF contigs are not grouped in coordinate order: " + contig_name);
            seen_contigs.emplace(contig_name, contigs.size());
            contigs.push_back(Contig{contig_name, {{static_cast<std::uint64_t>(offset), 0}}, 2000, 0, 0});
            current_contig = contig_name;
            last_start = 0;
        }
        last_start = start;
        auto& contig = contigs.back();
        while (static_cast<std::int64_t>(start) >
               static_cast<std::int64_t>(contig.blocks.size()) * contig.bin_width)
            contig.blocks.push_back({static_cast<std::uint64_t>(offset), 0});
        contig.longest_feature = std::max(contig.longest_feature, end - start + 1);
        ++contig.n_features;
        ++records;
        const double value = static_cast<double>(contig.longest_feature);
        const double delta = value - mean;
        mean += delta / static_cast<double>(records);
        m2 += delta * (value - mean);
    }
    for (std::size_t c = 0; c < contigs.size(); ++c) {
        auto& contig = contigs[c];
        const auto contig_end = c + 1 < contigs.size()
            ? contigs[c + 1].blocks.front().start : input_size;
        for (std::size_t i = 0; i < contig.blocks.size(); ++i) {
            contig.blocks[i].end = i + 1 < contig.blocks.size()
                ? contig.blocks[i + 1].start : contig_end;
            if (contig.blocks[i].end < contig.blocks[i].start)
                throw std::runtime_error("BAD_INPUT: invalid VCF block offsets");
        }
        optimize(contig);
    }
    const double variance = records <= 1 ? 0.0 : m2 / static_cast<double>(records - 1);
    const double stddev = std::sqrt(std::max(0.0, variance));
    // Publish the sidecar atomically.  A failed index build must never leave a
    // truncated `.idx` that downstream GATK/HTSlib interprets as complete.
    const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
    TemporaryPath temporary(index_path + ".tmp." + std::to_string(nonce));
    std::ofstream output(temporary.path(), std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    write_i32(output, 1480870228);  // TIDX
    write_i32(output, 1);           // LinearIndex
    write_i32(output, 3);           // AbstractIndex v3
    write_string(output, file_uri(input_path));
    write_i64(output, static_cast<std::int64_t>(input_size));
    write_i64(output, file_mtime_millis(input_path));
    write_string(output, "");
    write_i32(output, 0);
    write_i32(output, static_cast<std::int32_t>(dictionary.size() + 4));
    for (const auto& property : dictionary) {
        write_string(output, property.key);
        write_string(output, property.value);
    }
    write_string(output, "FEATURE_LENGTH_MEAN");
    write_string(output, java_double(mean));
    write_string(output, "FEATURE_LENGTH_STD_DEV");
    write_string(output, java_double(stddev));
    write_string(output, "MEAN_FEATURE_VARIANCE");
    write_string(output, java_double(variance));
    write_string(output, "FEATURE_COUNT");
    write_string(output, std::to_string(records));
    if (contigs.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: too many VCF contigs for Tribble index");
    write_i32(output, static_cast<std::int32_t>(contigs.size()));
    for (const auto& contig : contigs) {
        write_string(output, contig.name);
        write_i32(output, contig.bin_width);
        write_i32(output, static_cast<std::int32_t>(contig.blocks.size()));
        write_i32(output, contig.longest_feature);
        write_i32(output, 0);
        write_i32(output, contig.n_features);
        for (const auto& block : contig.blocks)
            write_i64(output, static_cast<std::int64_t>(block.start));
        write_i64(output, static_cast<std::int64_t>(contig.blocks.empty() ? 0 : contig.blocks.back().end));
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble index: " + index_path);
    output.close();
    std::error_code rename_error;
    std::filesystem::rename(temporary.path(), index_path, rename_error);
    if (rename_error) {
        // POSIX rename replaces an existing sidecar, while some filesystems
        // reject that operation.  Only the explicitly requested index target
        // is removed for the retry; the temporary guard still cleans up on
        // every failure path.
        std::error_code remove_error;
        std::filesystem::remove(index_path, remove_error);
        rename_error.clear();
        std::filesystem::rename(temporary.path(), index_path, rename_error);
    }
    if (rename_error)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot publish Tribble index: " + index_path);
    temporary.commit();
}

}  // namespace fastgatk::io
