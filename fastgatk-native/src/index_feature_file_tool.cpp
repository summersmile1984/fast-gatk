#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <cctype>
#include <charconv>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sys/stat.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string manifest;
    int threads = 1;
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string value_for(int& index, int argc, char** argv, const std::string& argument,
                     std::string_view name, const char* short_name = nullptr) {
    if (const auto value = inline_value(argument, name); !value.empty()) return value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument("missing value for " + std::string(name));
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-index-feature-file (GATK-compatible native contract-compatible)\n"
                         "  -I, --input FILE             feature file (BGZF VCF/GVCF/BED or BCF)\n"
                         "  -O, --output FILE            output .tbi/.csi/.idx (optional)\n"
                         "      --threads N              HTSlib decoder threads\n"
                         "      --output-manifest FILE   OutputManifest JSON\n"
                         "\n"
                         "Uncompressed VCF files produce an HTSJDK-compatible linear .idx;\n"
                         "other unsupported Tribble codecs fail closed and can be routed with --fallback.\n";
            std::exit(0);
        } else if (argument == "-I" || argument == "--input" ||
                   !inline_value(argument, "--input").empty()) {
            options.input = value_for(index, argc, argv, argument, "--input", "-I");
        } else if (argument == "-O" || argument == "--output" ||
                   !inline_value(argument, "--output").empty()) {
            options.output = value_for(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(value_for(index, argc, argv, argument, "--threads"));
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() ||
                   !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = value_for(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--seconds-between-progress-updates", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                     "--seconds-between-progress-updates");
                (void)value_for(index, argc, argv, argument, name);
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    if (!std::filesystem::is_regular_file(options.input))
        throw std::runtime_error("BAD_INPUT: feature file is not readable: " + options.input);
    return options;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char character : value) {
        if (character == '"' || character == '\\') out << '\\';
        if (character == '\n') out << "\\n";
        else if (character == '\r') out << "\\r";
        else if (character == '\t') out << "\\t";
        else out << character;
    }
    return out.str();
}

const char* format_name(enum htsExactFormat format) {
    switch (format) {
    case vcf: return "VCF";
    case bcf: return "BCF";
    case bed: return "BED";
    default: return "UNKNOWN";
    }
}

const char* compression_name(enum htsCompression compression) {
    switch (compression) {
    case bgzf: return "BGZF";
    case gzip: return "GZIP";
    case no_compression: return "NONE";
    default: return "OTHER";
    }
}

struct DetectedFormat {
    enum htsExactFormat format = unknown_format;
    enum htsCompression compression = no_compression;
};

DetectedFormat detect_format(const std::string& input) {
    htsFile* file = hts_open(input.c_str(), "r");
    if (!file) throw std::runtime_error("BAD_INPUT: HTSlib cannot open feature file: " + input);
    const htsFormat* detected = hts_get_format(file);
    DetectedFormat result;
    if (detected) {
        result.format = detected->format;
        result.compression = detected->compression;
    }
    if (hts_close(file) != 0)
        throw std::runtime_error("BAD_INPUT: HTSlib failed closing feature file: " + input);
    return result;
}

std::string default_index_path(const Options& options, const DetectedFormat& detected) {
    if (!options.output.empty()) return options.output;
    if (detected.format == bcf) return options.input + ".csi";
    if ((detected.format == vcf || detected.format == bed) && detected.compression == no_compression)
        return options.input + ".idx";
    return options.input + ".tbi";
}

constexpr int kTribbleMagic = 1480870228;
constexpr int kLinearIndexType = 1;
constexpr int kGvcfBinWidth = 128000;

struct LinearBlock {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct LinearContig {
    std::string name;
    std::vector<LinearBlock> blocks;
    int bin_width = kGvcfBinWidth;
    int longest_feature = 0;
    int n_features = 0;
};

struct LinearIndexStats {
    std::size_t records = 0;
    std::size_t contigs = 0;
    std::size_t blocks = 0;
    // DynamicIndexCreator publishes these values as v3 properties.  Keep the
    // summary separate from the linear layout so ordinary VCF can select the
    // same linear/interval-tree candidate as HTSJDK without carrying Java's
    // feature object graph into the native path.
    double feature_length_mean = 0.0;
    double feature_length_stddev = 0.0;
    double feature_length_variance = 0.0;
};

struct IntervalBlock {
    int start = 0;
    int end = 0;
    std::uint64_t file_start = 0;
    std::uint64_t file_end = 0;
};

struct IntervalContig {
    std::string name;
    std::vector<IntervalBlock> intervals;
};

struct FeatureSummary {
    std::size_t records = 0;
    std::uint64_t bases_seen = 0;
    int longest_feature = 0;
    // HTSJDK's RunningStat uses the usual online mean/M2 update.  The
    // DynamicIndexCreator pushes the running *maximum* feature length, not
    // the raw length, so retain that subtle behavior for property parity.
    double mean = 0.0;
    double m2 = 0.0;
    std::size_t stat_count = 0;
    std::vector<IntervalContig> interval_contigs;
};

void push_feature_stat(FeatureSummary& summary, double value) {
    ++summary.stat_count;
    const double delta = value - summary.mean;
    summary.mean += delta / static_cast<double>(summary.stat_count);
    const double delta2 = value - summary.mean;
    summary.m2 += delta * delta2;
}

double summary_variance(const FeatureSummary& summary) {
    // htsjdk.tribble.util.MathUtils.RunningStat reports the unbiased sample
    // variance (M2/(n-1)), not the population variance.  This is observable
    // in the BED fixture's FEATURE_LENGTH_STD_DEV and must be retained for
    // DynamicIndexCreator byte compatibility.
    return summary.stat_count <= 1 ? 0.0 :
        summary.m2 / static_cast<double>(summary.stat_count - 1);
}

double summary_stddev(const FeatureSummary& summary) {
    return std::sqrt(std::max(0.0, summary_variance(summary)));
}

std::string java_double_string(double value);

void optimize_linear_contig(LinearContig& contig) {
    // LinearIndex.optimize() adaptively merges adjacent bins until the
    // densest bin would contain more than MAX_FEATURES_PER_BIN (100) average
    // features, or the occupied-contig safety limit is reached.  This is
    // important for GVCFs: HTSJDK's advertised 128 kb starting bin is only a
    // starting point and the final .idx must carry the optimized width.
    constexpr double max_features_per_bin = 100.0;
    constexpr int max_bin_width = 1000000000;
    constexpr int max_occupied_bin_width = 1024000;
    while (contig.blocks.size() > 1) {
        std::uint64_t total_size = 0;
        double densest = -1.0;
        for (const auto& block : contig.blocks) {
            const std::uint64_t size = block.end - block.start;
            total_size += size;
        }
        const double average_feature_size = static_cast<double>(total_size) /
            static_cast<double>(std::max(1, contig.n_features));
        for (const auto& block : contig.blocks)
            densest = std::max(densest, static_cast<double>(block.end - block.start) /
                               average_feature_size);
        // HTSJDK's occupied-contig guard treats the configured maximum as a
        // terminal width (the next doubling is never attempted), so retain
        // the 1,024,000 bp layout observed in its serialized indices.
        const bool bad_width = contig.bin_width > max_bin_width || contig.bin_width < 0 ||
            (contig.n_features > 1 && contig.bin_width >= max_occupied_bin_width);
        if (densest > max_features_per_bin || bad_width) break;

        std::vector<LinearBlock> merged;
        merged.reserve((contig.blocks.size() + 1) / 2);
        for (std::size_t index = 0; index < contig.blocks.size(); index += 2) {
            LinearBlock block = contig.blocks[index];
            if (index + 1 < contig.blocks.size())
                block.end += contig.blocks[index + 1].end - contig.blocks[index + 1].start;
            merged.push_back(block);
        }
        // HTSJDK's optimize() deliberately returns the last layout before a
        // merge that would collapse an occupied chromosome to one block.  The
        // resulting leading empty block is observable in v3 bytes and is
        // retained for ordinary VCF compatibility (GVCF layouts generally
        // start with one occupied bin and never enter this branch).
        if (merged.size() == 1) break;
        contig.blocks.swap(merged);
        if (contig.bin_width > max_bin_width / 2) contig.bin_width = max_bin_width + 1;
        else contig.bin_width *= 2;
    }
}

void write_i32(std::ofstream& output, std::int32_t value) {
    const std::uint32_t bits = static_cast<std::uint32_t>(value);
    for (int shift = 0; shift < 32; shift += 8)
        output.put(static_cast<char>((bits >> shift) & 0xffU));
}

void write_i64(std::ofstream& output, std::int64_t value) {
    const std::uint64_t bits = static_cast<std::uint64_t>(value);
    for (int shift = 0; shift < 64; shift += 8)
        output.put(static_cast<char>((bits >> shift) & 0xffU));
}

void write_tribble_string(std::ofstream& output, const std::string& value) {
    if (value.find('\0') != std::string::npos)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: Tribble index strings cannot contain NUL");
    output.write(value.data(), static_cast<std::streamsize>(value.size()));
    output.put('\0');
}

std::string file_uri(const std::string& path) {
    const std::string absolute = std::filesystem::absolute(path).generic_string();
    std::ostringstream uri;
    uri << "file://";
    constexpr char hex[] = "0123456789ABCDEF";
    for (const unsigned char character : absolute) {
        const bool unreserved = (character >= 'a' && character <= 'z') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= '0' && character <= '9') || character == '-' ||
            character == '_' || character == '.' || character == '~' || character == '/';
        if (unreserved) uri << static_cast<char>(character);
        else uri << '%' << hex[(character >> 4) & 0xf] << hex[character & 0xf];
    }
    return uri.str();
}

std::int64_t file_mtime_millis(const std::string& path) {
    struct stat status {};
    if (::stat(path.c_str(), &status) != 0)
        throw std::runtime_error("BAD_INPUT: cannot stat feature file: " + path);
    return static_cast<std::int64_t>(status.st_mtim.tv_sec) * 1000 +
        static_cast<std::int64_t>(status.st_mtim.tv_nsec / 1000000);
}

std::vector<std::string_view> split_tabs(std::string_view line) {
    std::vector<std::string_view> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const std::size_t end = line.find('\t', begin);
        if (end == std::string_view::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, end - begin));
        begin = end + 1;
    }
    return fields;
}

int parse_positive_coordinate(std::string_view value, const char* field) {
    if (value.empty()) throw std::runtime_error(std::string("BAD_INPUT: missing VCF ") + field);
    std::string copy(value);
    std::size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(copy, &consumed);
    } catch (...) {
        throw std::runtime_error(std::string("BAD_INPUT: invalid VCF ") + field + " " + copy);
    }
    if (consumed != copy.size() || parsed < 1 || parsed > std::numeric_limits<int>::max())
        throw std::runtime_error(std::string("BAD_INPUT: invalid VCF ") + field + " " + copy);
    return static_cast<int>(parsed);
}

int parse_nonnegative_coordinate(std::string_view value, const char* field) {
    if (value.empty()) throw std::runtime_error(std::string("BAD_INPUT: missing BED ") + field);
    std::string copy(value);
    std::size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(copy, &consumed);
    } catch (...) {
        throw std::runtime_error(std::string("BAD_INPUT: invalid BED ") + field + " " + copy);
    }
    if (consumed != copy.size() || parsed < 0 || parsed >= std::numeric_limits<int>::max())
        throw std::runtime_error(std::string("BAD_INPUT: invalid BED ") + field + " " + copy);
    return static_cast<int>(parsed);
}

struct BedFeature {
    std::string contig;
    int start = 0;  // 1-based inclusive, matching BEDCodec's Feature view.
    int end = 0;    // 1-based inclusive (the BED end is exclusive).
};

std::optional<BedFeature> parse_bed_feature(const std::string& line) {
    if (line.empty() || line.front() == '#' || line.rfind("track", 0) == 0 ||
        line.rfind("browser", 0) == 0) return std::nullopt;
    const auto fields = split_tabs(line);
    if (fields.size() < 3)
        throw std::runtime_error("BAD_INPUT: malformed BED record while building Tribble index");
    const std::string contig(fields[0]);
    if (contig.empty()) throw std::runtime_error("BAD_INPUT: BED record has empty contig");
    const int start0 = parse_nonnegative_coordinate(fields[1], "start");
    const int end0 = parse_nonnegative_coordinate(fields[2], "end");
    if (end0 <= start0)
        throw std::runtime_error("BAD_INPUT: BED end must be greater than start");
    return BedFeature{contig, start0 + 1, end0};
}

int parse_end_coordinate(std::string_view info, int start, int ref_length) {
    int end = start + std::max(1, ref_length) - 1;
    std::size_t begin = 0;
    while (begin <= info.size()) {
        const std::size_t next = info.find(';', begin);
        const std::string_view token = info.substr(begin,
            next == std::string_view::npos ? info.size() - begin : next - begin);
        if (token.rfind("END=", 0) == 0) {
            const int candidate = parse_positive_coordinate(token.substr(4), "END");
            end = candidate;
            break;
        }
        if (next == std::string_view::npos) break;
        begin = next + 1;
    }
    if (end < start)
        throw std::runtime_error("BAD_INPUT: VCF END precedes POS");
    return end;
}

FeatureSummary summarize_vcf(const Options& options, std::uint64_t input_size) {
    std::ifstream input(options.input, std::ios::binary);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot read VCF: " + options.input);

    FeatureSummary summary;
    std::string current_contig;
    std::string last_contig;
    int last_start = 0;
    std::size_t records_in_interval = 0;
    std::string line;
    while (true) {
        const std::streampos position = input.tellg();
        if (!std::getline(input, line)) break;
        if (position < 0) throw std::runtime_error("BAD_INPUT: invalid VCF byte offset");
        const auto offset = static_cast<std::uint64_t>(position);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 5)
            throw std::runtime_error("BAD_INPUT: malformed VCF record while building Tribble index");
        const std::string contig(fields[0]);
        if (contig.empty()) throw std::runtime_error("BAD_INPUT: VCF record has empty contig");
        const int start = parse_positive_coordinate(fields[1], "POS");
        const int end = parse_end_coordinate(fields.size() > 7 ? fields[7] : std::string_view{},
                                             start, static_cast<int>(fields[3].size()));

        if (summary.records != 0) {
            // DynamicIndexCreator's basesSeen is deliberately based on the
            // previous feature's start, but resets to the new coordinate
            // when a contig changes or records are out of order.
            if (contig == last_contig && start >= last_start)
                summary.bases_seen += static_cast<std::uint64_t>(start - last_start);
            else
                summary.bases_seen += static_cast<std::uint64_t>(start);
        } else {
            summary.bases_seen = static_cast<std::uint64_t>(start);
        }
        if (contig == last_contig && start < last_start)
            throw std::runtime_error("BAD_INPUT: VCF records are not sorted by coordinate");

        const int length = std::max(1, end - start + 1);
        summary.longest_feature = std::max(summary.longest_feature, length);
        push_feature_stat(summary, static_cast<double>(summary.longest_feature));
        ++summary.records;
        last_contig = contig;
        last_start = start;

        if (contig != current_contig) {
            if (!current_contig.empty() && !summary.interval_contigs.empty() &&
                !summary.interval_contigs.back().intervals.empty()) {
                summary.interval_contigs.back().intervals.back().file_end = offset;
            }
            // The input must be grouped by contig for both HTSJDK index
            // creators.  Reject interleaving rather than silently producing a
            // query-unsound index.
            for (const auto& existing : summary.interval_contigs)
                if (existing.name == contig)
                    throw std::runtime_error("BAD_INPUT: VCF contigs are not grouped in coordinate order: " + contig);
            summary.interval_contigs.push_back(IntervalContig{contig, {}});
            current_contig = contig;
            records_in_interval = 0;
        }
        auto& intervals = summary.interval_contigs.back().intervals;
        if (intervals.empty() || records_in_interval >= 75) {
            if (!intervals.empty()) intervals.back().file_end = offset;
            intervals.push_back(IntervalBlock{start, end, offset, 0});
            records_in_interval = 0;
        }
        intervals.back().end = std::max(intervals.back().end, end);
        ++records_in_interval;
    }
    if (!summary.interval_contigs.empty() && !summary.interval_contigs.back().intervals.empty())
        summary.interval_contigs.back().intervals.back().file_end = input_size;
    return summary;
}

FeatureSummary summarize_bed(const Options& options, std::uint64_t input_size) {
    std::ifstream input(options.input, std::ios::binary);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot read BED: " + options.input);

    FeatureSummary summary;
    std::string current_contig;
    std::string last_contig;
    int last_start = 0;
    std::size_t records_in_interval = 0;
    std::string line;
    while (true) {
        const std::streampos position = input.tellg();
        if (!std::getline(input, line)) break;
        if (position < 0) throw std::runtime_error("BAD_INPUT: invalid BED byte offset");
        const auto offset = static_cast<std::uint64_t>(position);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto feature = parse_bed_feature(line);
        if (!feature) continue;

        if (summary.records != 0) {
            if (feature->contig == last_contig && feature->start >= last_start)
                summary.bases_seen += static_cast<std::uint64_t>(feature->start - last_start);
            else
                summary.bases_seen += static_cast<std::uint64_t>(feature->start);
        } else {
            summary.bases_seen = static_cast<std::uint64_t>(feature->start);
        }
        if (feature->contig == last_contig && feature->start < last_start)
            throw std::runtime_error("BAD_INPUT: BED records are not sorted by coordinate");

        const int length = feature->end - feature->start + 1;
        summary.longest_feature = std::max(summary.longest_feature, length);
        push_feature_stat(summary, static_cast<double>(summary.longest_feature));
        ++summary.records;
        last_contig = feature->contig;
        last_start = feature->start;

        if (feature->contig != current_contig) {
            if (!current_contig.empty() && !summary.interval_contigs.empty() &&
                !summary.interval_contigs.back().intervals.empty())
                summary.interval_contigs.back().intervals.back().file_end = offset;
            for (const auto& existing : summary.interval_contigs)
                if (existing.name == feature->contig)
                    throw std::runtime_error("BAD_INPUT: BED contigs are not grouped in coordinate order: " +
                                             feature->contig);
            summary.interval_contigs.push_back(IntervalContig{feature->contig, {}});
            current_contig = feature->contig;
            records_in_interval = 0;
        }
        auto& intervals = summary.interval_contigs.back().intervals;
        if (intervals.empty() || records_in_interval >= 75) {
            if (!intervals.empty()) intervals.back().file_end = offset;
            intervals.push_back(IntervalBlock{feature->start, feature->end, offset, 0});
            records_in_interval = 0;
        }
        intervals.back().end = std::max(intervals.back().end, feature->end);
        ++records_in_interval;
    }
    if (!summary.interval_contigs.empty() && !summary.interval_contigs.back().intervals.empty())
        summary.interval_contigs.back().intervals.back().file_end = input_size;
    return summary;
}

LinearIndexStats write_linear_vcf_index(const Options& options, const std::string& index_path,
                                        std::uint64_t input_size, int initial_bin_width,
                                        const FeatureSummary* summary = nullptr) {
    std::ifstream input(options.input, std::ios::binary);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot read VCF: " + options.input);

    std::vector<LinearContig> contigs;
    std::map<std::string, std::size_t> seen_contigs;
    std::string current_contig;
    int last_start = 0;
    std::uint64_t offset = 0;
    std::string line;
    while (true) {
        const std::streampos position = input.tellg();
        if (!std::getline(input, line)) break;
        if (position < 0) throw std::runtime_error("BAD_INPUT: invalid VCF byte offset");
        offset = static_cast<std::uint64_t>(position);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tabs(line);
        if (fields.size() < 5)
            throw std::runtime_error("BAD_INPUT: malformed VCF record while building Tribble index");
        const std::string contig(fields[0]);
        if (contig.empty()) throw std::runtime_error("BAD_INPUT: VCF record has empty contig");
        const int start = parse_positive_coordinate(fields[1], "POS");
        const int end = parse_end_coordinate(fields.size() > 7 ? fields[7] : std::string_view{},
                                             start, static_cast<int>(fields[3].size()));

        if (contig != current_contig) {
            if (seen_contigs.find(contig) != seen_contigs.end())
                throw std::runtime_error("BAD_INPUT: VCF contigs are not grouped in coordinate order: " + contig);
            seen_contigs.emplace(contig, contigs.size());
            contigs.push_back(LinearContig{contig, {{offset, 0}}, initial_bin_width, 0, 0});
            current_contig = contig;
            last_start = 0;
        }
        if (start < last_start)
            throw std::runtime_error("BAD_INPUT: VCF records are not sorted by coordinate");
        last_start = start;

        auto& current = contigs.back();
        while (static_cast<std::int64_t>(start) >
               static_cast<std::int64_t>(current.blocks.size()) * initial_bin_width)
            current.blocks.push_back({offset, 0});
        current.longest_feature = std::max(current.longest_feature, end - start + 1);
        ++current.n_features;
    }

    for (std::size_t contig_index = 0; contig_index < contigs.size(); ++contig_index) {
        auto& contig = contigs[contig_index];
        // LinearIndexCreator closes the previous chromosome when the first
        // feature of the next chromosome is observed.  Only the final
        // chromosome runs to EOF; using EOF for every contig would make
        // cross-contig seeks read unrelated records.
        const std::uint64_t contig_end = contig_index + 1 < contigs.size()
            ? contigs[contig_index + 1].blocks.front().start : input_size;
        for (std::size_t index = 0; index < contig.blocks.size(); ++index) {
            contig.blocks[index].end = index + 1 < contig.blocks.size()
                ? contig.blocks[index + 1].start : contig_end;
            if (contig.blocks[index].end < contig.blocks[index].start)
                throw std::runtime_error("BAD_INPUT: invalid VCF block offsets");
        }
        optimize_linear_contig(contig);
    }

    std::ofstream output(index_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    write_i32(output, kTribbleMagic);
    write_i32(output, kLinearIndexType);
    write_i32(output, 3);  // AbstractIndex current version.
    write_tribble_string(output, file_uri(options.input));
    write_i64(output, static_cast<std::int64_t>(input_size));
    write_i64(output, file_mtime_millis(options.input));
    write_tribble_string(output, "");  // Java LinearIndex leaves MD5 unset.
    write_i32(output, 0);  // flags
    const bool dynamic_properties = summary != nullptr && initial_bin_width == 2000;
    write_i32(output, dynamic_properties ? 4 : 0);
    if (dynamic_properties) {
        write_tribble_string(output, "FEATURE_LENGTH_MEAN");
        write_tribble_string(output, java_double_string(summary->mean));
        write_tribble_string(output, "FEATURE_LENGTH_STD_DEV");
        write_tribble_string(output, java_double_string(summary_stddev(*summary)));
        write_tribble_string(output, "MEAN_FEATURE_VARIANCE");
        write_tribble_string(output, java_double_string(summary_variance(*summary)));
        write_tribble_string(output, "FEATURE_COUNT");
        write_tribble_string(output, std::to_string(summary->records));
    }
    write_i32(output, static_cast<std::int32_t>(contigs.size()));
    std::size_t block_count = 0;
    std::size_t record_count = 0;
    for (const auto& contig : contigs) {
        write_tribble_string(output, contig.name);
        write_i32(output, contig.bin_width);
        write_i32(output, static_cast<std::int32_t>(contig.blocks.size()));
        write_i32(output, contig.longest_feature);
        write_i32(output, 0);  // OLD_V3_INDEX=false
        write_i32(output, contig.n_features);
        // LinearIndex$ChrIndex.write stores each block start, followed by the
        // end of the final block.  The reader reconstructs each block as
        // [previous_start, next_start), with the last value as its end.
        std::uint64_t previous_end = 0;
        for (const auto& block : contig.blocks) {
            write_i64(output, static_cast<std::int64_t>(block.start));
            previous_end = block.end;
        }
        write_i64(output, static_cast<std::int64_t>(previous_end));
        block_count += contig.blocks.size();
        record_count += static_cast<std::size_t>(contig.n_features);
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble index: " + index_path);
    LinearIndexStats result{record_count, contigs.size(), block_count};
    if (summary != nullptr) {
        result.feature_length_mean = summary->mean;
        result.feature_length_stddev = summary_stddev(*summary);
        result.feature_length_variance = summary_variance(*summary);
    }
    return result;
}

using TribbleProperties = std::vector<std::pair<std::string, std::string>>;

void write_tribble_header(std::ofstream& output, int index_type, const Options& options,
                          std::uint64_t input_size, const TribbleProperties& properties);

LinearIndexStats write_linear_bed_index(const Options& options, const std::string& index_path,
                                        std::uint64_t input_size, int initial_bin_width,
                                        const FeatureSummary* summary = nullptr) {
    std::ifstream input(options.input, std::ios::binary);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot read BED: " + options.input);

    std::vector<LinearContig> contigs;
    std::map<std::string, std::size_t> seen_contigs;
    std::string current_contig;
    int last_start = 0;
    std::string line;
    while (true) {
        const std::streampos position = input.tellg();
        if (!std::getline(input, line)) break;
        if (position < 0) throw std::runtime_error("BAD_INPUT: invalid BED byte offset");
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto feature = parse_bed_feature(line);
        if (!feature) continue;

        if (feature->contig != current_contig) {
            if (seen_contigs.find(feature->contig) != seen_contigs.end())
                throw std::runtime_error("BAD_INPUT: BED contigs are not grouped in coordinate order: " +
                                         feature->contig);
            seen_contigs.emplace(feature->contig, contigs.size());
            contigs.push_back(LinearContig{feature->contig, {{static_cast<std::uint64_t>(position), 0}},
                                           initial_bin_width, 0, 0});
            current_contig = feature->contig;
            last_start = 0;
        }
        if (feature->start < last_start)
            throw std::runtime_error("BAD_INPUT: BED records are not sorted by coordinate");
        last_start = feature->start;

        auto& current = contigs.back();
        while (static_cast<std::int64_t>(feature->start) >
               static_cast<std::int64_t>(current.blocks.size()) * initial_bin_width)
            current.blocks.push_back({static_cast<std::uint64_t>(position), 0});
        current.longest_feature = std::max(current.longest_feature,
                                           feature->end - feature->start + 1);
        ++current.n_features;
    }

    for (std::size_t contig_index = 0; contig_index < contigs.size(); ++contig_index) {
        auto& contig = contigs[contig_index];
        const std::uint64_t contig_end = contig_index + 1 < contigs.size()
            ? contigs[contig_index + 1].blocks.front().start : input_size;
        for (std::size_t index = 0; index < contig.blocks.size(); ++index) {
            contig.blocks[index].end = index + 1 < contig.blocks.size()
                ? contig.blocks[index + 1].start : contig_end;
            if (contig.blocks[index].end < contig.blocks[index].start)
                throw std::runtime_error("BAD_INPUT: invalid BED block offsets");
        }
        optimize_linear_contig(contig);
    }

    std::ofstream output(index_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    write_tribble_header(output, kLinearIndexType, options, input_size,
                          summary == nullptr ? TribbleProperties{} : TribbleProperties{
                              {"FEATURE_LENGTH_MEAN", java_double_string(summary->mean)},
                              {"FEATURE_LENGTH_STD_DEV", java_double_string(summary_stddev(*summary))},
                              {"MEAN_FEATURE_VARIANCE", java_double_string(summary_variance(*summary))},
                              {"FEATURE_COUNT", std::to_string(summary->records)}});
    write_i32(output, static_cast<std::int32_t>(contigs.size()));
    std::size_t block_count = 0;
    std::size_t record_count = 0;
    for (const auto& contig : contigs) {
        write_tribble_string(output, contig.name);
        write_i32(output, contig.bin_width);
        write_i32(output, static_cast<std::int32_t>(contig.blocks.size()));
        write_i32(output, contig.longest_feature);
        write_i32(output, 0);  // OLD_V3_INDEX=false
        write_i32(output, contig.n_features);
        std::uint64_t previous_end = 0;
        for (const auto& block : contig.blocks) {
            write_i64(output, static_cast<std::int64_t>(block.start));
            previous_end = block.end;
        }
        write_i64(output, static_cast<std::int64_t>(previous_end));
        block_count += contig.blocks.size();
        record_count += static_cast<std::size_t>(contig.n_features);
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble index: " + index_path);
    LinearIndexStats result{record_count, contigs.size(), block_count};
    if (summary != nullptr) {
        result.feature_length_mean = summary->mean;
        result.feature_length_stddev = summary_stddev(*summary);
        result.feature_length_variance = summary_variance(*summary);
    }
    return result;
}

std::string java_double_string(double value) {
    // Java's Double.toString is shortest-round-trip.  std::to_chars(general)
    // uses the same shortest-round-trip principle; normalize its exponent
    // spelling and Java's explicit .0 for integral values.
    if (value == 0.0) return "0.0";
    std::array<char, 128> buffer{};
    const auto converted = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                         value, std::chars_format::general);
    std::string text;
    if (converted.ec == std::errc{}) text.assign(buffer.data(), converted.ptr);
    else {
        std::ostringstream stream;
        stream << std::setprecision(17) << value;
        text = stream.str();
    }
    const auto exponent = text.find_first_of("eE");
    if (exponent == std::string::npos) {
        if (text.find('.') == std::string::npos) text += ".0";
    } else {
        text[exponent] = 'E';
        std::size_t digits = exponent + 1;
        if (digits < text.size() && (text[digits] == '+' || text[digits] == '-')) ++digits;
        while (digits + 1 < text.size() && text[digits] == '0')
            text.erase(digits, 1);
        if (text.find('.') == std::string::npos)
            text.insert(exponent, ".0");
    }
    return text;
}

void write_tribble_header(std::ofstream& output, int index_type, const Options& options,
                          std::uint64_t input_size, const TribbleProperties& properties) {
    write_i32(output, kTribbleMagic);
    write_i32(output, index_type);
    write_i32(output, 3);  // AbstractIndex current version.
    write_tribble_string(output, file_uri(options.input));
    write_i64(output, static_cast<std::int64_t>(input_size));
    write_i64(output, file_mtime_millis(options.input));
    write_tribble_string(output, "");  // Java leaves MD5 unset.
    write_i32(output, 0);  // flags
    write_i32(output, static_cast<std::int32_t>(properties.size()));
    for (const auto& property : properties) {
        write_tribble_string(output, property.first);
        write_tribble_string(output, property.second);
    }
}

LinearIndexStats write_interval_vcf_index(const Options& options, const std::string& index_path,
                                          std::uint64_t input_size, const FeatureSummary& summary) {
    std::ofstream output(index_path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write Tribble index: " + index_path);
    const TribbleProperties properties = {
        {"FEATURE_LENGTH_MEAN", java_double_string(summary.mean)},
        {"FEATURE_LENGTH_STD_DEV", java_double_string(summary_stddev(summary))},
        {"MEAN_FEATURE_VARIANCE", java_double_string(summary_variance(summary))},
        {"FEATURE_COUNT", std::to_string(summary.records)},
    };
    write_tribble_header(output, 2, options, input_size, properties);
    write_i32(output, static_cast<std::int32_t>(summary.interval_contigs.size()));
    std::size_t interval_count = 0;

    // HTSJDK's IntervalTreeIndex does not serialize intervals in coordinate
    // order.  IntervalTree.getIntervals() is a pre-order walk of the
    // red-black tree built by IntervalTree.insert(): equal starts are sent to
    // the left child, then the usual CLRS insert-fixup rotations are applied.
    // Reproduce that tree shape here so native dense indexes can be compared
    // byte-for-byte with DynamicIndexCreator, rather than only being query
    // compatible.  The min/max augmentation is used by HTSJDK for queries but
    // is not serialized; the topology and insertion comparator determine the
    // emitted order.
    struct TreeNode {
        IntervalBlock interval;
        int parent = -1;
        int left = -1;
        int right = -1;
        bool red = true;
    };

    auto emit_node = [&](const TreeNode& node) {
        if (node.interval.end < node.interval.start ||
            node.interval.file_end < node.interval.file_start ||
            node.interval.file_end - node.interval.file_start >
                static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: invalid interval block while writing Tribble index");
        }
        write_i32(output, node.interval.start);
        write_i32(output, node.interval.end);
        write_i64(output, static_cast<std::int64_t>(node.interval.file_start));
        write_i32(output, static_cast<std::int32_t>(
            node.interval.file_end - node.interval.file_start));
        ++interval_count;
    };

    for (const auto& contig : summary.interval_contigs) {
        write_tribble_string(output, contig.name);
        write_i32(output, static_cast<std::int32_t>(contig.intervals.size()));
        std::vector<TreeNode> nodes;
        nodes.reserve(contig.intervals.size());
        int root = -1;

        auto is_red = [&](const int index) {
            return index >= 0 && nodes[static_cast<std::size_t>(index)].red;
        };
        auto left_rotate = [&](const int x) {
            const int y = nodes[static_cast<std::size_t>(x)].right;
            if (y < 0) throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: malformed interval tree rotation");
            nodes[static_cast<std::size_t>(x)].right =
                nodes[static_cast<std::size_t>(y)].left;
            if (nodes[static_cast<std::size_t>(y)].left >= 0)
                nodes[static_cast<std::size_t>(nodes[static_cast<std::size_t>(y)].left)].parent = x;
            nodes[static_cast<std::size_t>(y)].parent =
                nodes[static_cast<std::size_t>(x)].parent;
            if (nodes[static_cast<std::size_t>(x)].parent < 0) {
                root = y;
            } else if (x == nodes[static_cast<std::size_t>(
                           nodes[static_cast<std::size_t>(x)].parent)].left) {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].left = y;
            } else {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].right = y;
            }
            nodes[static_cast<std::size_t>(y)].left = x;
            nodes[static_cast<std::size_t>(x)].parent = y;
        };
        auto right_rotate = [&](const int x) {
            const int y = nodes[static_cast<std::size_t>(x)].left;
            if (y < 0) throw std::runtime_error(
                "OUTPUT_CONTRACT_FAILURE: malformed interval tree rotation");
            nodes[static_cast<std::size_t>(x)].left =
                nodes[static_cast<std::size_t>(y)].right;
            if (nodes[static_cast<std::size_t>(y)].right >= 0)
                nodes[static_cast<std::size_t>(nodes[static_cast<std::size_t>(y)].right)].parent = x;
            nodes[static_cast<std::size_t>(y)].parent =
                nodes[static_cast<std::size_t>(x)].parent;
            if (nodes[static_cast<std::size_t>(x)].parent < 0) {
                root = y;
            } else if (x == nodes[static_cast<std::size_t>(
                           nodes[static_cast<std::size_t>(x)].parent)].right) {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].right = y;
            } else {
                nodes[static_cast<std::size_t>(
                    nodes[static_cast<std::size_t>(x)].parent)].left = y;
            }
            nodes[static_cast<std::size_t>(y)].right = x;
            nodes[static_cast<std::size_t>(x)].parent = y;
        };
        auto insert_fixup = [&](int z) {
            while (z != root && is_red(nodes[static_cast<std::size_t>(z)].parent)) {
                const int parent = nodes[static_cast<std::size_t>(z)].parent;
                const int grand = nodes[static_cast<std::size_t>(parent)].parent;
                if (parent == nodes[static_cast<std::size_t>(grand)].left) {
                    const int uncle = nodes[static_cast<std::size_t>(grand)].right;
                    if (is_red(uncle)) {
                        nodes[static_cast<std::size_t>(parent)].red = false;
                        nodes[static_cast<std::size_t>(uncle)].red = false;
                        nodes[static_cast<std::size_t>(grand)].red = true;
                        z = grand;
                    } else {
                        if (z == nodes[static_cast<std::size_t>(parent)].right) {
                            z = parent;
                            left_rotate(z);
                        }
                        const int fixed_parent = nodes[static_cast<std::size_t>(z)].parent;
                        const int fixed_grand = nodes[static_cast<std::size_t>(fixed_parent)].parent;
                        nodes[static_cast<std::size_t>(fixed_parent)].red = false;
                        nodes[static_cast<std::size_t>(fixed_grand)].red = true;
                        right_rotate(fixed_grand);
                    }
                } else {
                    const int uncle = nodes[static_cast<std::size_t>(grand)].left;
                    if (is_red(uncle)) {
                        nodes[static_cast<std::size_t>(parent)].red = false;
                        nodes[static_cast<std::size_t>(uncle)].red = false;
                        nodes[static_cast<std::size_t>(grand)].red = true;
                        z = grand;
                    } else {
                        if (z == nodes[static_cast<std::size_t>(parent)].left) {
                            z = parent;
                            right_rotate(z);
                        }
                        const int fixed_parent = nodes[static_cast<std::size_t>(z)].parent;
                        const int fixed_grand = nodes[static_cast<std::size_t>(fixed_parent)].parent;
                        nodes[static_cast<std::size_t>(fixed_parent)].red = false;
                        nodes[static_cast<std::size_t>(fixed_grand)].red = true;
                        left_rotate(fixed_grand);
                    }
                }
            }
            if (root >= 0) nodes[static_cast<std::size_t>(root)].red = false;
        };

        for (const auto& interval : contig.intervals) {
            int parent = -1;
            int cursor = root;
            while (cursor >= 0) {
                parent = cursor;
                // HTSJDK's treeInsert uses <= for the left branch, so equal
                // starts are deliberately not made stable on the right.
                cursor = interval.start <= nodes[static_cast<std::size_t>(cursor)].interval.start
                    ? nodes[static_cast<std::size_t>(cursor)].left
                    : nodes[static_cast<std::size_t>(cursor)].right;
            }
            const int inserted = static_cast<int>(nodes.size());
            nodes.push_back(TreeNode{interval, parent, -1, -1, true});
            if (parent < 0) root = inserted;
            else if (interval.start <= nodes[static_cast<std::size_t>(parent)].interval.start)
                nodes[static_cast<std::size_t>(parent)].left = inserted;
            else nodes[static_cast<std::size_t>(parent)].right = inserted;
            insert_fixup(inserted);
        }

        std::function<void(int)> preorder = [&](const int index) {
            if (index < 0) return;
            emit_node(nodes[static_cast<std::size_t>(index)]);
            preorder(nodes[static_cast<std::size_t>(index)].left);
            preorder(nodes[static_cast<std::size_t>(index)].right);
        };
        preorder(root);
    }
    output.flush();
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed writing Tribble interval index: " + index_path);
    return {summary.records, summary.interval_contigs.size(), interval_count,
            summary.mean, summary_stddev(summary), summary_variance(summary)};
}

int build_index(const Options& options, const DetectedFormat& detected, const std::string& index_path,
                std::string& index_type, LinearIndexStats* linear_stats = nullptr,
                std::uint64_t input_size = 0) {
    if (detected.format == bcf) {
        if (!index_path.empty() && (index_path.size() < 4 ||
            index_path.substr(index_path.size() - 4) != ".csi")) {
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: BCF index output must end in .csi");
        }
        index_type = "CSI";
        return bcf_index_build3(options.input.c_str(), index_path.c_str(), 14, options.threads);
    }

    if (detected.format != vcf && detected.format != bed)
        throw std::runtime_error("BACKEND_UNAVAILABLE: IndexFeatureFile supports BGZF VCF/GVCF/BED and BCF; use --fallback for this codec");
    if (detected.compression == no_compression &&
        (detected.format == vcf || detected.format == bed)) {
        const bool is_gvcf = options.input.size() >= 6 && options.input.ends_with(".g.vcf");
        if (is_gvcf && detected.format != vcf)
            throw std::runtime_error("BAD_INPUT: a BED path cannot use the .g.vcf suffix");
        const FeatureSummary summary = detected.format == vcf
            ? summarize_vcf(options, input_size) : summarize_bed(options, input_size);
        const auto write_linear = [&](int initial_bin_width) {
            return detected.format == vcf
                ? write_linear_vcf_index(options, index_path, input_size,
                                         initial_bin_width, &summary)
                : write_linear_bed_index(options, index_path, input_size,
                                          initial_bin_width, &summary);
        };
        if (is_gvcf) {
            index_type = "TRIBBLE_LINEAR_GVCF";
            if (linear_stats) *linear_stats = write_linear_vcf_index(
                options, index_path, input_size, kGvcfBinWidth, &summary);
        } else {
            // DynamicIndexCreator(FOR_SEEK_TIME) evaluates a 2000 bp linear
            // candidate against 75-features-per-interval.  Empty files use
            // the interval-tree representation (and still carry properties).
            const double density = summary.bases_seen == 0 ? 0.0 :
                static_cast<double>(summary.records) / static_cast<double>(summary.bases_seen);
            const double linear_score = 2000.0 * density *
                std::ceil(static_cast<double>(summary.longest_feature) / 2000.0);
            if (summary.records == 0 || !(linear_score < 75.0)) {
                index_type = "TRIBBLE_INTERVAL_TREE";
                if (linear_stats) *linear_stats = write_interval_vcf_index(
                    options, index_path, input_size, summary);
            } else {
                index_type = "TRIBBLE_LINEAR";
                if (linear_stats) *linear_stats = write_linear(2000);
            }
        }
        return 0;
    }
    if (detected.compression != bgzf)
        throw std::runtime_error("BACKEND_UNAVAILABLE: uncompressed Tribble dynamic/linear indices for this feature format are not implemented natively; use --fallback");
    if (index_path.size() < 4 || index_path.substr(index_path.size() - 4) != ".tbi")
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: VCF/BED index output must end in .tbi");
    index_type = "TABIX";
    const tbx_conf_t* configuration = detected.format == bed ? &tbx_conf_bed : &tbx_conf_vcf;
    return tbx_index_build3(options.input.c_str(), index_path.c_str(), 0, options.threads, configuration);
}

int main_impl(int argc, char** argv) {
    const Options options = parse_options(argc, argv);
    const DetectedFormat detected = detect_format(options.input);
    const std::string index_path = default_index_path(options, detected);
    if (index_path.empty()) throw std::invalid_argument("-O/--output cannot be empty");

    std::error_code error;
    const auto input_size = std::filesystem::file_size(options.input, error);
    if (error || input_size > std::numeric_limits<std::size_t>::max())
        throw std::runtime_error("BAD_INPUT: cannot determine feature file size: " + options.input);

    Kokkos::InitializationSettings settings;
    settings.set_num_threads(options.threads);
    Kokkos::initialize(settings);
    bool initialized = true;
    try {
        {
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host("index-feature-file-v1");
        host.records = 1;
        host.bytes = static_cast<std::size_t>(input_size);
        fastgatk::core::KernelPlan<ExecSpace> plan("index-feature-file");
        plan.begin_prepare(host);
        Kokkos::View<std::uint64_t> source_size("feature_file_size");
        Kokkos::View<std::uint64_t> observed_size("observed_feature_file_size");
        auto host_source = Kokkos::create_mirror_view(source_size);
        host_source() = input_size;
        Kokkos::deep_copy(source_size, host_source);
        fastgatk::core::DeviceBatch<ExecSpace> device(1);
        device.bind("source_size", source_size);
        device.bind("observed_size", observed_size);
        ExecSpace().fence();
        plan.end_prepare(device);
        plan.begin_execute();
        Kokkos::parallel_for("index_feature_file_metadata", Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int) { observed_size() = source_size(); });
        ExecSpace().fence();
        plan.end_execute();
        auto host_observed = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), observed_size);
        if (host_observed() != input_size)
            throw std::runtime_error("INTERNAL_ERROR: Kokkos metadata projection changed input size");

        std::string index_type;
        LinearIndexStats linear_stats;
        const int result = build_index(options, detected, index_path, index_type, &linear_stats, input_size);
        if (result != 0 || !std::filesystem::is_regular_file(index_path) ||
            std::filesystem::file_size(index_path, error) == 0 || error) {
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: HTSlib failed to create index " + index_path +
                                     " (return code " + std::to_string(result) + ")");
        }

        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            const auto& telemetry = plan.telemetry();
            manifest << "{\"schema_version\":1,\"tool\":\"IndexFeatureFile\",\"implementation\":\"fastgatk-index-feature-file\",\"status\":\"contract-compatible\",\"input\":\""
                     << json_escape(options.input) << "\",\"index\":\"" << json_escape(index_path)
                     << "\",\"format\":\"" << format_name(detected.format) << "\",\"compression\":\""
                     << compression_name(detected.compression) << "\",\"index_type\":\"" << index_type
                     << "\",\"linear_records\":" << linear_stats.records
                     << ",\"linear_contigs\":" << linear_stats.contigs
                     << ",\"linear_blocks\":" << linear_stats.blocks
                     << ",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\""
                     << ExecSpace::name() << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_prepare_seconds\":"
                     << telemetry.prepare_seconds << ",\"kernel_execute_seconds\":" << telemetry.execute_seconds
                     << ",\"input_bytes\":" << input_size << "}}\n";
        }
        std::cerr << "{\"tool\":\"IndexFeatureFile\",\"status\":\"contract-compatible\",\"format\":\""
                  << format_name(detected.format) << "\",\"compression\":\"" << compression_name(detected.compression)
                  << "\",\"index\":\"" << json_escape(index_path) << "\",\"index_type\":\"" << index_type
                  << "\",\"linear_records\":" << linear_stats.records
                  << ",\"linear_contigs\":" << linear_stats.contigs
                  << ",\"linear_blocks\":" << linear_stats.blocks
                  << ",\"kernel_execution_space\":\"" << ExecSpace::name() << "\",\"kernel_execute_seconds\":"
                  << plan.telemetry().execute_seconds << "}\n";
        }
        Kokkos::finalize();
        initialized = false;
        return 0;
    } catch (...) {
        if (initialized) Kokkos::finalize();
        throw;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return main_impl(argc, argv);
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
