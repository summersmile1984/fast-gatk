#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/output.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cctype>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::string reference;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> comments;
    std::string output;
    std::string manifest;
    bool allow_overlaps = false;
    bool reorder_input_by_first_variant = false;
    bool disable_sequence_dictionary_validation = false;
    bool quiet = false;
    bool create_index = true;
    int compression_level = 2;
};

struct IntervalInputStats {
    std::size_t files = 0;
    std::size_t records = 0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* long_name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, long_name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == long_name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto character : text) {
        if (character == '"' || character == '\\') escaped << '\\';
        if (character == '\n') escaped << "\\n";
        else if (character == '\r') escaped << "\\r";
        else if (character == '\t') escaped << "\\t";
        else escaped << character;
    }
    return escaped.str();
}

bool suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    return path.size() >= value.size() && path.compare(path.size() - value.size(), value.size(), value) == 0;
}

bool file_complete(const std::string& path) {
    std::error_code error;
    const bool regular = std::filesystem::is_regular_file(path, error);
    if (error || !regular) return false;
    const auto size = std::filesystem::file_size(path, error);
    return !error && size > 0;
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-gather-vcfs (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                 input VCF/BCF (repeatable)\n"
                         "  -L, --intervals REGION           interval/file selector (repeatable)\n"
                         "      --interval-set-rule RULE     UNION (default) or INTERSECTION\n"
                         "  -CO, --COMMENT TEXT              add GatherVcfs comment metadata\n"
                         "  -O, --output FILE                gathered VCF/BCF\n"
                         "      --reorder-input-by-first-variant, -RI\n"
                         "                                  reorder shards by first record\n"
                         "      --CREATE_INDEX BOOL          write a VCF index (default true)\n"
                         "      --COMPRESSION_LEVEL INT      BGZF compression level (default 2)\n"
                         "      --allow-overlaps             allow non-disjoint shard coordinates\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-I" || is_option(argument, "--input") || is_option(argument, "--INPUT")) {
            const char* name = is_option(argument, "--INPUT") ? "--INPUT" : "--input";
            options.inputs.push_back(require_value(index, argc, argv, argument, name, "-I"));
        }
        else if (argument == "-R" || is_option(argument, "--reference") ||
                 is_option(argument, "--REFERENCE_SEQUENCE")) {
            const char* name = is_option(argument, "--REFERENCE_SEQUENCE")
                ? "--REFERENCE_SEQUENCE" : "--reference";
            options.reference = require_value(index, argc, argv, argument, name, "-R");
        }
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                is_option(argument, "--region") ? "--region" :
                is_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (argument == "-isr" || is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule", "-isr");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (argument == "-CO" || is_option(argument, "--COMMENT"))
            options.comments.push_back(require_value(index, argc, argv, argument, "--COMMENT", "-CO"));
        else if (argument == "-O" || is_option(argument, "--output") || is_option(argument, "--OUTPUT")) {
            const char* name = is_option(argument, "--OUTPUT") ? "--OUTPUT" : "--output";
            options.output = require_value(index, argc, argv, argument, name, "-O");
        }
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--allow-overlaps")
            options.allow_overlaps = true;
        else if (argument == "-RI" || argument == "--reorder-input-by-first-variant" ||
                 argument.rfind("--reorder-input-by-first-variant=", 0) == 0 ||
                 argument == "--REORDER_INPUT_BY_FIRST_VARIANT" ||
                 argument.rfind("--REORDER_INPUT_BY_FIRST_VARIANT=", 0) == 0) {
            const char* name = argument.rfind("--REORDER_INPUT_BY_FIRST_VARIANT", 0) == 0
                ? "--REORDER_INPUT_BY_FIRST_VARIANT" : "--reorder-input-by-first-variant";
            options.reorder_input_by_first_variant = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        }
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0 ||
                 argument == "--CREATE_INDEX" || argument.rfind("--CREATE_INDEX=", 0) == 0) {
            const char* name = argument.rfind("--CREATE_INDEX", 0) == 0
                ? "--CREATE_INDEX" : "--create-output-variant-index";
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (is_option(argument, "--COMPRESSION_LEVEL")) {
            const auto value = require_value(index, argc, argv, argument, "--COMPRESSION_LEVEL");
            try {
                options.compression_level = std::stoi(value);
            } catch (...) {
                throw std::invalid_argument("invalid integer for --COMPRESSION_LEVEL: " + value);
            }
            if (options.compression_level < 0 || options.compression_level > 9)
                throw std::invalid_argument("--COMPRESSION_LEVEL must be between 0 and 9");
        } else if (argument == "--quiet" || argument == "--QUIET") {
            const char* name = argument == "--QUIET" ? "--QUIET" : "--quiet";
            options.quiet = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (argument == "--disable-sequence-dictionary-validation") {
            options.disable_sequence_dictionary_validation = true;
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--VERBOSITY")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--VERBOSITY", 0) == 0 ? "--VERBOSITY" : "--verbosity"));
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.inputs.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    return options;
}

#if FASTGATK_HAS_HTSLIB

void close_inputs(std::vector<htsFile*>& files, std::vector<bcf_hdr_t*>& headers) {
    for (auto* header : headers) bcf_hdr_destroy(header);
    for (auto* file : files) bcf_close(file);
    headers.clear();
    files.clear();
}

struct Region {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
};

struct FirstVariantCoordinate {
    std::string path;
    std::string contig;
    int position = std::numeric_limits<int>::max();
    int contig_order = std::numeric_limits<int>::max();
    bool has_record = false;
    std::size_t input_order = 0;
};

struct ContigDictionaryEntry {
    std::string name;
    std::int64_t length = -1;
};

std::vector<ContigDictionaryEntry> contig_dictionary(const bcf_hdr_t* header) {
    int count = 0;
    const char** names = bcf_hdr_seqnames(header, &count);
    if (names == nullptr && count != 0)
        throw std::runtime_error("BAD_INPUT: cannot enumerate VCF sequence dictionary");
    std::vector<ContigDictionaryEntry> result;
    result.reserve(static_cast<std::size_t>(std::max(0, count)));
    for (int index = 0; index < count; ++index) {
        if (names[index] == nullptr) continue;
        ContigDictionaryEntry entry;
        entry.name = names[index];
        const auto* hrec = bcf_hdr_get_hrec(header, BCF_HL_CTG, "ID", names[index], nullptr);
        if (hrec != nullptr) {
            for (int field = 0; field < hrec->nkeys; ++field) {
                if (hrec->keys[field] != nullptr && hrec->vals[field] != nullptr &&
                    std::string(hrec->keys[field]) == "length") {
                    try {
                        entry.length = std::stoll(hrec->vals[field]);
                    } catch (...) {
                        throw std::runtime_error("BAD_INPUT: invalid contig length for " + entry.name);
                    }
                }
            }
        }
        result.push_back(std::move(entry));
    }
    free(const_cast<char**>(names));
    return result;
}

void validate_sequence_dictionary(const bcf_hdr_t* expected, const bcf_hdr_t* observed,
                                  const std::string& path) {
    const auto lhs = contig_dictionary(expected);
    const auto rhs = contig_dictionary(observed);
    if (lhs.size() != rhs.size())
        throw std::runtime_error("BAD_INPUT: GatherVcfs sequence dictionary count differs in " + path);
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        if (lhs[index].name != rhs[index].name ||
            (lhs[index].length > 0 && rhs[index].length > 0 && lhs[index].length != rhs[index].length)) {
            throw std::runtime_error("BAD_INPUT: GatherVcfs sequence dictionary mismatch in " + path +
                                     " at contig " + std::to_string(index));
        }
    }
}

std::vector<std::string> reorder_inputs_by_first_variant(const std::vector<std::string>& paths) {
    if (paths.size() < 2) return paths;

    // Use the first shard's sequence dictionary as the ordering authority,
    // matching htsjdk's SAMSequenceDictionary comparator.  Header validation
    // still happens in the normal gather pass; this scan only decides shard
    // order and never emits records.
    htsFile* reference_file = bcf_open(paths.front().c_str(), "r");
    if (!reference_file) throw std::runtime_error("BAD_INPUT: cannot open VCF/BCF: " + paths.front());
    bcf_hdr_t* reference_header = bcf_hdr_read(reference_file);
    if (!reference_header) {
        bcf_close(reference_file);
        throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + paths.front());
    }
    std::map<std::string, int> contig_order;
    for (int rid = 0;; ++rid) {
        const char* name = bcf_hdr_id2name(reference_header, rid);
        if (!name) break;
        contig_order.emplace(name, rid);
    }
    bcf_hdr_destroy(reference_header);
    bcf_close(reference_file);

    std::vector<FirstVariantCoordinate> coordinates;
    coordinates.reserve(paths.size());
    for (std::size_t input_order = 0; input_order < paths.size(); ++input_order) {
        const auto& path = paths[input_order];
        htsFile* file = bcf_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open VCF/BCF: " + path);
        bcf_hdr_t* header = bcf_hdr_read(file);
        if (!header) {
            bcf_close(file);
            throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + path);
        }
        bcf1_t* record = bcf_init();
        if (!record) {
            bcf_hdr_destroy(header);
            bcf_close(file);
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate VCF record for reorder scan");
        }
        FirstVariantCoordinate coordinate;
        coordinate.path = path;
        coordinate.input_order = input_order;
        if (bcf_read(file, header, record) == 0) {
            bcf_unpack(record, BCF_UN_STR);
            if (record->rid < 0 || record->pos < 0) {
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(file);
                throw std::runtime_error("BAD_INPUT: GatherVcfs reorder scan encountered an invalid coordinate in " + path);
            }
            const char* contig = bcf_hdr_id2name(header, record->rid);
            if (!contig) {
                bcf_destroy(record);
                bcf_hdr_destroy(header);
                bcf_close(file);
                throw std::runtime_error("BAD_INPUT: GatherVcfs reorder scan encountered an unknown contig in " + path);
            }
            coordinate.contig = contig;
            coordinate.position = static_cast<int>(record->pos);
            coordinate.has_record = true;
            const auto found = contig_order.find(coordinate.contig);
            coordinate.contig_order = found == contig_order.end()
                ? std::numeric_limits<int>::max() : found->second;
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(file);
        coordinates.push_back(std::move(coordinate));
    }
    std::stable_sort(coordinates.begin(), coordinates.end(), [](const auto& left, const auto& right) {
        if (left.has_record != right.has_record) return left.has_record > right.has_record;
        if (left.contig_order != right.contig_order) return left.contig_order < right.contig_order;
        if (left.contig != right.contig) return left.contig < right.contig;
        if (left.position != right.position) return left.position < right.position;
        return left.input_order < right.input_order;
    });
    std::vector<std::string> ordered;
    ordered.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) ordered.push_back(coordinate.path);
    return ordered;
}

Region parse_region_token(const std::string& text, const bcf_hdr_t* header,
                          bool bed_coordinates = false) {
    Region region;
    const auto colon = text.find(':');
    const auto contig = colon == std::string::npos ? text : text.substr(0, colon);
    region.rid = bcf_hdr_name2id(header, contig.c_str());
    if (region.rid < 0)
        throw std::invalid_argument("BAD_INPUT: interval contig not in VCF header: " + contig);
    if (colon == std::string::npos) return region;
    const auto dash = text.find('-', colon + 1);
    const auto begin = std::stoll(text.substr(
        colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1));
    if (begin < 1 || begin > std::numeric_limits<int>::max())
        throw std::invalid_argument("invalid interval start");
    region.begin = static_cast<int>(bed_coordinates ? begin : begin - 1);
    if (dash != std::string::npos) {
        const auto end = std::stoll(text.substr(dash + 1));
        if (end < begin || end > std::numeric_limits<int>::max() ||
            (bed_coordinates && end == begin))
            throw std::invalid_argument("invalid interval end");
        region.end = static_cast<int>(bed_coordinates ? end : end);
    }
    return region;
}

bool has_suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    for (std::size_t index = 0; index < ending.size(); ++index) {
        const auto lhs = static_cast<char>(std::tolower(static_cast<unsigned char>(value[value.size() - ending.size() + index])));
        const auto rhs = static_cast<char>(std::tolower(static_cast<unsigned char>(ending[index])));
        if (lhs != rhs) return false;
    }
    return true;
}

void append_interval_input(const std::string& text, const bcf_hdr_t* header,
                           std::vector<Region>& regions, IntervalInputStats& stats) {
    std::error_code error;
    const bool regular = std::filesystem::is_regular_file(text, error) && !error;
    // GATK accepts interval files without requiring a particular suffix.  Do
    // not classify an obvious VCF/BCF as an interval file, but otherwise let
    // regular files use the interval-list parser (including temporary files
    // created by workflow engines with opaque names).
    const bool variant_file = has_suffix_ci(text, ".vcf") || has_suffix_ci(text, ".vcf.gz") ||
                              has_suffix_ci(text, ".bcf") || has_suffix_ci(text, ".bcf.gz");
    const bool interval_file = regular && !variant_file;
    if (!regular || !interval_file) {
        regions.push_back(parse_region_token(text, header));
        return;
    }
    ++stats.files;
    const bool bed = has_suffix_ci(text, ".bed") || has_suffix_ci(text, ".bed.gz");
    const auto lines = fastgatk::io::read_interval_lines(text);
    std::size_t records = 0;
    for (std::size_t line_number = 0; line_number < lines.size(); ++line_number) {
        const auto& line = lines[line_number];
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '@') continue;
        std::istringstream fields(line.substr(first));
        std::vector<std::string> tokens;
        for (std::string token; fields >> token;) tokens.push_back(std::move(token));
        if (!tokens.empty() && (tokens.front() == "track" || tokens.front() == "browser")) continue;
        // GATK interval files commonly use either Picard interval-list rows
        // (contig, start, end, strand, name) or one interval literal per
        // line (chr:start-end / chr).  Accept both forms; the former remains
        // one-based inclusive while literals use the same parser as -L.
        if (tokens.size() == 1) {
            regions.push_back(parse_region_token(tokens.front(), header));
            ++records;
            ++stats.records;
            continue;
        }
        if (tokens.size() < 3)
            throw std::invalid_argument("BAD_INPUT: malformed interval list line " +
                                        std::to_string(line_number + 1) + " in " + text);
        const auto rid = bcf_hdr_name2id(header, tokens[0].c_str());
        if (rid < 0)
            throw std::invalid_argument("BAD_INPUT: interval contig not in VCF header: " + tokens[0]);
        long long begin = 0;
        long long end = 0;
        try {
            begin = std::stoll(tokens[1]);
            end = std::stoll(tokens[2]);
        } catch (...) {
            throw std::invalid_argument("BAD_INPUT: malformed interval coordinates on line " +
                                        std::to_string(line_number + 1) + " in " + text);
        }
        if (bed) {
            if (begin < 0 || end <= begin || end > std::numeric_limits<int>::max())
            throw std::invalid_argument("BAD_INPUT: invalid BED interval on line " + std::to_string(line_number + 1));
        } else if (begin < 1 || end < begin || end > std::numeric_limits<int>::max()) {
            throw std::invalid_argument("BAD_INPUT: invalid interval-list coordinates on line " + std::to_string(line_number + 1));
        }
        Region region;
        region.rid = rid;
        region.begin = static_cast<int>(bed ? begin : begin - 1);
        region.end = static_cast<int>(end);
        regions.push_back(region);
        ++records;
        ++stats.records;
    }
    if (records == 0) throw std::invalid_argument("BAD_INPUT: interval list is empty: " + text);
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_gather_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.end <= region.begin) continue;
        if (!merged.empty() && merged.back().rid == region.rid && region.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, region.end);
        else merged.push_back(region);
    }
    regions.swap(merged);
}

void append_interval_input_with_rule(
    const std::string& text,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    IntervalInputStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_interval_input(text, header, incoming, stats);
    normalize_gather_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < regions[left].rid) { ++right; continue; }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end) {
            Region overlap = regions[left];
            overlap.begin = begin;
            overlap.end = end;
            intersection.push_back(overlap);
        }
        if (regions[left].end < incoming[right].end) ++left;
        else ++right;
    }
    regions.swap(intersection);
}

void validate_sample_headers(const bcf_hdr_t* expected, const bcf_hdr_t* actual,
                             const std::string& path) {
    if (expected->n[BCF_DT_SAMPLE] != actual->n[BCF_DT_SAMPLE])
        throw std::runtime_error("BAD_INPUT: GatherVcfs sample count differs in " + path);
    for (int index = 0; index < expected->n[BCF_DT_SAMPLE]; ++index) {
        const char* lhs = bcf_hdr_int2id(expected, BCF_DT_SAMPLE, index);
        const char* rhs = bcf_hdr_int2id(actual, BCF_DT_SAMPLE, index);
        if (lhs == nullptr || rhs == nullptr || std::string(lhs) != rhs)
            throw std::runtime_error("BAD_INPUT: GatherVcfs sample order differs in " + path);
    }
}

int record_end_exclusive(const bcf_hdr_t* header, const bcf1_t* record) {
    // VCF END is one-based inclusive, while HTSlib positions are zero-based.
    // Keep the interval contract in zero-based half-open coordinates so a
    // gVCF reference block cannot be accidentally gathered across a shard.
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, const_cast<bcf1_t*>(record),
                                           "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, static_cast<int>(values[0]));
    free(values);
    return end;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    const auto wall_begin = std::chrono::steady_clock::now();
    std::vector<htsFile*> inputs;
    std::vector<bcf_hdr_t*> headers;
    bcf_hdr_t* output_header = nullptr;
    htsFile* output = nullptr;
    bcf1_t* record = nullptr;
    std::uint64_t input_records = 0;
    std::uint64_t interval_skipped = 0;
    std::int32_t previous_rid = -1;
    std::int32_t previous_pos = -1;
    std::int32_t previous_end_exclusive = -1;
    IntervalInputStats interval_input_stats;
    try {
        // HTSlib 1.22 emits BCF 2.2 for native writes, whereas the pinned
        // GATK/htsjdk reader used by the replacement contract only accepts
        // the uncompressed BCF 2.1 form produced by GATK itself.  Do not
        // publish a file that downstream GATK cannot read; the dispatcher can
        // route this explicit boundary to the Java fallback instead.
        if (has_suffix_ci(options.output, ".bcf") || has_suffix_ci(options.output, ".bcf.gz"))
            throw std::runtime_error(
                "BACKEND_UNAVAILABLE: native GatherVcfs cannot emit GATK-readable BCF 2.1; "
                "use --fallback for BCF output");
        const auto input_paths = options.reorder_input_by_first_variant
            ? reorder_inputs_by_first_variant(options.inputs) : options.inputs;
        // Read and merge all headers before writing the output header. This
        // preserves FILTER/INFO/FORMAT definitions from later scatter shards.
        for (const auto& path : input_paths) {
            auto* input = bcf_open(path.c_str(), "r");
            if (!input) throw std::runtime_error("BAD_INPUT: cannot open VCF/BCF: " + path);
            auto* header = bcf_hdr_read(input);
            if (!header) {
                bcf_close(input);
                throw std::runtime_error("BAD_INPUT: cannot read VCF/BCF header: " + path);
            }
            if (!output_header) {
                output_header = bcf_hdr_dup(header);
                if (!output_header) {
                    bcf_hdr_destroy(header);
                    bcf_close(input);
                    throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate VCF header");
                }
            } else {
                validate_sample_headers(output_header, header, path);
                if (!options.disable_sequence_dictionary_validation)
                    validate_sequence_dictionary(output_header, header, path);
                if (!bcf_hdr_merge(output_header, header))
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot merge VCF headers");
            }
            inputs.push_back(input);
            headers.push_back(header);
        }
        bcf_hdr_append(output_header, "##source=fastgatk-gather-vcfs");
        bcf_hdr_append(output_header, "##fastgatk_gather_vcfs_status=prototype-stream-gather");
        for (const auto& comment : options.comments) {
            // Picard/GATK encodes -CO as a VCF meta-line.  Keep the user text
            // verbatim after the key; HTSlib performs the same header-line
            // validation as the native writer for embedded tabs/newlines.
            bcf_hdr_append(output_header, ("##GatherVcfs.comment=" + comment).c_str());
        }
        if (bcf_hdr_sync(output_header) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync GatherVcfs header");
        std::vector<Region> regions;
        regions.reserve(options.regions.size());
        bool first_interval_selector = true;
        for (const auto& text : options.regions) {
            append_interval_input_with_rule(
                text, output_header, regions, interval_input_stats,
                options.interval_set_rule, first_interval_selector);
            first_interval_selector = false;
        }
        normalize_gather_regions(regions);

        output = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
        if (!output) throw std::runtime_error("cannot open GatherVcfs output: " + options.output);
        if (suffix(options.output, ".gz") &&
            hts_set_opt(output, HTS_OPT_COMPRESSION_LEVEL, options.compression_level) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set GatherVcfs compression level");
        if (bcf_hdr_write(output, output_header) != 0)
            throw std::runtime_error("cannot write GatherVcfs header");
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");

        for (std::size_t file_index = 0; file_index < inputs.size(); ++file_index) {
            int read_status = 0;
            while ((read_status = bcf_read(inputs[file_index], headers[file_index], record)) == 0) {
                ++input_records;
                bcf_unpack(record, BCF_UN_ALL);
                if (record->rid < 0 || record->pos < 0)
                    throw std::runtime_error("BAD_INPUT: GatherVcfs encountered an invalid coordinate");
                const char* contig = bcf_hdr_id2name(headers[file_index], record->rid);
                if (!contig || bcf_hdr_name2id(output_header, contig) < 0)
                    throw std::runtime_error("BAD_INPUT: GatherVcfs contig is absent from merged header");
                const int output_rid = bcf_hdr_name2id(output_header, contig);
                const auto current_end_exclusive = record_end_exclusive(headers[file_index], record);
                if (!regions.empty() &&
                    !std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
                        // Interval selection is overlap-based.  This matters
                        // for gVCF reference blocks: a record whose POS is
                        // before -L may still cover the selected locus via
                        // INFO/END.  Coordinates are zero-based half-open in
                        // both Region and record_end_exclusive.
                        return output_rid == region.rid &&
                               current_end_exclusive > region.begin &&
                               record->pos < region.end;
                    })) {
                    ++interval_skipped;
                    bcf_clear(record);
                    continue;
                }
                if (!options.allow_overlaps &&
                    (output_rid < previous_rid ||
                     (output_rid == previous_rid &&
                      (record->pos < previous_pos || record->pos < previous_end_exclusive))))
                    throw std::runtime_error(
                        "BAD_INPUT: GatherVcfs shard records are not coordinate-sorted/disjoint; "
                        "overlapping intervals require --allow-overlaps");
                previous_rid = output_rid;
                previous_pos = record->pos;
                previous_end_exclusive = current_end_exclusive;
                if (bcf_translate(output_header, headers[file_index], record) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot translate gathered record header");
                if (bcf_write(output, output_header, record) != 0)
                    throw std::runtime_error("cannot write GatherVcfs record");
                bcf_clear(record);
            }
            if (read_status < -1)
                throw std::runtime_error("BAD_INPUT: GatherVcfs failed while reading shard " + input_paths[file_index]);
        }
        bcf_destroy(record); record = nullptr;
        bcf_close(output); output = nullptr;
        close_inputs(inputs, headers);

        std::string index_path;
        if (options.create_index && options.output != "-") {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build GatherVcfs index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const bool primary_complete = file_complete(options.output);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !index_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GatherVcfs output or index is missing/empty");
        std::error_code output_error;
        const auto output_bytes = std::filesystem::file_size(options.output, output_error);
        if (output_error || output_bytes == 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GatherVcfs output size is unavailable");
        const auto index_bytes = index_path.empty() ? std::uintmax_t{0}
                                                     : std::filesystem::file_size(index_path, output_error);
        if (!index_path.empty() && (output_error || index_bytes == 0))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: GatherVcfs index size is unavailable");
        const auto wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_begin).count();
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write GatherVcfs manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"GatherVcfs\",\"implementation\":\"fastgatk-gather-vcfs\",\"status\":\"prototype\"," 
                 << "\"execution_space\":\"Host\",\"determinism\":\"strict\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"stream_gather\":true,\"header_merge\":true,\"sample_order_validation\":true,\"coordinate_contract\":"
                 << (options.allow_overlaps ? "false" : "true")
                 << ",\"interval_subset\":" << (!regions.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_input_stats.files
                 << ",\"interval_list_records\":" << interval_input_stats.records
                 << ",\"reorder_input_by_first_variant\":"
                 << (options.reorder_input_by_first_variant ? "true" : "false")
                 << ",\"comments\":" << options.comments.size()
                 << ",\"compression_level\":" << options.compression_level
                 << ",\"sequence_dictionary_validation\":"
                 << (!options.disable_sequence_dictionary_validation ? "true" : "false")
                 << ",\"quiet\":" << (options.quiet ? "true" : "false")
                 << ",\"cloud_uri_staging\":false"
                 << ",\"vcf_index\":" << (index_path.empty() ? "false" : "true")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"vcf\",\"complete\":" << (primary_complete ? "true" : "false") << "}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":" + (index_complete ? "true}" : "false}"))
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size() << ",\"input_records\":" << input_records
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"output_bytes\":" << output_bytes
                 << ",\"index_bytes\":" << index_bytes
                 << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds << "}}\n";
        manifest.close();
        if (!manifest) throw std::runtime_error("cannot finalize GatherVcfs manifest");
        if (options.output != "-") {
            fastgatk::runtime::require_complete_output({
                std::filesystem::path(options.output), std::filesystem::path(index_path),
                std::filesystem::path(manifest_path), !index_path.empty(), true});
        }
        std::cout << "{\"tool\":\"GatherVcfs\",\"status\":\"prototype\",\"execution_space\":\"Host\",\"input_files\":"
                  << options.inputs.size() << ",\"input_records\":" << input_records
                  << ",\"interval_skipped\":" << interval_skipped
                  << ",\"interval_set_rule\":\"" << interval_set_rule_name(options.interval_set_rule) << "\""
                  << ",\"output_bytes\":" << output_bytes
                  << ",\"index_bytes\":" << index_bytes
                  << ",\"wall_seconds\":" << std::setprecision(12) << wall_seconds
                  << ",\"reorder_input_by_first_variant\":"
                  << (options.reorder_input_by_first_variant ? "true" : "false")
                  << ",\"comments\":" << options.comments.size()
                  << ",\"compression_level\":" << options.compression_level
                  << ",\"sequence_dictionary_validation\":"
                  << (!options.disable_sequence_dictionary_validation ? "true" : "false")
                  << ",\"quiet\":" << (options.quiet ? "true" : "false") << "}\n";
        bcf_hdr_destroy(output_header);
        return 0;
    } catch (...) {
        bcf_destroy(record);
        if (output) bcf_close(output);
        close_inputs(inputs, headers);
        if (output_header) bcf_hdr_destroy(output_header);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for GatherVcfs");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-gather-vcfs: " << error.what() << '\n';
        return 2;
    }
}
