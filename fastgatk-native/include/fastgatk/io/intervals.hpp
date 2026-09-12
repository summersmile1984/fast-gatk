#pragma once

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>
#include "fastgatk/io/hts_read_guard.hpp"

namespace fastgatk::io {

// Host-side interval representation shared by VCF tools.  Coordinates are
// zero-based, half-open; only the Host parser knows whether an input line was
// a literal/interval-list (1-based inclusive) or BED (0-based half-open).
struct IndexedInterval {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
    // Canonical contig name, retained even when a caller allows a selector
    // to be absent from the current header (for disjoint shard dictionaries).
    std::string contig;
};

struct IntervalFileStats {
    std::size_t files = 0;
    std::size_t records = 0;
};

inline bool has_suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    for (std::size_t index = 0; index < ending.size(); ++index) {
        const auto lhs = static_cast<char>(std::tolower(static_cast<unsigned char>(
            value[value.size() - ending.size() + index])));
        const auto rhs = static_cast<char>(std::tolower(static_cast<unsigned char>(ending[index])));
        if (lhs != rhs) return false;
    }
    return true;
}

inline IndexedInterval parse_interval_token(const std::string& text,
                                            const bcf_hdr_t* header,
                                            bool bed_coordinates = false,
                                            bool allow_unknown_contig = false) {
    IndexedInterval interval;
    const auto colon = text.find(':');
    const auto contig = colon == std::string::npos ? text : text.substr(0, colon);
    interval.rid = bcf_hdr_name2id(header, contig.c_str());
    interval.contig = contig;
    if (interval.rid < 0 && !allow_unknown_contig)
        throw std::invalid_argument("BAD_INPUT: interval contig not in VCF header: " + contig);
    if (colon == std::string::npos) return interval;
    const auto dash = text.find('-', colon + 1);
    const auto begin = std::stoll(text.substr(
        colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1));
    if (begin < (bed_coordinates ? 0 : 1) || begin > std::numeric_limits<int>::max())
        throw std::invalid_argument("invalid interval start");
    interval.begin = static_cast<int>(bed_coordinates ? begin : begin - 1);
    if (dash != std::string::npos) {
        const auto end = std::stoll(text.substr(dash + 1));
        if (end < begin || end > std::numeric_limits<int>::max() ||
            (bed_coordinates && end == begin))
            throw std::invalid_argument("invalid interval end");
        interval.end = static_cast<int>(end);
    }
    return interval;
}

inline std::vector<std::string> split_interval_line(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while (stream >> field) fields.push_back(field);
    return fields;
}

// Interval lists are frequently staged as bgzip/gzip files by Nextflow.  Use
// HTSlib's text reader for both plain and compressed inputs so all VCF tools
// share the same line semantics and do not silently reinterpret a .bed.gz as
// a 1-based Picard interval list.
inline std::vector<std::string> read_interval_lines(const std::string& path) {
    std::vector<std::string> lines;
    if (has_suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval list: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (fastgatk::io::read_text_line(file, &line, path) >= 0)
                lines.emplace_back(line.s == nullptr ? "" : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(file);
            throw;
        }
        const auto close_status = hts_close(file);
        free(line.s);
        if (close_status != 0) throw std::runtime_error("BAD_INPUT: failed reading interval list: " + path);
        return lines;
    }
    std::ifstream input(path);
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open interval list: " + path);
    for (std::string line; std::getline(input, line);) lines.push_back(std::move(line));
    return lines;
}

inline void append_interval_selector(const std::string& selector,
                                     const bcf_hdr_t* header,
                                     std::vector<IndexedInterval>& intervals,
                                     IntervalFileStats& stats,
                                     bool allow_unknown_contigs = false) {
    std::error_code error;
    const bool regular = std::filesystem::is_regular_file(selector, error) && !error;
    const bool variant_file = has_suffix_ci(selector, ".vcf") || has_suffix_ci(selector, ".vcf.gz") ||
                              has_suffix_ci(selector, ".bcf") || has_suffix_ci(selector, ".bcf.gz");
    if (!regular || variant_file) {
        intervals.push_back(parse_interval_token(selector, header, false, allow_unknown_contigs));
        return;
    }

    const bool bed = has_suffix_ci(selector, ".bed") || has_suffix_ci(selector, ".bed.gz");
    ++stats.files;
    std::size_t records = 0;
    const auto lines = read_interval_lines(selector);
    for (std::size_t line_number = 0; line_number < lines.size(); ++line_number) {
        const auto& line = lines[line_number];
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '@') continue;
        const auto fields = split_interval_line(line.substr(first));
        if (fields.empty()) continue;
        if (fields.front() == "track" || fields.front() == "browser") continue;
        try {
            if (bed) {
                if (fields.size() < 3)
                    throw std::invalid_argument("BED requires contig/start/end");
                const auto start = std::stoll(fields[1]);
                const auto end = std::stoll(fields[2]);
                if (start < 0 || end <= start || end > std::numeric_limits<int>::max())
                    throw std::invalid_argument("invalid BED coordinates");
                const auto interval = parse_interval_token(
                    fields[0] + ":" + std::to_string(start) + "-" + std::to_string(end),
                    header, true, allow_unknown_contigs);
                intervals.push_back(interval);
            } else {
                if (fields.size() == 2)
                    throw std::invalid_argument("interval list requires contig/start/end");
                const auto interval = fields.size() == 1
                    ? parse_interval_token(fields[0], header, false, allow_unknown_contigs)
                    : parse_interval_token(fields[0] + ":" + fields[1] + "-" + fields[2],
                                           header, false, allow_unknown_contigs);
                intervals.push_back(interval);
            }
        } catch (const std::exception& exception) {
            throw std::invalid_argument("BAD_INPUT: malformed interval list line " +
                                        std::to_string(line_number + 1) + " in " + selector +
                                        ": " + exception.what());
        }
        ++records;
    }
    if (records == 0) throw std::invalid_argument("BAD_INPUT: interval list is empty: " + selector);
    stats.records += records;
}

}  // namespace fastgatk::io
