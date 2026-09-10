#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include "fastgatk/hdf5_count_collection.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Interval {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::string name;
};

struct Annotation {
    Interval interval;
    double gc = std::numeric_limits<double>::quiet_NaN();
    double mappability = std::numeric_limits<double>::quiet_NaN();
    double segmental_duplication = std::numeric_limits<double>::quiet_NaN();
    bool has_gc = false;
    bool has_mappability = false;
    bool has_segmental_duplication = false;
};

struct CountTable {
    std::vector<Interval> intervals;
    std::vector<double> counts;
    std::map<std::string, std::int64_t> lengths;
    std::vector<std::string> sequence_order;
};

struct Options {
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    std::vector<std::string> count_inputs;
    std::string annotated_intervals;
    std::string output;
    std::string manifest;
    std::string merging_rule = "ALL";
    std::string set_rule = "UNION";
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    double minimum_gc = 0.1;
    double maximum_gc = 0.9;
    double minimum_mappability = 0.9;
    double maximum_mappability = 1.0;
    double minimum_segmental_duplication = 0.0;
    double maximum_segmental_duplication = 0.5;
    double low_count_threshold = 10.0;
    double low_count_percentage = 50.0;
    double extreme_minimum_percentile = 1.0;
    double extreme_maximum_percentile = 99.0;
    double extreme_percentage = 90.0;
    int threads = 1;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

bool suffix_ci(const std::string& value, const char* suffix);

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

double parse_double(int& index, int argc, char** argv, const std::string& argument, const char* name) {
    const auto text = require_value(index, argc, argv, argument, name);
    const auto value = std::stod(text);
    if (!std::isfinite(value)) throw std::invalid_argument(std::string(name) + " must be finite");
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-filter-intervals (GATK-compatible native prototype)\n"
                         "  -L, --intervals REGION           interval or interval-list/BED (repeatable)\n"
                         "      --annotated-intervals FILE   AnnotateIntervals TSV\n"
                         "  -XL, --exclude-intervals REGION  exclusion interval (repeatable)\n"
                         "  -O, --output FILE                filtered Picard interval-list\n"
                         "      --minimum-gc-content X       inclusive GC lower bound (default 0.1)\n"
                         "      --maximum-gc-content X       inclusive GC upper bound (default 0.9)\n"
                         "      --minimum-mappability X      inclusive mappability lower bound\n"
                         "      --maximum-mappability X      inclusive mappability upper bound\n"
                         "      --minimum-segmental-duplication-content X\n"
                         "      --maximum-segmental-duplication-content X\n"
                         "  -I, --input FILE                 TSV/HDF5 count collection (repeatable)\n"
                         "      --low-count-filter-count-threshold N\n"
                         "      --low-count-filter-percentage-of-samples X\n"
                         "      --extreme-count-filter-minimum-percentile X\n"
                         "      --extreme-count-filter-maximum-percentile X\n"
                         "      --extreme-count-filter-percentage-of-samples X\n"
                         "      --threads N                  Kokkos execution threads\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                has_option(argument, "--region") ? "--region" :
                has_option(argument, "--interval") ? "--interval" : "--intervals";
            options.regions.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals")) {
            options.exclusions.push_back(require_value(index, argc, argv, argument,
                "--exclude-intervals", "-XL"));
        } else if (has_option(argument, "--annotated-intervals")) {
            options.annotated_intervals = require_value(index, argc, argv, argument, "--annotated-intervals");
        } else if (argument == "-I" || has_option(argument, "--input")) {
            const auto path = require_value(index, argc, argv, argument, "--input", "-I");
            options.count_inputs.push_back(path);
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "-imr" || has_option(argument, "--interval-merging-rule")) {
            options.merging_rule = require_value(index, argc, argv, argument,
                "--interval-merging-rule", "-imr");
            if (options.merging_rule != "ALL" && options.merging_rule != "OVERLAPPING_ONLY")
                throw std::invalid_argument("--interval-merging-rule must be ALL or OVERLAPPING_ONLY");
        } else if (argument == "-isr" || has_option(argument, "--interval-set-rule")) {
            options.set_rule = require_value(index, argc, argv, argument,
                "--interval-set-rule", "-isr");
            if (options.set_rule != "UNION" && options.set_rule != "INTERSECTION")
                throw std::invalid_argument("--interval-set-rule must be UNION or INTERSECTION");
            // FilterIntervals validates the shared GATK interval argument
            // collection and only permits UNION.  INTERSECTION is advertised
            // by the generic argument help but rejected by this tool in
            // GATK 4.6.2.0; fail closed instead of producing a non-GATK set.
            if (options.set_rule == "INTERSECTION")
                throw std::invalid_argument("Interval set rule must be set to UNION.");
        } else if (argument == "-ip" || has_option(argument, "--interval-padding")) {
            options.interval_padding = std::stoll(require_value(index, argc, argv, argument,
                "--interval-padding", "-ip"));
            if (options.interval_padding < 0) throw std::invalid_argument("--interval-padding must be non-negative");
            // FilterIntervals delegates interval parsing to the copy-number
            // interval argument collection, which explicitly rejects any
            // padding.  Accepting a non-zero value here would silently select
            // different annotation rows from GATK, so fail closed.
            if (options.interval_padding != 0)
                throw std::invalid_argument("Interval padding must be set to 0.");
        } else if (argument == "-ixp" || has_option(argument, "--interval-exclusion-padding")) {
            options.exclusion_padding = std::stoll(require_value(index, argc, argv, argument,
                "--interval-exclusion-padding", "-ixp"));
            if (options.exclusion_padding < 0) throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
            if (options.exclusion_padding != 0)
                throw std::invalid_argument("Interval exclusion padding must be set to 0.");
        } else if (has_option(argument, "--minimum-gc-content")) {
            options.minimum_gc = parse_double(index, argc, argv, argument, "--minimum-gc-content");
        } else if (has_option(argument, "--maximum-gc-content")) {
            options.maximum_gc = parse_double(index, argc, argv, argument, "--maximum-gc-content");
        } else if (has_option(argument, "--minimum-mappability")) {
            options.minimum_mappability = parse_double(index, argc, argv, argument, "--minimum-mappability");
        } else if (has_option(argument, "--maximum-mappability")) {
            options.maximum_mappability = parse_double(index, argc, argv, argument, "--maximum-mappability");
        } else if (has_option(argument, "--minimum-segmental-duplication-content")) {
            options.minimum_segmental_duplication = parse_double(index, argc, argv, argument,
                "--minimum-segmental-duplication-content");
        } else if (has_option(argument, "--maximum-segmental-duplication-content")) {
            options.maximum_segmental_duplication = parse_double(index, argc, argv, argument,
                "--maximum-segmental-duplication-content");
        } else if (has_option(argument, "--low-count-filter-count-threshold")) {
            options.low_count_threshold = parse_double(index, argc, argv, argument,
                "--low-count-filter-count-threshold");
        } else if (has_option(argument, "--low-count-filter-percentage-of-samples")) {
            options.low_count_percentage = parse_double(index, argc, argv, argument,
                "--low-count-filter-percentage-of-samples");
        } else if (has_option(argument, "--extreme-count-filter-minimum-percentile")) {
            options.extreme_minimum_percentile = parse_double(index, argc, argv, argument,
                "--extreme-count-filter-minimum-percentile");
        } else if (has_option(argument, "--extreme-count-filter-maximum-percentile")) {
            options.extreme_maximum_percentile = parse_double(index, argc, argv, argument,
                "--extreme-count-filter-maximum-percentile");
        } else if (has_option(argument, "--extreme-count-filter-percentage-of-samples")) {
            options.extreme_percentage = parse_double(index, argc, argv, argument,
                "--extreme-count-filter-percentage-of-samples");
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
            if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
        } else if (argument == "--quiet" || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates") ||
                   has_option(argument, "--tmp-dir") || has_option(argument, "--java-options") ||
                   has_option(argument, "--gatk-config-file")) {
            if (argument != "--quiet") {
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    argument.rfind("--seconds-between-progress-updates", 0) == 0 ? "--seconds-between-progress-updates" :
                    argument.rfind("--tmp-dir", 0) == 0 ? "--tmp-dir" :
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--gatk-config-file");
            }
        } else if (has_option(argument, "--disable-sequence-dictionary-validation") ||
                   has_option(argument, "--use-jdk-deflater") || has_option(argument, "--use-jdk-inflater") ||
                   has_option(argument, "--lenient")) {
            // These flags do not alter the annotation-only algorithm.
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.regions.empty()) throw std::invalid_argument("-L/--intervals is required");
    if (options.annotated_intervals.empty() && options.count_inputs.empty())
        throw std::invalid_argument("--annotated-intervals or -I/--input is required for native FilterIntervals");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    // GATK 4.6.2.0's FilterIntervals validation requires the explicit
    // OVERLAPPING_ONLY mode (despite the generic help text's ALL default).
    if (options.merging_rule != "OVERLAPPING_ONLY")
        throw std::invalid_argument("Interval merging rule must be set to OVERLAPPING_ONLY.");
    if (options.minimum_gc > options.maximum_gc || options.minimum_mappability > options.maximum_mappability ||
        options.minimum_segmental_duplication > options.maximum_segmental_duplication ||
        options.low_count_threshold < 0.0 || options.low_count_percentage < 0.0 || options.low_count_percentage > 100.0 ||
        options.extreme_minimum_percentile < 0.0 || options.extreme_minimum_percentile > 100.0 ||
        options.extreme_maximum_percentile < 0.0 || options.extreme_maximum_percentile > 100.0 ||
        options.extreme_minimum_percentile > options.extreme_maximum_percentile ||
        options.extreme_percentage < 0.0 || options.extreme_percentage > 100.0)
        throw std::invalid_argument("filter minimum must not exceed maximum");
    return options;
}

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) fields.push_back(std::move(field));
    return fields;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (begin <= line.size()) {
        const auto end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

void remember_contig(std::vector<std::string>& order, const std::string& contig) {
    if (contig.empty()) return;
    if (std::find(order.begin(), order.end(), contig) == order.end()) order.push_back(contig);
}

bool suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    for (std::size_t i = 0; i < ending.size(); ++i) {
        const auto lhs = std::tolower(static_cast<unsigned char>(value[value.size() - ending.size() + i]));
        const auto rhs = std::tolower(static_cast<unsigned char>(ending[i]));
        if (lhs != rhs) return false;
    }
    return true;
}

std::vector<std::string> read_lines(const std::string& path) {
    std::vector<std::string> lines;
    if (suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        kstring_t line{0, 0, nullptr};
        while (hts_getline(file, '\n', &line) >= 0)
            lines.emplace_back(line.s == nullptr ? "" : std::string(line.s, line.l));
        free(line.s);
        if (hts_close(file) != 0) throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
    } else {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        for (std::string line; std::getline(input, line);) lines.push_back(std::move(line));
    }
    return lines;
}

Interval parse_token(const std::string& token) {
    const auto colon = token.find(':');
    const std::string contig = colon == std::string::npos ? token : token.substr(0, colon);
    if (contig.empty()) throw std::invalid_argument("invalid interval contig");
    if (colon == std::string::npos) return {contig, 1, std::numeric_limits<std::int64_t>::max(), {}};
    const auto dash = token.find('-', colon + 1);
    const auto start = std::stoll(token.substr(colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1));
    const auto end = dash == std::string::npos ? start : std::stoll(token.substr(dash + 1));
    if (start < 1 || end < start) throw std::invalid_argument("invalid interval coordinates: " + token);
    return {contig, start, end, {}};
}

std::vector<Interval> parse_selector(const std::string& selector,
                                     std::map<std::string, std::int64_t>& lengths,
                                     std::vector<std::string>& sequence_order) {
    std::error_code error;
    const bool file = std::filesystem::is_regular_file(selector, error) && !error;
    if (!file) return {parse_token(selector)};
    const bool bed = suffix_ci(selector, ".bed") || suffix_ci(selector, ".bed.gz");
    std::vector<Interval> result;
    for (const auto& line : read_lines(selector)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '@') {
            if (first != std::string::npos && line[first] == '@') {
                const auto fields = split_fields(line.substr(first));
                if (fields.size() >= 3 && fields[0] == "@SQ") {
                    std::string name;
                    std::int64_t length = 0;
                    for (std::size_t i = 1; i < fields.size(); ++i) {
                        if (fields[i].rfind("SN:", 0) == 0) name = fields[i].substr(3);
                        if (fields[i].rfind("LN:", 0) == 0) length = std::stoll(fields[i].substr(3));
                    }
                    if (!name.empty() && length > 0) {
                        lengths[name] = length;
                        remember_contig(sequence_order, name);
                    }
                }
            }
            continue;
        }
        const auto fields = split_fields(line.substr(first));
        if (fields.empty() || fields[0] == "track" || fields[0] == "browser") continue;
        if (bed) {
            if (fields.size() < 3) throw std::invalid_argument("BED requires contig/start/end");
            const auto start0 = std::stoll(fields[1]);
            const auto end0 = std::stoll(fields[2]);
            if (start0 < 0 || end0 <= start0) throw std::invalid_argument("invalid BED coordinates");
            result.push_back({fields[0], start0 + 1, end0, fields.size() >= 4 ? fields[3] : std::string{}});
        } else if (fields.size() == 1) {
            result.push_back(parse_token(fields[0]));
        } else if (fields.size() >= 3) {
            result.push_back({fields[0], std::stoll(fields[1]), std::stoll(fields[2]), fields.size() >= 5 ? fields[4] : std::string{}});
        } else {
            throw std::invalid_argument("interval-list row requires contig/start/end");
        }
    }
    if (result.empty()) throw std::invalid_argument("BAD_INPUT: interval file is empty: " + selector);
    return result;
}

std::map<std::string, std::size_t> make_contig_order(
    const std::vector<std::string>& sequence_order,
    const std::map<std::string, std::int64_t>& lengths) {
    std::map<std::string, std::size_t> order;
    std::size_t ordinal = 0;
    for (const auto& contig : sequence_order) {
        if (lengths.find(contig) != lengths.end() && order.emplace(contig, ordinal).second)
            ++ordinal;
    }
    // Inputs without a dictionary (for example a bare literal selector) still
    // need deterministic ordering.  Append their known contigs lexically
    // after the explicitly declared sequence dictionary.
    for (const auto& entry : lengths)
        if (order.emplace(entry.first, ordinal).second) ++ordinal;
    return order;
}

std::vector<Interval> normalize(std::vector<Interval> intervals, const std::string& rule,
                                const std::map<std::string, std::int64_t>& lengths,
                                const std::vector<std::string>& sequence_order) {
    const auto order = make_contig_order(sequence_order, lengths);
    const auto rank = [&](const std::string& contig) {
        const auto found = order.find(contig);
        return found == order.end() ? order.size() + 1 : found->second;
    };
    std::sort(intervals.begin(), intervals.end(), [&](const Interval& lhs, const Interval& rhs) {
        const auto lr = rank(lhs.contig), rr = rank(rhs.contig);
        return lr != rr ? lr < rr : (lhs.contig != rhs.contig ? lhs.contig < rhs.contig :
               (lhs.start != rhs.start ? lhs.start < rhs.start : lhs.end < rhs.end));
    });
    std::vector<Interval> merged;
    for (auto interval : intervals) {
        if (interval.end == std::numeric_limits<std::int64_t>::max()) {
            const auto length = lengths.find(interval.contig);
            if (length != lengths.end()) interval.end = length->second;
        }
        if (interval.end < interval.start) throw std::invalid_argument("BAD_INPUT: interval end precedes start");
        if (!merged.empty() && merged.back().contig == interval.contig &&
            (interval.start <= merged.back().end ||
             (rule == "ALL" && interval.start == merged.back().end + 1))) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(std::move(interval));
        }
    }
    return merged;
}

std::vector<Interval> intersect_intervals(
    const std::vector<Interval>& lhs,
    const std::vector<Interval>& rhs,
    const std::map<std::string, std::int64_t>& lengths,
    const std::vector<std::string>& sequence_order) {
    // normalize() orders known contigs by the sequence dictionary and only
    // then falls back to lexical order.  Keep exactly that order while doing
    // the two-pointer intersection; comparing contig strings directly is
    // wrong for dictionaries such as chr1, chr2, chr10.
    const auto order = make_contig_order(sequence_order, lengths);
    const auto rank = [&](const std::string& contig) {
        const auto found = order.find(contig);
        return found == order.end() ? order.size() + 1 : found->second;
    };
    const auto before = [&](const Interval& left, const Interval& right) {
        const auto left_rank = rank(left.contig);
        const auto right_rank = rank(right.contig);
        return left_rank != right_rank ? left_rank < right_rank : left.contig < right.contig;
    };
    std::vector<Interval> result;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < lhs.size() && right < rhs.size()) {
        if (lhs[left].contig != rhs[right].contig) {
            if (before(lhs[left], rhs[right])) ++left;
            else ++right;
            continue;
        }
        if (lhs[left].end < rhs[right].start) {
            ++left;
            continue;
        }
        if (rhs[right].end < lhs[left].start) {
            ++right;
            continue;
        }
        const auto start = std::max(lhs[left].start, rhs[right].start);
        const auto end = std::min(lhs[left].end, rhs[right].end);
        if (start <= end) {
            result.push_back({lhs[left].contig, start, end,
                              lhs[left].name.empty() ? rhs[right].name : lhs[left].name});
        }
        if (lhs[left].end < rhs[right].end) ++left;
        else if (rhs[right].end < lhs[left].end) ++right;
        else { ++left; ++right; }
    }
    return result;
}

std::string interval_key(const Interval& interval) {
    return interval.contig + "\t" + std::to_string(interval.start) + "\t" + std::to_string(interval.end);
}

struct AnnotatedTable {
    std::vector<std::string> headers;
    std::vector<Annotation> rows;
    std::map<std::string, std::int64_t> lengths;
    std::vector<std::string> sequence_order;
};

AnnotatedTable read_annotations(const std::string& path) {
    AnnotatedTable table;
    bool header_seen = false;
    std::size_t gc_column = std::numeric_limits<std::size_t>::max();
    std::size_t map_column = std::numeric_limits<std::size_t>::max();
    std::size_t dup_column = std::numeric_limits<std::size_t>::max();
    for (const auto& line : read_lines(path)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '@') {
            const auto fields = split_fields(line.substr(first));
            if (fields.size() >= 3 && fields[0] == "@SQ") {
                table.headers.push_back(line.substr(first));
                std::string name;
                std::int64_t length = 0;
                for (std::size_t i = 1; i < fields.size(); ++i) {
                    if (fields[i].rfind("SN:", 0) == 0) name = fields[i].substr(3);
                    if (fields[i].rfind("LN:", 0) == 0) length = std::stoll(fields[i].substr(3));
                }
                if (!name.empty() && length > 0) {
                    table.lengths[name] = length;
                    remember_contig(table.sequence_order, name);
                }
            }
            continue;
        }
        const auto fields = split_tab(line.substr(first));
        if (!header_seen) {
            if (fields.size() < 3 || fields[0] != "CONTIG" || fields[1] != "START" || fields[2] != "END")
                throw std::invalid_argument("BAD_INPUT: annotated intervals must start with CONTIG/START/END");
            for (std::size_t i = 3; i < fields.size(); ++i) {
                if (fields[i] == "GC_CONTENT") gc_column = i;
                else if (fields[i] == "MAPPABILITY") map_column = i;
                else if (fields[i] == "SEGMENTAL_DUPLICATION_CONTENT") dup_column = i;
            }
            if (gc_column == std::numeric_limits<std::size_t>::max() &&
                map_column == std::numeric_limits<std::size_t>::max() &&
                dup_column == std::numeric_limits<std::size_t>::max())
                throw std::invalid_argument("BAD_INPUT: annotated intervals contain no supported annotation columns");
            header_seen = true;
            continue;
        }
        if (fields.size() < 3) throw std::invalid_argument("BAD_INPUT: malformed annotated interval row");
        if (gc_column >= fields.size() && gc_column != std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("BAD_INPUT: GC_CONTENT column missing in annotated row");
        if (map_column >= fields.size() && map_column != std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("BAD_INPUT: MAPPABILITY column missing in annotated row");
        if (dup_column >= fields.size() && dup_column != std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("BAD_INPUT: SEGMENTAL_DUPLICATION_CONTENT column missing in annotated row");
        Annotation annotation;
        annotation.interval = {fields[0], std::stoll(fields[1]), std::stoll(fields[2]), {}};
        const auto parse_annotation = [&](std::size_t column, bool& present) {
            if (column == std::numeric_limits<std::size_t>::max()) return std::numeric_limits<double>::quiet_NaN();
            present = true;
            const auto value = fields[column];
            if (value == "NA" || value == "NaN" || value == ".") return std::numeric_limits<double>::quiet_NaN();
            return std::stod(value);
        };
        annotation.gc = parse_annotation(gc_column, annotation.has_gc);
        annotation.mappability = parse_annotation(map_column, annotation.has_mappability);
        annotation.segmental_duplication = parse_annotation(dup_column, annotation.has_segmental_duplication);
        if (annotation.interval.start < 1 || annotation.interval.end < annotation.interval.start)
            throw std::invalid_argument("BAD_INPUT: invalid annotated interval coordinates");
        table.rows.push_back(std::move(annotation));
    }
    if (!header_seen || table.rows.empty()) throw std::invalid_argument("BAD_INPUT: annotated intervals contain no rows");
    return table;
}

CountTable read_count_table(const std::string& path) {
    CountTable table;
    if (suffix_ci(path, ".h5") || suffix_ci(path, ".hdf5")) {
        if (!fastgatk::hdf5::available())
            throw std::invalid_argument("BACKEND_UNAVAILABLE: HDF5 count-based FilterIntervals requires HDF5 support; use --fallback");
        const auto collection = fastgatk::hdf5::read_simple_count_collection(path);
        table.intervals.reserve(collection.intervals.size());
        table.counts = collection.counts;
        for (const auto& interval : collection.intervals)
            table.intervals.push_back({interval.contig, interval.start, interval.end, {}});
        std::istringstream dictionary(collection.sequence_dictionary);
        for (std::string line; std::getline(dictionary, line);) {
            const auto fields = split_fields(line);
            if (fields.empty() || fields[0] != "@SQ") continue;
            std::string name;
            std::int64_t length = 0;
            for (std::size_t i = 1; i < fields.size(); ++i) {
                if (fields[i].rfind("SN:", 0) == 0) name = fields[i].substr(3);
                if (fields[i].rfind("LN:", 0) == 0) length = std::stoll(fields[i].substr(3));
            }
            if (!name.empty() && length > 0) {
                table.lengths[name] = length;
                remember_contig(table.sequence_order, name);
            }
        }
        if (table.intervals.empty() || table.intervals.size() != table.counts.size())
            throw std::invalid_argument("BAD_INPUT: HDF5 count collection contains no aligned intervals: " + path);
        for (std::size_t i = 0; i < table.intervals.size(); ++i) {
            if (table.intervals[i].start < 1 || table.intervals[i].end < table.intervals[i].start ||
                !std::isfinite(table.counts[i]) || table.counts[i] < 0.0)
                throw std::invalid_argument("BAD_INPUT: invalid HDF5 count row: " + path);
        }
        return table;
    }
    bool header_seen = false;
    for (const auto& line : read_lines(path)) {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        const auto text = line.substr(first);
        if (text[0] == '@') {
            const auto fields = split_fields(text);
            if (fields.size() >= 3 && fields[0] == "@SQ") {
                std::string name;
                std::int64_t length = 0;
                for (std::size_t i = 1; i < fields.size(); ++i) {
                    if (fields[i].rfind("SN:", 0) == 0) name = fields[i].substr(3);
                    if (fields[i].rfind("LN:", 0) == 0) length = std::stoll(fields[i].substr(3));
                }
                if (!name.empty() && length > 0) {
                    table.lengths[name] = length;
                    remember_contig(table.sequence_order, name);
                }
            }
            continue;
        }
        const auto fields = split_tab(text);
        if (!header_seen) {
            if (fields.size() != 4 || fields[0] != "CONTIG" || fields[1] != "START" ||
                fields[2] != "END" || fields[3] != "COUNT")
                throw std::invalid_argument("BAD_INPUT: count table must contain CONTIG/START/END/COUNT: " + path);
            header_seen = true;
            continue;
        }
        if (fields.size() != 4) throw std::invalid_argument("BAD_INPUT: malformed count row: " + path);
        Interval interval{fields[0], std::stoll(fields[1]), std::stoll(fields[2]), {}};
        const auto count = std::stod(fields[3]);
        if (interval.contig.empty() || interval.start < 1 || interval.end < interval.start ||
            !std::isfinite(count) || count < 0.0)
            throw std::invalid_argument("BAD_INPUT: invalid count row: " + path);
        table.intervals.push_back(std::move(interval));
        table.counts.push_back(count);
    }
    if (!header_seen || table.intervals.empty() || table.intervals.size() != table.counts.size())
        throw std::invalid_argument("BAD_INPUT: count table contains no rows: " + path);
    return table;
}

bool overlap(const Interval& lhs, const Interval& rhs) {
    return lhs.contig == rhs.contig && lhs.start <= rhs.end && rhs.start <= lhs.end;
}

double percentile(std::vector<double> values, double percent) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    if (percent <= 0.0) return values.front();
    if (percent >= 100.0) return values.back();
    // Apache Commons Math Percentile's legacy estimator (used by GATK 4.6):
    // position = p * (n + 1) / 100, clamped to the sample endpoints.
    const double position = (percent / 100.0) * static_cast<double>(values.size() + 1);
    if (position <= 1.0) return values.front();
    if (position >= static_cast<double>(values.size())) return values.back();
    const auto lower = static_cast<std::size_t>(std::floor(position)) - 1;
    const auto upper = lower + 1;
    return values[lower] + (position - std::floor(position)) * (values[upper] - values[lower]);
}

std::string json_escape(const std::string& text) {
    std::ostringstream out;
    for (char c : text) {
        if (c == '"' || c == '\\') out << '\\';
        if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else out << c;
    }
    return out.str();
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    const bool have_annotations = !options.annotated_intervals.empty();
    AnnotatedTable annotations;
    if (have_annotations) annotations = read_annotations(options.annotated_intervals);
    std::vector<CountTable> count_tables;
    count_tables.reserve(options.count_inputs.size());
    for (const auto& path : options.count_inputs) count_tables.push_back(read_count_table(path));
    std::vector<Interval> requested;
    std::vector<std::vector<Interval>> requested_sets;
    auto lengths = annotations.lengths;
    auto sequence_order = annotations.sequence_order;
    if (lengths.empty()) {
        for (const auto& table : count_tables) {
            for (const auto& contig : table.sequence_order) remember_contig(sequence_order, contig);
            for (const auto& entry : table.lengths) lengths.emplace(entry.first, entry.second);
        }
    } else {
        // If the annotation file has no dictionary, inherit one from the
        // count collection before appending selector-declared contigs.
        for (const auto& table : count_tables)
            for (const auto& contig : table.sequence_order) remember_contig(sequence_order, contig);
    }
    for (const auto& selector : options.regions) {
        auto intervals = parse_selector(selector, lengths, sequence_order);
        intervals = normalize(std::move(intervals), options.merging_rule, lengths, sequence_order);
        requested_sets.push_back(std::move(intervals));
    }
    std::vector<Interval> intervals;
    if (options.set_rule == "INTERSECTION" && !requested_sets.empty()) {
        intervals = requested_sets.front();
        for (std::size_t index = 1; index < requested_sets.size() && !intervals.empty(); ++index)
            intervals = intersect_intervals(intervals, requested_sets[index], lengths, sequence_order);
    } else {
        for (const auto& set : requested_sets)
            requested.insert(requested.end(), set.begin(), set.end());
        intervals = normalize(std::move(requested), options.merging_rule, lengths, sequence_order);
    }
    if (intervals.empty()) throw std::invalid_argument("BAD_INPUT: no intervals remain after normalization");
    if (options.interval_padding > 0) {
        for (auto& interval : intervals) {
            interval.start = std::max<std::int64_t>(1, interval.start - options.interval_padding);
            const auto found = lengths.find(interval.contig);
            const auto limit = found == lengths.end() ? interval.end + options.interval_padding : found->second;
            interval.end = std::min(limit, interval.end + options.interval_padding);
        }
        intervals = normalize(std::move(intervals), options.merging_rule, lengths, sequence_order);
    }
    std::vector<Interval> exclusions;
    for (const auto& selector : options.exclusions) {
        const auto parsed = parse_selector(selector, lengths, sequence_order);
        exclusions.insert(exclusions.end(), parsed.begin(), parsed.end());
    }
    exclusions = normalize(std::move(exclusions), options.merging_rule, lengths, sequence_order);
    if (options.exclusion_padding > 0) {
        for (auto& interval : exclusions) {
            interval.start = std::max<std::int64_t>(1, interval.start - options.exclusion_padding);
            const auto found = lengths.find(interval.contig);
            const auto limit = found == lengths.end() ? interval.end + options.exclusion_padding : found->second;
            interval.end = std::min(limit, interval.end + options.exclusion_padding);
        }
        exclusions = normalize(std::move(exclusions), options.merging_rule, lengths, sequence_order);
    }

    std::unordered_map<std::string, Annotation> by_key;
    by_key.reserve(annotations.rows.size());
    for (const auto& row : annotations.rows) by_key.emplace(interval_key(row.interval), row);
    std::vector<Annotation> selected;
    selected.reserve(intervals.size());
    for (const auto& interval : intervals) {
        if (have_annotations) {
            const auto found = by_key.find(interval_key(interval));
            if (found == by_key.end()) continue;
            selected.push_back(found->second);
        } else {
            selected.push_back(Annotation{interval});
        }
    }
    if (selected.empty()) throw std::invalid_argument("BAD_INPUT: no requested intervals are present in the selected annotation/count inputs");

    std::vector<std::vector<double>> count_matrix;
    if (!count_tables.empty()) {
        count_matrix.assign(count_tables.size(), std::vector<double>(selected.size(), 0.0));
        for (std::size_t sample = 0; sample < count_tables.size(); ++sample) {
            std::unordered_map<std::string, double> values;
            values.reserve(count_tables[sample].intervals.size());
            for (std::size_t row = 0; row < count_tables[sample].intervals.size(); ++row) {
                const auto inserted = values.emplace(interval_key(count_tables[sample].intervals[row]),
                                                      count_tables[sample].counts[row]);
                if (!inserted.second) throw std::invalid_argument("BAD_INPUT: duplicate interval in count table: " + options.count_inputs[sample]);
            }
            for (std::size_t interval = 0; interval < selected.size(); ++interval) {
                const auto found = values.find(interval_key(selected[interval].interval));
                if (found == values.end())
                    throw std::invalid_argument("BAD_INPUT: count table does not contain requested interval: " + options.count_inputs[sample]);
                count_matrix[sample][interval] = found->second;
            }
        }
    }

    const auto annotation_pass = [&](std::size_t index) {
        if (std::any_of(exclusions.begin(), exclusions.end(), [&](const Interval& exclusion) {
                return overlap(selected[index].interval, exclusion);
            })) return false;
        if (selected[index].has_gc &&
            !(std::isfinite(selected[index].gc) && selected[index].gc >= options.minimum_gc && selected[index].gc <= options.maximum_gc))
            return false;
        if (selected[index].has_mappability &&
            !(std::isfinite(selected[index].mappability) && selected[index].mappability >= options.minimum_mappability && selected[index].mappability <= options.maximum_mappability))
            return false;
        if (selected[index].has_segmental_duplication &&
            !(std::isfinite(selected[index].segmental_duplication) && selected[index].segmental_duplication >= options.minimum_segmental_duplication && selected[index].segmental_duplication <= options.maximum_segmental_duplication))
            return false;
        return true;
    };
    std::vector<std::uint8_t> count_failed(selected.size(), 0);
    if (!count_matrix.empty()) {
        const auto sample_count = count_matrix.size();
        std::vector<std::uint8_t> low_failed(selected.size(), 0);
        for (std::size_t interval = 0; interval < selected.size(); ++interval) {
            std::size_t low_count = 0;
            for (const auto& sample : count_matrix)
                if (sample[interval] < options.low_count_threshold) ++low_count;
            if (annotation_pass(interval) &&
                static_cast<double>(low_count) > options.low_count_percentage * static_cast<double>(sample_count) / 100.0)
                low_failed[interval] = 1;
        }
        std::vector<std::size_t> extreme_failures(selected.size(), 0);
        for (std::size_t sample = 0; sample < sample_count; ++sample) {
            std::vector<double> available;
            for (std::size_t interval = 0; interval < selected.size(); ++interval)
                if (annotation_pass(interval) && !low_failed[interval]) available.push_back(count_matrix[sample][interval]);
            const auto minimum = options.extreme_minimum_percentile == 0.0
                ? 0.0 : percentile(available, options.extreme_minimum_percentile);
            const auto maximum = options.extreme_maximum_percentile == 0.0
                ? 0.0 : percentile(available, options.extreme_maximum_percentile);
            for (std::size_t interval = 0; interval < selected.size(); ++interval) {
                if (!annotation_pass(interval) || low_failed[interval]) continue;
                if (!(minimum <= count_matrix[sample][interval] && count_matrix[sample][interval] <= maximum))
                    ++extreme_failures[interval];
            }
        }
        for (std::size_t interval = 0; interval < selected.size(); ++interval)
            count_failed[interval] = low_failed[interval] ||
                static_cast<double>(extreme_failures[interval]) > options.extreme_percentage * static_cast<double>(sample_count) / 100.0;
    }

    std::vector<std::uint8_t> keep(selected.size(), 0);
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(options.threads);
    Kokkos::initialize(settings);
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch filter_host("filter-intervals-mask-v1");
    filter_host.records = selected.size();
    filter_host.bytes = selected.size() * (3 * sizeof(double) + 4 * sizeof(std::uint8_t));
    fastgatk::core::KernelPlan<ExecSpace> filter_plan("filter-intervals-mask");
    filter_plan.begin_prepare(filter_host);
    {
        Kokkos::View<double*> gc("filter_intervals_gc", selected.size());
        Kokkos::View<double*> mappability("filter_intervals_mappability", selected.size());
        Kokkos::View<double*> duplication("filter_intervals_duplication", selected.size());
        Kokkos::View<std::uint8_t*> excluded("filter_intervals_excluded", selected.size());
        Kokkos::View<std::uint8_t*> count_failed_view("filter_intervals_count_failed", selected.size());
        Kokkos::View<std::uint8_t*> has_gc("filter_intervals_has_gc", selected.size());
        Kokkos::View<std::uint8_t*> has_mappability("filter_intervals_has_mappability", selected.size());
        Kokkos::View<std::uint8_t*> has_duplication("filter_intervals_has_duplication", selected.size());
        auto host_gc = Kokkos::create_mirror_view(gc);
        auto host_mappability = Kokkos::create_mirror_view(mappability);
        auto host_duplication = Kokkos::create_mirror_view(duplication);
        auto host_excluded = Kokkos::create_mirror_view(excluded);
        auto host_count_failed = Kokkos::create_mirror_view(count_failed_view);
        auto host_has_gc = Kokkos::create_mirror_view(has_gc);
        auto host_has_mappability = Kokkos::create_mirror_view(has_mappability);
        auto host_has_duplication = Kokkos::create_mirror_view(has_duplication);
        for (std::size_t i = 0; i < selected.size(); ++i) {
            host_gc(i) = selected[i].gc;
            host_mappability(i) = selected[i].mappability;
            host_duplication(i) = selected[i].segmental_duplication;
            host_has_gc(i) = selected[i].has_gc ? 1 : 0;
            host_has_mappability(i) = selected[i].has_mappability ? 1 : 0;
            host_has_duplication(i) = selected[i].has_segmental_duplication ? 1 : 0;
            host_count_failed(i) = count_failed[i];
            bool is_excluded = false;
            for (const auto& exclusion : exclusions) {
                if (overlap(selected[i].interval, exclusion)) {
                    is_excluded = true;
                    break;
                }
            }
            host_excluded(i) = is_excluded ? 1 : 0;
        }
        Kokkos::deep_copy(gc, host_gc);
        Kokkos::deep_copy(mappability, host_mappability);
        Kokkos::deep_copy(duplication, host_duplication);
        Kokkos::deep_copy(excluded, host_excluded);
        Kokkos::deep_copy(count_failed_view, host_count_failed);
        Kokkos::deep_copy(has_gc, host_has_gc);
        Kokkos::deep_copy(has_mappability, host_has_mappability);
        Kokkos::deep_copy(has_duplication, host_has_duplication);
        Kokkos::View<std::uint8_t*> device_keep("filter_intervals_keep", selected.size());
        fastgatk::core::DeviceBatch<ExecSpace> filter_device(selected.size());
        filter_device.bind("gc", gc);
        filter_device.bind("mappability", mappability);
        filter_device.bind("duplication", duplication);
        filter_device.bind("excluded", excluded);
        filter_device.bind("count_failed", count_failed_view);
        filter_device.bind("has_gc", has_gc);
        filter_device.bind("has_mappability", has_mappability);
        filter_device.bind("has_duplication", has_duplication);
        filter_device.bind("keep", device_keep);
        ExecSpace().fence();
        filter_plan.end_prepare(filter_device);
        const auto minimum_gc = options.minimum_gc;
        const auto maximum_gc = options.maximum_gc;
        const auto minimum_mappability = options.minimum_mappability;
        const auto maximum_mappability = options.maximum_mappability;
        const auto minimum_duplication = options.minimum_segmental_duplication;
        const auto maximum_duplication = options.maximum_segmental_duplication;
        filter_plan.begin_execute();
        Kokkos::parallel_for("filter_intervals_annotation_mask", Kokkos::RangePolicy<ExecSpace>(0, selected.size()),
            KOKKOS_LAMBDA(const std::size_t index) {
                bool pass = excluded(index) == 0 && count_failed_view(index) == 0;
                const auto gc_value = gc(index);
                const auto map_value = mappability(index);
                const auto dup_value = duplication(index);
                if (pass && has_gc(index)) pass = !Kokkos::isnan(gc_value) && gc_value >= minimum_gc && gc_value <= maximum_gc;
                if (pass && has_mappability(index)) pass = !Kokkos::isnan(map_value) && map_value >= minimum_mappability && map_value <= maximum_mappability;
                if (pass && has_duplication(index)) pass = !Kokkos::isnan(dup_value) && dup_value >= minimum_duplication && dup_value <= maximum_duplication;
                device_keep(index) = pass ? 1 : 0;
            });
        ExecSpace().fence();
        filter_plan.end_execute();
        auto host_keep = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_keep);
        for (std::size_t i = 0; i < selected.size(); ++i) keep[i] = host_keep(i);
    }
    Kokkos::finalize();

    // Apply the annotation predicates on Host as well to preserve optional-column presence.
    for (std::size_t i = 0; i < selected.size(); ++i) {
        bool pass = keep[i] != 0;
        if (pass && selected[i].has_gc)
            pass = std::isfinite(selected[i].gc) && selected[i].gc >= options.minimum_gc && selected[i].gc <= options.maximum_gc;
        if (pass && selected[i].has_mappability)
            pass = std::isfinite(selected[i].mappability) && selected[i].mappability >= options.minimum_mappability && selected[i].mappability <= options.maximum_mappability;
        if (pass && selected[i].has_segmental_duplication)
            pass = std::isfinite(selected[i].segmental_duplication) && selected[i].segmental_duplication >= options.minimum_segmental_duplication && selected[i].segmental_duplication <= options.maximum_segmental_duplication;
        keep[i] = pass ? 1 : 0;
    }

    std::map<std::string, std::size_t> contig_counts;
    for (std::size_t i = 0; i < selected.size(); ++i) if (keep[i]) ++contig_counts[selected[i].interval.contig];
    std::size_t kept = 0;
    for (const auto value : keep) kept += value != 0;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        if (keep[i] && contig_counts[selected[i].interval.contig] == 1) keep[i] = 0;
    }
    kept = 0;
    for (const auto value : keep) kept += value != 0;
    if (kept == 0) throw std::runtime_error("BAD_INPUT: filtering removed all intervals");

    std::ofstream output(options.output);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write filtered intervals: " + options.output);
    output << "@HD\tVN:1.6\n";
    for (const auto& header : annotations.headers) output << header << '\n';
    if (annotations.headers.empty()) {
        std::set<std::string> emitted;
        for (const auto& contig : sequence_order) {
            const auto found = lengths.find(contig);
            if (found != lengths.end()) {
                output << "@SQ\tSN:" << contig << "\tLN:" << found->second << '\n';
                emitted.insert(contig);
            }
        }
        for (const auto& [contig, length] : lengths) {
            if (emitted.insert(contig).second)
                output << "@SQ\tSN:" << contig << "\tLN:" << length << '\n';
        }
    }
    for (std::size_t i = 0; i < selected.size(); ++i)
        if (keep[i]) output << selected[i].interval.contig << '\t' << selected[i].interval.start << '\t'
                            << selected[i].interval.end << "\t+\t" << (selected[i].interval.name.empty() ? "interval_" + std::to_string(i) : selected[i].interval.name) << '\n';
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize filtered intervals");
    if (!options.manifest.empty()) {
        std::ofstream manifest(options.manifest);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
        manifest << "{\"schema_version\":1,\"tool\":\"FilterIntervals\",\"implementation\":\"fastgatk-filter-intervals\",\"status\":\"prototype\",\"input_intervals\":"
                 << selected.size() << ",\"output_intervals\":" << kept << ",\"execution_space\":\""
                 << Kokkos::DefaultExecutionSpace::name() << "\",\"interval_set_rule\":\""
                 << options.set_rule << "\",\"interval_merging_rule\":\""
                 << options.merging_rule << "\",\"interval_padding\":"
                 << options.interval_padding << ",\"interval_exclusion_padding\":"
                 << options.exclusion_padding << ",\"excluded_intervals\":"
                 << exclusions.size() << ",\"output\":\"" << json_escape(options.output)
                 << "\",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\"" << ExecSpace::name()
                 << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_batches\":1,\"kernel_records\":" << selected.size()
                 << ",\"kernel_prepare_seconds\":" << filter_plan.telemetry().prepare_seconds
                 << ",\"kernel_execute_seconds\":" << filter_plan.telemetry().execute_seconds << "}}\n";
    }
    std::cout << "{\"tool\":\"FilterIntervals\",\"status\":\"prototype\",\"input_intervals\":"
              << selected.size() << ",\"output_intervals\":" << kept << ",\"kernel_execution_space\":\"" << ExecSpace::name()
              << "\",\"kernel_records\":" << selected.size()
              << ",\"kernel_execute_seconds\":" << filter_plan.telemetry().execute_seconds
              << ",\"output\":\"" << json_escape(options.output)
              << "\",\"interval_set_rule\":\"" << options.set_rule
              << "\",\"interval_merging_rule\":\"" << options.merging_rule
              << "\",\"interval_padding\":" << options.interval_padding
              << ",\"interval_exclusion_padding\":" << options.exclusion_padding
              << ",\"excluded_intervals\":" << exclusions.size() << "}\n";
    return 0;
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
