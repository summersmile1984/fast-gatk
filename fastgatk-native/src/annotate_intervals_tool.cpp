#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include "fastgatk/io/hts_read_guard.hpp"

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Interval {
    std::string contig;
    std::int64_t start = 0; // one-based inclusive
    std::int64_t end = 0;   // one-based inclusive
};

struct TrackSegment {
    int tid = -1;
    std::int64_t start = 0; // zero-based inclusive
    std::int64_t end = 0;   // zero-based exclusive
    double score = 1.0;
};

struct Options {
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    std::string output;
    std::string manifest;
    std::string mappability_track;
    std::string segmental_duplication_track;
    // AnnotateIntervals is an interval walker.  Keep the generic GATK
    // selector controls on the native path instead of rejecting otherwise
    // ordinary CNV target lists.
    // CopyNumberArgumentValidationUtils requires this tool's selector
    // collection to be OVERLAPPING_ONLY (the shared engine default is ALL).
    std::string interval_merging_rule = "ALL";
    std::string interval_set_rule = "UNION";
    std::int64_t interval_padding = 0;
    std::int64_t exclusion_padding = 0;
    // GATK exposes this as a FeatureManager query-cache tuning knob.  It does
    // not alter annotation values, but accepting and recording it is required
    // for direct CLI replacement (including zero/negative values accepted by
    // the pinned 4.6.2.0 implementation).
    std::int64_t feature_query_lookahead = 1'000'000;
    int threads = 1;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-annotate-intervals (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE             reference FASTA (+ .fai)\n"
                         "  -L, --intervals REGION           interval or interval-list/BED (repeatable)\n"
                         "  -XL, --exclude-intervals REGION  subtract interval or interval-list/BED (repeatable)\n"
                         "      --interval-merging-rule RULE ALL (default) or OVERLAPPING_ONLY\n"
                         "      --interval-set-rule RULE      UNION (default) or INTERSECTION\n"
                         "      --interval-padding N          pad included intervals\n"
                         "      --interval-exclusion-padding N pad exclusions before subtraction\n"
                         "      --feature-query-lookahead N  feature-track query cache lookahead\n"
                         "  -O, --output FILE                annotated intervals TSV\n"
                         "      --mappability-track FILE     optional merged BED(.gz) track\n"
                         "      --segmental-duplication-track FILE\n"
                         "      --threads N                  Kokkos execution threads\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                has_option(argument, "--region") ? "--region" :
                has_option(argument, "--interval") ? "--interval" : "--intervals";
            options.regions.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals")) {
            options.exclusions.push_back(require_value(index, argc, argv, argument,
                                                        "--exclude-intervals", "-XL"));
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--mappability-track")) {
            options.mappability_track = require_value(index, argc, argv, argument, "--mappability-track");
        } else if (has_option(argument, "--segmental-duplication-track")) {
            options.segmental_duplication_track = require_value(index, argc, argv, argument,
                                                                "--segmental-duplication-track");
        } else if (has_option(argument, "--feature-query-lookahead")) {
            options.feature_query_lookahead = std::stoll(require_value(
                index, argc, argv, argument, "--feature-query-lookahead"));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (has_option(argument, "--interval-merging-rule")) {
            auto rule = require_value(index, argc, argv, argument, "--interval-merging-rule");
            std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char c) {
                return static_cast<char>(std::toupper(c));
            });
            if (rule != "ALL" && rule != "OVERLAPPING_ONLY")
                throw std::invalid_argument("invalid --interval-merging-rule: " + rule);
            options.interval_merging_rule = std::move(rule);
        } else if (has_option(argument, "--interval-set-rule")) {
            auto rule = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char c) {
                return static_cast<char>(std::toupper(c));
            });
            if (rule != "UNION" && rule != "INTERSECTION")
                throw std::invalid_argument("invalid --interval-set-rule: " + rule);
            options.interval_set_rule = std::move(rule);
        } else if (has_option(argument, "--interval-padding") || has_option(argument, "--interval-exclusion-padding")) {
            const bool exclusion = argument.rfind("--interval-exclusion-padding", 0) == 0;
            const char* name = exclusion ? "--interval-exclusion-padding" : "--interval-padding";
            const auto value = std::stoll(require_value(index, argc, argv, argument, name));
            if (value < 0) throw std::invalid_argument(std::string(name) + " must be non-negative");
            // CopyNumberArgumentValidationUtils rejects padding for
            // AnnotateIntervals.  Accept the common switches so the native
            // failure is visible at parse time, but never silently create a
            // result with semantics Java refuses.
            if (value != 0) {
                throw std::invalid_argument(exclusion
                    ? "Interval exclusion padding must be set to 0."
                    : "Interval padding must be set to 0.");
            }
            if (exclusion) options.exclusion_padding = value;
            else options.interval_padding = value;
        } else if (argument == "--quiet" || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates") ||
                   has_option(argument, "--java-options")) {
            if (argument != "--quiet") {
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.regions.empty()) throw std::invalid_argument("-L/--intervals is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    if (options.interval_merging_rule != "OVERLAPPING_ONLY")
        throw std::invalid_argument("Interval merging rule must be set to OVERLAPPING_ONLY.");
    if (options.interval_set_rule != "UNION")
        throw std::invalid_argument("native AnnotateIntervals requires --interval-set-rule UNION");
    // CopyNumberArgumentValidationUtils rejects padding for AnnotateIntervals
    // (the generic interval engine's non-zero defaults are not valid here).
    if (options.interval_padding != 0 || options.exclusion_padding != 0)
        throw std::invalid_argument("Interval padding must be set to 0.");
    return options;
}

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) fields.push_back(std::move(field));
    return fields;
}

bool suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    for (std::size_t i = 0; i < ending.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(value[value.size() - ending.size() + i])) !=
            std::tolower(static_cast<unsigned char>(ending[i]))) return false;
    }
    return true;
}

std::vector<std::string> read_lines(const std::string& path) {
    std::vector<std::string> lines;
    if (suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        kstring_t line{0, 0, nullptr};
        while (fastgatk::io::read_text_line(file, &line, path) >= 0)
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

int reference_tid(const faidx_t* fai, const std::string& contig) {
    for (int index = 0; index < faidx_nseq(fai); ++index)
        if (contig == faidx_iseq(fai, index)) return index;
    return -1;
}

std::vector<TrackSegment> read_bed_track(const std::string& path, const faidx_t* fai) {
    const auto lines = read_lines(path);
    std::vector<TrackSegment> segments;
    for (std::size_t line_number = 0; line_number < lines.size(); ++line_number) {
        const auto first = lines[line_number].find_first_not_of(" \t\r\n");
        if (first == std::string::npos || lines[line_number][first] == '#' ||
            lines[line_number].rfind("track", first) == first ||
            lines[line_number].rfind("browser", first) == first)
            continue;
        const auto fields = split_fields(lines[line_number].substr(first));
        if (fields.size() < 3)
            throw std::invalid_argument("BAD_INPUT: BED track requires contig/start/end at line " +
                                        std::to_string(line_number + 1) + ": " + path);
        const auto tid = reference_tid(fai, fields[0]);
        if (tid < 0) throw std::invalid_argument("BAD_INPUT: BED track contig absent from reference: " + fields[0]);
        std::size_t start_consumed = 0, end_consumed = 0;
        const auto start = std::stoll(fields[1], &start_consumed);
        const auto end = std::stoll(fields[2], &end_consumed);
        const auto length = static_cast<std::int64_t>(faidx_seq_len64(fai, fields[0].c_str()));
        if (start_consumed != fields[1].size() || end_consumed != fields[2].size() ||
            start < 0 || end <= start || end > length)
            throw std::invalid_argument("BAD_INPUT: invalid BED track coordinates at line " +
                                        std::to_string(line_number + 1) + ": " + path);
        // htsjdk's default BEDCodec uses StartOffset.ONE: getStart() is
        // raw_start + 1 and the annotator deliberately queries
        // getEnd() - 1.  In zero-based half-open storage this is therefore
        // [raw_start, raw_end - 1), which excludes the BED end base just as
        // the Java implementation does.  BED column 4 is a name; column 5
        // is the optional score.
        const auto effective_end = end - 1;
        double score = 1.0;
        if (fields.size() >= 5) {
            std::size_t score_consumed = 0;
            score = std::stod(fields[4], &score_consumed);
            if (score_consumed != fields[4].size() || !std::isfinite(score)) score = 1.0;
        }
        segments.push_back(TrackSegment{tid, start, effective_end, score});
    }
    std::sort(segments.begin(), segments.end(), [](const TrackSegment& lhs, const TrackSegment& rhs) {
        if (lhs.tid != rhs.tid) return lhs.tid < rhs.tid;
        if (lhs.start != rhs.start) return lhs.start < rhs.start;
        return lhs.end < rhs.end;
    });
    for (std::size_t index = 1; index < segments.size(); ++index) {
        if (segments[index].tid == segments[index - 1].tid &&
            segments[index].start < segments[index - 1].end)
            throw std::invalid_argument("BAD_INPUT: BED track intervals must be merged before AnnotateIntervals: " + path);
    }
    return segments;
}

Interval parse_token(const std::string& token, const faidx_t* fai) {
    const auto colon = token.find(':');
    const std::string contig = colon == std::string::npos ? token : token.substr(0, colon);
    if (contig.empty() || faidx_seq_len64(fai, contig.c_str()) < 0)
        throw std::invalid_argument("BAD_INPUT: interval contig is absent from reference: " + contig);
    const std::int64_t length = static_cast<std::int64_t>(faidx_seq_len64(fai, contig.c_str()));
    if (colon == std::string::npos) return {contig, 1, length};
    const auto dash = token.find('-', colon + 1);
    const auto start = std::stoll(token.substr(colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1));
    const auto end = dash == std::string::npos ? start : std::stoll(token.substr(dash + 1));
    if (start < 1 || end < start || end > length)
        throw std::invalid_argument("BAD_INPUT: interval is outside reference: " + token);
    return {contig, start, end};
}

std::vector<Interval> parse_selector(const std::string& selector, const faidx_t* fai) {
    std::error_code error;
    const bool file = std::filesystem::is_regular_file(selector, error) && !error;
    if (!file) return {parse_token(selector, fai)};
    const bool bed = suffix_ci(selector, ".bed") || suffix_ci(selector, ".bed.gz");
    std::vector<Interval> result;
    const auto lines = read_lines(selector);
    for (std::size_t line_number = 0; line_number < lines.size(); ++line_number) {
        const auto first = lines[line_number].find_first_not_of(" \t\r\n");
        if (first == std::string::npos || lines[line_number][first] == '#' || lines[line_number][first] == '@') continue;
        const auto fields = split_fields(lines[line_number].substr(first));
        if (fields.empty() || fields[0] == "track" || fields[0] == "browser") continue;
        try {
            if (bed) {
                if (fields.size() < 3) throw std::invalid_argument("BED requires contig/start/end");
                const auto start0 = std::stoll(fields[1]);
                const auto end0 = std::stoll(fields[2]);
                if (start0 < 0 || end0 <= start0) throw std::invalid_argument("invalid BED coordinates");
                result.push_back(parse_token(fields[0] + ":" + std::to_string(start0 + 1) + "-" + std::to_string(end0), fai));
            } else if (fields.size() == 1) {
                result.push_back(parse_token(fields[0], fai));
            } else if (fields.size() >= 3) {
                result.push_back(parse_token(fields[0] + ":" + fields[1] + "-" + fields[2], fai));
            } else {
                throw std::invalid_argument("interval list row requires contig/start/end");
            }
        } catch (const std::exception& exception) {
            throw std::invalid_argument("BAD_INPUT: malformed interval line " + std::to_string(line_number + 1) +
                                        " in " + selector + ": " + exception.what());
        }
    }
    if (result.empty()) throw std::invalid_argument("BAD_INPUT: interval file is empty: " + selector);
    return result;
}

std::vector<Interval> normalize_intervals(const std::vector<Interval>& input, const faidx_t* fai,
                                          const std::string& merging_rule = "ALL") {
    std::unordered_map<std::string, int> order;
    for (int i = 0; i < faidx_nseq(fai); ++i) order.emplace(faidx_iseq(fai, i), i);
    auto result = input;
    std::sort(result.begin(), result.end(), [&](const Interval& lhs, const Interval& rhs) {
        const auto li = order.at(lhs.contig), ri = order.at(rhs.contig);
        return li != ri ? li < ri : (lhs.start != rhs.start ? lhs.start < rhs.start : lhs.end < rhs.end);
    });
    std::vector<Interval> merged;
    for (const auto& interval : result) {
        const bool merge_abutting = merging_rule == "ALL";
        const bool overlaps = !merged.empty() && merged.back().contig == interval.contig &&
            (merge_abutting ? interval.start <= merged.back().end + 1
                             : interval.start <= merged.back().end);
        if (overlaps) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }
    return merged;
}

std::vector<Interval> pad_intervals(const std::vector<Interval>& input, std::int64_t padding,
                                    const faidx_t* fai) {
    if (padding == 0) return input;
    std::vector<Interval> result;
    result.reserve(input.size());
    for (const auto& interval : input) {
        const auto length = static_cast<std::int64_t>(faidx_seq_len64(fai, interval.contig.c_str()));
        result.push_back({interval.contig,
                          std::max<std::int64_t>(1, interval.start - padding),
                          std::min<std::int64_t>(length, interval.end + padding)});
    }
    return result;
}

std::vector<Interval> intersect_intervals(const std::vector<Interval>& lhs,
                                          const std::vector<Interval>& rhs,
                                          const faidx_t* fai) {
    const auto left = normalize_intervals(lhs, fai, "ALL");
    const auto right = normalize_intervals(rhs, fai, "ALL");
    std::vector<Interval> result;
    std::size_t li = 0, ri = 0;
    while (li < left.size() && ri < right.size()) {
        if (left[li].contig < right[ri].contig) { ++li; continue; }
        if (right[ri].contig < left[li].contig) { ++ri; continue; }
        const auto start = std::max(left[li].start, right[ri].start);
        const auto end = std::min(left[li].end, right[ri].end);
        if (start <= end) result.push_back({left[li].contig, start, end});
        if (left[li].end < right[ri].end) ++li;
        else ++ri;
    }
    return result;
}

std::vector<Interval> subtract_intervals(const std::vector<Interval>& input,
                                         const std::vector<Interval>& exclusions,
                                         const faidx_t* fai) {
    if (exclusions.empty()) return input;
    const auto normalized_exclusions = normalize_intervals(exclusions, fai, "ALL");
    std::vector<Interval> result;
    for (const auto& source : input) {
        std::int64_t cursor = source.start;
        for (const auto& exclusion : normalized_exclusions) {
            if (exclusion.contig != source.contig) continue;
            if (exclusion.start > source.end) break;
            if (exclusion.end < cursor) continue;
            if (exclusion.start > cursor)
                result.push_back({source.contig, cursor,
                                  std::min(source.end, exclusion.start - 1)});
            if (exclusion.end >= source.end) { cursor = source.end + 1; break; }
            cursor = std::max(cursor, exclusion.end + 1);
        }
        if (cursor <= source.end) result.push_back({source.contig, cursor, source.end});
    }
    return result;
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
    faidx_t* fai = fai_load(options.reference.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
    try {
        std::vector<Interval> requested;
        bool first_selector = true;
        for (const auto& selector : options.regions) {
            auto selected = pad_intervals(parse_selector(selector, fai), options.interval_padding, fai);
            selected = normalize_intervals(selected, fai, options.interval_merging_rule);
            if (options.interval_set_rule == "INTERSECTION" && !first_selector)
                requested = intersect_intervals(requested, selected, fai);
            else {
                requested.insert(requested.end(), selected.begin(), selected.end());
                if (options.interval_set_rule == "UNION")
                    requested = normalize_intervals(requested, fai, options.interval_merging_rule);
            }
            first_selector = false;
        }
        std::vector<Interval> requested_exclusions;
        for (const auto& selector : options.exclusions) {
            const auto selected = parse_selector(selector, fai);
            requested_exclusions.insert(requested_exclusions.end(), selected.begin(), selected.end());
        }
        const auto padded_exclusions = pad_intervals(requested_exclusions, options.exclusion_padding, fai);
        const auto intervals = normalize_intervals(
            subtract_intervals(requested, padded_exclusions, fai), fai, options.interval_merging_rule);
        if (intervals.empty()) throw std::invalid_argument("BAD_INPUT: no intervals remain after normalization");
        const auto mappability = options.mappability_track.empty()
            ? std::vector<TrackSegment>{} : read_bed_track(options.mappability_track, fai);
        const auto segmental_duplication = options.segmental_duplication_track.empty()
            ? std::vector<TrackSegment>{} : read_bed_track(options.segmental_duplication_track, fai);

        std::vector<std::string> sequences;
        sequences.reserve(intervals.size());
        std::vector<std::size_t> offsets{0};
        for (const auto& interval : intervals) {
            hts_pos_t length = 0;
            char* sequence = faidx_fetch_seq64(fai, interval.contig.c_str(), interval.start - 1, interval.end - 1, &length);
            if (!sequence || length != interval.end - interval.start + 1) {
                free(sequence);
                throw std::runtime_error("BAD_INPUT: failed to fetch reference interval: " + interval.contig);
            }
            sequences.emplace_back(sequence, static_cast<std::size_t>(length));
            offsets.push_back(offsets.back() + static_cast<std::size_t>(length));
            free(sequence);
        }
        std::vector<unsigned char> packed;
        packed.reserve(offsets.back());
        for (const auto& sequence : sequences)
            packed.insert(packed.end(), sequence.begin(), sequence.end());

        std::vector<int> interval_tids;
        interval_tids.reserve(intervals.size());
        for (const auto& interval : intervals) interval_tids.push_back(reference_tid(fai, interval.contig));

        // Tracks are sorted by (reference tid, start).  Keep a compact per-
        // contig index so every interval can binary-search the first possible
        // overlap instead of scanning all BED rows.  This is important for a
        // WGS interval list where the track can contain millions of segments.
        const auto contig_count = static_cast<std::size_t>(faidx_nseq(fai));
        auto build_track_offsets = [contig_count](const std::vector<TrackSegment>& track) {
            std::vector<std::size_t> offsets(contig_count + 1, 0);
            for (const auto& segment : track) {
                if (segment.tid >= 0 && static_cast<std::size_t>(segment.tid) < contig_count)
                    ++offsets[static_cast<std::size_t>(segment.tid) + 1];
            }
            for (std::size_t index = 1; index < offsets.size(); ++index)
                offsets[index] += offsets[index - 1];
            return offsets;
        };
        const auto mappability_offsets = build_track_offsets(mappability);
        const auto segmental_offsets = build_track_offsets(segmental_duplication);

        Kokkos::InitializationSettings settings;
        settings.set_num_threads(options.threads);
        Kokkos::initialize(settings);
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch annotation_host("annotate-intervals-gc-v1");
        annotation_host.records = intervals.size();
        annotation_host.bytes = packed.size() * sizeof(unsigned char) +
                                offsets.size() * sizeof(std::size_t) +
                                interval_tids.size() * sizeof(int) +
                                intervals.size() * 2 * sizeof(std::int64_t) +
                                (mappability_offsets.size() + segmental_offsets.size()) * sizeof(std::size_t) +
                                (mappability.size() + segmental_duplication.size()) * sizeof(TrackSegment);
        fastgatk::core::KernelPlan<ExecSpace> annotation_plan("annotate-intervals-gc-tracks");
        annotation_plan.begin_prepare(annotation_host);
        std::vector<std::uint64_t> gc_counts(intervals.size());
        std::vector<std::uint64_t> at_counts(intervals.size());
        std::vector<double> mappability_values(intervals.size(), 0.0);
        std::vector<double> segmental_duplication_values(intervals.size(), 0.0);
        {
            Kokkos::View<unsigned char*> device_bases("annotated_interval_bases", packed.size());
            Kokkos::View<std::size_t*> device_offsets("annotated_interval_offsets", offsets.size());
            Kokkos::View<int*> device_tids("annotated_interval_tids", intervals.size());
            Kokkos::View<std::int64_t*> device_starts("annotated_interval_starts", intervals.size());
            Kokkos::View<std::int64_t*> device_ends("annotated_interval_ends", intervals.size());
            Kokkos::View<std::uint64_t*> device_gc("annotated_interval_gc", intervals.size());
            Kokkos::View<std::uint64_t*> device_at("annotated_interval_at", intervals.size());
            Kokkos::View<int*> device_mappability_tids("mappability_tids", mappability.size());
            Kokkos::View<std::int64_t*> device_mappability_starts("mappability_starts", mappability.size());
            Kokkos::View<std::int64_t*> device_mappability_ends("mappability_ends", mappability.size());
            Kokkos::View<double*> device_mappability_scores("mappability_scores", mappability.size());
            Kokkos::View<std::size_t*> device_mappability_offsets("mappability_offsets", mappability_offsets.size());
            Kokkos::View<double*> device_mappability("mappability", intervals.size());
            Kokkos::View<int*> device_segmental_tids("segmental_tids", segmental_duplication.size());
            Kokkos::View<std::int64_t*> device_segmental_starts("segmental_starts", segmental_duplication.size());
            Kokkos::View<std::int64_t*> device_segmental_ends("segmental_ends", segmental_duplication.size());
            Kokkos::View<double*> device_segmental_scores("segmental_scores", segmental_duplication.size());
            Kokkos::View<std::size_t*> device_segmental_offsets("segmental_offsets", segmental_offsets.size());
            Kokkos::View<double*> device_segmental("segmental_duplication", intervals.size());
            auto host_bases = Kokkos::create_mirror_view(device_bases);
            auto host_offsets = Kokkos::create_mirror_view(device_offsets);
            auto host_tids = Kokkos::create_mirror_view(device_tids);
            auto host_starts = Kokkos::create_mirror_view(device_starts);
            auto host_ends = Kokkos::create_mirror_view(device_ends);
            auto host_mappability_tids = Kokkos::create_mirror_view(device_mappability_tids);
            auto host_mappability_starts = Kokkos::create_mirror_view(device_mappability_starts);
            auto host_mappability_ends = Kokkos::create_mirror_view(device_mappability_ends);
            auto host_mappability_scores = Kokkos::create_mirror_view(device_mappability_scores);
            auto host_mappability_offsets = Kokkos::create_mirror_view(device_mappability_offsets);
            auto host_segmental_tids = Kokkos::create_mirror_view(device_segmental_tids);
            auto host_segmental_starts = Kokkos::create_mirror_view(device_segmental_starts);
            auto host_segmental_ends = Kokkos::create_mirror_view(device_segmental_ends);
            auto host_segmental_scores = Kokkos::create_mirror_view(device_segmental_scores);
            auto host_segmental_offsets = Kokkos::create_mirror_view(device_segmental_offsets);
            for (std::size_t i = 0; i < packed.size(); ++i) host_bases(i) = packed[i];
            for (std::size_t i = 0; i < offsets.size(); ++i) host_offsets(i) = offsets[i];
            for (std::size_t i = 0; i < intervals.size(); ++i) {
                host_tids(i) = interval_tids[i];
                host_starts(i) = intervals[i].start - 1;
                host_ends(i) = intervals[i].end;
            }
            for (std::size_t i = 0; i < mappability.size(); ++i) {
                host_mappability_tids(i) = mappability[i].tid;
                host_mappability_starts(i) = mappability[i].start;
                host_mappability_ends(i) = mappability[i].end;
                host_mappability_scores(i) = mappability[i].score;
            }
            for (std::size_t i = 0; i < mappability_offsets.size(); ++i)
                host_mappability_offsets(i) = mappability_offsets[i];
            for (std::size_t i = 0; i < segmental_duplication.size(); ++i) {
                host_segmental_tids(i) = segmental_duplication[i].tid;
                host_segmental_starts(i) = segmental_duplication[i].start;
                host_segmental_ends(i) = segmental_duplication[i].end;
                host_segmental_scores(i) = segmental_duplication[i].score;
            }
            for (std::size_t i = 0; i < segmental_offsets.size(); ++i)
                host_segmental_offsets(i) = segmental_offsets[i];
            Kokkos::deep_copy(device_bases, host_bases);
            Kokkos::deep_copy(device_offsets, host_offsets);
            Kokkos::deep_copy(device_tids, host_tids);
            Kokkos::deep_copy(device_starts, host_starts);
            Kokkos::deep_copy(device_ends, host_ends);
            Kokkos::deep_copy(device_mappability_tids, host_mappability_tids);
            Kokkos::deep_copy(device_mappability_starts, host_mappability_starts);
            Kokkos::deep_copy(device_mappability_ends, host_mappability_ends);
            Kokkos::deep_copy(device_mappability_scores, host_mappability_scores);
            Kokkos::deep_copy(device_mappability_offsets, host_mappability_offsets);
            Kokkos::deep_copy(device_segmental_tids, host_segmental_tids);
            Kokkos::deep_copy(device_segmental_starts, host_segmental_starts);
            Kokkos::deep_copy(device_segmental_ends, host_segmental_ends);
            Kokkos::deep_copy(device_segmental_scores, host_segmental_scores);
            Kokkos::deep_copy(device_segmental_offsets, host_segmental_offsets);
            Kokkos::deep_copy(device_gc, std::uint64_t{0});
            Kokkos::deep_copy(device_at, std::uint64_t{0});
            fastgatk::core::DeviceBatch<ExecSpace> annotation_device(intervals.size());
            annotation_device.bind("bases", device_bases);
            annotation_device.bind("offsets", device_offsets);
            annotation_device.bind("tids", device_tids);
            annotation_device.bind("starts", device_starts);
            annotation_device.bind("ends", device_ends);
            annotation_device.bind("gc", device_gc);
            annotation_device.bind("at", device_at);
            annotation_device.bind("mappability_tids", device_mappability_tids);
            annotation_device.bind("mappability_starts", device_mappability_starts);
            annotation_device.bind("mappability_ends", device_mappability_ends);
            annotation_device.bind("mappability_scores", device_mappability_scores);
            annotation_device.bind("mappability_offsets", device_mappability_offsets);
            annotation_device.bind("mappability", device_mappability);
            annotation_device.bind("segmental_tids", device_segmental_tids);
            annotation_device.bind("segmental_starts", device_segmental_starts);
            annotation_device.bind("segmental_ends", device_segmental_ends);
            annotation_device.bind("segmental_scores", device_segmental_scores);
            annotation_device.bind("segmental_offsets", device_segmental_offsets);
            annotation_device.bind("segmental_duplication", device_segmental);
            ExecSpace().fence();
            annotation_plan.end_prepare(annotation_device);
            annotation_plan.begin_execute();
            Kokkos::parallel_for("annotate_intervals_gc_tracks", Kokkos::RangePolicy<ExecSpace>(0, intervals.size()),
                KOKKOS_LAMBDA(const std::size_t interval) {
                    std::uint64_t gc = 0, at = 0;
                    for (std::size_t position = device_offsets(interval); position < device_offsets(interval + 1); ++position) {
                        const auto base = device_bases(position);
                        if (base == 'G' || base == 'g' || base == 'C' || base == 'c') ++gc;
                        else if (base == 'A' || base == 'a' || base == 'T' || base == 't') ++at;
                    }
                    device_gc(interval) = gc;
                    device_at(interval) = at;
                    const auto start = device_starts(interval);
                    const auto end = device_ends(interval);
                    const auto length = end - start;
                    double map_sum = 0.0;
                    const auto tid = device_tids(interval);
                    const auto map_contig_begin = device_mappability_offsets(tid);
                    const auto map_contig_end = device_mappability_offsets(tid + 1);
                    std::size_t map_begin = map_contig_begin;
                    std::size_t map_end = map_contig_end;
                    while (map_begin < map_end) {
                        const auto mid = map_begin + (map_end - map_begin) / 2;
                        if (device_mappability_ends(mid) <= start) map_begin = mid + 1;
                        else map_end = mid;
                    }
                    for (std::size_t track = map_begin; track < map_contig_end &&
                         device_mappability_starts(track) < end; ++track) {
                        const auto overlap_start = start > device_mappability_starts(track)
                            ? start : device_mappability_starts(track);
                        const auto overlap_end = end < device_mappability_ends(track)
                            ? end : device_mappability_ends(track);
                        if (overlap_end > overlap_start)
                            map_sum += static_cast<double>(overlap_end - overlap_start) * device_mappability_scores(track);
                    }
                    device_mappability(interval) = length > 0 ? map_sum / static_cast<double>(length) : 0.0;
                    double segmental_sum = 0.0;
                    const auto seg_contig_begin = device_segmental_offsets(tid);
                    const auto seg_contig_end = device_segmental_offsets(tid + 1);
                    std::size_t seg_begin = seg_contig_begin;
                    std::size_t seg_end = seg_contig_end;
                    while (seg_begin < seg_end) {
                        const auto mid = seg_begin + (seg_end - seg_begin) / 2;
                        if (device_segmental_ends(mid) <= start) seg_begin = mid + 1;
                        else seg_end = mid;
                    }
                    for (std::size_t track = seg_begin; track < seg_contig_end &&
                         device_segmental_starts(track) < end; ++track) {
                        const auto overlap_start = start > device_segmental_starts(track)
                            ? start : device_segmental_starts(track);
                        const auto overlap_end = end < device_segmental_ends(track)
                            ? end : device_segmental_ends(track);
                        if (overlap_end > overlap_start)
                            segmental_sum += static_cast<double>(overlap_end - overlap_start) * device_segmental_scores(track);
                    }
                    device_segmental(interval) = length > 0 ? segmental_sum / static_cast<double>(length) : 0.0;
                });
            ExecSpace().fence();
            annotation_plan.end_execute();
            auto host_gc = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_gc);
            auto host_at = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_at);
            auto host_mappability = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_mappability);
            auto host_segmental = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_segmental);
            for (std::size_t i = 0; i < intervals.size(); ++i) {
                gc_counts[i] = host_gc(i);
                at_counts[i] = host_at(i);
                mappability_values[i] = host_mappability(i);
                segmental_duplication_values[i] = host_segmental(i);
            }
        }
        Kokkos::finalize();

        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write annotated intervals: " + options.output);
        output << "@HD\tVN:1.6\n";
        for (int i = 0; i < faidx_nseq(fai); ++i)
            output << "@SQ\tSN:" << faidx_iseq(fai, i) << "\tLN:" << faidx_seq_len64(fai, faidx_iseq(fai, i)) << '\n';
        output << "CONTIG\tSTART\tEND\tGC_CONTENT";
        if (!options.mappability_track.empty()) output << "\tMAPPABILITY";
        if (!options.segmental_duplication_track.empty()) output << "\tSEGMENTAL_DUPLICATION_CONTENT";
        output << '\n' << std::fixed << std::setprecision(6);
        for (std::size_t i = 0; i < intervals.size(); ++i) {
            const auto denominator = gc_counts[i] + at_counts[i];
            output << intervals[i].contig << '\t' << intervals[i].start << '\t' << intervals[i].end << '\t';
            if (denominator == 0) output << "NaN";
            else output << static_cast<double>(gc_counts[i]) / static_cast<double>(denominator);
            if (!options.mappability_track.empty()) {
                output << '\t';
                if (std::isfinite(mappability_values[i])) output << mappability_values[i];
                else output << "NaN";
            }
            if (!options.segmental_duplication_track.empty()) {
                output << '\t';
                if (std::isfinite(segmental_duplication_values[i])) output << segmental_duplication_values[i];
                else output << "NaN";
            }
            output << '\n';
        }
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize annotated intervals: " + options.output);
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"AnnotateIntervals\",\"implementation\":\"fastgatk-annotate-intervals\",\"status\":\"prototype\",\"reference\":\""
                     << json_escape(options.reference) << "\",\"intervals\":" << intervals.size()
                     << ",\"excluded_intervals\":" << requested_exclusions.size()
                     << ",\"interval_merging_rule\":\"" << options.interval_merging_rule
                     << "\",\"interval_set_rule\":\"" << options.interval_set_rule
                     << "\",\"interval_padding\":" << options.interval_padding
                     << ",\"interval_exclusion_padding\":" << options.exclusion_padding
                     << ",\"mappability_segments\":" << mappability.size()
                     << ",\"segmental_duplication_segments\":" << segmental_duplication.size()
                     << ",\"feature_query_lookahead\":" << options.feature_query_lookahead
                     << ",\"mappability_track\":" << (options.mappability_track.empty() ? "null" :
                         ("\"" + json_escape(options.mappability_track) + "\""))
                     << ",\"segmental_duplication_track\":" << (options.segmental_duplication_track.empty() ? "null" :
                         ("\"" + json_escape(options.segmental_duplication_track) + "\""))
                     << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\",\"output\":\""
                     << json_escape(options.output) << "\",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\"" << ExecSpace::name() << "\",\"kernel_execution_policy\":\"RangePolicy\",\"track_index_strategy\":\"per-contig-binary-search\",\"kernel_batches\":1,\"kernel_records\":" << intervals.size()
                     << ",\"kernel_prepare_seconds\":" << annotation_plan.telemetry().prepare_seconds
                     << ",\"feature_query_lookahead\":" << options.feature_query_lookahead
                     << ",\"kernel_execute_seconds\":" << annotation_plan.telemetry().execute_seconds << "}}\n";
        }
        std::cout << "{\"tool\":\"AnnotateIntervals\",\"status\":\"prototype\",\"intervals\":"
                  << intervals.size() << ",\"kernel_execution_space\":\"" << ExecSpace::name()
                  << "\",\"excluded_intervals\":" << requested_exclusions.size()
                  << ",\"interval_merging_rule\":\"" << options.interval_merging_rule
                  << "\",\"interval_set_rule\":\"" << options.interval_set_rule
                  << "\",\"kernel_records\":" << intervals.size()
                  << ",\"mappability\":" << (!options.mappability_track.empty() ? "true" : "false")
                  << ",\"segmental_duplication\":" << (!options.segmental_duplication_track.empty() ? "true" : "false")
                  << ",\"mappability_segments\":" << mappability.size()
                  << ",\"segmental_duplication_segments\":" << segmental_duplication.size()
                  << ",\"feature_query_lookahead\":" << options.feature_query_lookahead
                  << ",\"kernel_execute_seconds\":" << annotation_plan.telemetry().execute_seconds
                  << ",\"output\":\"" << json_escape(options.output) << "\"}\n";
        fai_destroy(fai);
        return 0;
    } catch (...) {
        fai_destroy(fai);
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
