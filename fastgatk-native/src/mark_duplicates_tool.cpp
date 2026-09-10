#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/output.hpp"
#include "optional_boolean.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#if FASTGATK_HAS_HTSLIB
#include <htslib/sam.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string reference;
    std::string metrics;
    std::string manifest;
    std::string tmp_dir;
    std::string read_name_regex;
    std::string tagging_policy = "DontTag";
    std::size_t max_records_in_memory = 0;
    std::int64_t optical_duplicate_pixel_distance = 100;
    bool clear_duplicates = false;
    bool remove_duplicates = false;
    bool remove_sequencing_duplicates = false;
    // Keep the completed first-pass metadata runs in tmp_dir so an interrupted
    // second pass can be resumed without rereading/rehashing the input.
    bool resume_spill = false;
    // Picard/GATK MarkDuplicates 4.6.2.0 defaults CREATE_INDEX=false.  Keep
    // index creation opt-in so a direct replacement does not add an
    // unexpected sidecar to every coordinate-sorted BAM output.
    bool create_index = false;
    bool assume_sorted = false;
    // Picard MarkDuplicates defaults to ADD_PG_TAG_TO_READS=true.  Keep the
    // native default aligned with that behavior; callers that need a header
    // without a program record can pass --add-pg-tag=false.
    bool add_pg_tag = true;
    // Picard's default PROGRAM_RECORD_ID is MarkDuplicates; HTSlib will
    // preserve the caller-supplied override for chained @PG provenance.
    std::string pg_id = "MarkDuplicates";
    std::string pg_name = "MarkDuplicates";
    std::string pg_version = "fastgatk-native";
    std::string pg_command_line;
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
            std::cout << "fastgatk-mark-duplicates (GATK/Picard-compatible native prototype)\n"
                         "  -I, --input FILE                 coordinate-sorted SAM/BAM/CRAM\n"
                         "  -O, --output FILE                output SAM/BAM/CRAM\n"
                         "  --metrics-file FILE              duplicate metrics\n"
                         "  --remove-duplicates              omit all duplicate records\n"
                         "  --remove-sequencing-duplicates   omit only optical duplicate records\n"
                         "  --max-records-in-memory N       fail closed above this group bound\n"
                         "  --tmp-dir DIR                   scratch directory contract\n"
                         "  --resume-spill                  resume a completed input-pass checkpoint\n"
                         "  --optical-duplicate-pixel-distance N\n"
                         "  --tagging-policy DontTag|OpticalOnly|All   optional DT duplicate tag\n"
                         "  --add-pg-tag[=true|false]       add a SAM @PG record (default true)\n"
                         "  --program-record-id ID          @PG ID\n"
                         "  --program-group-name NAME       @PG PN\n"
                         "  --program-group-version VERSION @PG VN\n"
                         "  --program-group-command-line CL @PG CL\n"
                         "  -R, --reference FILE             CRAM reference (optional)\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-I" || is_option(argument, "--input"))
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (is_option(argument, "--metrics-file") || is_option(argument, "--output-metrics"))
            options.metrics = require_value(index, argc, argv, argument,
                argument.rfind("--output-metrics", 0) == 0 ? "--output-metrics" : "--metrics-file");
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (is_option(argument, "--tmp-dir") || is_option(argument, "--TMP_DIR"))
            options.tmp_dir = require_value(index, argc, argv, argument,
                argument.rfind("--TMP_DIR", 0) == 0 ? "--TMP_DIR" : "--tmp-dir");
        else if (is_option(argument, "--max-records-in-memory") || is_option(argument, "--MAX_RECORDS_IN_RAM")) {
            const auto value = require_value(index, argc, argv, argument,
                argument.rfind("--MAX_RECORDS_IN_RAM", 0) == 0 ? "--MAX_RECORDS_IN_RAM" : "--max-records-in-memory");
            try {
                options.max_records_in_memory = static_cast<std::size_t>(std::stoull(value));
            } catch (const std::exception&) {
                throw std::invalid_argument("invalid max-records-in-memory: " + value);
            }
            if (options.max_records_in_memory == 0)
                throw std::invalid_argument("max-records-in-memory must be positive");
        }
        else if (is_option(argument, "--read-name-regex"))
            options.read_name_regex = require_value(index, argc, argv, argument, "--read-name-regex");
        else if (is_option(argument, "--tagging-policy"))
            options.tagging_policy = require_value(index, argc, argv, argument, "--tagging-policy");
        else if (argument == "--add-pg-tag" || argument.rfind("--add-pg-tag=", 0) == 0) {
            options.add_pg_tag = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--add-pg-tag");
        }
        else if (is_option(argument, "--program-record-id"))
            options.pg_id = require_value(index, argc, argv, argument, "--program-record-id");
        else if (is_option(argument, "--program-group-name"))
            options.pg_name = require_value(index, argc, argv, argument, "--program-group-name");
        else if (is_option(argument, "--program-group-version"))
            options.pg_version = require_value(index, argc, argv, argument, "--program-group-version");
        else if (is_option(argument, "--program-group-command-line"))
            options.pg_command_line = require_value(index, argc, argv, argument, "--program-group-command-line");
        else if (is_option(argument, "--optical-duplicate-pixel-distance")) {
            const auto value = require_value(index, argc, argv, argument, "--optical-duplicate-pixel-distance");
            try {
                options.optical_duplicate_pixel_distance = std::stoll(value);
            } catch (const std::exception&) {
                throw std::invalid_argument("invalid optical-duplicate-pixel-distance: " + value);
            }
            if (options.optical_duplicate_pixel_distance < 0)
                throw std::invalid_argument("optical-duplicate-pixel-distance must be non-negative");
        }
        else if (argument == "--remove-duplicates")
            options.remove_duplicates = true;
        else if (argument == "--remove-sequencing-duplicates")
            options.remove_sequencing_duplicates = true;
        else if (argument == "--resume-spill")
            options.resume_spill = true;
        else if (argument == "--assume-sorted" || argument == "--assume-sort-order")
            options.assume_sorted = true;
        else if (argument == "--create-output-bam-index" ||
                 argument.rfind("--create-output-bam-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-bam-index");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--use-jdk-deflater" || argument == "--use-jdk-inflater") {
        } else if (argument == "--clear-duplicates") {
            options.clear_duplicates = true;
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.resume_spill && options.tmp_dir.empty())
        throw std::invalid_argument("--resume-spill requires --tmp-dir");
    if (options.tagging_policy != "DontTag" && options.tagging_policy != "OpticalOnly" &&
        options.tagging_policy != "All" && options.tagging_policy != "donttag" &&
        options.tagging_policy != "opticalonly" && options.tagging_policy != "all")
        throw std::invalid_argument("unsupported --tagging-policy: " + options.tagging_policy);
    if (options.add_pg_tag && options.pg_id.empty())
        throw std::invalid_argument("--program-record-id must not be empty when --add-pg-tag is enabled");
    if (options.metrics.empty()) options.metrics = options.output + ".metrics.txt";
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct DuplicateKey {
    int tid = -1;
    std::int64_t pos = -1;
    int mtid = -1;
    std::int64_t mpos = -1;
    int orientation = 0;
    std::string library;

    bool operator==(const DuplicateKey& other) const {
        return std::tie(tid, pos, mtid, mpos, orientation, library) ==
               std::tie(other.tid, other.pos, other.mtid, other.mpos, other.orientation, other.library);
    }
};

bool key_less(const DuplicateKey& left, const DuplicateKey& right) {
    return std::tie(left.tid, left.pos, left.mtid, left.mpos, left.orientation, left.library) <
           std::tie(right.tid, right.pos, right.mtid, right.mpos, right.orientation, right.library);
}

struct DuplicateKeyLess {
    bool operator()(const DuplicateKey& left, const DuplicateKey& right) const {
        return key_less(left, right);
    }
};

std::int64_t five_prime_position(const bam1_t* record) {
    if (record->core.tid < 0 || record->core.pos < 0) return -1;
    // Picard's ReadEnds uses SAMRecord.getUnclippedStart/End.  HTSlib's
    // bam_endpos() intentionally excludes terminal soft/hard clips, so add
    // the terminal clip on reverse reads and subtract the leading clip on
    // forward reads before forming a duplicate key.
    std::int64_t leading_clip = 0;
    std::int64_t trailing_clip = 0;
    if (record->core.n_cigar > 0) {
        const auto* cigar = bam_get_cigar(record);
        std::size_t index = 0;
        while (index < record->core.n_cigar) {
            const auto operation = bam_cigar_op(cigar[index]);
            if (operation != BAM_CSOFT_CLIP && operation != BAM_CHARD_CLIP) break;
            leading_clip += bam_cigar_oplen(cigar[index++]);
        }
        index = record->core.n_cigar;
        while (index > 0) {
            const auto operation = bam_cigar_op(cigar[index - 1]);
            if (operation != BAM_CSOFT_CLIP && operation != BAM_CHARD_CLIP) break;
            trailing_clip += bam_cigar_oplen(cigar[--index]);
        }
    }
    if ((record->core.flag & BAM_FREVERSE) != 0)
        return bam_endpos(record) - 1 + trailing_clip;
    return record->core.pos - leading_clip;
}

std::int64_t mate_five_prime_position(const bam1_t* record) {
    if (record == nullptr || record->core.mtid < 0 || record->core.mpos < 0) return -1;
    auto position = static_cast<std::int64_t>(record->core.mpos);
    if ((record->core.flag & BAM_FMREVERSE) == 0) return position;
    // The mate's unclipped 5' coordinate is its alignment end when the mate
    // is reverse-strand.  MC is emitted by modern aligners and lets us derive
    // that coordinate without loading a second record; when it is absent we
    // retain the SAM leftmost position as a conservative fallback.
    const auto* tag = bam_aux_get(record, "MC");
    if (tag != nullptr && *tag == 'Z') {
        const auto* text = bam_aux2Z(tag);
        if (text != nullptr) {
            std::int64_t reference_span = 0;
            const char* cursor = text;
            while (*cursor != '\0') {
                char* end = nullptr;
                const auto length = std::strtoll(cursor, &end, 10);
                if (end == cursor || length < 0 || *end == '\0') break;
                const char operation = *end++;
                if (operation == 'M' || operation == 'D' || operation == 'N' ||
                    operation == '=' || operation == 'X') reference_span += length;
                cursor = end;
            }
            if (reference_span > 0) position += reference_span - 1;
        }
    }
    return position;
}

std::string read_group(const bam1_t* record) {
    const auto* tag = bam_aux_get(record, "RG");
    if (tag == nullptr || *tag != 'Z') return {};
    const auto* value = bam_aux2Z(tag);
    return value == nullptr ? std::string{} : std::string(value);
}

std::string library_for_read_group(sam_hdr_t* header, const std::string& rg) {
    if (header == nullptr || rg.empty()) return {};
    kstring_t value = {0, 0, nullptr};
    const int status = sam_hdr_find_tag_id(header, "RG", "ID", rg.c_str(), "LB", &value);
    std::string result;
    if (status == 0 && value.s != nullptr) result = value.s;
    free(value.s);
    return result;
}

DuplicateKey duplicate_key(const bam1_t* record, sam_hdr_t* header) {
    const bool paired = (record->core.flag & BAM_FPAIRED) != 0 && record->core.mtid >= 0;
    auto own_tid = record->core.tid;
    auto own_position = five_prime_position(record);
    auto own_reverse = (record->core.flag & BAM_FREVERSE) != 0;
    auto mate_tid = paired ? record->core.mtid : -1;
    auto mate_position = paired ? mate_five_prime_position(record) : -1;
    auto mate_reverse = paired && ((record->core.flag & BAM_FMREVERSE) != 0);
    // Canonicalize the two ends so read1/read2 form one duplicate fragment.
    // The orientation bits follow the sorted endpoint order and therefore are
    // identical when the same fragment is encountered from either mate.
    if (paired && std::tie(mate_tid, mate_position, mate_reverse) <
                  std::tie(own_tid, own_position, own_reverse)) {
        std::swap(own_tid, mate_tid);
        std::swap(own_position, mate_position);
        std::swap(own_reverse, mate_reverse);
    }
    const int orientation = static_cast<int>(own_reverse ? 1 : 0) |
                            static_cast<int>(mate_reverse ? 2 : 0);
    return {own_tid, own_position, paired ? mate_tid : -1,
            paired ? mate_position : -1, orientation,
            library_for_read_group(header, read_group(record))};
}

struct OpticalCoordinate {
    std::int64_t tile = 0;
    std::int64_t x = 0;
    std::int64_t y = 0;
};

bool parse_decimal(const std::string& value, std::int64_t& result) {
    if (value.empty()) return false;
    try {
        std::size_t consumed = 0;
        result = std::stoll(value, &consumed);
        return consumed == value.size();
    } catch (const std::exception&) {
        return false;
    }
}

bool parse_optical_coordinate(const char* qname, const std::string& expression,
                              OpticalCoordinate& coordinate) {
    if (qname == nullptr) return false;
    const std::string name(qname);
    if (!expression.empty()) {
        try {
            std::smatch matches;
            if (!std::regex_match(name, matches, std::regex(expression)) || matches.size() < 4)
                return false;
            return parse_decimal(matches[1].str(), coordinate.tile) &&
                   parse_decimal(matches[2].str(), coordinate.x) &&
                   parse_decimal(matches[3].str(), coordinate.y);
        } catch (const std::regex_error&) {
            return false;
        }
    }
    std::array<std::string, 3> suffixes{};
    std::size_t end = name.size();
    for (int index = 2; index >= 0; --index) {
        const auto begin = name.rfind(':', end == 0 ? std::string::npos : end - 1);
        if (begin == std::string::npos) return false;
        suffixes[static_cast<std::size_t>(index)] = name.substr(begin + 1, end - begin - 1);
        end = begin;
    }
    return parse_decimal(suffixes[0], coordinate.tile) &&
           parse_decimal(suffixes[1], coordinate.x) && parse_decimal(suffixes[2], coordinate.y);
}

bool is_optical_duplicate_name(const std::string& representative, const std::string& duplicate,
                               const Options& options) {
    OpticalCoordinate left, right;
    if (!parse_optical_coordinate(representative.c_str(), options.read_name_regex, left) ||
        !parse_optical_coordinate(duplicate.c_str(), options.read_name_regex, right) ||
        left.tile != right.tile)
        return false;
    const auto dx = left.x - right.x;
    const auto dy = left.y - right.y;
    const auto distance_squared = dx * dx + dy * dy;
    const auto limit = options.optical_duplicate_pixel_distance;
    return distance_squared <= limit * limit;
}

std::int16_t quality_score(const bam1_t* record) {
    const auto* qualities = bam_get_qual(record);
    std::int32_t score = 0;
    for (int index = 0; index < record->core.l_qseq; ++index) {
        // htsjdk's byte[] is signed.  Keep that detail so the reserved 0xff
        // missing-quality value and out-of-range byte values are excluded in
        // exactly the same way as DuplicateScoringStrategy.
        const auto quality = static_cast<std::int8_t>(qualities[index]);
        if (quality >= 15) score += quality;
    }
    score = std::min<std::int32_t>(score, 16383);
    if ((record->core.flag & BAM_FQCFAIL) != 0) score -= 16384;
    // Picard's default DuplicateScoringStrategy is SUM_OF_BASE_QUALITIES.
    // Mapping quality is deliberately not part of this score: including it
    // changes the representative when two duplicate reads have the same
    // base-quality sum but different MAPQ values.  The Java return type is a
    // short, including the intentional vendor-QC negative range.
    return static_cast<std::int16_t>(score);
}

std::uint64_t estimate_library_size(std::uint64_t read_pairs_examined,
                                    std::uint64_t read_pair_duplicates) {
    if (read_pairs_examined == 0 || read_pair_duplicates == 0 ||
        read_pair_duplicates >= read_pairs_examined)
        return 0;
    const double observed = static_cast<double>(read_pairs_examined);
    const double duplicate_fraction = static_cast<double>(read_pair_duplicates) / observed;
    auto expected_duplicate_fraction = [observed](double library_size) {
        const double ratio = observed / library_size;
        return 1.0 - (1.0 / ratio) * (-std::expm1(-ratio));
    };
    double lower = 1.0;
    double upper = std::max(2.0, observed);
    while (expected_duplicate_fraction(upper) > duplicate_fraction && upper < 1.0e18)
        upper *= 2.0;
    for (int iteration = 0; iteration < 120; ++iteration) {
        const double middle = lower + (upper - lower) * 0.5;
        if (expected_duplicate_fraction(middle) > duplicate_fraction) lower = middle;
        else upper = middle;
    }
    const auto estimate = static_cast<std::uint64_t>(std::llround(upper));
    return estimate == 0 ? 1 : estimate;
}

std::string fragment_name(const bam1_t* record, std::size_t order) {
    const char* qname = bam_get_qname(record);
    std::string name = qname == nullptr ? std::string{} : std::string(qname);
    if ((record->core.flag & BAM_FPAIRED) == 0)
        return "__unpaired_" + std::to_string(order);
    {
        if (name.size() > 2 && name.compare(name.size() - 2, 2, "/1") == 0)
            name.resize(name.size() - 2);
        else if (name.size() > 2 && name.compare(name.size() - 2, 2, "/2") == 0)
            name.resize(name.size() - 2);
    }
    if (name.empty()) name = "__record_" + std::to_string(order);
    return name;
}

std::string pair_identifier(const bam1_t* record) {
    const char* qname = bam_get_qname(record);
    std::string name = qname == nullptr ? std::string{} : std::string(qname);
    if (name.size() > 2 &&
        (name.compare(name.size() - 2, 2, "/1") == 0 ||
         name.compare(name.size() - 2, 2, "/2") == 0))
        name.resize(name.size() - 2);
    return name;
}

std::string pair_lookup_identifier(const bam1_t* record) {
    // Picard's mate cache is keyed by read-group ID plus read name.  Keeping
    // the RG in this transient key prevents equal qnames from different
    // libraries from being joined into one fragment.
    return read_group(record) + '\t' + pair_identifier(record);
}

struct FragmentSummary {
    // Picard stores DuplicateScoringStrategy's result as a signed short.  A
    // paired fragment receives the Java short sum of its two read scores, so
    // retaining the signed value here is important for vendor-QC and very
    // long/high-quality reads (unsigned accumulation changes the winner).
    std::int16_t score = 0;
    std::uint64_t record_count = 0;
    std::uint64_t first_order = std::numeric_limits<std::uint64_t>::max();
    std::string optical_name;
    bool duplicate = false;
    bool optical_duplicate = false;
};

// MarkDuplicates' ReadEnds builder pairs the two alignments by read name and
// uses each alignment's own unclipped five-prime coordinate.  A SAM/BAM record
// only carries the mate's leftmost coordinate, so deriving a pair key from
// mate fields alone is wrong when MC is absent (and for clipped mates even
// when it is present).  Keep the small first-pass observation needed to form
// the Java-equivalent key once both mapped mates have been seen.
struct PairObservation {
    std::string name;
    DuplicateKey fallback_key;
    int tid = -1;
    std::int64_t coordinate = -1;
    bool reverse = false;
    std::int16_t score = 0;
    std::uint64_t first_order = 0;
    std::string optical_name;
};

DuplicateKey paired_duplicate_key(const PairObservation& left,
                                  const PairObservation& right) {
    auto first_tid = left.tid;
    auto first_coordinate = left.coordinate;
    auto first_reverse = left.reverse;
    auto second_tid = right.tid;
    auto second_coordinate = right.coordinate;
    auto second_reverse = right.reverse;
    if (std::tie(second_tid, second_coordinate, second_reverse) <
        std::tie(first_tid, first_coordinate, first_reverse)) {
        std::swap(first_tid, second_tid);
        std::swap(first_coordinate, second_coordinate);
        std::swap(first_reverse, second_reverse);
    }
    const int orientation = static_cast<int>(first_reverse ? 1 : 0) |
                            static_cast<int>(second_reverse ? 2 : 0);
    return {first_tid, first_coordinate, second_tid, second_coordinate,
            orientation, left.fallback_key.library.empty()
                ? right.fallback_key.library : left.fallback_key.library};
}

struct DuplicateGroupSummary {
    bool paired = false;
    std::uint64_t record_count = 0;
    std::map<std::string, FragmentSummary> fragments;
    std::string representative;
};

// ReadEndsMDComparator uses physical location (tile, x, y) as the
// deterministic tie-break after the duplicate score and only then falls back
// to the read's file index.  A map keyed by fragment name has a different
// iteration order, so selecting the lexicographically smallest name is not
// Picard-compatible for equal-score pairs (e.g. x=999 versus x=888).
std::tuple<std::int64_t, std::int64_t, std::int64_t, std::uint64_t>
representative_tie_key(const FragmentSummary& fragment, const Options& options) {
    OpticalCoordinate coordinate{-1, -1, -1};
    (void)parse_optical_coordinate(fragment.optical_name.c_str(), options.read_name_regex, coordinate);
    return {coordinate.tile, coordinate.x, coordinate.y, fragment.first_order};
}

// Picard emits one DuplicationMetrics row per library (LB).  Keep the
// counters in Host memory, not in the duplicate-key spill, because they are
// small and are independent of the second-pass alignment payload.  An empty
// library is retained as an empty key and rendered as "*" in the metrics
// table, matching the aggregate row convention used by the prototype.
struct LibraryMetrics {
    std::uint64_t unpaired_reads_examined = 0;
    std::uint64_t read_pairs_examined = 0;
    std::uint64_t secondary_or_supplementary = 0;
    std::uint64_t unmapped_reads = 0;
    std::uint64_t unpaired_read_duplicates = 0;
    std::uint64_t read_pair_duplicates = 0;
    std::uint64_t read_pair_optical_duplicates = 0;
    std::uint64_t duplicate_records = 0;
};

// Spill format is deliberately private and versioned by the manifest.  It
// stores only duplicate-key metadata; alignment records remain in the input
// BAM/CRAM and are streamed again during the second pass.
struct SpillRun {
    std::string path;
    std::vector<std::uint64_t> offsets;
};

// The alignment payload stays in the input BAM/CRAM; a checkpoint therefore
// only needs the deterministic duplicate summaries and the first-pass
// counters.  The format is text (with quoted strings) so an interrupted
// Nextflow/SLURM retry can be inspected and rejected without trusting a
// partially-written binary control file.
struct SpillCheckpoint {
    std::string input;
    std::string output;
    std::string reference;
    std::string tagging_policy;
    std::string pg_id;
    std::string pg_name;
    std::string pg_version;
    std::string pg_command_line;
    std::uint64_t input_size = 0;
    std::int64_t input_mtime = 0;
    std::size_t external_group_record_limit = 0;
    std::uint64_t input_records = 0;
    std::uint64_t unmapped_records = 0;
    std::uint64_t secondary_or_supplementary = 0;
    bool clear_duplicates = false;
    bool remove_duplicates = false;
    bool remove_sequencing_duplicates = false;
    bool add_pg_tag = true;
    std::vector<std::string> runs;
    std::map<std::string, LibraryMetrics> library_metrics;
    std::map<std::string, DuplicateKey> pair_keys;
};

std::filesystem::path checkpoint_path(const Options& options) {
    const auto identity = options.input + "\n" + options.output + "\n" +
        options.reference + "\n" + options.tagging_policy + "\n" +
        (options.clear_duplicates ? "1" : "0") + (options.remove_duplicates ? "1" : "0") +
        (options.remove_sequencing_duplicates ? "1" : "0") +
        (options.add_pg_tag ? "1" : "0") + options.pg_id + "\n" + options.pg_name +
        "\n" + options.pg_version + "\n" + options.pg_command_line;
    const auto hash = std::hash<std::string>{}(identity);
    std::ostringstream name;
    name << "fastgatk-markduplicates-checkpoint-" << std::hex << hash << ".txt";
    return std::filesystem::path(options.tmp_dir) / name.str();
}

template <typename T>
void write_pod(std::ofstream& stream, const T& value) {
    stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot write MarkDuplicates spill run");
}

template <typename T>
void read_pod(std::ifstream& stream, T& value) {
    stream.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!stream) throw std::runtime_error("BAD_INPUT: truncated MarkDuplicates spill run");
}

void write_string(std::ofstream& stream, const std::string& value) {
    const auto size = static_cast<std::uint64_t>(value.size());
    write_pod(stream, size);
    if (size != 0) stream.write(value.data(), static_cast<std::streamsize>(size));
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot write MarkDuplicates spill string");
}

void read_string(std::ifstream& stream, std::string& value) {
    std::uint64_t size = 0;
    read_pod(stream, size);
    if (size > (1ULL << 30)) throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates spill string length");
    value.resize(static_cast<std::size_t>(size));
    if (size != 0) stream.read(value.data(), static_cast<std::streamsize>(size));
    if (!stream) throw std::runtime_error("BAD_INPUT: truncated MarkDuplicates spill string");
}

void write_spill_entry(std::ofstream& stream, const DuplicateKey& key,
                       const DuplicateGroupSummary& summary) {
    write_pod(stream, key.tid);
    write_pod(stream, key.pos);
    write_pod(stream, key.mtid);
    write_pod(stream, key.mpos);
    write_pod(stream, key.orientation);
    write_string(stream, key.library);
    const std::uint8_t paired = summary.paired ? 1 : 0;
    write_pod(stream, paired);
    write_pod(stream, summary.record_count);
    write_string(stream, summary.representative);
    const auto count = static_cast<std::uint64_t>(summary.fragments.size());
    write_pod(stream, count);
    for (const auto& [name, fragment] : summary.fragments) {
        write_string(stream, name);
        write_pod(stream, fragment.score);
        write_pod(stream, fragment.record_count);
        write_pod(stream, fragment.first_order);
        write_string(stream, fragment.optical_name);
        const std::uint8_t duplicate = fragment.duplicate ? 1 : 0;
        const std::uint8_t optical = fragment.optical_duplicate ? 1 : 0;
        write_pod(stream, duplicate);
        write_pod(stream, optical);
    }
}

void read_spill_entry(std::ifstream& stream, DuplicateKey& key,
                      DuplicateGroupSummary& summary) {
    read_pod(stream, key.tid);
    read_pod(stream, key.pos);
    read_pod(stream, key.mtid);
    read_pod(stream, key.mpos);
    read_pod(stream, key.orientation);
    read_string(stream, key.library);
    std::uint8_t paired = 0;
    read_pod(stream, paired);
    summary = {};
    summary.paired = paired != 0;
    read_pod(stream, summary.record_count);
    read_string(stream, summary.representative);
    std::uint64_t count = 0;
    read_pod(stream, count);
    if (count > (1ULL << 30)) throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates spill group size");
    for (std::uint64_t index = 0; index < count; ++index) {
        std::string name;
        FragmentSummary fragment;
        read_string(stream, name);
        read_pod(stream, fragment.score);
        read_pod(stream, fragment.record_count);
        read_pod(stream, fragment.first_order);
        read_string(stream, fragment.optical_name);
        std::uint8_t duplicate = 0;
        std::uint8_t optical = 0;
        read_pod(stream, duplicate);
        read_pod(stream, optical);
        fragment.duplicate = duplicate != 0;
        fragment.optical_duplicate = optical != 0;
        summary.fragments.emplace(std::move(name), std::move(fragment));
    }
}

void write_checkpoint(const std::filesystem::path& path, const Options& options,
                      std::size_t external_group_record_limit,
                      std::uint64_t input_records, std::uint64_t unmapped_records,
                      std::uint64_t secondary_or_supplementary,
                      const std::map<std::string, LibraryMetrics>& library_metrics,
                      const std::vector<SpillRun>& runs,
                      const std::map<std::string, DuplicateKey>& pair_keys) {
    std::error_code stat_error;
    const auto input_size = std::filesystem::file_size(options.input, stat_error);
    if (stat_error)
        throw std::runtime_error("BAD_INPUT: cannot stat MarkDuplicates input for checkpoint: " + options.input);
    const auto input_mtime = std::filesystem::last_write_time(options.input, stat_error)
        .time_since_epoch().count();
    if (stat_error)
        throw std::runtime_error("BAD_INPUT: cannot timestamp MarkDuplicates input for checkpoint: " + options.input);
    const auto temporary = path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::trunc);
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot write MarkDuplicates checkpoint");
    // V3 adds repaired paired-key metadata needed by resumed no-MC inputs;
    // V2 added fragment first-order metadata to spill entries and preserved
    // REMOVE_SEQUENCING_DUPLICATES in the restart signature.
    stream << "FASTGATK_MARKDUPLICATES_CHECKPOINT_V3\n"
           << "input " << std::quoted(options.input) << "\n"
           << "output " << std::quoted(options.output) << "\n"
           << "reference " << std::quoted(options.reference) << "\n"
           << "tagging_policy " << std::quoted(options.tagging_policy) << "\n"
           << "pg_id " << std::quoted(options.pg_id) << "\n"
           << "pg_name " << std::quoted(options.pg_name) << "\n"
           << "pg_version " << std::quoted(options.pg_version) << "\n"
           << "pg_command_line " << std::quoted(options.pg_command_line) << "\n"
           << "input_size " << input_size << "\n"
           << "input_mtime " << input_mtime << "\n"
           << "external_group_record_limit " << external_group_record_limit << "\n"
           << "input_records " << input_records << "\n"
           << "unmapped_records " << unmapped_records << "\n"
           << "secondary_or_supplementary " << secondary_or_supplementary << "\n"
           << "clear_duplicates " << (options.clear_duplicates ? 1 : 0) << "\n"
           << "remove_duplicates " << (options.remove_duplicates ? 1 : 0) << "\n"
           << "remove_sequencing_duplicates " << (options.remove_sequencing_duplicates ? 1 : 0) << "\n"
           << "add_pg_tag " << (options.add_pg_tag ? 1 : 0) << "\n"
           << "runs " << runs.size() << "\n";
    for (const auto& run : runs) stream << "run " << std::quoted(run.path) << "\n";
    stream << "library_metrics " << library_metrics.size() << "\n";
    for (const auto& [library, metrics] : library_metrics) {
        stream << "library " << std::quoted(library) << ' '
               << metrics.unpaired_reads_examined << ' '
               << metrics.read_pairs_examined << ' '
               << metrics.secondary_or_supplementary << ' '
               << metrics.unmapped_reads << ' '
               << metrics.unpaired_read_duplicates << ' '
               << metrics.read_pair_duplicates << ' '
               << metrics.read_pair_optical_duplicates << ' '
               << metrics.duplicate_records << '\n';
    }
    stream << "pair_keys " << pair_keys.size() << "\n";
    for (const auto& [name, key] : pair_keys) {
        stream << "pair_key " << std::quoted(name) << ' '
               << key.tid << ' ' << key.pos << ' ' << key.mtid << ' '
               << key.mpos << ' ' << key.orientation << ' '
               << std::quoted(key.library) << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot finalize MarkDuplicates checkpoint");
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) {
        std::filesystem::remove(path, error);
        error.clear();
        std::filesystem::rename(temporary, path, error);
    }
    if (error) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot install MarkDuplicates checkpoint: " + path.string());
}

SpillCheckpoint read_checkpoint(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open MarkDuplicates checkpoint: " + path.string());
    std::string version;
    if (!(stream >> version) || version != "FASTGATK_MARKDUPLICATES_CHECKPOINT_V3")
        throw std::runtime_error("BAD_INPUT: unsupported MarkDuplicates checkpoint version");
    SpillCheckpoint checkpoint;
    std::string key;
    std::size_t expected_runs = 0;
    std::size_t expected_libraries = 0;
    std::size_t expected_pair_keys = 0;
    while (stream >> key) {
        if (key == "input") stream >> std::quoted(checkpoint.input);
        else if (key == "output") stream >> std::quoted(checkpoint.output);
        else if (key == "reference") stream >> std::quoted(checkpoint.reference);
        else if (key == "tagging_policy") stream >> std::quoted(checkpoint.tagging_policy);
        else if (key == "pg_id") stream >> std::quoted(checkpoint.pg_id);
        else if (key == "pg_name") stream >> std::quoted(checkpoint.pg_name);
        else if (key == "pg_version") stream >> std::quoted(checkpoint.pg_version);
        else if (key == "pg_command_line") stream >> std::quoted(checkpoint.pg_command_line);
        else if (key == "input_size") stream >> checkpoint.input_size;
        else if (key == "input_mtime") stream >> checkpoint.input_mtime;
        else if (key == "external_group_record_limit") stream >> checkpoint.external_group_record_limit;
        else if (key == "input_records") stream >> checkpoint.input_records;
        else if (key == "unmapped_records") stream >> checkpoint.unmapped_records;
        else if (key == "secondary_or_supplementary") stream >> checkpoint.secondary_or_supplementary;
        else if (key == "clear_duplicates") { int value = 0; stream >> value; checkpoint.clear_duplicates = value != 0; }
        else if (key == "remove_duplicates") { int value = 0; stream >> value; checkpoint.remove_duplicates = value != 0; }
        else if (key == "remove_sequencing_duplicates") { int value = 0; stream >> value; checkpoint.remove_sequencing_duplicates = value != 0; }
        else if (key == "add_pg_tag") { int value = 1; stream >> value; checkpoint.add_pg_tag = value != 0; }
        else if (key == "runs") {
            stream >> expected_runs;
            if (expected_runs > (1ULL << 20)) throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates checkpoint run count");
            checkpoint.runs.reserve(expected_runs);
        } else if (key == "run") {
            std::string value;
            stream >> std::quoted(value);
            checkpoint.runs.push_back(std::move(value));
        } else if (key == "library_metrics") {
            stream >> expected_libraries;
            if (expected_libraries > (1ULL << 20)) throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates checkpoint library count");
        } else if (key == "library") {
            std::string name;
            LibraryMetrics metrics;
            stream >> std::quoted(name)
                   >> metrics.unpaired_reads_examined
                   >> metrics.read_pairs_examined
                   >> metrics.secondary_or_supplementary
                   >> metrics.unmapped_reads
                   >> metrics.unpaired_read_duplicates
                   >> metrics.read_pair_duplicates
                   >> metrics.read_pair_optical_duplicates
                   >> metrics.duplicate_records;
            checkpoint.library_metrics.emplace(std::move(name), metrics);
        } else if (key == "pair_keys") {
            stream >> expected_pair_keys;
            if (expected_pair_keys > (1ULL << 30))
                throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates checkpoint pair-key count");
            checkpoint.pair_keys.clear();
        } else if (key == "pair_key") {
            std::string name;
            DuplicateKey value;
            stream >> std::quoted(name) >> value.tid >> value.pos >> value.mtid
                   >> value.mpos >> value.orientation >> std::quoted(value.library);
            checkpoint.pair_keys.emplace(std::move(name), std::move(value));
        } else {
            throw std::runtime_error("BAD_INPUT: unknown MarkDuplicates checkpoint field: " + key);
        }
        if (!stream) throw std::runtime_error("BAD_INPUT: truncated MarkDuplicates checkpoint");
    }
    if (checkpoint.input.empty() || checkpoint.output.empty() || checkpoint.tagging_policy.empty() ||
        checkpoint.external_group_record_limit == 0 || checkpoint.runs.size() != expected_runs ||
        checkpoint.library_metrics.size() != expected_libraries ||
        checkpoint.pair_keys.size() != expected_pair_keys)
        throw std::runtime_error("BAD_INPUT: incomplete MarkDuplicates checkpoint");
    for (const auto& run : checkpoint.runs) {
        if (!file_complete(run))
            throw std::runtime_error("BAD_INPUT: missing or empty MarkDuplicates spill run: " + run);
    }
    return checkpoint;
}

SpillRun load_spill_run(const std::string& path) {
    SpillRun run;
    run.path = path;
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open MarkDuplicates spill run: " + path);
    while (true) {
        const auto offset = stream.tellg();
        if (offset < 0) throw std::runtime_error("BAD_INPUT: invalid MarkDuplicates spill run offset: " + path);
        const int next = stream.peek();
        if (next == std::char_traits<char>::eof()) {
            if (stream.eof()) break;
            throw std::runtime_error("BAD_INPUT: cannot scan MarkDuplicates spill run: " + path);
        }
        DuplicateKey key;
        DuplicateGroupSummary summary;
        run.offsets.push_back(static_cast<std::uint64_t>(offset));
        read_spill_entry(stream, key, summary);
    }
    if (run.offsets.empty()) throw std::runtime_error("BAD_INPUT: empty MarkDuplicates spill run: " + path);
    return run;
}

SpillRun spill_groups(const std::map<DuplicateKey, DuplicateGroupSummary, DuplicateKeyLess>& groups,
                      const std::string& directory, std::size_t run_number) {
    SpillRun run;
    run.path = (std::filesystem::path(directory) /
                ("fastgatk-markduplicates-run-" + std::to_string(run_number) + ".bin")).string();
    std::ofstream stream(run.path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot create MarkDuplicates spill run: " + run.path);
    for (const auto& [key, summary] : groups) {
        const auto offset = static_cast<std::uint64_t>(stream.tellp());
        run.offsets.push_back(offset);
        write_spill_entry(stream, key, summary);
    }
    stream.flush();
    if (!stream || !file_complete(run.path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: MarkDuplicates spill run is incomplete");
    return run;
}

void merge_group_summary(DuplicateGroupSummary& destination,
                         const DuplicateGroupSummary& source) {
    destination.paired = destination.paired || source.paired;
    destination.record_count += source.record_count;
    for (const auto& [name, fragment] : source.fragments) {
        auto& target = destination.fragments[name];
        target.score = static_cast<std::int16_t>(
            static_cast<std::int32_t>(target.score) + static_cast<std::int32_t>(fragment.score));
        target.record_count += fragment.record_count;
        target.first_order = std::min(target.first_order, fragment.first_order);
        if (target.optical_name.empty()) target.optical_name = fragment.optical_name;
        target.duplicate = target.duplicate || fragment.duplicate;
        target.optical_duplicate = target.optical_duplicate || fragment.optical_duplicate;
    }
}

bool lookup_spill_run(const SpillRun& run, const DuplicateKey& wanted,
                      DuplicateGroupSummary& result) {
    if (run.offsets.empty()) return false;
    std::ifstream stream(run.path, std::ios::binary);
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open MarkDuplicates spill run: " + run.path);
    std::size_t low = 0;
    std::size_t high = run.offsets.size();
    while (low < high) {
        const auto middle = low + (high - low) / 2;
        stream.clear();
        stream.seekg(static_cast<std::streamoff>(run.offsets[middle]));
        if (!stream) throw std::runtime_error("BAD_INPUT: cannot seek MarkDuplicates spill run");
        DuplicateKey key;
        DuplicateGroupSummary summary;
        read_spill_entry(stream, key, summary);
        if (key_less(key, wanted)) {
            low = middle + 1;
        } else if (key_less(wanted, key)) {
            high = middle;
        } else {
            result = std::move(summary);
            return true;
        }
    }
    return false;
}

bool lookup_spills(const std::vector<SpillRun>& runs, const DuplicateKey& wanted,
                   DuplicateGroupSummary& result) {
    bool found = false;
    result = {};
    for (const auto& run : runs) {
        DuplicateGroupSummary candidate;
        if (lookup_spill_run(run, wanted, candidate)) {
            if (!found) result = std::move(candidate);
            else merge_group_summary(result, candidate);
            found = true;
        }
    }
    return found;
}

void cleanup_spills(const std::vector<SpillRun>& runs) {
    for (const auto& run : runs) {
        std::error_code error;
        std::filesystem::remove(run.path, error);
    }
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    samFile* input = nullptr;
    samFile* second_input = nullptr;
    samFile* output = nullptr;
    sam_hdr_t* header = nullptr;
    sam_hdr_t* second_header = nullptr;
    bam1_t* record = nullptr;
    std::map<DuplicateKey, DuplicateGroupSummary, DuplicateKeyLess> groups;
    // Keys for paired records whose mate fields are insufficient to recreate
    // the Java ReadEnds pair (notably records without MC).  This map is only
    // populated when the independently-computed pair key differs from both
    // record-local fallbacks; ordinary MC-bearing input stays streaming.
    std::map<std::string, DuplicateKey> pair_keys;
    std::map<std::string, PairObservation> pending_pairs;
    std::vector<SpillRun> spill_runs;
    const bool external_spill = !options.tmp_dir.empty();
    std::size_t spill_run_number = 0;
    std::uint64_t input_records = 0, output_records = 0, duplicate_records = 0;
    std::uint64_t optical_duplicate_records = 0;
    std::uint64_t paired_examined = 0, unpaired_examined = 0, paired_duplicates = 0, unpaired_duplicates = 0;
    std::uint64_t paired_duplicate_records = 0, optical_duplicate_pairs = 0;
    std::uint64_t secondary_or_supplementary = 0, unmapped_records = 0;
    std::map<std::string, LibraryMetrics> library_metrics;
    std::filesystem::path checkpoint;
    bool resumed = false;
    // Picard writes a duplicate-set histogram for paired fragments.  The
    // first element is the number of sets of a given size, followed by the
    // optical and non-optical subsets.
    std::map<std::size_t, std::array<std::uint64_t, 3>> duplicate_set_histogram;
    try {
        input = sam_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open SAM/BAM/CRAM: " + options.input);
        if (!options.reference.empty() && hts_set_fai_filename(input, options.reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure input reference");
        header = sam_hdr_read(input);
        if (!header) throw std::runtime_error("BAD_INPUT: cannot read SAM/BAM header");
        if (!options.tmp_dir.empty()) {
            std::error_code error;
            std::filesystem::create_directories(options.tmp_dir, error);
            if (error || !std::filesystem::is_directory(options.tmp_dir))
                throw std::runtime_error("RESOURCE_EXHAUSTED: cannot create MarkDuplicates tmp directory: " + options.tmp_dir);
        }
        const auto safe_budget = resources.safe_memory_budget_bytes();
        const std::size_t external_group_record_limit = options.max_records_in_memory != 0
            ? options.max_records_in_memory
            : (safe_budget == 0 ? 100000 : std::max<std::size_t>(1,
                static_cast<std::size_t>(safe_budget / (4 * (sizeof(bam1_t) + 256)))));
        const std::size_t external_key_limit = options.max_records_in_memory != 0
            ? options.max_records_in_memory : 100000;
        kstring_t sort_order = {0, 0, nullptr};
        const bool coordinate_sorted = options.assume_sorted ||
            (sam_hdr_find_tag_hd(header, "SO", &sort_order) == 0 &&
             sort_order.s != nullptr && std::strcmp(sort_order.s, "coordinate") == 0);
        free(sort_order.s);
        if (!coordinate_sorted) throw std::runtime_error("BAD_INPUT: MarkDuplicates requires coordinate-sorted input or --assume-sorted");
        if (options.add_pg_tag) {
            const auto command_line = options.pg_command_line.empty()
                ? std::string("fastgatk MarkDuplicates") : options.pg_command_line;
            if (sam_hdr_add_pg(header, options.pg_name.c_str(),
                               "ID", options.pg_id.c_str(),
                               "VN", options.pg_version.c_str(),
                               "CL", command_line.c_str(), NULL) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot add MarkDuplicates @PG header");
        }
        // External-spill runs publish a checkpoint after their first pass.  A
        // normal invocation never consumes an old checkpoint implicitly;
        // callers must opt in with --resume-spill so stale retry state cannot
        // silently change the requested execution.
        checkpoint = external_spill ? checkpoint_path(options) : std::filesystem::path{};
        std::error_code checkpoint_error;
        const bool resume_requested = options.resume_spill && !checkpoint.empty() &&
            std::filesystem::is_regular_file(checkpoint, checkpoint_error) && !checkpoint_error;
        if (options.resume_spill && !resume_requested)
            throw std::runtime_error("BAD_INPUT: --resume-spill requested but no MarkDuplicates checkpoint exists");
        auto add_fragment = [&](const DuplicateKey& key, const std::string& fragment,
                                std::int16_t score, std::uint64_t record_count,
                                std::uint64_t first_order, const std::string& optical_name,
                                bool paired) {
            auto& summary = groups[key];
            summary.paired = summary.paired || paired;
            summary.record_count += record_count;
            auto& fragment_summary = summary.fragments[fragment];
            fragment_summary.record_count += record_count;
            fragment_summary.score = static_cast<std::int16_t>(
                static_cast<std::int32_t>(fragment_summary.score) +
                static_cast<std::int32_t>(score));
            fragment_summary.first_order = std::min(fragment_summary.first_order, first_order);
            if (fragment_summary.optical_name.empty()) fragment_summary.optical_name = optical_name;
            if (options.max_records_in_memory != 0 &&
                summary.record_count > options.max_records_in_memory && !external_spill)
                throw std::runtime_error("RESOURCE_EXHAUSTED: duplicate group exceeds --max-records-in-memory; use --tmp-dir for external spill or explicit Picard fallback");
            const auto budget = resources.safe_memory_budget_bytes();
            if (budget != 0 && summary.record_count * (sizeof(bam1_t) + 256) > budget / 4 && !external_spill)
                throw std::runtime_error("RESOURCE_EXHAUSTED: duplicate group exceeds safe memory budget; use --tmp-dir for external spill or explicit Picard fallback");
            if (external_spill && summary.record_count >= external_group_record_limit) {
                std::map<DuplicateKey, DuplicateGroupSummary, DuplicateKeyLess> chunk;
                chunk.emplace(key, std::move(summary));
                groups.erase(key);
                spill_runs.push_back(spill_groups(chunk, options.tmp_dir, spill_run_number++));
            }
            if (external_spill && groups.size() >= external_key_limit) {
                spill_runs.push_back(spill_groups(groups, options.tmp_dir, spill_run_number++));
                groups.clear();
            }
        };
        if (resume_requested) {
            const auto saved = read_checkpoint(checkpoint);
            if (saved.input != options.input || saved.output != options.output ||
                saved.reference != options.reference || saved.tagging_policy != options.tagging_policy ||
                saved.pg_id != options.pg_id || saved.pg_name != options.pg_name ||
                saved.pg_version != options.pg_version || saved.pg_command_line != options.pg_command_line ||
                saved.clear_duplicates != options.clear_duplicates ||
                saved.remove_duplicates != options.remove_duplicates ||
                saved.remove_sequencing_duplicates != options.remove_sequencing_duplicates ||
                saved.add_pg_tag != options.add_pg_tag ||
                saved.external_group_record_limit != external_group_record_limit)
                throw std::runtime_error("BAD_INPUT: MarkDuplicates checkpoint does not match input/output/options/budget");
            std::error_code stat_error;
            const auto current_size = std::filesystem::file_size(options.input, stat_error);
            const auto current_mtime = stat_error ? std::int64_t{0} :
                std::filesystem::last_write_time(options.input, stat_error).time_since_epoch().count();
            if (stat_error || current_size != saved.input_size || current_mtime != saved.input_mtime)
                throw std::runtime_error("BAD_INPUT: MarkDuplicates input changed since spill checkpoint");
            for (const auto& path : saved.runs) spill_runs.push_back(load_spill_run(path));
            input_records = saved.input_records;
            unmapped_records = saved.unmapped_records;
            secondary_or_supplementary = saved.secondary_or_supplementary;
            library_metrics = saved.library_metrics;
            pair_keys = saved.pair_keys;
            spill_run_number = spill_runs.size();
            resumed = true;
            if (sam_close(input) != 0)
                throw std::runtime_error("cannot close input while resuming MarkDuplicates");
            input = nullptr;
        } else {
            record = bam_init1();
            if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed");
            // First pass: collect only duplicate-key/fragment metadata.  This
            // is deliberately separate from output so a pair whose mates are
            // far apart in coordinate order is still assigned one consensus.
            int last_tid = -1;
            std::int64_t last_pos = -1;
            bool have_last_coordinate = false;
            while (sam_read1(input, header, record) >= 0) {
                ++input_records;
                if (!options.assume_sorted && record->core.tid >= 0 && record->core.pos >= 0) {
                    if (have_last_coordinate &&
                        (record->core.tid < last_tid ||
                         (record->core.tid == last_tid && record->core.pos < last_pos)))
                        throw std::runtime_error("BAD_INPUT: MarkDuplicates input is not coordinate sorted");
                    last_tid = record->core.tid;
                    last_pos = record->core.pos;
                    have_last_coordinate = true;
                }
                if ((record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0) {
                    auto& metrics = library_metrics[library_for_read_group(header, read_group(record))];
                    if ((record->core.flag & BAM_FUNMAP) != 0) ++unmapped_records;
                    if ((record->core.flag & BAM_FUNMAP) != 0) ++metrics.unmapped_reads;
                    if ((record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0)
                        ++secondary_or_supplementary;
                    if ((record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0)
                        ++metrics.secondary_or_supplementary;
                    continue;
                }
                const bool both_mapped_pair =
                    (record->core.flag & BAM_FPAIRED) != 0 &&
                    (record->core.flag & BAM_FMUNMAP) == 0 &&
                    record->core.mtid >= 0;
                if (both_mapped_pair) {
                    const auto name = pair_identifier(record);
                    const auto lookup_name = pair_lookup_identifier(record);
                    const auto pending = pending_pairs.find(lookup_name);
                    PairObservation observation;
                    observation.name = name;
                    observation.fallback_key = duplicate_key(record, header);
                    observation.tid = record->core.tid;
                    observation.coordinate = five_prime_position(record);
                    observation.reverse = (record->core.flag & BAM_FREVERSE) != 0;
                    observation.score = quality_score(record);
                    observation.first_order = input_records;
                    if (const auto* qname = bam_get_qname(record); qname != nullptr)
                        observation.optical_name = qname;
                    if (pending == pending_pairs.end()) {
                        pending_pairs.emplace(lookup_name, std::move(observation));
                        const auto budget = resources.safe_memory_budget_bytes();
                        if (budget != 0 && pending_pairs.size() * (sizeof(PairObservation) + 96) > budget / 4)
                            throw std::runtime_error("RESOURCE_EXHAUSTED: paired observation map exceeds safe memory budget");
                    } else {
                        const auto first = std::move(pending->second);
                        pending_pairs.erase(pending);
                        const auto key = paired_duplicate_key(first, observation);
                        const auto fragment = name.empty()
                            ? fragment_name(record, static_cast<std::size_t>(input_records)) : name;
                        add_fragment(key, fragment,
                                     static_cast<std::int16_t>(
                                         static_cast<std::int32_t>(first.score) +
                                         static_cast<std::int32_t>(observation.score)),
                                     2, std::min(first.first_order, observation.first_order),
                                     first.optical_name.empty() ? observation.optical_name : first.optical_name,
                                     true);
                        // Only records for which the old record-local key was
                        // insufficient need a qname lookup during pass two.
                        // This keeps the usual MC-bearing path streaming while
                        // making no-MC/clipped mates exact on the bounded map.
                        if (!(first.fallback_key == key && observation.fallback_key == key)) {
                            pair_keys[lookup_name] = key;
                            const auto budget = resources.safe_memory_budget_bytes();
                            if (budget != 0 && pair_keys.size() * (sizeof(DuplicateKey) + 96) > budget / 4)
                                throw std::runtime_error("RESOURCE_EXHAUSTED: paired-key repair map exceeds safe memory budget");
                        }
                    }
                } else {
                    add_fragment(duplicate_key(record, header),
                                 fragment_name(record, static_cast<std::size_t>(input_records)),
                                 quality_score(record), 1, input_records,
                                 bam_get_qname(record) == nullptr ? std::string{} : std::string(bam_get_qname(record)),
                                 false);
                }
            }
            // An input can contain a mapped read whose mate record is absent.
            // Picard retains that read as a fragment using its record-local
            // ReadEnds key; do the same after the paired observations close.
            for (const auto& [name, observation] : pending_pairs) {
                add_fragment(observation.fallback_key, observation.name,
                             observation.score, 1, observation.first_order,
                             observation.optical_name, true);
                pair_keys[name] = observation.fallback_key;
            }
            pending_pairs.clear();
            bam_destroy1(record); record = nullptr;
            if (sam_close(input) != 0) throw std::runtime_error("cannot close MarkDuplicates first pass input");
            input = nullptr;

            if (external_spill && !groups.empty()) {
                spill_runs.push_back(spill_groups(groups, options.tmp_dir, spill_run_number++));
                groups.clear();
            }
            if (!checkpoint.empty())
            write_checkpoint(checkpoint, options, external_group_record_limit,
                                 input_records, unmapped_records,
                                 secondary_or_supplementary, library_metrics, spill_runs,
                                 pair_keys);
            // Test-only fault injection used by the restart contract verifier.
            // It models a preemption/OOM after the durable first-pass state has
            // been published, leaving the checkpoint and spill runs available
            // for a later --resume-spill invocation.
            if (!checkpoint.empty()) {
                const char* fail_after_checkpoint = std::getenv(
                    "FASTGATK_MARKDUPLICATES_FAIL_AFTER_CHECKPOINT");
                if (fail_after_checkpoint != nullptr &&
                    std::string(fail_after_checkpoint) == "1")
                    throw std::runtime_error(
                        "RESOURCE_EXHAUSTED: injected MarkDuplicates checkpoint interruption");
            }
        }

        // Decide representatives and compute metrics before the streaming
        // second pass.  No alignment record is retained in memory.  The same
        // accounting function is used for in-memory and spilled groups so
        // per-library metrics cannot diverge between resource modes.
        const auto account_group = [&](const DuplicateKey& key,
                                       DuplicateGroupSummary& summary) {
            if (summary.fragments.empty()) return;
            auto representative = summary.fragments.begin();
            for (auto it = std::next(summary.fragments.begin()); it != summary.fragments.end(); ++it) {
                if (it->second.score > representative->second.score ||
                    (it->second.score == representative->second.score &&
                     representative_tie_key(it->second, options) <
                         representative_tie_key(representative->second, options)))
                    representative = it;
            }
            summary.representative = representative->first;
            auto& library = library_metrics[key.library];
            std::size_t optical_count = 0;
            if (summary.paired) {
                paired_examined += summary.fragments.size();
                library.read_pairs_examined += summary.fragments.size();
            } else {
                unpaired_examined += summary.fragments.size();
                library.unpaired_reads_examined += summary.fragments.size();
            }
            for (auto& [name, fragment] : summary.fragments) {
                fragment.duplicate = name != summary.representative;
                if (!fragment.duplicate) continue;
                if (summary.paired) {
                    ++paired_duplicates;
                    paired_duplicate_records += fragment.record_count;
                    ++library.read_pair_duplicates;
                } else {
                    ++unpaired_duplicates;
                    ++library.unpaired_read_duplicates;
                }
                library.duplicate_records += summary.paired ? fragment.record_count : 1;
                // Picard reports optical duplicate clusters for paired
                // fragments.  A coordinate-colliding unpaired duplicate may
                // have a sequencer-style name, but it remains an ordinary
                // library duplicate and must not be tagged DT:Z:SQ.
                fragment.optical_duplicate = summary.paired &&
                    is_optical_duplicate_name(
                        representative->second.optical_name, fragment.optical_name, options);
                if (fragment.optical_duplicate) {
                    ++optical_count;
                    if (summary.paired) {
                        ++optical_duplicate_pairs;
                        ++library.read_pair_optical_duplicates;
                    }
                    optical_duplicate_records += summary.paired ? fragment.record_count : 1;
                }
            }
            if (summary.paired && summary.fragments.size() > 1) {
                bool has_duplicate = false;
                for (const auto& [name, fragment] : summary.fragments) {
                    (void)name;
                    has_duplicate = has_duplicate || fragment.duplicate;
                }
                if (has_duplicate) {
                    auto& histogram = duplicate_set_histogram[summary.fragments.size()];
                    ++histogram[0];
                    // Picard's optical histogram bin includes the
                    // representative (optical_count + 1), while the
                    // non-optical histogram contains the remainder of the
                    // duplicate set (set_size - optical_count).  They are
                    // emitted as three overlaid histograms, so their bins
                    // need not share the all_sets key.
                    if (optical_count != 0)
                        ++duplicate_set_histogram[optical_count + 1][1];
                    ++duplicate_set_histogram[summary.fragments.size() - optical_count][2];
                }
            }
        };

        for (auto& [key, summary] : groups) {
            account_group(key, summary);
        }

        std::set<DuplicateKey, DuplicateKeyLess> external_seen_keys;
        const auto external_fragment_decision = [&](const DuplicateKey& key,
                                                     const std::string& fragment_name_value) {
            DuplicateGroupSummary summary;
            if (!lookup_spills(spill_runs, key, summary))
                throw std::runtime_error("INTERNAL_ERROR: duplicate key missing in spilled metadata");
            if (summary.fragments.empty())
                throw std::runtime_error("INTERNAL_ERROR: spilled duplicate group has no fragments");
            auto representative = summary.fragments.begin();
            for (auto it = std::next(summary.fragments.begin()); it != summary.fragments.end(); ++it) {
                if (it->second.score > representative->second.score ||
                    (it->second.score == representative->second.score &&
                     representative_tie_key(it->second, options) <
                         representative_tie_key(representative->second, options)))
                    representative = it;
            }
            if (external_seen_keys.insert(key).second)
                account_group(key, summary);
            const auto fragment = summary.fragments.find(fragment_name_value);
            if (fragment == summary.fragments.end())
                throw std::runtime_error("INTERNAL_ERROR: spilled duplicate fragment missing in second pass");
            const bool duplicate = fragment->first != representative->first;
            const bool optical = duplicate && summary.paired &&
                is_optical_duplicate_name(
                    representative->second.optical_name, fragment->second.optical_name, options);
            return std::pair<bool, bool>{duplicate, optical};
        };

        second_input = sam_open(options.input.c_str(), "r");
        if (!second_input) throw std::runtime_error("BAD_INPUT: cannot reopen MarkDuplicates input: " + options.input);
        if (!options.reference.empty() && hts_set_fai_filename(second_input, options.reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure second-pass input reference");
        second_header = sam_hdr_read(second_input);
        if (!second_header) throw std::runtime_error("BAD_INPUT: cannot read second-pass SAM/BAM header");
        const char* output_mode = suffix(options.output, ".cram") ? "wc" :
            suffix(options.output, ".sam") ? "w" : "wb";
        output = sam_open(options.output.c_str(), output_mode);
        if (!output) throw std::runtime_error("cannot open MarkDuplicates output: " + options.output);
        if (!options.reference.empty() && hts_set_fai_filename(output, options.reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure output reference");
        if (sam_hdr_write(output, header) < 0) throw std::runtime_error("cannot write MarkDuplicates header");
        record = bam_init1();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed for second pass");
        std::uint64_t second_order = 0;
        while (sam_read1(second_input, second_header, record) >= 0) {
            ++second_order;
            // Picard MarkDuplicates runs with CLEAR_DT=true by default and
            // starts each record from a clean duplicate state.  This applies
            // to representatives, secondary/supplementary records and reads
            // that were already marked by an upstream pass; otherwise a
            // stale BAM_FDUP/DT would survive on a representative or on a
            // record excluded from the current duplicate-key pass.  Keep the
            // historical --clear-duplicates switch as an idempotent alias,
            // but make the default match the pinned Picard behavior too.
            record->core.flag &= ~BAM_FDUP;
            auto* duplicate_tag = bam_aux_get(record, "DT");
            if (duplicate_tag != nullptr) bam_aux_del(record, duplicate_tag);
            bool duplicate = false;
            bool optical = false;
            if ((record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0) {
                const auto identifier = pair_lookup_identifier(record);
                const auto repaired = pair_keys.find(identifier);
                const auto key = repaired == pair_keys.end()
                    ? duplicate_key(record, second_header) : repaired->second;
                const auto fragment = fragment_name(record, static_cast<std::size_t>(second_order));
                if (external_spill) {
                    const auto decision = external_fragment_decision(key, fragment);
                    duplicate = decision.first;
                    optical = decision.second;
                } else {
                    const auto group_it = groups.find(key);
                    if (group_it == groups.end()) throw std::runtime_error("INTERNAL_ERROR: duplicate key missing in second pass");
                    const auto fragment_it = group_it->second.fragments.find(fragment);
                    if (fragment_it == group_it->second.fragments.end())
                        throw std::runtime_error("INTERNAL_ERROR: duplicate fragment missing in second pass");
                    duplicate = fragment_it->second.duplicate;
                    optical = fragment_it->second.optical_duplicate;
                }
            }
            if (duplicate) {
                record->core.flag |= BAM_FDUP;
                ++duplicate_records;
                const bool tagging_all = options.tagging_policy == "All" || options.tagging_policy == "all";
                const bool tagging_optical_only = options.tagging_policy == "OpticalOnly" ||
                    options.tagging_policy == "opticalonly";
                // Picard's OpticalOnly policy tags only the optical members
                // of a duplicate set (DT:Z:SQ); ordinary library duplicates
                // intentionally carry no DT tag.  Keep the decision at the
                // output boundary after optical classification so All and
                // OpticalOnly share exactly the same duplicate accounting.
                if (tagging_all || (tagging_optical_only && optical)) {
                    const char* tag = optical ? "SQ" : "LB";
                    if (bam_aux_update_str(record, "DT", 3, tag) != 0)
                        throw std::runtime_error("cannot tag MarkDuplicates record");
                }
            }
            // Picard's ADD_PG_TAG_TO_READS also attaches the newly-created
            // program ID to every emitted alignment (not just the @PG header).
            // Replace an upstream PG value while retaining it in the header's
            // provenance chain; this is what htsjdk's SAMRecordProgramGroup
            // writer does for a MarkDuplicates pass.
            if (options.add_pg_tag) {
                if (options.pg_id.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: MarkDuplicates @PG ID is too long");
                if (bam_aux_update_str(record, "PG",
                                       static_cast<int>(options.pg_id.size() + 1),
                                       options.pg_id.c_str()) != 0)
                    throw std::runtime_error("cannot tag MarkDuplicates record with PG");
            }
            const bool omit_record = duplicate &&
                (options.remove_duplicates ||
                 (options.remove_sequencing_duplicates && optical));
            if (!omit_record) {
                if (sam_write1(output, header, record) < 0)
                    throw std::runtime_error("cannot write MarkDuplicates record");
                ++output_records;
            }
        }
        bam_destroy1(record); record = nullptr;
        if (sam_close(second_input) != 0) throw std::runtime_error("cannot close MarkDuplicates second-pass input");
        second_input = nullptr;
        if (sam_close(output) != 0) throw std::runtime_error("cannot close MarkDuplicates output");
        output = nullptr;

        // Picard excludes optical duplicate pairs from both the examined and
        // duplicate populations before solving the Lander-Waterman estimate.
        const auto estimated_library_size = estimate_library_size(
            paired_examined - std::min(paired_examined, optical_duplicate_pairs),
            paired_duplicates - std::min(paired_duplicates, optical_duplicate_pairs));
        std::ofstream metrics(options.metrics);
        if (!metrics) throw std::runtime_error("cannot write MarkDuplicates metrics: " + options.metrics);
        metrics << "## FASTGATK-MARKDUPLICATES v3\n"
                << "## METRICS CLASS picard.sam.DuplicationMetrics\n"
                << "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n";
        for (const auto& [library_name, library] : library_metrics) {
            const auto examined_records = library.unpaired_reads_examined +
                2 * library.read_pairs_examined;
            const auto optical = std::min(library.read_pairs_examined,
                                           library.read_pair_optical_duplicates);
            const auto estimated = estimate_library_size(
                library.read_pairs_examined - optical,
                library.read_pair_duplicates - std::min(library.read_pair_duplicates, optical));
            const auto percent = examined_records == 0 ? 0.0 :
                static_cast<double>(library.duplicate_records) / examined_records;
            metrics << (library_name.empty() ? "Unknown Library" : library_name) << '\t'
                    << library.unpaired_reads_examined << '\t'
                    << library.read_pairs_examined << '\t'
                    << library.secondary_or_supplementary << '\t'
                    << library.unmapped_reads << '\t'
                    << library.unpaired_read_duplicates << '\t'
                    << library.read_pair_duplicates << '\t'
                    << library.read_pair_optical_duplicates << '\t'
                    << percent << '\t';
            // Picard leaves ESTIMATED_LIBRARY_SIZE blank when the input has
            // no paired duplicate population large enough for its estimator.
            if (estimated != 0) metrics << estimated;
            metrics << "\n";
        }
        if (!duplicate_set_histogram.empty()) {
            metrics << "\n## HISTOGRAM\tjava.lang.Double\n"
                    << "set_size\tall_sets\toptical_sets\tnon_optical_sets\n";
            for (const auto& [set_size, histogram] : duplicate_set_histogram)
                metrics << static_cast<double>(set_size) << '\t'
                        << histogram[0] << '\t' << histogram[1] << '\t'
                        << histogram[2] << '\n';
        }
        metrics.flush();
        if (!metrics || !file_complete(options.metrics)) throw std::runtime_error("cannot finalize MarkDuplicates metrics");
        std::string index_path;
        if (options.create_index && (suffix(options.output, ".bam") || suffix(options.output, ".cram"))) {
            index_path = options.output + (suffix(options.output, ".cram") ? ".crai" : ".bai");
            if (sam_index_build3(options.output.c_str(), index_path.c_str(), 0, 0) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build MarkDuplicates index");
        }
        const bool primary_complete = file_complete(options.output);
        const bool metrics_complete = file_complete(options.metrics);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !metrics_complete || !index_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: MarkDuplicates output/metrics/index incomplete");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write MarkDuplicates manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"MarkDuplicates\",\"implementation\":\"fastgatk-mark-duplicates\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"alignment\","
                 << "\"compatibility\":{\"duplicate_key\":true,\"paired_fragment_consensus\":true,\"library_partition\":true,\"optical_duplicate_detection\":true,\"mark_or_remove\":true,\"tagging_policy\":true,\"optical_only_tagging\":"
                 << ((options.tagging_policy == "OpticalOnly" || options.tagging_policy == "opticalonly") ? "true" : "false")
                 << ",\"metrics\":true,\"bounded_group\":"
                 << (!external_spill ? "true" : "false")
                 << ",\"global_key_pass\":true,\"two_pass_duplicate_marking\":true,\"chunked_group_spill\":true,\"spill_runs\":"
                 << (!spill_runs.empty() ? "true" : "false") << ",\"tmp_dir\":\""
                 << json_escape(options.tmp_dir) << "\",\"max_records_in_memory\":" << options.max_records_in_memory << ",\"bam_index\":"
                 << (index_path.empty() ? "false" : "true") << ",\"per_library_metrics\":true,\"pg_line\":"
                 << (options.add_pg_tag ? "true" : "false")
                 << ",\"remove_sequencing_duplicates\":"
                 << (options.remove_sequencing_duplicates ? "true" : "false")
                 << ",\"spill_checkpoint\":" << (!checkpoint.empty() ? "true" : "false")
                 << ",\"resumed\":" << (resumed ? "true" : "false")
                 << ",\"bit_identical_to_picard\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"alignment\",\"complete\":true},{\"path\":\""
                 << json_escape(options.metrics) << "\",\"kind\":\"duplicate-metrics\",\"complete\":true}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"alignment-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << input_records << ",\"output_records\":" << output_records
                 << ",\"duplicate_records\":" << duplicate_records
                 << ",\"optical_duplicate_records\":" << optical_duplicate_records
                 << ",\"estimated_library_size\":" << estimated_library_size
                 << ",\"library_count\":" << library_metrics.size()
                 << ",\"external_group_record_limit\":" << external_group_record_limit
                 << ",\"spill_run_count\":" << spill_runs.size()
                 << ",\"external_spill\":" << (external_spill ? "true" : "false")
                 << ",\"spill_checkpoint\":" << (!checkpoint.empty() ? "true" : "false")
                 << ",\"resumed\":" << (resumed ? "true" : "false")
                 << ",\"pg_line\":" << (options.add_pg_tag ? "true" : "false") << "}}\n";
        manifest.flush();
        if (!manifest) throw std::runtime_error("cannot finalize MarkDuplicates manifest");
        manifest.close();
        fastgatk::runtime::require_complete_output({
            std::filesystem::path(options.output), std::filesystem::path(index_path),
            std::filesystem::path(manifest_path), !index_path.empty(), true});
        std::cout << "{\"tool\":\"MarkDuplicates\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"output_records\":" << output_records
                  << ",\"duplicate_records\":" << duplicate_records
                  << ",\"optical_duplicate_records\":" << optical_duplicate_records
                  << ",\"estimated_library_size\":" << estimated_library_size
                  << ",\"spill_run_count\":" << spill_runs.size()
                  << ",\"resumed\":" << (resumed ? "true" : "false") << "}\n";
        // A successful output is self-contained; remove restart state after
        // all output contracts have passed.  Checkpoints are retained only
        // when an interrupted second pass reaches the catch block below.
        cleanup_spills(spill_runs);
        if (!checkpoint.empty()) {
            std::error_code checkpoint_error;
            std::filesystem::remove(checkpoint, checkpoint_error);
        }
        sam_hdr_destroy(header);
        sam_hdr_destroy(second_header);
        return 0;
    } catch (...) {
        std::error_code checkpoint_error;
        const bool preserve_checkpoint = !checkpoint.empty() &&
            std::filesystem::is_regular_file(checkpoint, checkpoint_error) && !checkpoint_error;
        // If the first pass completed, keep the checkpoint and all validated
        // runs so a retry can resume the second pass.  Any failure before the
        // checkpoint is atomically published still cleans up its partial runs.
        if (!preserve_checkpoint) cleanup_spills(spill_runs);
        bam_destroy1(record);
        if (output) sam_close(output);
        if (second_input) sam_close(second_input);
        if (input) sam_close(input);
        if (header) sam_hdr_destroy(header);
        if (second_header) sam_hdr_destroy(second_header);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for MarkDuplicates");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-mark-duplicates: " << error.what() << '\n';
        return 2;
    }
}
