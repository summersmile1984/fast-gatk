#include "fastgatk/io/bam_htsjdk.hpp"
#include "fastgatk/runtime/resource.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/sam.h>
#include "fastgatk/io/hts_read_guard.hpp"
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string reference;
    std::string manifest;
    std::string tmp_dir;
    std::string sort_order = "coordinate";
    std::size_t max_records_in_memory = 0;
    int compression_level = 2;
    std::string validation_stringency = "SILENT";
    // Picard SortSam 4.6.2.0 defaults CREATE_INDEX=false.  Keep index
    // creation opt-in so a SAM output or a coordinate BAM has the same
    // default side-effect surface as the direct replacement.
    bool create_index = false;
    bool resume_spill = false;
    bool add_pg_tag = false;
    std::string pg_id = "fastgatk-sort-sam";
    std::string pg_name = "SortSam";
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

std::size_t parse_size(const std::string& value, const char* option) {
    std::size_t consumed = 0;
    const auto parsed = std::stoull(value, &consumed);
    if (consumed != value.size() || parsed == 0) throw std::invalid_argument(std::string("invalid ") + option);
    return static_cast<std::size_t>(parsed);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-sort-sam (GATK/Picard-compatible native prototype)\n"
                         "  -I, --input FILE                 input SAM/BAM/CRAM\n"
                         "      --INPUT FILE                  Picard alias for -I\n"
                         "  -O, --output FILE                sorted SAM/BAM/CRAM\n"
                         "      --OUTPUT FILE                 Picard alias for -O\n"
                         "  -R, --reference FILE             CRAM reference (optional)\n"
                         "      --REFERENCE_SEQUENCE FILE    Picard alias for -R\n"
                         "  --sort-order coordinate|queryname|duplicate (-SO/--SORT_ORDER)\n"
                         "  --max-records-in-memory N (--MAX_RECORDS_IN_RAM)\n"
                         "  --tmp-dir DIR (--TMP_DIR)         spill directory\n"
                         "  --COMPRESSION_LEVEL 0..9          BGZF compression level\n"
                         "  --VALIDATION_STRINGENCY MODE      STRICT/LENIENT/SILENT\n"
                         "  --CREATE_INDEX[=true|false]       create BAM/CRAM index\n"
                         "  --resume-spill                   resume a completed input-pass checkpoint\n"
                         "  --add-pg-tag[=true|false]         add a SAM @PG record\n"
                         "  --program-record-id ID            @PG ID\n"
                         "  --program-group-name NAME         @PG PN\n"
                         "  --program-group-version VERSION   @PG VN\n"
                         "  --program-group-command-line CL   @PG CL\n"
                         "      --output-manifest FILE      OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-I" || is_option(argument, "--input") || is_option(argument, "--INPUT")) {
            const char* name = is_option(argument, "--INPUT") ? "--INPUT" : "--input";
            options.input = require_value(index, argc, argv, argument, name, "-I");
        }
        else if (argument == "-O" || is_option(argument, "--output") || is_option(argument, "--OUTPUT")) {
            const char* name = is_option(argument, "--OUTPUT") ? "--OUTPUT" : "--output";
            options.output = require_value(index, argc, argv, argument, name, "-O");
        }
        else if (argument == "-R" || is_option(argument, "--reference") ||
                 is_option(argument, "--REFERENCE_SEQUENCE")) {
            const char* name = is_option(argument, "--REFERENCE_SEQUENCE")
                ? "--REFERENCE_SEQUENCE" : "--reference";
            options.reference = require_value(index, argc, argv, argument, name, "-R");
        }
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "-SO" || is_option(argument, "--sort-order") || is_option(argument, "--SORT_ORDER")) {
            const char* name = argument == "-SO" ? "-SO" :
                (is_option(argument, "--SORT_ORDER") ? "--SORT_ORDER" : "--sort-order");
            options.sort_order = require_value(index, argc, argv, argument, name, "-SO");
            std::transform(options.sort_order.begin(), options.sort_order.end(), options.sort_order.begin(),
                           [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
        }
        else if (is_option(argument, "--max-records-in-memory") || is_option(argument, "--MAX_RECORDS_IN_RAM")) {
            const char* name = is_option(argument, "--MAX_RECORDS_IN_RAM")
                ? "--MAX_RECORDS_IN_RAM" : "--max-records-in-memory";
            options.max_records_in_memory = parse_size(require_value(index, argc, argv, argument, name), name);
        }
        else if (is_option(argument, "--tmp-dir") || is_option(argument, "--TMP_DIR")) {
            const char* name = is_option(argument, "--TMP_DIR") ? "--TMP_DIR" : "--tmp-dir";
            options.tmp_dir = require_value(index, argc, argv, argument, name);
        }
        else if (is_option(argument, "--COMPRESSION_LEVEL")) {
            const auto value = require_value(index, argc, argv, argument, "--COMPRESSION_LEVEL");
            try { options.compression_level = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid integer for --COMPRESSION_LEVEL: " + value); }
            if (options.compression_level < 0 || options.compression_level > 9)
                throw std::invalid_argument("--COMPRESSION_LEVEL must be between 0 and 9");
        }
        else if (is_option(argument, "--VALIDATION_STRINGENCY")) {
            options.validation_stringency = require_value(index, argc, argv, argument, "--VALIDATION_STRINGENCY");
            std::transform(options.validation_stringency.begin(), options.validation_stringency.end(),
                           options.validation_stringency.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (options.validation_stringency != "STRICT" && options.validation_stringency != "LENIENT" &&
                options.validation_stringency != "SILENT")
                throw std::invalid_argument("--VALIDATION_STRINGENCY must be STRICT, LENIENT, or SILENT");
        }
        else if (argument == "--resume-spill")
            options.resume_spill = true;
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
        else if (argument == "--CREATE_INDEX" || argument.rfind("--CREATE_INDEX=", 0) == 0 ||
                 argument == "--create-output-bam-index" ||
                 argument.rfind("--create-output-bam-index=", 0) == 0 ||
                 argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            const char* name = argument.rfind("--CREATE_INDEX", 0) == 0
                ? "--CREATE_INDEX"
                : (argument.rfind("--create-output-variant-index", 0) == 0
                    ? "--create-output-variant-index" : "--create-output-bam-index");
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument == "--QUIET" ||
                   argument.rfind("--quiet=", 0) == 0 || argument.rfind("--QUIET=", 0) == 0) {
            const char* name = argument.rfind("--QUIET", 0) == 0 ? "--QUIET" : "--quiet";
            (void)fastgatk::native::parse_optional_boolean(index, argc, argv, argument, name);
        } else if (argument == "--disable-sequence-dictionary-validation" ||
                   argument == "--use-jdk-deflater" || argument == "--use-jdk-inflater" ||
                   argument == "--USE_JDK_DEFLATER" || argument == "--USE_JDK_INFLATER") {
            // Accepted compatibility flags. HTSlib is the selected codec backend.
            if (argument != "--disable-sequence-dictionary-validation") {
                (void)fastgatk::native::parse_optional_boolean(index, argc, argv, argument,
                    argument.rfind("--USE_", 0) == 0 ?
                        (argument.rfind("--USE_JDK_DEFLATER", 0) == 0 ? "--USE_JDK_DEFLATER" : "--USE_JDK_INFLATER") :
                        (argument.rfind("--use-jdk-deflater", 0) == 0 ? "--use-jdk-deflater" : "--use-jdk-inflater"));
            }
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
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.resume_spill && options.tmp_dir.empty())
        throw std::invalid_argument("--resume-spill requires --tmp-dir");
    if (options.add_pg_tag && options.pg_id.empty())
        throw std::invalid_argument("--program-record-id must not be empty when --add-pg-tag is enabled");
    if (options.sort_order != "coordinate" && options.sort_order != "queryname" &&
        options.sort_order != "duplicate")
        throw std::invalid_argument("--sort-order must be coordinate, queryname, or duplicate");
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct RecordLess {
    std::string sort_order;
    const std::map<std::string, std::string>* read_group_libraries = nullptr;

    static int compare_int(const int left, const int right) {
        return left < right ? -1 : (left > right ? 1 : 0);
    }

    static int compare_queryname(const bam1_t* left, const bam1_t* right) {
        const int name_compare = std::strcmp(bam_get_qname(left), bam_get_qname(right));
        if (name_compare != 0) return name_compare;
        const bool left_paired = (left->core.flag & BAM_FPAIRED) != 0;
        const bool right_paired = (right->core.flag & BAM_FPAIRED) != 0;
        if (left_paired != right_paired) return left_paired ? -1 : 1;
        if (left_paired) {
            const bool left_first = (left->core.flag & BAM_FREAD1) != 0;
            const bool right_second = (right->core.flag & BAM_FREAD2) != 0;
            if (left_first && right_second) return -1;
            const bool left_second = (left->core.flag & BAM_FREAD2) != 0;
            const bool right_first = (right->core.flag & BAM_FREAD1) != 0;
            if (left_second && right_first) return 1;
        }
        const bool left_reverse = (left->core.flag & BAM_FREVERSE) != 0;
        const bool right_reverse = (right->core.flag & BAM_FREVERSE) != 0;
        if (left_reverse != right_reverse) return left_reverse ? 1 : -1;
        const bool left_secondary = (left->core.flag & BAM_FSECONDARY) != 0;
        const bool right_secondary = (right->core.flag & BAM_FSECONDARY) != 0;
        if (left_secondary != right_secondary) return left_secondary ? 1 : -1;
        const bool left_supplementary = (left->core.flag & BAM_FSUPPLEMENTARY) != 0;
        const bool right_supplementary = (right->core.flag & BAM_FSUPPLEMENTARY) != 0;
        if (left_supplementary != right_supplementary)
            return left_supplementary ? 1 : -1;
        const auto* left_hi = bam_aux_get(left, "HI");
        const auto* right_hi = bam_aux_get(right, "HI");
        if (left_hi == nullptr && right_hi != nullptr) return -1;
        if (left_hi != nullptr && right_hi == nullptr) return 1;
        if (left_hi != nullptr && right_hi != nullptr) {
            const int hi_compare = compare_int(bam_aux2i(left_hi), bam_aux2i(right_hi));
            if (hi_compare != 0) return hi_compare;
        }
        return 0;
    }

    struct CigarShape {
        std::int64_t reference_length = 0;
        std::int64_t leading_clip = 0;
        std::int64_t trailing_clip = 0;
    };

    static CigarShape cigar_shape(const bam1_t* record) {
        CigarShape shape;
        const auto* cigar = bam_get_cigar(record);
        for (std::size_t index = 0; index < record->core.n_cigar; ++index) {
            const auto operation = bam_cigar_op(cigar[index]);
            const auto length = static_cast<std::int64_t>(bam_cigar_oplen(cigar[index]));
            if (bam_cigar_type(operation) & 2) shape.reference_length += length;
        }
        std::size_t index = 0;
        while (index < record->core.n_cigar) {
            const auto operation = bam_cigar_op(cigar[index]);
            if (operation != BAM_CSOFT_CLIP && operation != BAM_CHARD_CLIP) break;
            shape.leading_clip += bam_cigar_oplen(cigar[index++]);
        }
        index = record->core.n_cigar;
        while (index > 0) {
            const auto operation = bam_cigar_op(cigar[index - 1]);
            if (operation != BAM_CSOFT_CLIP && operation != BAM_CHARD_CLIP) break;
            shape.trailing_clip += bam_cigar_oplen(cigar[--index]);
        }
        return shape;
    }

    static CigarShape cigar_shape(const char* text) {
        if (text == nullptr || *text == '\0' || std::strcmp(text, "*") == 0)
            throw std::runtime_error("BAD_INPUT: duplicate SortSam requires a valid MC tag for mapped mates");
        CigarShape shape;
        bool leading_clips = true;
        const char* cursor = text;
        while (*cursor != '\0') {
            char* end = nullptr;
            const auto length = std::strtoll(cursor, &end, 10);
            if (end == cursor || length <= 0 || *end == '\0')
                throw std::runtime_error("BAD_INPUT: malformed MC tag for duplicate SortSam");
            const char operation = *end++;
            if (std::strchr("MIDNSHP=X", operation) == nullptr)
                throw std::runtime_error("BAD_INPUT: malformed MC tag for duplicate SortSam");
            if (operation == 'M' || operation == 'D' || operation == 'N' ||
                operation == '=' || operation == 'X')
                shape.reference_length += length;
            const bool clipping = operation == 'S' || operation == 'H';
            if (leading_clips && clipping)
                shape.leading_clip += length;
            else
                leading_clips = false;
            if (clipping) shape.trailing_clip += length;
            else shape.trailing_clip = 0;
            cursor = end;
        }
        return shape;
    }

    static bool paired_and_both_mapped(const bam1_t* record) {
        return (record->core.flag & BAM_FPAIRED) != 0 &&
               (record->core.flag & BAM_FUNMAP) == 0 &&
               (record->core.flag & BAM_FMUNMAP) == 0;
    }

    static std::int64_t read_coordinate(const bam1_t* record) {
        const auto shape = cigar_shape(record);
        if ((record->core.flag & BAM_FREVERSE) != 0)
            return static_cast<std::int64_t>(record->core.pos) +
                   shape.reference_length - 1 + shape.trailing_clip;
        return static_cast<std::int64_t>(record->core.pos) - shape.leading_clip;
    }

    static CigarShape mate_cigar_shape(const bam1_t* record) {
        const auto* value = bam_aux_get(record, "MC");
        if (value == nullptr || *value != 'Z')
            throw std::runtime_error("BAD_INPUT: duplicate SortSam requires MC tags for mapped mates");
        return cigar_shape(bam_aux2Z(value));
    }

    static std::int64_t mate_coordinate(const bam1_t* record) {
        if (!paired_and_both_mapped(record)) return -1;
        const auto shape = mate_cigar_shape(record);
        if ((record->core.flag & BAM_FMREVERSE) != 0)
            return static_cast<std::int64_t>(record->core.mpos) +
                   shape.reference_length - 1 + shape.trailing_clip;
        return static_cast<std::int64_t>(record->core.mpos) - shape.leading_clip;
    }

    static int orientation(const bam1_t* record) {
        if (!paired_and_both_mapped(record))
            return (record->core.flag & BAM_FREVERSE) != 0 ? 5 : 2;
        const bool read_reverse = (record->core.flag & BAM_FREVERSE) != 0;
        const bool mate_reverse = (record->core.flag & BAM_FMREVERSE) != 0;
        if (read_reverse) return mate_reverse ? 4 : 3;
        return mate_reverse ? 1 : 0;
    }

    std::string library(const bam1_t* record) const {
        const auto* value = bam_aux_get(record, "RG");
        if (value != nullptr && *value == 'Z' && read_group_libraries != nullptr) {
            const auto* rg = bam_aux2Z(value);
            const auto found = rg == nullptr ? read_group_libraries->end()
                                             : read_group_libraries->find(rg);
            if (found != read_group_libraries->end() && !found->second.empty())
                return found->second;
        }
        return "Unknown Library";
    }

    static bool has_mapped_end(const bam1_t* record) {
        return (record->core.flag & BAM_FUNMAP) == 0 ||
               ((record->core.flag & BAM_FPAIRED) != 0 &&
                (record->core.flag & BAM_FMUNMAP) == 0);
    }

    static bool has_unmapped_end(const bam1_t* record) {
        return (record->core.flag & BAM_FUNMAP) != 0 ||
               ((record->core.flag & BAM_FPAIRED) != 0 &&
                (record->core.flag & BAM_FMUNMAP) != 0);
    }

    static std::int16_t duplicate_score(const bam1_t* record) {
        std::int32_t score = 0;
        if ((record->core.flag & BAM_FUNMAP) == 0)
            score = static_cast<std::int32_t>(std::min<std::int64_t>(
                cigar_shape(record).reference_length, 16383));
        if (paired_and_both_mapped(record)) {
            score += static_cast<std::int32_t>(std::min<std::int64_t>(
                mate_cigar_shape(record).reference_length, 16383));
        }
        auto result = static_cast<std::int16_t>(score);
        if ((record->core.flag & BAM_FQCFAIL) != 0)
            result = static_cast<std::int16_t>(static_cast<std::int32_t>(result) - 16384);
        return result;
    }

    static std::string canonical_name(const bam1_t* record) {
        const char* name = bam_get_qname(record);
        const auto* value = bam_aux_get(record, "RG");
        if (value == nullptr || *value != 'Z') return name == nullptr ? std::string{} : name;
        const auto* rg = bam_aux2Z(value);
        return std::string(rg == nullptr ? "" : rg) + ":" +
               (name == nullptr ? std::string{} : std::string(name));
    }

    int compare_duplicate(const bam1_t* left, const bam1_t* right) const {
        const auto left_library = library(left);
        const auto right_library = library(right);
        if (left_library != right_library) return left_library < right_library ? -1 : 1;
        const int left_tid = left->core.tid;
        const int right_tid = right->core.tid;
        if (left_tid != right_tid) {
            if (left_tid < 0) return 1;
            if (right_tid < 0) return -1;
            return left_tid < right_tid ? -1 : 1;
        }
        const auto left_coordinate = read_coordinate(left);
        const auto right_coordinate = read_coordinate(right);
        if (left_coordinate != right_coordinate) return left_coordinate < right_coordinate ? -1 : 1;
        const int left_orientation = orientation(left);
        const int right_orientation = orientation(right);
        if (left_orientation != right_orientation) return left_orientation < right_orientation ? -1 : 1;
        if (paired_and_both_mapped(left) && paired_and_both_mapped(right)) {
            if (left->core.mtid != right->core.mtid)
                return left->core.mtid < right->core.mtid ? -1 : 1;
            const auto left_mate_coordinate = mate_coordinate(left);
            const auto right_mate_coordinate = mate_coordinate(right);
            if (left_mate_coordinate != right_mate_coordinate)
                return left_mate_coordinate < right_mate_coordinate ? -1 : 1;
        }
        const bool left_has_mapped = has_mapped_end(left);
        const bool right_has_mapped = has_mapped_end(right);
        if (left_has_mapped != right_has_mapped) return left_has_mapped ? -1 : 1;
        const bool left_paired = (left->core.flag & BAM_FPAIRED) != 0;
        const bool right_paired = (right->core.flag & BAM_FPAIRED) != 0;
        if (left_paired == right_paired) {
            const bool left_has_unmapped = has_unmapped_end(left);
            const bool right_has_unmapped = has_unmapped_end(right);
            if (left_has_unmapped != right_has_unmapped) return left_has_unmapped ? 1 : -1;
        } else {
            return left_paired ? -1 : 1;
        }
        if (left_paired != right_paired) return left_paired ? -1 : 1;
        const auto left_score = duplicate_score(left);
        const auto right_score = duplicate_score(right);
        if (left_score != right_score) return left_score > right_score ? -1 : 1;
        const auto left_canonical = canonical_name(left);
        const auto right_canonical = canonical_name(right);
        if (left_canonical != right_canonical) return left_canonical < right_canonical ? -1 : 1;
        const int name_compare = std::strcmp(bam_get_qname(left), bam_get_qname(right));
        if (name_compare != 0) return name_compare;
        if (left_paired && right_paired) {
            const int left_end = (left->core.flag & BAM_FREAD1) != 0 ? 0 : 1;
            const int right_end = (right->core.flag & BAM_FREAD1) != 0 ? 0 : 1;
            if (left_end != right_end) return left_end < right_end ? -1 : 1;
        }
        return 0;
    }

    bool operator()(const bam1_t* left, const bam1_t* right) const {
        if (sort_order == "queryname") {
            return compare_queryname(left, right) < 0;
        } else if (sort_order == "duplicate") {
            return compare_duplicate(left, right) < 0;
        } else {
            const bool left_unmapped = left->core.tid < 0 || left->core.pos < 0;
            const bool right_unmapped = right->core.tid < 0 || right->core.pos < 0;
            if (left_unmapped != right_unmapped) return !left_unmapped;
            if (!left_unmapped) {
                if (left->core.tid != right->core.tid) return left->core.tid < right->core.tid;
                if (left->core.pos != right->core.pos) return left->core.pos < right->core.pos;
            }
            const bool left_reverse = (left->core.flag & BAM_FREVERSE) != 0;
            const bool right_reverse = (right->core.flag & BAM_FREVERSE) != 0;
            if (left_reverse != right_reverse) return !left_reverse;
            const int name_compare = std::strcmp(bam_get_qname(left), bam_get_qname(right));
            if (name_compare != 0) return name_compare < 0;
            if (left->core.flag != right->core.flag) return left->core.flag < right->core.flag;
            if (left->core.qual != right->core.qual) return left->core.qual < right->core.qual;
            if (left->core.mtid != right->core.mtid) return left->core.mtid < right->core.mtid;
            if (left->core.mpos != right->core.mpos) return left->core.mpos < right->core.mpos;
            return left->core.isize < right->core.isize;
        }
    }
};

std::size_t record_bytes(const bam1_t* record) {
    return sizeof(bam1_t) + static_cast<std::size_t>(record->l_data) + 64;
}

void destroy_records(std::vector<bam1_t*>& records) {
    for (auto* record : records) bam_destroy1(record);
    records.clear();
}

std::filesystem::path spill_path(const std::filesystem::path& directory, std::size_t sequence) {
    return directory / ("fastgatk-sort-run-" + std::to_string(sequence) + ".bam");
}

struct SpillCheckpoint {
    std::string input;
    std::string output;
    std::string reference;
    std::string sort_order;
    std::uint64_t input_size = 0;
    std::int64_t input_mtime = 0;
    std::size_t record_limit = 0;
    std::uint64_t input_records = 0;
    std::uint64_t spilled_bytes = 0;
    std::vector<std::filesystem::path> runs;
};

std::filesystem::path checkpoint_path(const Options& options) {
    // A per-input/output name prevents two concurrent Nextflow tasks sharing
    // one scratch root from consuming each other's checkpoint.  The hash is
    // only a filename partition; all semantic fields are still validated on
    // read before a checkpoint is trusted.
    const auto identity = options.input + "\n" + options.output + "\n" + options.sort_order;
    const auto hash = std::hash<std::string>{}(identity);
    std::ostringstream name;
    name << "fastgatk-sort-checkpoint-" << std::hex << hash << ".txt";
    return std::filesystem::path(options.tmp_dir) / name.str();
}

void write_checkpoint(const std::filesystem::path& path, const Options& options,
                      std::size_t record_limit, std::uint64_t input_records,
                      std::uint64_t spilled_bytes,
                      const std::vector<std::filesystem::path>& runs) {
    const auto temporary = path.string() + ".tmp";
    std::error_code stat_error;
    const auto input_size = std::filesystem::file_size(options.input, stat_error);
    if (stat_error) throw std::runtime_error("BAD_INPUT: cannot stat SortSam input for checkpoint: " + options.input);
    const auto input_mtime = std::filesystem::last_write_time(options.input, stat_error)
        .time_since_epoch().count();
    if (stat_error) throw std::runtime_error("BAD_INPUT: cannot timestamp SortSam input for checkpoint: " + options.input);
    std::ofstream stream(temporary, std::ios::trunc);
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot write SortSam checkpoint");
    stream << "FASTGATK_SORT_CHECKPOINT_V1\n"
           << "input " << std::quoted(options.input) << "\n"
           << "output " << std::quoted(options.output) << "\n"
           << "reference " << std::quoted(options.reference) << "\n"
           << "sort_order " << std::quoted(options.sort_order) << "\n"
           << "input_size " << input_size << "\n"
           << "input_mtime " << input_mtime << "\n"
           << "record_limit " << record_limit << "\n"
           << "input_records " << input_records << "\n"
           << "spilled_bytes " << spilled_bytes << "\n"
           << "runs " << runs.size() << "\n";
    for (const auto& run : runs) stream << "run " << std::quoted(run.string()) << "\n";
    stream.flush();
    if (!stream) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot finalize SortSam checkpoint");
    std::error_code error;
    std::filesystem::rename(temporary, path, error);
    if (error) {
        std::filesystem::remove(path, error);
        error.clear();
        std::filesystem::rename(temporary, path, error);
    }
    if (error) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot install SortSam checkpoint: " + path.string());
}

SpillCheckpoint read_checkpoint(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("BAD_INPUT: cannot open SortSam checkpoint: " + path.string());
    std::string version;
    if (!(stream >> version) || version != "FASTGATK_SORT_CHECKPOINT_V1")
        throw std::runtime_error("BAD_INPUT: unsupported SortSam checkpoint version");
    SpillCheckpoint checkpoint;
    std::string key;
    while (stream >> key) {
        if (key == "input") stream >> std::quoted(checkpoint.input);
        else if (key == "output") stream >> std::quoted(checkpoint.output);
        else if (key == "reference") stream >> std::quoted(checkpoint.reference);
        else if (key == "sort_order") stream >> std::quoted(checkpoint.sort_order);
        else if (key == "input_size") stream >> checkpoint.input_size;
        else if (key == "input_mtime") stream >> checkpoint.input_mtime;
        else if (key == "record_limit") stream >> checkpoint.record_limit;
        else if (key == "input_records") stream >> checkpoint.input_records;
        else if (key == "spilled_bytes") stream >> checkpoint.spilled_bytes;
        else if (key == "runs") {
            std::size_t count = 0;
            stream >> count;
            if (count > (1ULL << 20)) throw std::runtime_error("BAD_INPUT: invalid SortSam checkpoint run count");
            checkpoint.runs.reserve(count);
        } else if (key == "run") {
            std::string value;
            stream >> std::quoted(value);
            checkpoint.runs.emplace_back(value);
        } else {
            throw std::runtime_error("BAD_INPUT: unknown SortSam checkpoint field: " + key);
        }
        if (!stream) throw std::runtime_error("BAD_INPUT: truncated SortSam checkpoint");
    }
    if (checkpoint.input.empty() || checkpoint.output.empty() ||
        (checkpoint.sort_order != "coordinate" && checkpoint.sort_order != "queryname" &&
         checkpoint.sort_order != "duplicate") ||
        checkpoint.record_limit == 0)
        throw std::runtime_error("BAD_INPUT: incomplete SortSam checkpoint");
    if (checkpoint.input_size == 0)
        throw std::runtime_error("BAD_INPUT: SortSam checkpoint has no input size");
    if (checkpoint.runs.empty())
        throw std::runtime_error("BAD_INPUT: SortSam checkpoint contains no spill runs");
    for (const auto& run : checkpoint.runs) {
        if (run.empty()) continue;
        if (!file_complete(run.string()))
            throw std::runtime_error("BAD_INPUT: missing or empty SortSam spill run: " + run.string());
    }
    return checkpoint;
}

void write_run(const std::filesystem::path& path, const sam_hdr_t* header,
              std::vector<bam1_t*>& records, const RecordLess& less) {
    std::sort(records.begin(), records.end(), less);
    samFile* output = sam_open(path.c_str(), "wb");
    if (!output) throw std::runtime_error("cannot open spill run: " + path.string());
    if (sam_hdr_write(output, header) < 0) {
        sam_close(output);
        throw std::runtime_error("cannot write spill run header: " + path.string());
    }
    for (const auto* record : records) {
        if (sam_write1(output, header, record) < 0) {
            sam_close(output);
            throw std::runtime_error("cannot write spill run: " + path.string());
        }
    }
    if (sam_close(output) != 0) throw std::runtime_error("cannot close spill run: " + path.string());
}

struct MergeRun {
    samFile* file = nullptr;
    sam_hdr_t* header = nullptr;
    bam1_t* record = nullptr;
};

struct HeapNode {
    std::size_t run = 0;
    bam1_t* record = nullptr;
};

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    samFile* input = nullptr;
    samFile* output = nullptr;
    sam_hdr_t* header = nullptr;
    bam1_t* record = nullptr;
    std::vector<bam1_t*> records;
    std::vector<std::filesystem::path> runs;
    std::vector<MergeRun> merge_runs;
    std::uint64_t input_records = 0;
    std::uint64_t output_records = 0;
    std::uint64_t staged_bytes = 0;
    std::uint64_t spilled_bytes = 0;
    std::size_t run_sequence = 0;
    std::filesystem::path checkpoint;
    bool resumed = false;
    const bool coordinate = options.sort_order == "coordinate";
    std::map<std::string, std::string> read_group_libraries;
    const RecordLess less{options.sort_order, &read_group_libraries};
    try {
        input = sam_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open SAM/BAM/CRAM: " + options.input);
        if (!options.reference.empty() && hts_set_fai_filename(input, options.reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure input reference");
        header = sam_hdr_read(input);
        if (!header) throw std::runtime_error("BAD_INPUT: cannot read SAM/BAM header");
        if (options.sort_order == "duplicate") {
            const char* text = sam_hdr_str(header);
            std::istringstream lines(text == nullptr ? "" : text);
            std::string line;
            while (std::getline(lines, line)) {
                if (line.rfind("@RG\t", 0) != 0) continue;
                std::string id, library;
                std::istringstream fields(line.substr(4));
                std::string field;
                while (std::getline(fields, field, '\t')) {
                    if (field.rfind("ID:", 0) == 0) id = field.substr(3);
                    else if (field.rfind("LB:", 0) == 0) library = field.substr(3);
                }
                if (!id.empty()) read_group_libraries[id] = library;
            }
        }
        if (sam_hdr_update_line(header, "HD", NULL, NULL, "SO", options.sort_order.c_str(), NULL) < 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update sort-order header");
        if (options.add_pg_tag) {
            const auto command_line = options.pg_command_line.empty()
                ? std::string("fastgatk SortSam") : options.pg_command_line;
            if (sam_hdr_add_pg(header, options.pg_name.c_str(),
                               "ID", options.pg_id.c_str(),
                               "VN", options.pg_version.c_str(),
                               "CL", command_line.c_str(), NULL) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot add SortSam @PG header");
        }

        const auto snapshot_budget = resources.safe_memory_budget_bytes();
        const std::size_t default_limit = snapshot_budget == 0
            ? 100000 : std::max<std::size_t>(1, static_cast<std::size_t>(snapshot_budget / 4096));
        const std::size_t record_limit = options.max_records_in_memory == 0
            ? default_limit : options.max_records_in_memory;
        std::filesystem::path spill_directory;
        if (!options.tmp_dir.empty()) spill_directory = options.tmp_dir;
        else if (!resources.scratch_directory.empty()) spill_directory = resources.scratch_directory;
        else spill_directory = std::filesystem::temp_directory_path();
        std::error_code directory_error;
        std::filesystem::create_directories(spill_directory, directory_error);
        if (directory_error || !std::filesystem::is_directory(spill_directory))
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot create sort tmp directory: " + spill_directory.string());

        checkpoint = options.resume_spill ? checkpoint_path(options) : std::filesystem::path{};
        std::error_code checkpoint_error;
        const bool resume_requested = options.resume_spill &&
            std::filesystem::is_regular_file(checkpoint, checkpoint_error) && !checkpoint_error;
        if (resume_requested) {
            const auto saved = read_checkpoint(checkpoint);
            if (saved.input != options.input || saved.output != options.output ||
                saved.reference != options.reference || saved.sort_order != options.sort_order ||
                saved.record_limit != record_limit)
                throw std::runtime_error("BAD_INPUT: SortSam checkpoint does not match input/output/reference/sort/budget");
            std::error_code stat_error;
            const auto current_size = std::filesystem::file_size(options.input, stat_error);
            const auto current_mtime = stat_error ? std::int64_t{0} :
                std::filesystem::last_write_time(options.input, stat_error).time_since_epoch().count();
            if (stat_error || current_size != saved.input_size || current_mtime != saved.input_mtime)
                throw std::runtime_error("BAD_INPUT: SortSam input changed since spill checkpoint");
            runs = saved.runs;
            run_sequence = runs.size();
            input_records = saved.input_records;
            spilled_bytes = saved.spilled_bytes;
            resumed = true;
            if (sam_close(input) != 0) throw std::runtime_error("cannot close input while resuming SortSam");
            input = nullptr;
        } else {
            record = bam_init1();
            if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed");
            while (fastgatk::io::read_alignment_record(input, header, record,
                                                     options.input) >= 0) {
                ++input_records;
                auto* copy = bam_dup1(record);
                if (!copy) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot stage BAM record");
                records.push_back(copy);
                staged_bytes += record_bytes(copy);
                if (records.size() >= record_limit ||
                    (snapshot_budget != 0 && staged_bytes * 2 >= snapshot_budget)) {
                    const auto path = spill_path(spill_directory, run_sequence++);
                    write_run(path, header, records, less);
                    spilled_bytes += staged_bytes;
                    runs.push_back(path);
                    destroy_records(records);
                    staged_bytes = 0;
                }
            }
            bam_destroy1(record); record = nullptr;
            if (!records.empty() || runs.empty()) {
                if (!records.empty()) {
                    const auto path = spill_path(spill_directory, run_sequence++);
                    write_run(path, header, records, less);
                    spilled_bytes += staged_bytes;
                    runs.push_back(path);
                    destroy_records(records);
                    staged_bytes = 0;
                } else {
                    // Empty input still produces a valid header-only output.
                    runs.emplace_back();
                }
            }
            if (sam_close(input) != 0) throw std::runtime_error("cannot close input");
            input = nullptr;
            if (options.resume_spill) {
                write_checkpoint(checkpoint, options, record_limit, input_records,
                                 spilled_bytes, runs);
            }
        }

        const char* output_mode = suffix(options.output, ".cram") ? "wc" :
            suffix(options.output, ".sam") ? "w" : "wb";
        output = sam_open(options.output.c_str(), output_mode);
        if (!output) throw std::runtime_error("cannot open SortSam output: " + options.output);
        if (!suffix(options.output, ".sam") &&
            hts_set_opt(output, HTS_OPT_COMPRESSION_LEVEL, options.compression_level) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set SortSam compression level");
        if (!options.reference.empty() && hts_set_fai_filename(output, options.reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure output reference");
        // htsjdk SAMTextHeaderCodec write semantics: @HD VN 1.6, SO reflects
        // the sorted output order, DT/PT timestamps re-emitted in the
        // process time zone.  ("duplicate" is a native-only order with no
        // SAM SO code; keep the input declaration for it.)
        const char* sort_code = (options.sort_order == "coordinate" ||
                                 options.sort_order == "queryname")
                                    ? options.sort_order.c_str() : nullptr;
        sam_hdr_t* output_header = fastgatk::io::htsjdk_header(header, sort_code);
        if (output_header == nullptr)
            throw std::runtime_error("INTERNAL_ERROR: cannot normalize SortSam header");
        const int header_status = sam_hdr_write(output, output_header);
        sam_hdr_destroy(output_header);
        if (header_status < 0) throw std::runtime_error("cannot write SortSam header");

        auto write_record = [&](bam1_t* value) {
            // Records flow through SortSam untouched; htsjdk then writes
            // binary attributes byte-identically (no value-based re-encode).
            if (sam_write1(output, header, value) < 0)
                throw std::runtime_error("cannot write SortSam record");
            ++output_records;
        };
        if (runs.size() == 1 && runs.front().empty()) {
            // header-only input
        } else {
            for (const auto& path : runs) {
                MergeRun state;
                state.file = sam_open(path.c_str(), "rb");
                if (!state.file) throw std::runtime_error("cannot reopen spill run: " + path.string());
                state.header = sam_hdr_read(state.file);
                state.record = bam_init1();
                if (!state.header || !state.record || sam_read1(state.file, state.header, state.record) < 0) {
                    if (state.record) bam_destroy1(state.record);
                    if (state.header) sam_hdr_destroy(state.header);
                    sam_close(state.file);
                    throw std::runtime_error("cannot initialize spill merge: " + path.string());
                }
                merge_runs.push_back(state);
            }
            auto heap_less = [&](const HeapNode& left, const HeapNode& right) {
                return less(right.record, left.record);
            };
            std::priority_queue<HeapNode, std::vector<HeapNode>, decltype(heap_less)> heap(heap_less);
            for (std::size_t index = 0; index < merge_runs.size(); ++index)
                heap.push({index, merge_runs[index].record});
            while (!heap.empty()) {
                const auto node = heap.top();
                heap.pop();
                write_record(node.record);
                auto& state = merge_runs[node.run];
                const int status = fastgatk::io::read_alignment_record(
                    state.file, state.header, state.record, runs[node.run].string());
                if (status >= 0) heap.push({node.run, state.record});
            }
        }
        if (sam_close(output) != 0) throw std::runtime_error("cannot close SortSam output");
        output = nullptr;
        for (auto& state : merge_runs) {
            bam_destroy1(state.record);
            sam_hdr_destroy(state.header);
            sam_close(state.file);
        }
        merge_runs.clear();
        for (const auto& path : runs) if (!path.empty()) std::filesystem::remove(path, directory_error);
        runs.clear();
        if (!checkpoint.empty()) std::filesystem::remove(checkpoint, directory_error);

        std::string index_path;
        if (options.create_index && coordinate && (suffix(options.output, ".bam") || suffix(options.output, ".cram"))) {
            index_path = options.output + (suffix(options.output, ".cram") ? ".crai" : ".bai");
            if (sam_index_build3(options.output.c_str(), index_path.c_str(), 0, 0) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build SortSam index");
        }
        const bool primary_complete = file_complete(options.output);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !index_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: SortSam output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write SortSam manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"SortSam\",\"implementation\":\"fastgatk-sort-sam\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"alignment\","
                 << "\"compatibility\":{\"external_memory_sort\":true,\"spill_runs\":true,\"sort_order\":\""
                 << json_escape(options.sort_order) << "\",\"bam_index\":" << (index_path.empty() ? "false" : "true")
                 << ",\"spill_checkpoint\":true,\"resumed\":" << (resumed ? "true" : "false")
                 << ",\"pg_line\":" << (options.add_pg_tag ? "true" : "false")
                 << ",\"compression_level\":" << options.compression_level
                 << ",\"validation_stringency\":\"" << options.validation_stringency << "\""
                 << ",\"bit_identical_to_picard\":false},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"alignment\",\"complete\":true}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"alignment-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << input_records << ",\"output_records\":" << output_records
                 << ",\"spill_runs\":" << run_sequence << ",\"spilled_bytes\":" << spilled_bytes
                 << ",\"record_limit\":" << record_limit
                 << ",\"spill_checkpoint\":" << (!checkpoint.empty() ? "true" : "false")
                  << ",\"resumed\":" << (resumed ? "true" : "false")
                  << ",\"pg_line\":" << (options.add_pg_tag ? "true" : "false")
                  << ",\"compression_level\":" << options.compression_level
                  << ",\"validation_stringency\":\"" << options.validation_stringency << "\"}}\n";
        if (!manifest) throw std::runtime_error("cannot finalize SortSam manifest");
        std::cout << "{\"tool\":\"SortSam\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"output_records\":" << output_records
                  << ",\"spill_runs\":" << run_sequence
                  << ",\"resumed\":" << (resumed ? "true" : "false") << "}\n";
        sam_hdr_destroy(header);
        return 0;
    } catch (...) {
        bam_destroy1(record);
        destroy_records(records);
        if (output) sam_close(output);
        if (input) sam_close(input);
        for (auto& state : merge_runs) {
            bam_destroy1(state.record);
            sam_hdr_destroy(state.header);
            sam_close(state.file);
        }
        std::error_code checkpoint_error;
        const bool preserve_checkpoint = !checkpoint.empty() &&
            std::filesystem::is_regular_file(checkpoint, checkpoint_error) && !checkpoint_error;
        if (!preserve_checkpoint) {
            for (const auto& path : runs) {
                std::error_code ignored;
                if (!path.empty()) std::filesystem::remove(path, ignored);
            }
        }
        if (header) sam_hdr_destroy(header);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for SortSam");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-sort-sam: " << error.what() << '\n';
        return 2;
    }
}
