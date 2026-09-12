#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include "fastgatk/io/hts_read_guard.hpp"
#include <zlib.h>

#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/runtime/pipeline.hpp"
#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::vector<std::string> inputs;
    std::vector<std::string> intervals;
    std::vector<std::string> exclusions;
    std::string reference;
    std::string output;
    std::string manifest;
    int min_median_mapq = 50;
    int min_base_quality = 20;
    int max_depth = 200;
    int threads = 1;
    std::size_t batch_records = 4096;
};

struct Interval {
    std::string contig;
    int tid = -1;
    std::int64_t begin = 0;
    std::int64_t end = 0;
};

struct Observation {
    std::int32_t locus = -1;
    std::int32_t base = -1;
    std::int32_t mapq = 0;
    std::uint8_t f1r2 = 0;
    std::uint8_t indel = 0;
};

// A decoded item is deliberately a bounded observation slice.  The Host
// reader still owns BAM/CRAM parsing and locus identity; the Kokkos stage only
// consumes flat numeric fields.  Keeping the global locus id in the slice
// avoids string/object access in the execution space.
struct F1R2Batch {
    std::vector<Observation> observations;
};

struct F1R2BatchResult {
    std::vector<std::int32_t> loci;
    std::vector<std::uint64_t> counts;
    std::vector<std::uint64_t> f1r2_counts;
    std::vector<std::uint64_t> observation_counts;
    std::vector<std::uint64_t> indel_counts;
    std::size_t kernel_batches = 0;
    std::size_t kernel_observations = 0;
    double kernel_prepare_seconds = 0.0;
    double kernel_execute_seconds = 0.0;
};

struct Locus {
    std::string sample;
    std::string contig;
    int tid = -1;
    std::int64_t position = -1;
};

using LocusKey = std::tuple<std::string, std::string, std::int64_t>;

struct SampleOutput {
    std::string sample;
    std::array<std::vector<std::uint64_t>, 64> ref_hist{};
    std::array<std::vector<std::uint64_t>, 384> alt_hist{};
    struct AltRow {
        std::string context;
        int ref_count = 0;
        int alt_count = 0;
        int ref_f1r2 = 0;
        int alt_f1r2 = 0;
        int depth = 0;
        char alt = 'N';
    };
    std::vector<AltRow> alt_rows;
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
    if ((argument == name || (short_name != nullptr && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-collect-f1r2-counts (GATK-compatible native path)\n"
                         "  -R, --reference FILE             reference FASTA (+ .fai)\n"
                         "  -I/--I, --input FILE              input BAM/CRAM (repeatable)\n"
                         "  -L, --intervals REGION           interval (repeatable)\n"
                         "  -XL, --exclude-intervals REGION  subtract interval (repeatable)\n"
                         "  -O, --output FILE                output .tar.gz\n"
                         "      --f1r2-median-mq INT         minimum median mapping quality\n"
                         "      --f1r2-min-bq INT            minimum base quality (strict >)\n"
                         "      --f1r2-max-depth INT         histogram depth cap\n"
                         "      --threads INT                 Kokkos execution threads\n"
                         "      --batch-records INT           bounded observation batch size\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-I" || has_option(argument, "--I") || has_option(argument, "--input")) {
            const char* name = argument.rfind("--input", 0) == 0 ? "--input" : "--I";
            options.inputs.push_back(require_value(index, argc, argv, argument, name, "-I"));
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                (argument.rfind("--interval", 0) == 0 ? "--interval" : "--region");
            options.intervals.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals") ||
                   has_option(argument, "--exclude-interval") || has_option(argument, "--exclude-region")) {
            const char* name = argument == "-XL" ? "--exclude-intervals" :
                (argument.rfind("--exclude-interval", 0) == 0 ? "--exclude-interval" :
                 argument.rfind("--exclude-region", 0) == 0 ? "--exclude-region" : "--exclude-intervals");
            options.exclusions.push_back(require_value(index, argc, argv, argument, name, "-XL"));
        } else if (argument == "-O" || has_option(argument, "--O") || has_option(argument, "--output")) {
            const char* name = argument.rfind("--output", 0) == 0 ? "--output" : "--O";
            options.output = require_value(index, argc, argv, argument, name, "-O");
        } else if (has_option(argument, "--f1r2-median-mq")) {
            options.min_median_mapq = std::stoi(require_value(index, argc, argv, argument, "--f1r2-median-mq"));
        } else if (has_option(argument, "--f1r2-min-bq")) {
            options.min_base_quality = std::stoi(require_value(index, argc, argv, argument, "--f1r2-min-bq"));
        } else if (has_option(argument, "--f1r2-max-depth")) {
            options.max_depth = std::stoi(require_value(index, argc, argv, argument, "--f1r2-max-depth"));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--batch-records")) {
            const auto value = require_value(index, argc, argv, argument, "--batch-records");
            try {
                options.batch_records = static_cast<std::size_t>(std::stoull(value));
            } catch (const std::exception&) {
                throw std::invalid_argument("invalid --batch-records: " + value);
            }
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Accepted GATK launcher/tool flags.
        } else if (has_option(argument, "--java-options") || has_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.inputs.empty()) throw std::invalid_argument("at least one -I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.output.size() < 7 || options.output.substr(options.output.size() - 7) != ".tar.gz")
        throw std::invalid_argument("-O/--output must end in .tar.gz");
    if (options.min_median_mapq < 0 || options.min_base_quality < 0 ||
        options.max_depth < 1 || options.threads < 1 || options.batch_records < 1)
        throw std::invalid_argument("F1R2 parameters must be non-negative and max-depth/threads positive");
    return options;
}

int base_index(char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

char base_char(int index) {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    return index >= 0 && index < 4 ? bases[static_cast<std::size_t>(index)] : 'N';
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char c : value) {
        if (c == '"' || c == '\\') out << '\\';
        if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else out << c;
    }
    return out.str();
}

std::string url_encode(const std::string& value) {
    std::ostringstream out;
    out << std::uppercase << std::hex;
    for (const unsigned char c : value) {
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
            (c >= '0' && c <= '9') || c == '-' || c == '_' || c == '.' || c == '~')
            out << static_cast<char>(c);
        else out << '%' << std::setw(2) << std::setfill('0') << static_cast<int>(c)
                 << std::setfill(' ');
    }
    return out.str();
}

std::vector<Interval> parse_intervals(const std::vector<std::string>& raw,
                                      const bam_hdr_t* header) {
    std::vector<Interval> result;
    for (const auto& value : raw) {
        const auto colon = value.find(':');
        const std::string contig = colon == std::string::npos ? value : value.substr(0, colon);
        int tid = -1;
        for (int i = 0; i < header->n_targets; ++i)
            if (contig == header->target_name[i]) { tid = i; break; }
        if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in BAM header: " + contig);
        std::int64_t begin = 0;
        std::int64_t end = header->target_len[tid];
        if (colon != std::string::npos) {
            const auto coordinates = value.substr(colon + 1);
            const auto dash = coordinates.find('-');
            const auto start_text = coordinates.substr(0, dash);
            const auto end_text = dash == std::string::npos ? start_text : coordinates.substr(dash + 1);
            try {
                begin = std::stoll(start_text) - 1;
                end = std::stoll(end_text);
            } catch (const std::exception&) {
                throw std::runtime_error("BAD_INPUT: malformed interval: " + value);
            }
        }
        if (begin < 0 || end <= begin || end > header->target_len[tid])
            throw std::runtime_error("BAD_INPUT: interval outside contig: " + value);
        result.push_back(Interval{contig, tid, begin, end});
    }
    std::sort(result.begin(), result.end(), [](const Interval& left, const Interval& right) {
        return std::tie(left.tid, left.begin, left.end) < std::tie(right.tid, right.begin, right.end);
    });
    std::vector<Interval> merged;
    for (const auto& interval : result) {
        if (!merged.empty() && merged.back().tid == interval.tid &&
            interval.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, interval.end);
        else merged.push_back(interval);
    }
    return merged;
}

std::vector<Interval> subtract_intervals(const std::vector<Interval>& input,
                                         const std::vector<Interval>& exclusions) {
    if (exclusions.empty()) return input;
    std::vector<Interval> result;
    for (const auto& source : input) {
        std::int64_t cursor = source.begin;
        for (const auto& exclusion : exclusions) {
            if (exclusion.tid != source.tid) continue;
            if (exclusion.end <= cursor) continue;
            if (exclusion.begin >= source.end) break;
            if (exclusion.begin > cursor)
                result.push_back(Interval{source.contig, source.tid, cursor,
                                          std::min(source.end, exclusion.begin)});
            if (exclusion.end >= source.end) { cursor = source.end; break; }
            cursor = std::max(cursor, exclusion.end);
        }
        if (cursor < source.end)
            result.push_back(Interval{source.contig, source.tid, cursor, source.end});
    }
    return result;
}

bool selected(const std::vector<Interval>& intervals, const bam_hdr_t* header,
              int tid, std::int64_t position) {
    if (intervals.empty()) return true;
    const auto* contig = tid >= 0 && tid < header->n_targets ? header->target_name[tid] : nullptr;
    if (contig == nullptr) return false;
    return std::any_of(intervals.begin(), intervals.end(), [&](const Interval& interval) {
        return interval.contig == contig && position >= interval.begin && position < interval.end;
    });
}

std::unordered_map<std::string, std::string> read_group_samples(bam_hdr_t* header) {
    std::unordered_map<std::string, std::string> result;
    const char* raw = sam_hdr_str(header);
    if (raw == nullptr) return result;
    std::istringstream lines(raw);
    for (std::string line; std::getline(lines, line);) {
        if (line.rfind("@RG\t", 0) != 0) continue;
        std::string id;
        std::string sample;
        std::istringstream fields(line);
        for (std::string field; std::getline(fields, field, '\t');) {
            if (field.rfind("ID:", 0) == 0) id = field.substr(3);
            else if (field.rfind("SM:", 0) == 0) sample = field.substr(3);
        }
        if (!id.empty() && !sample.empty()) result[id] = sample;
    }
    return result;
}

bool good_cigar(const bam1_t* record) {
    if (record->core.n_cigar == 0) return false;
    std::int64_t reference_span = 0;
    const auto* cigar = bam_get_cigar(record);
    for (std::uint32_t index = 0; index < record->core.n_cigar; ++index) {
        const auto op = bam_cigar_op(cigar[index]);
        const auto length = bam_cigar_oplen(cigar[index]);
        if (length == 0 || op > BAM_CPAD) return false;
        if (bam_cigar_type(op) & 2) reference_span += length;
    }
    return reference_span > 0;
}

bool keep_read(const bam1_t* record) {
    constexpr std::uint16_t rejected_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY |
                                             BAM_FDUP | BAM_FQCFAIL;
    if ((record->core.flag & rejected_flags) != 0) return false;
    if (record->core.qual < 20 || record->core.qual == 255) return false;
    if (record->core.l_qseq < 30) return false;
    return good_cigar(record);
}

void append_observation(std::vector<Locus>& loci,
                        std::map<LocusKey, std::size_t>& locus_index,
                        std::vector<Observation>& observations,
                        const std::string& sample, const bam_hdr_t* header,
                        int tid, std::int64_t position, int base, int mapq,
                        bool f1r2, bool indel, const std::vector<Interval>& intervals) {
    if (!selected(intervals, header, tid, position)) return;
    const std::string contig(header->target_name[tid]);
    // TIDs are local to each BAM/CRAM header.  Use the canonical contig name
    // in the aggregation key so the same sample/locus is merged even when
    // multiple inputs list contigs in different orders.
    const auto key = std::make_tuple(sample, contig, position);
    auto found = locus_index.find(key);
    std::size_t locus = 0;
    if (found == locus_index.end()) {
        locus = loci.size();
        locus_index.emplace(key, locus);
        loci.push_back(Locus{sample, contig, tid, position});
    } else {
        locus = found->second;
    }
    observations.push_back(Observation{static_cast<std::int32_t>(locus), base, mapq,
                                       static_cast<std::uint8_t>(f1r2),
                                       static_cast<std::uint8_t>(indel)});
}

struct IngestResult {
    std::size_t records_seen = 0;
    bool indexed_traversal = false;
    std::size_t iterator_intervals = 0;
};

[[maybe_unused]] IngestResult ingest_bam(const std::string& path, const Options& options,
                        const std::vector<Interval>& intervals,
                        std::vector<Locus>& loci,
                        std::map<LocusKey, std::size_t>& locus_index,
                        std::vector<Observation>& observations,
                        std::vector<std::string>& discovered_samples) {
    samFile* input = sam_open(path.c_str(), "r");
    if (input == nullptr) throw std::runtime_error("BAD_INPUT: cannot open BAM/CRAM: " + path);
    bam_hdr_t* header = sam_hdr_read(input);
    if (header == nullptr) {
        sam_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read BAM/CRAM header: " + path);
    }
    const auto rg_samples = read_group_samples(header);
    std::string default_sample = rg_samples.empty() ? "UNKNOWN" : rg_samples.begin()->second;
    if (std::find(discovered_samples.begin(), discovered_samples.end(), default_sample) == discovered_samples.end())
        discovered_samples.push_back(default_sample);
    bam1_t* record = bam_init1();
    if (record == nullptr) {
        bam_hdr_destroy(header);
        sam_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed");
    }
    hts_idx_t* index = nullptr;
    if (!intervals.empty()) index = sam_index_load(input, path.c_str());
    IngestResult result;
    result.indexed_traversal = index != nullptr && !intervals.empty();
    auto process_record = [&]() {
        ++result.records_seen;
        if (!keep_read(record)) return;
        std::string sample = default_sample;
        const auto* rg = bam_aux_get(record, "RG");
        if (rg != nullptr && *rg == 'Z' && bam_aux2Z(rg) != nullptr) {
            const auto found = rg_samples.find(bam_aux2Z(rg));
            if (found != rg_samples.end()) sample = found->second;
        }
        if (std::find(discovered_samples.begin(), discovered_samples.end(), sample) == discovered_samples.end())
            discovered_samples.push_back(sample);
        const bool f1r2 = (((record->core.flag & BAM_FREVERSE) != 0) &&
                           ((record->core.flag & BAM_FREAD1) == 0)) ||
                          (((record->core.flag & BAM_FREVERSE) == 0) &&
                           ((record->core.flag & BAM_FREAD1) != 0));
        const auto* sequence = bam_get_seq(record);
        const auto* qualities = bam_get_qual(record);
        const auto* cigar = bam_get_cigar(record);
        const char* bases = "=ACMGRSVTWYHKDBN";
        std::int64_t reference_position = record->core.pos;
        std::uint32_t read_position = 0;
        for (std::uint32_t index = 0; index < record->core.n_cigar; ++index) {
            const auto op = bam_cigar_op(cigar[index]);
            const auto length = bam_cigar_oplen(cigar[index]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (std::uint32_t offset = 0; offset < length; ++offset) {
                    const auto base_quality = qualities[read_position + offset] == 0xff ? 0 :
                        static_cast<int>(qualities[read_position + offset]);
                    if (base_quality > options.min_base_quality) {
                        const int base = base_index(bases[bam_seqi(sequence, read_position + offset) & 15]);
                        append_observation(loci, locus_index, observations, sample, header,
                                           record->core.tid, reference_position + offset, base,
                                           record->core.qual, f1r2, false, intervals);
                    }
                }
                reference_position += length;
                read_position += length;
            } else if (op == BAM_CINS) {
                if (reference_position > record->core.pos)
                    append_observation(loci, locus_index, observations, sample, header,
                                       record->core.tid, reference_position - 1, -1,
                                       record->core.qual, f1r2, true, intervals);
                read_position += length;
            } else if (op == BAM_CDEL) {
                for (std::uint32_t offset = 0; offset < length; ++offset)
                    append_observation(loci, locus_index, observations, sample, header,
                                       record->core.tid, reference_position + offset, -1,
                                       record->core.qual, f1r2, true, intervals);
                reference_position += length;
            } else if (op == BAM_CREF_SKIP) {
                reference_position += length;
            } else if (op == BAM_CSOFT_CLIP) {
                read_position += length;
            } else if (op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                // no sequence/reference cursor movement
            } else {
                throw std::runtime_error("BAD_INPUT: unsupported CIGAR operator in " + path);
            }
        }
    };
    try {
        if (result.indexed_traversal) {
            for (const auto& interval : intervals) {
                const auto tid = sam_hdr_name2tid(header, interval.contig.c_str());
                if (tid < 0) continue;
                hts_itr_t* iterator = sam_itr_queryi(index, tid, interval.begin, interval.end);
                if (iterator == nullptr) continue;
                ++result.iterator_intervals;
                int status = 0;
                while ((status = sam_itr_next(input, iterator, record)) >= 0) process_record();
                hts_itr_destroy(iterator);
                if (status < -1)
                    throw std::runtime_error("BAD_INPUT: indexed BAM/CRAM iterator failed: " + path);
            }
        } else {
            while (fastgatk::io::read_alignment_record(input, header, record, path) >= 0)
                process_record();
        }
    } catch (...) {
        bam_destroy1(record);
        hts_idx_destroy(index);
        bam_hdr_destroy(header);
        sam_close(input);
        throw;
    }
    bam_destroy1(record);
    bam_hdr_destroy(header);
    hts_idx_destroy(index);
    if (sam_close(input) != 0) throw std::runtime_error("BAD_INPUT: failed reading BAM/CRAM: " + path);
    return result;
}

// Streaming Host decoder used by the bounded pipeline. It deliberately keeps
// only one HtsReader batch plus the current F1R2 observation slice alive;
// locus metadata and mapping-quality vectors are compact aggregate state.
class F1R2StreamDecoder {
public:
    F1R2StreamDecoder(const Options& options,
                      const std::vector<Interval>& intervals,
                      std::vector<Locus>& loci,
                      std::map<LocusKey, std::size_t>& locus_index,
                      std::vector<std::vector<int>>& mapping_qualities,
                      std::vector<std::string>& discovered_samples)
        : options_(options), intervals_(intervals), loci_(loci),
          locus_index_(locus_index), mapping_qualities_(mapping_qualities),
          discovered_samples_(discovered_samples) {
        // A ReadBatch may contain records that produce multiple observations;
        // retaining a small record batch bounds worst-case overshoot while
        // preserving HTSlib indexed traversal.
        const auto reader_batch_records = std::size_t{1};
        for (const auto& input : options_.inputs) {
            auto reader = std::make_unique<fastgatk::io::HtsReader>(
                input, options_.reference, options_.intervals, reader_batch_records);
            const auto& samples = reader->header().read_group_samples;
            const auto default_sample = samples.empty() ? std::string("UNKNOWN") : samples.begin()->second;
            if (std::find(discovered_samples_.begin(), discovered_samples_.end(), default_sample) ==
                discovered_samples_.end())
                discovered_samples_.push_back(default_sample);
            readers_.push_back(std::move(reader));
        }
    }

    F1R2StreamDecoder(const F1R2StreamDecoder&) = delete;
    F1R2StreamDecoder& operator=(const F1R2StreamDecoder&) = delete;

    std::optional<F1R2Batch> next() {
        // Keep draining filtered/empty records iteratively.  A recursive
        // retry here would make a BAM containing a long run of rejected
        // records consume one host stack frame per record, which defeats the
        // bounded-streaming contract on pathological inputs.
        while (active_reader_ < readers_.size()) {
            F1R2Batch result;
            result.observations.reserve(options_.batch_records);
            while (result.observations.size() < options_.batch_records &&
                   active_reader_ < readers_.size()) {
                auto& reader = *readers_[active_reader_];
                if (raw_record_ >= raw_batch_.records()) {
                    raw_batch_.clear();
                    raw_record_ = 0;
                    if (!reader.next(raw_batch_)) {
                        ++active_reader_;
                        continue;
                    }
                }
                process_record(reader, raw_record_, result);
                ++raw_record_;
            }
            if (!result.observations.empty()) return result;
        }
        return std::nullopt;
    }

    std::size_t records_seen() const noexcept {
        std::size_t total = 0;
        for (const auto& reader : readers_) total += reader->records_read();
        return total;
    }

    std::size_t indexed_inputs() const noexcept {
        return static_cast<std::size_t>(std::count_if(readers_.begin(), readers_.end(),
            [](const auto& reader) { return reader->indexed(); }));
    }

    std::size_t sequential_inputs() const noexcept {
        return readers_.size() - indexed_inputs();
    }

    std::size_t iterator_intervals() const noexcept {
        std::size_t total = 0;
        for (const auto& reader : readers_)
            if (reader->indexed()) total += reader->intervals().size();
        return total;
    }

private:
    static bool good_read(const fastgatk::io::ReadBatch& batch, std::size_t record) {
        if (record >= batch.records() || batch.flags.size() != batch.records() ||
            batch.mapq.size() != batch.records() || batch.positions.size() != batch.records() ||
            batch.offsets.size() != batch.records() + 1 ||
            batch.cigar_offsets.size() != batch.records() + 1)
            return false;
        constexpr std::uint16_t rejected_flags = 0x4U | 0x100U | 0x200U | 0x400U | 0x800U;
        if ((batch.flags[record] & rejected_flags) != 0 || batch.mapq[record] < 20 ||
            batch.mapq[record] == 255 || batch.offsets[record + 1] - batch.offsets[record] < 30)
            return false;
        const auto begin = batch.cigar_offsets[record];
        const auto end = batch.cigar_offsets[record + 1];
        if (begin >= end || end > batch.cigar_ops.size()) return false;
        std::uint64_t reference_span = 0;
        for (std::size_t index = begin; index < end; ++index) {
            const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
            if (!operation.valid()) return false;
            if (operation.consumes_reference()) reference_span += operation.length;
        }
        return reference_span > 0;
    }

    bool selected(const fastgatk::io::HeaderSummary& header, int tid,
                  std::int64_t position) const {
        if (intervals_.empty()) return true;
        if (tid < 0 || static_cast<std::size_t>(tid) >= header.contigs.size()) return false;
        const auto& contig = header.contigs[static_cast<std::size_t>(tid)];
        return std::any_of(intervals_.begin(), intervals_.end(), [&](const auto& interval) {
            return interval.contig == contig && position >= interval.begin && position < interval.end;
        });
    }

    std::size_t ensure_locus(const std::string& sample, const std::string& contig,
                             int tid, std::int64_t position, int mapq) {
        const auto key = std::make_tuple(sample, contig, position);
        auto found = locus_index_.find(key);
        std::size_t locus = 0;
        if (found == locus_index_.end()) {
            locus = loci_.size();
            locus_index_.emplace(key, locus);
            loci_.push_back(Locus{sample, contig, tid, position});
            mapping_qualities_.emplace_back();
        } else {
            locus = found->second;
        }
        mapping_qualities_[locus].push_back(mapq);
        return locus;
    }

    void append(F1R2Batch& result, const fastgatk::io::HeaderSummary& header,
                const std::string& sample, int tid, std::int64_t position, int base,
                int mapq, bool f1r2, bool indel) {
        if (!selected(header, tid, position)) return;
        const auto& contig = header.contigs[static_cast<std::size_t>(tid)];
        const auto locus = ensure_locus(sample, contig, tid, position, mapq);
        result.observations.push_back(Observation{
            static_cast<std::int32_t>(locus), base, mapq,
            static_cast<std::uint8_t>(f1r2), static_cast<std::uint8_t>(indel)});
    }

    void process_record(const fastgatk::io::HtsReader& reader, std::size_t record,
                        F1R2Batch& result) {
        const auto& batch = raw_batch_;
        if (!good_read(batch, record)) return;
        const auto tid = batch.tids[record];
        if (tid < 0 || static_cast<std::size_t>(tid) >= reader.header().contigs.size()) return;
        const auto read_begin = batch.offsets[record];
        const auto read_end = batch.offsets[record + 1];
        const auto rg_begin = batch.read_group_offsets.size() == batch.records() + 1
            ? batch.read_group_offsets[record] : 0;
        const auto rg_end = batch.read_group_offsets.size() == batch.records() + 1
            ? batch.read_group_offsets[record + 1] : 0;
        const std::string rg(batch.read_groups.begin() + static_cast<std::ptrdiff_t>(rg_begin),
                             batch.read_groups.begin() + static_cast<std::ptrdiff_t>(rg_end));
        const auto& samples = reader.header().read_group_samples;
        const auto sample_it = samples.find(rg);
        const std::string sample = sample_it == samples.end()
            ? (samples.empty() ? "UNKNOWN" : samples.begin()->second) : sample_it->second;
        if (std::find(discovered_samples_.begin(), discovered_samples_.end(), sample) ==
            discovered_samples_.end())
            discovered_samples_.push_back(sample);
        const auto flags = batch.flags[record];
        const bool f1r2 = (((flags & 0x10U) != 0) && ((flags & 0x40U) == 0)) ||
                          (((flags & 0x10U) == 0) && ((flags & 0x40U) != 0));
        std::int64_t reference_position = batch.positions[record];
        std::uint32_t read_position = 0;
        const auto cigar_begin = batch.cigar_offsets[record];
        const auto cigar_end = batch.cigar_offsets[record + 1];
        for (std::size_t index = cigar_begin; index < cigar_end; ++index) {
            const auto operation = fastgatk::io::CigarOp::unpack(batch.cigar_ops[index]);
            const auto length = operation.length;
            if (operation.projects_base()) {
                if (read_position > read_end - read_begin ||
                    length > read_end - read_begin - read_position)
                    throw std::runtime_error("BAD_INPUT: CIGAR consumes beyond read sequence");
                for (std::uint32_t offset = 0; offset < length; ++offset) {
                    const auto quality = batch.qualities[read_begin + read_position + offset];
                    if (quality <= static_cast<std::uint8_t>(options_.min_base_quality)) continue;
                    const auto base = base_index(static_cast<char>(
                        batch.bases[read_begin + read_position + offset]));
                    append(result, reader.header(), sample, tid, reference_position + offset,
                           base, batch.mapq[record], f1r2, false);
                }
                reference_position += length;
                read_position += length;
            } else if (operation.code == fastgatk::io::CigarOpCode::Insertion) {
                if (read_position > read_end - read_begin ||
                    length > read_end - read_begin - read_position)
                    throw std::runtime_error("BAD_INPUT: CIGAR consumes beyond read sequence");
                if (reference_position > batch.positions[record])
                    append(result, reader.header(), sample, tid, reference_position - 1, -1,
                           batch.mapq[record], f1r2, true);
                read_position += length;
            } else if (operation.code == fastgatk::io::CigarOpCode::Deletion) {
                for (std::uint32_t offset = 0; offset < length; ++offset)
                    append(result, reader.header(), sample, tid, reference_position + offset, -1,
                           batch.mapq[record], f1r2, true);
                reference_position += length;
            } else if (operation.code == fastgatk::io::CigarOpCode::ReferenceSkip) {
                reference_position += length;
            } else if (operation.code == fastgatk::io::CigarOpCode::SoftClip) {
                if (read_position > read_end - read_begin ||
                    length > read_end - read_begin - read_position)
                    throw std::runtime_error("BAD_INPUT: CIGAR consumes beyond read sequence");
                read_position += length;
            } else if (operation.code != fastgatk::io::CigarOpCode::HardClip &&
                       operation.code != fastgatk::io::CigarOpCode::Padding) {
                throw std::runtime_error("BAD_INPUT: unsupported CIGAR operator");
            }
        }
        if (read_position > read_end - read_begin)
            throw std::runtime_error("BAD_INPUT: CIGAR read cursor overflow");
    }

    const Options& options_;
    const std::vector<Interval>& intervals_;
    std::vector<Locus>& loci_;
    std::map<LocusKey, std::size_t>& locus_index_;
    std::vector<std::vector<int>>& mapping_qualities_;
    std::vector<std::string>& discovered_samples_;
    std::vector<std::unique_ptr<fastgatk::io::HtsReader>> readers_;
    std::size_t active_reader_ = 0;
    fastgatk::io::ReadBatch raw_batch_;
    std::size_t raw_record_ = 0;
};

std::vector<std::string> all_contexts() {
    std::vector<std::string> result;
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    result.reserve(64);
    for (const char left : bases)
        for (const char middle : bases)
            for (const char right : bases)
                result.emplace_back(std::string{left, middle, right});
    return result;
}

std::size_t alt_histogram_index(const std::vector<std::string>& contexts,
                                const std::string& context, int alt, bool f1r2) {
    const auto context_it = std::find(contexts.begin(), contexts.end(), context);
    if (context_it == contexts.end()) throw std::logic_error("context missing from histogram index");
    const std::size_t context_index = static_cast<std::size_t>(context_it - contexts.begin());
    int rank = 0;
    const int middle = base_index(context[1]);
    for (int candidate = 0; candidate < 4; ++candidate) {
        if (candidate == middle) continue;
        if (candidate == alt) break;
        ++rank;
    }
    return context_index * 6 + static_cast<std::size_t>(rank * 2 + (f1r2 ? 0 : 1));
}

std::string histogram_content(const SampleOutput& output, const std::vector<std::string>& contexts,
                              bool reference_histogram, int max_depth) {
    std::ostringstream out;
    out << "## htsjdk.samtools.metrics.StringHeader\n# " << output.sample
        << "\n\n\n## HISTOGRAM\tjava.lang.Integer\n";
    out << "depth";
    if (reference_histogram) {
        for (const auto& context : contexts) out << '\t' << context;
    } else {
        for (const auto& context : contexts) {
            const int middle = base_index(context[1]);
            for (int alt = 0; alt < 4; ++alt) {
                if (alt == middle) continue;
                out << '\t' << context << '_' << base_char(alt) << "_F1R2"
                    << '\t' << context << '_' << base_char(alt) << "_F2R1";
            }
        }
    }
    out << '\n';
    for (int depth = 1; depth <= max_depth; ++depth) {
        out << depth;
        if (reference_histogram) {
            for (std::size_t index = 0; index < contexts.size(); ++index)
                out << '\t' << output.ref_hist[index][static_cast<std::size_t>(depth)];
        } else {
            for (std::size_t index = 0; index < output.alt_hist.size(); ++index)
                out << '\t' << output.alt_hist[index][static_cast<std::size_t>(depth)];
        }
        out << '\n';
    }
    return out.str();
}

std::string alt_table_content(const SampleOutput& output) {
    std::ostringstream out;
    out << "#<METADATA>SAMPLE=" << output.sample << '\n'
        << "context\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt\n";
    for (const auto& row : output.alt_rows) {
        out << row.context << '\t' << row.ref_count << '\t' << row.alt_count << '\t'
            << row.ref_f1r2 << '\t' << row.alt_f1r2 << '\t' << row.depth << '\t' << row.alt << '\n';
    }
    return out.str();
}

void append_tar_field(std::array<unsigned char, 512>& header, std::size_t offset,
                      std::size_t length, const std::string& value) {
    const auto count = std::min(length, value.size());
    std::copy_n(value.data(), count, header.data() + offset);
}

void append_tar_octal(std::array<unsigned char, 512>& header, std::size_t offset,
                      std::size_t length, std::uint64_t value) {
    std::ostringstream octal;
    octal << std::oct << value;
    const auto text = octal.str();
    std::fill(header.begin() + offset, header.begin() + offset + length, '0');
    const auto begin = text.size() < length - 1 ? 0 : text.size() - (length - 1);
    std::copy(text.begin() + begin, text.end(),
              header.begin() + offset + (length - 1) - (text.size() - begin));
    header[offset + length - 1] = '\0';
}

void write_member(gzFile stream, const std::string& member, const std::string& content) {
    std::array<unsigned char, 512> header{};
    append_tar_field(header, 0, 100, member);
    append_tar_field(header, 100, 8, "0000644\0");
    append_tar_field(header, 257, 6, "ustar\0");
    append_tar_field(header, 263, 2, "00");
    append_tar_octal(header, 124, 12, content.size());
    std::fill(header.begin() + 148, header.begin() + 156, ' ');
    header[156] = '0';
    std::uint32_t checksum = 0;
    for (const auto byte : header) checksum += byte;
    append_tar_octal(header, 148, 8, checksum);
    auto write = [&](const void* pointer, std::size_t bytes) {
        const auto* data = static_cast<const unsigned char*>(pointer);
        while (bytes != 0) {
            const auto chunk = static_cast<unsigned int>(std::min<std::size_t>(bytes, std::numeric_limits<unsigned int>::max()));
            if (gzwrite(stream, data, chunk) != static_cast<int>(chunk))
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write F1R2 archive");
            data += chunk;
            bytes -= chunk;
        }
    };
    write(header.data(), header.size());
    write(content.data(), content.size());
    const std::array<unsigned char, 512> zeros{};
    const auto remainder = content.size() % 512;
    if (remainder != 0) write(zeros.data(), 512 - remainder);
}

void write_tar_gz(const std::string& output,
                  const std::vector<std::pair<std::string, std::string>>& members) {
    gzFile stream = gzopen(output.c_str(), "wb");
    if (stream == nullptr) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create archive: " + output);
    try {
        for (const auto& member : members) write_member(stream, member.first, member.second);
        const std::array<unsigned char, 1024> zeros{};
        if (gzwrite(stream, zeros.data(), static_cast<unsigned int>(zeros.size())) !=
            static_cast<int>(zeros.size()))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize archive");
        if (gzclose(stream) != Z_OK) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize archive");
    } catch (...) {
        gzclose(stream);
        throw;
    }
}

std::string manifest_path_for(const Options& options) {
    return options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    faidx_t* reference = nullptr;
    try {
        const auto options = parse_options(argc, argv);
        reference = fai_load(options.reference.c_str());
        if (reference == nullptr)
            throw std::runtime_error("BAD_INPUT: cannot load reference FASTA/index: " + options.reference);

        samFile* header_input = sam_open(options.inputs.front().c_str(), "r");
        if (header_input == nullptr) throw std::runtime_error("BAD_INPUT: cannot open input BAM/CRAM");
        bam_hdr_t* first_header = sam_hdr_read(header_input);
        if (first_header == nullptr) {
            sam_close(header_input);
            throw std::runtime_error("BAD_INPUT: cannot read input BAM/CRAM header");
        }
        const auto include_intervals = parse_intervals(options.intervals, first_header);
        const auto excluded_intervals = parse_intervals(options.exclusions, first_header);
        if (options.intervals.empty() && !excluded_intervals.empty())
            throw std::runtime_error("UNSUPPORTED_PARAMETER: --exclude-intervals without -L/--intervals requires full-reference subtraction");
        const auto intervals = subtract_intervals(include_intervals, excluded_intervals);
        if (!options.intervals.empty() && intervals.empty())
            throw std::runtime_error("BAD_INPUT: no intervals remain after applying --exclude-intervals");
        bam_hdr_destroy(first_header);
        sam_close(header_input);

        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;

        std::vector<Locus> loci;
        std::map<LocusKey, std::size_t> locus_index;
        std::vector<std::vector<int>> mapping_qualities;
        std::vector<std::string> discovered_samples;
        F1R2StreamDecoder decoder(options, intervals, loci, locus_index,
                                  mapping_qualities, discovered_samples);
        std::vector<std::uint64_t> counts;
        std::vector<std::uint64_t> f1r2_counts;
        std::vector<std::uint64_t> observation_counts;
        std::vector<std::uint64_t> indel_counts;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        std::uint64_t kernel_batches = 0;
        std::uint64_t kernel_observations = 0;
        double kernel_prepare_seconds = 0.0;
        double kernel_execute_seconds = 0.0;
        fastgatk::runtime::ThreeStagePipeline<F1R2Batch, F1R2BatchResult, F1R2BatchResult>::Metrics pipeline_metrics;
        {
            using Pipeline = fastgatk::runtime::ThreeStagePipeline<F1R2Batch, F1R2BatchResult, F1R2BatchResult>;
            const auto safe_budget = resources.safe_memory_budget_bytes();
            const auto decoded_capacity = safe_budget > 0
                ? std::max<std::uint64_t>(safe_budget / 4, 1)
                : std::max<std::uint64_t>(64ULL * 1024ULL * 1024ULL,
                    static_cast<std::uint64_t>(options.batch_records) * sizeof(Observation) * 2ULL);
            const auto result_bytes = static_cast<std::uint64_t>(options.batch_records) *
                (sizeof(std::int32_t) + sizeof(std::uint64_t) * 10ULL) +
                sizeof(F1R2BatchResult);
            const auto result_floor = safe_budget > 0
                ? std::max<std::uint64_t>(safe_budget / 4, 4096ULL)
                : 64ULL * 1024ULL * 1024ULL;
            const auto result_capacity = std::max<std::uint64_t>(result_bytes * 2ULL, result_floor);
            Pipeline pipeline(
                Pipeline::Limits{decoded_capacity, result_capacity, result_capacity},
                [&decoder]() -> std::optional<F1R2Batch> {
                    return decoder.next();
                },
                [&](F1R2Batch batch) -> std::optional<F1R2BatchResult> {
                    if (batch.observations.empty()) return std::nullopt;
                    std::unordered_map<std::int32_t, std::size_t> local_index;
                    local_index.reserve(batch.observations.size());
                    std::vector<std::int32_t> local_loci;
                    local_loci.reserve(batch.observations.size());
                    for (const auto& observation : batch.observations) {
                        if (local_index.emplace(observation.locus, local_loci.size()).second)
                            local_loci.push_back(observation.locus);
                    }
                    Kokkos::View<std::int32_t*> device_locus("f1r2_observation_locus", batch.observations.size());
                    Kokkos::View<std::int32_t*> device_base("f1r2_observation_base", batch.observations.size());
                    Kokkos::View<std::uint8_t*> device_f1r2("f1r2_observation_orientation", batch.observations.size());
                    Kokkos::View<std::uint8_t*> device_indel("f1r2_observation_indel", batch.observations.size());
                    Kokkos::View<std::uint64_t*> device_counts("f1r2_base_counts", local_loci.size() * 4);
                    Kokkos::View<std::uint64_t*> device_f1r2_counts("f1r2_f1r2_counts", local_loci.size() * 4);
                    Kokkos::View<std::uint64_t*> device_observation_counts("f1r2_observation_counts", local_loci.size());
                    Kokkos::View<std::uint64_t*> device_indel_counts("f1r2_indel_counts", local_loci.size());
                    auto host_locus = Kokkos::create_mirror_view(device_locus);
                    auto host_base = Kokkos::create_mirror_view(device_base);
                    auto host_f1r2 = Kokkos::create_mirror_view(device_f1r2);
                    auto host_indel = Kokkos::create_mirror_view(device_indel);
                    for (std::size_t index = 0; index < batch.observations.size(); ++index) {
                        host_locus(index) = static_cast<std::int32_t>(
                            local_index.at(batch.observations[index].locus));
                        host_base(index) = batch.observations[index].base;
                        host_f1r2(index) = batch.observations[index].f1r2;
                        host_indel(index) = batch.observations[index].indel;
                    }
                    Kokkos::deep_copy(device_locus, host_locus);
                    Kokkos::deep_copy(device_base, host_base);
                    Kokkos::deep_copy(device_f1r2, host_f1r2);
                    Kokkos::deep_copy(device_indel, host_indel);
                    Kokkos::deep_copy(device_counts, std::uint64_t{0});
                    Kokkos::deep_copy(device_f1r2_counts, std::uint64_t{0});
                    Kokkos::deep_copy(device_observation_counts, std::uint64_t{0});
                    Kokkos::deep_copy(device_indel_counts, std::uint64_t{0});
                    fastgatk::core::HostBatch count_host("collect-f1r2-counts-v1");
                    count_host.records = batch.observations.size();
                    count_host.bytes = batch.observations.size() *
                        (sizeof(std::int32_t) * 2 + sizeof(std::uint8_t) * 2);
                    fastgatk::core::KernelPlan<ExecSpace> count_plan("collect-f1r2-counts");
                    count_plan.begin_prepare(count_host);
                    fastgatk::core::DeviceBatch<ExecSpace> count_device(batch.observations.size());
                    count_device.bind("locus", device_locus);
                    count_device.bind("base", device_base);
                    count_device.bind("f1r2", device_f1r2);
                    count_device.bind("indel", device_indel);
                    count_device.bind("counts", device_counts);
                    count_device.bind("f1r2_counts", device_f1r2_counts);
                    count_device.bind("observation_counts", device_observation_counts);
                    count_device.bind("indel_counts", device_indel_counts);
                    ExecSpace().fence();
                    count_plan.end_prepare(count_device);
                    count_plan.begin_execute();
                    Kokkos::parallel_for("collect_f1r2_counts",
                        Kokkos::RangePolicy<ExecSpace>(0, batch.observations.size()),
                        KOKKOS_LAMBDA(const std::size_t index) {
                        const auto locus = static_cast<std::size_t>(device_locus(index));
                        Kokkos::atomic_add(&device_observation_counts(locus), std::uint64_t{1});
                        if (device_indel(index) != 0)
                            Kokkos::atomic_add(&device_indel_counts(locus), std::uint64_t{1});
                        const auto base = device_base(index);
                        if (base >= 0 && base < 4) {
                            Kokkos::atomic_add(&device_counts(
                                locus * 4 + static_cast<std::size_t>(base)), std::uint64_t{1});
                            if (device_f1r2(index) != 0)
                                Kokkos::atomic_add(&device_f1r2_counts(
                                    locus * 4 + static_cast<std::size_t>(base)), std::uint64_t{1});
                        }
                    });
                    ExecSpace().fence();
                    count_plan.end_execute();
                    auto host_counts = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace{}, device_counts);
                    auto host_f1r2_counts = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace{}, device_f1r2_counts);
                    auto host_observation_counts = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace{}, device_observation_counts);
                    auto host_indel_counts = Kokkos::create_mirror_view_and_copy(
                        Kokkos::HostSpace{}, device_indel_counts);
                    F1R2BatchResult result;
                    result.loci = std::move(local_loci);
                    result.counts.resize(result.loci.size() * 4);
                    result.f1r2_counts.resize(result.loci.size() * 4);
                    result.observation_counts.resize(result.loci.size());
                    result.indel_counts.resize(result.loci.size());
                    for (std::size_t index = 0; index < result.counts.size(); ++index) {
                        result.counts[index] = host_counts(index);
                        result.f1r2_counts[index] = host_f1r2_counts(index);
                    }
                    for (std::size_t index = 0; index < result.loci.size(); ++index) {
                        result.observation_counts[index] = host_observation_counts(index);
                        result.indel_counts[index] = host_indel_counts(index);
                    }
                    result.kernel_batches = 1;
                    result.kernel_observations = batch.observations.size();
                    result.kernel_prepare_seconds = count_plan.telemetry().prepare_seconds;
                    result.kernel_execute_seconds = count_plan.telemetry().execute_seconds;
                    return result;
                },
                [](F1R2BatchResult result) -> std::optional<F1R2BatchResult> {
                    return result;
                },
                [&](F1R2BatchResult result) {
                    for (std::size_t index = 0; index < result.loci.size(); ++index) {
                        const auto global_locus = static_cast<std::size_t>(result.loci[index]);
                        if (global_locus >= observation_counts.size()) {
                            const auto new_loci = global_locus + 1;
                            counts.resize(new_loci * 4, 0);
                            f1r2_counts.resize(new_loci * 4, 0);
                            observation_counts.resize(new_loci, 0);
                            indel_counts.resize(new_loci, 0);
                        }
                        for (std::size_t base = 0; base < 4; ++base) {
                            counts[global_locus * 4 + base] += result.counts[index * 4 + base];
                            f1r2_counts[global_locus * 4 + base] +=
                                result.f1r2_counts[index * 4 + base];
                        }
                        observation_counts[global_locus] += result.observation_counts[index];
                        indel_counts[global_locus] += result.indel_counts[index];
                    }
                    kernel_batches += result.kernel_batches;
                    kernel_observations += result.kernel_observations;
                    kernel_prepare_seconds += result.kernel_prepare_seconds;
                    kernel_execute_seconds += result.kernel_execute_seconds;
                },
                [](const F1R2Batch& batch) {
                    return std::max<std::uint64_t>(1, batch.observations.size() * sizeof(Observation));
                },
                [](const F1R2BatchResult& result) {
                    return static_cast<std::uint64_t>(sizeof(F1R2BatchResult) +
                        result.loci.size() * sizeof(std::int32_t) +
                        result.counts.size() * sizeof(std::uint64_t) +
                        result.f1r2_counts.size() * sizeof(std::uint64_t) +
                        result.observation_counts.size() * sizeof(std::uint64_t) +
                        result.indel_counts.size() * sizeof(std::uint64_t));
                },
                [](const F1R2BatchResult& result) {
                    return static_cast<std::uint64_t>(sizeof(F1R2BatchResult) +
                        result.loci.size() * sizeof(std::int32_t) +
                        result.counts.size() * sizeof(std::uint64_t) +
                        result.f1r2_counts.size() * sizeof(std::uint64_t) +
                        result.observation_counts.size() * sizeof(std::uint64_t) +
                        result.indel_counts.size() * sizeof(std::uint64_t));
                });
            pipeline_metrics = pipeline.run();
        }

        const std::size_t records_seen = decoder.records_seen();
        const std::size_t indexed_inputs = decoder.indexed_inputs();
        const std::size_t sequential_inputs = decoder.sequential_inputs();
        const std::size_t iterator_intervals = decoder.iterator_intervals();
        const std::size_t locus_count = loci.size();
        counts.resize(locus_count * 4, 0);
        f1r2_counts.resize(locus_count * 4, 0);
        observation_counts.resize(locus_count, 0);
        indel_counts.resize(locus_count, 0);
        const auto contexts = all_contexts();
        std::map<std::string, SampleOutput> outputs;
        for (const auto& sample : discovered_samples) {
            SampleOutput output;
            output.sample = sample;
            for (auto& histogram : output.ref_hist) histogram.assign(static_cast<std::size_t>(options.max_depth + 1), 0);
            for (auto& histogram : output.alt_hist) histogram.assign(static_cast<std::size_t>(options.max_depth + 1), 0);
            outputs.emplace(sample, std::move(output));
        }
        std::unordered_map<std::string, std::string> reference_cache;
        std::size_t accepted_loci = 0;
        std::size_t reference_loci = 0;
        std::size_t alternate_loci = 0;
        for (std::size_t locus = 0; locus < loci.size(); ++locus) {
            auto& qualities = mapping_qualities[locus];
            if (qualities.empty()) continue;
            std::sort(qualities.begin(), qualities.end());
            const double median = qualities.size() % 2 == 1
                ? qualities[qualities.size() / 2]
                : 0.5 * (qualities[qualities.size() / 2 - 1] + qualities[qualities.size() / 2]);
            const auto depth = counts[locus * 4] + counts[locus * 4 + 1] +
                               counts[locus * 4 + 2] + counts[locus * 4 + 3];
            const bool bad_indel = indel_counts[locus] > depth / 100 ||
                                   (depth == 0 && observation_counts[locus] > 0);
            if (depth == 0 || bad_indel || median < options.min_median_mapq) continue;
            auto sequence_it = reference_cache.find(loci[locus].contig);
            std::string context;
            if (sequence_it == reference_cache.end()) {
                hts_pos_t fetched_length = 0;
                char* sequence = faidx_fetch_seq64(reference, loci[locus].contig.c_str(), 0,
                                                   faidx_seq_len64(reference, loci[locus].contig.c_str()) - 1,
                                                   &fetched_length);
                if (sequence == nullptr) throw std::runtime_error("BAD_INPUT: reference contig not found: " + loci[locus].contig);
                std::string materialized(sequence, sequence + fetched_length);
                free(sequence);
                for (char& base : materialized) base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
                sequence_it = reference_cache.emplace(loci[locus].contig, std::move(materialized)).first;
            }
            const auto& sequence = sequence_it->second;
            if (loci[locus].position < 1 || loci[locus].position + 1 >= static_cast<std::int64_t>(sequence.size())) continue;
            context = sequence.substr(static_cast<std::size_t>(loci[locus].position - 1), 3);
            if (std::any_of(context.begin(), context.end(), [](char base) { return base_index(base) < 0; })) continue;
            const int ref = base_index(context[1]);
            int alt = -1;
            std::uint64_t alt_count = 0;
            for (int base = 0; base < 4; ++base) {
                if (base == ref) continue;
                if (counts[locus * 4 + static_cast<std::size_t>(base)] > alt_count) {
                    alt = base;
                    alt_count = counts[locus * 4 + static_cast<std::size_t>(base)];
                }
            }
            auto& output = outputs.at(loci[locus].sample);
            const auto context_it = std::find(contexts.begin(), contexts.end(), context);
            if (context_it == contexts.end()) continue;
            const auto context_index = static_cast<std::size_t>(context_it - contexts.begin());
            if (alt < 0 || alt_count == 0) {
                output.ref_hist[context_index][static_cast<std::size_t>(std::min<std::uint64_t>(depth, options.max_depth))]++;
                ++reference_loci;
            } else {
                ++alternate_loci;
                const auto ref_count = counts[locus * 4 + static_cast<std::size_t>(ref)];
                const auto ref_f1r2 = f1r2_counts[locus * 4 + static_cast<std::size_t>(ref)];
                const auto alt_f1r2 = f1r2_counts[locus * 4 + static_cast<std::size_t>(alt)];
                if (alt_count == 1) {
                    const auto index_f1r2 = alt_histogram_index(contexts, context, alt, alt_f1r2 == 1);
                    output.alt_hist[index_f1r2][static_cast<std::size_t>(std::min<std::uint64_t>(depth, options.max_depth))]++;
                } else {
                    output.alt_rows.push_back(SampleOutput::AltRow{
                        context, static_cast<int>(ref_count), static_cast<int>(alt_count),
                        static_cast<int>(ref_f1r2), static_cast<int>(alt_f1r2),
                        static_cast<int>(depth), base_char(alt)});
                }
            }
            ++accepted_loci;
        }

        std::vector<std::pair<std::string, std::string>> members;
        for (const auto& [sample, output] : outputs) {
            const auto encoded = url_encode(sample);
            members.emplace_back(encoded + ".ref_histogram", histogram_content(output, contexts, true, options.max_depth));
            members.emplace_back(encoded + ".alt_histogram", histogram_content(output, contexts, false, options.max_depth));
            members.emplace_back(encoded + ".alt_table", alt_table_content(output));
        }
        if (members.empty()) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: no sample outputs were produced");
        write_tar_gz(options.output, members);
        if (!std::filesystem::is_regular_file(options.output) || std::filesystem::file_size(options.output) == 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: output archive is missing or empty");
        std::ofstream manifest(manifest_path_for(options));
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write output manifest");
        manifest << "{\"schema_version\":1,\"tool\":\"CollectF1R2Counts\","
                 << "\"implementation\":\"fastgatk-collect-f1r2-kokkos-v1\",\"status\":\"prototype\","
                 << "\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\","
                 << "\"primary_output_kind\":\"collect-f1r2-counts-tar-gz\","
                 << "\"compatibility\":{\"standard_tar_members\":true,\"bam_pileup\":true,"
                 << "\"mutect2_default_filters\":true,\"indexed_interval_iterator\":true,"
                 << "\"sequential_fallback\":true,\"exclude_intervals\":true,\"bit_identical_to_gatk\":false},"
                 << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
                 << "\",\"kind\":\"collect-f1r2-counts-tar-gz\",\"complete\":true}],"
                 << "\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_files\":" << options.inputs.size()
                 << ",\"records_seen\":" << records_seen
                 << ",\"observations\":" << kernel_observations
                 << ",\"loci\":" << loci.size()
                 << ",\"accepted_loci\":" << accepted_loci
                 << ",\"reference_loci\":" << reference_loci
                 << ",\"alternate_loci\":" << alternate_loci
                 << ",\"samples\":" << outputs.size()
                 << ",\"indexed_inputs\":" << indexed_inputs
                 << ",\"sequential_inputs\":" << sequential_inputs
                 << ",\"iterator_intervals\":" << iterator_intervals
                 << ",\"excluded_intervals\":" << excluded_intervals.size()
                 << ",\"indexed_traversal\":" << (indexed_inputs > 0 ? "true" : "false")
                 << ",\"min_median_mapq\":" << options.min_median_mapq
                 << ",\"min_base_quality\":" << options.min_base_quality
                 << ",\"max_depth\":" << options.max_depth
                 << ",\"batch_records\":" << options.batch_records
                 << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"kernel_execution_space\":\"" << ExecSpace::name() << "\""
                 << ",\"kernel_execution_policy\":\"RangePolicy\""
                 << ",\"kernel_batches\":" << kernel_batches
                 << ",\"kernel_observations\":" << kernel_observations
                 << ",\"kernel_prepare_seconds\":" << kernel_prepare_seconds
                 << ",\"kernel_execute_seconds\":" << kernel_execute_seconds
                 << ",\"pipeline_lifecycle\":\"Host decode->bounded queue->Kokkos compute->encode->sink\""
                 << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                 << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                 << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items
                 << ",\"pipeline_decoded_bytes\":" << pipeline_metrics.decoded_bytes
                 << ",\"pipeline_computed_bytes\":" << pipeline_metrics.computed_bytes
                 << ",\"pipeline_encoded_bytes\":" << pipeline_metrics.encoded_bytes
                 << ",\"pipeline_peak_decoded_bytes\":" << pipeline_metrics.peak_decoded_bytes
                 << ",\"pipeline_peak_computed_bytes\":" << pipeline_metrics.peak_computed_bytes
                 << ",\"pipeline_peak_encoded_bytes\":" << pipeline_metrics.peak_encoded_bytes << "}}\n";
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize output manifest");
        std::cout << "{\"tool\":\"CollectF1R2Counts\",\"status\":\"prototype\","
                  << "\"samples\":" << outputs.size() << ",\"records_seen\":" << records_seen
                  << ",\"observations\":" << kernel_observations << ",\"loci\":" << loci.size()
                  << ",\"accepted_loci\":" << accepted_loci
                  << ",\"indexed_inputs\":" << indexed_inputs
                  << ",\"sequential_inputs\":" << sequential_inputs
                  << ",\"iterator_intervals\":" << iterator_intervals
                  << ",\"kernel_batches\":" << kernel_batches
                  << ",\"kernel_observations\":" << kernel_observations
                  << ",\"kernel_execution_space\":\"" << ExecSpace::name() << "\""
                  << ",\"pipeline_decoded_items\":" << pipeline_metrics.decoded_items
                  << ",\"pipeline_computed_items\":" << pipeline_metrics.computed_items
                  << ",\"pipeline_encoded_items\":" << pipeline_metrics.encoded_items << "}\n";
        fai_destroy(reference);
        reference = nullptr;
        Kokkos::finalize();
        initialized = false;
        return 0;
    } catch (const std::exception& error) {
        if (reference != nullptr) fai_destroy(reference);
        if (initialized) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
