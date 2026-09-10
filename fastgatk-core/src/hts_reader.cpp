#include "fastgatk/io/hts_reader.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#if FASTGATK_HAS_HTSLIB
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#endif

namespace fastgatk::io {

void ReadBatch::clear() {
    offsets.clear();
    bases.clear();
    qualities.clear();
    positions.clear();
    tids.clear();
    mapq.clear();
    cigar_offsets.clear();
    cigar_ops.clear();
    flags.clear();
    mate_tids.clear();
    mate_positions.clear();
    template_lengths.clear();
    name_offsets.clear();
    names.clear();
    read_group_offsets.clear();
    read_groups.clear();
    original_alignment_offsets.clear();
    original_alignments.clear();
    original_alignment_present.clear();
    mate_contig_offsets.clear();
    mate_contigs.clear();
    mate_contig_present.clear();
    sample_ids.clear();
    source_records.clear();
    source_input_ids.clear();
    insertion_quality_offsets.clear();
    insertion_qualities.clear();
    deletion_quality_offsets.clear();
    deletion_qualities.clear();
    flow_tp_offsets.clear();
    flow_tp.clear();
    flow_t0_offsets.clear();
    flow_t0_phred.clear();
    flow_order_offsets.clear();
    flow_orders.clear();
    flow_max_hmer.clear();
}

std::size_t ReadBatch::records() const {
    return positions.size();
}

std::size_t ReadBatch::bases_count() const {
    return bases.size();
}

std::uint64_t ReadBatch::bytes() const noexcept {
    std::uint64_t result = 0;
    const auto add = [&result](std::size_t count, std::size_t width = 1) {
        const auto bytes = static_cast<std::uint64_t>(count) * width;
        result = bytes > std::numeric_limits<std::uint64_t>::max() - result
            ? std::numeric_limits<std::uint64_t>::max() : result + bytes;
    };
    add(offsets.size(), sizeof(std::uint32_t));
    add(bases.size()); add(qualities.size());
    add(positions.size(), sizeof(std::int32_t));
    add(tids.size(), sizeof(std::int32_t)); add(mapq.size());
    add(cigar_offsets.size(), sizeof(std::uint32_t));
    add(cigar_ops.size(), sizeof(std::uint32_t));
    add(flags.size(), sizeof(std::uint16_t));
    add(mate_tids.size(), sizeof(std::int32_t));
    add(mate_positions.size(), sizeof(std::int32_t));
    add(template_lengths.size(), sizeof(std::int32_t));
    add(name_offsets.size(), sizeof(std::uint32_t)); add(names.size());
    add(read_group_offsets.size(), sizeof(std::uint32_t)); add(read_groups.size());
    add(original_alignment_offsets.size(), sizeof(std::uint32_t)); add(original_alignments.size());
    add(original_alignment_present.size());
    add(mate_contig_offsets.size(), sizeof(std::uint32_t)); add(mate_contigs.size());
    add(mate_contig_present.size());
    add(sample_ids.size(), sizeof(std::int32_t));
    add(source_records.size(), sizeof(std::uint32_t));
    add(source_input_ids.size(), sizeof(std::uint32_t));
    add(insertion_quality_offsets.size(), sizeof(std::uint32_t));
    add(insertion_qualities.size());
    add(deletion_quality_offsets.size(), sizeof(std::uint32_t));
    add(deletion_qualities.size());
    add(flow_tp_offsets.size(), sizeof(std::uint32_t)); add(flow_tp.size());
    add(flow_t0_offsets.size(), sizeof(std::uint32_t)); add(flow_t0_phred.size());
    add(flow_order_offsets.size(), sizeof(std::uint32_t)); add(flow_orders.size());
    add(flow_max_hmer.size(), sizeof(std::uint16_t));
    return result;
}

bool ReadBatch::has_cigar() const noexcept {
    return !cigar_offsets.empty();
}

bool ReadBatch::cigar_layout_valid() const noexcept {
    if (cigar_offsets.empty()) return true;
    if (cigar_offsets.size() != records() + 1 || cigar_offsets.front() != 0 ||
        cigar_offsets.back() != cigar_ops.size() || offsets.size() != records() + 1 ||
        offsets.empty() || offsets.back() != bases.size()) return false;
    for (std::size_t i = 1; i < offsets.size(); ++i)
        if (offsets[i] < offsets[i - 1]) return false;
    for (std::size_t i = 1; i < cigar_offsets.size(); ++i)
        if (cigar_offsets[i] < cigar_offsets[i - 1]) return false;
    for (std::size_t record = 0; record < records(); ++record) {
        std::uint64_t read_consumed = 0;
        for (std::size_t i = cigar_offsets[record]; i < cigar_offsets[record + 1]; ++i) {
            const auto operation = CigarOp::unpack(cigar_ops[i]);
            if (!operation.valid()) return false;
            if (operation.consumes_read()) read_consumed += operation.length;
        }
        if (read_consumed != offsets[record + 1] - offsets[record]) {
            // Unmapped BAM records commonly carry sequence but no CIGAR.  They
            // are retained for compatibility and are excluded by projection;
            // mapped records must have a CIGAR that consumes exactly l_qseq.
            const bool unmapped = flags.size() == records() && (flags[record] & 0x4U) != 0;
            if (!(unmapped && read_consumed == 0)) return false;
        }
    }
    return true;
}

bool ReadBatch::cigar_record_layout_valid(const std::size_t record) const noexcept {
    if (cigar_offsets.empty()) return true;
    if (record >= records() || cigar_offsets.size() != records() + 1 ||
        offsets.size() != records() + 1 || offsets.front() > offsets.back() ||
        offsets.back() != bases.size() || cigar_offsets.front() > cigar_offsets.back() ||
        cigar_offsets.back() != cigar_ops.size()) return false;
    if (offsets[record] > offsets[record + 1] ||
        cigar_offsets[record] > cigar_offsets[record + 1]) return false;
    std::uint64_t read_consumed = 0;
    for (std::size_t i = cigar_offsets[record]; i < cigar_offsets[record + 1]; ++i) {
        const auto operation = CigarOp::unpack(cigar_ops[i]);
        if (!operation.valid()) return false;
        if (operation.consumes_read()) read_consumed += operation.length;
    }
    if (read_consumed == offsets[record + 1] - offsets[record]) return true;
    const bool unmapped = flags.size() == records() && (flags[record] & 0x4U) != 0;
    return unmapped && read_consumed == 0;
}

bool ReadBatch::indel_quality_layout_valid() const noexcept {
    const auto valid = [this](const std::vector<std::uint32_t>& offsets,
                              const std::vector<std::uint8_t>& payload) {
        if (offsets.empty()) return payload.empty();
        if (offsets.size() != records() + 1 || offsets.front() != 0 ||
            offsets.back() != payload.size()) return false;
        for (std::size_t i = 1; i < offsets.size(); ++i)
            if (offsets[i] < offsets[i - 1]) return false;
        return true;
    };
    return valid(insertion_quality_offsets, insertion_qualities) &&
           valid(deletion_quality_offsets, deletion_qualities);
}

namespace {

bool checked_add(std::int64_t& value, std::uint32_t increment) noexcept {
    if (increment > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max() - value))
        return false;
    value += static_cast<std::int64_t>(increment);
    return true;
}

bool base_layout_valid(const ReadBatch& batch, std::size_t record,
                       std::size_t& begin, std::size_t& end) noexcept {
    if (record >= batch.records() || batch.offsets.size() != batch.records() + 1 ||
        batch.offsets.front() > batch.offsets.back() ||
        batch.offsets.back() != batch.bases.size() ||
        batch.qualities.size() != batch.bases.size()) return false;
    begin = batch.offsets[record];
    end = batch.offsets[record + 1];
    return begin <= end && end <= batch.bases.size();
}

}  // namespace

bool project_read_offset(const ReadBatch& batch, std::size_t record,
                         std::size_t read_offset, ReadProjection& projection) noexcept {
    std::size_t read_begin = 0, read_end = 0;
    if (!base_layout_valid(batch, record, read_begin, read_end) ||
        read_offset >= read_end - read_begin || record >= batch.positions.size() ||
        batch.positions[record] < 0) return false;

    projection.read_offset = static_cast<std::uint32_t>(read_offset);
    projection.operation = CigarOpCode::Match;
    projection.reference_position = batch.positions[record];
    if (!batch.has_cigar()) {
        return checked_add(projection.reference_position,
                           static_cast<std::uint32_t>(read_offset));
    }
    if (!batch.cigar_record_layout_valid(record))
        return false;

    std::size_t cigar_begin = batch.cigar_offsets[record];
    const std::size_t cigar_end = batch.cigar_offsets[record + 1];
    std::uint64_t read_cursor = 0;
    for (std::size_t i = cigar_begin; i < cigar_end; ++i) {
        const CigarOp operation = CigarOp::unpack(batch.cigar_ops[i]);
        if (!operation.valid()) return false;
        if (operation.projects_base()) {
            if (read_offset >= read_cursor &&
                read_offset - read_cursor < operation.length) {
                if (!checked_add(projection.reference_position,
                                 static_cast<std::uint32_t>(read_offset - read_cursor)))
                    return false;
                projection.operation = operation.code;
                return true;
            }
            read_cursor += operation.length;
            if (!checked_add(projection.reference_position, operation.length)) return false;
        } else {
            if (operation.consumes_read()) read_cursor += operation.length;
            if (operation.consumes_reference() &&
                !checked_add(projection.reference_position, operation.length)) return false;
        }
        if (read_cursor > read_end - read_begin) return false;
    }
    return false;
}

std::int64_t reference_end(const ReadBatch& batch, std::size_t record) noexcept {
    std::size_t read_begin = 0, read_end = 0;
    if (!base_layout_valid(batch, record, read_begin, read_end) ||
        record >= batch.positions.size() || batch.positions[record] < 0) return -1;
    std::int64_t end = batch.positions[record];
    if (!batch.has_cigar())
        return checked_add(end, static_cast<std::uint32_t>(read_end - read_begin)) ? end : -1;
    if (!batch.cigar_record_layout_valid(record)) return -1;
    for (std::size_t i = batch.cigar_offsets[record]; i < batch.cigar_offsets[record + 1]; ++i) {
        const CigarOp operation = CigarOp::unpack(batch.cigar_ops[i]);
        if (!operation.valid()) return -1;
        if (operation.consumes_reference() && !checked_add(end, operation.length)) return -1;
    }
    return end;
}

#if FASTGATK_HAS_HTSLIB

std::string simple_count_sequence_dictionary(const sam_hdr_t* header) {
    // HDF5SimpleCountCollection.write does not serialize the complete BAM
    // header.  It constructs a fresh SAMFileHeader from the sequence
    // dictionary, whose current HTSJDK codec representation is a VN:1.6 HD
    // line followed by the original SQ records (including AS/M5/UR/SP tags),
    // with no SO or RG lines.  Preserve those SQ tags from HTSlib's normalized
    // header instead of reducing the dictionary to name/length pairs.
    std::string result = "@HD\tVN:1.6\n";
    const auto* raw = sam_hdr_str(const_cast<sam_hdr_t*>(header));
    const auto raw_length = sam_hdr_length(const_cast<sam_hdr_t*>(header));
    if (raw != nullptr && raw_length != 0) {
        std::size_t begin = 0;
        while (begin < raw_length) {
            auto end = begin;
            while (end < raw_length && raw[end] != '\n') ++end;
            std::string line(raw + begin, raw + end);
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.rfind("@SQ\t", 0) == 0) {
                result += line;
                result.push_back('\n');
            }
            begin = end < raw_length ? end + 1 : raw_length;
        }
    }
    // sam_hdr_str() should always expose SQ lines for a valid BAM/CRAM/SAM
    // header, but retain a deterministic fallback for synthetic HTSlib
    // headers in tests and unusual readers.
    if (result == "@HD\tVN:1.6\n") {
        for (int index = 0; index < header->n_targets; ++index) {
            result += "@SQ\tSN:";
            result += header->target_name[index];
            result += "\tLN:";
            result += std::to_string(header->target_len[index]);
            result.push_back('\n');
        }
    }
    return result;
}

struct HtsReader::Impl {
    htsFile* file = nullptr;
    bam_hdr_t* header = nullptr;
    bam1_t* record = nullptr;
    hts_idx_t* index = nullptr;
    hts_itr_t* iterator = nullptr;
    HeaderSummary summary;
    std::string backend = "HTSlib 1.22.1";
    std::string source_path;
    std::size_t batch_records = 1024;
    std::vector<HtsInterval> requested_intervals;
    std::vector<HtsInterval> intervals;
    std::vector<HtsInterval> exclusions;
    std::size_t interval_files = 0;
    std::size_t interval_records = 0;
    std::size_t records_read = 0;
    std::size_t active_interval = 0;
    bool indexed = false;
    bool include_unmapped = false;
    std::string selected_sample;
    std::map<std::string, std::string> flow_orders_by_rg;
    std::map<std::string, std::uint16_t> flow_max_hmer_by_rg;

    Impl(const std::string& path, const std::string& reference,
         const std::vector<std::string>& regions, std::size_t batch_size,
         const std::string& sample,
         const std::vector<std::string>& exclusion_regions,
         std::int64_t interval_padding,
         std::int64_t exclusion_padding,
         HtsIntervalSetRule interval_set_rule,
         bool include_unmapped_records)
        : source_path(path), batch_records(batch_size),
          include_unmapped(include_unmapped_records), selected_sample(sample) {
        if (batch_records == 0) throw std::invalid_argument("HTS batch size must be positive");
        if (interval_padding < 0 || exclusion_padding < 0)
            throw std::invalid_argument("HTS interval padding must be non-negative");
        file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open HTS file: " + path);
        if (!reference.empty() && hts_set_fai_filename(file, reference.c_str()) != 0)
            throw std::runtime_error("BAD_INPUT: cannot configure reference for CRAM: " + reference);
        header = sam_hdr_read(file);
        if (!header) throw std::runtime_error("BAD_INPUT: cannot read SAM/BAM/CRAM header: " + path);
        record = bam_init1();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bam_init1 failed");
        for (int i = 0; i < header->n_targets; ++i) {
            summary.contigs.emplace_back(header->target_name[i]);
            summary.contig_lengths.push_back(header->target_len[i]);
            std::string assembly;
            kstring_t value{0, 0, nullptr};
            if (sam_hdr_find_tag_pos(header, "SQ", i, "AS", &value) == 0 && value.s)
                assembly.assign(value.s, value.l);
            if (value.s) std::free(value.s);
            summary.contig_assemblies.push_back(std::move(assembly));
        }
        summary.sequence_dictionary = simple_count_sequence_dictionary(header);
        parse_flow_read_groups();
        if (!selected_sample.empty()) {
            const auto sample_it = std::find(summary.samples.begin(), summary.samples.end(),
                                             selected_sample);
            if (sample_it == summary.samples.end())
                throw std::runtime_error("BAD_INPUT: requested sample is absent from @RG SM: " +
                                         selected_sample);
            // Writers expose one VCF/GVCF sample column.  Put the selected
            // sample first while retaining the remaining header order for
            // diagnostics and downstream metadata consumers.
            std::rotate(summary.samples.begin(), sample_it, sample_it + 1);
        }
        std::vector<std::vector<HtsInterval>> interval_sets;
        for (const auto& region : regions) {
            if (region.empty()) continue;
            std::vector<HtsInterval> parsed;
            parse_selector(region, parsed);
            if (!parsed.empty()) interval_sets.push_back(std::move(parsed));
        }
        if (interval_set_rule == HtsIntervalSetRule::Intersection && !interval_sets.empty()) {
            intervals = std::move(interval_sets.front());
            normalize_intervals(intervals);
            for (std::size_t set_index = 1; set_index < interval_sets.size(); ++set_index) {
                normalize_intervals(interval_sets[set_index]);
                intervals = intersect_intervals(intervals, interval_sets[set_index]);
                if (intervals.empty()) break;
            }
        } else {
            for (auto& parsed : interval_sets)
                intervals.insert(intervals.end(), parsed.begin(), parsed.end());
        }
        for (const auto& region : exclusion_regions) {
            if (region.empty()) continue;
            parse_selector(region, exclusions);
        }
        const auto apply_padding = [this](std::vector<HtsInterval>& values, std::int64_t padding) {
            if (padding == 0) return;
            for (auto& interval : values) {
                interval.start = std::max<std::int64_t>(0, interval.start - padding);
                const auto max_length = static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max());
                interval.end = interval.end > max_length - padding
                    ? max_length : interval.end + padding;
                interval.end = std::min<std::int64_t>(interval.end,
                    header->target_len[interval.tid]);
            }
        };
        apply_padding(intervals, interval_padding);
        apply_padding(exclusions, exclusion_padding);
        requested_intervals = intervals;
        normalize_intervals(intervals);
        normalize_intervals(exclusions);
        // Indexed iteration is opportunistic: BAM/CRAM/SAM inputs without an
        // index retain the historical sequential filtering behavior.  The
        // index is particularly important for region-level streaming, where
        // reopening one reader per bounded tile must not rescan the file.
        if (!intervals.empty()) {
            indexed = has_index();
        }
    }

    bool has_index() noexcept {
        if (index != nullptr) return true;
        index = sam_index_load(file, source_path.c_str());
        return index != nullptr;
    }

    void parse_flow_read_groups() {
        const char* text = sam_hdr_str(header);
        if (text == nullptr) return;
        std::istringstream lines(text);
        std::string line;
        while (std::getline(lines, line)) {
            if (line.rfind("@RG", 0) != 0) continue;
            std::string id;
            std::string flow_order;
            std::string sample;
            std::uint16_t max_hmer = 12;
            std::istringstream fields(line);
            std::string field;
            while (std::getline(fields, field, '\t')) {
                if (field.size() > 3 && field.rfind("ID:", 0) == 0) id = field.substr(3);
                else if (field.size() > 3 && field.rfind("FO:", 0) == 0) flow_order = field.substr(3);
                else if (field.size() > 3 && field.rfind("SM:", 0) == 0) sample = field.substr(3);
                else if (field.size() > 3 && field.rfind("mc:", 0) == 0) {
                    try {
                        const auto parsed = std::stoul(field.substr(3));
                        if (parsed > 0 && parsed < 256) max_hmer = static_cast<std::uint16_t>(parsed);
                    } catch (...) {
                        throw std::runtime_error("BAD_INPUT: malformed @RG mc flow tag");
                    }
                }
            }
            if (!id.empty() && !flow_order.empty()) {
                flow_orders_by_rg[id] = flow_order;
                flow_max_hmer_by_rg[id] = max_hmer;
            }
            if (!id.empty() && !sample.empty()) {
                const auto existing = summary.read_group_samples.find(id);
                if (existing != summary.read_group_samples.end() && existing->second != sample)
                    throw std::runtime_error("BAD_INPUT: @RG ID maps to multiple SM samples: " + id);
                summary.read_group_samples[id] = sample;
            }
            if (!sample.empty() && std::find(summary.samples.begin(), summary.samples.end(), sample) ==
                    summary.samples.end()) {
                summary.samples.push_back(sample);
            }
        }
    }

    ~Impl() {
        if (iterator) hts_itr_destroy(iterator);
        if (index) hts_idx_destroy(index);
        if (record) bam_destroy1(record);
        if (header) bam_hdr_destroy(header);
        if (file) hts_close(file);
    }

    bool next_record() {
        if (!indexed) {
            const auto status = sam_read1(file, header, record);
            if (status >= 0) ++records_read;
            return status >= 0;
        }
        while (active_interval < intervals.size()) {
            if (iterator == nullptr) {
                const auto& interval = intervals[active_interval];
                iterator = sam_itr_queryi(index, interval.tid,
                                          static_cast<hts_pos_t>(interval.start),
                                          static_cast<hts_pos_t>(interval.end));
                if (iterator == nullptr) {
                    ++active_interval;
                    continue;
                }
            }
            const auto status = sam_itr_next(file, iterator, record);
            if (status >= 0) {
                ++records_read;
                return true;
            }
            hts_itr_destroy(iterator);
            iterator = nullptr;
            ++active_interval;
        }
        return false;
    }

    static std::string trim(std::string value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return {};
        const auto last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    }

    static bool suffix_ci(const std::string& value, const char* ending) {
        const std::string suffix(ending);
        if (value.size() < suffix.size()) return false;
        const auto offset = value.size() - suffix.size();
        for (std::size_t i = 0; i < suffix.size(); ++i) {
            if (std::tolower(static_cast<unsigned char>(value[offset + i])) !=
                std::tolower(static_cast<unsigned char>(suffix[i]))) return false;
        }
        return true;
    }

    int find_tid(const std::string& contig) const {
        for (int i = 0; i < header->n_targets; ++i)
            if (contig == header->target_name[i]) return i;
        return -1;
    }

    static std::int64_t parse_coordinate(const std::string& text) {
        if (text.empty()) throw std::runtime_error("BAD_INPUT: empty interval coordinate");
        std::size_t consumed = 0;
        const auto value = std::stoll(text, &consumed);
        if (consumed != text.size() || value < 1)
            throw std::runtime_error("BAD_INPUT: invalid interval coordinate: " + text);
        return value;
    }

    void append_literal(const std::string& raw, std::vector<HtsInterval>& target) {
        const auto region = trim(raw);
        const auto colon = region.find(':');
        const auto contig = region.substr(0, colon);
        const auto tid = find_tid(contig);
        if (tid < 0) throw std::runtime_error("BAD_INPUT: region contig not in header: " + contig);
        if (colon == std::string::npos) {
            target.push_back(HtsInterval{tid, 0, header->target_len[tid]});
            return;
        }
        const auto coordinates = region.substr(colon + 1);
        const auto dash = coordinates.find('-');
        const auto start = parse_coordinate(coordinates.substr(0, dash));
        const auto end = dash == std::string::npos
            ? start : parse_coordinate(coordinates.substr(dash + 1));
        if (end < start || end > header->target_len[tid])
            throw std::runtime_error("BAD_INPUT: interval outside contig: " + region);
        target.push_back(HtsInterval{tid, start - 1, end});
    }

    void append_file(const std::string& path, std::vector<HtsInterval>& target) {
        htsFile* input = hts_open(path.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        const bool bed = suffix_ci(path, ".bed") || suffix_ci(path, ".bed.gz");
        std::size_t parsed = 0;
        kstring_t raw_line{0, 0, nullptr};
        try {
        while (hts_getline(input, '\n', &raw_line) >= 0) {
            std::string line = trim(raw_line.s == nullptr ? std::string{} :
                                    std::string(raw_line.s, raw_line.l));
            if (line.empty() || line[0] == '#' || line[0] == '@') continue;
            std::istringstream fields(line);
            std::string first;
            fields >> first;
            if (first.empty()) continue;
            if (first == "track" || first == "browser") continue;
            if (first.find(':') != std::string::npos) {
                append_literal(first, target);
                ++parsed;
                continue;
            }
            std::string second, third;
            if (!(fields >> second >> third))
                throw std::runtime_error("BAD_INPUT: malformed interval file line: " + line);
            const auto tid = find_tid(first);
            if (tid < 0) throw std::runtime_error("BAD_INPUT: interval contig not in header: " + first);
            std::size_t consumed_start = 0, consumed_end = 0;
            const auto start = std::stoll(second, &consumed_start);
            const auto end = std::stoll(third, &consumed_end);
            if (consumed_start != second.size() || consumed_end != third.size())
                throw std::runtime_error("BAD_INPUT: malformed interval coordinates: " + line);
            const auto begin = bed ? start : start - 1;
            const auto finish = end;
            if (begin < 0 || finish <= begin || finish > header->target_len[tid])
                throw std::runtime_error("BAD_INPUT: interval outside contig: " + line);
            target.push_back(HtsInterval{tid, begin, finish});
            ++parsed;
        }
        } catch (...) {
            free(raw_line.s);
            hts_close(input);
            throw;
        }
        const auto close_status = hts_close(input);
        free(raw_line.s);
        if (close_status != 0)
            throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
        if (parsed == 0) throw std::runtime_error("BAD_INPUT: interval file contains no intervals: " + path);
        ++interval_files;
        interval_records += parsed;
    }

    void parse_selector(const std::string& selector,
                        std::vector<HtsInterval>& target) {
        std::error_code error;
        if (std::filesystem::is_regular_file(selector, error) && !error) append_file(selector, target);
        else append_literal(selector, target);
    }

    static void normalize_intervals(std::vector<HtsInterval>& values) {
        std::sort(values.begin(), values.end(), [](const auto& left, const auto& right) {
            if (left.tid != right.tid) return left.tid < right.tid;
            if (left.start != right.start) return left.start < right.start;
            return left.end < right.end;
        });
        std::vector<HtsInterval> merged;
        for (const auto& interval : values) {
            if (!merged.empty() && merged.back().tid == interval.tid &&
                interval.start <= merged.back().end) {
                merged.back().end = std::max(merged.back().end, interval.end);
            } else merged.push_back(interval);
        }
        values.swap(merged);
    }

    static std::vector<HtsInterval> intersect_intervals(
            const std::vector<HtsInterval>& left,
            const std::vector<HtsInterval>& right) {
        std::vector<HtsInterval> result;
        for (const auto& first : left) {
            for (const auto& second : right) {
                if (first.tid != second.tid) continue;
                const auto start = std::max(first.start, second.start);
                const auto end = std::min(first.end, second.end);
                if (start < end) result.push_back(HtsInterval{first.tid, start, end});
            }
        }
        normalize_intervals(result);
        return result;
    }

    bool selected(const bam1_t* record) const {
        // Most callers apply a tool-specific read filter that excludes
        // unmapped records.  Smoke/read-metrics adapters can opt in to the
        // GATK unbounded traversal behavior explicitly; interval traversals
        // still cannot select an unmapped record.
        if (record->core.tid < 0 || record->core.pos < 0)
            return include_unmapped && intervals.empty();
        std::int64_t span = 0;
        const std::uint32_t* cigar = bam_get_cigar(record);
        for (std::uint32_t i = 0; i < record->core.n_cigar; ++i) {
            const auto operation = CigarOp::unpack(cigar[i]);
            if (operation.valid() && operation.consumes_reference())
                span += operation.length;
        }
        // Unmapped/legacy records without a reference-consuming CIGAR retain
        // the historical qseq overlap behavior for region filtering.
        if (span == 0) span = record->core.l_qseq;
        const auto begin = static_cast<std::int64_t>(record->core.pos);
        const auto end = begin + span;
        const auto overlaps = [&](const std::vector<HtsInterval>& values) {
            return std::any_of(values.begin(), values.end(), [&](const auto& interval) {
                return record->core.tid == interval.tid && begin < interval.end && end > interval.start;
            });
        };
        if (!intervals.empty() && !overlaps(intervals)) return false;
        // Exclusions are intentionally not applied at whole-read selection:
        // a read crossing an excluded locus can still provide evidence for an
        // included locus.  The calling pipeline consumes exclusion_intervals
        // per projected base, preserving GATK traversal semantics.
        return true;
    }

    static void append_string(std::vector<std::uint8_t>& payload,
                              std::vector<std::uint32_t>& offsets,
                              const char* value) {
        const std::size_t length = value == nullptr ? 0 : std::strlen(value);
        if (payload.size() > std::numeric_limits<std::uint32_t>::max() - length)
            throw std::runtime_error("RESOURCE_EXHAUSTED: HTS metadata exceeds uint32 offsets");
        if (value != nullptr)
            payload.insert(payload.end(), value, value + length);
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    }

    static void append_tp_tag(const bam1_t* record,
                              std::vector<std::int8_t>& payload,
                              std::vector<std::uint32_t>& offsets) {
        const auto* aux = bam_aux_get(record, "tp");
        if (aux == nullptr) {
            offsets.push_back(static_cast<std::uint32_t>(payload.size()));
            return;
        }
        if (*aux != 'B') throw std::runtime_error("BAD_INPUT: flow tp tag must be a BAM B array");
        const auto subtype = static_cast<char>(aux[1]);
        const auto length = bam_auxB_len(aux);
        if (payload.size() > std::numeric_limits<std::uint32_t>::max() - length)
            throw std::runtime_error("RESOURCE_EXHAUSTED: flow tp offsets exceed uint32 range");
        for (std::uint32_t index = 0; index < length; ++index) {
            auto value = bam_auxB2i(aux, index);
            if ((subtype == 'C' || subtype == 'c') && subtype == 'C' && value > 127) value -= 256;
            if (value < -128 || value > 127)
                throw std::runtime_error("BAD_INPUT: flow tp value does not fit signed byte");
            payload.push_back(static_cast<std::int8_t>(value));
        }
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    }

    static void append_t0_tag(const bam1_t* record,
                              std::vector<std::uint8_t>& payload,
                              std::vector<std::uint32_t>& offsets) {
        const auto* aux = bam_aux_get(record, "t0");
        if (aux == nullptr) {
            offsets.push_back(static_cast<std::uint32_t>(payload.size()));
            return;
        }
        if (*aux != 'Z') throw std::runtime_error("BAD_INPUT: flow t0 tag must be a SAM Z string");
        const char* value = bam_aux2Z(aux);
        if (value == nullptr) throw std::runtime_error("BAD_INPUT: malformed flow t0 tag");
        const auto length = std::strlen(value);
        if (payload.size() > std::numeric_limits<std::uint32_t>::max() - length)
            throw std::runtime_error("RESOURCE_EXHAUSTED: flow t0 offsets exceed uint32 range");
        for (std::size_t index = 0; index < length; ++index) {
            const int phred = static_cast<unsigned char>(value[index]) - 33;
            if (phred < 0 || phred > 255)
                throw std::runtime_error("BAD_INPUT: flow t0 FASTQ quality is out of range");
            payload.push_back(static_cast<std::uint8_t>(phred));
        }
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    }

    static void append_indel_quality_tag(const bam1_t* record, const char* tag,
                                         std::vector<std::uint8_t>& payload,
                                         std::vector<std::uint32_t>& offsets) {
        const auto* aux = bam_aux_get(record, tag);
        if (aux == nullptr) {
            offsets.push_back(static_cast<std::uint32_t>(payload.size()));
            return;
        }
        if (*aux != 'Z')
            throw std::runtime_error(std::string("BAD_INPUT: ") + tag +
                                     " indel quality tag must be a SAM Z string");
        const char* value = bam_aux2Z(aux);
        if (value == nullptr)
            throw std::runtime_error(std::string("BAD_INPUT: malformed ") + tag +
                                     " indel quality tag");
        const auto length = std::strlen(value);
        if (payload.size() > std::numeric_limits<std::uint32_t>::max() - length)
            throw std::runtime_error("RESOURCE_EXHAUSTED: indel quality offsets exceed uint32 range");
        for (std::size_t index = 0; index < length; ++index) {
            const int phred = static_cast<unsigned char>(value[index]) - 33;
            if (phred < 0 || phred > 255)
                throw std::runtime_error(std::string("BAD_INPUT: ") + tag +
                                         " indel quality is outside FASTQ byte range");
            payload.push_back(static_cast<std::uint8_t>(phred));
        }
        offsets.push_back(static_cast<std::uint32_t>(payload.size()));
    }

    bool next(ReadBatch& batch) {
        batch.clear();
        batch.offsets.push_back(0);
        batch.cigar_offsets.push_back(0);
        batch.name_offsets.push_back(0);
        batch.read_group_offsets.push_back(0);
        batch.original_alignment_offsets.push_back(0);
        batch.mate_contig_offsets.push_back(0);
        batch.insertion_quality_offsets.push_back(0);
        batch.deletion_quality_offsets.push_back(0);
        batch.flow_tp_offsets.push_back(0);
        batch.flow_t0_offsets.push_back(0);
        batch.flow_order_offsets.push_back(0);
        const char* bases = "=ACMGRSVTWYHKDBN";
        while (batch.records() < batch_records) {
            if (!next_record()) break;
            if (!selected(record)) continue;
            if (!selected_sample.empty()) {
                const auto* rg = bam_aux_get(record, "RG");
                const char* rg_name = rg != nullptr && *rg == 'Z' ? bam_aux2Z(rg) : nullptr;
                const auto sample_it = rg_name == nullptr
                    ? summary.read_group_samples.end()
                    : summary.read_group_samples.find(rg_name);
                if (sample_it == summary.read_group_samples.end() ||
                    sample_it->second != selected_sample)
                    continue;
            }
            const int length = record->core.l_qseq;
            const std::uint8_t* sequence = bam_get_seq(record);
            const std::uint8_t* qualities = bam_get_qual(record);
            for (int i = 0; i < length; ++i) {
                const auto code = bam_seqi(sequence, i);
                const char base = bases[code & 15];
                batch.bases.push_back(static_cast<std::uint8_t>(base));
                batch.qualities.push_back(qualities[i] == 0xff ? 0 : qualities[i]);
            }
            batch.offsets.push_back(static_cast<std::uint32_t>(batch.bases.size()));
            const std::uint32_t* cigar = bam_get_cigar(record);
            if (batch.cigar_ops.size() > std::numeric_limits<std::uint32_t>::max() - record->core.n_cigar)
                throw std::runtime_error("RESOURCE_EXHAUSTED: HTS CIGAR offsets exceed uint32 range");
            for (std::uint32_t i = 0; i < record->core.n_cigar; ++i)
                batch.cigar_ops.push_back(cigar[i]);
            batch.cigar_offsets.push_back(static_cast<std::uint32_t>(batch.cigar_ops.size()));
            batch.positions.push_back(record->core.pos);
            batch.tids.push_back(record->core.tid);
            batch.mapq.push_back(record->core.qual);
            batch.flags.push_back(record->core.flag);
            batch.mate_tids.push_back(record->core.mtid);
            batch.mate_positions.push_back(record->core.mpos);
            batch.template_lengths.push_back(record->core.isize);
            append_string(batch.names, batch.name_offsets, bam_get_qname(record));
            const std::uint8_t* rg = bam_aux_get(record, "RG");
            append_string(batch.read_groups, batch.read_group_offsets,
                          rg != nullptr && *rg == 'Z' ? bam_aux2Z(rg) : nullptr);
            const std::uint8_t* oa = bam_aux_get(record, "OA");
            batch.original_alignment_present.push_back(oa != nullptr ? 1U : 0U);
            append_string(batch.original_alignments, batch.original_alignment_offsets,
                          oa != nullptr && *oa == 'Z' ? bam_aux2Z(oa) : nullptr);
            const std::uint8_t* xm = bam_aux_get(record, "XM");
            batch.mate_contig_present.push_back(xm != nullptr ? 1U : 0U);
            append_string(batch.mate_contigs, batch.mate_contig_offsets,
                          xm != nullptr && *xm == 'Z' ? bam_aux2Z(xm) : nullptr);
            append_indel_quality_tag(record, "BI", batch.insertion_qualities,
                                     batch.insertion_quality_offsets);
            append_indel_quality_tag(record, "BD", batch.deletion_qualities,
                                     batch.deletion_quality_offsets);
            append_tp_tag(record, batch.flow_tp, batch.flow_tp_offsets);
            append_t0_tag(record, batch.flow_t0_phred, batch.flow_t0_offsets);
            std::string rg_name = rg != nullptr && *rg == 'Z' && bam_aux2Z(rg) != nullptr
                ? bam_aux2Z(rg) : std::string{};
            const auto flow = flow_orders_by_rg.find(rg_name);
            if (flow == flow_orders_by_rg.end()) {
                batch.flow_order_offsets.push_back(static_cast<std::uint32_t>(batch.flow_orders.size()));
                batch.flow_max_hmer.push_back(0);
            } else {
                if (batch.flow_orders.size() > std::numeric_limits<std::uint32_t>::max() - flow->second.size())
                    throw std::runtime_error("RESOURCE_EXHAUSTED: flow order offsets exceed uint32 range");
                batch.flow_orders.insert(batch.flow_orders.end(), flow->second.begin(), flow->second.end());
                batch.flow_order_offsets.push_back(static_cast<std::uint32_t>(batch.flow_orders.size()));
                batch.flow_max_hmer.push_back(flow_max_hmer_by_rg[rg_name]);
            }
        }
        return !batch.positions.empty();
    }
};

#else

struct HtsReader::Impl {
    HeaderSummary summary;
    std::string backend = "unavailable (HTSlib not compiled)";
};

#endif

HtsReader::HtsReader(const std::string& path, const std::string& reference,
                     const std::string& region, std::size_t batch_records,
                     const std::string& sample, bool include_unmapped)
#if FASTGATK_HAS_HTSLIB
    : HtsReader(path, reference, region.empty() ? std::vector<std::string>{}
                                              : std::vector<std::string>{region}, batch_records,
                sample, {}, 0, 0, HtsIntervalSetRule::Union, include_unmapped) {}
#else
    : impl_(std::make_unique<Impl>()) {
    (void)path; (void)reference; (void)region; (void)batch_records; (void)sample;
    (void)include_unmapped;
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with FAST_GATK_HTSLIB_ROOT to read BAM/CRAM");
}
#endif

HtsReader::~HtsReader() = default;
HtsReader::HtsReader(HtsReader&&) noexcept = default;
HtsReader& HtsReader::operator=(HtsReader&&) noexcept = default;

bool HtsReader::next(ReadBatch& batch) {
#if FASTGATK_HAS_HTSLIB
    return impl_->next(batch);
#else
    (void)batch;
    throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib reader is not compiled");
#endif
}

void HtsReader::set_batch_records(std::size_t batch_records) {
#if FASTGATK_HAS_HTSLIB
    if (batch_records == 0) throw std::invalid_argument("HTS batch size must be positive");
    impl_->batch_records = batch_records;
#else
    (void)batch_records;
    throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib reader is not compiled");
#endif
}

std::size_t HtsReader::batch_records() const noexcept {
#if FASTGATK_HAS_HTSLIB
    return impl_->batch_records;
#else
    return 0;
#endif
}

HtsReader::HtsReader(const std::string& path, const std::string& reference,
                     const std::vector<std::string>& regions, std::size_t batch_records,
                     const std::string& sample,
                     const std::vector<std::string>& exclusions,
                     std::int64_t interval_padding,
                     std::int64_t exclusion_padding,
                     HtsIntervalSetRule interval_set_rule,
                     bool include_unmapped)
#if FASTGATK_HAS_HTSLIB
    : impl_(std::make_unique<Impl>(path, reference, regions, batch_records, sample,
                                   exclusions, interval_padding, exclusion_padding,
                                   interval_set_rule, include_unmapped)) {}
#else
    : impl_(std::make_unique<Impl>()) {
    (void)path; (void)reference; (void)regions; (void)batch_records; (void)sample;
    (void)exclusions; (void)interval_padding; (void)exclusion_padding;
    (void)interval_set_rule; (void)include_unmapped;
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with FAST_GATK_HTSLIB_ROOT to read BAM/CRAM");
}
#endif

bool HtsReader::compiled() const {
#if FASTGATK_HAS_HTSLIB
    return true;
#else
    return false;
#endif
}

bool HtsReader::indexed() const noexcept {
#if FASTGATK_HAS_HTSLIB
    return impl_->indexed;
#else
    return false;
#endif
}

bool HtsReader::has_index() const noexcept {
#if FASTGATK_HAS_HTSLIB
    return impl_->has_index();
#else
    return false;
#endif
}

std::size_t HtsReader::records_read() const noexcept {
#if FASTGATK_HAS_HTSLIB
    return impl_->records_read;
#else
    return 0;
#endif
}

const std::string& HtsReader::backend_description() const { return impl_->backend; }
const HeaderSummary& HtsReader::header() const { return impl_->summary; }
const std::vector<HtsInterval>& HtsReader::requested_intervals() const {
#if FASTGATK_HAS_HTSLIB
    return impl_->requested_intervals;
#else
    static const std::vector<HtsInterval> empty;
    return empty;
#endif
}
const std::vector<HtsInterval>& HtsReader::intervals() const {
#if FASTGATK_HAS_HTSLIB
    return impl_->intervals;
#else
    static const std::vector<HtsInterval> empty;
    return empty;
#endif
}
const std::vector<HtsInterval>& HtsReader::exclusion_intervals() const {
#if FASTGATK_HAS_HTSLIB
    return impl_->exclusions;
#else
    static const std::vector<HtsInterval> empty;
    return empty;
#endif
}
std::size_t HtsReader::interval_file_inputs() const {
#if FASTGATK_HAS_HTSLIB
    return impl_->interval_files;
#else
    return 0;
#endif
}
std::size_t HtsReader::interval_file_records() const {
#if FASTGATK_HAS_HTSLIB
    return impl_->interval_records;
#else
    return 0;
#endif
}

}  // namespace fastgatk::io
