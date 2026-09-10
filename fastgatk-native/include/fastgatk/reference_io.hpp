#pragma once

#include <htslib/faidx.h>
#include <htslib/hts.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fastgatk::reference {

struct SequenceInfo {
    std::string name;
    std::int64_t length = 0;
    std::string md5;
};

struct FastaInfo {
    std::string path;
    std::vector<SequenceInfo> sequences;
};

inline std::string dictionary_path(const std::string& fasta) {
    std::filesystem::path path(fasta);
    path.replace_extension(".dict");
    return path.string();
}

inline std::unordered_map<std::string, SequenceInfo> read_dictionary(const std::string& fasta) {
    std::unordered_map<std::string, SequenceInfo> result;
    std::ifstream input(dictionary_path(fasta));
    if (!input) return result;
    for (std::string line; std::getline(input, line);) {
        if (line.rfind("@SQ", 0) != 0) continue;
        SequenceInfo sequence;
        std::istringstream fields(line);
        std::string field;
        fields >> field;
        while (fields >> field) {
            if (field.rfind("SN:", 0) == 0) sequence.name = field.substr(3);
            else if (field.rfind("LN:", 0) == 0) sequence.length = std::stoll(field.substr(3));
            else if (field.rfind("M5:", 0) == 0) sequence.md5 = field.substr(3);
        }
        if (!sequence.name.empty()) result.emplace(sequence.name, std::move(sequence));
    }
    return result;
}

inline std::string md5_for_sequence(faidx_t* fai, const std::string& contig,
                                    std::int64_t length) {
    hts_pos_t fetched = 0;
    char* sequence = faidx_fetch_seq64(fai, contig.c_str(), 0, length - 1, &fetched);
    if (!sequence || fetched != length) {
        free(sequence);
        throw std::runtime_error("BAD_INPUT: failed to fetch reference sequence: " + contig);
    }
    hts_md5_context* context = hts_md5_init();
    if (!context) {
        free(sequence);
        throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib MD5 context allocation failed");
    }
    hts_md5_update(context, sequence, static_cast<unsigned long>(fetched));
    unsigned char digest[16]{};
    hts_md5_final(digest, context);
    hts_md5_destroy(context);
    free(sequence);
    char hex[33]{};
    hts_md5_hex(hex, digest);
    return std::string(hex);
}

inline std::string fetch_sequence(faidx_t* fai, const std::string& contig,
                                  std::int64_t start, std::int64_t end) {
    hts_pos_t fetched = 0;
    char* sequence = faidx_fetch_seq64(fai, contig.c_str(), start - 1, end - 1, &fetched);
    if (!sequence || fetched != end - start + 1) {
        free(sequence);
        throw std::runtime_error("BAD_INPUT: failed to fetch reference interval: " + contig);
    }
    std::string result(sequence, static_cast<std::size_t>(fetched));
    free(sequence);
    return result;
}

inline FastaInfo load_fasta_info(const std::string& path, std::string_view md5_mode) {
    faidx_t* fai = fai_load(path.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + path);
    try {
        const auto dictionary = read_dictionary(path);
        FastaInfo result;
        result.path = path;
        result.sequences.reserve(static_cast<std::size_t>(faidx_nseq(fai)));
        for (int index = 0; index < faidx_nseq(fai); ++index) {
            const char* name = faidx_iseq(fai, index);
            const auto length = static_cast<std::int64_t>(faidx_seq_len64(fai, name));
            SequenceInfo info{name, length, {}};
            const auto entry = dictionary.find(info.name);
            if (md5_mode == "USE_DICT") {
                if (entry == dictionary.end() || entry->second.md5.empty()) {
                    throw std::invalid_argument("BAD_INPUT: running in USE_DICT mode, but MD5 is missing for sequence " + info.name);
                }
                info.md5 = entry->second.md5;
            } else if (md5_mode == "RECALCULATE_IF_MISSING" && entry != dictionary.end() && !entry->second.md5.empty()) {
                info.md5 = entry->second.md5;
            } else {
                info.md5 = md5_for_sequence(fai, info.name, info.length);
            }
            result.sequences.push_back(std::move(info));
        }
        fai_destroy(fai);
        return result;
    } catch (...) {
        fai_destroy(fai);
        throw;
    }
}

inline std::string basename(const std::string& path) {
    return std::filesystem::path(path).filename().string();
}

}  // namespace fastgatk::reference
