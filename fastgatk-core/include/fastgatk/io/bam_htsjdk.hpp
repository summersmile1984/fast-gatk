#pragma once

// htsjdk BAM serialization conformance for native BAM outputs.
//
// Verified against the pinned htsjdk (via the GATK 4.6.2.0 jar) with
// controlled read->modify->write experiments:
//
//  * Integer aux scalars are re-encoded to the smallest type that fits the
//    VALUE, signed types preferred at each width:
//      [-128,127] -> 'c'   [0,255] -> 'C'   [-32768,32767] -> 's'
//      [0,65535] -> 'S'   int32 rest -> 'i'   uint32-only -> 'I'
//    Other entry kinds (A/Z/H/f/B) pass through byte-identically; B-array
//    subtypes are preserved as stored (htsjdk would normalize by Java array
//    type when it constructs the value itself).
//  * SAMRecord#setAttribute semantics: a tag that already exists is replaced
//    IN PLACE; a new tag is inserted before the first existing tag that is
//    greater under the little-endian 16-bit tag code ((tag[1]<<8)|tag[0]),
//    or appended at the end when none is greater.  The pre-existing tag
//    order is never otherwise changed.
//  * Header text: @HD VN is always emitted as the htsjdk write version
//    1.6; @HD field order is preserved; SO is replaced when the caller sets
//    a sort order; @HD is created when absent.  DT (@RG) / PT (@PG)
//    timestamps are parsed and re-emitted as yyyy-MM-dd'T'HH:mm:ss±HHmm in
//    the process time zone; date-only and timezone-less values parse as
//    UTC.  Unparseable values pass through unchanged.

#include <htslib/sam.h>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

namespace fastgatk::io {

namespace bam_htsjdk_detail {

inline std::uint16_t tag_key(const std::uint8_t* tag) {
    return static_cast<std::uint16_t>((static_cast<std::uint16_t>(tag[1]) << 8) |
                                      static_cast<std::uint16_t>(tag[0]));
}

// Full byte size of one aux entry starting at p (3-byte tag+type header).
inline std::size_t entry_size(const std::uint8_t* p, const std::uint8_t* end) {
    const std::uint8_t type = p[2];
    switch (type) {
        case 'A': case 'c': case 'C': return 4;
        case 's': case 'S': return 5;
        case 'i': case 'I': case 'f': return 7;
        case 'Z': case 'H': {
            const void* nul = std::memchr(p + 3, 0, static_cast<std::size_t>(end - p - 3));
            return nul ? static_cast<std::size_t>(static_cast<const std::uint8_t*>(nul) - p) + 1
                       : static_cast<std::size_t>(end - p);
        }
        case 'B': {
            if (end - p < 8) return static_cast<std::size_t>(end - p);
            std::int32_t count = 0;
            std::memcpy(&count, p + 4, 4);
            const std::uint8_t sub = p[3];
            std::size_t width = 0;
            switch (sub) {
                case 'c': case 'C': width = 1; break;
                case 's': case 'S': width = 2; break;
                case 'i': case 'I': case 'f': width = 4; break;
                default: return static_cast<std::size_t>(end - p);
            }
            return 8 + static_cast<std::size_t>(count) * width;
        }
        default: return static_cast<std::size_t>(end - p);
    }
}

// htsjdk's value -> BAM integer type rule.
inline char integer_type(std::int64_t value) {
    if (value >= -128 && value <= 127) return 'c';
    if (value >= 0 && value <= 255) return 'C';
    if (value >= -32768 && value <= 32767) return 's';
    if (value >= 0 && value <= 65535) return 'S';
    if (value >= -2147483647LL - 1 && value <= 2147483647LL) return 'i';
    return 'I';
}

inline void encode_integer(std::vector<std::uint8_t>& out, const std::uint8_t* tag,
                           std::int64_t value) {
    const char type = integer_type(value);
    out.push_back(tag[0]);
    out.push_back(tag[1]);
    out.push_back(static_cast<std::uint8_t>(type));
    switch (type) {
        case 'c': case 'C': {
            out.push_back(static_cast<std::uint8_t>(value & 0xFF));
            break;
        }
        case 's': case 'S': {
            const std::uint16_t raw = static_cast<std::uint16_t>(value & 0xFFFF);
            out.push_back(static_cast<std::uint8_t>(raw & 0xFF));
            out.push_back(static_cast<std::uint8_t>((raw >> 8) & 0xFF));
            break;
        }
        case 'i': {
            const std::uint32_t raw = static_cast<std::uint32_t>(value & 0xFFFFFFFFLL);
            for (int shift = 0; shift < 32; shift += 8)
                out.push_back(static_cast<std::uint8_t>((raw >> shift) & 0xFF));
            break;
        }
        default: {  // 'I'
            const std::uint32_t raw = static_cast<std::uint32_t>(value & 0xFFFFFFFFLL);
            for (int shift = 0; shift < 32; shift += 8)
                out.push_back(static_cast<std::uint8_t>((raw >> shift) & 0xFF));
            break;
        }
    }
}

inline std::int64_t decode_integer(const std::uint8_t* payload, std::uint8_t type) {
    switch (type) {
        case 'c': return static_cast<std::int8_t>(payload[0]);
        case 'C': return payload[0];
        case 's': {
            std::int16_t v = 0;
            std::memcpy(&v, payload, 2);
            return v;
        }
        case 'S': {
            std::uint16_t v = 0;
            std::memcpy(&v, payload, 2);
            return v;
        }
        case 'i': {
            std::int32_t v = 0;
            std::memcpy(&v, payload, 4);
            return v;
        }
        default: {  // 'I'
            std::uint32_t v = 0;
            std::memcpy(&v, payload, 4);
            return v;
        }
    }
}

// Rewrite the aux block of `record` through `transform` (a callable mapping
// an entry list to an output buffer) and shrink-fit the record.
inline void replace_aux(bam1_t* record, const std::vector<std::uint8_t>& aux) {
    const std::uint8_t* data = record->data;
    const std::size_t prefix = static_cast<std::size_t>(bam_get_aux(record) - data);
    const std::size_t needed = prefix + aux.size();
    if (needed > record->m_data) {
        record->data = static_cast<std::uint8_t*>(std::realloc(record->data, needed));
        record->m_data = static_cast<int>(needed);
    }
    if (!aux.empty())
        std::memcpy(record->data + prefix, aux.data(), aux.size());
    record->l_data = static_cast<int>(needed);
}

}  // namespace bam_htsjdk_detail

// Re-encode the aux block the way htsjdk's BAMRecordCodec does on write:
// integer scalars normalized by value, all other entries copied verbatim,
// entry order preserved.  The record never grows.
inline void normalize_bam_aux(bam1_t* record) {
    namespace d = bam_htsjdk_detail;
    const std::uint8_t* begin = bam_get_aux(record);
    const std::uint8_t* end = record->data + record->l_data;
    std::vector<std::uint8_t> out;
    out.reserve(static_cast<std::size_t>(end - begin));
    for (const std::uint8_t* p = begin; p + 3 <= end;) {
        const std::size_t size = d::entry_size(p, end);
        const std::uint8_t type = p[2];
        if (type == 'c' || type == 'C' || type == 's' || type == 'S' ||
            type == 'i' || type == 'I') {
            d::encode_integer(out, p, d::decode_integer(p + 3, type));
        } else {
            out.insert(out.end(), p, p + size);
        }
        p += size;
    }
    d::replace_aux(record, out);
}

// SAMRecord#setAttribute("TAG", "value") semantics for Z values: replace in
// place when the tag exists, otherwise insert before the first existing tag
// that is greater under the little-endian 16-bit tag code, appending when
// none is greater.  The rewritten block is fully normalized.
inline void set_bam_aux_z(bam1_t* record, const char tag[2], const char* value) {
    namespace d = bam_htsjdk_detail;
    const std::uint8_t* begin = bam_get_aux(record);
    const std::uint8_t* end = record->data + record->l_data;
    const std::uint16_t key =
        static_cast<std::uint16_t>((static_cast<std::uint16_t>(tag[1]) << 8) |
                                   static_cast<std::uint16_t>(tag[0]));
    std::vector<std::uint8_t> out;
    out.reserve(static_cast<std::size_t>(end - begin) + std::strlen(value) + 4);
    bool written = false;
    for (const std::uint8_t* p = begin; p + 3 <= end;) {
        const std::size_t size = d::entry_size(p, end);
        const std::uint16_t entry_key = d::tag_key(p);
        if (!written && entry_key == key) {
            // Replace the value in place; the tag keeps its position.
            out.push_back(p[0]);
            out.push_back(p[1]);
            out.push_back('Z');
            out.insert(out.end(), value, value + std::strlen(value) + 1);
            written = true;
        } else if (!written && entry_key > key) {
            // Insert before the first greater tag.
            out.push_back(static_cast<std::uint8_t>(tag[0]));
            out.push_back(static_cast<std::uint8_t>(tag[1]));
            out.push_back('Z');
            out.insert(out.end(), value, value + std::strlen(value) + 1);
            written = true;
            // fall through: the current entry still has to be copied
            const std::uint8_t type = p[2];
            if (type == 'c' || type == 'C' || type == 's' || type == 'S' ||
                type == 'i' || type == 'I') {
                d::encode_integer(out, p, d::decode_integer(p + 3, type));
            } else {
                out.insert(out.end(), p, p + size);
            }
        } else {
            const std::uint8_t type = p[2];
            if (type == 'c' || type == 'C' || type == 's' || type == 'S' ||
                type == 'i' || type == 'I') {
                d::encode_integer(out, p, d::decode_integer(p + 3, type));
            } else {
                out.insert(out.end(), p, p + size);
            }
        }
        p += size;
    }
    if (!written) {
        out.push_back(static_cast<std::uint8_t>(tag[0]));
        out.push_back(static_cast<std::uint8_t>(tag[1]));
        out.push_back('Z');
        out.insert(out.end(), value, value + std::strlen(value) + 1);
    }
    d::replace_aux(record, out);
}

// ---- header text ------------------------------------------------------

namespace bam_htsjdk_detail {

// yyyy-MM-dd'T'HH:mm:ss with optional ±HHmm (or ±HH:MM); date-only and
// timezone-less forms parse as UTC.  Returns epoch seconds on success.
inline bool parse_htsjdk_date(const std::string& value, std::time_t* epoch) {
    std::tm tm{};
    std::size_t used = 0;
    if (std::sscanf(value.c_str(), "%4d-%2d-%2dT%2d:%2d:%2d",
                    &tm.tm_year, &tm.tm_mon, &tm.tm_mday,
                    &tm.tm_hour, &tm.tm_min, &tm.tm_sec) == 6) {
        used = 19;  // yyyy-MM-ddTHH:mm:ss
    } else if (std::sscanf(value.c_str(), "%4d-%2d-%2d",
                           &tm.tm_year, &tm.tm_mon, &tm.tm_mday) == 3) {
        used = 10;  // yyyy-MM-dd
    } else {
        return false;
    }
    if (value.size() != used && value.size() != used + 5 &&
        value.size() != used + 6)
        return false;
    tm.tm_year -= 1900;
    tm.tm_mon -= 1;
    tm.tm_isdst = 0;
    std::time_t base = timegm(&tm);
    if (value.size() > used) {
        // ±HHmm or ±HH:MM
        const char sign = value[used];
        if (sign != '+' && sign != '-') return false;
        int hours = 0, minutes = 0;
        const char* tail = value.c_str() + used + 1;
        if (value.size() == used + 6) {  // ±HH:MM
            if (std::sscanf(tail, "%2d:%2d", &hours, &minutes) != 2) return false;
        } else if (std::sscanf(tail, "%2d%2d", &hours, &minutes) != 2) {
            return false;
        }
        const int offset = hours * 3600 + minutes * 60;
        base += (sign == '+') ? -offset : offset;
    }
    *epoch = base;
    return true;
}

inline std::string format_htsjdk_date(std::time_t epoch) {
    std::tm local{};
    localtime_r(&epoch, &local);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", &local);
    return buffer;
}

inline bool is_date_field(const std::string& line_type, const std::string& key) {
    return (line_type == "@RG" && key == "DT") ||
           (line_type == "@PG" && key == "PT");
}

}  // namespace bam_htsjdk_detail

// htsjdk SAMTextHeaderCodec write semantics (see the file header).
// sort_order == nullptr keeps the input SO value.
inline std::string htsjdk_header_text(sam_hdr_t* header,
                                      const char* sort_order) {
    namespace d = bam_htsjdk_detail;
    const char* raw = sam_hdr_str(header);
    std::string text = raw ? raw : "";
    std::string out;
    out.reserve(text.size() + 32);
    bool have_hd = false;
    std::size_t start = 0;
    while (start <= text.size()) {
        std::size_t nl = text.find('\n', start);
        if (nl == std::string::npos) nl = text.size();
        std::string line = text.substr(start, nl - start);
        start = nl + 1;
        if (line.empty()) continue;
        const std::string type = line.substr(0, 3);
        if (type == "@HD") {
            have_hd = true;
            // Preserve field order; rewrite VN always, SO when requested.
            std::string rebuilt = "@HD";
            std::size_t p = 3;
            if (p < line.size() && line[p] == '\t') ++p;
            while (p < line.size()) {
                std::size_t tab = line.find('\t', p);
                if (tab == std::string::npos) tab = line.size();
                std::string field = line.substr(p, tab - p);
                p = tab + 1;
                if (field.empty()) continue;
                if (field.rfind("VN:", 0) == 0) field = "VN:1.6";
                else if (sort_order != nullptr && field.rfind("SO:", 0) == 0)
                    field = std::string("SO:") + sort_order;
                rebuilt += "\t" + field;
            }
            if (sort_order != nullptr && rebuilt.find("\tSO:") == std::string::npos)
                rebuilt += std::string("\tSO:") + sort_order;
            out += rebuilt + "\n";
        } else if (type == "@RG" || type == "@PG") {
            std::string rebuilt = type;
            std::size_t p = 3;
            if (p < line.size() && line[p] == '\t') ++p;
            while (p < line.size()) {
                std::size_t tab = line.find('\t', p);
                if (tab == std::string::npos) tab = line.size();
                std::string field = line.substr(p, tab - p);
                p = tab + 1;
                if (field.empty()) continue;
                if (field.size() > 3 && field[2] == ':' &&
                    d::is_date_field(type, field.substr(0, 2))) {
                    std::time_t epoch = 0;
                    if (d::parse_htsjdk_date(field.substr(3), &epoch))
                        field = field.substr(0, 3) + d::format_htsjdk_date(epoch);
                }
                rebuilt += "\t" + field;
            }
            out += rebuilt + "\n";
        } else {
            out += line + "\n";
        }
        if (nl == text.size()) break;
    }
    if (!have_hd) {
        std::string hd = "@HD\tVN:1.6";
        if (sort_order != nullptr) hd += std::string("\tSO:") + sort_order;
        out = hd + "\n" + out;
    }
    return out;
}

// Parse a normalized header text back into an sam_hdr_t for writing.
inline sam_hdr_t* htsjdk_header(sam_hdr_t* header, const char* sort_order) {
    const std::string text = htsjdk_header_text(header, sort_order);
    return sam_hdr_parse(text.size(), text.c_str());
}

}  // namespace fastgatk::io
