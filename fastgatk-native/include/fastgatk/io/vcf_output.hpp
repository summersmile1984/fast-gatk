#pragma once

// Shared variant (VCF/BCF) output opening for the native tools.
//
// The pinned GATK/htsjdk reader accepts BCF only as a *plain, uncompressed*
// BCF 2.1 file:
//   * htsjdk's BCF2Codec rejects any minor version other than 1, and
//   * htsjdk has no BCF-over-BGZF codec, so a BGZF-wrapped BCF fails to open
//     with "no suitable codecs found" whatever version it holds.
// The vendored HTSlib is patched accordingly (third_party/HTSLIB_BCF21.patch):
// it writes BCF 2.1 by default and, when asked for
// htsFormat{format=bcf, compression=no_compression}, writes it as a plain
// byte stream instead of putting it inside a BGZF container.
//
// Calling bcf_open(path, "w") - the obvious spelling - does NOT produce BCF:
// bcf_open() is a macro for hts_open(), and mode "w" means *text*, so a
// ".bcf" output would silently receive VCF text.  Use open_variant_output()
// instead of hand-rolling the mode string.

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include <htslib/hts.h>
#include <htslib/vcf.h>

namespace fastgatk::io {

// Case-insensitive suffix test.  Deliberately case-insensitive, unlike the
// per-tool `suffix()` helpers: the file name only selects a container here and
// htsjdk/htsjdk-based tooling is case-insensitive about ".bcf".
inline bool variant_output_has_suffix(const std::string& path, const char* ending) {
    const std::string value(ending);
    if (path.size() < value.size()) return false;
    const std::size_t offset = path.size() - value.size();
    for (std::size_t index = 0; index < value.size(); ++index) {
        const auto lhs = static_cast<char>(std::tolower(static_cast<unsigned char>(path[offset + index])));
        const auto rhs = static_cast<char>(std::tolower(static_cast<unsigned char>(value[index])));
        if (lhs != rhs) return false;
    }
    return true;
}

// True when `path` names a BCF output (.bcf, .bcf.gz or .bcf.bgz).
inline bool is_bcf_output_path(const std::string& path) {
    return variant_output_has_suffix(path, ".bcf") ||
           variant_output_has_suffix(path, ".bcf.gz") ||
           variant_output_has_suffix(path, ".bcf.bgz");
}

// hts_open() mode string for an output path, matching open_variant_output().
// Kept for callers that need the mode rather than a handle.
inline const char* variant_output_mode(const std::string& path) {
    if (is_bcf_output_path(path))
        return variant_output_has_suffix(path, ".bcf") ? "w" : "wbz";
    return variant_output_has_suffix(path, ".gz") || variant_output_has_suffix(path, ".bgz")
        ? "wz" : "w";
}

// Open `path` for VCF/BCF writing with the container the name implies:
//   .bcf                  -> plain (non-BGZF) BCF 2.1, the form GATK reads
//   .bcf.gz / .bcf.bgz    -> BCF 2.1 inside BGZF
//   .vcf.gz / .vcf.bgz    -> BGZF-compressed VCF text
//   anything else         -> plain VCF text
// Returns nullptr on failure, exactly like bcf_open()/hts_open().
inline htsFile* open_variant_output(const std::string& path) {
    if (is_bcf_output_path(path) && variant_output_has_suffix(path, ".bcf")) {
        htsFormat format;
        std::memset(&format, 0, sizeof(format));
        format.category = variant_data;
        format.format = bcf;
        format.compression = no_compression;
        format.compression_level = -1;
        format.version.major = 2;
        format.version.minor = 1;
        return hts_open_format(path.c_str(), "w", &format);
    }
    return bcf_open(path.c_str(), variant_output_mode(path));
}

// BCF outputs are written without an index sidecar: GATK indexes BCF with a
// Tribble index holding BCF byte offsets, which neither the text scanner nor
// tabix can produce, and publishing an index built from a bogus parse would be
// worse than publishing none.  Emits the tool's own warning, matching
// GatherVcfs.
inline void warn_bcf_output_has_no_index(const char* tool, const std::string& path) {
    std::cerr << tool << ": warning: --CREATE_INDEX is not supported for BCF output; "
              << path << " is written without an index sidecar\n";
}

// htsjdk VCFEncoder.formatVCFDouble(): values below 0.01 render as "%.3e",
// values below one as "%.3f", and all others as "%.2f" (a magnitude below
// 1e-20 collapses to the literal "0.00").  Java's Formatter rounds HALF_UP
// on the double's shortest decimal representation, while printf-style
// iostreams round the exact binary expansion half-to-even: the AF half-way
// value 0.5295 becomes "0.530" through GATK but "0.529" through
// std::fixed.  Rebuild the Java result from std::to_chars' shortest
// round-trip digits so INFO and FORMAT doubles are byte-equal to htsjdk.
inline std::string format_vcf_double(const double value) {
    if (std::isnan(value) || std::isinf(value)) return std::to_string(value);
    if (value < 0.01 && std::abs(value) < 1.0e-20) return "0.00";
    const bool negative = std::signbit(value);
    const double magnitude = std::abs(value);
    char buffer[64];
    const auto converted =
        std::to_chars(buffer, buffer + sizeof(buffer), magnitude);
    const std::string text(buffer, converted.ptr);
    // Normalize to an integer digit string plus a power of ten:
    // magnitude = <digits> * 10^power.
    const auto exponent_marker = text.find('e');
    const std::string mantissa = text.substr(0, exponent_marker);
    const int explicit_exponent = exponent_marker == std::string::npos
        ? 0 : std::stoi(text.substr(exponent_marker + 1));
    const auto point = mantissa.find('.');
    const std::string integer = point == std::string::npos
        ? mantissa : mantissa.substr(0, point);
    const std::string fraction = point == std::string::npos
        ? std::string() : mantissa.substr(point + 1);
    std::string digits = integer + fraction;
    int power = explicit_exponent - static_cast<int>(fraction.size());
    // Leading zeros must not inflate the significant-digit count used for
    // the scientific exponent below; the integer value (and thus `power`)
    // is unchanged by stripping them.
    const auto first_digit = digits.find_first_not_of('0');
    if (first_digit == std::string::npos) digits = "0";
    else if (first_digit > 0) digits = digits.substr(first_digit);
    // round_half_up(<digits> * 10^shift) as an integer digit string.
    const auto round_scaled = [&digits](const int shift) {
        if (shift >= 0)
            return digits + std::string(static_cast<std::size_t>(shift), '0');
        const auto cut = static_cast<std::size_t>(-shift);
        if (cut >= digits.size()) {
            const std::string padded =
                std::string(cut - digits.size(), '0') + digits;
            return padded.front() >= '5' ? std::string("1") : std::string("0");
        }
        std::string result = digits.substr(0, digits.size() - cut);
        if (digits[digits.size() - cut] >= '5') {
            std::size_t index = result.size();
            while (index > 0 && result[index - 1] == '9') result[--index] = '0';
            if (index > 0) ++result[index - 1];
            else result.insert(result.begin(), '1');
        }
        return result;
    };
    std::string rendered;
    if (value < 0.01) {
        // "%.3e": one integer digit plus three decimals of the mantissa.
        const int exponent = static_cast<int>(digits.size()) - 1 + power;
        std::string scaled = round_scaled(4 - static_cast<int>(digits.size()));
        int final_exponent = exponent;
        if (scaled.size() > 4U) {  // rounding carried (9.9999 -> 10.000)
            ++final_exponent;
            scaled = scaled.substr(0, 4U);
        }
        rendered.assign(1, scaled.front());
        rendered += '.';
        rendered += scaled.substr(1, 3);
        rendered += 'e';
        rendered += final_exponent < 0 ? '-' : '+';
        const int printed = final_exponent < 0 ? -final_exponent : final_exponent;
        if (printed < 10) rendered += '0';
        rendered += std::to_string(printed);
    } else {
        const int decimals = value < 1.0 ? 3 : 2;
        const std::string scaled =
            round_scaled(power + decimals);
        if (scaled.size() <= static_cast<std::size_t>(decimals)) {
            rendered = "0.";
            rendered += std::string(
                static_cast<std::size_t>(decimals) - scaled.size(), '0');
            rendered += scaled;
        } else {
            rendered = scaled.substr(0, scaled.size() -
                                         static_cast<std::size_t>(decimals));
            rendered += '.';
            rendered += scaled.substr(scaled.size() -
                                      static_cast<std::size_t>(decimals));
        }
    }
    if (negative) rendered.insert(rendered.begin(), '-');
    return rendered;
}

}  // namespace fastgatk::io
