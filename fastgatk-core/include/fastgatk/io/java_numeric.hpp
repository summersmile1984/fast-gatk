#pragma once

#include <array>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fastgatk::io {

inline std::string format_java_shortest(std::string text) {
    const bool negative = !text.empty() && text.front() == '-';
    if (negative) text.erase(text.begin());

    int exponent = 0;
    const auto exponent_position = text.find_first_of("eE");
    if (exponent_position != std::string::npos) {
        exponent = std::stoi(text.substr(exponent_position + 1));
        text.erase(exponent_position);
    }
    const auto point = text.find('.');
    int decimal_position = point == std::string::npos
        ? static_cast<int>(text.size()) : static_cast<int>(point);
    std::string digits = text;
    if (point != std::string::npos) digits.erase(point, 1);
    while (!digits.empty() && digits.front() == '0') {
        digits.erase(digits.begin());
        --decimal_position;
    }
    if (digits.empty()) return negative ? "-0.0" : "0.0";
    decimal_position += exponent;

    // Java uses fixed notation for -3 <= scientific exponent < 7.
    const int scientific_exponent = decimal_position - 1;
    const bool scientific = scientific_exponent < -3 || scientific_exponent >= 7;
    std::string result;
    if (!scientific) {
        if (negative) result.push_back('-');
        if (decimal_position <= 0) {
            result += "0.";
            result.append(static_cast<std::size_t>(-decimal_position), '0');
            result += digits;
        } else if (decimal_position >= static_cast<int>(digits.size())) {
            result += digits;
            result.append(static_cast<std::size_t>(decimal_position - digits.size()), '0');
            result += ".0";
        } else {
            result.append(digits.data(), static_cast<std::size_t>(decimal_position));
            result.push_back('.');
            result.append(digits.data() + decimal_position,
                          digits.size() - static_cast<std::size_t>(decimal_position));
        }
        if (result.find('.') == std::string::npos) result += ".0";
        return result;
    }

    if (negative) result.push_back('-');
    result.push_back(digits.front());
    if (digits.size() > 1) {
        result.push_back('.');
        result.append(digits.data() + 1, digits.size() - 1);
    } else {
        result += ".0";
    }
    result.push_back('E');
    // Double.toString omits a plus sign for positive exponents (e.g. 1.0E7).
    result += scientific_exponent < 0 ? "-" : "";
    result += std::to_string(std::abs(scientific_exponent));
    return result;
}

// Java's Double.toString contract used by htsjdk/GATK table writers.  The
// standard C++ shortest-round-trip conversion supplies the selected digits;
// this adapter applies Java's fixed/scientific cut-over, exponent spelling,
// and explicit .0/negative-zero rules at the file boundary.
inline std::string java_double(double value) {
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return value < 0.0 ? "-Infinity" : "Infinity";
    if (value == 0.0) return std::signbit(value) ? "-0.0" : "0.0";
    std::array<char, 128> buffer{};
    const auto converted = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value,
                                         std::chars_format::general);
    if (converted.ec != std::errc{}) {
        std::ostringstream fallback;
        fallback << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
        return fallback.str();
    }
    return format_java_shortest(std::string(buffer.data(), converted.ptr));
}

// Float.toString follows the same textual cut-over but must first select the
// shortest decimal that round-trips to binary32; promoting to double before
// conversion would expose the float's binary residue (for example 0.2f).
inline std::string java_float(float value) {
    if (std::isnan(value)) return "NaN";
    if (std::isinf(value)) return value < 0.0f ? "-Infinity" : "Infinity";
    if (value == 0.0f) return std::signbit(value) ? "-0.0" : "0.0";
    std::array<char, 128> buffer{};
    const auto converted = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value,
                                         std::chars_format::general);
    if (converted.ec != std::errc{}) {
        std::ostringstream fallback;
        fallback << std::setprecision(std::numeric_limits<float>::max_digits10) << value;
        return fallback.str();
    }
    return format_java_shortest(std::string(buffer.data(), converted.ptr));
}

}  // namespace fastgatk::io
