#pragma once

#include <stdexcept>
#include <string>

// GATK/Barclay optional booleans are intentionally kept in a tiny header so
// every native command applies the same fail-closed parsing rule.  A bare
// switch means true; only true/false/1/0 may follow inline or as the next
// token.  The helper has no dependency on HTSlib or Kokkos.
namespace fastgatk::native {

inline bool is_optional_boolean_literal(const std::string& value) {
    return value == "true" || value == "false" || value == "1" || value == "0";
}

inline bool parse_optional_boolean(int& index, int argc, char** argv,
                                   const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    const bool inline_form = argument.rfind(prefix, 0) == 0;
    auto value = inline_form ? argument.substr(prefix.size()) : std::string{};
    // `--flag=` is an explicitly supplied, but empty, value.  Treat it as an
    // invalid literal instead of silently converting it to the bare-switch
    // default; this matches Barclay's fail-closed argument binding.
    if (inline_form && value.empty())
        throw std::invalid_argument(std::string("invalid boolean for ") + name + ": ");
    if (value.empty() && argument == name && index + 1 < argc) {
        const std::string next(argv[index + 1]);
        // There are no positional arguments in these tools.  Consume a
        // non-option token so an invalid literal reports the boolean error
        // at the switch that owns it instead of becoming a later unknown
        // argument; a following option still means a bare true switch.
        if (!next.empty() && next.front() != '-') value = argv[++index];
    }
    if (value.empty()) return true;
    if (!is_optional_boolean_literal(value))
        throw std::invalid_argument(std::string("invalid boolean for ") + name + ": " + value);
    return value == "true" || value == "1";
}

}  // namespace fastgatk::native
