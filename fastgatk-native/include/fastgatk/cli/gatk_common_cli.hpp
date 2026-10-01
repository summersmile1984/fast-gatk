#pragma once

#include "fastgatk/cli/gatk_cli_catalog.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace fastgatk::cli {

inline std::string option_name(std::string_view token) {
    const auto equals = token.find('=');
    return std::string(equals == std::string_view::npos ? token : token.substr(0, equals));
}

inline bool is_boolean_token(std::string_view token) {
    return token == "true" || token == "false" || token == "1" || token == "0" ||
           token == "TRUE" || token == "FALSE";
}

// Consume a GATK-public option that this tool's dedicated parser did not
// handle.  Optional booleans accept a following true/false token or an
// embedded '='; value options require a following token or '=' form.
inline bool consume(int& index, int argc, char** argv, std::string_view tool) {
    if (index < 1 || index >= argc) return false;
    const std::string token(argv[index]);
    if (token.size() < 2 || token[0] != '-') return false;
    const std::string name = option_name(token);
    const auto& sets = options_for(tool);
    const bool is_flag = sets.flags.count(name) != 0;
    const bool is_value = sets.values.count(name) != 0;
    if (!is_flag && !is_value) return false;
    if (is_flag) {
        if (token.find('=') != std::string::npos) {
            const auto value = token.substr(token.find('=') + 1);
            if (!is_boolean_token(value))
                throw std::invalid_argument("boolean option expects true or false: " + token);
            return true;
        }
        if (index + 1 < argc && is_boolean_token(argv[index + 1]))
            ++index;
        return true;
    }
    if (token.find('=') != std::string::npos) {
        if (token.size() == name.size() + 1)
            throw std::invalid_argument("empty value for option " + name);
        return true;
    }
    if (index + 1 >= argc || argv[index + 1][0] == '-')
        throw std::invalid_argument("missing value for option " + name);
    ++index;
    return true;
}

}  // namespace fastgatk::cli
