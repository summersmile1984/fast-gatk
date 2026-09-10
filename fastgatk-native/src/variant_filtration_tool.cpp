#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/io/intervals.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/tribble_index.hpp"
#include "optional_boolean.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string output;
    std::string manifest;
    std::string reference;
    std::vector<std::string> regions;
    fastgatk::io::HtsIntervalSetRule interval_set_rule =
        fastgatk::io::HtsIntervalSetRule::Union;
    std::vector<std::string> expressions;
    std::vector<std::string> names;
    std::vector<std::string> genotype_expressions;
    std::vector<std::string> genotype_names;
    std::vector<std::string> masks;
    std::string mask_name = "Mask";
    std::string mask_description;
    int mask_extension = 0;
    bool filter_records_not_in_mask = false;
    int cluster_size = 0;
    int cluster_window_size = 0;
    bool missing_fails = false;
    bool invert_filter_expression = false;
    bool invert_genotype_filter_expression = false;
    bool apply_allele_specific_filters = false;
    bool invalidate_previous_filters = false;
    bool set_filtered_genotypes_to_nocall = false;
    bool create_index = true;
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
            std::cout << "fastgatk-variant-filtration (GATK-compatible native prototype)\n"
                         "  -V, --variant FILE             input VCF/BCF\n"
                         "  -O, --output FILE              filtered VCF/BCF\n"
                         "  --filter-expression EXPR       JEXL subset: QUAL/INFO/FORMAT and Allele comparisons, &&/||/!\n"
                         "  --filter-name NAME             FILTER label (repeatable)\n"
                         "  --genotype-filter-expression EXPR  genotype JEXL predicate\n"
                         "  --genotype-filter-name NAME       genotype FT label (repeatable)\n"
                         "  --mask FILE                       mask VCF/BCF (repeatable)\n"
                         "  --mask-name NAME                  mask FILTER label\n"
                         "  --mask-description TEXT           mask FILTER header description\n"
                         "  --mask-extension N                extend mask intervals by N bases\n"
                         "  --filter-not-in-mask              filter records outside the mask\n"
                         "  --interval-set-rule RULE          UNION (default) or INTERSECTION\n"
                         "  --cluster-size N                  clustered-event count threshold\n"
                         "  --cluster-window-size N           clustered-event window in bases\n"
                         "  --apply-allele-specific-filters   evaluate AS_* Number=A annotations per ALT\n"
                         "  --invalidate-previous-filters [BOOL]  clear existing site FILTER labels\n"
                         "  --set-filtered-genotype-to-no-call[=BOOL]  set filtered GT alleles to ./.\n"
                         "  --missing-values-evaluate-as-failing\n"
                         "      --output-manifest FILE    OutputManifest JSON\n";
            std::exit(0);
        }
        if (argument == "-V" || is_option(argument, "--variant"))
            options.input = require_value(index, argc, argv, argument, "--variant", "-V");
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (argument == "-R" || is_option(argument, "--reference"))
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        else if (argument == "-L" || is_option(argument, "--intervals") ||
                 is_option(argument, "--interval") || is_option(argument, "--region"))
            options.regions.push_back(require_value(
                index, argc, argv, argument,
                argument == "-L" ? "--intervals" :
                is_option(argument, "--region") ? "--region" :
                is_option(argument, "--intervals") ? "--intervals" : "--interval", "-L"));
        else if (is_option(argument, "--interval-set-rule")) {
            auto value = require_value(index, argc, argv, argument, "--interval-set-rule");
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
            if (value == "UNION") options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Union;
            else if (value == "INTERSECTION")
                options.interval_set_rule = fastgatk::io::HtsIntervalSetRule::Intersection;
            else throw std::invalid_argument("invalid --interval-set-rule: " + value);
        }
        else if (is_option(argument, "--filter-expression"))
            options.expressions.push_back(require_value(index, argc, argv, argument, "--filter-expression"));
        else if (is_option(argument, "--filter-name"))
            options.names.push_back(require_value(index, argc, argv, argument, "--filter-name"));
        else if (is_option(argument, "--genotype-filter-expression"))
            options.genotype_expressions.push_back(require_value(index, argc, argv, argument,
                                                                 "--genotype-filter-expression"));
        else if (is_option(argument, "--genotype-filter-name"))
            options.genotype_names.push_back(require_value(index, argc, argv, argument,
                                                            "--genotype-filter-name"));
        else if (is_option(argument, "--mask"))
            options.masks.push_back(require_value(index, argc, argv, argument, "--mask"));
        else if (is_option(argument, "--mask-name"))
            options.mask_name = require_value(index, argc, argv, argument, "--mask-name");
        else if (is_option(argument, "--mask-description"))
            options.mask_description = require_value(index, argc, argv, argument, "--mask-description");
        else if (is_option(argument, "--mask-extension")) {
            const auto value = require_value(index, argc, argv, argument, "--mask-extension");
            try { options.mask_extension = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --mask-extension: " + value); }
            if (options.mask_extension < 0) throw std::invalid_argument("--mask-extension must be non-negative");
        }
        else if (is_option(argument, "--cluster-size")) {
            const auto value = require_value(index, argc, argv, argument, "--cluster-size");
            try { options.cluster_size = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --cluster-size: " + value); }
            if (options.cluster_size < 0) throw std::invalid_argument("--cluster-size must be non-negative");
        }
        else if (is_option(argument, "--cluster-window-size")) {
            const auto value = require_value(index, argc, argv, argument, "--cluster-window-size");
            try { options.cluster_window_size = std::stoi(value); }
            catch (...) { throw std::invalid_argument("invalid --cluster-window-size: " + value); }
            if (options.cluster_window_size < 0)
                throw std::invalid_argument("--cluster-window-size must be non-negative");
        }
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--missing-values-evaluate-as-failing")
            options.missing_fails = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--missing-values-evaluate-as-failing");
        else if (argument.rfind("--missing-values-evaluate-as-failing=", 0) == 0)
            // Barclay binds this particular GATK argument as a separated
            // Boolean token; an embedded '=' is rejected as an option name.
            throw std::invalid_argument(
                "invalid option spelling for --missing-values-evaluate-as-failing; use a separated Boolean value");
        else if (argument == "--filter-not-in-mask") options.filter_records_not_in_mask = true;
        else if (argument == "--invert-filter-expression" || argument == "--invfilter")
            options.invert_filter_expression = true;
        else if (argument == "--invert-genotype-filter-expression" || argument == "--invG-filter")
            options.invert_genotype_filter_expression = true;
        else if (argument == "--apply-allele-specific-filters")
            options.apply_allele_specific_filters = true;
        else if (argument == "--invalidate-previous-filters" ||
                 argument.rfind("--invalidate-previous-filters=", 0) == 0)
            options.invalidate_previous_filters = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--invalidate-previous-filters");
        else if (argument == "--set-filtered-genotype-to-no-call")
            options.set_filtered_genotypes_to_nocall = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--set-filtered-genotype-to-no-call");
        else if (argument.rfind("--set-filtered-genotype-to-no-call=", 0) == 0)
            // Barclay rejects this Boolean's embedded '=' spelling; keep the
            // same separated-token boundary as GATK 4.6.2.0.
            throw std::invalid_argument(
                "invalid option spelling for --set-filtered-genotype-to-no-call; use a separated Boolean value");
        else if (argument == "--create-output-variant-index" ||
                 argument.rfind("--create-output-variant-index=", 0) == 0) {
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        } else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility flags.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" : "--verbosity");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-V/--variant is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.expressions.empty() && options.genotype_expressions.empty() && options.masks.empty() &&
        options.cluster_size == 0 && !options.set_filtered_genotypes_to_nocall)
        throw std::invalid_argument("at least one site/genotype filter expression or --mask is required");
    if (options.expressions.size() != options.names.size())
        throw std::invalid_argument("each --filter-expression requires one --filter-name");
    if (options.genotype_expressions.size() != options.genotype_names.size())
        throw std::invalid_argument("each --genotype-filter-expression requires one --genotype-filter-name");
    if ((options.cluster_size == 0) != (options.cluster_window_size == 0))
        throw std::invalid_argument("--cluster-size and --cluster-window-size must be provided together");
    if (options.cluster_size == 1)
        throw std::invalid_argument("--cluster-size must be at least 2 when enabled");
    if (options.filter_records_not_in_mask && options.masks.empty())
        throw std::invalid_argument("--filter-not-in-mask requires --mask");
    if (!options.mask_description.empty() && options.masks.empty())
        throw std::invalid_argument("--mask-description requires --mask");
    return options;
}

#if FASTGATK_HAS_HTSLIB

struct NumericExpr;

struct Expr {
    enum class Kind { Predicate, And, Or, Not } kind = Kind::Predicate;
    std::string field;
    std::string sample;
    std::string method;
    std::string op;
    double threshold = 0.0;
    bool string_compare = false;
    // JEXL's =~/!~ operators use the same string operand surface as ==/!=,
    // but evaluate the RHS as a regular expression.  Keep this explicit so
    // an expression such as vc.getID() !~ "decoy.*" cannot be mistaken for
    // a numeric comparison.
    bool regex_compare = false;
    bool boolean_rhs = false;
    std::string string_threshold;
    bool string_method = false;
    std::string string_method_name;
    std::string string_method_argument;
    bool null_compare = false;
    // Fixed index for an INFO vector expression such as AF[0] or
    // vc.getAttribute("AF")[0].  Keep this separate from allele_index,
    // which is the dynamic alternate-allele selector used by AS_* filters.
    int field_index = -1;
    int allele_index = -1;
    bool reference_allele = false;
    std::shared_ptr<NumericExpr> numeric_left;
    std::shared_ptr<NumericExpr> numeric_right;
    std::shared_ptr<Expr> left;
    std::shared_ptr<Expr> right;
};

// Numeric JEXL operands are kept as a small Host-side AST.  The VCF record
// remains owned by HTSlib and all field loads still go through read_field;
// this node only adds deterministic arithmetic around an existing scalar
// field/INFO element.  Keeping the tree separate from Expr avoids changing
// the boolean/regex/missing-value semantics already covered by the oracle.
struct NumericExpr {
    enum class Kind { Literal, Field, Add, Subtract, Multiply, Divide, Modulo, Negate };
    Kind kind = Kind::Field;
    double literal = 0.0;
    std::string field;
    std::shared_ptr<NumericExpr> left;
    std::shared_ptr<NumericExpr> right;
};

struct EvalResult {
    bool value = false;
    bool present = true;
};

std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool wrapped_by_parentheses(const std::string& value) {
    if (value.size() < 2 || value.front() != '(' || value.back() != ')') return false;
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const char character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') ++depth;
        else if (character == ')' && --depth == 0 && index + 1 != value.size()) return false;
    }
    return depth == 0;
}

std::size_t find_top_level(const std::string& value, const std::string& token) {
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    for (std::size_t index = 0; index + token.size() <= value.size(); ++index) {
        const char character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') { ++depth; continue; }
        if (character == ')') { --depth; continue; }
        if (depth == 0 && value.compare(index, token.size(), token) == 0) return index;
    }
    return std::string::npos;
}

// Return the right-most top-level arithmetic operator.  Selecting the
// right-most operator gives the usual left-associative parse for operators at
// the same precedence (for example `DP - 1 - 2`) while the precedence split
// below keeps multiplication/division tighter than addition/subtraction.
std::pair<std::size_t, char> find_top_level_arithmetic(const std::string& value) {
    int depth = 0;
    bool quoted = false;
    char quote = 0;
    std::size_t additive = std::string::npos;
    std::size_t multiplicative = std::string::npos;
    char additive_op = 0;
    char multiplicative_op = 0;
    for (std::size_t index = 0; index < value.size(); ++index) {
        const char character = value[index];
        if ((character == '\'' || character == '"') && (!quoted || character == quote)) {
            if (quoted) quoted = false;
            else { quoted = true; quote = character; }
            continue;
        }
        if (quoted) continue;
        if (character == '(') { ++depth; continue; }
        if (character == ')') { --depth; continue; }
        if (depth != 0 || (character != '+' && character != '-' &&
                           character != '*' && character != '/' && character != '%'))
            continue;
        // A leading sign (or a sign immediately following another operator
        // or an opening parenthesis) is unary and belongs to the operand.
        if ((character == '+' || character == '-') &&
            (index == 0 || value[index - 1] == '(' || value[index - 1] == '+' ||
             value[index - 1] == '-' || value[index - 1] == '*' ||
             value[index - 1] == '/' || value[index - 1] == '%' ||
             ((value[index - 1] == 'e' || value[index - 1] == 'E') && index > 1 &&
              std::isdigit(static_cast<unsigned char>(value[index - 2])))))
            continue;
        if (character == '+' || character == '-') {
            additive = index;
            additive_op = character;
        } else {
            multiplicative = index;
            multiplicative_op = character;
        }
    }
    if (additive != std::string::npos) return {additive, additive_op};
    if (multiplicative != std::string::npos) return {multiplicative, multiplicative_op};
    return {std::string::npos, 0};
}

std::shared_ptr<NumericExpr> parse_numeric_expression(const std::string& source) {
    auto value = trim_copy(source);
    if (value.empty()) throw std::invalid_argument("empty numeric expression");
    while (wrapped_by_parentheses(value)) value = trim_copy(value.substr(1, value.size() - 2));
    const auto split = find_top_level_arithmetic(value);
    if (split.first != std::string::npos) {
        auto node = std::make_shared<NumericExpr>();
        if (split.second == '+') node->kind = NumericExpr::Kind::Add;
        else if (split.second == '-') node->kind = NumericExpr::Kind::Subtract;
        else if (split.second == '*') node->kind = NumericExpr::Kind::Multiply;
        else if (split.second == '/') node->kind = NumericExpr::Kind::Divide;
        else node->kind = NumericExpr::Kind::Modulo;
        node->left = parse_numeric_expression(value.substr(0, split.first));
        node->right = parse_numeric_expression(value.substr(split.first + 1));
        return node;
    }
    if (!value.empty() && value.front() == '-') {
        auto node = std::make_shared<NumericExpr>();
        node->kind = NumericExpr::Kind::Negate;
        node->left = parse_numeric_expression(value.substr(1));
        return node;
    }
    std::size_t consumed = 0;
    try {
        const auto literal = std::stod(value, &consumed);
        if (consumed == value.size()) {
            auto node = std::make_shared<NumericExpr>();
            node->kind = NumericExpr::Kind::Literal;
            node->literal = literal;
            return node;
        }
    } catch (...) {
        // The operand is a field/method expression; read_field will apply the
        // normal HTSlib missing-value policy at evaluation time.
    }
    auto node = std::make_shared<NumericExpr>();
    node->kind = NumericExpr::Kind::Field;
    node->field = value;
    return node;
}

bool numeric_expression_has_operator(const NumericExpr& expression) {
    return expression.kind != NumericExpr::Kind::Literal &&
           expression.kind != NumericExpr::Kind::Field;
}

std::string unquote(std::string value) {
    value = trim_copy(std::move(value));
    if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                              (value.front() == '\'' && value.back() == '\'')))
        return value.substr(1, value.size() - 2);
    return value;
}

bool parse_call(const std::string& raw, std::string& name, std::string& argument) {
    const auto value = trim_copy(raw);
    const auto open = value.find('(');
    if (open == std::string::npos || value.empty() || value.back() != ')' || open == 0) return false;
    name = trim_copy(value.substr(0, open));
    argument = unquote(value.substr(open + 1, value.size() - open - 2));
    return true;
}

bool parse_string_method(const std::string& raw, std::string& field,
                         std::string& method, std::string& argument) {
    const auto value = trim_copy(raw);
    struct Method { const char* token; const char* name; };
    constexpr Method methods[] = {
        {".contains(", "contains"}, {".startsWith(", "startsWith"},
        {".endsWith(", "endsWith"}, {".matches(", "matches"},
    };
    for (const auto& candidate : methods) {
        const auto position = value.rfind(candidate.token);
        if (position == std::string::npos || value.empty() || value.back() != ')') continue;
        const auto argument_begin = position + std::strlen(candidate.token);
        if (argument_begin >= value.size() - 1) return false;
        field = trim_copy(value.substr(0, position));
        method = candidate.name;
        argument = unquote(value.substr(argument_begin, value.size() - argument_begin - 1));
        return !field.empty();
    }
    return false;
}

// Lower the common HTSJDK Allele accessors used by GATK JEXL.  Keeping the
// allele selection explicit in the expression node is important for
// multiallelic records: getAlternateAllele(i) is zero-based over ALT alleles,
// while the VCF record stores the reference at index zero.
bool parse_allele_method(const std::string& raw, bool& reference_allele,
                         int& allele_index, std::string& method) {
    const auto value = trim_copy(raw);
    const std::string reference_prefixes[] = {"vc.getReference().", "getReference()."};
    for (const auto& prefix : reference_prefixes) {
        if (value.rfind(prefix, 0) != 0) continue;
        const auto suffix = value.substr(prefix.size());
        if (suffix == "isSymbolic()" || suffix == "isReference()" ||
            suffix == "isNoCall()" || suffix == "isCalled()" ||
            suffix == "isBreakpoint()" || suffix == "isSingleBreakend()" ||
            suffix == "isNonReference()" || suffix == "isNonRefAllele()" ||
            suffix == "length()") {
            reference_allele = true;
            allele_index = -1;
            method = suffix;
            return true;
        }
    }
    const std::string alternate_prefix = "vc.getAlternateAllele(";
    const std::string alternate_prefix_short = "getAlternateAllele(";
    const std::string* selected_prefix = nullptr;
    if (value.rfind(alternate_prefix, 0) == 0) selected_prefix = &alternate_prefix;
    else if (value.rfind(alternate_prefix_short, 0) == 0) selected_prefix = &alternate_prefix_short;
    if (selected_prefix == nullptr) return false;
    const auto close = value.find(").", selected_prefix->size());
    if (close == std::string::npos || close + 2 >= value.size()) return false;
    const auto index_text = value.substr(selected_prefix->size(), close - selected_prefix->size());
    try {
        std::size_t consumed = 0;
        const auto index = std::stoi(index_text, &consumed);
        if (consumed != index_text.size() || index < 0) return false;
        const auto suffix = value.substr(close + 2);
        if (suffix != "isSymbolic()" && suffix != "isReference()" &&
            suffix != "isNoCall()" && suffix != "isCalled()" &&
            suffix != "isBreakpoint()" && suffix != "isSingleBreakend()" &&
            suffix != "isNonReference()" && suffix != "isNonRefAllele()" &&
            suffix != "length()") return false;
        reference_allele = false;
        allele_index = index;
        method = suffix;
        return true;
    } catch (...) {
        return false;
    }
}

bool indexed_format_method(const std::string& method, const char* name, int& index) {
    const std::string prefix = std::string(name) + "()[";
    if (method.rfind(prefix, 0) == 0 && method.size() > prefix.size() && method.back() == ']') {
        try {
            std::size_t consumed = 0;
            index = std::stoi(method.substr(prefix.size(), method.size() - prefix.size() - 1), &consumed);
            return consumed == method.size() - prefix.size() - 1 && index >= 0;
        } catch (...) {
            return false;
        }
    }
    // HTSJDK/JEXL also exposes collection elements with a dot suffix, e.g.
    // `vc.getGenotype("S").getAD().1`.  Keep the zero-based index explicit
    // and share the same reader with bracket notation.
    const std::string dot_prefix = std::string(name) + "().";
    if (method.rfind(dot_prefix, 0) != 0 || method.size() <= dot_prefix.size()) return false;
    try {
        std::size_t consumed = 0;
        index = std::stoi(method.substr(dot_prefix.size()), &consumed);
        return consumed == method.size() - dot_prefix.size() && index >= 0;
    } catch (...) {
        return false;
    }
}

// Parse the common HTSJDK/JEXL indexed INFO spellings.  The base expression
// is returned without the bracket suffix so the existing HTSlib readers can
// fetch the vector and select the requested element deterministically.
bool parse_indexed_info_field(const std::string& raw, std::string& base, int& index) {
    const auto value = trim_copy(raw);
    // Explicit genotype FORMAT accessors (for example
    // vc.getGenotype("S1").getAD()[1]) are lowered by the GENOTYPE reader;
    // do not mistake their bracket suffix for an INFO-vector index.
    if (value.find("getGenotype(") != std::string::npos) return false;
    if (value.size() < 4) return false;
    std::string index_text;
    std::size_t base_end = std::string::npos;
    if (value.back() == ']') {
        const auto open = value.rfind('[');
        if (open == std::string::npos || open == 0 || open + 1 >= value.size() - 1) return false;
        index_text = value.substr(open + 1, value.size() - open - 2);
        base_end = open;
    } else if (value.back() == ')') {
        // HTSJDK's JEXL collection accessor is vc.getAttribute("TAG").get(i).
        const auto token = value.rfind(".get(");
        if (token == std::string::npos || token == 0 || token + 5 >= value.size() - 1) return false;
        index_text = value.substr(token + 5, value.size() - token - 6);
        base_end = token;
    } else {
        // JEXL's list accessor shorthand is `.N` (for example
        // `vc.getAttribute("AF").1`).  Only accept an all-digit suffix so
        // ordinary method calls and dotted field names are not lowered.
        const auto dot = value.rfind('.');
        if (dot == std::string::npos || dot == 0 || dot + 1 >= value.size()) return false;
        index_text = value.substr(dot + 1);
        if (!std::all_of(index_text.begin(), index_text.end(),
                         [](unsigned char character) { return std::isdigit(character) != 0; }))
            return false;
        base_end = dot;
    }
    try {
        std::size_t consumed = 0;
        const auto parsed = std::stoi(index_text, &consumed);
        if (consumed != index_text.size() || parsed < 0) return false;
        base = trim_copy(value.substr(0, base_end));
        if (base.empty()) return false;
        index = parsed;
        return true;
    } catch (...) {
        return false;
    }
}

bool allele_is_breakpoint(const std::string& allele) {
    return allele.size() > 1 && allele.find_first_of("[]") != std::string::npos;
}

bool allele_is_single_breakend(const std::string& allele) {
    return allele.size() > 1 && (allele.front() == '.' || allele.back() == '.');
}

bool allele_is_symbolic(const std::string& allele) {
    // HTSJDK treats the VCF spanning-deletion allele `*` as symbolic.  It is
    // one character long, so keep it outside the structural-allele length
    // guard used for breakends and symbolic `<...>` alleles.
    return allele == "*" || (allele.size() > 1 &&
        (allele.front() == '<' || allele.back() == '>' ||
         allele_is_breakpoint(allele) || allele_is_single_breakend(allele)));
}

std::string variant_type(bcf1_t* record) {
    const auto types = bcf_get_variant_types(record);
    const bool snp = (types & VCF_SNP) != 0;
    const bool indel = (types & VCF_INDEL) != 0;
    const bool mnp = (types & VCF_MNP) != 0;
    if (types == VCF_REF) return "NO_VARIATION";
    if (snp && indel) return "MIXED";
    if (snp) return "SNP";
    if (indel) return "INDEL";
    if (mnp) return "MNP";
    for (int allele = 1; allele < record->n_allele; ++allele)
        if (record->d.allele[allele] != nullptr &&
            allele_is_symbolic(record->d.allele[allele])) return "SYMBOLIC";
    return "OTHER";
}

bool transition_state(bcf1_t* record, bool want_transition) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2 || record->d.allele == nullptr || record->d.allele[0] == nullptr ||
        std::strlen(record->d.allele[0]) != 1)
        return false;
    const char reference = static_cast<char>(std::toupper(record->d.allele[0][0]));
    bool seen = false;
    for (int index = 1; index < record->n_allele; ++index) {
        const char* alternate = record->d.allele[index];
        if (alternate == nullptr || std::strlen(alternate) != 1) return false;
        const char allele = static_cast<char>(std::toupper(alternate[0]));
        const bool transition = (reference == 'A' && allele == 'G') ||
                                (reference == 'G' && allele == 'A') ||
                                (reference == 'C' && allele == 'T') ||
                                (reference == 'T' && allele == 'C');
        if (transition != want_transition) return false;
        seen = true;
    }
    return seen;
}

bool filtered_record(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    if (record->d.n_flt == 0) return false;
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && std::strcmp(name, "PASS") != 0 && std::strcmp(name, ".") != 0) return true;
    }
    return false;
}

bool filter_set_contains(const bcf_hdr_t* header, bcf1_t* record, const std::string& wanted) {
    if (wanted.empty()) return false;
    bcf_unpack(record, BCF_UN_FLT);
    for (int index = 0; index < record->d.n_flt; ++index) {
        const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
        if (name != nullptr && wanted == name) return true;
    }
    return false;
}

void sort_filter_labels(const bcf_hdr_t* header, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_FLT);
    if (record->d.n_flt <= 1) return;
    std::vector<int> ids(record->d.flt, record->d.flt + record->d.n_flt);
    std::stable_sort(ids.begin(), ids.end(), [&](int left, int right) {
        const auto* left_name = bcf_hdr_int2id(header, BCF_DT_ID, left);
        const auto* right_name = bcf_hdr_int2id(header, BCF_DT_ID, right);
        if (left_name == nullptr || right_name == nullptr)
            return left < right;
        if (std::strcmp(left_name, right_name) != 0)
            return std::strcmp(left_name, right_name) < 0;
        return left < right;
    });
    if (std::equal(ids.begin(), ids.end(), record->d.flt)) return;
    if (bcf_update_filter(header, record, ids.data(), static_cast<int>(ids.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sort FILTER labels");
}

// Match VariantContext's sample-aware polymorphism semantics: only called
// non-reference alleles make a record polymorphic.  This deliberately ignores
// unused ALT alleles and no-call entries, which is important for filtered
// gVCF/reference-confidence records.
bool polymorphic_in_samples(const bcf_hdr_t* header, bcf1_t* record) {
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || record->n_sample <= 0 || genotype_count <= 0 ||
        genotype_count % record->n_sample != 0) {
        free(genotypes);
        return false;
    }
    const int ploidy = genotype_count / record->n_sample;
    bool polymorphic = false;
    for (int sample = 0; sample < record->n_sample && !polymorphic; ++sample) {
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto encoded = genotypes[sample * ploidy + allele];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) continue;
            if (bcf_gt_allele(encoded) > 0) {
                polymorphic = true;
                break;
            }
        }
    }
    free(genotypes);
    return polymorphic;
}

struct FiltrationGenotypeCounts {
    int called_chromosomes = 0;
    int no_call_count = 0;
    int hom_ref_count = 0;
    int het_count = 0;
    int hom_var_count = 0;
    bool has_genotypes = false;
};

std::string canonical_genotype_filter_status(std::string status) {
    if (status.empty() || status == "." || status == "PASS") return "PASS";
    std::vector<std::string> labels;
    std::size_t begin = 0;
    while (begin <= status.size()) {
        const auto end = status.find(';', begin);
        const auto label = status.substr(begin, end == std::string::npos
                                                   ? std::string::npos : end - begin);
        if (!label.empty() && label != "." && label != "PASS") labels.push_back(label);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    std::sort(labels.begin(), labels.end());
    labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
    if (labels.empty()) return "PASS";
    std::ostringstream output;
    for (std::size_t index = 0; index < labels.size(); ++index) {
        if (index != 0) output << ';';
        output << labels[index];
    }
    return output.str();
}

bool genotype_counts(const bcf_hdr_t* header, bcf1_t* record,
                     FiltrationGenotypeCounts& counts) {
    int32_t* genotypes = nullptr;
    int genotype_count = 0;
    const auto length = bcf_get_genotypes(header, record, &genotypes, &genotype_count);
    if (length <= 0 || record->n_sample <= 0 || genotype_count <= 0 ||
        genotype_count % record->n_sample != 0) {
        free(genotypes);
        return false;
    }
    counts.has_genotypes = true;
    const int ploidy = genotype_count / record->n_sample;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        bool called = true;
        bool any_alt = false;
        bool all_same = true;
        int first_allele = -1;
        for (int copy = 0; copy < ploidy; ++copy) {
            const auto encoded = genotypes[sample * ploidy + copy];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                called = false;
                continue;
            }
            ++counts.called_chromosomes;
            const int allele = bcf_gt_allele(encoded);
            if (first_allele < 0) first_allele = allele;
            else if (allele != first_allele) all_same = false;
            if (allele > 0) any_alt = true;
        }
        if (!called) ++counts.no_call_count;
        else if (all_same && first_allele == 0) ++counts.hom_ref_count;
        else if (all_same && any_alt) ++counts.hom_var_count;
        else ++counts.het_count;
    }
    free(genotypes);
    return true;
}

std::shared_ptr<Expr> parse_expression(const std::string& source, bool genotype_context = false) {
    auto value = trim_copy(source);
    while (wrapped_by_parentheses(value)) value = trim_copy(value.substr(1, value.size() - 2));
    if (value.empty()) throw std::invalid_argument("empty filter expression");
    const auto or_position = find_top_level(value, "||");
    if (or_position != std::string::npos) {
        auto node = std::make_shared<Expr>();
        node->kind = Expr::Kind::Or;
        node->left = parse_expression(value.substr(0, or_position), genotype_context);
        node->right = parse_expression(value.substr(or_position + 2), genotype_context);
        return node;
    }
    const auto and_position = find_top_level(value, "&&");
    if (and_position != std::string::npos) {
        auto node = std::make_shared<Expr>();
        node->kind = Expr::Kind::And;
        node->left = parse_expression(value.substr(0, and_position), genotype_context);
        node->right = parse_expression(value.substr(and_position + 2), genotype_context);
        return node;
    }
    if (value.front() == '!') {
        auto node = std::make_shared<Expr>();
        node->kind = Expr::Kind::Not;
        node->left = parse_expression(value.substr(1), genotype_context);
        return node;
    }

    auto node = std::make_shared<Expr>();
    node->kind = Expr::Kind::Predicate;
    const std::vector<std::string> operators{"<=", ">=", "==", "!=", "=~", "!~", "<", ">"};
    std::size_t operator_position = std::string::npos;
    std::string selected_operator;
    for (const auto& op : operators) {
        const auto position = find_top_level(value, op);
        if (position != std::string::npos && (operator_position == std::string::npos || position < operator_position)) {
            operator_position = position;
            selected_operator = op;
        }
    }
    std::string string_method_field;
    std::string string_method_name;
    std::string string_method_argument;
    if (parse_string_method(value, string_method_field, string_method_name, string_method_argument)) {
        node->field = "STRING_METHOD";
        node->string_method = true;
        node->string_method_name = std::move(string_method_name);
        node->string_method_argument = std::move(string_method_argument);
        node->sample = std::move(string_method_field);
    } else if (operator_position != std::string::npos) {
        node->field = trim_copy(value.substr(0, operator_position));
        const auto rhs = trim_copy(value.substr(operator_position + selected_operator.size()));
        if (rhs == "null") {
            if (selected_operator != "==" && selected_operator != "!=")
                throw std::invalid_argument("unsupported null filter expression: " + source);
            node->null_compare = true;
        } else if ((selected_operator == "==" || selected_operator == "!=") &&
                   rhs.rfind("VariantContext.Type.", 0) == 0) {
            node->string_compare = true;
            node->string_threshold = rhs.substr(std::strlen("VariantContext.Type."));
        } else if ((selected_operator == "==" || selected_operator == "!=") &&
                   (rhs == "true" || rhs == "false")) {
            // JEXL exposes Java boolean methods as first-class values.  The
            // evaluator represents booleans as deterministic 0/1 values for
            // comparison, while retaining this marker to reject nonsensical
            // ordering operators on a boolean RHS.
            node->boolean_rhs = true;
            node->threshold = rhs == "true" ? 1.0 : 0.0;
        } else if (selected_operator == "=~" || selected_operator == "!~") {
            if (rhs.empty())
                throw std::invalid_argument("empty regex filter expression: " + source);
            node->string_compare = true;
            node->regex_compare = true;
            node->string_threshold = unquote(rhs);
            // Accept both Java/JEXL quoted patterns and slash-delimited
            // regex literals.  The latter is useful when porting expressions
            // from command lines that use /pattern/ notation.
            if (node->string_threshold.size() >= 2 &&
                node->string_threshold.front() == '/' &&
                node->string_threshold.back() == '/')
                node->string_threshold = node->string_threshold.substr(
                    1, node->string_threshold.size() - 2);
        } else if (rhs.size() >= 2 && ((rhs.front() == '"' && rhs.back() == '"') ||
                                (rhs.front() == '\'' && rhs.back() == '\''))) {
            if (selected_operator != "==" && selected_operator != "!=")
                throw std::invalid_argument("unsupported string filter expression: " + source);
            node->string_compare = true;
            node->string_threshold = unquote(rhs);
        } else {
            // Keep simple numeric literals on the legacy scalar path, but
            // lower compound operands (for example `DP / 2` or
            // `QUAL + 5`) into the deterministic numeric AST.
            auto rhs_numeric = parse_numeric_expression(rhs);
            if (numeric_expression_has_operator(*rhs_numeric)) {
                node->numeric_right = std::move(rhs_numeric);
            } else {
                std::size_t consumed = 0;
                try { node->threshold = std::stod(rhs, &consumed); }
                catch (...) { throw std::invalid_argument("unsupported non-numeric filter expression: " + source); }
                if (consumed != rhs.size()) throw std::invalid_argument("unsupported filter expression: " + source);
            }
        }
        node->op = selected_operator;
    } else {
        bool reference_allele = false;
        int allele_index = -1;
        std::string allele_method;
        if (parse_allele_method(value, reference_allele, allele_index, allele_method)) {
            node->field = "ALLELE_METHOD";
            node->method = std::move(allele_method);
            node->reference_allele = reference_allele;
            node->allele_index = allele_index;
        } else {
        std::string method;
        std::string argument;
        if (parse_call(value, method, argument)) {
            if (method.rfind("vc.", 0) == 0) method.erase(0, 3);
            const std::set<std::string> site_methods{
                "isSNP", "isIndel", "isMNP", "isTransition", "isTransversion", "isMixed",
                "isSymbolic", "isVariant", "isBiallelic", "isMultiallelic", "isFiltered",
                "isPass", "isNotFiltered", "isNoVariation", "isPolymorphicInSamples",
                "isMonomorphicInSamples", "hasGenotypes", "hasAlternateAllele", "hasAttribute"
            };
            if (site_methods.count(method) != 0) {
                if (method == "hasAttribute" && argument.empty())
                    throw std::invalid_argument("hasAttribute requires a tag: " + source);
                node->field = "SITE_METHOD";
                node->method = method;
                node->sample = argument;
            } else {
                node->field = value;
            }
        } else {
            node->field = value;
        }
        }
    }

    // Compound arithmetic on the left-hand side must be recognized before
    // field canonicalization (which intentionally handles one field at a
    // time).  The evaluator marks this predicate as ARITHMETIC and resolves
    // each leaf through the same read_field/allele-index contract.
    if (node->kind == Expr::Kind::Predicate && node->numeric_left == nullptr &&
        node->field != "GENOTYPE" && node->field != "SITE_METHOD" &&
        node->field != "ALLELE_METHOD" && node->field != "STRING_METHOD") {
        auto lhs_numeric = parse_numeric_expression(node->field);
        if (numeric_expression_has_operator(*lhs_numeric)) {
            node->numeric_left = std::move(lhs_numeric);
            node->field = "ARITHMETIC";
        }
    }

    // INFO vectors use a fixed zero-based element index.  Do this before
    // canonicalizing vc.getAttribute("TAG") so both TAG[0] and the explicit
    // VariantContext spelling share the same HTSlib reader.  Genotype-context
    // AD[1]/PL[1] is intentionally handled below by the FORMAT reader.
    if (!genotype_context) {
        std::string base;
        int index = -1;
        if (parse_indexed_info_field(node->field, base, index)) {
            node->field = std::move(base);
            node->field_index = index;
        }
    }

    // Canonicalize the common GATK JEXL forms into a compact predicate.
    const std::string attribute_prefix = "vc.getAttribute(";
    if (node->field.rfind(attribute_prefix, 0) == 0 && node->field.back() == ')') {
        node->field = unquote(node->field.substr(attribute_prefix.size(), node->field.size() - attribute_prefix.size() - 1));
    }
    if (node->field == "vc.getStart()" || node->field == "getStart()") node->field = "START";
    if (node->field == "vc.getEnd()" || node->field == "getEnd()") node->field = "END";
    if (node->field == "vc.getType()" || node->field == "getType()") node->field = "TYPE";
    if (node->field == "vc.getNAlleles()" || node->field == "getNAlleles()") node->field = "NALLELES";
    if (node->field == "vc.getNSamples()" || node->field == "getNSamples()") node->field = "NSAMPLES";
    if (node->field == "vc.getCalledChrCount()" || node->field == "getCalledChrCount()") node->field = "CALLED_CHR_COUNT";
    if (node->field == "vc.getNoCallCount()" || node->field == "getNoCallCount()") node->field = "NO_CALL_COUNT";
    if (node->field == "vc.getHomRefCount()" || node->field == "getHomRefCount()") node->field = "HOM_REF_COUNT";
    if (node->field == "vc.getHetCount()" || node->field == "getHetCount()") node->field = "HET_COUNT";
    if (node->field == "vc.getHomVarCount()" || node->field == "getHomVarCount()") node->field = "HOM_VAR_COUNT";
    if (node->field == "vc.getAlleles().size()" || node->field == "getAlleles().size()")
        node->field = "ALLELES_SIZE";
    if (node->field == "vc.getFilters().size()" || node->field == "getFilters().size()")
        node->field = "FILTERS_SIZE";
    if (node->field == "vc.getGenotypes().size()" || node->field == "getGenotypes().size()")
        node->field = "GENOTYPES_SIZE";
    if (node->field == "vc.getFilters().isEmpty()" || node->field == "getFilters().isEmpty()") {
        node->field = "SITE_METHOD";
        node->method = "filtersIsEmpty";
    }
    if (node->field == "vc.getGenotypes().isEmpty()" || node->field == "getGenotypes().isEmpty()") {
        node->field = "SITE_METHOD";
        node->method = "genotypesIsEmpty";
    }
    // The no-comparison parser path already lowers vc.isSNP()/isPass() and
    // friends.  Do the same lowering when a caller compares the boolean
    // result explicitly (for example `vc.isSNP() == true`); GATK accepts
    // both spellings and they must share one evaluator implementation.
    if (node->field != "SITE_METHOD" && node->field != "GENOTYPE" &&
        node->field != "ALLELE_METHOD" && node->field != "STRING_METHOD") {
        std::string site_method;
        std::string site_argument;
        if (parse_call(node->field, site_method, site_argument)) {
            if (site_method.rfind("vc.", 0) == 0) site_method.erase(0, 3);
            const std::set<std::string> comparable_site_methods{
                "isSNP", "isIndel", "isMNP", "isTransition", "isTransversion", "isMixed",
                "isSymbolic", "isVariant", "isBiallelic", "isMultiallelic", "isFiltered",
                "isPass", "isNotFiltered", "isNoVariation", "isPolymorphicInSamples",
                "isMonomorphicInSamples", "hasGenotypes", "hasAlternateAllele", "hasAttribute"
            };
            if (comparable_site_methods.count(site_method) != 0 &&
                (site_method == "hasAlternateAllele" || site_method == "hasAttribute"
                     ? !site_argument.empty() : site_argument.empty())) {
                node->field = "SITE_METHOD";
                node->method = std::move(site_method);
                node->sample = std::move(site_argument);
            }
        }
    }
    const std::string genotype_prefix = "vc.getGenotype(";
    if (node->field.rfind(genotype_prefix, 0) == 0) {
        const auto close = node->field.find(')', genotype_prefix.size());
        if (close == std::string::npos) throw std::invalid_argument("malformed genotype expression: " + source);
        node->sample = unquote(node->field.substr(genotype_prefix.size(), close - genotype_prefix.size()));
        auto suffix = node->field.substr(close + 1);
        if (!suffix.empty() && suffix.front() == '.') suffix.erase(0, 1);
        node->method = suffix;
        node->field = "GENOTYPE";
    }
    // The same allele accessors can be used on the left-hand side of a
    // comparison, for example getAlternateAllele(0).length() > 1.
    if (node->field != "ALLELE_METHOD") {
        bool reference_allele = false;
        int allele_index = -1;
        std::string allele_method;
        if (parse_allele_method(node->field, reference_allele, allele_index, allele_method)) {
            node->field = "ALLELE_METHOD";
            node->method = std::move(allele_method);
            node->reference_allele = reference_allele;
            node->allele_index = allele_index;
        }
    }
    // The genotype-filter expression surface in GATK is intentionally
    // compact: `GQ < 20`, `DP == 0`, `AD[1] == 0`, and `isHet()` are all
    // evaluated in the context of each sample.  Lower these spellings to
    // the same GENOTYPE reader used by the explicit
    // vc.getGenotype("sample") form.  A later evaluator call supplies the
    // current sample index, so no sample name is required here.
    if (genotype_context && (node->field == "GQ" || node->field == "DP" || node->field == "MIN_DP")) {
        node->method = node->field == "GQ" ? "getGQ()" :
                       node->field == "MIN_DP" ? "getMIN_DP()" : "getDP()";
        node->field = "GENOTYPE";
    } else if (genotype_context &&
               (node->field.rfind("AD[", 0) == 0 || node->field.rfind("PL[", 0) == 0 ||
                node->field.rfind("AD.", 0) == 0 || node->field.rfind("PL.", 0) == 0)) {
        const bool bracket = node->field.front() == 'A' && node->field.size() > 2 && node->field[2] == '[';
        const auto tag = node->field.rfind("AD", 0) == 0 ? "getAD" : "getPL";
        std::string index_text;
        if (bracket) {
            const auto close = node->field.find(']');
            if (close == std::string::npos || close + 1 != node->field.size())
                throw std::invalid_argument("invalid compact genotype vector field: " + source);
            index_text = node->field.substr(3, close - 3);
        } else {
            index_text = node->field.substr(3);
        }
        const auto compact = std::string(tag) + (bracket ? "()[" : "().") + index_text + (bracket ? "]" : "");
        int ignored = -1;
        if (!indexed_format_method(compact, tag, ignored))
            throw std::invalid_argument("invalid compact genotype vector field: " + source);
        node->method = compact;
        node->field = "GENOTYPE";
    }
    if (genotype_context && node->field != "GENOTYPE" && node->field != "SITE_METHOD" &&
        node->field != "ALLELE_METHOD" && node->field != "STRING_METHOD" &&
        node->op.empty()) {
        std::string compact_method;
        std::string compact_argument;
        if (parse_call(node->field, compact_method, compact_argument)) {
            const std::set<std::string> genotype_methods{
                "isAvailable", "g.isAvailable", "isNoCall", "g.isNoCall",
                "isCalled", "g.isCalled", "isHet", "g.isHet", "isHom", "g.isHom",
                "isHomRef", "g.isHomRef", "isHomVar", "g.isHomVar",
                "hasDP", "g.hasDP", "hasGQ", "g.hasGQ", "hasAD", "g.hasAD",
                "hasPL", "g.hasPL"
            };
            if (compact_argument.empty() && genotype_methods.count(compact_method) != 0) {
                node->field = "GENOTYPE";
                node->method = compact_method + "()";
            }
        }
    }
    if (genotype_context && node->field != "GENOTYPE") {
        const std::set<std::string> compact_boolean_fields{
            "isAvailable", "isNoCall", "isCalled", "isHet", "isHom",
            "isHomRef", "isHomVar", "hasDP", "hasGQ", "hasAD", "hasPL"
        };
        if (compact_boolean_fields.count(node->field) != 0) {
            node->method = node->field + "()";
            node->field = "GENOTYPE";
        }
    }
    const std::vector<std::string> fields{"QUAL", "DP", "QD", "TLOD", "AF", "SOR", "FS", "MQ",
                                          "START", "END", "TYPE", "NALLELES", "NSAMPLES",
                                          "CALLED_CHR_COUNT", "NO_CALL_COUNT", "HOM_REF_COUNT",
                                          "HET_COUNT", "HOM_VAR_COUNT",
                                          "ALLELES_SIZE", "FILTERS_SIZE",
                                          "GENOTYPES_SIZE",
                                          "SITE_METHOD", "STRING_METHOD", "GENOTYPE", "ALLELE_METHOD",
                                          "ARITHMETIC"};
    const bool allele_specific_field = node->field.rfind("AS_", 0) == 0 && node->field.size() > 3;
    const bool indexed_info_field = node->field_index >= 0;
    if (std::find(fields.begin(), fields.end(), node->field) == fields.end() &&
        !node->null_compare && !node->string_compare && !allele_specific_field && !indexed_info_field)
        throw std::invalid_argument("unsupported filter field in expression: " + source);
    if (node->field == "GENOTYPE" && node->op.empty() &&
        node->method != "isHet()" && node->method != "isHom()" && node->method != "isHomRef()" &&
        node->method != "isHomVar()" && node->method != "isNoCall()" &&
        node->method != "isCalled()" && node->method != "isAvailable()" &&
        node->method != "hasDP()" && node->method != "hasGQ()" &&
        node->method != "hasAD()" && node->method != "hasPL()")
        throw std::invalid_argument("unsupported genotype predicate: " + source);
    if (node->field == "SITE_METHOD" && !node->op.empty() && !node->boolean_rhs)
        throw std::invalid_argument("site methods require a boolean RHS: " + source);
    if (node->field == "ALLELE_METHOD" && node->op.empty() &&
        node->method != "isSymbolic()" && node->method != "isReference()" &&
        node->method != "isNoCall()" && node->method != "isCalled()" &&
        node->method != "isBreakpoint()" && node->method != "isSingleBreakend()" &&
        node->method != "isNonReference()" && node->method != "isNonRefAllele()")
        throw std::invalid_argument("unsupported boolean allele method: " + source);
    if (node->field == "ALLELE_METHOD" && !node->op.empty() && node->method != "length()")
        throw std::invalid_argument("only allele length supports numeric comparison: " + source);
    if (node->string_compare && node->field == "GENOTYPE")
        throw std::invalid_argument("unsupported string genotype filter field: " + source);
    if (node->field == "GENOTYPE" && !node->op.empty()) {
        int ignored = -1;
        const bool indexed_ad = indexed_format_method(node->method, "getAD", ignored);
        const bool indexed_pl = indexed_format_method(node->method, "getPL", ignored);
        if (node->method != "getDP()" && node->method != "getGQ()" &&
            node->method != "getMIN_DP()" &&
            node->method != "getPhredScaledQual()" && node->method != "getPloidy()" &&
            node->method != "isAvailable()" && node->method != "isNoCall()" &&
            node->method != "isCalled()" && node->method != "isHet()" &&
            node->method != "isHom()" && node->method != "isHomRef()" &&
            node->method != "isHomVar()" && node->method != "hasDP()" &&
            node->method != "hasGQ()" && node->method != "hasAD()" &&
            node->method != "hasPL()" &&
            !indexed_ad && !indexed_pl)
            throw std::invalid_argument("unsupported genotype field: " + source);
        if ((node->method.rfind("getAD()[", 0) == 0 || node->method.rfind("getPL()[", 0) == 0) &&
            !indexed_ad && !indexed_pl)
            throw std::invalid_argument("invalid indexed genotype field: " + source);
    }
    return node;
}

bool uses_allele_specific_field(const Expr& expression) {
    if (expression.kind != Expr::Kind::Predicate) {
        return (expression.left && uses_allele_specific_field(*expression.left)) ||
               (expression.right && uses_allele_specific_field(*expression.right));
    }
    // A bare AS_* operand is evaluated once per alternate allele.  An
    // explicitly indexed AS_* vector element is a scalar site predicate and
    // must not unexpectedly require --apply-allele-specific-filters.
    return expression.field_index < 0 && expression.field.rfind("AS_", 0) == 0;
}

bool compare(double value, const std::string& op, double threshold) {
    if (op == "<") return value < threshold;
    if (op == "<=") return value <= threshold;
    if (op == ">") return value > threshold;
    if (op == ">=") return value >= threshold;
    if (op == "==") return value == threshold;
    if (op == "!=") return value != threshold;
    return false;
}

bool read_field(const bcf_hdr_t* header, bcf1_t* record, const std::string& field,
                double& value, int allele_index = -1, int field_index = -1) {
    if (field == "QUAL") {
        if (bcf_float_is_missing(record->qual) || bcf_float_is_vector_end(record->qual)) return false;
        value = record->qual;
        return true;
    }
    if (field == "START") {
        if (record->pos < 0) return false;
        value = static_cast<double>(record->pos + 1);
        return true;
    }
    if (field == "END") {
        int32_t* end = nullptr;
        int count = 0;
        const auto length = bcf_get_info_int32(header, record, "END", &end, &count);
        if (length > 0 && count > 0 && end[0] != bcf_int32_missing && end[0] != bcf_int32_vector_end) {
            value = end[0];
            free(end);
            return true;
        }
        free(end);
        bcf_unpack(record, BCF_UN_STR);
        if (record->pos < 0 || record->d.allele == nullptr || record->d.allele[0] == nullptr) return false;
        value = static_cast<double>(record->pos + std::strlen(record->d.allele[0]));
        return true;
    }
    if (field == "NALLELES") {
        value = static_cast<double>(record->n_allele);
        return true;
    }
    if (field == "ALLELES_SIZE") {
        value = static_cast<double>(record->n_allele);
        return true;
    }
    if (field == "FILTERS_SIZE") {
        bcf_unpack(record, BCF_UN_FLT);
        value = static_cast<double>(record->d.n_flt);
        return true;
    }
    if (field == "vc.getAlleles().size()" || field == "getAlleles().size()") {
        value = static_cast<double>(record->n_allele);
        return true;
    }
    if (field == "vc.getFilters().size()" || field == "getFilters().size()") {
        bcf_unpack(record, BCF_UN_FLT);
        value = static_cast<double>(record->d.n_flt);
        return true;
    }
    if (field == "vc.getGenotypes().size()" || field == "getGenotypes().size()") {
        value = static_cast<double>(record->n_sample);
        return true;
    }
    if (field == "GENOTYPES_SIZE") {
        value = static_cast<double>(record->n_sample);
        return true;
    }
    if (field == "NSAMPLES") {
        value = static_cast<double>(record->n_sample);
        return true;
    }
    if (field == "CALLED_CHR_COUNT" || field == "NO_CALL_COUNT" ||
        field == "HOM_REF_COUNT" || field == "HET_COUNT" || field == "HOM_VAR_COUNT") {
        FiltrationGenotypeCounts counts;
        if (!genotype_counts(header, record, counts)) return false;
        if (field == "CALLED_CHR_COUNT") value = counts.called_chromosomes;
        else if (field == "NO_CALL_COUNT") value = counts.no_call_count;
        else if (field == "HOM_REF_COUNT") value = counts.hom_ref_count;
        else if (field == "HET_COUNT") value = counts.het_count;
        else value = counts.hom_var_count;
        return true;
    }
    float* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_float(header, record, field.c_str(), &values, &count);
    const auto selected_float = field_index >= 0 ? field_index : (allele_index >= 0 ? allele_index : 0);
    if (length > selected_float && count > selected_float &&
        !bcf_float_is_missing(values[selected_float]) &&
        !bcf_float_is_vector_end(values[selected_float])) {
        value = values[selected_float];
        free(values);
        return std::isfinite(value);
    }
    free(values);
    int32_t* integer_values = nullptr;
    int integer_count = 0;
    const auto integer_length = bcf_get_info_int32(header, record, field.c_str(), &integer_values, &integer_count);
    const auto selected_integer = field_index >= 0 ? field_index : (allele_index >= 0 ? allele_index : 0);
    if (integer_length <= selected_integer || integer_count <= selected_integer ||
        integer_values[selected_integer] == bcf_int32_missing ||
        integer_values[selected_integer] == bcf_int32_vector_end) {
        free(integer_values);
        return false;
    }
    value = integer_values[selected_integer];
    free(integer_values);
    return true;
}

bool canonicalize_numeric_leaf(std::string raw, std::string& field, int& field_index) {
    field = trim_copy(std::move(raw));
    field_index = -1;
    std::string indexed_base;
    int indexed = -1;
    if (parse_indexed_info_field(field, indexed_base, indexed)) {
        field = std::move(indexed_base);
        field_index = indexed;
    }
    const std::string attribute_prefix = "vc.getAttribute(";
    if (field.rfind(attribute_prefix, 0) == 0 && field.back() == ')')
        field = unquote(field.substr(attribute_prefix.size(), field.size() - attribute_prefix.size() - 1));
    if (field == "vc.getStart()" || field == "getStart()") field = "START";
    else if (field == "vc.getEnd()" || field == "getEnd()") field = "END";
    else if (field == "vc.getNAlleles()" || field == "getNAlleles()") field = "NALLELES";
    else if (field == "vc.getNSamples()" || field == "getNSamples()") field = "NSAMPLES";
    else if (field == "vc.getCalledChrCount()" || field == "getCalledChrCount()") field = "CALLED_CHR_COUNT";
    else if (field == "vc.getNoCallCount()" || field == "getNoCallCount()") field = "NO_CALL_COUNT";
    else if (field == "vc.getHomRefCount()" || field == "getHomRefCount()") field = "HOM_REF_COUNT";
    else if (field == "vc.getHetCount()" || field == "getHetCount()") field = "HET_COUNT";
    else if (field == "vc.getHomVarCount()" || field == "getHomVarCount()") field = "HOM_VAR_COUNT";
    else if (field == "vc.getAlleles().size()" || field == "getAlleles().size()") field = "ALLELES_SIZE";
    else if (field == "vc.getFilters().size()" || field == "getFilters().size()") field = "FILTERS_SIZE";
    else if (field == "vc.getGenotypes().size()" || field == "getGenotypes().size()") field = "GENOTYPES_SIZE";
    return !field.empty();
}

bool read_allele_method(bcf1_t* record, const Expr& expression,
                        double& value, bool& boolean_value);
bool read_genotype(const bcf_hdr_t* header, bcf1_t* record, const Expr& expression,
                   double& value, bool& boolean_value, int sample_override = -1);

bool evaluate_numeric_expression(const NumericExpr& expression,
                                const bcf_hdr_t* header, bcf1_t* record,
                                double& value, int allele_index = -1,
                                int sample_index = -1) {
    if (expression.kind == NumericExpr::Kind::Literal) {
        value = expression.literal;
        return true;
    }
    if (expression.kind == NumericExpr::Kind::Field) {
        // Compact genotype fields are valid inside a genotype-filter
        // arithmetic operand (for example `GQ + 1 >= 51`).  Resolve them via
        // the same FORMAT reader used by ordinary genotype predicates before
        // attempting INFO canonicalization.
        if (sample_index >= 0) {
            std::string method;
            if (expression.field == "GQ") method = "getGQ()";
            else if (expression.field == "DP") method = "getDP()";
            else if (expression.field == "MIN_DP") method = "getMIN_DP()";
            else if (expression.field.rfind("AD[", 0) == 0 || expression.field.rfind("AD.", 0) == 0 ||
                     expression.field.rfind("PL[", 0) == 0 || expression.field.rfind("PL.", 0) == 0) {
                const bool ad = expression.field.rfind("AD", 0) == 0;
                const auto tag = ad ? "getAD" : "getPL";
                const bool bracket = expression.field.size() > 3 && expression.field[2] == '[';
                std::string index_text;
                if (bracket && expression.field.back() == ']')
                    index_text = expression.field.substr(3, expression.field.size() - 4);
                else if (!bracket)
                    index_text = expression.field.substr(3);
                try {
                    std::size_t consumed = 0;
                    const auto index = std::stoi(index_text, &consumed);
                    if (consumed == index_text.size() && index >= 0)
                        method = std::string(tag) + (bracket ? "()[" : "().") + index_text + (bracket ? "]" : "");
                } catch (...) {
                }
            }
            if (!method.empty()) {
                Expr genotype;
                genotype.field = "GENOTYPE";
                genotype.method = std::move(method);
                // read_genotype's scalar FORMAT accessors are reached from
                // the comparison path; a non-empty marker keeps the helper
                // from interpreting getGQ()/getDP() as a bare boolean
                // predicate before returning the numeric value.
                genotype.op = "==";
                bool ignored = false;
                return read_genotype(header, record, genotype, value, ignored, sample_index);
            }
        }
        std::string field;
        int field_index = -1;
        if (!canonicalize_numeric_leaf(expression.field, field, field_index)) return false;
        bool reference_allele = false;
        int selected_allele = -1;
        std::string method;
        if (parse_allele_method(field, reference_allele, selected_allele, method)) {
            Expr allele;
            allele.field = "ALLELE_METHOD";
            allele.method = method;
            allele.reference_allele = reference_allele;
            allele.allele_index = selected_allele;
            bool ignored = false;
            return read_allele_method(record, allele, value, ignored);
        }
        return read_field(header, record, field, value, allele_index, field_index);
    }
    if (expression.kind == NumericExpr::Kind::Negate) {
        if (!expression.left || !evaluate_numeric_expression(*expression.left, header, record, value, allele_index, sample_index))
            return false;
        value = -value;
        return true;
    }
    if (!expression.left || !expression.right) return false;
    double left = 0.0;
    double right = 0.0;
    if (!evaluate_numeric_expression(*expression.left, header, record, left, allele_index, sample_index) ||
        !evaluate_numeric_expression(*expression.right, header, record, right, allele_index, sample_index))
        return false;
    switch (expression.kind) {
        case NumericExpr::Kind::Add: value = left + right; break;
        case NumericExpr::Kind::Subtract: value = left - right; break;
        case NumericExpr::Kind::Multiply: value = left * right; break;
        case NumericExpr::Kind::Divide: value = left / right; break;
        case NumericExpr::Kind::Modulo: value = std::fmod(left, right); break;
        default: return false;
    }
    return true;
}

bool read_string_field(const bcf_hdr_t* header, bcf1_t* record,
                       const std::string& field, std::string& value) {
    if (field == "TYPE") {
        value = variant_type(record);
        return true;
    }
    if (field == "vc.getID()" || field == "getID()") {
        bcf_unpack(record, BCF_UN_STR);
        if (record->d.id == nullptr || std::strcmp(record->d.id, ".") == 0) return false;
        value = record->d.id;
        return true;
    }
    if (field == "vc.getContig()" || field == "getContig()") {
        const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
        if (contig == nullptr) return false;
        value = contig;
        return true;
    }
    if (field == "vc.getFilters()" || field == "getFilters()") {
        bcf_unpack(record, BCF_UN_FLT);
        value.clear();
        for (int index = 0; index < record->d.n_flt; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
            if (name == nullptr || std::strcmp(name, "PASS") == 0 || std::strcmp(name, ".") == 0) continue;
            if (!value.empty()) value.push_back(';');
            value += name;
        }
        return true;
    }
    if (field == "vc.getReference().getBaseString()" || field == "getReference().getBaseString()") {
        bcf_unpack(record, BCF_UN_STR);
        if (record->d.allele == nullptr || record->n_allele == 0 || record->d.allele[0] == nullptr)
            return false;
        value = record->d.allele[0];
        return true;
    }
    const std::string alternate_prefix = "vc.getAlternateAllele(";
    const std::string alternate_suffix = ").getBaseString()";
    if (field.rfind(alternate_prefix, 0) == 0 &&
        field.size() > alternate_prefix.size() + alternate_suffix.size() &&
        field.compare(field.size() - alternate_suffix.size(), alternate_suffix.size(), alternate_suffix) == 0) {
        try {
            const auto index = std::stoi(field.substr(alternate_prefix.size(),
                field.size() - alternate_prefix.size() - alternate_suffix.size()));
            bcf_unpack(record, BCF_UN_STR);
            if (index < 0 || index + 1 >= record->n_allele || record->d.allele[index + 1] == nullptr)
                return false;
            value = record->d.allele[index + 1];
            return true;
        } catch (...) { return false; }
    }
    std::string attribute = field;
    const std::string prefix = "vc.getAttribute(";
    if (attribute.rfind(prefix, 0) == 0 && attribute.back() == ')')
        attribute = unquote(attribute.substr(prefix.size(), attribute.size() - prefix.size() - 1));
    char* text = nullptr;
    int text_size = 0;
    const auto length = bcf_get_info_string(header, record, attribute.c_str(), &text, &text_size);
    if (length < 0 || text == nullptr) {
        free(text);
        int32_t* integers = nullptr;
        int integer_count = 0;
        const auto integer_length = bcf_get_info_int32(header, record, attribute.c_str(), &integers, &integer_count);
        if (integer_length > 0 && integer_count > 0 && integers[0] != bcf_int32_missing &&
            integers[0] != bcf_int32_vector_end) {
            value = std::to_string(integers[0]);
            free(integers);
            return true;
        }
        free(integers);
        return false;
    }
    value.assign(text, static_cast<std::size_t>(length));
    free(text);
    return true;
}

bool read_allele_method(bcf1_t* record, const Expr& expression,
                        double& value, bool& boolean_value) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->d.allele == nullptr || record->n_allele <= 0) return false;
    const int index = expression.reference_allele ? 0 : expression.allele_index + 1;
    if (index < 0 || index >= record->n_allele || record->d.allele[index] == nullptr) return false;
    const std::string allele(record->d.allele[index]);
    if (expression.method == "length()") {
        const bool symbolic = allele_is_symbolic(allele);
        value = static_cast<double>(allele == "." || symbolic || allele == "*" ? 0 : allele.size());
        return allele != ".";
    }
    if (expression.method == "isSymbolic()") {
        boolean_value = allele_is_symbolic(allele);
        return true;
    }
    if (expression.method == "isBreakpoint()") {
        boolean_value = allele_is_breakpoint(allele);
        return true;
    }
    if (expression.method == "isSingleBreakend()") {
        boolean_value = allele_is_single_breakend(allele);
        return true;
    }
    if (expression.method == "isReference()") {
        boolean_value = expression.reference_allele;
        return true;
    }
    if (expression.method == "isNoCall()") {
        boolean_value = allele == ".";
        return true;
    }
    if (expression.method == "isCalled()") {
        boolean_value = allele != ".";
        return true;
    }
    if (expression.method == "isNonReference()") {
        boolean_value = !expression.reference_allele && allele != ".";
        return true;
    }
    if (expression.method == "isNonRefAllele()") {
        boolean_value = allele == "<NON_REF>";
        return true;
    }
    return false;
}

bool read_site_method(const bcf_hdr_t* header, bcf1_t* record, const Expr& expression,
                      bool& value) {
    const auto& method = expression.method;
    if (method == "isSNP") value = variant_type(record) == "SNP";
    else if (method == "isIndel") value = variant_type(record) == "INDEL";
    else if (method == "isMNP") value = variant_type(record) == "MNP";
    else if (method == "isTransition") value = transition_state(record, true);
    else if (method == "isTransversion") value = transition_state(record, false);
    else if (method == "isSymbolic") value = variant_type(record) == "SYMBOLIC";
    else if (method == "isVariant") value = record->n_allele > 1 && variant_type(record) != "NO_VARIATION";
    else if (method == "isBiallelic") value = record->n_allele == 2;
    else if (method == "isMultiallelic") value = record->n_allele > 2;
    else if (method == "isFiltered") value = filtered_record(header, record);
    else if (method == "isPass" || method == "isNotFiltered") value = !filtered_record(header, record);
    else if (method == "isNoVariation") value = variant_type(record) == "NO_VARIATION";
    else if (method == "isPolymorphicInSamples") value = polymorphic_in_samples(header, record);
    else if (method == "isMonomorphicInSamples") value = !polymorphic_in_samples(header, record);
    else if (method == "hasGenotypes") {
        FiltrationGenotypeCounts counts;
        value = genotype_counts(header, record, counts) && counts.has_genotypes;
    }
    else if (method == "hasAlternateAllele") {
        try {
            const auto index = std::stoi(expression.sample);
            value = index >= 0 && index + 1 < record->n_allele;
        } catch (...) { value = false; }
    }
    else if (method == "isMixed") {
        std::set<std::string> categories;
        bcf_unpack(record, BCF_UN_STR);
        if (record->n_allele > 1 && record->d.allele != nullptr && record->d.allele[0] != nullptr) {
            const auto reference_length = std::strlen(record->d.allele[0]);
            for (int index = 1; index < record->n_allele; ++index) {
                const char* alternate = record->d.allele[index];
                if (alternate == nullptr) continue;
                if (alternate[0] == '<' || alternate[0] == '*') categories.insert("SYMBOLIC");
                else if (std::strlen(alternate) == reference_length && reference_length > 1) categories.insert("MNP");
                else if (std::strlen(alternate) == reference_length) categories.insert("SNP");
                else categories.insert("INDEL");
            }
        }
        value = categories.size() > 1;
    } else if (method == "hasAttribute") {
        const auto id = bcf_hdr_id2int(header, BCF_DT_ID, expression.sample.c_str());
        value = false;
        if (id >= 0) {
            bcf_unpack(record, BCF_UN_INFO);
            for (int index = 0; index < record->n_info; ++index)
                if (record->d.info[index].key == id && record->d.info[index].vptr != nullptr) {
                    value = true;
                    break;
                }
        }
    } else if (method == "filtersIsEmpty") {
        bcf_unpack(record, BCF_UN_FLT);
        value = true;
        for (int index = 0; index < record->d.n_flt; ++index) {
            const auto* name = bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[index]);
            if (name != nullptr && std::strcmp(name, "PASS") != 0 && std::strcmp(name, ".") != 0) {
                value = false;
                break;
            }
        }
    } else if (method == "genotypesIsEmpty") {
        value = record->n_sample == 0;
    } else return false;
    return true;
}

bool read_genotype(const bcf_hdr_t* header, bcf1_t* record, const Expr& expression,
                   double& value, bool& boolean_value, int sample_override) {
    // GATK permits both an explicitly named genotype accessor
    // (vc.getGenotype("S1").getGQ()) and the compact per-genotype form
    // (GQ < 20, AD[1] == 0, isHet()).  The latter is evaluated once for each
    // sample by apply_genotype_filters; keep the sample selection at this
    // boundary instead of baking a sample into the parsed expression.
    const int sample_index = expression.sample.empty()
        ? sample_override
        : bcf_hdr_id2int(header, BCF_DT_SAMPLE, expression.sample.c_str());
    if (sample_index < 0 || sample_index >= bcf_hdr_nsamples(header)) return false;
    if (!expression.op.empty() &&
        (expression.method == "isAvailable()" || expression.method == "isNoCall()" ||
         expression.method == "isCalled()" || expression.method == "isHet()" ||
         expression.method == "isHom()" || expression.method == "isHomRef()" ||
         expression.method == "isHomVar()" || expression.method == "hasDP()" ||
         expression.method == "hasGQ()" || expression.method == "hasAD()" ||
         expression.method == "hasPL()")) {
        // The compact GATK spelling uses numeric predicates such as
        // `isHet == 1`.  Reuse the predicate reader and expose its boolean
        // result as the numeric 0/1 value expected by the comparison node.
        Expr predicate = expression;
        predicate.op.clear();
        predicate.string_compare = false;
        predicate.null_compare = false;
        double ignored = 0.0;
        bool predicate_value = false;
        if (!read_genotype(header, record, predicate, ignored, predicate_value, sample_override))
            return false;
        value = predicate_value ? 1.0 : 0.0;
        boolean_value = predicate_value;
        return true;
    }
    if (expression.op.empty()) {
        if (expression.method == "hasDP()" || expression.method == "hasGQ()" ||
            expression.method == "hasAD()" || expression.method == "hasPL()") {
            const char* tag = expression.method == "hasDP()" ? "DP" :
                              expression.method == "hasGQ()" ? "GQ" :
                              expression.method == "hasAD()" ? "AD" : "PL";
            int32_t* values = nullptr;
            int count = 0;
            const auto length = bcf_get_format_int32(header, record, tag, &values, &count);
            const int sample_count = bcf_hdr_nsamples(header);
            bool available = false;
            if (length > 0 && sample_count > 0 && count >= sample_count && count % sample_count == 0) {
                const int width = count / sample_count;
                if (width > 0) {
                    const auto observed = values[sample_index * width];
                    available = observed != bcf_int32_missing && observed != bcf_int32_vector_end;
                }
            }
            free(values);
            boolean_value = available;
            return true;
        }
        int32_t* genotypes = nullptr;
        int count = 0;
        const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
        const int sample_count = bcf_hdr_nsamples(header);
        if (length <= 0 || sample_count <= 0 || count < sample_count ||
            count % sample_count != 0) {
            free(genotypes);
            return false;
        }
        const int ploidy = count / sample_count;
        if (ploidy <= 0 || count <= sample_index * ploidy) {
            free(genotypes);
            return false;
        }
        if (expression.method == "isAvailable()" || expression.method == "isAvailable") {
            boolean_value = true;
            free(genotypes);
            return true;
        }
        bool no_call = false;
        bool all_same = true;
        bool all_ref = true;
        bool all_alt = true;
        int first_allele = -1;
        for (int allele_index = 0; allele_index < ploidy; ++allele_index) {
            const auto encoded = genotypes[sample_index * ploidy + allele_index];
            if (encoded == bcf_int32_vector_end || bcf_gt_is_missing(encoded)) {
                no_call = true;
                continue;
            }
            const int allele = bcf_gt_allele(encoded);
            if (allele_index == 0) first_allele = allele;
            else if (allele != first_allele) all_same = false;
            all_ref = all_ref && allele == 0;
            all_alt = all_alt && allele > 0;
        }
        if (expression.method == "isNoCall()") boolean_value = no_call;
        else if (expression.method == "isCalled()") boolean_value = !no_call;
        else if (no_call) { free(genotypes); return false; }
        else if (expression.method == "isHet()") boolean_value = !all_same;
        else if (expression.method == "isHom()") boolean_value = all_same;
        else if (expression.method == "isHomRef()") boolean_value = all_ref;
        else if (expression.method == "isHomVar()") boolean_value = all_same && all_alt;
        else { free(genotypes); return false; }
        free(genotypes);
        return true;
    }
    if (expression.method == "getDP()" || expression.method == "getMIN_DP()") {
        int32_t* values = nullptr;
        int count = 0;
        const auto* tag = expression.method == "getMIN_DP()" ? "MIN_DP" : "DP";
        const auto length = bcf_get_format_int32(header, record, tag, &values, &count);
        if (length <= sample_index) { free(values); return false; }
        if (values[sample_index] == bcf_int32_missing || values[sample_index] == bcf_int32_vector_end) { free(values); return false; }
        value = values[sample_index];
        free(values);
        return true;
    }
    if (expression.method == "getGQ()" || expression.method == "getPhredScaledQual()") {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(header, record, "GQ", &values, &count);
        if (length <= sample_index) { free(values); return false; }
        if (values[sample_index] == bcf_int32_missing || values[sample_index] == bcf_int32_vector_end) { free(values); return false; }
        value = values[sample_index];
        free(values);
        return true;
    }
    if (expression.method == "getPloidy()") {
        int32_t* genotypes = nullptr;
        int count = 0;
        const auto length = bcf_get_genotypes(header, record, &genotypes, &count);
        const int sample_count = bcf_hdr_nsamples(header);
        if (length <= 0 || sample_count <= 0 || count < sample_count || count % sample_count != 0) {
            free(genotypes);
            return false;
        }
        value = static_cast<double>(count / sample_count);
        free(genotypes);
        return true;
    }
    int requested_index = -1;
    std::string vector_method;
    if (indexed_format_method(expression.method, "getAD", requested_index)) vector_method = "AD";
    else if (indexed_format_method(expression.method, "getPL", requested_index)) vector_method = "PL";
    if (!vector_method.empty()) {
        int32_t* values = nullptr;
        int count = 0;
        const auto length = bcf_get_format_int32(header, record, vector_method.c_str(), &values, &count);
        const int sample_count = bcf_hdr_nsamples(header);
        if (length <= 0 || sample_count <= 0 || count < sample_count || count % sample_count != 0) {
            free(values);
            return false;
        }
        const int width = count / sample_count;
        if (requested_index < 0 || requested_index >= width) {
            free(values);
            return false;
        }
        const auto observed = values[sample_index * width + requested_index];
        if (observed == bcf_int32_missing || observed == bcf_int32_vector_end) {
            free(values);
            return false;
        }
        value = observed;
        free(values);
        return true;
    }
    return false;
}

EvalResult evaluate(const Expr& expression, const bcf_hdr_t* header, bcf1_t* record,
                    bool missing_fails, int allele_index = -1, int sample_index = -1) {
    if (expression.kind == Expr::Kind::And) {
        const auto left = evaluate(*expression.left, header, record, missing_fails, allele_index, sample_index);
        if (!left.value) return {false, left.present};
        const auto right = evaluate(*expression.right, header, record, missing_fails, allele_index, sample_index);
        return {right.value, left.present && right.present};
    }
    if (expression.kind == Expr::Kind::Or) {
        const auto left = evaluate(*expression.left, header, record, missing_fails, allele_index, sample_index);
        if (left.value) return {true, left.present};
        const auto right = evaluate(*expression.right, header, record, missing_fails, allele_index, sample_index);
        return {right.value, left.present && right.present};
    }
    if (expression.kind == Expr::Kind::Not) {
        const auto child = evaluate(*expression.left, header, record, missing_fails, allele_index, sample_index);
        return {!child.value, child.present};
    }
    bool present = false;
    bool boolean_value = false;
    double value = 0.0;
    double threshold = expression.threshold;
    if (expression.numeric_left != nullptr) {
        present = evaluate_numeric_expression(*expression.numeric_left, header, record,
                                              value, allele_index, sample_index);
    } else if (expression.field == "GENOTYPE")
        present = read_genotype(header, record, expression, value, boolean_value, sample_index);
    else if (expression.field == "ALLELE_METHOD")
        present = read_allele_method(record, expression, value, boolean_value);
    else if (expression.field == "SITE_METHOD") {
        present = read_site_method(header, record, expression, boolean_value);
        // Site methods are Java booleans.  When explicitly compared to a
        // boolean literal, expose the same 0/1 representation used by the
        // genotype predicate path so the ordinary comparator remains the
        // single source of ordering/equality semantics.
        if (present && !expression.op.empty()) value = boolean_value ? 1.0 : 0.0;
    }
    else if (expression.field == "STRING_METHOD") {
        std::string observed;
        present = read_string_field(header, record, expression.sample, observed);
        if (present) {
            const auto& argument = expression.string_method_argument;
            if (expression.string_method_name == "contains" &&
                (expression.sample == "vc.getFilters()" || expression.sample == "getFilters()"))
                boolean_value = filter_set_contains(header, record, argument);
            else if (expression.string_method_name == "contains") boolean_value = observed.find(argument) != std::string::npos;
            else if (expression.string_method_name == "startsWith") boolean_value = observed.rfind(argument, 0) == 0;
            else if (expression.string_method_name == "endsWith")
                boolean_value = observed.size() >= argument.size() &&
                    observed.compare(observed.size() - argument.size(), argument.size(), argument) == 0;
            else if (expression.string_method_name == "matches") {
                try { boolean_value = std::regex_search(observed, std::regex(argument)); }
                catch (const std::regex_error&) {
                    throw std::invalid_argument("invalid regex in VariantFiltration expression");
                }
            }
        }
    }
    else if (expression.string_compare) {
        std::string observed;
        present = read_string_field(header, record, expression.field, observed);
        if (present) {
            if (expression.regex_compare) {
                try {
                    const bool matched = std::regex_search(
                        observed, std::regex(expression.string_threshold));
                    boolean_value = expression.op == "=~" ? matched : !matched;
                } catch (const std::regex_error&) {
                    throw std::invalid_argument("invalid regex in VariantFiltration expression");
                }
            } else {
                boolean_value = expression.op == "==" ? observed == expression.string_threshold
                                                        : observed != expression.string_threshold;
            }
        }
    } else if (expression.numeric_left == nullptr) {
        present = read_field(header, record, expression.field, value, allele_index,
                             expression.field_index);
    }
    if (expression.null_compare) {
        // An indexed numeric element is absent when the requested vector
        // slot cannot be read; do not fall back to element zero and thereby
        // turn AF[1] == null into a false positive.
        if (!present && expression.field_index >= 0)
            return {expression.op == "==", true};
        if (!present) {
            std::string observed;
            present = read_string_field(header, record, expression.field, observed);
        }
        return {expression.op == "==" ? !present : present, true};
    }
    // JEXL coerces a missing String operand to null for regex matching; in
    // GATK, `missing !~ "pattern"` therefore evaluates true while
    // `missing =~ "pattern"` evaluates false.  Keep this distinct from the
    // configurable numeric missing-values policy used by ordinary fields.
    if (!present && expression.string_compare && expression.regex_compare)
        return {expression.op == "!~", true};
    if (!present) return {missing_fails, false};
    if (expression.op.empty()) return {boolean_value, true};
    if (expression.string_compare) return {boolean_value, true};
    if (expression.numeric_right != nullptr &&
        !evaluate_numeric_expression(*expression.numeric_right, header, record, threshold,
                                     allele_index, sample_index))
        return {missing_fails, false};
    return {compare(value, expression.op, threshold), true};
}

struct GenotypeFilterResult {
    std::uint64_t filtered = 0;
    std::uint64_t set_to_no_call = 0;
};

GenotypeFilterResult apply_genotype_filters(
    const bcf_hdr_t* input_header, const bcf_hdr_t* output_header, bcf1_t* record,
    const std::vector<std::shared_ptr<Expr>>& rules,
    const std::vector<std::string>& names, bool missing_fails,
    bool invert_expression, bool set_filtered_to_nocall) {
    if (record->n_sample <= 0 || (rules.empty() && !set_filtered_to_nocall)) return {};
    std::vector<std::string> statuses(static_cast<std::size_t>(record->n_sample), "PASS");
    // VariantFiltration's Java getGenotypeFilters() starts with any existing
    // genotype FT labels and appends newly matched labels.  Preserve those
    // labels instead of resetting every sample to PASS; this matters when a
    // pipeline applies multiple filtration stages to the same VCF.
    char** existing_filters = nullptr;
    int existing_count = 0;
    const auto existing_length = bcf_get_format_string(
        input_header, record, "FT", &existing_filters, &existing_count);
    if (existing_length > 0 && existing_filters != nullptr &&
        existing_count >= record->n_sample) {
        for (int sample = 0; sample < record->n_sample; ++sample) {
            const auto* observed = existing_filters[sample];
            if (observed != nullptr && *observed != '\0' && std::strcmp(observed, ".") != 0 &&
                std::strcmp(observed, "PASS") != 0)
                statuses[static_cast<std::size_t>(sample)] = observed;
        }
    }
    if (existing_filters != nullptr) {
        if (existing_count > 0) free(existing_filters[0]);
        free(existing_filters);
    }
    GenotypeFilterResult result;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        const char* sample_name = bcf_hdr_int2id(input_header, BCF_DT_SAMPLE, sample);
        if (sample_name == nullptr) continue;
        for (std::size_t rule_index = 0; rule_index < rules.size(); ++rule_index) {
            const auto& rule = *rules[rule_index];
            // The current parser accepts the explicit GATK form
            // vc.getGenotype("SAMPLE").(...).  A rule is evaluated only for
            // its named sample; an empty sample is reserved for a future
            // all-sample expression lowering.
            if (rule.kind == Expr::Kind::Predicate && !rule.sample.empty() &&
                rule.sample != sample_name)
                continue;
            const auto eval_result = evaluate(rule, input_header, record, missing_fails, -1, sample);
            const bool failed = invert_expression ? !eval_result.value : eval_result.value;
            if (!failed) continue;
            auto& status = statuses[static_cast<std::size_t>(sample)];
            if (status == "PASS" || status.empty())
                status = names[rule_index];
            else
                status += ';' + names[rule_index];
            ++result.filtered;
        }
        if (statuses[static_cast<std::size_t>(sample)].empty())
            statuses[static_cast<std::size_t>(sample)] = "PASS";
        statuses[static_cast<std::size_t>(sample)] = canonical_genotype_filter_status(
            statuses[static_cast<std::size_t>(sample)]);
    }
    // With only --set-filtered-genotype-to-no-call, GATK normalizes an
    // existing FORMAT/FT field (including '.' -> PASS) while preserving the
    // absence of FT when the input did not declare it.
    if (!rules.empty() || (set_filtered_to_nocall && existing_length > 0)) {
        std::vector<const char*> values;
        values.reserve(statuses.size());
        for (const auto& status : statuses) values.push_back(status.c_str());
        if (bcf_update_format_string(output_header, record, "FT", values.data(),
                                     static_cast<int>(values.size())) != 0)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update genotype FT");
    }

    // GATK's --set-filtered-genotype-to-no-call changes only called
    // genotypes that have a non-PASS FT status (including a pre-existing FT),
    // and turns every allele in that genotype into an unphased no-call.  Keep
    // the operation after FT materialization so existing chained filters are
    // honored even when this invocation has no new genotype expression.
    if (set_filtered_to_nocall) {
        int32_t* encoded = nullptr;
        int count = 0;
        const auto length = bcf_get_genotypes(input_header, record, &encoded, &count);
        if (length > 0 && encoded != nullptr && count > 0 &&
            count % record->n_sample == 0) {
            std::vector<int32_t> genotypes(encoded, encoded + count);
            const int ploidy = count / record->n_sample;
            for (int sample = 0; sample < record->n_sample; ++sample) {
                const auto& status = statuses[static_cast<std::size_t>(sample)];
                const bool filtered = !status.empty() && status != "PASS" && status != ".";
                if (!filtered) continue;
                bool called = false;
                for (int allele = 0; allele < ploidy; ++allele) {
                    const auto value = genotypes[static_cast<std::size_t>(sample * ploidy + allele)];
                    if (value != bcf_int32_vector_end && !bcf_gt_is_missing(value)) {
                        called = true;
                        break;
                    }
                }
                if (!called) continue;
                for (int allele = 0; allele < ploidy; ++allele)
                    genotypes[static_cast<std::size_t>(sample * ploidy + allele)] = bcf_gt_missing;
                ++result.set_to_no_call;
            }
            if (result.set_to_no_call > 0 &&
                bcf_update_genotypes(output_header, record, genotypes.data(), count) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot set filtered GT to no-call");
        }
        free(encoded);
    }
    return result;
}

void recompute_chromosome_counts(const bcf_hdr_t* header, bcf1_t* record) {
    if (record->n_sample <= 0 || record->n_allele <= 1) return;
    int32_t* encoded = nullptr;
    int count = 0;
    const auto length = bcf_get_genotypes(header, record, &encoded, &count);
    if (length <= 0 || encoded == nullptr || count <= 0 || count % record->n_sample != 0) {
        free(encoded);
        return;
    }
    const int ploidy = count / record->n_sample;
    std::vector<int32_t> ac(static_cast<std::size_t>(record->n_allele - 1), 0);
    int32_t an = 0;
    for (int sample = 0; sample < record->n_sample; ++sample) {
        for (int allele = 0; allele < ploidy; ++allele) {
            const auto value = encoded[sample * ploidy + allele];
            if (value == bcf_int32_vector_end || bcf_gt_is_missing(value)) continue;
            const int index = bcf_gt_allele(value);
            if (index < 0 || index >= record->n_allele) continue;
            ++an;
            if (index > 0) ++ac[static_cast<std::size_t>(index - 1)];
        }
    }
    free(encoded);
    if (bcf_update_info_int32(header, record, "AC", ac.data(), static_cast<int>(ac.size())) != 0 ||
        bcf_update_info_int32(header, record, "AN", &an, 1) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update AC/AN");
    std::vector<float> af(ac.size(), 0.0F);
    for (std::size_t index = 0; index < ac.size(); ++index)
        af[index] = an == 0 ? 0.0F : static_cast<float>(ac[index]) / static_cast<float>(an);
    if (bcf_update_info_float(header, record, "AF", af.data(), static_cast<int>(af.size())) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update AF");
}

std::uint64_t apply_allele_specific_filters(
    const bcf_hdr_t* input_header, const bcf_hdr_t* output_header, bcf1_t* record,
    const std::vector<std::shared_ptr<Expr>>& rules,
    const std::vector<std::string>& names, bool missing_fails, bool invert_expression) {
    if (record->n_allele <= 1 || rules.empty()) return 0;
    std::vector<std::string> statuses(static_cast<std::size_t>(record->n_allele - 1), "PASS");
    std::uint64_t filtered = 0;
    for (int allele = 0; allele < record->n_allele - 1; ++allele) {
        for (std::size_t rule_index = 0; rule_index < rules.size(); ++rule_index) {
            if (!uses_allele_specific_field(*rules[rule_index])) continue;
            const auto result = evaluate(*rules[rule_index], input_header, record,
                                          missing_fails, allele);
            const bool failed = invert_expression ? !result.value : result.value;
            if (!failed) continue;
            auto& status = statuses[static_cast<std::size_t>(allele)];
            if (status == "PASS" || status.empty()) status = names[rule_index];
            else status += ';' + names[rule_index];
            ++filtered;
        }
    }
    std::ostringstream encoded;
    for (std::size_t index = 0; index < statuses.size(); ++index) {
        if (index != 0) encoded << ',';
        encoded << statuses[index];
    }
    const auto value = encoded.str();
    if (bcf_update_info_string(output_header, record, "AS_FilterStatus", value.c_str()) != 0)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update AS_FilterStatus");
    return filtered;
}

struct MaskInterval {
    std::string contig;
    int begin = 0;
    int end = 0;
};

struct Region {
    int rid = -1;
    int begin = 0;
    int end = std::numeric_limits<int>::max();
};

void append_region_selector(const std::string& selector, const bcf_hdr_t* header,
                            std::vector<Region>& regions,
                            fastgatk::io::IntervalFileStats& stats) {
    std::vector<fastgatk::io::IndexedInterval> parsed;
    fastgatk::io::append_interval_selector(selector, header, parsed, stats);
    for (const auto& interval : parsed)
        regions.push_back(Region{interval.rid, interval.begin, interval.end});
}

const char* interval_set_rule_name(fastgatk::io::HtsIntervalSetRule rule) {
    return rule == fastgatk::io::HtsIntervalSetRule::Intersection
        ? "INTERSECTION" : "UNION";
}

void normalize_filter_regions(std::vector<Region>& regions) {
    std::sort(regions.begin(), regions.end(), [](const Region& left, const Region& right) {
        if (left.rid != right.rid) return left.rid < right.rid;
        if (left.begin != right.begin) return left.begin < right.begin;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(regions.size());
    for (const auto& region : regions) {
        if (region.end <= region.begin) continue;
        if (!merged.empty() && merged.back().rid == region.rid && region.begin <= merged.back().end)
            merged.back().end = std::max(merged.back().end, region.end);
        else merged.push_back(region);
    }
    regions.swap(merged);
}

void append_filter_region_selector_with_rule(
    const std::string& selector,
    const bcf_hdr_t* header,
    std::vector<Region>& regions,
    fastgatk::io::IntervalFileStats& stats,
    fastgatk::io::HtsIntervalSetRule rule,
    bool first_selector) {
    std::vector<Region> incoming;
    append_region_selector(selector, header, incoming, stats);
    normalize_filter_regions(incoming);
    if (rule == fastgatk::io::HtsIntervalSetRule::Union || first_selector) {
        regions.insert(regions.end(), incoming.begin(), incoming.end());
        return;
    }
    std::vector<Region> intersection;
    std::size_t left = 0;
    std::size_t right = 0;
    while (left < regions.size() && right < incoming.size()) {
        if (regions[left].rid < incoming[right].rid) { ++left; continue; }
        if (incoming[right].rid < regions[left].rid) { ++right; continue; }
        const int begin = std::max(regions[left].begin, incoming[right].begin);
        const int end = std::min(regions[left].end, incoming[right].end);
        if (begin < end) intersection.push_back(Region{regions[left].rid, begin, end});
        if (regions[left].end < incoming[right].end) ++left;
        else ++right;
    }
    regions.swap(intersection);
}

int record_end_exclusive(const bcf_hdr_t* header, bcf1_t* record) {
    int end = static_cast<int>(record->pos + std::max<hts_pos_t>(1, record->rlen));
    int32_t* values = nullptr;
    int count = 0;
    const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
    if (length > 0 && count > 0 && values[0] > 0)
        end = std::max(end, values[0]);
    free(values);
    return end;
}

bool in_regions(const bcf_hdr_t* header, bcf1_t* record,
                const std::vector<Region>& regions) {
    if (regions.empty()) return true;
    const auto end = record_end_exclusive(header, record);
    return std::any_of(regions.begin(), regions.end(), [&](const Region& region) {
        return record->rid == region.rid && end > region.begin && record->pos < region.end;
    });
}

std::vector<MaskInterval> load_mask_intervals(const std::vector<std::string>& paths) {
    std::vector<MaskInterval> intervals;
    for (const auto& path : paths) {
        htsFile* input = bcf_open(path.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open VariantFiltration mask: " + path);
        bcf_hdr_t* header = bcf_hdr_read(input);
        if (!header) {
            bcf_close(input);
            throw std::runtime_error("BAD_INPUT: cannot read VariantFiltration mask header: " + path);
        }
        bcf1_t* record = bcf_init();
        if (!record) {
            bcf_hdr_destroy(header);
            bcf_close(input);
            throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate VariantFiltration mask record");
        }
        while (bcf_read(input, header, record) == 0) {
            bcf_unpack(record, BCF_UN_ALL);
            const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
            if (contig != nullptr && record->pos >= 0) {
                const auto span = static_cast<int>(std::max<hts_pos_t>(1, record->rlen));
                int end = static_cast<int>(record->pos) + span;
                int32_t* values = nullptr;
                int count = 0;
                const auto length = bcf_get_info_int32(header, record, "END", &values, &count);
                if (length > 0 && count > 0 && values[0] > 0)
                    end = std::max(end, static_cast<int>(values[0]));
                free(values);
                intervals.push_back({contig, static_cast<int>(record->pos), end});
            }
            bcf_clear(record);
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        bcf_close(input);
    }
    return intervals;
}

bool overlaps_mask(const bcf_hdr_t* header, const bcf1_t* record,
                   const std::vector<MaskInterval>& intervals, int extension) {
    if (intervals.empty() || record->rid < 0 || record->pos < 0) return false;
    const auto* contig = bcf_hdr_id2name(header, record->rid);
    if (contig == nullptr) return false;
    const int begin = static_cast<int>(record->pos) - extension;
    // Apply the same record-span contract used by interval selection: a gVCF
    // block's INFO/END extends the masked record beyond HTSlib's rlen (which
    // is usually one for symbolic <NON_REF> records).
    const int end = record_end_exclusive(header, const_cast<bcf1_t*>(record)) + extension;
    return std::any_of(intervals.begin(), intervals.end(), [&](const MaskInterval& interval) {
        return interval.contig == contig && begin < interval.end && end > interval.begin;
    });
}

using ClusterPosition = std::pair<std::string, int>;

std::set<ClusterPosition> load_cluster_positions(const std::string& path,
                                                  int cluster_size,
                                                  int cluster_window_size) {
    std::map<std::string, std::vector<int>> by_contig;
    htsFile* input = bcf_open(path.c_str(), "r");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open VariantFiltration input for clustering: " + path);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (!header) {
        bcf_close(input);
        throw std::runtime_error("BAD_INPUT: cannot read VariantFiltration clustering header: " + path);
    }
    bcf1_t* record = bcf_init();
    if (!record) {
        bcf_hdr_destroy(header);
        bcf_close(input);
        throw std::runtime_error("RESOURCE_EXHAUSTED: cannot allocate VariantFiltration clustering record");
    }
    while (bcf_read(input, header, record) == 0) {
        bcf_unpack(record, BCF_UN_STR);
        const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(header, record->rid) : nullptr;
        if (contig != nullptr && record->pos >= 0)
            by_contig[contig].push_back(static_cast<int>(record->pos));
        bcf_clear(record);
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    bcf_close(input);

    std::set<ClusterPosition> clustered;
    for (auto& [contig, positions] : by_contig) {
        std::sort(positions.begin(), positions.end());
        if (positions.size() < static_cast<std::size_t>(cluster_size)) continue;
        std::vector<int> difference(positions.size() + 1, 0);
        std::size_t left = 0;
        for (std::size_t right = 0; right < positions.size(); ++right) {
            while (left <= right && positions[right] - positions[left] > cluster_window_size)
                ++left;
            if (right + 1 - left >= static_cast<std::size_t>(cluster_size)) {
                ++difference[left];
                --difference[right + 1];
            }
        }
        int active = 0;
        for (std::size_t index = 0; index < positions.size(); ++index) {
            active += difference[index];
            if (active > 0) clustered.emplace(contig, positions[index]);
        }
    }
    return clustered;
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    htsFile* input = nullptr;
    htsFile* output = nullptr;
    bcf_hdr_t* input_header = nullptr;
    bcf_hdr_t* output_header = nullptr;
    bcf1_t* record = nullptr;
    std::uint64_t input_records = 0;
    std::uint64_t filtered_records = 0;
    std::uint64_t genotype_filtered = 0;
    std::uint64_t genotype_set_to_nocall = 0;
    std::uint64_t masked_records = 0;
    std::uint64_t clustered_records = 0;
    std::uint64_t allele_filtered = 0;
    std::uint64_t interval_skipped = 0;
    fastgatk::io::IntervalFileStats interval_file_stats;
    try {
        std::vector<std::shared_ptr<Expr>> rules;
        rules.reserve(options.expressions.size());
        bool allele_specific_rules_present = false;
        for (const auto& expression : options.expressions) {
            auto parsed = parse_expression(expression);
            const bool allele_specific = uses_allele_specific_field(*parsed);
            if (allele_specific) allele_specific_rules_present = true;
            if (allele_specific && !options.apply_allele_specific_filters)
                throw std::invalid_argument(
                    "UNSUPPORTED_PARAMETER: AS_* expression requires --apply-allele-specific-filters");
            rules.push_back(std::move(parsed));
        }
        std::vector<std::shared_ptr<Expr>> genotype_rules;
        genotype_rules.reserve(options.genotype_expressions.size());
        for (const auto& expression : options.genotype_expressions)
            genotype_rules.push_back(parse_expression(expression, true));
        input = bcf_open(options.input.c_str(), "r");
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open variant input: " + options.input);
        input_header = bcf_hdr_read(input);
        if (!input_header) throw std::runtime_error("BAD_INPUT: cannot read variant header");
        std::vector<Region> regions;
        regions.reserve(options.regions.size());
        bool first_selector = true;
        for (const auto& text : options.regions) {
            append_filter_region_selector_with_rule(
                text, input_header, regions, interval_file_stats,
                options.interval_set_rule, first_selector);
            first_selector = false;
        }
        normalize_filter_regions(regions);
        const auto mask_intervals = load_mask_intervals(options.masks);
        const auto clustered_positions = options.cluster_size == 0
            ? std::set<ClusterPosition>{}
            : load_cluster_positions(options.input, options.cluster_size, options.cluster_window_size);
        output_header = bcf_hdr_dup(input_header);
        if (!output_header) throw std::runtime_error("RESOURCE_EXHAUSTED: cannot duplicate variant header");
        for (const auto& name : options.names) {
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, name.c_str()) < 0)
                bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=fastgatk VariantFiltration rule>", name.c_str());
        }
        if (!options.masks.empty() && bcf_hdr_id2int(output_header, BCF_DT_ID, options.mask_name.c_str()) < 0) {
            const std::string description = options.mask_description.empty()
                ? (options.filter_records_not_in_mask
                    ? "Doesn't overlap a user-input mask"
                    : "Overlaps a user-input mask")
                : options.mask_description;
            // bcf_hdr_printf uses a printf format string, so keep the user
            // supplied description as data.  This also preserves literal '%'
            // characters in a GATK-compatible --mask-description value.
            bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=%s>",
                           options.mask_name.c_str(), description.c_str());
        }
        if (options.apply_allele_specific_filters && allele_specific_rules_present && bcf_hdr_id2int(
                output_header, BCF_DT_ID, "AS_FilterStatus") < 0)
            bcf_hdr_append(output_header,
                           "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=Allele-specific filter status>");
        const std::string cluster_name = "ClusteredEvents";
        if (!clustered_positions.empty() && bcf_hdr_id2int(output_header, BCF_DT_ID, cluster_name.c_str()) < 0)
            bcf_hdr_printf(output_header, "##FILTER=<ID=%s,Description=fastgatk VariantFiltration cluster>",
                           cluster_name.c_str());
        bcf_hdr_append(output_header, "##source=fastgatk-variant-filtration");
        bcf_hdr_append(output_header, "##fastgatk_variant_filtration_status=prototype-jexl-subset");
        if (!genotype_rules.empty() && bcf_hdr_id2int(output_header, BCF_DT_ID, "FT") < 0)
            bcf_hdr_append(output_header, "##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype-level filter>");
        if (options.set_filtered_genotypes_to_nocall) {
            // GATK advertises these standard annotations whenever a filtered
            // called genotype may be replaced with no-call and counts change.
            // Preserve existing declarations for chained VCF stages.
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AC") < 0)
                bcf_hdr_append(output_header,
                               "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count in genotypes>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AN") < 0)
                bcf_hdr_append(output_header,
                               "##INFO=<ID=AN,Number=1,Type=Integer,Description=Total number of alleles in called genotypes>");
            if (bcf_hdr_id2int(output_header, BCF_DT_ID, "AF") < 0)
                bcf_hdr_append(output_header,
                               "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>");
        }
        if (bcf_hdr_sync(output_header) != 0) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot sync output header");
        output = bcf_open(options.output.c_str(), suffix(options.output, ".gz") ? "wz" : "w");
        if (!output) throw std::runtime_error("cannot open VariantFiltration output: " + options.output);
        if (bcf_hdr_write(output, output_header) != 0) throw std::runtime_error("cannot write VariantFiltration header");
        record = bcf_init();
        if (!record) throw std::runtime_error("RESOURCE_EXHAUSTED: bcf_init failed");
        while (bcf_read(input, input_header, record) == 0) {
            ++input_records;
            bcf_unpack(record, BCF_UN_ALL);
            if (!in_regions(input_header, record, regions)) {
                ++interval_skipped;
                bcf_clear(record);
                continue;
            }
            // GATK's invalidate-previous-filters only resets the variant-level
            // FILTER set. Genotype FT values are deliberately left intact;
            // chained site/genotype filtering relies on that distinction.
            if (options.invalidate_previous_filters &&
                bcf_update_filter(output_header, record, nullptr, 0) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot invalidate previous FILTER values");
            bool any_failed = false;
            for (std::size_t index = 0; index < rules.size(); ++index) {
                // AS_* predicates are per-ALT annotations.  They must not
                // leak into the site FILTER column (which has no ALT index).
                if (uses_allele_specific_field(*rules[index])) continue;
                const auto result = evaluate(*rules[index], input_header, record, options.missing_fails);
                const bool failed = options.invert_filter_expression ? !result.value : result.value;
                if (!failed) continue;
                // GATK appends a new site FILTER to an already filtered record;
                // it clears only PASS/no-filter state.  Clearing unconditionally
                // loses provenance when VariantFiltration is chained.
                if (!any_failed && record->d.n_flt == 0 &&
                    bcf_update_filter(output_header, record, nullptr, 0) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear PASS filter");
                any_failed = true;
                const int filter_id = bcf_hdr_id2int(output_header, BCF_DT_ID, options.names[index].c_str());
                if (filter_id < 0 || bcf_add_filter(output_header, record, filter_id) < 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update FILTER column");
            }
            // GATK's addMaskIfCoversVariant uses
            // `maskVariants.isEmpty() == filterRecordsNotInMask`: normal mode
            // filters an overlapping record, while --filter-not-in-mask
            // filters a record that has no overlap.  Keep the overlap query
            // unchanged so extension and END-span semantics are shared by
            // both modes.  Since mask_overlap is the inverse of
            // maskVariants.isEmpty(), the filter condition is `!=` here.
            const bool mask_overlap = overlaps_mask(
                input_header, record, mask_intervals, options.mask_extension);
            if (mask_overlap != options.filter_records_not_in_mask) {
                if (!any_failed && record->d.n_flt == 0 &&
                    bcf_update_filter(output_header, record, nullptr, 0) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear PASS filter for mask");
                any_failed = true;
                const int mask_id = bcf_hdr_id2int(output_header, BCF_DT_ID, options.mask_name.c_str());
                if (mask_id < 0 || bcf_add_filter(output_header, record, mask_id) < 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update mask FILTER column");
                ++masked_records;
            }
            const auto* contig = record->rid >= 0 ? bcf_hdr_id2name(input_header, record->rid) : nullptr;
            if (contig != nullptr && clustered_positions.contains({contig, static_cast<int>(record->pos)})) {
                if (!any_failed && record->d.n_flt == 0 &&
                    bcf_update_filter(output_header, record, nullptr, 0) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot clear PASS filter for cluster");
                any_failed = true;
                const int cluster_id = bcf_hdr_id2int(output_header, BCF_DT_ID, cluster_name.c_str());
                if (cluster_id < 0 || bcf_add_filter(output_header, record, cluster_id) < 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot update cluster FILTER column");
                ++clustered_records;
            }
            if (any_failed) ++filtered_records;
            if (options.apply_allele_specific_filters && allele_specific_rules_present)
                allele_filtered += apply_allele_specific_filters(
                    input_header, output_header, record, rules, options.names,
                    options.missing_fails, options.invert_filter_expression);
            const auto genotype_result = apply_genotype_filters(
                input_header, output_header, record, genotype_rules,
                options.genotype_names, options.missing_fails,
                options.invert_genotype_filter_expression,
                options.set_filtered_genotypes_to_nocall);
            genotype_filtered += genotype_result.filtered;
            genotype_set_to_nocall += genotype_result.set_to_no_call;
            if (genotype_result.set_to_no_call > 0)
                recompute_chromosome_counts(output_header, record);
            sort_filter_labels(output_header, record);
            if (bcf_translate(output_header, input_header, record) != 0 ||
                bcf_write(output, output_header, record) != 0)
                throw std::runtime_error("cannot write VariantFiltration record");
            bcf_clear(record);
        }
        bcf_destroy(record); record = nullptr;
        bcf_close(output); output = nullptr;
        bcf_close(input); input = nullptr;
        std::string index_path;
        if (options.create_index) {
            if (suffix(options.output, ".gz")) {
                index_path = options.output + ".tbi";
                if (tbx_index_build3(options.output.c_str(), index_path.c_str(), 0, 0, &tbx_conf_vcf) != 0)
                    throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot build VariantFiltration index");
            } else {
                index_path = options.output + ".idx";
                fastgatk::io::write_uncompressed_vcf_tribble_index(options.output, index_path);
            }
        }
        const bool primary_complete = file_complete(options.output);
        const bool index_complete = index_path.empty() || file_complete(index_path);
        if (!primary_complete || !index_complete)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: VariantFiltration output or index is missing/empty");
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("cannot write VariantFiltration manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"VariantFiltration\",\"implementation\":\"fastgatk-variant-filtration\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"vcf\","
                 << "\"compatibility\":{\"jexl_boolean_subset\":true,\"allele_methods\":true,\"genotype_predicates\":true,\"missing_values_failing\":"
                 << (options.missing_fails ? "true" : "false") << ",\"genotype_ft\":"
                 << (!genotype_rules.empty() ? "true" : "false")
                 << ",\"set_filtered_genotype_to_no_call\":"
                 << (options.set_filtered_genotypes_to_nocall ? "true" : "false")
                 << ",\"jexl_info_vector_indexing\":true,\"jexl_boolean_comparisons\":true,\"jexl_regex_operators\":true,\"jexl_arithmetic\":true"
                 << ",\"invert_filter_expression\":" << (options.invert_filter_expression ? "true" : "false")
                 << ",\"invert_genotype_filter_expression\":"
                 << (options.invert_genotype_filter_expression ? "true" : "false")
                 << ",\"allele_specific_filters\":"
                 << (options.apply_allele_specific_filters && allele_specific_rules_present ? "true" : "false")
                 << ",\"mask_filter\":" << (!options.masks.empty() ? "true" : "false")
                 << ",\"filter_not_in_mask\":"
                 << (options.filter_records_not_in_mask ? "true" : "false")
                 << ",\"mask_description\":\"" << json_escape(options.mask_description) << "\""
                 << ",\"mask_extension\":" << options.mask_extension
                 << ",\"cluster_filter\":" << (options.cluster_size != 0 ? "true" : "false")
                 << ",\"cluster_size\":" << options.cluster_size
                 << ",\"cluster_window_size\":" << options.cluster_window_size
                 << ",\"invalidate_previous_filters\":"
                 << (options.invalidate_previous_filters ? "true" : "false")
                 << ",\"interval_subset\":" << (!regions.empty() ? "true" : "false")
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"vcf_index\":"
                 << (index_path.empty() ? "false" : "true")
                 << ",\"bit_identical_to_gatk\":false},\"outputs\":[{\"path\":\"" << json_escape(options.output)
                 << "\",\"kind\":\"vcf\",\"complete\":true}"
                 << (index_path.empty() ? "" : ",{\"path\":\"" + json_escape(index_path) + "\",\"kind\":\"vcf-index\",\"complete\":true}")
                 << "],\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"input_records\":" << input_records << ",\"filtered_records\":" << filtered_records
                 << ",\"genotype_filtered\":" << genotype_filtered << ",\"masked_records\":" << masked_records
                 << ",\"genotype_set_to_no_call\":" << genotype_set_to_nocall
                 << ",\"clustered_records\":" << clustered_records
                 << ",\"allele_filtered\":" << allele_filtered
                 << ",\"invalidate_previous_filters\":"
                 << (options.invalidate_previous_filters ? "true" : "false")
                 << ",\"interval_skipped\":" << interval_skipped
                 << ",\"interval_set_rule\":\""
                 << interval_set_rule_name(options.interval_set_rule) << "\""
                 << ",\"interval_list_inputs\":" << interval_file_stats.files
                 << ",\"interval_list_records\":" << interval_file_stats.records << "}}\n";
        std::cout << "{\"tool\":\"VariantFiltration\",\"status\":\"prototype\",\"input_records\":"
                  << input_records << ",\"filtered_records\":" << filtered_records
                  << ",\"genotype_filtered\":" << genotype_filtered << ",\"masked_records\":" << masked_records
                  << ",\"genotype_set_to_no_call\":" << genotype_set_to_nocall
                  << ",\"clustered_records\":" << clustered_records
                  << ",\"filter_not_in_mask\":"
                  << (options.filter_records_not_in_mask ? "true" : "false")
                  << ",\"allele_filtered\":" << allele_filtered
                  << ",\"invalidate_previous_filters\":"
                  << (options.invalidate_previous_filters ? "true" : "false")
                  << ",\"interval_skipped\":" << interval_skipped
                  << ",\"interval_set_rule\":\""
                  << interval_set_rule_name(options.interval_set_rule) << "\"}\n";
        bcf_hdr_destroy(output_header);
        return 0;
    } catch (...) {
        bcf_destroy(record);
        if (output) bcf_close(output);
        if (input) bcf_close(input);
        if (output_header) bcf_hdr_destroy(output_header);
        if (input_header) bcf_hdr_destroy(input_header);
        throw;
    }
}

#else
int run_tool(const Options&, const fastgatk::runtime::ResourceSnapshot&) {
    throw std::runtime_error("BACKEND_UNAVAILABLE: build with HTSlib for VariantFiltration");
}
#endif

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-variant-filtration: " << error.what() << '\n';
        return 2;
    }
}
