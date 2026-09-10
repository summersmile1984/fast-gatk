#include "fastgatk/runtime/output.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <vector>

namespace fastgatk::runtime {
namespace {

constexpr std::size_t kMaxManifestBytes = 64U * 1024U * 1024U;

bool regular_nonempty(const std::filesystem::path& path, std::uintmax_t& bytes,
                      std::string& reason, const char* label) {
    std::error_code error;
    if (path.empty()) {
        reason = std::string(label) + " path is empty";
        return false;
    }
    if (!std::filesystem::is_regular_file(path, error) || error) {
        reason = std::string(label) + " is missing or not a regular file: " + path.string();
        return false;
    }
    bytes = std::filesystem::file_size(path, error);
    if (error || bytes == 0) {
        reason = std::string(label) + " is empty or cannot be sized: " + path.string();
        return false;
    }
    return true;
}

bool has_boolean_field(const std::string& text, const std::string& field,
                       const bool expected) {
    const std::string key = "\"" + field + "\"";
    std::size_t cursor = 0;
    while ((cursor = text.find(key, cursor)) != std::string::npos) {
        cursor += key.size();
        const auto colon = text.find(':', cursor);
        if (colon == std::string::npos) return false;
        std::size_t value = colon + 1;
        while (value < text.size() &&
               std::isspace(static_cast<unsigned char>(text[value]))) ++value;
        const auto token = expected ? "true" : "false";
        if (text.compare(value, std::char_traits<char>::length(token), token) == 0)
            return true;
    }
    return false;
}

bool has_nonempty_array_field(const std::string& text, const std::string& field) {
    const std::string key = "\"" + field + "\"";
    const auto key_position = text.find(key);
    if (key_position == std::string::npos) return false;
    const auto colon = text.find(':', key_position + key.size());
    if (colon == std::string::npos) return false;
    const auto open = text.find('[', colon + 1);
    if (open == std::string::npos) return false;
    std::size_t first = open + 1;
    while (first < text.size() &&
           std::isspace(static_cast<unsigned char>(text[first]))) ++first;
    return first < text.size() && text[first] != ']';
}

bool has_numeric_schema_version(const std::string& text) {
    const std::string key = "\"schema_version\"";
    const auto key_position = text.find(key);
    if (key_position == std::string::npos) return false;
    const auto colon = text.find(':', key_position + key.size());
    if (colon == std::string::npos) return false;
    std::size_t cursor = colon + 1;
    while (cursor < text.size() &&
           std::isspace(static_cast<unsigned char>(text[cursor]))) ++cursor;
    if (cursor >= text.size() || !std::isdigit(static_cast<unsigned char>(text[cursor])))
        return false;
    while (cursor < text.size() &&
           std::isdigit(static_cast<unsigned char>(text[cursor]))) ++cursor;
    return true;
}

bool outputs_are_complete(const std::string& text) {
    const std::string key = "\"outputs\"";
    const auto key_position = text.find(key);
    if (key_position == std::string::npos) return false;
    const auto colon = text.find(':', key_position + key.size());
    if (colon == std::string::npos) return false;
    const auto open = text.find('[', colon + 1);
    if (open == std::string::npos) return false;

    std::size_t object_start = std::string::npos;
    int object_depth = 0;
    bool in_string = false;
    bool escaped = false;
    std::size_t objects = 0;
    for (std::size_t cursor = open + 1; cursor < text.size(); ++cursor) {
        const char character = text[cursor];
        if (in_string) {
            if (escaped) escaped = false;
            else if (character == '\\') escaped = true;
            else if (character == '"') in_string = false;
            continue;
        }
        if (character == '"') {
            in_string = true;
            continue;
        }
        if (character == '{') {
            if (object_depth == 0) object_start = cursor;
            ++object_depth;
            continue;
        }
        if (character == '}' && object_depth > 0) {
            --object_depth;
            if (object_depth == 0 && object_start != std::string::npos) {
                const auto object = text.substr(object_start, cursor - object_start + 1);
                if (has_boolean_field(object, "complete", false) ||
                    !has_boolean_field(object, "complete", true))
                    return false;
                ++objects;
                object_start = std::string::npos;
            }
            continue;
        }
        if (character == ']' && object_depth == 0)
            return objects != 0;
    }
    return false;
}

bool has_temporary_residue(const OutputBundle& bundle) {
    const auto has_residue = [](const std::filesystem::path& path) {
        if (path.empty()) return false;
        const std::string base = path.string();
        std::error_code error;
        for (const auto& suffix : {std::string{".tmp"}, std::string{".partial"},
                                   std::string{".incomplete"}}) {
            if (std::filesystem::exists(base + suffix, error) && !error) return true;
            error.clear();
        }
        return false;
    };
    return has_residue(bundle.primary) || has_residue(bundle.index) ||
           has_residue(bundle.manifest);
}

bool read_manifest(const std::filesystem::path& path, std::string& text,
                   std::uintmax_t& bytes, std::string& reason) {
    if (!regular_nonempty(path, bytes, reason, "manifest")) return false;
    if (bytes > kMaxManifestBytes) {
        reason = "manifest exceeds 64 MiB safety bound: " + path.string();
        return false;
    }
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        reason = "manifest cannot be opened: " + path.string();
        return false;
    }
    text.assign(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
    if (!stream.good() && !stream.eof()) {
        reason = "manifest read failed: " + path.string();
        return false;
    }
    auto begin = text.find_first_not_of(" \t\r\n");
    auto end = text.find_last_not_of(" \t\r\n");
    if (begin == std::string::npos || text[begin] != '{' || text[end] != '}') {
        reason = "manifest is not a JSON object: " + path.string();
        return false;
    }
    if (!has_numeric_schema_version(text) || !has_nonempty_array_field(text, "outputs") ||
        !outputs_are_complete(text)) {
        reason = "manifest is missing numeric schema_version, non-empty outputs, or complete outputs: " + path.string();
        return false;
    }
    if (has_boolean_field(text, "complete", false) ||
        !has_boolean_field(text, "complete", true)) {
        reason = "manifest does not mark every output complete: " + path.string();
        return false;
    }
    return true;
}

std::string validation_message(const OutputValidation& validation) {
    return validation.reason.empty() ? "output bundle is incomplete" : validation.reason;
}

void ensure_distinct_destinations(const OutputBundle& bundle) {
    const auto same = [](const std::filesystem::path& left,
                         const std::filesystem::path& right) {
        return !left.empty() && !right.empty() && left.lexically_normal() == right.lexically_normal();
    };
    if (bundle.primary.empty() || bundle.manifest.empty() ||
        same(bundle.primary, bundle.manifest) || same(bundle.primary, bundle.index) ||
        same(bundle.index, bundle.manifest))
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: output bundle paths must be distinct and non-empty");
}

}  // namespace

OutputValidation validate_output_bundle(const OutputBundle& bundle) {
    OutputValidation validation;
    if (bundle.primary.empty() || bundle.manifest.empty() ||
        (bundle.require_index && bundle.index.empty())) {
        validation.reason = "primary, manifest, and required index paths must be present";
        return validation;
    }
    if (bundle.reject_temporary_residue && has_temporary_residue(bundle)) {
        validation.temporary_residue = true;
        validation.reason = "temporary output residue is present";
        return validation;
    }
    if (!regular_nonempty(bundle.primary, validation.primary_bytes,
                          validation.reason, "primary output")) return validation;
    if (!bundle.index.empty() &&
        !regular_nonempty(bundle.index, validation.index_bytes,
                          validation.reason, "index")) return validation;
    std::string manifest;
    if (!read_manifest(bundle.manifest, manifest, validation.manifest_bytes,
                       validation.reason)) return validation;
    validation.complete = true;
    return validation;
}

void require_complete_output(const OutputBundle& bundle) {
    const auto validation = validate_output_bundle(bundle);
    if (!validation.complete)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: " + validation_message(validation));
}

void publish_output_bundle(const OutputBundle& staging,
                           const OutputBundle& destination) {
    ensure_distinct_destinations(destination);
    if (staging.primary.empty() || staging.manifest.empty() ||
        (destination.require_index && staging.index.empty()))
        throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: staging bundle is missing a required path");

    auto staging_check = staging;
    staging_check.require_index = destination.require_index;
    // Staging names conventionally end in .tmp/.partial; only the destination
    // contract rejects sibling residue.  The files themselves are still
    // required to be regular, non-empty, and covered by a complete manifest.
    staging_check.reject_temporary_residue = false;
    require_complete_output(staging_check);

    const std::vector<std::pair<std::filesystem::path, std::filesystem::path>> files{
        {staging.primary, destination.primary},
        {staging.index, destination.index},
        {staging.manifest, destination.manifest},
    };
    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> moves;
    moves.reserve(files.size());
    std::error_code error;
    for (const auto& [from, to] : files) {
        if (from.empty() && to.empty()) continue;
        if (from.empty() || to.empty())
            throw std::invalid_argument("OUTPUT_CONTRACT_FAILURE: staging/destination sidecar mismatch");
        if (std::filesystem::exists(to, error) && !error)
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: destination already exists: " + to.string());
        error.clear();
        moves.emplace_back(from, to);
    }

    std::vector<std::pair<std::filesystem::path, std::filesystem::path>> published;
    try {
        for (const auto& [from, to] : moves) {
            std::filesystem::rename(from, to);
            published.emplace_back(from, to);
        }
        require_complete_output(destination);
    } catch (const std::exception& failure) {
        bool rollback_failed = false;
        for (auto it = published.rbegin(); it != published.rend(); ++it) {
            std::error_code rollback_error;
            std::filesystem::rename(it->second, it->first, rollback_error);
            rollback_failed = rollback_failed || static_cast<bool>(rollback_error);
        }
        std::ostringstream message;
        message << "OUTPUT_CONTRACT_FAILURE: atomic output publish failed: " << failure.what();
        if (rollback_failed) message << "; rollback also failed";
        throw std::runtime_error(message.str());
    }
}

}  // namespace fastgatk::runtime
