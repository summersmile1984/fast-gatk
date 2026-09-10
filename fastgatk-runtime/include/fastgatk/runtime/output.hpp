#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace fastgatk::runtime {

// The output bundle is the smallest publish unit shared by native tools:
// primary payload, optional/required index, and a completion manifest.  The
// runtime does not interpret tool-specific JSON fields, but it does require a
// bounded, complete manifest envelope before a bundle is considered visible.
struct OutputBundle {
    std::filesystem::path primary;
    std::filesystem::path index;
    std::filesystem::path manifest;
    bool require_index = false;
    bool reject_temporary_residue = true;
};

struct OutputValidation {
    bool complete = false;
    bool temporary_residue = false;
    std::string reason;
    std::uintmax_t primary_bytes = 0;
    std::uintmax_t index_bytes = 0;
    std::uintmax_t manifest_bytes = 0;
};

// Read-only validation.  It never creates, removes, or renames a file, so it
// is safe to use before a manifest is published or by a gather/resume path.
// A non-empty index path is always checked; require_index additionally makes
// an empty index path an error (for tools that promise indexed output).
OutputValidation validate_output_bundle(const OutputBundle& bundle);

// Throw OUTPUT_CONTRACT_FAILURE when the bundle is not complete.  This is the
// common fail-closed boundary for callers that are about to report success.
void require_complete_output(const OutputBundle& bundle);

// Publish a fully validated staging bundle.  All staging files must be on the
// same filesystem as their destinations and destinations must not already
// exist.  The manifest is renamed last as the completion/commit marker.  If a
// rename fails, already-published files are moved back to their staging paths;
// no destination is intentionally left visible after a failed commit.
void publish_output_bundle(const OutputBundle& staging,
                           const OutputBundle& destination);

}  // namespace fastgatk::runtime
