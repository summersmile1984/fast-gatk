#include "fastgatk/runtime/output.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_file(const std::filesystem::path& path, const std::string& text) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("cannot create smoke file: " + path.string());
    stream << text;
    if (!stream) throw std::runtime_error("cannot write smoke file: " + path.string());
}

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto root = std::filesystem::temp_directory_path() /
        ("fastgatk-runtime-output-smoke-" + std::to_string(stamp));
    try {
        std::filesystem::create_directories(root);
        const auto final_primary = root / "calls.vcf";
        const auto final_index = root / "calls.vcf.tbi";
        const auto final_manifest = root / "calls.manifest.json";
        const auto staging_primary = root / "calls.vcf.fastgatk-stage";
        const auto staging_index = root / "calls.vcf.tbi.fastgatk-stage";
        const auto staging_manifest = root / "calls.manifest.json.fastgatk-stage";
        const fastgatk::runtime::OutputBundle destination{
            final_primary, final_index, final_manifest, true, true};
        const fastgatk::runtime::OutputBundle staging{
            staging_primary, staging_index, staging_manifest, true, false};
        const auto manifest = std::string(
            "{\"schema_version\":1,\"tool\":\"runtime-smoke\",\"outputs\":["
            "{\"path\":\"calls.vcf\",\"complete\":true},"
            "{\"path\":\"calls.vcf.tbi\",\"complete\":true}]}");

        // Complete staging is published only after all three artifacts pass
        // validation.  The manifest is moved last and acts as the commit mark.
        write_file(staging_primary, "##fileformat=VCFv4.3\n");
        write_file(staging_index, "tabix-index\n");
        write_file(staging_manifest, manifest);
        fastgatk::runtime::publish_output_bundle(staging, destination);
        const auto complete = fastgatk::runtime::validate_output_bundle(destination);
        require(complete.complete && complete.primary_bytes > 0 && complete.index_bytes > 0,
                "complete output bundle was not published");
        require(!std::filesystem::exists(staging_primary) &&
                    !std::filesystem::exists(staging_index) &&
                    !std::filesystem::exists(staging_manifest),
                "staging files remained after successful publish");

        // Every item in outputs[] must carry complete:true.  A manifest that
        // merely has one completed item must not hide a missing completion
        // marker for another sidecar.
        write_file(final_manifest,
                   "{\"schema_version\":1,\"outputs\":[{\"path\":\"calls.vcf\"}]}");
        const auto incomplete_manifest = fastgatk::runtime::validate_output_bundle(destination);
        require(!incomplete_manifest.complete &&
                    incomplete_manifest.reason.find("schema_version") != std::string::npos,
                "manifest with incomplete outputs[] was accepted");
        write_file(final_manifest, manifest);
        require(fastgatk::runtime::validate_output_bundle(destination).complete,
                "bundle did not recover after incomplete manifest replacement");

        // Missing index is an incomplete result even when primary and manifest
        // are present.  The throwing boundary uses a stable failure class.
        std::filesystem::remove(final_index);
        const auto missing_index = fastgatk::runtime::validate_output_bundle(destination);
        require(!missing_index.complete && missing_index.reason.find("index") != std::string::npos,
                "missing index was accepted");
        bool missing_index_threw = false;
        try {
            fastgatk::runtime::require_complete_output(destination);
        } catch (const std::runtime_error& error) {
            missing_index_threw = std::string(error.what()).find("OUTPUT_CONTRACT_FAILURE") !=
                std::string::npos;
        }
        require(missing_index_threw, "missing index did not fail closed");
        write_file(final_index, "tabix-index\n");

        // A leftover sibling temp file invalidates the final bundle.  This is
        // the state a killed writer must not expose to gather/resume.
        const auto temporary = std::filesystem::path(final_primary.string() + ".tmp");
        write_file(temporary, "partial\n");
        const auto residue = fastgatk::runtime::validate_output_bundle(destination);
        require(!residue.complete && residue.temporary_residue,
                "temporary output residue was accepted");
        std::filesystem::remove(temporary);
        require(fastgatk::runtime::validate_output_bundle(destination).complete,
                "bundle did not recover after temporary residue cleanup");

        // A failed commit must leave no destination artifact.  Validation is
        // performed before the first rename, so the staged primary remains
        // available for retry and the final path stays absent.
        std::filesystem::remove(final_primary);
        std::filesystem::remove(final_index);
        std::filesystem::remove(final_manifest);
        const auto failed_staging_primary = root / "failed.vcf.fastgatk-stage";
        const auto failed_staging_index = root / "failed.vcf.tbi.fastgatk-stage";
        const auto failed_staging_manifest = root / "failed.manifest.json.fastgatk-stage";
        write_file(failed_staging_primary, "partial candidate\n");
        write_file(failed_staging_manifest, manifest);
        const fastgatk::runtime::OutputBundle failed_staging{
            failed_staging_primary, failed_staging_index, failed_staging_manifest, true, false};
        bool publish_threw = false;
        try {
            fastgatk::runtime::publish_output_bundle(failed_staging, destination);
        } catch (const std::runtime_error& error) {
            publish_threw = std::string(error.what()).find("OUTPUT_CONTRACT_FAILURE") !=
                std::string::npos;
        }
        require(publish_threw, "incomplete staging bundle was published");
        require(!std::filesystem::exists(final_primary) &&
                    !std::filesystem::exists(final_index) &&
                    !std::filesystem::exists(final_manifest) &&
                    std::filesystem::exists(failed_staging_primary),
                "failed publish left a visible or unrecoverable destination");

        std::cout << "{\"status\":\"pass\",\"complete_publish\":true,"
                     "\"missing_index_fail_closed\":true,\"temporary_residue_fail_closed\":true,"
                     "\"failed_publish_atomic\":true}\n";
        std::error_code cleanup_error;
        std::filesystem::remove_all(root, cleanup_error);
        return 0;
    } catch (const std::exception& error) {
        std::error_code cleanup_error;
        std::filesystem::remove_all(root, cleanup_error);
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
