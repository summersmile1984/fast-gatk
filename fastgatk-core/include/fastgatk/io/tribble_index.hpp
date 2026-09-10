#pragma once

#include <string>

namespace fastgatk::io {

// Write the HTSJDK/Tribble LinearIndex v3 used by VariantContextWriter for
// an ordinary (uncompressed) VCF.  The caller chooses this path only for
// plain text output; BGZF VCF continues to use the HTSlib Tabix writer at the
// tool boundary.  The function is deliberately independent of HTSlib so the
// same implementation is available to every native variant tool.
void write_uncompressed_vcf_tribble_index(const std::string& input_path,
                                           const std::string& index_path);

}  // namespace fastgatk::io
