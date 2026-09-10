#include "fastgatk/calling/pipeline.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_trim(const std::string& label, const std::string& cigar,
                const std::size_t reference_length,
                const std::size_t haplotype_length,
                const std::size_t trim_start, const std::size_t trim_end,
                const std::size_t expected_begin, const std::size_t expected_end) {
    const auto actual = fastgatk::calling::trim_haplotype_by_reference_cigar(
        cigar, reference_length, haplotype_length, trim_start, trim_end);
    if (!actual || actual->sequence_begin != expected_begin ||
        actual->sequence_end != expected_end) {
        throw std::runtime_error("GATK Haplotype.trim mismatch for " + label);
    }
}

void check_rejected(const std::string& label, const std::string& cigar,
                    const std::size_t reference_length,
                    const std::size_t haplotype_length,
                    const std::size_t trim_start, const std::size_t trim_end) {
    if (fastgatk::calling::trim_haplotype_by_reference_cigar(
            cigar, reference_length, haplotype_length, trim_start, trim_end)) {
        throw std::runtime_error("GATK Haplotype.trim should reject " + label);
    }
}

void check_trimmed_population_order() {
    using Entry = fastgatk::calling::TrimmedHaplotypePopulationEntry;
    // This is the collision in AssemblyResultSet.trimDownHaplotypes(): a
    // reference source arriving after an equal-base non-reference source
    // must replace it, and calculateOriginalByTrimmedHaplotypes() then uses
    // Haplotype.SIZE_AND_BASE_ORDER (length, then bases).
    const std::vector<Entry> haplotypes{
        {"GGAA", false}, {"AT", false}, {"CCC", false},
        {"GGAA", true}, {"AA", false}, {"CC", true},
    };
    const std::vector<std::size_t> expected{4, 1, 5, 2, 3};
    const auto actual = fastgatk::calling::normalized_trimmed_haplotype_indices(haplotypes);
    if (actual != expected)
        throw std::runtime_error("GATK AssemblyResultSet.trimTo population order mismatch");
}

}  // namespace

int main() {
    try {
        // Direct port of HaplotypeUnitTest.
        // makeTrimmingDataWithLeadingAndTrailingInsertions().  The source
        // intervals are inclusive; the Host interface uses half-open spans.
        check_trim("internal insertion", "2M2I4M", 6, 8, 1, 6, 1, 8);
        check_trim("leading insertion", "2M2I4M", 6, 8, 2, 6, 4, 8);
        check_trim("later leading insertion", "2M2I4M", 6, 8, 3, 6, 5, 8);
        check_trim("internal insertion before right boundary", "2M2I4M", 6, 8, 0, 3, 0, 5);
        check_trim("trailing insertion", "2M2I4M", 6, 8, 0, 2, 0, 2);

        // Direct port of AlignmentUtilsUnitTest's deletion boundary matrix.
        // GATK returns null when either retained endpoint lies in a deletion.
        check_trim("deletion full span", "2M2D2M", 6, 4, 0, 6, 0, 4);
        check_trim("deletion suffix", "2M2D2M", 6, 4, 4, 6, 2, 4);
        check_trim("deletion left prefix", "2M2D2M", 6, 4, 0, 2, 0, 2);
        check_rejected("deletion left boundary", "2M2D2M", 6, 4, 2, 6);
        check_rejected("deletion right boundary", "2M2D2M", 6, 4, 0, 4);
        check_rejected("deletion internal endpoint", "2M2D2M", 6, 4, 0, 3);

        // AlignmentUtilsUnitTest also pins terminal insertions: they are not
        // swallowed by a full-reference trim, unlike a genuine internal I.
        check_trim("leading terminal insertion", "2I4M", 4, 6, 0, 4, 2, 6);
        check_trim("trailing terminal insertion", "4M2I", 4, 6, 0, 4, 0, 4);
        check_trimmed_population_order();

        std::cout << "{\"status\":\"pass\",\"source\":\"GATK HaplotypeUnitTest, AlignmentUtilsUnitTest, and AssemblyResultSet.trimTo\"}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
