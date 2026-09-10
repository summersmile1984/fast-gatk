#include "fastgatk/io/flow_codec.hpp"
#include "fastgatk/kernels/pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Case {
    const char* name;
    const char* read;
    const char* haplotype;
    const char* flow_order;
    std::uint64_t expected_bits;
};

fastgatk::pairhmm::KokkosBatchResult score(const Case& item) {
    fastgatk::io::FlowReadTags tags;
    tags.bases.assign(item.read, item.read + std::strlen(item.read));
    tags.qualities.assign(tags.bases.size(), 30);
    tags.tp.assign(tags.bases.size(), 0);
    tags.flow_order = item.flow_order;
    tags.max_hmer = 12;
    tags.filling_value = 0.001;
    tags.keep_boundary_flows = true;
    const auto decoded = fastgatk::io::decode_flow_read(tags);
    const auto encoded = fastgatk::io::encode_flow_key(
        std::vector<std::uint8_t>(item.haplotype,
                                  item.haplotype + std::strlen(item.haplotype)),
        item.flow_order);

    fastgatk::pairhmm::FlowPairHmmRead read;
    read.key = decoded.key;
    read.flow_order = decoded.flow_order;
    read.probabilities = decoded.probabilities;
    read.insertion_gop.assign(read.key.size(), 40);
    read.deletion_gop.assign(read.key.size(), 40);
    read.gap_continuation.assign(read.key.size(), 10);
    fastgatk::pairhmm::FlowPairHmmHaplotype haplotype;
    haplotype.key = encoded.key;
    haplotype.flow_order = encoded.flow_order;
    return fastgatk::pairhmm::compute_kokkos_flow(
        {read}, {haplotype}, {{0, 0}}, 1);
}

void verify_flow_alignment_engine() {
    // GATK's normal `FlowBased` engine is not FlowBasedHMM: it selects the
    // highest direct, cycle-aligned placement.  Of the two legal placements
    // below, the second has four P=0.5 cells and must win over four P=0.01
    // cells in the first.
    fastgatk::pairhmm::FlowPairHmmRead read;
    read.key = {0, 1, 0, 1};
    read.flow_order = {'A', 'C', 'G', 'T'};
    read.probabilities.assign(read.key.size() * 256, 0.001);
    const auto set_probability = [&](const std::size_t flow, const std::size_t hmer,
                                     const double probability) {
        read.probabilities[flow * 256 + hmer] = probability;
    };
    set_probability(0, 0, 0.01); set_probability(0, 1, 0.5);
    set_probability(1, 1, 0.01); set_probability(1, 0, 0.5);
    set_probability(2, 0, 0.01); set_probability(2, 1, 0.5);
    set_probability(3, 1, 0.01); set_probability(3, 0, 0.5);
    fastgatk::pairhmm::FlowPairHmmHaplotype haplotype;
    haplotype.key = {0, 1, 0, 1, 1, 0, 1, 0};
    haplotype.flow_order = {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
    const auto result = fastgatk::pairhmm::compute_kokkos_flow_alignment(
        {read}, {haplotype}, {{0, 0}}, 1);
    const auto expected = 4.0 * std::log10(0.5);
    if (result.likelihoods.size() != 1 || result.execution_policy != "RangePolicy" ||
        std::abs(result.likelihoods.front() - expected) > 1e-12)
        throw std::runtime_error("GATK FlowBased direct-alignment mismatch");

    // The Java loop does not score partial-read placements.
    haplotype.key.resize(3);
    haplotype.flow_order.resize(3);
    const auto no_placement = fastgatk::pairhmm::compute_kokkos_flow_alignment(
        {read}, {haplotype}, {{0, 0}}, 1);
    if (no_placement.likelihoods.size() != 1 ||
        !std::isinf(no_placement.likelihoods.front()) ||
        no_placement.likelihoods.front() >= 0.0)
        throw std::runtime_error("GATK FlowBased direct-alignment accepted partial read");
}

}  // namespace

int main() {
    try {
        Kokkos::initialize();
        const std::vector<Case> cases{
            {"match", "ACGTACGT", "ACGTACGT", "ACGT", 0xbfefb6c8b0925800ULL},
            {"mismatch", "ACGTACGT", "AGGTACGT", "ACGT", 0xc0119a41e76fc8c0ULL},
            {"homopolymer", "AAACCCGGGTTT", "AACCCGGTTT", "ACGT", 0xc01aa25ea7075500ULL},
            {"phase", "TACGTACG", "TACGTTACG", "TGCA", 0xc011fb604fc6fcc0ULL},
        };
        double max_delta = 0.0;
        for (const auto& item : cases) {
            const auto result = score(item);
            if (result.likelihoods.size() != 1 || !std::isfinite(result.likelihoods[0]))
                throw std::runtime_error(std::string("invalid flow result for ") + item.name);
            const auto observed_bits = std::bit_cast<std::uint64_t>(result.likelihoods[0]);
            const auto expected_value = std::bit_cast<double>(item.expected_bits);
            const auto delta = std::abs(result.likelihoods[0] - expected_value);
            max_delta = std::max(max_delta, delta);
            if (observed_bits != item.expected_bits)
                throw std::runtime_error(std::string("GATK flow likelihood mismatch for ") + item.name +
                                         " observed_bits=0x" +
                                         [&] { std::ostringstream stream; stream << std::hex << observed_bits; return stream.str(); }() +
                                         " expected_bits=0x" +
                                         [&] { std::ostringstream stream; stream << std::hex << item.expected_bits; return stream.str(); }());
        }
        verify_flow_alignment_engine();
        Kokkos::finalize();
        std::cout << std::setprecision(17)
                  << "{\"status\":\"pass\",\"suite\":\"flow-pairhmm-gatk-oracle\",\"max_abs_delta\":"
                  << max_delta << "}\n";
        return 0;
    } catch (const std::exception& error) {
        if (Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
