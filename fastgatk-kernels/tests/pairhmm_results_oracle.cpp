#include "fastgatk/kernels/pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using fastgatk::pairhmm::PairHmmHaplotype;
using fastgatk::pairhmm::PairHmmPrecision;
using fastgatk::pairhmm::PairHmmRead;
using fastgatk::pairhmm::PairHmmRequest;

std::vector<std::uint8_t> decode_fastq_qualities(const std::string& qualities) {
    std::vector<std::uint8_t> decoded;
    decoded.reserve(qualities.size());
    for (const auto value : qualities) {
        if (static_cast<unsigned char>(value) < 33U)
            throw std::runtime_error("PairHMM results file contains an invalid FASTQ quality");
        decoded.push_back(static_cast<std::uint8_t>(
            static_cast<unsigned char>(value) - 33U));
    }
    return decoded;
}

struct ParsedResults {
    std::vector<PairHmmRead> reads;
    std::vector<PairHmmHaplotype> haplotypes;
    std::vector<PairHmmRequest> requests;
    std::vector<double> expected;
};

ParsedResults parse_results(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open PairHMM results file: " + path);
    ParsedResults parsed;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::istringstream fields(line);
        std::string haplotype;
        std::string read;
        std::string read_quality;
        std::string insertion_quality;
        std::string deletion_quality;
        std::string gap_continuation;
        std::string expected;
        if (!(fields >> haplotype >> read >> read_quality >> insertion_quality >>
              deletion_quality >> gap_continuation >> expected)) {
            throw std::runtime_error("malformed GATK PairHMM results row");
        }
        PairHmmRead read_record;
        read_record.bases.assign(read.begin(), read.end());
        read_record.qualities = decode_fastq_qualities(read_quality);
        read_record.insertion_gop = decode_fastq_qualities(insertion_quality);
        read_record.deletion_gop = decode_fastq_qualities(deletion_quality);
        read_record.gap_continuation = decode_fastq_qualities(gap_continuation);
        if (read_record.bases.empty() || read_record.bases.size() != read_record.qualities.size() ||
            read_record.bases.size() != read_record.insertion_gop.size() ||
            read_record.bases.size() != read_record.deletion_gop.size() ||
            read_record.bases.size() != read_record.gap_continuation.size())
            throw std::runtime_error("GATK PairHMM results row has inconsistent read lengths");
        PairHmmHaplotype haplotype_record;
        haplotype_record.bases.assign(haplotype.begin(), haplotype.end());
        if (haplotype_record.bases.empty())
            throw std::runtime_error("GATK PairHMM results row has an empty haplotype");
        const auto read_id = static_cast<std::uint32_t>(parsed.reads.size());
        const auto haplotype_id = static_cast<std::uint32_t>(parsed.haplotypes.size());
        parsed.reads.push_back(std::move(read_record));
        parsed.haplotypes.push_back(std::move(haplotype_record));
        parsed.requests.push_back(PairHmmRequest{read_id, haplotype_id});
        parsed.expected.push_back(std::stod(expected));
    }
    if (parsed.requests.empty())
        throw std::runtime_error("GATK PairHMM results file contains no rows");
    return parsed;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        if (argc != 2) {
            std::cerr << "usage: fastgatk-pairhmm-results-oracle GATK_PAIR_HMM_RESULTS_FILE\n";
            return 2;
        }
        const auto parsed = parse_results(argv[1]);
        Kokkos::initialize();
        initialized = true;
        const auto result = fastgatk::pairhmm::compute_kokkos_bucketed(
            parsed.reads, parsed.haplotypes, parsed.requests, 1, PairHmmPrecision::Float64);
        const auto float_result = fastgatk::pairhmm::compute_kokkos_bucketed(
            parsed.reads, parsed.haplotypes, parsed.requests, 1, PairHmmPrecision::Float32);
        if (result.likelihoods.size() != parsed.expected.size())
            throw std::runtime_error("PairHMM result row count changed");
        if (float_result.likelihoods.size() != parsed.expected.size())
            throw std::runtime_error("Float32 PairHMM result row count changed");
        constexpr double kPrintedScientificTolerance = 5.0e-6;
        double max_abs_delta = 0.0;
        double sum_abs_delta = 0.0;
        double float_max_abs_delta = 0.0;
        double float_sum_abs_delta = 0.0;
        std::size_t within_tolerance = 0;
        std::size_t float_within_tolerance = 0;
        for (std::size_t index = 0; index < parsed.expected.size(); ++index) {
            const auto delta = std::abs(result.likelihoods[index] - parsed.expected[index]);
            const auto float_delta = std::abs(float_result.likelihoods[index] - parsed.expected[index]);
            max_abs_delta = std::max(max_abs_delta, delta);
            sum_abs_delta += delta;
            float_max_abs_delta = std::max(float_max_abs_delta, float_delta);
            float_sum_abs_delta += float_delta;
            if (delta <= kPrintedScientificTolerance) ++within_tolerance;
            if (float_delta <= kPrintedScientificTolerance) ++float_within_tolerance;
        }
        // PairHMM's GATK debug writer serializes `%e` with six digits after the
        // decimal point.  This is a decimal-observation bound, not a relaxed
        // model oracle: it is the only error permitted when replaying those
        // exact read/haplotype/quality rows through the shared Kokkos kernel.
        if (within_tolerance != parsed.expected.size())
            throw std::runtime_error("Kokkos PairHMM differs from a GATK debug row beyond the "
                                     "fixed six-decimal serialization bound");
        // GATK's default native PairHMM is GKL's Float32 path with its own
        // low-scaled-sum Context<double> recomputation.  Verify that public
        // default separately: otherwise a Float32 underflow can be hidden by
        // the opt-in Float64/Java compatibility replay above.
        if (float_within_tolerance != parsed.expected.size())
            throw std::runtime_error("Kokkos GKL Float32 PairHMM differs from a GATK debug row beyond the "
                                     "fixed six-decimal serialization bound");
        Kokkos::finalize();
        initialized = false;
        std::cout << std::setprecision(12)
                  << "{\"status\":\"pass\",\"rows\":" << parsed.expected.size()
                  << ",\"max_abs_delta\":" << max_abs_delta
                  << ",\"mean_abs_delta\":" << sum_abs_delta / parsed.expected.size()
                  << ",\"float32_max_abs_delta\":" << float_max_abs_delta
                  << ",\"float32_mean_abs_delta\":" << float_sum_abs_delta / parsed.expected.size()
                  << ",\"float32_within_decimal_bound\":" << float_within_tolerance
                  << ",\"decimal_bound\":" << kPrintedScientificTolerance << "}\n";
        return 0;
    } catch (const std::exception& error) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
