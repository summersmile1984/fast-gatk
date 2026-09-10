#include "pairhmm_kokkos.hpp"

#include <Kokkos_Core.hpp>

#include <bit>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using fastgatk::pairhmm::PairIndexBatch;
using fastgatk::pairhmm::PairInput;

namespace {

struct Options {
    std::size_t records = 512;
    std::size_t read = 150;
    std::size_t hap = 160;
    int threads = 1;
    int iterations = 10;
    std::uint64_t seed = 42;
    std::string workload = "independent";
    std::string values_out;
    std::string sums_out;
};

struct SplitMix64 {
    std::uint64_t state;
    std::uint64_t next() {
        state += 0x9e3779b97f4a7c15ULL;
        std::uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
};

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        auto value = [&](const char* name) -> std::string {
            const std::string prefix = std::string(name) + "=";
            return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
        };
        if (!value("--pairs").empty()) options.records = std::stoull(value("--pairs"));
        else if (!value("--read-len").empty()) options.read = std::stoull(value("--read-len"));
        else if (!value("--hap-len").empty()) options.hap = std::stoull(value("--hap-len"));
        else if (!value("--threads").empty()) options.threads = std::stoi(value("--threads"));
        else if (!value("--iterations").empty()) options.iterations = std::stoi(value("--iterations"));
        else if (!value("--seed").empty()) options.seed = std::stoull(value("--seed"));
        else if (!value("--workload").empty()) options.workload = value("--workload");
        else if (!value("--values-out").empty()) options.values_out = value("--values-out");
        else if (!value("--sums-out").empty()) options.sums_out = value("--sums-out");
        else if (argument == "--help") {
            std::cout << "pairhmm-kokkos --pairs=N --read-len=N --hap-len=N --threads=N "
                         "--iterations=N --seed=N --workload=independent|matrix8\n";
            std::exit(0);
        } else throw std::runtime_error("unknown option: " + argument);
    }
    if (options.threads < 1 || options.records < 1 || options.read < 1 || options.hap < 1)
        throw std::runtime_error("sizes and thread count must be positive");
    return options;
}

std::vector<PairInput> make_records(const Options& options) {
    SplitMix64 rng{options.seed};
    static constexpr std::uint8_t bases[] = {'A', 'C', 'G', 'T'};
    std::vector<PairInput> records(options.records);
    for (PairInput& input : records) {
        input.haplotype.resize(options.hap);
        input.read.resize(options.read);
        input.read_qual.resize(options.read);
        input.insertion_gop.resize(options.read);
        input.deletion_gop.resize(options.read);
        input.gap_continuation.resize(options.read);
        for (std::uint8_t& base : input.haplotype) base = bases[rng.next() % 4];
        for (std::size_t i = 0; i < options.read; ++i) {
            input.read[i] = (rng.next() % 10 == 0) ? bases[rng.next() % 4]
                                                   : input.haplotype[i % options.hap];
            input.read_qual[i] = static_cast<std::uint8_t>(25 + rng.next() % 16);
            input.insertion_gop[i] = static_cast<std::uint8_t>(35 + rng.next() % 8);
            input.deletion_gop[i] = static_cast<std::uint8_t>(35 + rng.next() % 8);
            input.gap_continuation[i] = static_cast<std::uint8_t>(10 + rng.next() % 10);
        }
    }
    return records;
}

PairIndexBatch make_pairs(const Options& options) {
    PairIndexBatch pairs;
    if (options.workload == "independent") {
        pairs.read_ids.resize(options.records);
        pairs.haplotype_ids.resize(options.records);
        for (std::size_t i = 0; i < options.records; ++i) {
            pairs.read_ids[i] = static_cast<std::uint32_t>(i);
            pairs.haplotype_ids[i] = static_cast<std::uint32_t>(i);
        }
        return pairs;
    }
    if (options.workload != "matrix8")
        throw std::runtime_error("unknown workload: " + options.workload);
    constexpr std::size_t block = 8;
    for (std::size_t base = 0; base < options.records; base += block) {
        const std::size_t count = std::min(block, options.records - base);
        for (std::size_t read = 0; read < count; ++read) {
            for (std::size_t hap = 0; hap < count; ++hap) {
                pairs.read_ids.push_back(static_cast<std::uint32_t>(base + read));
                pairs.haplotype_ids.push_back(static_cast<std::uint32_t>(base + hap));
            }
        }
    }
    return pairs;
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        const Options options = parse(argc, argv);
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(options.threads);
        Kokkos::initialize(settings);
        initialized = true;

        const auto records = make_records(options);
        const PairIndexBatch pairs = make_pairs(options);
        const auto result = fastgatk::pairhmm::compute_kokkos(records, pairs, options.iterations);
        std::cout << "# kokkos_backend=" << fastgatk::pairhmm::kokkos_backend_description()
                  << " threads=" << options.threads << '\n';
        std::cout << "{\"mode\":\"kokkos\",\"workload\":\"" << options.workload
                  << "\",\"records\":" << options.records
                  << ",\"computed_pairs\":" << pairs.read_ids.size()
                  << ",\"read_len\":" << options.read << ",\"hap_len\":" << options.hap
                  << ",\"threads\":" << options.threads << ",\"iterations\":" << options.iterations
                  << ",\"execution_space\":\"" << result.execution_space << "\""
                  << ",\"api\":\"Kokkos::Experimental::simd<double>\""
                  << ",\"simd_width\":" << result.simd_width
                  << ",\"prepare_seconds\":" << std::setprecision(12) << result.prepare_seconds
                  << ",\"seconds\":" << result.seconds
                  << ",\"pairs_per_second\":" << result.pairs_per_second
                  << ",\"checksum\":" << result.checksum << "}\n";
        if (!options.values_out.empty()) {
            std::ofstream output(options.values_out);
            if (!output) throw std::runtime_error("cannot open values output: " + options.values_out);
            for (double value : result.likelihoods)
                output << std::hex << std::bit_cast<std::uint64_t>(value) << '\n';
        }
        if (!options.sums_out.empty()) {
            std::ofstream output(options.sums_out);
            if (!output) throw std::runtime_error("cannot open sums output: " + options.sums_out);
            for (double value : result.scaled_sums)
                output << std::hex << std::bit_cast<std::uint64_t>(value) << '\n';
        }
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        return 2;
    }
}
