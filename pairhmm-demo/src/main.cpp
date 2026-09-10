#include "pairhmm.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using fastgatk::pairhmm::PairInput;

namespace {
struct Options { std::string mode = "both"; std::string backend = "auto"; std::size_t pairs = 256; std::size_t read = 150; std::size_t hap = 160; int threads = 1; int iterations = 5; std::uint64_t seed = 42; std::string values_out; };

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
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a(argv[i]);
        auto value = [&](const char* name) -> std::string {
            const std::string p = std::string(name) + "=";
            if (a.rfind(p, 0) != 0) return {};
            return a.substr(p.size());
        };
        if (!value("--mode").empty()) o.mode = value("--mode");
        else if (!value("--backend").empty()) o.backend = value("--backend");
        else if (!value("--pairs").empty()) o.pairs = std::stoull(value("--pairs"));
        else if (!value("--read-len").empty()) o.read = std::stoull(value("--read-len"));
        else if (!value("--hap-len").empty()) o.hap = std::stoull(value("--hap-len"));
        else if (!value("--threads").empty()) o.threads = std::stoi(value("--threads"));
        else if (!value("--iterations").empty()) o.iterations = std::stoi(value("--iterations"));
        else if (!value("--seed").empty()) o.seed = std::stoull(value("--seed"));
        else if (!value("--values-out").empty()) o.values_out = value("--values-out");
        else if (a == "--help") { std::cout << "pairhmm-demo --mode=scalar|simd|both --backend=auto|avx2|avx512 --pairs=N --read-len=N --hap-len=N --threads=N --iterations=N --seed=N\n"; std::exit(0); }
        else throw std::runtime_error("unknown option: " + a);
    }
    return o;
}

fastgatk::pairhmm::SimdBackend parse_backend(const std::string& backend) {
    if (backend == "auto") return fastgatk::pairhmm::SimdBackend::Auto;
    if (backend == "avx2") return fastgatk::pairhmm::SimdBackend::Avx2;
    if (backend == "avx512") return fastgatk::pairhmm::SimdBackend::Avx512;
    throw std::runtime_error("unknown backend: " + backend);
}

std::vector<PairInput> make_inputs(const Options& o) {
    SplitMix64 rng{o.seed};
    static constexpr std::uint8_t bases[] = {'A','C','G','T'};
    std::vector<PairInput> out(o.pairs);
    for (auto& p : out) {
        p.haplotype.resize(o.hap); p.read.resize(o.read);
        p.read_qual.resize(o.read); p.insertion_gop.resize(o.read); p.deletion_gop.resize(o.read); p.gap_continuation.resize(o.read);
        for (auto& b : p.haplotype) b = bases[rng.next() % 4];
        for (std::size_t i = 0; i < o.read; ++i) {
            p.read[i] = (rng.next() % 10 == 0) ? bases[rng.next() % 4] : p.haplotype[i % o.hap];
            p.read_qual[i] = static_cast<std::uint8_t>(25 + rng.next() % 16);
            p.insertion_gop[i] = static_cast<std::uint8_t>(35 + rng.next() % 8);
            p.deletion_gop[i] = static_cast<std::uint8_t>(35 + rng.next() % 8);
            p.gap_continuation[i] = static_cast<std::uint8_t>(10 + rng.next() % 10);
        }
    }
    return out;
}

void print_result(const char* name, const fastgatk::pairhmm::BatchResult& r, const Options& o) {
    std::cout << "{\"mode\":\"" << name << "\",\"pairs\":" << o.pairs << ",\"read_len\":" << o.read
              << ",\"hap_len\":" << o.hap << ",\"threads\":" << o.threads << ",\"iterations\":" << o.iterations
              << ",\"seconds\":" << std::setprecision(12) << r.seconds << ",\"pairs_per_second\":" << r.pairs_per_second
              << ",\"checksum\":" << r.checksum << "}\n";
}
}

int main(int argc, char** argv) {
    try {
        const Options o = parse(argc, argv);
        const auto inputs = make_inputs(o);
        std::cout << "# cpu_simd=" << fastgatk::pairhmm::cpu_simd_description() << " threads=" << o.threads << "\n";
        fastgatk::pairhmm::BatchResult scalar, simd;
        if (o.mode == "scalar" || o.mode == "both") { scalar = fastgatk::pairhmm::compute_scalar(inputs, o.threads, o.iterations); print_result("scalar", scalar, o); }
        if (o.mode == "simd" || o.mode == "both") {
            simd = fastgatk::pairhmm::compute_simd(inputs, o.threads, o.iterations, parse_backend(o.backend));
            const std::string simd_name = "simd-" + o.backend;
            print_result(simd_name.c_str(), simd, o);
        }
        if (!o.values_out.empty()) {
            const auto& values = (o.mode == "scalar") ? scalar.likelihoods : simd.likelihoods;
            std::ofstream out(o.values_out);
            if (!out) throw std::runtime_error("cannot open values output: " + o.values_out);
            for (double x : values) out << std::hex << std::bit_cast<std::uint64_t>(x) << '\n';
        }
        if (o.mode == "both") {
            double max_abs = 0.0; std::size_t differing = 0;
            for (std::size_t i = 0; i < inputs.size(); ++i) { const double d = std::abs(scalar.likelihoods[i] - simd.likelihoods[i]); max_abs = std::max(max_abs, d); if (std::bit_cast<std::uint64_t>(scalar.likelihoods[i]) != std::bit_cast<std::uint64_t>(simd.likelihoods[i])) ++differing; }
            std::cout << "{\"compatibility\":{\"max_abs\":" << std::setprecision(12) << max_abs << ",\"bit_different\":" << differing << "}}\n";
        }
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << '\n'; return 2; }
}
