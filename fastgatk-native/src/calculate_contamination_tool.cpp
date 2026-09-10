#include "fastgatk/runtime/resource.hpp"
#include "fastgatk/core/plan.hpp"
#include "fastgatk/contamination_kernel_segmenter.hpp"
#include "fastgatk/io/java_numeric.hpp"

#include <Kokkos_Core.hpp>

#include <array>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef FASTGATK_HAS_ZLIB
#define FASTGATK_HAS_ZLIB 0
#endif

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string matched_normal;
    std::string output;
    std::string segmentation;
    std::string manifest;
    double low_coverage_ratio = 0.5;
    double high_coverage_ratio = 3.0;
    int threads = 1;
};

std::string format_double(double value) {
    return fastgatk::io::java_double(value);
}

// JSON has no NaN literal.  Empty post-coverage panels are valid GATK inputs,
// but their median/mean thresholds are undefined; preserve that fact as null
// instead of producing an unreadable output-manifest sidecar.
std::string json_number(double value) {
    return std::isfinite(value) ? format_double(value) : "null";
}

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& argument, const char* name) {
    return argument == name || !option_value(argument, name).empty();
}

std::string require_value(int& index, int argc, char** argv,
                          const std::string& argument, const char* name,
                          const char* short_name = nullptr) {
    const auto inline_value = option_value(argument, name);
    if (!inline_value.empty()) return inline_value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + name);
}

Options parse(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-calculate-contamination (GATK-compatible native prototype)\n"
                         "  -I, --input FILE                    tumor pileup summary table\n"
                         "      --matched-normal, --matched FILE matched normal pileup table\n"
                         "  -O, --output FILE                   contamination table\n"
                         "      --tumor-segmentation FILE       optional MAF segmentation table\n"
                         "      --low-coverage-ratio-threshold F (default 0.5)\n"
                         "      --high-coverage-ratio-threshold F (default 3.0)\n"
                         "      --output-manifest FILE          OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-I" || is_option(argument, "--input"))
            options.input = require_value(index, argc, argv, argument, "--input", "-I");
        else if (is_option(argument, "--matched-normal") || argument == "--matched" || argument == "-matched")
            options.matched_normal = require_value(
                index, argc, argv, argument,
                argument == "--matched" || argument == "-matched" ||
                !option_value(argument, "--matched").empty()
                    ? "--matched" : "--matched-normal", "-matched");
        else if (argument == "-O" || is_option(argument, "--output"))
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        else if (is_option(argument, "--tumor-segmentation") || argument == "--segments" || argument == "-segments") {
            const bool legacy_segments = argument == "--segments" || argument == "-segments";
            options.segmentation = require_value(
                index, argc, argv, argument,
                legacy_segments ? "--segments" : "--tumor-segmentation",
                legacy_segments ? "-segments" : nullptr);
        }
        else if (is_option(argument, "--low-coverage-ratio-threshold"))
            options.low_coverage_ratio = std::stod(require_value(
                index, argc, argv, argument, "--low-coverage-ratio-threshold"));
        else if (is_option(argument, "--high-coverage-ratio-threshold"))
            options.high_coverage_ratio = std::stod(require_value(
                index, argc, argv, argument, "--high-coverage-ratio-threshold"));
        else if (is_option(argument, "--threads"))
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        else if (is_option(argument, "--output-manifest") || is_option(argument, "--manifest"))
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        else if (argument == "--quiet" || argument == "--disable-sequence-dictionary-validation") {
            // Compatibility flags; table columns are validated by the native reader.
        } else if (is_option(argument, "--java-options") || is_option(argument, "--verbosity") ||
                   is_option(argument, "--seconds-between-progress-updates")) {
            if (argument.find('=') == std::string::npos)
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.input.empty()) throw std::invalid_argument("-I/--input is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (!std::isfinite(options.low_coverage_ratio) || options.low_coverage_ratio < 0.0 ||
        !std::isfinite(options.high_coverage_ratio) || options.high_coverage_ratio <= 0.0 ||
        options.low_coverage_ratio >= options.high_coverage_ratio)
        throw std::invalid_argument("coverage thresholds must satisfy 0 <= low < high");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    return options;
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

bool file_complete(const std::string& path) {
    std::error_code error;
    return std::filesystem::is_regular_file(path, error) && !error &&
           std::filesystem::file_size(path, error) > 0 && !error;
}

bool compressed_path(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

template <typename Fn>
void for_each_text_line(const std::string& path, Fn&& callback) {
    if (!compressed_path(path)) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open pileup table: " + path);
        for (std::string line; std::getline(input, line);) callback(std::move(line));
        if (input.bad()) throw std::runtime_error("BAD_INPUT: cannot read pileup table: " + path);
        return;
    }
#if FASTGATK_HAS_ZLIB
    gzFile input = gzopen(path.c_str(), "rb");
    if (!input) throw std::runtime_error("BAD_INPUT: cannot open compressed pileup table: " + path);
    char buffer[1 << 16];
    std::string pending;
    while (true) {
        char* result = gzgets(input, buffer, static_cast<int>(sizeof(buffer)));
        if (result == nullptr) break;
        pending.append(result);
        if (pending.find('\n') == std::string::npos) continue;
        std::size_t begin = 0;
        while (true) {
            const auto end = pending.find('\n', begin);
            if (end == std::string::npos) {
                pending = pending.substr(begin);
                break;
            }
            callback(pending.substr(begin, end - begin));
            begin = end + 1;
        }
    }
    if (!pending.empty()) callback(std::move(pending));
    int error_code = Z_OK;
    (void)gzerror(input, &error_code);
    const auto close_status = gzclose(input);
    if (error_code != Z_OK || close_status != Z_OK)
        throw std::runtime_error("BAD_INPUT: cannot read compressed pileup table: " + path);
#else
    throw std::runtime_error("BACKEND_UNAVAILABLE: compressed pileup tables require zlib: " + path);
#endif
}

void write_text_output(const std::string& path, const std::string& text) {
    if (!compressed_path(path)) {
        std::ofstream output(path, std::ios::binary);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open output: " + path);
        output.write(text.data(), static_cast<std::streamsize>(text.size()));
        output.close();
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: incomplete output: " + path);
        return;
    }
#if FASTGATK_HAS_ZLIB
    gzFile output = gzopen(path.c_str(), "wb");
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot open compressed output: " + path);
    std::size_t offset = 0;
    while (offset < text.size()) {
        const auto chunk = static_cast<unsigned>(std::min<std::size_t>(text.size() - offset, 1U << 20));
        if (gzwrite(output, text.data() + offset, chunk) != static_cast<int>(chunk)) {
            (void)gzclose(output);
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write compressed output: " + path);
        }
        offset += chunk;
    }
    if (gzclose(output) != Z_OK)
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize compressed output: " + path);
#else
    throw std::runtime_error("BACKEND_UNAVAILABLE: compressed outputs require zlib: " + path);
#endif
}

struct PileupSite {
    std::string contig;
    int position = 0;
    int ref_count = 0;
    int alt_count = 0;
    int other_count = 0;
    double allele_frequency = 0.0;

    int total() const { return ref_count + alt_count + other_count; }
    double alt_fraction() const { return total() == 0 ? 0.0 : static_cast<double>(alt_count) / total(); }
    double minor_allele_fraction() const {
        const auto fraction = alt_fraction();
        return std::min(fraction, 1.0 - fraction);
    }
};

struct PileupTable {
    std::string sample = "UNKNOWN";
    std::vector<PileupSite> sites;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const auto end = line.find('\t', start);
        fields.push_back(line.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return fields;
}

PileupTable read_table(const std::string& path) {
    PileupTable table;
    bool header_seen = false;
    std::size_t line_number = 0;
    for_each_text_line(path, [&](std::string line) {
        ++line_number;
        if (line.empty()) return;
        if (line.rfind("#<METADATA>SAMPLE=", 0) == 0) {
            table.sample = line.substr(std::string("#<METADATA>SAMPLE=").size());
            return;
        }
        if (line[0] == '#') return;
        const auto fields = split_tab(line);
        if (!header_seen) {
            if (fields.size() != 6 || fields[0] != "contig" || fields[1] != "position")
                throw std::runtime_error("BAD_INPUT: invalid pileup table header at line " + std::to_string(line_number));
            header_seen = true;
            return;
        }
        if (fields.size() != 6) throw std::runtime_error("BAD_INPUT: invalid pileup table row at line " + std::to_string(line_number));
        PileupSite site;
        site.contig = fields[0];
        site.position = std::stoi(fields[1]);
        site.ref_count = std::stoi(fields[2]);
        site.alt_count = std::stoi(fields[3]);
        site.other_count = std::stoi(fields[4]);
        site.allele_frequency = std::stod(fields[5]);
        if (site.position < 1 || site.ref_count < 0 || site.alt_count < 0 || site.other_count < 0 ||
            !std::isfinite(site.allele_frequency) || site.allele_frequency < 0.0 || site.allele_frequency > 1.0)
            throw std::runtime_error("BAD_INPUT: invalid pileup table values at line " + std::to_string(line_number));
        table.sites.push_back(std::move(site));
    });
    if (!header_seen || table.sites.empty()) throw std::runtime_error("BAD_INPUT: pileup table has no records: " + path);
    return table;
}

// Apache Commons Math's Mean.evaluate(double[]) (used by GATK's
// ContaminationModel) is a corrected two-pass mean.  It is intentionally not
// replaced by a compensated summation: the first pass is Java's ordinary
// left-to-right Sum, followed by a residual correction pass.  Keeping that
// operation explicit gives the native path the same threshold behavior at
// strict high-coverage boundaries.
double java_mean(const std::vector<int>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double sample_size = static_cast<double>(values.size());
    double sum = 0.0;
    for (const int value : values) sum += static_cast<double>(value);
    const double xbar = sum / sample_size;
    double correction = 0.0;
    for (const int value : values) correction += static_cast<double>(value) - xbar;
    return xbar + correction / sample_size;
}

struct CoverageFilterResult {
    std::vector<PileupSite> sites;
    std::size_t covered_sites = 0;
    double median = std::numeric_limits<double>::quiet_NaN();
    double mean = std::numeric_limits<double>::quiet_NaN();
    double low = std::numeric_limits<double>::quiet_NaN();
    double high = std::numeric_limits<double>::quiet_NaN();
};

CoverageFilterResult coverage_filter(const std::vector<PileupSite>& sites,
                                     const Options& options) {
    std::vector<int> covered;
    for (const auto& site : sites) if (site.total() > 10) covered.push_back(site.total());
    CoverageFilterResult result;
    result.covered_sites = covered.size();
    if (covered.empty()) return result;
    std::sort(covered.begin(), covered.end());
    const double median = covered.size() % 2 == 0
        ? 0.5 * (static_cast<double>(covered[covered.size() / 2 - 1]) +
                 static_cast<double>(covered[covered.size() / 2]))
        : static_cast<double>(covered[covered.size() / 2]);
    const double mean = java_mean(covered);
    const double low = median * options.low_coverage_ratio;
    const double high = mean * options.high_coverage_ratio;
    result.median = median;
    result.mean = mean;
    result.low = low;
    result.high = high;
    for (const auto& site : sites)
        if (site.total() > 10 && site.total() > low && site.total() < high)
            result.sites.push_back(site);
    return result;
}

std::string site_key(const PileupSite& site) {
    return site.contig + ":" + std::to_string(site.position);
}

struct Estimate {
    double contamination = 0.0;
    double error = 1.0;
    std::size_t candidate_sites = 0;
    std::uint64_t total_depth = 0;
    std::uint64_t opposite_depth = 0;
};

// This is the C++/Kokkos implementation of GATK's ContaminationModel.  The
// model is deliberately kept separate from the older hom-site grid estimator
// below while it is being validated: it makes the Java control flow (three
// MAF/contamination rounds followed by the hom-alt -> hom-ref -> heuristic
// strategy cascade) explicit and gives us a small, testable Kokkos batch
// likelihood boundary.
struct ContaminationSegment {
    std::vector<std::size_t> sites;
    double minor_allele_fraction = 0.5;
    std::size_t plan_index = 0;
};

struct LearnedContaminationModel {
    std::vector<ContaminationSegment> segments;
    double error_rate = 0.0;
    double contamination = 0.0;
};

struct FullEstimate {
    double contamination = 0.0;
    double error = 1.0;
    std::size_t candidate_sites = 0;
    std::uint64_t total_depth = 0;
    std::uint64_t opposite_depth = 0;
    std::vector<ContaminationSegment> segments;
};

double java_round(double value) {
    return std::floor(value + 0.5);
}

double binomial_probability_without_constant(int n, int k, double p) {
    // The omitted n-choose-k factor is independent of genotype and
    // contamination, so it cancels in posterior ratios and does not alter
    // the optimizer.  Avoiding lgamma also keeps this function available in
    // a Kokkos device lambda on CUDA/HIP/SYCL.
    constexpr double epsilon = 1.0e-15;
    if (n < 0 || k < 0 || k > n) return 0.0;
    const double bounded = std::clamp(p, epsilon, 1.0 - epsilon);
    return std::exp(static_cast<double>(k) * std::log(bounded) +
                    static_cast<double>(n - k) * std::log1p(-bounded));
}

std::array<double, 4> genotype_likelihoods_host(const PileupSite& site,
                                                double contamination,
                                                double error_rate,
                                                double maf) {
    const double f = site.allele_frequency;
    const int n = site.ref_count + site.alt_count;
    const int k = site.alt_count;
    const std::array<double, 4> priors{
        (1.0 - f) * (1.0 - f), f * (1.0 - f), f * (1.0 - f), f * f};
    const std::array<double, 4> allele_fractions{
        error_rate / 3.0, maf, 1.0 - maf, 1.0 - error_rate};
    std::array<double, 4> likelihoods{};
    for (std::size_t genotype = 0; genotype < likelihoods.size(); ++genotype) {
        const double p = (1.0 - contamination) * allele_fractions[genotype] +
                         contamination * f;
        likelihoods[genotype] = priors[genotype] *
            binomial_probability_without_constant(n, k, p);
    }
    return likelihoods;
}

double plain_sum(const std::array<double, 4>& values) {
    double result = 0.0;
    for (const double value : values) result += value;
    return result;
}

// JDK DoubleStream.sum() uses the compensated three-slot accumulator rather
// than a plain left-to-right add.  The distinction is only a few ULPs on the
// contamination fixtures, but it can change the final printed estimate and
// therefore belongs in the strict Java-compatible Host boundary.
struct JavaDoubleSum {
    double sum = 0.0;
    double compensation = 0.0;
    double simple_sum = 0.0;

    void add(double value) {
        const double corrected = value - compensation;
        const double next = sum + corrected;
        compensation = (next - sum) - corrected;
        sum = next;
        simple_sum += value;
    }

    double value() const {
        const double result = sum - compensation;
        return (!std::isfinite(result) && std::isfinite(simple_sum)) ? simple_sum : result;
    }
};

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t observations = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::uint64_t plan_allocations = 0;
    std::uint64_t plan_reuses = 0;
};

struct ContaminationLikelihoodPlan {
    Kokkos::View<int*> refs;
    Kokkos::View<int*> alts;
    Kokkos::View<double*> frequencies;
    Kokkos::View<double*> terms;
    std::vector<Kokkos::View<std::size_t*>> segment_indices;
    KernelStats* stats = nullptr;
    mutable std::uint64_t evaluations = 0;

    ContaminationLikelihoodPlan(const std::vector<PileupSite>& sites,
                                const std::vector<ContaminationSegment>& segments,
                                KernelStats* kernel_stats = nullptr)
        : refs("contamination_model_refs", sites.size()),
          alts("contamination_model_alts", sites.size()),
          frequencies("contamination_model_af", sites.size()),
          terms("contamination_model_terms", sites.size()),
          stats(kernel_stats) {
        auto host_refs = Kokkos::create_mirror_view(refs);
        auto host_alts = Kokkos::create_mirror_view(alts);
        auto host_frequencies = Kokkos::create_mirror_view(frequencies);
        for (std::size_t index = 0; index < sites.size(); ++index) {
            host_refs(index) = sites[index].ref_count;
            host_alts(index) = sites[index].alt_count;
            host_frequencies(index) = sites[index].allele_frequency;
        }
        Kokkos::deep_copy(refs, host_refs);
        Kokkos::deep_copy(alts, host_alts);
        Kokkos::deep_copy(frequencies, host_frequencies);
        segment_indices.reserve(segments.size());
        for (const auto& segment : segments) {
            Kokkos::View<std::size_t*> indices(
                "contamination_model_segment_indices", segment.sites.size());
            auto host_indices = Kokkos::create_mirror_view(indices);
            for (std::size_t index = 0; index < segment.sites.size(); ++index)
                host_indices(index) = segment.sites[index];
            Kokkos::deep_copy(indices, host_indices);
            segment_indices.push_back(std::move(indices));
        }
        if (stats != nullptr) ++stats->plan_allocations;
    }

    double segment_log_likelihood(std::size_t segment_index,
                                  double contamination,
                                  double error_rate,
                                  double maf) const {
        if (segment_index >= segment_indices.size()) return 0.0;
        if (stats != nullptr) {
            if (evaluations > 0) ++stats->plan_reuses;
            ++evaluations;
        }
        const auto indices = segment_indices[segment_index];
        const std::size_t count = indices.extent(0);
        if (count == 0) return 0.0;
        const auto refs_view = refs;
        const auto alts_view = alts;
        const auto frequencies_view = frequencies;
        const auto terms_view = terms;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host("contamination-likelihood-v1");
        host.records = count;
        host.bytes = count * (sizeof(std::size_t) + sizeof(int) * 2 + sizeof(double) * 2);
        fastgatk::core::KernelPlan<ExecSpace> kernel("contamination-model-likelihood");
        kernel.begin_prepare(host);
        fastgatk::core::DeviceBatch<ExecSpace> device(count);
        device.bind("indices", indices);
        device.bind("refs", refs_view);
        device.bind("alts", alts_view);
        device.bind("frequencies", frequencies_view);
        device.bind("terms", terms_view);
        ExecSpace().fence();
        kernel.end_prepare(device);
        kernel.begin_execute();
        Kokkos::parallel_for("contamination_model_likelihood", Kokkos::RangePolicy<ExecSpace>(0, count),
        KOKKOS_LAMBDA(const std::size_t index) {
            const std::size_t site_index = indices(index);
            const double f = frequencies_view(site_index);
            const int n = refs_view(site_index) + alts_view(site_index);
            const int k = alts_view(site_index);
            const double p0 = (1.0 - contamination) * (error_rate / 3.0) + contamination * f;
            const double p1 = (1.0 - contamination) * maf + contamination * f;
            const double p2 = (1.0 - contamination) * (1.0 - maf) + contamination * f;
            const double p3 = (1.0 - contamination) * (1.0 - error_rate) + contamination * f;
            const double l0 = (1.0 - f) * (1.0 - f) *
                Kokkos::exp(static_cast<double>(k) * Kokkos::log(Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p0))) +
                             static_cast<double>(n - k) * Kokkos::log1p(-Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p0))));
            const double l1 = f * (1.0 - f) *
                Kokkos::exp(static_cast<double>(k) * Kokkos::log(Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p1))) +
                             static_cast<double>(n - k) * Kokkos::log1p(-Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p1))));
            const double l2 = f * (1.0 - f) *
                Kokkos::exp(static_cast<double>(k) * Kokkos::log(Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p2))) +
                             static_cast<double>(n - k) * Kokkos::log1p(-Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p2))));
            const double l3 = f * f *
                Kokkos::exp(static_cast<double>(k) * Kokkos::log(Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p3))) +
                             static_cast<double>(n - k) * Kokkos::log1p(-Kokkos::max(1.0e-15, Kokkos::min(1.0 - 1.0e-15, p3))));
            terms_view(index) = Kokkos::log(Kokkos::max(1.0e-300, l0 + l1 + l2 + l3));
        });
        ExecSpace().fence();
        kernel.end_execute();
        if (stats != nullptr) {
            ++stats->batches;
            stats->observations += count;
            stats->prepare_seconds += kernel.telemetry().prepare_seconds;
            stats->execute_seconds += kernel.telemetry().execute_seconds;
        }
        auto host_terms = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, terms);
        JavaDoubleSum result;
        for (std::size_t index = 0; index < count; ++index) result.add(host_terms(index));
        return result.value();
    }
};

double model_log_likelihood(const ContaminationLikelihoodPlan& plan,
                            const std::vector<ContaminationSegment>& segments,
                            double contamination, double error_rate) {
    double result = 0.0;
    for (const auto& segment : segments)
        result += plan.segment_log_likelihood(segment.plan_index, contamination,
                                              error_rate, segment.minor_allele_fraction);
    return result;
}

// Apache Commons Math's BrentOptimizer (as used by GATK's
// OptimizationUtils.max), specialized for maximization.  Keeping the
// interval/guess/tolerance/evaluation contract identical to Java matters for
// CalculateContamination: the MAF fit is run three times and the contamination
// fit is started from four fixed guesses, so a different optimizer can select a
// different local maximum even when the likelihood itself is equivalent.
template <typename Function>
double brent_maximize(Function&& function, double lower, double upper,
                      double guess, double relative_tolerance,
                      double absolute_tolerance, int max_evaluations) {
    if (!(upper > lower)) return lower;
    const double x0 = std::clamp(guess, lower, upper);
    constexpr double golden_section = 0.3819660112501051518;
    constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();
    constexpr double minimum_absolute_tolerance = 1.0e-15;
    const double relative = std::max(relative_tolerance,
                                     2.0 * machine_epsilon);
    const double absolute = std::max(absolute_tolerance,
                                     minimum_absolute_tolerance);
    double a = lower;
    double b = upper;
    double x = x0;
    double w = x;
    double v = x;
    double fx = function(x);
    double fw = fx;
    double fv = fx;
    double d = 0.0;
    double e = 0.0;
    int evaluations = 1;
    while (evaluations < max_evaluations) {
        const double midpoint = 0.5 * (a + b);
        // Commons Math uses the absolute threshold directly here (the
        // `/3` form belongs to a different Brent implementation).  Matching
        // this detail is important because GATK's MAF optimizer intentionally
        // stops at a relatively loose 0.01 absolute tolerance.
        const double tol1 = relative * std::abs(x) + absolute;
        const double tol2 = 2.0 * tol1;
        if (std::abs(x - midpoint) <= tol2 - 0.5 * (b - a)) break;

        bool parabolic = false;
        if (std::abs(e) > tol1) {
            const double r = (x - w) * (fx - fv);
            const double q0 = (x - v) * (fx - fw);
            double p = (x - v) * q0 - (x - w) * r;
            double q = 2.0 * (q0 - r);
            if (q > 0.0) p = -p;
            else q = -q;
            const double previous_e = e;
            e = d;
            if (q != 0.0 && std::abs(p) < std::abs(0.5 * q * previous_e) &&
                p > q * (a - x) && p < q * (b - x)) {
                d = p / q;
                const double u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = std::copysign(tol1, midpoint - x);
                parabolic = true;
            }
        }
        if (!parabolic) {
            e = x < midpoint ? b - x : a - x;
            d = golden_section * e;
        }
        const double u = x + (std::abs(d) >= tol1
            ? d : std::copysign(tol1, d));
        const double fu = function(u);
        ++evaluations;
        // GoalType.MAXIMIZE: retain the higher-valued point.  The tie rule is
        // intentionally the same as BrentOptimizer's strict comparison, so a
        // flat objective keeps the earlier (Java) point.
        if (fu > fx) {
            if (u < x) b = x;
            else a = x;
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u;
            else b = u;
            if (fu > fw || w == x) {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if (fu > fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return x;
}

std::vector<ContaminationSegment> build_contamination_segments(
        const std::vector<PileupSite>& sites,
        fastgatk::contamination::SegmenterTelemetry* segmenter_stats = nullptr) {
    std::map<std::string, std::vector<std::size_t>> by_contig;
    // ContaminationSegmenter uses Collectors.groupingBy(), whose default
    // HashMap iteration order is observable in the segmentation sidecar.  A
    // lexical std::map would put "10" before "2" and a first-seen vector
    // would not reproduce Java's bucket order.  Retain insertion positions so
    // the small Java HashMap bucket layout (including head insertion on
    // collisions) can be reproduced below while still using the map for
    // deterministic per-contig lookup.
    std::vector<std::string> contig_order;
    std::unordered_map<std::string, bool> seen_contigs;
    contig_order.reserve(sites.size());
    for (std::size_t index = 0; index < sites.size(); ++index) {
        by_contig[sites[index].contig].push_back(index);
        if (seen_contigs.emplace(sites[index].contig, true).second)
            contig_order.push_back(sites[index].contig);
    }
    std::unordered_map<std::string, std::size_t> first_seen;
    first_seen.reserve(contig_order.size() * 2 + 1);
    for (std::size_t index = 0; index < contig_order.size(); ++index)
        first_seen.emplace(contig_order[index], index);
    std::size_t hash_capacity = 16;
    while (contig_order.size() > (hash_capacity * 3) / 4)
        hash_capacity *= 2;
    const auto java_hash_bucket = [hash_capacity](const std::string& key) {
        std::uint32_t hash = 0;
        for (const unsigned char character : key)
            hash = hash * 31U + static_cast<std::uint32_t>(character);
        hash ^= hash >> 16U;
        return static_cast<std::size_t>(hash & static_cast<std::uint32_t>(hash_capacity - 1));
    };
    std::sort(contig_order.begin(), contig_order.end(), [&](const auto& left, const auto& right) {
        const auto left_bucket = java_hash_bucket(left);
        const auto right_bucket = java_hash_bucket(right);
        if (left_bucket != right_bucket) return left_bucket < right_bucket;
        // HashMap.putVal inserts a colliding node at the bucket head.
        return first_seen.at(left) > first_seen.at(right);
    });
    std::vector<ContaminationSegment> segments;
    for (const auto& contig : contig_order) {
        auto& indices = by_contig[contig];
        std::sort(indices.begin(), indices.end(), [&](std::size_t left, std::size_t right) {
            if (sites[left].position != sites[right].position)
                return sites[left].position < sites[right].position;
            return left < right;
        });
        std::vector<std::size_t> hets;
        for (const auto index : indices) {
            const double fraction = sites[index].alt_fraction();
            if (fraction >= 0.1 && fraction <= 0.9) hets.push_back(index);
        }
        if (hets.empty()) continue;
        std::vector<double> het_fractions;
        het_fractions.reserve(hets.size());
        for (const auto index : hets) het_fractions.push_back(sites[index].alt_fraction());
        const auto changepoints = fastgatk::contamination::find_changepoints(
            het_fractions, 10, 100, 50, 1.0, 1.0, segmenter_stats);
        std::vector<std::size_t> starts{0};
        for (const auto changepoint : changepoints)
            starts.push_back(changepoint + 1);
        std::vector<std::size_t> ends = changepoints;
        ends.push_back(hets.size() - 1);
        for (std::size_t block = 0; block < starts.size(); ++block) {
            const auto first_het = hets[starts[block]];
            const auto last_het = hets[ends[block]];
            ContaminationSegment segment;
            for (const auto index : indices)
                if (sites[index].position >= sites[first_het].position &&
                    sites[index].position <= sites[last_het].position)
                    segment.sites.push_back(index);
            if (!segment.sites.empty()) {
                segment.plan_index = segments.size();
                segments.push_back(std::move(segment));
            }
        }
    }
    return segments;
}

double model_error_rate(const std::vector<PileupSite>& sites) {
    std::uint64_t total = 0;
    std::uint64_t other = 0;
    for (const auto& site : sites) {
        total += static_cast<std::uint64_t>(site.total());
        other += static_cast<std::uint64_t>(site.other_count);
    }
    return total == 0 ? std::numeric_limits<double>::quiet_NaN()
                      : 1.5 * static_cast<double>(other) / static_cast<double>(total);
}

LearnedContaminationModel learn_contamination_model(const std::vector<PileupSite>& sites,
                                                    KernelStats* kernel_stats = nullptr,
                                                    fastgatk::contamination::SegmenterTelemetry* segmenter_stats = nullptr) {
    LearnedContaminationModel model;
    model.segments = build_contamination_segments(sites, segmenter_stats);
    model.error_rate = model_error_rate(sites);
    if (model.segments.empty()) {
        if (!sites.empty()) {
            ContaminationSegment all;
            all.sites.resize(sites.size());
            std::iota(all.sites.begin(), all.sites.end(), 0);
            all.plan_index = 0;
            model.segments.push_back(std::move(all));
        }
        return model;
    }
    ContaminationLikelihoodPlan likelihood_plan(sites, model.segments, kernel_stats);
    for (int iteration = 0; iteration < 3; ++iteration) {
        for (auto& segment : model.segments) {
            // GATK uses OptimizationUtils.max(objective, 0.1, 0.5, 0.4,
            // 0.01, 0.01, 20) for each segment.
            segment.minor_allele_fraction = brent_maximize(
                [&](double maf) {
                    return likelihood_plan.segment_log_likelihood(
                        segment.plan_index, model.contamination, model.error_rate, maf);
                }, 0.1, 0.5, 0.4, 0.01, 0.01, 20);
        }
        std::vector<ContaminationSegment> selected;
        std::size_t total_sites = 0;
        for (const auto& segment : model.segments) total_sites += segment.sites.size();
        std::size_t selected_sites = 0;
        for (const auto& segment : model.segments)
            if (segment.minor_allele_fraction > 0.0) selected_sites += segment.sites.size();
        for (double threshold = 0.40; threshold > 0.0; threshold -= 0.04) {
            selected.clear();
            selected_sites = 0;
            for (const auto& segment : model.segments)
                if (segment.minor_allele_fraction > threshold) {
                    selected.push_back(segment);
                    selected_sites += segment.sites.size();
                }
            if (total_sites > 0 && static_cast<double>(selected_sites) /
                    static_cast<double>(total_sites) > 0.25) break;
        }
        if (selected.empty()) selected = model.segments;
        // ContaminationModel evaluates the same Brent objective from four
        // fixed initial guesses and keeps the optimum with the largest
        // objective value.  Do not collapse this to one global grid search:
        // the likelihood can be locally multimodal for small panels.
        constexpr std::array<double, 4> initial_guesses{0.02, 0.05, 0.1, 0.2};
        double best_contamination = 0.0;
        double best_value = -std::numeric_limits<double>::infinity();
        for (const double initial : initial_guesses) {
            const double candidate = brent_maximize(
                [&](double contamination) {
                    return model_log_likelihood(likelihood_plan, selected,
                                                contamination, model.error_rate);
                }, 0.0, 0.5, initial, 1.0e-4, 1.0e-4, 30);
            const double value = model_log_likelihood(
                likelihood_plan, selected, candidate, model.error_rate);
            if (value > best_value) {
                best_value = value;
                best_contamination = candidate;
            }
        }
        model.contamination = best_contamination;
    }
    return model;
}

bool model_genotype_probability(const PileupSite& site, double contamination,
                                double error_rate, double maf, int genotype) {
    const auto likelihoods = genotype_likelihoods_host(site, contamination, error_rate, maf);
    const double total = plain_sum(likelihoods);
    return total > 0.0 && likelihoods[static_cast<std::size_t>(genotype)] / total > 0.5;
}

double percentile90(std::vector<double> values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double rank = 0.90 * static_cast<double>(values.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(rank));
    const auto upper = static_cast<std::size_t>(std::ceil(rank));
    return values[lower] + (rank - static_cast<double>(lower)) *
           (values[upper] - values[lower]);
}

FullEstimate estimate_contamination_full(const std::vector<PileupSite>& tumor,
                                          const std::optional<std::vector<PileupSite>>& matched,
                                          KernelStats* kernel_stats = nullptr,
                                          fastgatk::contamination::SegmenterTelemetry* segmenter_stats = nullptr) {
    const auto& genotyping_sites = matched ? *matched : tumor;
    const auto model = learn_contamination_model(genotyping_sites, kernel_stats, segmenter_stats);
    FullEstimate result;
    result.segments = model.segments;
    std::vector<const PileupSite*> homs;
    auto add_genotype_sites = [&](int genotype, double min_maf) {
        std::vector<const PileupSite*> selected;
        for (const auto& segment : model.segments) {
            if (!(segment.minor_allele_fraction > min_maf)) continue;
            for (const auto index : segment.sites) {
                const auto& site = genotyping_sites[index];
                if (model_genotype_probability(site, model.contamination,
                                               model.error_rate,
                                               segment.minor_allele_fraction,
                                               genotype))
                    selected.push_back(&site);
            }
        }
        for (const auto& site : tumor) {
            if (std::any_of(selected.begin(), selected.end(), [&](const auto* candidate) {
                    return candidate->contig == site.contig && candidate->position == site.position;
                })) homs.push_back(&site);
        }
    };
    auto calculate_strategy = [&](int strategy, double min_maf) -> std::pair<double, double> {
        homs.clear();
        if (strategy == 3 || strategy == 0) add_genotype_sites(strategy, min_maf);
        else {
            std::vector<const PileupSite*> candidates;
            for (const auto& site : tumor)
                if (site.alt_fraction() < 0.15) candidates.push_back(&site);
            std::vector<double> fractions;
            for (const auto* site : candidates) fractions.push_back(site->alt_fraction());
            const double threshold = std::max(0.1, percentile90(std::move(fractions)));
            for (const auto* site : candidates)
                if (site->alt_fraction() <= threshold) homs.push_back(site);
        }
        const bool hom_alt = strategy == 3;
        const double tumor_error = model_error_rate(tumor);
        double total_depth = 0.0;
        double opposite_depth = 0.0;
        JavaDoubleSum weighted;
        JavaDoubleSum coefficient_one;
        JavaDoubleSum coefficient_two;
        for (const auto* site : homs) {
            const double depth = static_cast<double>(site->total());
            const double frequency = hom_alt ? (1.0 - site->allele_frequency)
                                             : site->allele_frequency;
            const double opposite = static_cast<double>(hom_alt ? site->ref_count : site->alt_count);
            total_depth += depth;
            opposite_depth += opposite;
            weighted.add(depth * frequency);
            coefficient_one.add(frequency * depth);
            coefficient_two.add(frequency * (1.0 - frequency) * depth * depth);
        }
        const double weighted_value = weighted.value();
        const double coefficient_one_value = coefficient_one.value();
        const double coefficient_two_value = coefficient_two.value();
        if (!(weighted_value > 0.0) || homs.empty()) return {std::numeric_limits<double>::quiet_NaN(), 1.0};
        const double error_depth = java_round(total_depth * tumor_error / 3.0);
        const double contamination_estimate = std::max(0.0, opposite_depth - error_depth) / weighted_value;
        const auto error_function = [&](double contamination) {
            return std::sqrt(std::max(0.0, coefficient_one_value * contamination * (1.0 - contamination) +
                                      coefficient_two_value * contamination * contamination)) / weighted_value;
        };
        // Port MathUtils.binarySearchFindZero literally.  In particular it
        // reevaluates bottom/top/mid in that order and stops at a 1e-6
        // bracket; using a fixed iteration count changes the reported error by
        // several 1e-7 on the GATK contamination fixtures.
        auto java_signum = [](double value) {
            if (std::isnan(value)) return std::numeric_limits<double>::quiet_NaN();
            if (value > 0.0) return 1.0;
            if (value < 0.0) return -1.0;
            return value;  // preserve Java's +/-0.0/zero sign semantics
        };
        auto binary_search_find_zero = [&](auto&& function, double lower, double upper,
                                            double precision) -> std::optional<double> {
            double bottom = lower;
            double top = upper;
            while (top - bottom > precision) {
                const double middle = (bottom + top) / 2.0;
                const double bottom_value = function(bottom);
                const double top_value = function(top);
                const double middle_value = function(middle);
                if (java_signum(bottom_value) == java_signum(top_value))
                    return std::nullopt;
                if (java_signum(bottom_value) == java_signum(middle_value))
                    bottom = middle;
                else
                    top = middle;
            }
            return (bottom + top) / 2.0;
        };
        auto error_bound = [&](bool upper) {
            const double lower = upper ? contamination_estimate : 0.0;
            const double high = upper ? 1.0 : contamination_estimate;
            const auto root = binary_search_find_zero(
                [&](double value) {
                    return upper ? value - error_function(value) - contamination_estimate
                                  : value + error_function(value) - contamination_estimate;
                }, lower, high, 1.0e-6);
            return root.has_value()
                ? std::max(std::abs(*root - contamination_estimate), 0.0)
                : error_function(contamination_estimate);
        };
        return {std::min(contamination_estimate, 1.0),
                std::max(error_bound(true), error_bound(false))};
    };
    for (double min_maf = 0.40; min_maf >= 0.0; min_maf -= 0.04) {
        const int strategy = min_maf > 0.25 ? 3 : (min_maf > 0.20 ? 0 : 1);
        const auto estimate = calculate_strategy(strategy, min_maf);
        if (std::isfinite(estimate.first) && estimate.second < estimate.first * 0.2 + 0.001) {
            result.contamination = estimate.first;
            result.error = estimate.second;
            break;
        }
        result.contamination = estimate.first;
        result.error = estimate.second;
    }
    if (!std::isfinite(result.contamination)) {
        result.contamination = 0.0;
        result.error = 1.0;
    }
    result.candidate_sites = homs.size();
    for (const auto* site : homs) {
        result.total_depth += static_cast<std::uint64_t>(site->total());
        result.opposite_depth += static_cast<std::uint64_t>(site->alt_count);
    }
    return result;
}

[[maybe_unused]] Estimate estimate_contamination(const std::vector<PileupSite>& tumor,
                                const std::optional<std::vector<PileupSite>>& matched,
                                const Options& options) {
    (void)options;
    std::unordered_map<std::string, PileupSite> matched_by_key;
    if (matched) {
        matched_by_key.reserve(matched->size() * 2 + 1);
        for (const auto& site : *matched) matched_by_key.emplace(site_key(site), site);
    }
    std::vector<double> opposite(tumor.size(), 0.0);
    std::vector<double> denominator(tumor.size(), 0.0);
    std::vector<std::uint8_t> candidate(tumor.size(), 0);
    std::vector<double> observed_counts;
    std::vector<double> total_counts;
    std::vector<double> candidate_af;
    std::vector<std::uint8_t> observe_ref;
    Kokkos::View<int*> ref_counts("contam_ref", tumor.size());
    Kokkos::View<int*> alt_counts("contam_alt", tumor.size());
    Kokkos::View<int*> other_counts("contam_other", tumor.size());
    Kokkos::View<double*> allele_frequencies("contam_af", tumor.size());
    auto host_ref = Kokkos::create_mirror_view(ref_counts);
    auto host_alt = Kokkos::create_mirror_view(alt_counts);
    auto host_other = Kokkos::create_mirror_view(other_counts);
    auto host_af = Kokkos::create_mirror_view(allele_frequencies);
    for (std::size_t index = 0; index < tumor.size(); ++index) {
        host_ref(index) = tumor[index].ref_count;
        host_alt(index) = tumor[index].alt_count;
        host_other(index) = tumor[index].other_count;
        host_af(index) = tumor[index].allele_frequency;
    }
    Kokkos::View<double*> device_opposite("contam_opposite", tumor.size());
    Kokkos::View<double*> device_denominator("contam_denominator", tumor.size());
    Kokkos::View<std::uint8_t*> device_candidate("contam_candidate", tumor.size());
    Kokkos::deep_copy(ref_counts, host_ref);
    Kokkos::deep_copy(alt_counts, host_alt);
    Kokkos::deep_copy(other_counts, host_other);
    Kokkos::deep_copy(allele_frequencies, host_af);
    const bool matched_mode = matched.has_value();
    Kokkos::parallel_for("calculate_contamination_site_kernel", Kokkos::RangePolicy<>(0, tumor.size()),
        KOKKOS_LAMBDA(const std::size_t index) {
            const int depth = ref_counts(index) + alt_counts(index) + other_counts(index);
            const double af = allele_frequencies(index);
            bool use = false;
            double observed = 0.0;
            double expected = 0.0;
            if (depth > 10) {
                // In tumor-only mode, low-alt sites are conservative hom-ref
                // candidates.  In matched mode the host marks candidates below.
                if (!matched_mode && (static_cast<double>(alt_counts(index)) / depth) < 0.20) {
                    use = true;
                    observed = alt_counts(index);
                    expected = depth * af;
                }
            }
            device_candidate(index) = use ? 1 : 0;
            device_opposite(index) = observed;
            device_denominator(index) = expected;
        });
    Kokkos::fence();
    auto host_opposite = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_opposite);
    auto host_denominator = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_denominator);
    auto host_candidate = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_candidate);
    double numerator = 0.0;
    double denominator_sum = 0.0;
    double error_bases = 0.0;
    std::uint64_t total_depth = 0;
    std::uint64_t opposite_depth = 0;
    for (std::size_t index = 0; index < tumor.size(); ++index) {
        const auto& tumor_site = tumor[index];
        const auto matched_site = matched_by_key.find(site_key(tumor_site));
        const bool hom_alt = matched && matched_site != matched_by_key.end() &&
                             matched_site->second.alt_fraction() > 0.80;
        if (matched && !hom_alt) continue;
        double observed = 0.0;
        double expected = 0.0;
        if (matched) {
            observed = tumor_site.ref_count;
            expected = tumor_site.total() * (1.0 - tumor_site.allele_frequency);
        } else if (host_candidate(index)) {
            observed = host_opposite(index);
            expected = host_denominator(index);
        } else continue;
        if (expected <= 0.0) continue;
        total_depth += static_cast<std::uint64_t>(tumor_site.total());
        opposite_depth += static_cast<std::uint64_t>(std::max(0.0, observed));
        numerator += observed;
        denominator_sum += expected;
        error_bases += tumor_site.other_count;
        candidate[index] = 1;
        observed_counts.push_back(observed);
        total_counts.push_back(static_cast<double>(tumor_site.total()));
        candidate_af.push_back(std::clamp(tumor_site.allele_frequency, 1.0e-9, 1.0 - 1.0e-9));
        observe_ref.push_back(matched ? 1 : 0);
    }
    Estimate result;
    result.candidate_sites = std::count(candidate.begin(), candidate.end(), static_cast<std::uint8_t>(1));
    result.total_depth = total_depth;
    result.opposite_depth = opposite_depth;
    if (denominator_sum <= 0.0 || result.candidate_sites == 0) return result;
    // Evaluate the binomial contamination likelihood on a fixed grid.  A
    // matched-normal hom-alt site observes REF reads with p=c; a tumor-only
    // low-ALT site observes ALT reads with p=c*population_AF.  The fixed grid
    // makes reductions deterministic while the Kokkos kernel is portable to
    // Serial/OpenMP/CUDA/HIP.
    constexpr int GRID_SIZE = 2001;
    constexpr double GRID_STEP = 1.0 / static_cast<double>(GRID_SIZE - 1);
    constexpr double EPS = 1.0e-12;
    Kokkos::View<double*> device_observed("contam_likelihood_observed", observed_counts.size());
    Kokkos::View<double*> device_total("contam_likelihood_total", total_counts.size());
    Kokkos::View<double*> device_af("contam_likelihood_af", candidate_af.size());
    Kokkos::View<std::uint8_t*> device_ref("contam_likelihood_ref", observe_ref.size());
    auto host_likelihood_observed = Kokkos::create_mirror_view(device_observed);
    auto host_likelihood_total = Kokkos::create_mirror_view(device_total);
    auto host_likelihood_af = Kokkos::create_mirror_view(device_af);
    auto host_likelihood_ref = Kokkos::create_mirror_view(device_ref);
    for (std::size_t index = 0; index < observed_counts.size(); ++index) {
        host_likelihood_observed(index) = observed_counts[index];
        host_likelihood_total(index) = total_counts[index];
        host_likelihood_af(index) = candidate_af[index];
        host_likelihood_ref(index) = observe_ref[index];
    }
    Kokkos::deep_copy(device_observed, host_likelihood_observed);
    Kokkos::deep_copy(device_total, host_likelihood_total);
    Kokkos::deep_copy(device_af, host_likelihood_af);
    Kokkos::deep_copy(device_ref, host_likelihood_ref);
    Kokkos::View<double*> device_likelihood("contam_likelihood_grid", GRID_SIZE);
    const std::size_t candidate_count = observed_counts.size();
    Kokkos::parallel_for("calculate_contamination_likelihood_grid", GRID_SIZE,
        KOKKOS_LAMBDA(const int grid_index) {
            const double contamination = static_cast<double>(grid_index) * GRID_STEP;
            double likelihood = 0.0;
            for (std::size_t site = 0; site < candidate_count; ++site) {
                const double p_raw = device_ref(site) != 0
                    ? contamination
                    : contamination * device_af(site);
                const double p = Kokkos::min(1.0 - EPS, Kokkos::max(EPS, p_raw));
                const double k = device_observed(site);
                const double n = device_total(site);
                likelihood += k * Kokkos::log(p) + (n - k) * Kokkos::log(1.0 - p);
            }
            device_likelihood(grid_index) = likelihood;
        });
    Kokkos::fence();
    auto host_likelihood = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, device_likelihood);
    int best_grid = 0;
    for (int grid_index = 1; grid_index < GRID_SIZE; ++grid_index)
        if (host_likelihood(grid_index) > host_likelihood(best_grid)) best_grid = grid_index;
    result.contamination = static_cast<double>(best_grid) * GRID_STEP;
    double information = 0.0;
    for (std::size_t site = 0; site < observed_counts.size(); ++site) {
        const double p_raw = observe_ref[site] != 0
            ? result.contamination : result.contamination * candidate_af[site];
        const double p = std::clamp(p_raw, EPS, 1.0 - EPS);
        const double derivative = observe_ref[site] != 0 ? 1.0 : candidate_af[site];
        information += total_counts[site] * derivative * derivative / (p * (1.0 - p));
    }
    result.error = std::clamp(std::sqrt(1.0 / std::max(1.0, information)), 0.0, 1.0);
    return result;
}

[[maybe_unused]] void write_segmentation(const std::string& path, const std::string& sample,
                        const std::vector<PileupSite>& sites) {
    std::map<std::string, std::vector<double>> values;
    std::map<std::string, std::pair<int, int>> bounds;
    for (const auto& site : sites) {
        values[site.contig].push_back(site.minor_allele_fraction());
        auto& bound = bounds[site.contig];
        if (bound.first == 0 || site.position < bound.first) bound.first = site.position;
        bound.second = std::max(bound.second, site.position);
    }
    std::ostringstream output;
    output << "#<METADATA>SAMPLE=" << sample << '\n'
           << "contig\tstart\tend\tminor_allele_fraction\n";
    for (auto& entry : values) {
        auto& fractions = entry.second;
        std::sort(fractions.begin(), fractions.end());
        const double median = fractions.size() % 2 == 0
            ? 0.5 * (fractions[fractions.size() / 2 - 1] + fractions[fractions.size() / 2])
            : fractions[fractions.size() / 2];
        const auto bound = bounds[entry.first];
        output << entry.first << '\t' << bound.first << '\t' << bound.second << '\t'
               << format_double(median) << '\n';
    }
    write_text_output(path, output.str());
    if (!file_complete(path)) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: segmentation is incomplete");
}

void write_model_segmentation(const std::string& path, const std::string& sample,
                              const std::vector<ContaminationSegment>& segments,
                              const std::vector<PileupSite>& sites) {
    std::ostringstream output;
    output << "#<METADATA>SAMPLE=" << sample << '\n'
           << "contig\tstart\tend\tminor_allele_fraction\n";
    for (const auto& segment : segments) {
        if (segment.sites.empty()) continue;
        const auto& first = sites[segment.sites.front()];
        const auto& last = sites[segment.sites.back()];
        output << first.contig << '\t' << first.position << '\t' << last.position << '\t'
               << format_double(segment.minor_allele_fraction) << '\n';
    }
    write_text_output(path, output.str());
    if (!file_complete(path))
        throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: segmentation is incomplete");
}

int run_tool(const Options& options, const fastgatk::runtime::ResourceSnapshot& resources) {
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
    Kokkos::initialize(settings);
    bool initialized = true;
    try {
        const auto tumor_raw = read_table(options.input);
        const auto tumor_filter = coverage_filter(tumor_raw.sites, options);
        const auto& tumor = tumor_filter.sites;
        std::optional<std::vector<PileupSite>> matched;
        std::string matched_sample;
        CoverageFilterResult matched_filter;
        if (!options.matched_normal.empty()) {
            const auto matched_raw = read_table(options.matched_normal);
            matched_sample = matched_raw.sample;
            matched_filter = coverage_filter(matched_raw.sites, options);
            matched = matched_filter.sites;
        }
        // GATK deliberately treats an empty post-coverage panel as a valid
        // estimate boundary.  ContaminationModel then falls through its
        // hom-site strategies and emits the conservative 0.0 +/- 1.0 record;
        // it also writes a header-only MAF sidecar when requested.  Rejecting
        // the input here made native workflows fail exactly where the Java
        // tool produces a usable (if uninformative) downstream table.
        KernelStats kernel_stats;
        fastgatk::contamination::SegmenterTelemetry segmenter_stats;
        const auto estimate = estimate_contamination_full(
            tumor, matched, &kernel_stats, &segmenter_stats);
        std::ostringstream output;
        output << "sample\tcontamination\terror\n"
               << tumor_raw.sample << '\t' << format_double(estimate.contamination) << '\t'
               << format_double(estimate.error) << '\n';
        write_text_output(options.output, output.str());
        if (!file_complete(options.output))
            throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: contamination table is incomplete");
        if (!options.segmentation.empty()) {
            // GATK learns the genotyping model from the matched normal, but
            // writes the optional segmentation from the tumor model.  Keep
            // those two responsibilities separate instead of leaking normal
            // coordinates/MAFs into the tumor sidecar.
            if (matched) {
                const auto tumor_model = learn_contamination_model(
                    tumor, &kernel_stats, &segmenter_stats);
                write_model_segmentation(options.segmentation, tumor_raw.sample,
                                         tumor_model.segments, tumor);
            } else {
                write_model_segmentation(options.segmentation, tumor_raw.sample,
                                         estimate.segments, tumor);
            }
        }
        const auto manifest_path = options.manifest.empty() ? options.output + ".manifest.json" : options.manifest;
        std::ofstream manifest(manifest_path);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create manifest: " + manifest_path);
        manifest << "{\"schema_version\":1,\"tool\":\"CalculateContamination\","
                 << "\"implementation\":\"fastgatk-calculate-contamination\",\"status\":\"prototype\","
                 << "\"primary_output\":\"" << json_escape(options.output) << "\",\"primary_output_kind\":\"table\","
                 << "\"compatibility\":{\"gatk_table_input\":true,\"coverage_filter\":true,"
                 << "\"empty_post_coverage_result\":true,"
                 << "\"matched_normal_hom_site_estimate\":" << (matched ? "true" : "false")
                 << ",\"tumor_only_estimate\":" << (matched ? "false" : "true")
                 << ",\"kokkos_site_kernel\":true,\"kernel_segmenter\":true,"
                 << "\"segment_coordinates_gatk_oracle\":true,\"bit_identical_to_gatk\":false,"
                 << "\"full_contamination_model\":true,\"commons_math_brent_optimizer\":true,"
                 << "\"commons_math_binary_search_error\":true,\"strict_numeric_oracle\":true,"
                 << "\"compressed_io\":"
                 << (compressed_path(options.input) || compressed_path(options.output) ? "true" : "false")
                 << "},\"outputs\":[{\"path\":\""
                 << json_escape(options.output) << "\",\"kind\":\"table\",\"complete\":true}],"
                 << "\"telemetry\":{\"resources\":" << resources.to_json()
                 << ",\"sample\":\"" << json_escape(tumor_raw.sample) << "\",\"matched_sample\":\""
                 << json_escape(matched_sample) << "\",\"raw_sites\":" << tumor_raw.sites.size()
                 << ",\"covered_sites\":" << tumor_filter.covered_sites
                 << ",\"filtered_sites\":" << tumor.size()
                 << ",\"empty_filtered_tumor\":" << (tumor.empty() ? "true" : "false")
                 << ",\"empty_filtered_matched\":"
                 << (matched && matched->empty() ? "true" : "false")
                 << ",\"coverage_median\":" << json_number(tumor_filter.median)
                 << ",\"coverage_mean\":" << json_number(tumor_filter.mean)
                 << ",\"coverage_low_threshold\":" << json_number(tumor_filter.low)
                 << ",\"coverage_high_threshold\":" << json_number(tumor_filter.high)
                 << ",\"candidate_sites\":" << estimate.candidate_sites
                 << ",\"total_depth\":" << estimate.total_depth << ",\"opposite_depth\":" << estimate.opposite_depth
                 << ",\"segmentation\":" << (!options.segmentation.empty() ? "true" : "false")
                 << ",\"input_compressed\":" << (compressed_path(options.input) ? "true" : "false")
                 << ",\"output_compressed\":" << (compressed_path(options.output) ? "true" : "false")
                 << ",\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"optimizer\":\"CommonsMathBrentOptimizer\",\"error_search\":\"MathUtils.binarySearchFindZero\""
                 << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
                 << ",\"kernel_execution_policy\":\"RangePolicy\""
                 << ",\"kernel_batches\":" << kernel_stats.batches
                 << ",\"kernel_observations\":" << kernel_stats.observations
                 << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                 << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds
                 << ",\"likelihood_plan_allocations\":" << kernel_stats.plan_allocations
                 << ",\"likelihood_plan_reuses\":" << kernel_stats.plan_reuses
                 << ",\"segmenter_kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\""
                 << ",\"segmenter_kernel_execution_space\":\""
                 << json_escape(segmenter_stats.execution_space)
                 << "\",\"segmenter_kernel_execution_policy\":\""
                 << segmenter_stats.execution_policy
                 << "\",\"segmenter_kernel_batches\":" << segmenter_stats.batches
                 << ",\"segmenter_kernel_observations\":" << segmenter_stats.observations
                 << ",\"segmenter_kernel_prepare_seconds\":" << segmenter_stats.prepare_seconds
                 << ",\"segmenter_kernel_execute_seconds\":" << segmenter_stats.execute_seconds << "}}\n";
        std::cout << "{\"tool\":\"CalculateContamination\",\"status\":\"prototype\","
                  << "\"sample\":\"" << json_escape(tumor_raw.sample) << "\",\"contamination\":"
                  << estimate.contamination << ",\"error\":" << estimate.error
                  << ",\"optimizer\":\"CommonsMathBrentOptimizer\",\"strict_numeric_oracle\":true"
                  << ",\"kernel_batches\":" << kernel_stats.batches
                  << ",\"kernel_observations\":" << kernel_stats.observations
                  << ",\"likelihood_plan_allocations\":" << kernel_stats.plan_allocations
                  << ",\"likelihood_plan_reuses\":" << kernel_stats.plan_reuses
                  << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                  << "\",\"segmenter_kernel_batches\":" << segmenter_stats.batches
                  << ",\"segmenter_kernel_observations\":" << segmenter_stats.observations
                  << "}\n";
        Kokkos::finalize();
        initialized = false;
        return 0;
    } catch (...) {
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        throw;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_tool(parse(argc, argv), fastgatk::runtime::ResourceSnapshot::probe());
    } catch (const std::exception& error) {
        std::cerr << "fastgatk-calculate-contamination: " << error.what() << '\n';
        return 2;
    }
}
