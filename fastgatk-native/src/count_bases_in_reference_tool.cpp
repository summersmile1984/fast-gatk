#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include "fastgatk/io/hts_read_guard.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct Region {
    std::string contig;
    hts_pos_t start = 1;
    hts_pos_t end = 0;
};

struct Options {
    std::string reference;
    std::vector<std::string> selectors;
    std::string output;
    std::string manifest;
    int threads = 1;
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string value_for(int& index, int argc, char** argv, const std::string& argument,
                     std::string_view name, const char* short_name = nullptr) {
    if (const auto value = inline_value(argument, name); !value.empty()) return value;
    if ((argument == name || (short_name && argument == short_name)) && index + 1 < argc)
        return argv[++index];
    throw std::invalid_argument("missing value for " + std::string(name));
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-count-bases-in-reference (GATK-compatible native)\n"
                         "  -R, --reference FILE       reference FASTA (+ .fai)\n"
                         "  -L, --intervals REGION     optional contig or contig:start-end (repeatable)\n"
                         "  -O, --output FILE          optional copy of the base-count report\n"
                         "      --threads N            Kokkos execution threads\n"
                         "      --output-manifest FILE OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || argument == "--reference" ||
                   !inline_value(argument, "--reference").empty()) {
            options.reference = value_for(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || argument == "--intervals" || argument == "--interval" ||
                   argument == "--region" || !inline_value(argument, "--intervals").empty() ||
                   !inline_value(argument, "--interval").empty() ||
                   !inline_value(argument, "--region").empty()) {
            const std::string_view name = argument == "-L" || argument == "--intervals" ? "--intervals" :
                (argument == "--interval" || !inline_value(argument, "--interval").empty() ? "--interval" :
                 (argument == "--region" || !inline_value(argument, "--region").empty() ? "--region" : "--intervals"));
            options.selectors.push_back(value_for(index, argc, argv, argument, name, "-L"));
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(value_for(index, argc, argv, argument, "--threads"));
        } else if (argument == "-O" || argument == "--output" ||
                   !inline_value(argument, "--output").empty()) {
            options.output = value_for(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() ||
                   !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = value_for(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--seconds-between-progress-updates", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                     "--seconds-between-progress-updates");
                (void)value_for(index, argc, argv, argument, name);
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    return options;
}

Region parse_region(const std::string& selector, const faidx_t* fai) {
    const auto colon = selector.find(':');
    const std::string contig = colon == std::string::npos ? selector : selector.substr(0, colon);
    if (contig.empty()) throw std::invalid_argument("empty interval contig");
    const auto length = faidx_seq_len64(fai, contig.c_str());
    if (length < 0) throw std::invalid_argument("interval contig is absent from reference: " + contig);
    if (colon == std::string::npos) return {contig, 1, length};
    const auto dash = selector.find('-', colon + 1);
    const auto start_text = selector.substr(colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1);
    const auto end_text = dash == std::string::npos ? start_text : selector.substr(dash + 1);
    if (start_text.empty() || end_text.empty()) throw std::invalid_argument("malformed interval: " + selector);
    const auto start = static_cast<hts_pos_t>(std::stoll(start_text));
    const auto end = static_cast<hts_pos_t>(std::stoll(end_text));
    if (start < 1 || end < start || end > length)
        throw std::invalid_argument("interval is outside reference: " + selector);
    return {contig, start, end};
}

std::string trim_copy(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

bool suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    const auto offset = value.size() - ending.size();
    for (std::size_t index = 0; index < ending.size(); ++index) {
        if (std::tolower(static_cast<unsigned char>(value[offset + index])) !=
            std::tolower(static_cast<unsigned char>(ending[index]))) return false;
    }
    return true;
}

void append_interval_file(const std::string& path, const faidx_t* fai,
                          std::vector<Region>& output) {
    const bool bed = suffix_ci(path, ".bed") || suffix_ci(path, ".bed.gz");
    std::size_t parsed = 0;
    const auto process_line = [&](std::string line) {
        line = trim_copy(std::move(line));
        if (line.empty() || line.front() == '#' || line.front() == '@' ||
            line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) return;
        std::istringstream fields(line);
        std::string first;
        fields >> first;
        if (first.find(':') != std::string::npos) {
            output.push_back(parse_region(first, fai));
            ++parsed;
            return;
        }
        std::string second, third;
        if (!(fields >> second >> third))
            throw std::invalid_argument("malformed interval file line: " + line);
        const auto length = faidx_seq_len64(fai, first.c_str());
        if (length < 0) throw std::invalid_argument("interval contig is absent from reference: " + first);
        std::size_t consumed_start = 0, consumed_end = 0;
        const auto start = std::stoll(second, &consumed_start);
        const auto end = std::stoll(third, &consumed_end);
        if (consumed_start != second.size() || consumed_end != third.size())
            throw std::invalid_argument("malformed interval coordinates: " + line);
        // BED is zero-based half-open; Picard/GATK interval-list rows are
        // one-based inclusive.  Keep Region in the latter representation,
        // matching parse_region and the faidx_fetch_seq64 conversion below.
        const auto begin = bed ? start + 1 : start;
        const auto finish = end;
        if (begin < 1 || finish < begin || finish > length)
            throw std::invalid_argument("interval is outside reference: " + line);
        output.push_back({first, static_cast<hts_pos_t>(begin),
                          static_cast<hts_pos_t>(finish)});
        ++parsed;
    };
    if (suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::invalid_argument("cannot open compressed interval file: " + path);
        kstring_t line{0, 0, nullptr};
        try {
            while (fastgatk::io::read_text_line(file, &line, path) >= 0)
                process_line(line.s == nullptr ? std::string{} : std::string(line.s, line.l));
        } catch (...) {
            free(line.s);
            hts_close(file);
            throw;
        }
        const auto close_status = hts_close(file);
        free(line.s);
        if (close_status != 0)
            throw std::invalid_argument("failed reading interval file: " + path);
    } else {
        std::ifstream stream(path);
        if (!stream) throw std::invalid_argument("cannot open interval file: " + path);
        for (std::string line; std::getline(stream, line);)
            process_line(std::move(line));
        if (stream.bad()) throw std::invalid_argument("failed reading interval file: " + path);
    }
    if (parsed == 0) throw std::invalid_argument("interval file contains no intervals: " + path);
}

std::vector<Region> selected_regions(const Options& options, const faidx_t* fai) {
    if (options.selectors.empty()) {
        std::vector<Region> result;
        result.reserve(static_cast<std::size_t>(faidx_nseq(fai)));
        for (int index = 0; index < faidx_nseq(fai); ++index) {
            const char* contig = faidx_iseq(fai, index);
            result.push_back({contig, 1, faidx_seq_len64(fai, contig)});
        }
        return result;
    }
    std::vector<Region> result;
    result.reserve(options.selectors.size());
    for (const auto& selector : options.selectors) {
        std::error_code error;
        if (std::filesystem::is_regular_file(selector, error) && !error)
            append_interval_file(selector, fai, result);
        else
            result.push_back(parse_region(selector, fai));
    }
    // GATK's interval traversal is set-like: overlapping selectors must not
    // count the same reference base more than once.  Merge after parsing so
    // literal, BED, and interval-list selectors share one deterministic path.
    std::sort(result.begin(), result.end(), [](const Region& left, const Region& right) {
        if (left.contig != right.contig) return left.contig < right.contig;
        if (left.start != right.start) return left.start < right.start;
        return left.end < right.end;
    });
    std::vector<Region> merged;
    merged.reserve(result.size());
    for (const auto& region : result) {
        if (!merged.empty() && merged.back().contig == region.contig &&
            region.start <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, region.end);
        } else {
            merged.push_back(region);
        }
    }
    result.swap(merged);
    return result;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char character : value) {
        if (character == '"' || character == '\\') out << '\\';
        if (character == '\n') out << "\\n";
        else if (character == '\r') out << "\\r";
        else if (character == '\t') out << "\\t";
        else out << character;
    }
    return out.str();
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    faidx_t* fai = fai_load(options.reference.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
    try {
        const auto regions = selected_regions(options, fai);
        std::vector<unsigned char> bases;
        for (const auto& region : regions) {
            hts_pos_t fetched = 0;
            char* sequence = faidx_fetch_seq64(fai, region.contig.c_str(), region.start - 1,
                                               region.end - 1, &fetched);
            if (!sequence || fetched != region.end - region.start + 1) {
                free(sequence);
                throw std::runtime_error("BAD_INPUT: failed to fetch reference interval: " + region.contig);
            }
            bases.insert(bases.end(), sequence, sequence + fetched);
            free(sequence);
        }

        Kokkos::InitializationSettings settings;
        settings.set_num_threads(options.threads);
        Kokkos::initialize(settings);
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        fastgatk::core::HostBatch host("count-bases-in-reference-v1");
        host.records = bases.size();
        host.bytes = bases.size() * sizeof(unsigned char);
        fastgatk::core::KernelPlan<ExecSpace> plan("count-bases-in-reference");
        plan.begin_prepare(host);
        std::array<std::uint64_t, 256> counts{};
        {
            Kokkos::View<unsigned char*> device_bases("reference_bases", bases.size());
            Kokkos::View<std::uint64_t*> device_counts("reference_base_counts", 256);
            auto host_bases = Kokkos::create_mirror_view(device_bases);
            for (std::size_t i = 0; i < bases.size(); ++i) host_bases(i) = bases[i];
            Kokkos::deep_copy(device_bases, host_bases);
            Kokkos::deep_copy(device_counts, std::uint64_t{0});
            fastgatk::core::DeviceBatch<ExecSpace> device(bases.size());
            device.bind("bases", device_bases);
            device.bind("counts", device_counts);
            ExecSpace().fence();
            plan.end_prepare(device);
            plan.begin_execute();
            Kokkos::parallel_for("count_bases_in_reference", Kokkos::RangePolicy<ExecSpace>(0, bases.size()),
                KOKKOS_LAMBDA(const std::size_t index) {
                    Kokkos::atomic_fetch_add(&device_counts(static_cast<unsigned char>(device_bases(index))),
                                             std::uint64_t{1});
                });
            ExecSpace().fence();
            plan.end_execute();
            auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_counts);
            for (std::size_t i = 0; i < counts.size(); ++i) counts[i] = host_counts(i);
        }
        Kokkos::finalize();

        std::ostringstream report;
        for (std::size_t i = 0; i < counts.size(); ++i)
            if (counts[i] > 0) report << static_cast<char>(i) << " : " << counts[i] << '\n';
        const auto text = report.str();
        std::cout << text;
        if (!options.output.empty()) {
            std::ofstream output(options.output);
            if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write output: " + options.output);
            output << text;
            if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize output: " + options.output);
        }
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            const auto& telemetry = plan.telemetry();
            manifest << "{\"schema_version\":1,\"tool\":\"CountBasesInReference\",\"implementation\":\"fastgatk-count-bases-in-reference\",\"status\":\"contract-compatible\",\"reference\":\""
                     << json_escape(options.reference) << "\",\"regions\":" << regions.size()
                     << ",\"bases\":" << bases.size() << ",\"output\":\"" << json_escape(options.output)
                     << "\",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\""
                     << ExecSpace::name() << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_batches\":1,\"kernel_records\":"
                     << bases.size() << ",\"kernel_prepare_seconds\":" << telemetry.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << telemetry.execute_seconds << "}}\n";
        }
        std::cerr << "{\"tool\":\"CountBasesInReference\",\"status\":\"contract-compatible\",\"regions\":"
                  << regions.size() << ",\"bases\":" << bases.size() << ",\"kernel_execution_space\":\""
                  << ExecSpace::name() << "\",\"kernel_records\":" << bases.size()
                  << ",\"kernel_execute_seconds\":" << plan.telemetry().execute_seconds << "}\n";
        fai_destroy(fai);
        return 0;
    } catch (...) {
        fai_destroy(fai);
        throw;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return main_impl(argc, argv);
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
