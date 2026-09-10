#include "fastgatk/core/plan.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Interval {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
};

struct Options {
    std::string reference;
    std::vector<std::string> regions;
    std::vector<std::string> exclusions;
    std::string output;
    std::string manifest;
    std::int64_t padding = 250;
    std::int64_t exclusion_padding = 0;
    std::int64_t bin_length = 1000;
    // GATK 4.6.2.0 validates the shared interval collection and requires
    // OVERLAPPING_ONLY for this tool, even though the generic argument
    // default is ALL. Keep omission on the same fail-closed path as ALL.
    std::string interval_merging_rule = "ALL";
    int threads = 1;
};

struct KernelStats {
    std::uint64_t batches = 0;
    std::uint64_t records = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

std::string option_value(const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

bool has_option(const std::string& argument, const char* name) {
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

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-preprocess-intervals (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE             reference FASTA (+ .fai)\n"
                         "  -L, --intervals REGION           interval or interval-list/BED (repeatable; optional for WGS)\n"
                         "  -XL, --exclude-intervals REGION  subtract interval or interval-list/BED (repeatable)\n"
                         "      --padding N                  padding on each side (default 250)\n"
                         "      --interval-exclusion-padding N  pad exclusions before subtraction\n"
                         "      --bin-length N               bin size, 0 disables binning (default 1000)\n"
                         "  -O, --output FILE                Picard interval-list output\n"
                         "      --threads N                  Kokkos execution threads\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || has_option(argument, "--reference")) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-L" || has_option(argument, "--intervals") ||
                   has_option(argument, "--interval") || has_option(argument, "--region")) {
            const char* name = argument == "-L" ? "--intervals" :
                has_option(argument, "--region") ? "--region" :
                has_option(argument, "--interval") ? "--interval" : "--intervals";
            options.regions.push_back(require_value(index, argc, argv, argument, name, "-L"));
        } else if (argument == "-XL" || has_option(argument, "--exclude-intervals")) {
            const char* name = argument == "-XL" ? "--exclude-intervals" : "--exclude-intervals";
            options.exclusions.push_back(require_value(index, argc, argv, argument, name, "-XL"));
        } else if (argument == "-O" || has_option(argument, "--output")) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (has_option(argument, "--padding")) {
            options.padding = std::stoll(require_value(index, argc, argv, argument, "--padding"));
        } else if (has_option(argument, "--bin-length")) {
            options.bin_length = std::stoll(require_value(index, argc, argv, argument, "--bin-length"));
        } else if (has_option(argument, "--threads")) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (has_option(argument, "--output-manifest") || has_option(argument, "--manifest")) {
            options.manifest = require_value(index, argc, argv, argument,
                argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest");
        } else if (has_option(argument, "--interval-merging-rule")) {
            options.interval_merging_rule = require_value(index, argc, argv, argument, "--interval-merging-rule");
        } else if (has_option(argument, "--interval-set-rule")) {
            const auto rule = require_value(index, argc, argv, argument, "--interval-set-rule");
            if (rule != "UNION")
                throw std::invalid_argument("native PreprocessIntervals requires --interval-set-rule UNION");
        } else if (has_option(argument, "--interval-padding") || has_option(argument, "--interval-exclusion-padding")) {
            const bool exclusion = argument.rfind("--interval-exclusion-padding", 0) == 0;
            const char* name = exclusion ? "--interval-exclusion-padding" : "--interval-padding";
            const auto value = std::stoll(require_value(index, argc, argv, argument, name));
            if (value < 0) throw std::invalid_argument(std::string(name) + " must be non-negative");
            if (exclusion) {
                // CopyNumberArgumentValidationUtils deliberately rejects a
                // non-zero exclusion padding for PreprocessIntervals.  Keep
                // the option visible for command-line compatibility while
                // matching that fail-closed behavior exactly.
                if (value != 0)
                    throw std::invalid_argument("native PreprocessIntervals requires --interval-exclusion-padding 0");
                options.exclusion_padding = value;
            } else if (value != 0) {
                throw std::invalid_argument("native PreprocessIntervals requires --interval-padding 0; use --padding for target padding");
            }
        } else if (argument == "--quiet" || has_option(argument, "--verbosity") ||
                   has_option(argument, "--seconds-between-progress-updates") ||
                   has_option(argument, "--java-options")) {
            if (argument != "--quiet") {
                (void)require_value(index, argc, argv, argument,
                    argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    argument.rfind("--verbosity", 0) == 0 ? "--verbosity" :
                    "--seconds-between-progress-updates");
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.padding < 0) throw std::invalid_argument("--padding must be non-negative");
    if (options.exclusion_padding < 0) throw std::invalid_argument("--interval-exclusion-padding must be non-negative");
    if (options.bin_length < 0) throw std::invalid_argument("--bin-length must be non-negative");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    if (options.interval_merging_rule != "OVERLAPPING_ONLY")
        throw std::invalid_argument("Interval merging rule must be set to OVERLAPPING_ONLY.");
    return options;
}

std::vector<std::string> split_fields(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) fields.push_back(std::move(field));
    return fields;
}

bool suffix_ci(const std::string& value, const char* suffix) {
    const std::string ending(suffix);
    if (value.size() < ending.size()) return false;
    for (std::size_t i = 0; i < ending.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(value[value.size() - ending.size() + i])) !=
            std::tolower(static_cast<unsigned char>(ending[i]))) return false;
    }
    return true;
}

std::vector<std::string> read_lines(const std::string& path) {
    std::vector<std::string> lines;
    if (suffix_ci(path, ".gz")) {
        htsFile* file = hts_open(path.c_str(), "r");
        if (!file) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        kstring_t line{0, 0, nullptr};
        while (hts_getline(file, '\n', &line) >= 0)
            lines.emplace_back(line.s == nullptr ? "" : std::string(line.s, line.l));
        free(line.s);
        if (hts_close(file) != 0) throw std::runtime_error("BAD_INPUT: failed reading interval file: " + path);
    } else {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("BAD_INPUT: cannot open interval file: " + path);
        for (std::string line; std::getline(input, line);) lines.push_back(std::move(line));
    }
    return lines;
}

Interval parse_token(const std::string& token, const faidx_t* fai) {
    const auto colon = token.find(':');
    const std::string contig = colon == std::string::npos ? token : token.substr(0, colon);
    const auto raw_length = faidx_seq_len64(fai, contig.c_str());
    if (contig.empty() || raw_length < 0)
        throw std::invalid_argument("BAD_INPUT: interval contig is absent from reference: " + contig);
    const std::int64_t length = static_cast<std::int64_t>(raw_length);
    if (colon == std::string::npos) return {contig, 1, length};
    const auto dash = token.find('-', colon + 1);
    const auto start = std::stoll(token.substr(colon + 1, dash == std::string::npos ? std::string::npos : dash - colon - 1));
    const auto end = dash == std::string::npos ? start : std::stoll(token.substr(dash + 1));
    if (start < 1 || end < start || end > length)
        throw std::invalid_argument("BAD_INPUT: interval is outside reference: " + token);
    return {contig, start, end};
}

std::vector<Interval> parse_selector(const std::string& selector, const faidx_t* fai) {
    std::error_code error;
    const bool file = std::filesystem::is_regular_file(selector, error) && !error;
    if (!file) return {parse_token(selector, fai)};
    const bool bed = suffix_ci(selector, ".bed") || suffix_ci(selector, ".bed.gz");
    std::vector<Interval> result;
    for (std::size_t line_number = 0; const auto& line : read_lines(selector)) {
        ++line_number;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#' || line[first] == '@') continue;
        const auto fields = split_fields(line.substr(first));
        if (fields.empty() || fields[0] == "track" || fields[0] == "browser") continue;
        try {
            if (bed) {
                if (fields.size() < 3) throw std::invalid_argument("BED requires contig/start/end");
                const auto start0 = std::stoll(fields[1]);
                const auto end0 = std::stoll(fields[2]);
                if (start0 < 0 || end0 <= start0) throw std::invalid_argument("invalid BED coordinates");
                result.push_back(parse_token(fields[0] + ":" + std::to_string(start0 + 1) + "-" + std::to_string(end0), fai));
            } else if (fields.size() == 1) {
                result.push_back(parse_token(fields[0], fai));
            } else if (fields.size() >= 3) {
                result.push_back(parse_token(fields[0] + ":" + fields[1] + "-" + fields[2], fai));
            } else {
                throw std::invalid_argument("interval list row requires contig/start/end");
            }
        } catch (const std::exception& exception) {
            throw std::invalid_argument("BAD_INPUT: malformed interval line " + std::to_string(line_number) +
                                        " in " + selector + ": " + exception.what());
        }
    }
    if (result.empty()) throw std::invalid_argument("BAD_INPUT: interval file is empty: " + selector);
    return result;
}

std::vector<Interval> normalize_intervals(const std::vector<Interval>& input, const faidx_t* fai) {
    std::unordered_map<std::string, int> order;
    for (int i = 0; i < faidx_nseq(fai); ++i) order.emplace(faidx_iseq(fai, i), i);
    auto result = input;
    std::sort(result.begin(), result.end(), [&](const Interval& lhs, const Interval& rhs) {
        const auto li = order.at(lhs.contig), ri = order.at(rhs.contig);
        return li != ri ? li < ri : (lhs.start != rhs.start ? lhs.start < rhs.start : lhs.end < rhs.end);
    });
    std::vector<Interval> merged;
    for (const auto& interval : result) {
        if (!merged.empty() && merged.back().contig == interval.contig && interval.start <= merged.back().end) {
            merged.back().end = std::max(merged.back().end, interval.end);
        } else {
            merged.push_back(interval);
        }
    }
    return merged;
}

std::vector<Interval> pad_intervals(const std::vector<Interval>& input, std::int64_t padding, const faidx_t* fai) {
    std::vector<Interval> padded;
    padded.reserve(input.size());
    for (const auto& interval : input) {
        const auto length = static_cast<std::int64_t>(faidx_seq_len64(fai, interval.contig.c_str()));
        padded.push_back({interval.contig, std::max<std::int64_t>(1, interval.start - padding),
                          std::min<std::int64_t>(length, interval.end + padding)});
    }
    for (std::size_t i = 0; i + 1 < padded.size(); ++i) {
        if (padded[i].contig != padded[i + 1].contig || padded[i].end < padded[i + 1].start) continue;
        const auto midpoint = (input[i].end + input[i + 1].start) / 2;
        padded[i].end = midpoint;
        padded[i + 1].start = midpoint + 1;
    }
    return padded;
}

std::vector<Interval> subtract_intervals(const std::vector<Interval>& input,
                                         const std::vector<Interval>& exclusions) {
    if (exclusions.empty()) return input;
    std::map<std::string, std::vector<Interval>> by_contig;
    for (const auto& exclusion : exclusions) by_contig[exclusion.contig].push_back(exclusion);
    for (auto& [contig, values] : by_contig) {
        std::sort(values.begin(), values.end(), [](const Interval& lhs, const Interval& rhs) {
            return lhs.start != rhs.start ? lhs.start < rhs.start : lhs.end < rhs.end;
        });
        std::vector<Interval> merged;
        for (const auto& value : values) {
            // OVERLAPPING_ONLY is the default GATK interval merge rule.  Keep
            // adjacent bases distinct (start == end + 1) while coalescing
            // genuine overlap, so subtraction never consumes an unrequested
            // boundary base.
            if (!merged.empty() && value.start <= merged.back().end)
                merged.back().end = std::max(merged.back().end, value.end);
            else
                merged.push_back(value);
        }
        values = std::move(merged);
    }

    std::vector<Interval> result;
    for (const auto& source : input) {
        std::int64_t cursor = source.start;
        const auto found = by_contig.find(source.contig);
        if (found != by_contig.end()) {
            for (const auto& exclusion : found->second) {
                if (exclusion.end < cursor) continue;
                if (exclusion.start > source.end) break;
                if (exclusion.start > cursor)
                    result.push_back({source.contig, cursor,
                                      std::min(source.end, exclusion.start - 1)});
                if (exclusion.end >= source.end) {
                    cursor = source.end + 1;
                    break;
                }
                cursor = std::max(cursor, exclusion.end + 1);
            }
        }
        if (cursor <= source.end) result.push_back({source.contig, cursor, source.end});
    }
    return result;
}

std::vector<Interval> generate_bins(const std::vector<Interval>& padded, std::int64_t bin_length) {
    if (bin_length == 0) return padded;
    std::vector<Interval> bins;
    for (const auto& interval : padded) {
        for (std::int64_t start = interval.start; start <= interval.end; start += bin_length) {
            bins.push_back({interval.contig, start, std::min(interval.end, start + bin_length - 1)});
        }
    }
    return bins;
}

std::vector<Interval> filter_n_only(const std::vector<Interval>& bins, const faidx_t* fai,
                                    int threads, KernelStats& stats) {
    std::vector<std::uint8_t> keep(bins.size(), 0);
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(threads);
    Kokkos::initialize(settings);
    try {
        constexpr std::size_t batch_size = 4096;
        using ExecSpace = Kokkos::DefaultExecutionSpace;
        for (std::size_t batch_start = 0; batch_start < bins.size(); batch_start += batch_size) {
            const auto batch_end = std::min(bins.size(), batch_start + batch_size);
            std::vector<unsigned char> packed;
            std::vector<std::size_t> offsets{0};
            packed.reserve((batch_end - batch_start) * 1024);
            for (std::size_t i = batch_start; i < batch_end; ++i) {
                hts_pos_t length = 0;
                char* sequence = faidx_fetch_seq64(fai, bins[i].contig.c_str(), bins[i].start - 1, bins[i].end - 1, &length);
                if (!sequence || length != bins[i].end - bins[i].start + 1) {
                    free(sequence);
                    throw std::runtime_error("BAD_INPUT: failed to fetch reference bin: " + bins[i].contig);
                }
                packed.insert(packed.end(), sequence, sequence + length);
                offsets.push_back(offsets.back() + static_cast<std::size_t>(length));
                free(sequence);
            }
            fastgatk::core::HostBatch filter_host("preprocess-intervals-filter-n-v1");
            filter_host.records = batch_end - batch_start;
            filter_host.bytes = packed.size() * sizeof(unsigned char) + offsets.size() * sizeof(std::size_t);
            fastgatk::core::KernelPlan<ExecSpace> filter_plan("preprocess-intervals-filter-n");
            filter_plan.begin_prepare(filter_host);
            Kokkos::View<unsigned char*> device_bases("preprocess_interval_bases", packed.size());
            Kokkos::View<std::size_t*> device_offsets("preprocess_interval_offsets", offsets.size());
            Kokkos::View<std::uint8_t*> device_keep("preprocess_interval_keep", batch_end - batch_start);
            auto host_bases = Kokkos::create_mirror_view(device_bases);
            auto host_offsets = Kokkos::create_mirror_view(device_offsets);
            for (std::size_t i = 0; i < packed.size(); ++i) host_bases(i) = packed[i];
            for (std::size_t i = 0; i < offsets.size(); ++i) host_offsets(i) = offsets[i];
            Kokkos::deep_copy(device_bases, host_bases);
            Kokkos::deep_copy(device_offsets, host_offsets);
            fastgatk::core::DeviceBatch<ExecSpace> filter_device(batch_end - batch_start);
            filter_device.bind("bases", device_bases);
            filter_device.bind("offsets", device_offsets);
            filter_device.bind("keep", device_keep);
            ExecSpace().fence();
            filter_plan.end_prepare(filter_device);
            filter_plan.begin_execute();
            Kokkos::parallel_for("preprocess_intervals_filter_n", Kokkos::RangePolicy<ExecSpace>(0, batch_end - batch_start),
                KOKKOS_LAMBDA(const std::size_t index) {
                    bool all_n = true;
                    for (std::size_t position = device_offsets(index); position < device_offsets(index + 1); ++position) {
                        const auto base = device_bases(position);
                        if (base != 'N' && base != 'n') { all_n = false; break; }
                    }
                    device_keep(index) = all_n ? 0 : 1;
                });
            ExecSpace().fence();
            filter_plan.end_execute();
            ++stats.batches;
            stats.records += batch_end - batch_start;
            stats.prepare_seconds += filter_plan.telemetry().prepare_seconds;
            stats.execute_seconds += filter_plan.telemetry().execute_seconds;
            auto host_keep = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), device_keep);
            for (std::size_t i = 0; i < batch_end - batch_start; ++i) keep[batch_start + i] = host_keep(i);
        }
        Kokkos::finalize();
    } catch (...) {
        Kokkos::finalize();
        throw;
    }
    std::vector<Interval> result;
    result.reserve(bins.size());
    for (std::size_t i = 0; i < bins.size(); ++i) if (keep[i]) result.push_back(bins[i]);
    return result;
}

std::string json_escape(const std::string& text) {
    std::ostringstream out;
    for (char c : text) {
        if (c == '"' || c == '\\') out << '\\';
        if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else out << c;
    }
    return out.str();
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    faidx_t* fai = fai_load(options.reference.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
    try {
        std::vector<Interval> requested;
        if (options.regions.empty()) {
            for (int i = 0; i < faidx_nseq(fai); ++i) {
                const auto contig = std::string(faidx_iseq(fai, i));
                requested.push_back({contig, 1, static_cast<std::int64_t>(faidx_seq_len64(fai, contig.c_str()))});
            }
        } else {
            for (const auto& selector : options.regions) {
                const auto intervals = parse_selector(selector, fai);
                requested.insert(requested.end(), intervals.begin(), intervals.end());
            }
        }
        const auto normalized = normalize_intervals(requested, fai);
        if (normalized.empty()) throw std::invalid_argument("BAD_INPUT: no intervals remain after normalization");
        std::vector<Interval> requested_exclusions;
        for (const auto& selector : options.exclusions) {
            const auto parsed = parse_selector(selector, fai);
            requested_exclusions.insert(requested_exclusions.end(), parsed.begin(), parsed.end());
        }
        const auto normalized_exclusions = normalize_intervals(requested_exclusions, fai);
        const auto padded_exclusions = pad_intervals(normalized_exclusions, options.exclusion_padding, fai);
        const auto after_exclusion = subtract_intervals(normalized, padded_exclusions);
        if (after_exclusion.empty())
            throw std::invalid_argument("BAD_INPUT: no intervals remain after applying --exclude-intervals");
        const auto padded = pad_intervals(after_exclusion, options.padding, fai);
        const auto unfiltered_bins = generate_bins(padded, options.bin_length);
        KernelStats kernel_stats;
        const auto bins = filter_n_only(unfiltered_bins, fai, options.threads, kernel_stats);

        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write interval list: " + options.output);
        output << "@HD\tVN:1.6\tSO:coordinate\n";
        for (int i = 0; i < faidx_nseq(fai); ++i) {
            const auto contig = std::string(faidx_iseq(fai, i));
            output << "@SQ\tSN:" << contig << "\tLN:" << faidx_seq_len64(fai, contig.c_str()) << '\n';
        }
        for (const auto& interval : bins)
            output << interval.contig << '\t' << interval.start << '\t' << interval.end << "\t+\t.\n";
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot finalize interval list: " + options.output);
        if (!options.manifest.empty()) {
            std::ofstream manifest(options.manifest);
            if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
            manifest << "{\"schema_version\":1,\"tool\":\"PreprocessIntervals\",\"implementation\":\"fastgatk-preprocess-intervals\",\"status\":\"prototype\",\"reference\":\""
                     << json_escape(options.reference) << "\",\"input_intervals\":" << normalized.size()
                     << ",\"excluded_intervals\":" << normalized_exclusions.size()
                     << ",\"intervals_after_exclusion\":" << after_exclusion.size()
                     << ",\"unfiltered_bins\":" << unfiltered_bins.size() << ",\"output_intervals\":" << bins.size()
                     << ",\"padding\":" << options.padding << ",\"exclusion_padding\":" << options.exclusion_padding
                     << ",\"bin_length\":" << options.bin_length
                     << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\",\"output\":\""
                     << json_escape(options.output) << "\",\"telemetry\":{\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                     << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_batches\":" << kernel_stats.batches
                     << ",\"kernel_records\":" << kernel_stats.records
                     << ",\"kernel_prepare_seconds\":" << kernel_stats.prepare_seconds
                     << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds << "}}\n";
        }
        std::cout << "{\"tool\":\"PreprocessIntervals\",\"status\":\"prototype\",\"input_intervals\":"
                  << normalized.size() << ",\"unfiltered_bins\":" << unfiltered_bins.size()
                  << ",\"excluded_intervals\":" << normalized_exclusions.size()
                  << ",\"intervals_after_exclusion\":" << after_exclusion.size()
                  << ",\"output_intervals\":" << bins.size() << ",\"kernel_execution_space\":\"" << Kokkos::DefaultExecutionSpace::name()
                  << "\",\"kernel_batches\":" << kernel_stats.batches << ",\"kernel_records\":" << kernel_stats.records
                  << ",\"kernel_execute_seconds\":" << kernel_stats.execute_seconds << ",\"output\":\""
                  << json_escape(options.output) << "\"}\n";
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
