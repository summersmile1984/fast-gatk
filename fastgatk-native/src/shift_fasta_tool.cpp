#include "fastgatk/core/plan.hpp"
#include "fastgatk/reference_io.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::string reference;
    std::string output;
    std::string shift_back_output;
    std::string interval_file_name;
    std::vector<int> offsets;
    int line_width = 60;
    int threads = 1;
    std::string manifest;
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

void append_offsets(const std::string& value, std::vector<int>& offsets) {
    std::string token;
    std::istringstream stream(value);
    while (std::getline(stream, token, ',')) {
        if (token.empty()) throw std::invalid_argument("BAD_INPUT: empty shift offset");
        offsets.push_back(std::stoi(token));
    }
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument(argv[index]);
        if (argument == "--help" || argument == "-h") {
            std::cout << "fastgatk-shift-fasta (GATK-compatible native prototype)\n"
                         "  -R, --reference FILE          input circular FASTA (+ .fai)\n"
                         "  -O, --output FILE             shifted FASTA\n"
                         "      --shift-back-output FILE  UCSC chain output\n"
                         "      --shift-offset-list LIST  comma-separated offsets\n"
                         "      --interval-file-name BASE write .intervals/.shifted.intervals\n"
                         "      --line-width N            bases per sequence line (default 60)\n"
                         "      --output-manifest FILE    OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || argument == "--reference" || !inline_value(argument, "--reference").empty()) {
            options.reference = value_for(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-O" || argument == "--output" || !inline_value(argument, "--output").empty()) {
            options.output = value_for(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--shift-back-output" || !inline_value(argument, "--shift-back-output").empty()) {
            options.shift_back_output = value_for(index, argc, argv, argument, "--shift-back-output");
        } else if (argument == "--shift-offset-list" || !inline_value(argument, "--shift-offset-list").empty()) {
            append_offsets(value_for(index, argc, argv, argument, "--shift-offset-list"), options.offsets);
        } else if (argument == "--interval-file-name" || !inline_value(argument, "--interval-file-name").empty()) {
            options.interval_file_name = value_for(index, argc, argv, argument, "--interval-file-name");
        } else if (argument == "--line-width" || !inline_value(argument, "--line-width").empty()) {
            options.line_width = std::stoi(value_for(index, argc, argv, argument, "--line-width"));
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(value_for(index, argc, argv, argument, "--threads"));
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
    if (options.output.empty()) throw std::invalid_argument("-O/--output is required");
    if (options.shift_back_output.empty()) throw std::invalid_argument("--shift-back-output is required");
    if (options.line_width < 1) throw std::invalid_argument("--line-width must be positive");
    if (options.threads < 1) throw std::invalid_argument("--threads must be positive");
    return options;
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

struct KernelStats {
    std::size_t records = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

std::string rotate_kokkos(const std::string& sequence, int offset, int threads, KernelStats& stats) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    const std::size_t length = sequence.size();
    const std::size_t split = static_cast<std::size_t>(offset);
    Kokkos::View<unsigned char*> input("shift_fasta_input", length);
    Kokkos::View<unsigned char*> output("shift_fasta_output", length);
    auto host = Kokkos::create_mirror_view(input);
    for (std::size_t index = 0; index < length; ++index)
        host(index) = static_cast<unsigned char>(std::toupper(static_cast<unsigned char>(sequence[index])));
    fastgatk::core::HostBatch batch("shift-fasta-v1");
    batch.records = length;
    batch.bytes = length;
    fastgatk::core::KernelPlan<ExecSpace> plan("shift-fasta-rotate");
    plan.begin_prepare(batch);
    Kokkos::deep_copy(input, host);
    ExecSpace().fence();
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(length);
    device_batch.bind("input", input);
    device_batch.bind("output", output);
    plan.end_prepare(device_batch);
    plan.begin_execute();
    Kokkos::parallel_for("shift_fasta_rotate", Kokkos::RangePolicy<ExecSpace>(0, length),
        KOKKOS_LAMBDA(const std::size_t index) {
            const std::size_t source = (index + split) % length;
            output(index) = input(source);
        });
    ExecSpace().fence();
    plan.end_execute();
    auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
    std::string rotated(length, '\0');
    for (std::size_t index = 0; index < length; ++index) rotated[index] = static_cast<char>(result(index));
    stats.records += length;
    stats.prepare_seconds += plan.telemetry().prepare_seconds;
    stats.execute_seconds += plan.telemetry().execute_seconds;
    (void)threads;
    return rotated;
}

std::string md5_for_bytes(const std::string& sequence) {
    hts_md5_context* context = hts_md5_init();
    if (!context) throw std::runtime_error("BACKEND_UNAVAILABLE: HTSlib MD5 context allocation failed");
    hts_md5_update(context, sequence.data(), static_cast<unsigned long>(sequence.size()));
    unsigned char digest[16]{};
    hts_md5_final(digest, context);
    hts_md5_destroy(context);
    char hex[33]{};
    hts_md5_hex(hex, digest);
    return std::string(hex);
}

void write_dictionary(const std::string& fasta, const std::vector<std::pair<std::string, std::pair<std::int64_t, std::string>>>& sequences) {
    std::ofstream dictionary(fastgatk::reference::dictionary_path(fasta));
    if (!dictionary) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write sequence dictionary");
    dictionary << "@HD\tVN:1.6\n";
    for (const auto& [name, info] : sequences)
        dictionary << "@SQ\tSN:" << name << "\tLN:" << info.first << "\tM5:" << info.second << "\n";
}

void write_sequence(std::ofstream& output, const std::string& name, const std::string& sequence, int line_width) {
    output << ">" << name << "\n";
    for (std::size_t offset = 0; offset < sequence.size(); offset += static_cast<std::size_t>(line_width))
        output.write(sequence.data() + static_cast<std::streamoff>(offset),
                     static_cast<std::streamsize>(std::min<std::size_t>(line_width, sequence.size() - offset))) << '\n';
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    faidx_t* fai = fai_load(options.reference.c_str());
    if (!fai) throw std::runtime_error("BACKEND_UNAVAILABLE: reference FASTA requires a readable .fai: " + options.reference);
    try {
        const int count = faidx_nseq(fai);
        std::int64_t reference_length = 0;
        for (int tid = 0; tid < count; ++tid) {
            const std::string name = faidx_iseq(fai, tid);
            const auto contig_length = faidx_seq_len64(fai, name.c_str());
            if (contig_length < 0 || reference_length > std::numeric_limits<std::int64_t>::max() - contig_length)
                throw std::invalid_argument("BAD_INPUT: reference length is too long");
            reference_length += contig_length;
        }
        if (reference_length > std::numeric_limits<int>::max())
            throw std::invalid_argument("BAD_INPUT: reference length is too long");
        if (!options.offsets.empty() && static_cast<int>(options.offsets.size()) != count)
            throw std::invalid_argument("BAD_INPUT: shift offset list size must equal reference contig count");
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(options.threads);
        Kokkos::initialize(settings);
        bool initialized = true;
        KernelStats stats;
        std::ofstream output(options.output, std::ios::binary);
        std::ofstream chain(options.shift_back_output, std::ios::binary);
        if (!output || !chain) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create ShiftFasta output");
        std::ofstream regular_intervals;
        std::ofstream shifted_intervals;
        if (!options.interval_file_name.empty()) {
            regular_intervals.open(options.interval_file_name + ".intervals");
            shifted_intervals.open(options.interval_file_name + ".shifted.intervals");
            if (!regular_intervals || !shifted_intervals)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot create ShiftFasta interval sidecars");
        }
        int chain_id = 1;
        int shifted_contigs = 0;
        std::vector<std::pair<std::string, std::pair<std::int64_t, std::string>>> dictionary_entries;
        try {
            for (int tid = 0; tid < count; ++tid) {
                const std::string name = faidx_iseq(fai, tid);
                const auto length64 = static_cast<std::int64_t>(faidx_seq_len64(fai, name.c_str()));
                if (length64 > std::numeric_limits<int>::max()) throw std::invalid_argument("BAD_INPUT: reference length is too long");
                const int length = static_cast<int>(length64);
                const int offset = options.offsets.empty() ? length / 2 : options.offsets[static_cast<std::size_t>(tid)];
                if (offset <= 0 || offset >= length) continue;
                hts_pos_t fetched = 0;
                char* raw = faidx_fetch_seq64(fai, name.c_str(), 0, length - 1, &fetched);
                if (!raw || fetched != length) {
                    free(raw);
                    throw std::runtime_error("BAD_INPUT: failed to fetch reference contig: " + name);
                }
                const std::string sequence(raw, static_cast<std::size_t>(fetched));
                free(raw);
                const auto shifted = rotate_kokkos(sequence, offset, options.threads, stats);
                write_sequence(output, name, shifted, options.line_width);
                dictionary_entries.push_back({name, {length, md5_for_bytes(shifted)}});
                const int shift_back_offset = length - offset;
                chain << "chain\t" << shift_back_offset << "\t" << name << "\t" << length
                      << "\t+\t0\t" << shift_back_offset << "\t" << name << "\t" << length
                      << "\t+\t" << offset << "\t" << length << "\t" << chain_id++ << "\n"
                      << shift_back_offset << "\n\n";
                chain << "chain\t" << offset - 1 << "\t" << name << "\t" << length
                      << "\t+\t" << shift_back_offset << "\t" << length << "\t" << name << "\t" << length
                      << "\t+\t0\t" << offset << "\t" << chain_id++ << "\n"
                      << offset << "\n\n";
                if (regular_intervals) {
                    const int interval_start = offset / 2;
                    const int interval_end = interval_start + length / 2 - 1;
                    regular_intervals << name << ":" << interval_start << "-" << interval_end << "\n";
                    shifted_intervals << name << ":" << interval_start << "-" << interval_end + length % 2 << "\n";
                }
                ++shifted_contigs;
            }
            output.flush();
            chain.flush();
            if (!output || !chain) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed finalizing ShiftFasta output");
            if (fai_build(options.output.c_str()) != 0)
                throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: failed creating shifted FASTA index");
            write_dictionary(options.output, dictionary_entries);
            if (!options.manifest.empty()) {
                std::ofstream manifest(options.manifest);
                if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest");
                manifest << "{\"schema_version\":1,\"tool\":\"ShiftFasta\",\"implementation\":\"fastgatk-shift-fasta\",\"status\":\"prototype\",\"execution_space\":\""
                         << Kokkos::DefaultExecutionSpace::name() << "\",\"reference\":\"" << json_escape(options.reference)
                         << "\",\"shifted_contigs\":" << shifted_contigs << ",\"records\":" << stats.records
                         << ",\"prepare_seconds\":" << stats.prepare_seconds << ",\"execute_seconds\":" << stats.execute_seconds
                         << ",\"output\":\"" << json_escape(options.output) << "\",\"shift_back_output\":\""
                         << json_escape(options.shift_back_output) << "\",\"dictionary\":\""
                         << json_escape(fastgatk::reference::dictionary_path(options.output)) << "\"}\n";
            }
            Kokkos::finalize();
            initialized = false;
        } catch (...) {
            if (initialized) Kokkos::finalize();
            throw;
        }
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
