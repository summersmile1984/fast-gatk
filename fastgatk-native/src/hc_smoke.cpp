#include "fastgatk/core/plan.hpp"
#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/runtime/resource.hpp"

#include <Kokkos_Core.hpp>

#include <cstdint>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_ZLIB
#include <zlib.h>
#endif

namespace {

struct Options {
    std::string input;
    std::string reference;
    std::string region;
    std::string output;
    std::string manifest;
    std::string telemetry;
    std::size_t batch_records = 1024;
    int threads = 1;
    std::vector<std::string> compatibility_options;
};

std::string value(const std::string& arg, const char* name) {
    const std::string prefix = std::string(name) + "=";
    return arg.rfind(prefix, 0) == 0 ? arg.substr(prefix.size()) : std::string{};
}

bool is_option(const std::string& arg, const char* name) {
    return arg == name || !value(arg, name).empty();
}

std::string require_value(int& index, int argc, char** argv, const std::string& arg,
                          const char* long_name, const char* short_name = nullptr) {
    if (const auto inline_value = value(arg, long_name); !inline_value.empty()) return inline_value;
    if (short_name && arg == short_name && index + 1 < argc) return argv[++index];
    if (arg == long_name && index + 1 < argc) return argv[++index];
    throw std::invalid_argument(std::string("missing value for ") + long_name);
}

bool is_integer(const std::string& text) {
    if (text.empty()) return false;
    std::size_t offset = text[0] == '-' ? 1 : 0;
    if (offset == text.size()) return false;
    return std::all_of(text.begin() + static_cast<std::ptrdiff_t>(offset), text.end(),
                       [](unsigned char ch) { return std::isdigit(ch) != 0; });
}

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help" || arg == "-h") {
            std::cout << "fastgatk-hc-smoke [GATK-compatible HaplotypeCaller options]\n"
                         "  -I, --input FILE                 input SAM/BAM/CRAM\n"
                         "  -R, --reference FILE             reference FASTA (CRAM)\n"
                         "  -L, --intervals REGION           contig:start-end\n"
                         "  -O, --output FILE                VCF or JSON summary\n"
                         "      --output-manifest FILE       OutputManifest JSON\n"
                         "      --telemetry FILE             JSON telemetry summary\n"
                         "      --batch-records N             streaming batch size\n"
                         "      --threads N                  Kokkos worker count\n"
                         "Compatibility aliases include -I/-R/-L/-O, --interval, and\n"
                         "--native-pair-hmm-threads. This build is a smoke adapter and emits\n"
                         "a valid empty VCF; it is not yet a biological variant caller.\n";
            std::exit(0);
        } else if (arg == "--version") {
            std::cout << "fastgatk-hc-smoke 0.2 (compatibility smoke adapter)\n";
            std::exit(0);
        } else if (is_option(arg, "--input")) options.input = require_value(i, argc, argv, arg, "--input", "-I");
        else if (arg == "-I") options.input = require_value(i, argc, argv, arg, "--input", "-I");
        else if (is_option(arg, "--reference")) options.reference = require_value(i, argc, argv, arg, "--reference", "-R");
        else if (arg == "-R") options.reference = require_value(i, argc, argv, arg, "--reference", "-R");
        else if (is_option(arg, "--region")) options.region = require_value(i, argc, argv, arg, "--region");
        else if (is_option(arg, "--intervals")) options.region = require_value(i, argc, argv, arg, "--intervals", "-L");
        else if (arg == "-L") options.region = require_value(i, argc, argv, arg, "--intervals", "-L");
        else if (is_option(arg, "--interval")) options.region = require_value(i, argc, argv, arg, "--interval");
        else if (is_option(arg, "--output")) options.output = require_value(i, argc, argv, arg, "--output", "-O");
        else if (arg == "-O") options.output = require_value(i, argc, argv, arg, "--output", "-O");
        else if (is_option(arg, "--output-manifest")) options.manifest = require_value(i, argc, argv, arg, "--output-manifest");
        else if (is_option(arg, "--manifest")) options.manifest = require_value(i, argc, argv, arg, "--manifest");
        else if (is_option(arg, "--telemetry")) options.telemetry = require_value(i, argc, argv, arg, "--telemetry");
        else if (is_option(arg, "--batch-records")) options.batch_records = std::stoull(require_value(i, argc, argv, arg, "--batch-records"));
        else if (is_option(arg, "--threads")) options.threads = std::stoi(require_value(i, argc, argv, arg, "--threads"));
        else if (is_option(arg, "--native-pair-hmm-threads")) {
            const auto text = require_value(i, argc, argv, arg, "--native-pair-hmm-threads");
            if (!is_integer(text)) throw std::invalid_argument("invalid --native-pair-hmm-threads");
            options.threads = std::stoi(text);
            options.compatibility_options.push_back("--native-pair-hmm-threads=" + text);
        } else if (arg == "--native-pair-hmm-use-double-precision" ||
                   arg == "--disable-read-filter" || arg == "--disable-tool-default-read-filters" ||
                   arg == "--create-output-variant-index" || arg == "--add-output-vcf-command-line" ||
                   arg == "--disable-sequence-dictionary-validation" || arg == "--quiet") {
            std::string setting = "true";
            if (i + 1 < argc && argv[i + 1][0] != '-') setting = argv[++i];
            options.compatibility_options.push_back(arg + "=" + setting);
        } else if (is_option(arg, "--read-filter") || is_option(arg, "--tmp-dir") ||
                   is_option(arg, "--verbosity") || is_option(arg, "--java-options") ||
                   is_option(arg, "--seconds-between-progress-updates") ||
                   is_option(arg, "--interval-set-rule") || is_option(arg, "--interval-padding") ||
                   is_option(arg, "--bam-output") || is_option(arg, "--bam-writer-type") ||
                   is_option(arg, "--native-pair-hmm-log-performance")) {
            options.compatibility_options.push_back(arg + "=" + require_value(i, argc, argv, arg,
                arg.substr(0, arg.find('=' )).c_str()));
        } else if (arg == "--gvcf" || arg == "-ERC") {
            // Preserve the command-line contract but make the smoke limitation visible.
            options.compatibility_options.push_back(arg);
        } else throw std::runtime_error("unknown option: " + arg);
    }
    if (options.input.empty()) throw std::invalid_argument("--input is required");
    if (options.batch_records == 0 || options.threads < 1) throw std::invalid_argument("invalid batch/thread count");
    return options;
}

std::string json_escape(const std::string& text) {
    std::ostringstream escaped;
    for (const auto ch : text) {
        switch (ch) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default: escaped << ch; break;
        }
    }
    return escaped.str();
}

struct BatchStats {
    std::uint64_t reads = 0;
    std::uint64_t bases = 0;
    std::uint64_t quality_sum = 0;
    std::uint64_t mapq_sum = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

BatchStats process_batch(const fastgatk::io::ReadBatch& input) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::View<std::uint8_t*> bases("hc_bases", input.bases.size());
    Kokkos::View<std::uint8_t*> qualities("hc_qualities", input.qualities.size());
    Kokkos::View<std::uint8_t*> mapq("hc_mapq", input.mapq.size());
    Kokkos::View<std::uint32_t*> offsets("hc_offsets", input.offsets.size());
    auto bases_host = Kokkos::create_mirror_view(bases);
    auto qualities_host = Kokkos::create_mirror_view(qualities);
    auto mapq_host = Kokkos::create_mirror_view(mapq);
    auto offsets_host = Kokkos::create_mirror_view(offsets);
    for (std::size_t i = 0; i < input.bases.size(); ++i) {
        bases_host(i) = input.bases[i]; qualities_host(i) = input.qualities[i];
    }
    for (std::size_t i = 0; i < input.mapq.size(); ++i) mapq_host(i) = input.mapq[i];
    for (std::size_t i = 0; i < input.offsets.size(); ++i) offsets_host(i) = input.offsets[i];
    Kokkos::deep_copy(bases, bases_host); Kokkos::deep_copy(qualities, qualities_host);
    Kokkos::deep_copy(mapq, mapq_host); Kokkos::deep_copy(offsets, offsets_host);

    fastgatk::core::HostBatch host_batch("hc-smoke-read-batch-v1");
    host_batch.records = input.records();
    host_batch.bytes = input.bases.size() + input.qualities.size() + input.mapq.size() +
                       input.offsets.size() * sizeof(std::uint32_t);
    fastgatk::core::KernelPlan<ExecSpace> plan("hc-smoke-pileup");
    plan.begin_prepare(host_batch);
    fastgatk::core::DeviceBatch<ExecSpace> device_batch(input.records());
    device_batch.bind("bases", bases); device_batch.bind("qualities", qualities);
    device_batch.bind("mapq", mapq); device_batch.bind("offsets", offsets);
    Kokkos::fence();
    plan.end_prepare(device_batch);
    plan.begin_execute();
    std::uint64_t bases_sum = 0, quality_sum = 0, mapq_sum = 0;
    Kokkos::parallel_reduce("hc_smoke_bases", Kokkos::RangePolicy<ExecSpace>(0, bases.extent(0)),
        KOKKOS_LAMBDA(const std::size_t i, std::uint64_t& value) { value += bases(i) != static_cast<std::uint8_t>('N'); }, bases_sum);
    Kokkos::parallel_reduce("hc_smoke_quality", Kokkos::RangePolicy<ExecSpace>(0, qualities.extent(0)),
        KOKKOS_LAMBDA(const std::size_t i, std::uint64_t& value) { value += qualities(i); }, quality_sum);
    Kokkos::parallel_reduce("hc_smoke_mapq", Kokkos::RangePolicy<ExecSpace>(0, mapq.extent(0)),
        KOKKOS_LAMBDA(const std::size_t i, std::uint64_t& value) { value += mapq(i); }, mapq_sum);
    Kokkos::fence();
    plan.end_execute();
    return BatchStats{input.records(), bases_sum, quality_sum, mapq_sum,
                      plan.telemetry().prepare_seconds, plan.telemetry().execute_seconds};
}

std::string make_json(const Options& options, const fastgatk::io::HtsReader& reader,
                      const BatchStats& stats, std::size_t batches,
                      double prepare_seconds, double execute_seconds) {
    std::ostringstream output;
    output << "{\"schema_version\":1,\"tool\":\"HaplotypeCaller\",\"implementation\":\"fastgatk-hc-smoke\",\"status\":\"smoke\""
           << ",\"htslib_backend\":\"" << json_escape(reader.backend_description()) << "\""
           << ",\"execution_space\":\"" << Kokkos::DefaultExecutionSpace::name() << "\""
           << ",\"input\":\"" << json_escape(options.input) << "\""
           << ",\"reference\":\"" << json_escape(options.reference) << "\""
           << ",\"interval\":\"" << json_escape(options.region) << "\""
           << ",\"batches\":" << batches << ",\"reads\":" << stats.reads
           << ",\"bases\":" << stats.bases << ",\"quality_sum\":" << stats.quality_sum
           << ",\"mapq_sum\":" << stats.mapq_sum
           << ",\"prepare_seconds\":" << std::setprecision(12) << prepare_seconds
           << ",\"execute_seconds\":" << execute_seconds << "}\n";
    return output.str();
}

bool has_suffix(const std::string& path, const char* suffix) {
    const std::string ending(suffix);
    return path.size() >= ending.size() &&
           path.compare(path.size() - ending.size(), ending.size(), ending) == 0;
}

bool is_vcf_output(const std::string& path) {
    return has_suffix(path, ".vcf") || has_suffix(path, ".vcf.gz") ||
           has_suffix(path, ".g.vcf") || has_suffix(path, ".g.vcf.gz");
}

std::string make_vcf(const fastgatk::io::HtsReader& reader, const BatchStats& stats) {
    std::ostringstream vcf;
    vcf << "##fileformat=VCFv4.2\n"
        << "##source=fastgatk-hc-smoke\n"
        << "##fastgatk_status=smoke\n"
        << "##fastgatk_reads=" << stats.reads << "\n"
        << "##fastgatk_bases=" << stats.bases << "\n";
    const auto& header = reader.header();
    for (std::size_t i = 0; i < header.contigs.size(); ++i)
        vcf << "##contig=<ID=" << header.contigs[i] << ",length=" << header.contig_lengths[i] << ">\n";
    vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    return vcf.str();
}

void write_artifact(const std::string& path, const std::string& text) {
    if (has_suffix(path, ".gz")) {
#if FASTGATK_HAS_ZLIB
        gzFile file = gzopen(path.c_str(), "wb");
        if (!file) throw std::runtime_error("cannot open compressed output: " + path);
        const int written = gzwrite(file, text.data(), static_cast<unsigned int>(text.size()));
        const int close_status = gzclose(file);
        if (written != static_cast<int>(text.size()) || close_status != Z_OK)
            throw std::runtime_error("cannot write compressed output: " + path);
#else
        throw std::runtime_error("compressed output requires zlib: " + path);
#endif
    } else {
        std::ofstream file(path);
        if (!file) throw std::runtime_error("cannot open output: " + path);
        file << text;
        if (!file) throw std::runtime_error("cannot write output: " + path);
    }
}

std::string make_manifest(const Options& options, const fastgatk::io::HtsReader& reader,
                          const BatchStats& stats, std::size_t batches,
                          double prepare_seconds, double execute_seconds) {
    const std::string summary = make_json(options, reader, stats, batches, prepare_seconds, execute_seconds);
    std::ostringstream manifest;
    manifest << "{\"schema_version\":1,\"tool\":\"HaplotypeCaller\","
             << "\"implementation\":\"fastgatk-hc-smoke\",\"status\":\"smoke\","
             << "\"primary_output\":\"" << json_escape(options.output) << "\","
             << "\"primary_output_kind\":\"" << (is_vcf_output(options.output) ? "vcf" : "json") << "\","
             << "\"input\":\"" << json_escape(options.input) << "\","
             << "\"reference\":\"" << json_escape(options.reference) << "\","
             << "\"interval\":\"" << json_escape(options.region) << "\","
             << "\"outputs\":[{\"path\":\"" << json_escape(options.output)
             << "\",\"kind\":\"" << (is_vcf_output(options.output) ? "vcf" : "json")
             << "\",\"complete\":true}],"
             << "\"compatibility\":{\"gatk_parameter_aliases\":true,"
             << "\"biological_variant_calling\":false,\"bit_identical_to_gatk\":false},"
             << "\"telemetry\":" << summary << "}\n";
    return manifest.str();
}

}  // namespace

int main(int argc, char** argv) {
    bool initialized = false;
    try {
        Options options = parse(argc, argv);
        const auto resources = fastgatk::runtime::ResourceSnapshot::probe();
        Kokkos::InitializationSettings settings;
        settings.set_num_threads(resources.effective_threads(static_cast<std::size_t>(options.threads)));
        Kokkos::initialize(settings);
        initialized = true;
        fastgatk::io::HtsReader reader(options.input, options.reference, options.region,
                                       options.batch_records, {}, true);
        fastgatk::io::ReadBatch batch;
        BatchStats totals;
        std::size_t batches = 0;
        double prepare_seconds = 0.0, execute_seconds = 0.0;
        while (reader.next(batch)) {
            const auto stats = process_batch(batch);
            totals.reads += stats.reads; totals.bases += stats.bases;
            totals.quality_sum += stats.quality_sum; totals.mapq_sum += stats.mapq_sum;
            ++batches;
            // The per-batch plan is intentionally visible in this smoke adapter;
            // production HC will pool it across AssemblyRegions.
            prepare_seconds += stats.prepare_seconds;
            execute_seconds += stats.execute_seconds;
        }
        const std::string json = make_json(options, reader, totals, batches, prepare_seconds, execute_seconds);
        if (!options.output.empty()) {
            if (is_vcf_output(options.output)) write_artifact(options.output, make_vcf(reader, totals));
            else write_artifact(options.output, json);
        }
        if (!options.telemetry.empty()) write_artifact(options.telemetry, json);
        if (options.manifest.empty() && is_vcf_output(options.output)) options.manifest = options.output + ".manifest.json";
        if (!options.manifest.empty()) write_artifact(options.manifest,
            make_manifest(options, reader, totals, batches, prepare_seconds, execute_seconds));
        std::cout << json;
        Kokkos::finalize();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        if (initialized && Kokkos::is_initialized()) Kokkos::finalize();
        return 2;
    }
}
