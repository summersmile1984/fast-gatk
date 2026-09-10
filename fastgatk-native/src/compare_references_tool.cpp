#include "fastgatk/core/plan.hpp"
#include "fastgatk/reference_io.hpp"

#include <Kokkos_Core.hpp>

#include <htslib/faidx.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Options {
    std::string reference;
    std::vector<std::string> references_to_compare;
    std::string output;
    std::string manifest;
    std::string md5_mode = "RECALCULATE_IF_MISSING";
    std::string base_mode = "NO_BASE_COMPARISON";
    std::string base_output;
    bool display_by_name = false;
    bool only_differing = false;
    int threads = 1;
};

std::string inline_value(const std::string& argument, std::string_view name) {
    const std::string prefix = std::string(name) + "=";
    return argument.rfind(prefix, 0) == 0 ? argument.substr(prefix.size()) : std::string{};
}

std::string require_value(int& index, int argc, char** argv, const std::string& argument,
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
            std::cout << "fastgatk-compare-references (GATK-compatible native contract-compatible)\n"
                         "  -R, --reference FILE             primary reference FASTA (+ .fai)\n"
                         "  -refcomp, --references-to-compare FILE  comparison FASTA (repeatable)\n"
                         "  -O, --output FILE                reference comparison TSV\n"
                         "      --md5-calculation-mode MODE  USE_DICT|RECALCULATE_IF_MISSING|ALWAYS_RECALCULATE\n"
                         "      --display-sequences-by-name  print name-keyed analysis\n"
                         "      --display-only-differing-sequences\n"
                         "      --base-comparison MODE       NO_BASE_COMPARISON|FIND_SNPS_ONLY|FULL_ALIGNMENT\n"
                         "      --base-comparison-output DIR output directory for base comparison\n"
                         "      --threads N                  Kokkos execution threads\n"
                         "      --output-manifest FILE       OutputManifest JSON\n";
            std::exit(0);
        } else if (argument == "-R" || !inline_value(argument, "--reference").empty()) {
            options.reference = require_value(index, argc, argv, argument, "--reference", "-R");
        } else if (argument == "-refcomp" || argument == "--refcomp" || argument == "--references-to-compare" ||
                   !inline_value(argument, "--references-to-compare").empty() ||
                   !inline_value(argument, "--refcomp").empty()) {
            const auto name = argument.rfind("--refcomp", 0) == 0 ? "--refcomp" : "--references-to-compare";
            options.references_to_compare.push_back(require_value(index, argc, argv, argument, name, "-refcomp"));
        } else if (argument == "-O" || !inline_value(argument, "--output").empty()) {
            options.output = require_value(index, argc, argv, argument, "--output", "-O");
        } else if (argument == "--md5-calculation-mode" || !inline_value(argument, "--md5-calculation-mode").empty() ||
                   argument == "-md5-calculation-mode") {
            options.md5_mode = require_value(index, argc, argv, argument, "--md5-calculation-mode", "-md5-calculation-mode");
        } else if (argument == "--display-sequences-by-name") {
            options.display_by_name = true;
        } else if (argument == "--display-only-differing-sequences") {
            options.only_differing = true;
        } else if (argument == "--base-comparison" || !inline_value(argument, "--base-comparison").empty()) {
            options.base_mode = require_value(index, argc, argv, argument, "--base-comparison");
        } else if (argument == "--base-comparison-output" || !inline_value(argument, "--base-comparison-output").empty()) {
            options.base_output = require_value(index, argc, argv, argument, "--base-comparison-output");
        } else if (argument == "--threads" || !inline_value(argument, "--threads").empty()) {
            options.threads = std::stoi(require_value(index, argc, argv, argument, "--threads"));
        } else if (argument == "--output-manifest" || argument == "--manifest" ||
                   !inline_value(argument, "--output-manifest").empty() || !inline_value(argument, "--manifest").empty()) {
            const auto name = argument.rfind("--manifest", 0) == 0 ? "--manifest" : "--output-manifest";
            options.manifest = require_value(index, argc, argv, argument, name);
        } else if (argument == "--quiet" || argument.rfind("--verbosity", 0) == 0 ||
                   argument.rfind("--java-options", 0) == 0 || argument.rfind("--seconds-between-progress-updates", 0) == 0) {
            if (argument != "--quiet") {
                const auto name = argument.rfind("--java-options", 0) == 0 ? "--java-options" :
                    (argument.rfind("--verbosity", 0) == 0 ? "--verbosity" : "--seconds-between-progress-updates");
                (void)require_value(index, argc, argv, argument, name);
            }
        } else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    if (options.reference.empty()) throw std::invalid_argument("-R/--reference is required");
    if (options.references_to_compare.empty()) throw std::invalid_argument("-refcomp/--references-to-compare is required");
    if (options.md5_mode != "USE_DICT" && options.md5_mode != "RECALCULATE_IF_MISSING" && options.md5_mode != "ALWAYS_RECALCULATE")
        throw std::invalid_argument("unsupported --md5-calculation-mode: " + options.md5_mode);
    if (options.base_mode != "NO_BASE_COMPARISON" && options.base_mode != "FIND_SNPS_ONLY" && options.base_mode != "FULL_ALIGNMENT")
        throw std::invalid_argument("unsupported --base-comparison: " + options.base_mode);
    if (options.base_mode != "NO_BASE_COMPARISON") {
        if (options.references_to_compare.size() != 1)
            throw std::invalid_argument("base comparison requires exactly two references");
        if (options.base_output.empty()) throw std::invalid_argument("--base-comparison-output is required for base comparison");
        if (!std::filesystem::is_directory(options.base_output))
            throw std::invalid_argument("base comparison output directory does not exist: " + options.base_output);
    }
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

struct Row {
    std::string md5;
    std::int64_t length = 0;
    std::vector<std::string> names;
};

struct Table {
    std::vector<Row> rows;
    std::vector<fastgatk::reference::FastaInfo> references;
    std::unordered_map<std::string, std::size_t> row_index;
};

Table build_table(const Options& options) {
    Table table;
    table.references.reserve(options.references_to_compare.size() + 1);
    table.references.push_back(fastgatk::reference::load_fasta_info(options.reference, options.md5_mode));
    for (const auto& path : options.references_to_compare)
        table.references.push_back(fastgatk::reference::load_fasta_info(path, options.md5_mode));
    for (const auto& reference : table.references) {
        for (const auto& sequence : reference.sequences) {
            auto [iterator, inserted] = table.row_index.emplace(sequence.md5, table.rows.size());
            if (inserted) table.rows.push_back({sequence.md5, sequence.length,
                                                std::vector<std::string>(table.references.size(), "---")});
            auto& row = table.rows[iterator->second];
            row.names[&reference - table.references.data()] = sequence.name;
        }
    }
    return table;
}

std::string table_text(const Table& table) {
    std::ostringstream out;
    out << "MD5\tLength";
    for (const auto& reference : table.references)
        out << '\t' << fastgatk::reference::basename(reference.path);
    out << '\n';
    for (const auto& row : table.rows) {
        out << row.md5 << '\t' << row.length;
        for (const auto& name : row.names) out << '\t' << name;
        out << '\n';
    }
    return out.str();
}

enum class Status { Exact, DifferentNames, DifferentSequence, DifferentPresence, Superset, Subset };

std::unordered_map<std::string, fastgatk::reference::SequenceInfo> sequence_map(const fastgatk::reference::FastaInfo& reference) {
    std::unordered_map<std::string, fastgatk::reference::SequenceInfo> result;
    for (const auto& sequence : reference.sequences) result.emplace(sequence.name, sequence);
    return result;
}

std::vector<Status> analyze_pair(const Table& table, std::size_t left, std::size_t right) {
    bool exact = true;
    bool differing_names = false;
    for (const auto& row : table.rows) {
        const auto& a = row.names[left];
        const auto& b = row.names[right];
        if (a != b) exact = false;
        if (a != "---" && b != "---" && a != b) differing_names = true;
    }
    const auto left_map = sequence_map(table.references[left]);
    const auto right_map = sequence_map(table.references[right]);
    bool differing_sequence = false;
    bool differing_presence = false;
    bool left_subset = true;
    bool right_subset = true;
    for (const auto& [name, sequence] : left_map) {
        const auto it = right_map.find(name);
        if (it == right_map.end()) {
            left_subset = false;
        } else if (it->second.md5 != sequence.md5) {
            differing_sequence = true;
        }
    }
    for (const auto& [name, sequence] : right_map) {
        if (left_map.find(name) == left_map.end()) {
            right_subset = false;
        }
    }
    // A missing table cell is a true presence difference only when the
    // missing name is absent from the other dictionary.  If the same name is
    // present under another MD5, Java ReferenceSequenceTable reports
    // DIFFER_IN_SEQUENCE rather than also reporting DIFFER_IN_SEQUENCES_PRESENT.
    for (const auto& row : table.rows) {
        const auto& a = row.names[left];
        const auto& b = row.names[right];
        if (a == "---" && b != "---" && left_map.find(b) == left_map.end()) differing_presence = true;
        if (b == "---" && a != "---" && right_map.find(a) == right_map.end()) differing_presence = true;
    }
    std::vector<Status> statuses;
    if (exact) statuses.push_back(Status::Exact);
    else {
        if (differing_names) statuses.push_back(Status::DifferentNames);
        if (differing_sequence) statuses.push_back(Status::DifferentSequence);
        if (differing_presence) {
            if (!differing_names && left_subset && !right_subset) statuses.push_back(Status::Subset);
            else if (!differing_names && right_subset && !left_subset) statuses.push_back(Status::Superset);
            else statuses.push_back(Status::DifferentPresence);
        }
    }
    if (statuses.empty()) statuses.push_back(Status::DifferentPresence);
    return statuses;
}

const char* status_name(Status status) {
    switch (status) {
        case Status::Exact: return "EXACT_MATCH";
        case Status::DifferentNames: return "DIFFER_IN_SEQUENCE_NAMES";
        case Status::DifferentSequence: return "DIFFER_IN_SEQUENCE";
        case Status::DifferentPresence: return "DIFFER_IN_SEQUENCES_PRESENT";
        case Status::Superset: return "SUPERSET";
        case Status::Subset: return "SUBSET";
    }
    return "UNKNOWN";
}

std::string display_by_name(const Table& table, bool only_differing) {
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (const auto& reference : table.references)
        for (const auto& sequence : reference.sequences)
            if (seen.insert(sequence.name).second) names.push_back(sequence.name);
    std::ostringstream out;
    out << "*********************************************************\nName \tMD5 \tReference\n";
    for (const auto& name : names) {
        std::vector<const Row*> rows;
        for (const auto& row : table.rows)
            if (std::find(row.names.begin(), row.names.end(), name) != row.names.end()) rows.push_back(&row);
        if (only_differing && rows.size() <= 1) continue;
        out << name;
        for (const auto* row : rows) {
            out << "\n\t" << row->md5 << '\t';
            for (std::size_t index = 0; index < row->names.size(); ++index)
                if (row->names[index] == name) out << fastgatk::reference::basename(table.references[index].path) << '\t';
        }
        out << '\n';
    }
    return out.str();
}

struct BaseComparisonResult {
    std::string output;
    std::size_t compared_bases = 0;
    std::size_t mismatches = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
};

BaseComparisonResult find_snps(const Table& table, const Options& options) {
    const auto& reference1 = table.references[0];
    const auto& reference2 = table.references[1];
    faidx_t* fai1 = fai_load(reference1.path.c_str());
    faidx_t* fai2 = fai_load(reference2.path.c_str());
    if (!fai1 || !fai2) {
        if (fai1) fai_destroy(fai1);
        if (fai2) fai_destroy(fai2);
        throw std::runtime_error("BACKEND_UNAVAILABLE: base comparison requires readable FASTA indexes");
    }
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(options.threads);
    Kokkos::initialize(settings);
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    fastgatk::core::HostBatch host("compare-references-snps-v1");
    BaseComparisonResult result;
    std::ostringstream rows;
    rows << "Sequence Name\tPosition\t" << fastgatk::reference::basename(reference1.path)
         << '\t' << fastgatk::reference::basename(reference2.path) << '\n';
    std::size_t compared = 0;
    for (const auto& sequence : reference1.sequences) {
        const auto it = std::find_if(reference2.sequences.begin(), reference2.sequences.end(),
                                     [&](const auto& candidate) { return candidate.name == sequence.name; });
        if (it == reference2.sequences.end() || it->length != sequence.length || it->md5 == sequence.md5) continue;
        const auto bases1 = fastgatk::reference::fetch_sequence(fai1, sequence.name, 1, sequence.length);
        const auto bases2 = fastgatk::reference::fetch_sequence(fai2, sequence.name, 1, sequence.length);
        const std::size_t length = bases1.size();
        host.records += length;
        host.bytes += length * 2;
        std::vector<unsigned char> mismatches(length);
        {
            Kokkos::View<char*> d1("reference1_bases", length);
            Kokkos::View<char*> d2("reference2_bases", length);
            Kokkos::View<unsigned char*> dm("reference_mismatch", length);
            auto h1 = Kokkos::create_mirror_view(d1);
            auto h2 = Kokkos::create_mirror_view(d2);
            for (std::size_t index = 0; index < length; ++index) { h1(index) = bases1[index]; h2(index) = bases2[index]; }
            Kokkos::deep_copy(d1, h1);
            Kokkos::deep_copy(d2, h2);
            fastgatk::core::DeviceBatch<ExecSpace> device(length);
            device.bind("ref1", d1); device.bind("ref2", d2); device.bind("mismatch", dm);
            fastgatk::core::KernelPlan<ExecSpace> plan("compare-references-snps");
            plan.begin_prepare(host);
            ExecSpace().fence();
            plan.end_prepare(device);
            plan.begin_execute();
            Kokkos::parallel_for("compare_references_snps", Kokkos::RangePolicy<ExecSpace>(0, length),
                KOKKOS_LAMBDA(const std::size_t index) { dm(index) = d1(index) == d2(index) ? 0 : 1; });
            ExecSpace().fence();
            plan.end_execute();
            auto hm = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dm);
            for (std::size_t index = 0; index < length; ++index) mismatches[index] = hm(index);
            result.prepare_seconds += plan.telemetry().prepare_seconds;
            result.execute_seconds += plan.telemetry().execute_seconds;
        }
        for (std::size_t index = 0; index < length; ++index) {
            if (mismatches[index] == 0) continue;
            rows << sequence.name << '\t' << (index + 1) << '\t' << bases1[index] << '\t' << bases2[index] << '\n';
            ++result.mismatches;
        }
        compared += length;
    }
    Kokkos::finalize();
    fai_destroy(fai1);
    fai_destroy(fai2);
    result.compared_bases = compared;
    const auto filename = fastgatk::reference::basename(reference1.path) + "_" +
        fastgatk::reference::basename(reference2.path) + "_snps.tsv";
    result.output = (std::filesystem::path(options.base_output) / filename).string();
    std::ofstream output(result.output);
    if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write SNP table: " + result.output);
    output << rows.str();
    return result;
}

int main_impl(int argc, char** argv) {
    const auto options = parse_options(argc, argv);
    if (options.base_mode == "FULL_ALIGNMENT")
        throw std::invalid_argument("FULL_ALIGNMENT requires the external MUMmer backend; use --fallback");
    const auto table = build_table(options);
    const auto table_output = table_text(table);
    if (options.output.empty()) std::cout << table_output;
    else {
        std::ofstream output(options.output);
        if (!output) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write table: " + options.output);
        output << table_output;
    }
    const auto references = table.references.size();
    for (std::size_t left = 0; left < references; ++left) {
        for (std::size_t right = left + 1; right < references; ++right) {
            const auto statuses = analyze_pair(table, left, right);
            std::cout << "*********************************************************\nREFERENCE PAIR: "
                      << fastgatk::reference::basename(table.references[left].path) << ", "
                      << fastgatk::reference::basename(table.references[right].path) << "\nStatus:\n";
            for (const auto status : statuses) std::cout << '\t' << status_name(status) << '\n';
        }
    }
    if (options.display_by_name) std::cout << display_by_name(table, options.only_differing);
    BaseComparisonResult base_result;
    if (options.base_mode == "FIND_SNPS_ONLY") {
        base_result = find_snps(table, options);
        std::cerr << "{\"base_comparison\":\"FIND_SNPS_ONLY\",\"output\":\""
                  << json_escape(base_result.output) << "\",\"compared_bases\":" << base_result.compared_bases
                  << ",\"mismatches\":" << base_result.mismatches << "}\n";
    }
    if (!options.manifest.empty()) {
        std::ofstream manifest(options.manifest);
        if (!manifest) throw std::runtime_error("OUTPUT_CONTRACT_FAILURE: cannot write manifest: " + options.manifest);
        manifest << "{\"schema_version\":1,\"tool\":\"CompareReferences\",\"implementation\":\"fastgatk-compare-references\",\"status\":\"contract-compatible\",\"references\":"
                 << references << ",\"rows\":" << table.rows.size() << ",\"md5_calculation_mode\":\"" << options.md5_mode
                 << "\",\"base_comparison\":\"" << options.base_mode << "\",\"output\":\"" << json_escape(options.output)
                 << "\",\"telemetry\":{\"host_reference_io\":true,\"kernel_lifecycle\":\"HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect\",\"kernel_execution_space\":\""
                 << Kokkos::DefaultExecutionSpace::name() << "\",\"kernel_execution_policy\":\"RangePolicy\",\"kernel_records\":"
                 << base_result.compared_bases << ",\"kernel_prepare_seconds\":" << base_result.prepare_seconds
                 << ",\"kernel_execute_seconds\":" << base_result.execute_seconds << "}}\n";
    }
    std::cerr << "{\"tool\":\"CompareReferences\",\"status\":\"contract-compatible\",\"references\":" << references
              << ",\"rows\":" << table.rows.size() << ",\"base_comparison\":\"" << options.base_mode << "\"}\n";
    return 0;
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
