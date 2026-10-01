// fastgatk-funcotator — GATK-compatible Funcotator (VCF output).
//
// Mirrors GATK's Funcotator on the gencode + simpleXSV surface:
//   * data-source discovery <ds>/<Name>/<hg19|hg38>/<name>.config with the
//     required key validation of DataSourceUtils (gencode mandatory),
//   * gencode funcotations (GTF + transcript FASTA) with CANONICAL/ALL
//     transcript selection and the 22-field Gencode_<version>_* layout,
//   * simpleXSV funcotations keyed by GENE_NAME with the
//     <name>_<column> field naming,
//   * FUNCOTATION VCF rendering: Number=A, per-ALT `[block#block]`,
//     fields `|`-joined, values sanitized per FuncotatorUtils,
//   * --annotation-default / --annotation-override, and the
//     "Funcotator Version" header used for the re-annotation refusal.
//
// Determinism: everything is Host-side and ordered; outputs are byte-stable.

#include "optional_boolean.hpp"

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::invalid_argument(message);
}

struct Options {
    std::string variants;
    std::string output;
    std::string reference;
    std::string ref_version;                 // hg19 | hg38
    std::vector<std::string> data_sources;   // repeatable
    std::string output_format = "VCF";
    std::string transcript_selection_mode = "CANONICAL";
    std::vector<std::string> annotation_defaults;
    std::vector<std::string> annotation_overrides;
};

bool is_option(const std::string& argument, const char* name) {
    return argument == name ||
        argument.rfind(std::string(name) + "=", 0) == 0;
}

std::string option_value(int& index, int argc, char** argv,
                         const std::string& argument, const char* name) {
    const std::string prefix = std::string(name) + "=";
    if (argument.rfind(prefix, 0) == 0) return argument.substr(prefix.size());
    if (index + 1 >= argc) fail(std::string("missing value for ") + name);
    return argv[++index];
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            std::cout <<
                "fastgatk-funcotator (GATK-compatible native prototype)\n"
                "usage: fastgatk-funcotator -V <calls.vcf> -O <out.vcf> -R <ref.fa>\n"
                "  --ref-version hg19|hg38        reference version (required)\n"
                "  --data-sources-path DIR        data source root (repeatable)\n"
                "  --output-file-format VCF       output format (VCF)\n"
                "  --transcript-selection-mode M  CANONICAL|ALL (default CANONICAL)\n"
                "  --annotation-default K:V       default funcotation value\n"
                "  --annotation-override K:V      override funcotation value\n";
            std::exit(0);
        }
        if (is_option(argument, "-V") || is_option(argument, "--variant"))
            options.variants = option_value(index, argc, argv, argument,
                argument.rfind("-V", 0) == 0 ? "-V" : "--variant");
        else if (is_option(argument, "-O") || is_option(argument, "--output"))
            options.output = option_value(index, argc, argv, argument,
                argument.rfind("-O", 0) == 0 ? "-O" : "--output");
        else if (is_option(argument, "-R") || is_option(argument, "--reference"))
            options.reference = option_value(index, argc, argv, argument,
                argument.rfind("-R", 0) == 0 ? "-R" : "--reference");
        else if (is_option(argument, "--ref-version"))
            options.ref_version = option_value(index, argc, argv, argument, "--ref-version");
        else if (is_option(argument, "--data-sources-path"))
            options.data_sources.push_back(
                option_value(index, argc, argv, argument, "--data-sources-path"));
        else if (is_option(argument, "--output-file-format"))
            options.output_format = option_value(index, argc, argv, argument,
                "--output-file-format");
        else if (is_option(argument, "--transcript-selection-mode"))
            options.transcript_selection_mode = option_value(index, argc, argv, argument,
                "--transcript-selection-mode");
        else if (is_option(argument, "--annotation-default"))
            options.annotation_defaults.push_back(
                option_value(index, argc, argv, argument, "--annotation-default"));
        else if (is_option(argument, "--annotation-override"))
            options.annotation_overrides.push_back(
                option_value(index, argc, argv, argument, "--annotation-override"));
        else
            fail("unknown option: " + argument);
    }
    if (options.variants.empty()) fail("--variant is required");
    if (options.output.empty()) fail("--output is required");
    if (options.reference.empty()) fail("--reference is required");
    if (options.ref_version != "hg19" && options.ref_version != "hg38")
        fail("--ref-version must be hg19 or hg38");
    if (options.data_sources.empty()) fail("--data-sources-path is required");
    if (options.output_format != "VCF")
        fail("unsupported --output-file-format (supported: VCF)");
    if (options.transcript_selection_mode != "CANONICAL" &&
        options.transcript_selection_mode != "ALL")
        fail("unsupported --transcript-selection-mode (supported: CANONICAL, ALL)");
    return options;
}

// ---------------------------------------------------------------- FASTA ----

struct Fasta {
    faidx_t* index = nullptr;
    explicit Fasta(const std::string& path) {
        index = fai_load(path.c_str());
        if (index == nullptr) fail("cannot load reference index: " + path);
    }
    ~Fasta() { if (index != nullptr) fai_destroy(index); }
    std::string fetch(const std::string& contig, std::int64_t start_1based,
                      std::int64_t end_1based) const {
        int length = 0;
        char* bases = faidx_fetch_seq(index, contig.c_str(),
            static_cast<int>(start_1based - 1), static_cast<int>(end_1based - 1),
            &length);
        if (bases == nullptr) return {};
        std::string sequence(bases, static_cast<std::size_t>(length));
        free(bases);
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),
            [](unsigned char character) {
                return static_cast<char>(std::toupper(character));
            });
        return sequence;
    }
};

// --------------------------------------------------------------- config ----

std::map<std::string, std::string> load_config(const std::string& path) {
    std::map<std::string, std::string> values;
    std::ifstream input(path);
    if (!input) fail("cannot read config: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto equals = line.find('=');
        if (equals == std::string::npos) continue;
        std::string key = line.substr(0, equals);
        std::string value = line.substr(equals + 1);
        auto trim = [](std::string& text) {
            while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front())))
                text.erase(text.begin());
            while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())))
                text.pop_back();
        };
        trim(key);
        trim(value);
        values[key] = value;
    }
    return values;
}

struct DataSource {
    std::string name;
    std::string version;
    std::string type;
    std::string src_file;
    std::string ncbi_build;
    std::string gencode_fasta;
    std::string xsv_key;
    std::size_t xsv_key_column = 0;
    std::string xsv_delimiter = "\t";
    std::vector<std::pair<std::string, std::string>> xsv_rows;  // key -> raw line
    std::vector<std::string> xsv_header;
    std::string directory;
};

std::string config_value(const std::map<std::string, std::string>& config,
                         const std::string& key, bool required,
                         const std::string& context) {
    const auto found = config.find(key);
    if (found == config.end() || found->second.empty()) {
        if (required) fail("missing required config key '" + key + "' in " + context);
        return {};
    }
    return found->second;
}

std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    while (start <= text.size()) {
        const auto next = text.find(delimiter, start);
        parts.push_back(text.substr(start, next == std::string::npos
            ? std::string::npos : next - start));
        if (next == std::string::npos) break;
        start = next + 1;
    }
    return parts;
}

std::vector<DataSource> discover_data_sources(const Options& options) {
    std::vector<DataSource> sources;
    bool found_gencode = false;
    for (const auto& root : options.data_sources) {
        const std::filesystem::path root_path(root);
        if (!std::filesystem::is_directory(root_path)) fail("data source path is not a directory: " + root);
        for (const auto& entry : std::filesystem::directory_iterator(root_path)) {
            if (!entry.is_directory()) continue;
            const auto version_dir = entry.path() / options.ref_version;
            if (!std::filesystem::is_directory(version_dir)) continue;
            for (const auto& file : std::filesystem::directory_iterator(version_dir)) {
                if (file.path().extension() != ".config") continue;
                const auto config = load_config(file.path().string());
                DataSource source;
                source.directory = version_dir.string();
                source.name = config_value(config, "name", true, file.path().string());
                source.version = config_value(config, "version", true, file.path().string());
                source.type = config_value(config, "type", true, file.path().string());
                source.src_file = config_value(config, "src_file", true, file.path().string());
                config_value(config, "origin_location", true, file.path().string());
                config_value(config, "preprocessing_script", true, file.path().string());
                if (source.type == "gencode") {
                    source.ncbi_build = config_value(config, "ncbi_build_version", false,
                                                     file.path().string());
                    source.gencode_fasta = config_value(config, "gencode_fasta_path", true,
                                                        file.path().string());
                    found_gencode = true;
                } else if (source.type == "simpleXSV") {
                    source.xsv_key = config_value(config, "xsv_key", true, file.path().string());
                    source.xsv_key_column = static_cast<std::size_t>(std::stoul(
                        config_value(config, "xsv_key_column", true, file.path().string())));
                    std::string delimiter = config_value(config, "xsv_delimiter", true,
                                                         file.path().string());
                    if (delimiter == "\\t") delimiter = "\t";
                    source.xsv_delimiter = delimiter;
                    config_value(config, "xsv_permissive_cols", true, file.path().string());
                } else {
                    fail("unsupported data source type (supported: gencode, simpleXSV): " +
                         source.type);
                }
                sources.push_back(std::move(source));
            }
        }
    }
    if (!found_gencode)
        fail("No Gencode data source was specified.  Funcotator requires a gencode "
             "data source for transcript-based annotation.");
    // datasourceComparator: gencode first, then by name (stable).
    std::stable_sort(sources.begin(), sources.end(),
        [](const DataSource& left, const DataSource& right) {
            const bool left_gencode = left.type == "gencode";
            const bool right_gencode = right.type == "gencode";
            if (left_gencode != right_gencode) return left_gencode;
            return left.name < right.name;
        });
    return sources;
}

// ------------------------------------------------------------------ GTF ----

struct Exon {
    std::int64_t start = 0;   // 1-based inclusive
    std::int64_t end = 0;
    int number = 0;
};

struct Transcript {
    std::string id;
    std::string gene_name;
    std::string gene_id;
    std::string strand;
    std::int64_t start = 0;
    std::int64_t end = 0;
    std::vector<Exon> exons;
    std::vector<std::pair<std::int64_t, std::int64_t>> coding;
    bool appris_principal = false;
    bool basic = false;
};

std::string gtf_attribute(const std::string& attributes, const std::string& key) {
    const auto marker = attributes.find(key + " \"");
    if (marker == std::string::npos) return {};
    const auto begin = marker + key.size() + 2;
    const auto end = attributes.find('"', begin);
    return attributes.substr(begin, end - begin);
}

std::vector<Transcript> load_gtf(const std::string& path) {
    std::map<std::string, Transcript> transcripts;
    std::ifstream input(path);
    if (!input) fail("cannot read GTF: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split(line, '\t');
        if (fields.size() < 9) continue;
        const std::string& feature = fields[2];
        const std::string& attributes = fields[8];
        const std::string transcript_id = gtf_attribute(attributes, "transcript_id");
        if (transcript_id.empty()) continue;
        auto& transcript = transcripts[transcript_id];
        transcript.id = transcript_id;
        transcript.gene_name = gtf_attribute(attributes, "gene_name");
        transcript.gene_id = gtf_attribute(attributes, "gene_id");
        transcript.strand = fields[6];
        const auto start = std::strtoll(fields[3].c_str(), nullptr, 10);
        const auto end = std::strtoll(fields[4].c_str(), nullptr, 10);
        if (transcript.start == 0 || start < transcript.start) transcript.start = start;
        if (end > transcript.end) transcript.end = end;
        if (attributes.find("appris_principal") != std::string::npos)
            transcript.appris_principal = true;
        if (attributes.find("tag \"basic\"") != std::string::npos) transcript.basic = true;
        if (feature == "exon") {
            Exon exon;
            exon.start = start;
            exon.end = end;
            const std::string number = gtf_attribute(attributes, "exon_number");
            exon.number = number.empty() ? 0 : std::stoi(number);
            transcript.exons.push_back(exon);
        } else if (feature == "CDS") {
            transcript.coding.emplace_back(start, end);
        }
    }
    std::vector<Transcript> ordered;
    ordered.reserve(transcripts.size());
    for (auto& [id, transcript] : transcripts) {
        std::sort(transcript.exons.begin(), transcript.exons.end(),
            [](const Exon& left, const Exon& right) {
                return left.start < right.start;
            });
        std::sort(transcript.coding.begin(), transcript.coding.end());
        ordered.push_back(std::move(transcript));
    }
    // GTF order is preserved for ties by the map iteration only loosely;
    // keep transcript-id order for determinism, then let selection rules
    // decide preference.
    std::sort(ordered.begin(), ordered.end(),
        [](const Transcript& left, const Transcript& right) {
            return left.id < right.id;
        });
    return ordered;
}

// ------------------------------------------------------------ XSV sources --

void load_xsv(DataSource& source) {
    const auto path = std::filesystem::path(source.directory) / source.src_file;
    std::ifstream input(path);
    if (!input) fail("cannot read XSV source: " + path.string());
    std::string line;
    bool first = true;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto columns = split(line, source.xsv_delimiter.front());
        if (first) {
            source.xsv_header = columns;
            first = false;
            continue;
        }
        if (source.xsv_key_column >= columns.size()) continue;
        source.xsv_rows.emplace_back(columns[source.xsv_key_column], line);
    }
}

// ------------------------------------------------------- funcotation core --

std::string sanitize(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char character : value) {
        if (character == ',' || character == ';' || character == '=' ||
            character == '\t' || character == ' ' || character == '|' ||
            character == '\n' || character == '#') {
            out << "_%";
            const char* digits = "0123456789ABCDEF";
            out << digits[character >> 4] << digits[character & 0xF] << '_';
        } else {
            out << static_cast<char>(character);
        }
    }
    return out.str();
}

// Standard genetic code translation of a 3-base codon (uppercase).
std::string translate_codon(const std::string& codon) {
    static const char* bases = "TCAG";
    static const char* table =
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    if (codon.size() != 3) return "?";
    int index = 0;
    for (int position = 0; position < 3; ++position) {
        const auto found = std::char_traits<char>::find(bases, 4, codon[position]);
        if (found == nullptr) return "?";
        index = index * 4 + static_cast<int>(found - bases);
    }
    const char amino = table[index];
    return amino == '*' ? "*" : std::string(1, amino);
}

std::string reverse_complement(const std::string& sequence) {
    std::string result;
    result.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        switch (*iterator) {
            case 'A': result.push_back('T'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break;
        }
    }
    return result;
}

std::string to_lower(const std::string& text) {
    std::string result = text;
    std::transform(result.begin(), result.end(), result.begin(),
        [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return result;
}

struct FuncotationBlock {
    std::vector<std::pair<std::string, std::string>> fields;
};

std::string render_block(const FuncotationBlock& block) {
    // FuncotatorUtils.renderSanitizedFuncotationForVcf joins the *sanitized*
    // field values with raw '|' separators.
    std::string out;
    for (std::size_t index = 0; index < block.fields.size(); ++index) {
        if (index != 0) out += '|';
        out += sanitize(block.fields[index].second);
    }
    return out;
}

std::string variant_type(const std::string& ref, const std::string& alt) {
    if (ref.size() == 1 && alt.size() == 1) return "SNP";
    if (ref.size() == alt.size()) {
        if (ref.size() == 2) return "DNP";
        if (ref.size() == 3) return "TNP";
        return "ONP";
    }
    return ref.size() < alt.size() ? "INS" : "DEL";
}

std::string genome_change(const std::string& contig, std::int64_t start,
                          const std::string& ref, const std::string& alt) {
    std::ostringstream out;
    out << "g." << contig << ':' << start;
    if (ref.size() == 1 && alt.size() == 1) {
        out << ref << '>' << alt;
    } else if (ref.size() < alt.size()) {
        out << '_' << (start + static_cast<std::int64_t>(ref.size()) - 1)
            << "ins" << alt.substr(ref.size());
    } else if (ref.size() > alt.size()) {
        out << '_' << (start + static_cast<std::int64_t>(ref.size()) - 1)
            << "del" << ref.substr(alt.size());
    } else {
        out << '_' << (start + static_cast<std::int64_t>(ref.size()) - 1)
            << "ins" << alt.substr(1) << "del" << ref.substr(1);
    }
    return out.str();
}

// Java Double.toString semantics: shortest round-trip decimal.
std::string double_text(const double value) {
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    return std::string(buffer, result.ptr);
}

double gc_content(const Fasta& reference, const std::string& contig,
                  std::int64_t start, const std::string& ref, const std::string& alt,
                  int window = 200) {
    const int leading = (ref.size() == alt.size()) ? window : window - 1;
    const std::int64_t fetch_start = start - leading;
    const std::int64_t fetch_end = start + static_cast<std::int64_t>(ref.size()) - 1 + window;
    // getBases clamps at the contig boundary instead of failing.
    const std::int64_t clamped_start = std::max<std::int64_t>(1, fetch_start);
    std::string bases = reference.fetch(contig, clamped_start, fetch_end);
    std::size_t gc = 0;
    for (const char base : bases)
        if (base == 'G' || base == 'C') ++gc;
    return bases.empty() ? 0.0
        : static_cast<double>(gc) / static_cast<double>(bases.size());
}

std::string reference_context(const Fasta& reference, const std::string& contig,
                              std::int64_t start, const std::string& ref,
                              const std::string& alt, std::string& context_ref,
                              std::string& context_alt) {
    // GATK's GencodeFuncotation referenceContext is 21 bases for a SNP:
    // 10 leading + ref + 10 trailing; for indels the window is built around
    // the padding base with the alt replacing the event body.
    const std::int64_t context_start = start - 10;
    std::string leading = context_start >= 1
        ? reference.fetch(contig, context_start, start - 1) : std::string();
    while (leading.size() < 10) leading.insert(leading.begin(), 'N');
    const std::string trailing = reference.fetch(
        contig, start + static_cast<std::int64_t>(ref.size()),
        start + static_cast<std::int64_t>(ref.size()) + 9);
    context_ref = leading + ref + trailing;
    context_alt = leading + alt + trailing;
    return context_ref;
}

struct TranscriptContext {
    bool overlaps = false;
    bool in_coding = false;
    bool in_exon = false;
    bool on_forward = true;
    int exon_number = 0;
    std::int64_t transcript_pos = 0;   // 1-based cDNA position of the ref start
    std::int64_t coding_pos = 0;       // 1-based CDS position of the ref start
    std::string exon_sequence_before;  // transcript bases before the variant
};

bool in_interval(const std::pair<std::int64_t, std::int64_t>& span,
                 std::int64_t position) {
    return position >= span.first && position <= span.second;
}

TranscriptContext classify_position(const Transcript& transcript,
                                    std::int64_t start, std::int64_t end) {
    TranscriptContext context;
    context.on_forward = transcript.strand != "-";
    for (const auto& exon : transcript.exons) {
        if (end >= exon.start && start <= exon.end) {
            context.overlaps = true;
            context.in_exon = true;
            context.exon_number = exon.number;
            break;
        }
    }
    for (const auto& span : transcript.coding) {
        if (end >= span.first && start <= span.second) {
            context.overlaps = true;
            context.in_coding = true;
            break;
        }
    }
    if (!context.overlaps && start <= transcript.end && end >= transcript.start)
        context.overlaps = true;   // intronic
    return context;
}

std::string trim_alt_for_change(const std::string& ref, const std::string& alt,
                                std::int64_t& change_position) {
    // Minimum-representation style trimming for change strings.
    std::string local_ref = ref;
    std::string local_alt = alt;
    std::size_t prefix = 0;
    while (local_ref.size() - prefix > 1 && local_alt.size() - prefix > 1 &&
           local_ref[prefix] == local_alt[prefix])
        ++prefix;
    change_position += static_cast<std::int64_t>(prefix);
    local_ref = local_ref.substr(prefix);
    local_alt = local_alt.substr(prefix);
    while (local_ref.size() > 1 && local_alt.size() > 1 &&
           local_ref.back() == local_alt.back()) {
        local_ref.pop_back();
        local_alt.pop_back();
    }
    return local_ref + ">" + local_alt;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);

        // Refuse already-annotated inputs (Funcotator.onTraversalStart).
        {
            htsFile* probe = bcf_open(options.variants.c_str(), "r");
            if (probe == nullptr) fail("cannot open input VCF: " + options.variants);
            bcf_hdr_t* probe_header = bcf_hdr_read(probe);
            if (probe_header == nullptr) {
                hts_close(probe);
                fail("cannot read input VCF header");
            }
            for (int line = 0; line < probe_header->nhrec; ++line) {
                const bcf_hrec_t* record = probe_header->hrec[line];
                for (int item = 0; record != nullptr && item < record->nkeys; ++item) {
                    if (std::string(record->keys[item]) == "ID" &&
                        record->vals[item] != nullptr &&
                        std::string(record->vals[item]) == "Funcotator Version") {
                        bcf_hdr_destroy(probe_header);
                        hts_close(probe);
                        fail("The input VCF is already annotated by Funcotator.");
                    }
                }
            }
            bcf_hdr_destroy(probe_header);
            hts_close(probe);
        }

        auto data_sources = discover_data_sources(options);
        for (auto& source : data_sources)
            if (source.type == "simpleXSV") load_xsv(source);

        std::vector<Transcript> transcripts;
        for (const auto& source : data_sources)
            if (source.type == "gencode")
                for (auto&& transcript : load_gtf(
                        (std::filesystem::path(source.directory) / source.src_file).string()))
                    transcripts.push_back(std::move(transcript));

        Fasta reference(options.reference);
        std::optional<Fasta> transcripts_fasta;
        for (const auto& source : data_sources)
            if (source.type == "gencode" && !source.gencode_fasta.empty()) {
                transcripts_fasta.emplace(
                    (std::filesystem::path(source.directory) / source.gencode_fasta).string());
                break;
            }

        // Header: input metadata + Funcotator lines.
        htsFile* input = bcf_open(options.variants.c_str(), "r");
        if (input == nullptr) fail("cannot open input VCF: " + options.variants);
        bcf_hdr_t* header = bcf_hdr_read(input);
        if (header == nullptr) {
            hts_close(input);
            fail("cannot read input VCF header");
        }
        kstring_t raw_header = {0, 0, nullptr};
        if (bcf_hdr_format(header, 0, &raw_header) < 0 || raw_header.s == nullptr) {
            free(raw_header.s);
            fail("cannot format input header");
        }
        const std::string header_text(raw_header.s, raw_header.l);
        free(raw_header.s);

        // Field list: gencode fields (data-source order) then XSV fields then
        // default/override keys (Funcotator field ordering: data sources in
        // comparator order, gencode first).
        std::vector<std::string> field_names;
        const auto gencode_field = [](const std::string& name, const std::string& version,
                                      const std::string& field) {
            return name + "_" + version + "_" + field;
        };
        std::string gencode_name = "Gencode";
        std::string gencode_version = "Mini";
        for (const auto& source : data_sources) {
            if (source.type != "gencode") continue;
            gencode_name = source.name;
            gencode_version = source.version;
            for (const std::string field : {"hugoSymbol", "ncbiBuild", "chromosome",
                     "start", "end", "variantClassification",
                     "secondaryVariantClassification", "variantType", "refAllele",
                     "tumorSeqAllele1", "tumorSeqAllele2", "genomeChange",
                     "annotationTranscript", "transcriptStrand", "transcriptExon",
                     "transcriptPos", "cDnaChange", "codonChange", "proteinChange",
                     "gcContent", "referenceContext", "otherTranscripts"})
                field_names.push_back(gencode_field(source.name, source.version, field));
        }
        for (const auto& source : data_sources) {
            if (source.type != "simpleXSV") continue;
            for (std::size_t column = 0; column < source.xsv_header.size(); ++column) {
                if (column == source.xsv_key_column) continue;
                std::string name = source.xsv_header[column];
                std::replace(name.begin(), name.end(), ' ', '_');
                field_names.push_back(source.name + "_" + name);
            }
        }
        std::map<std::string, std::string> defaults;
        std::map<std::string, std::string> overrides;
        for (const auto& entry : options.annotation_defaults) {
            const auto separator = entry.find(':');
            if (separator == std::string::npos) fail("invalid --annotation-default: " + entry);
            defaults[entry.substr(0, separator)] = entry.substr(separator + 1);
        }
        for (const auto& entry : options.annotation_overrides) {
            const auto separator = entry.find(':');
            if (separator == std::string::npos) fail("invalid --annotation-override: " + entry);
            overrides[entry.substr(0, separator)] = entry.substr(separator + 1);
        }
        for (const auto& [key, value] : defaults) {
            (void)value;
            if (std::find(field_names.begin(), field_names.end(), key) == field_names.end())
                field_names.push_back(key);
        }
        for (const auto& [key, value] : overrides) {
            (void)value;
            if (std::find(field_names.begin(), field_names.end(), key) == field_names.end())
                field_names.push_back(key);
        }

        // Funcotator Version line: <tool> | <name> <version> <mode> | ...
        std::string data_source_info;
        for (const auto& source : data_sources) {
            if (!data_source_info.empty()) data_source_info += " | ";
            data_source_info += source.name + " " + source.version;
            if (source.type == "gencode")
                data_source_info += " " + options.transcript_selection_mode;
        }

        std::string funcotation_description =
            "Functional annotation from the Funcotator tool.  Funcotation fields are: ";
        for (std::size_t index = 0; index < field_names.size(); ++index) {
            if (index != 0) funcotation_description += '|';
            funcotation_description += field_names[index];
        }

        // htsjdk VCFEncoder ordering: fileformat first, then FILTER, FORMAT
        // and INFO in sorted ID order, then contig lines, then remaining
        // metadata in input order.  The implicit PASS definition is dropped
        // and Funcotator adds its own ##source line.
        std::map<std::string, std::string> header_format_lines;
        std::map<std::string, std::string> header_info_lines;
        std::vector<std::string> header_filter_lines;
        std::vector<std::string> header_contig_lines;
        std::vector<std::string> header_other_lines;
        {
            std::size_t start_line = 0;
            while (start_line < header_text.size()) {
                const auto newline = header_text.find('\n', start_line);
                std::string line = header_text.substr(start_line,
                    newline == std::string::npos ? std::string::npos : newline - start_line);
                start_line = newline == std::string::npos ? header_text.size() : newline + 1;
                if (line.empty() || line.rfind("#CHROM", 0) == 0 ||
                    line.rfind("##fileformat=", 0) == 0)
                    continue;
                const auto id_marker = line.find("=<ID=");
                std::string id;
                if (id_marker != std::string::npos) {
                    const auto id_end = line.find(',', id_marker + 5);
                    id = line.substr(id_marker + 5, id_end == std::string::npos
                        ? std::string::npos : id_end - (id_marker + 5));
                }
                if (line.rfind("##FILTER=", 0) == 0) {
                    if (id != "PASS") header_filter_lines.push_back(line);
                } else if (line.rfind("##FORMAT=", 0) == 0) {
                    header_format_lines[id] = line;
                } else if (line.rfind("##INFO=", 0) == 0) {
                    header_info_lines[id] = line;
                } else if (line.rfind("##contig=", 0) == 0) {
                    header_contig_lines.push_back(line);
                } else {
                    header_other_lines.push_back(line);
                }
            }
        }
        header_info_lines["FUNCOTATION"] =
            "##INFO=<ID=FUNCOTATION,Number=A,Type=String,Description=\"" +
            funcotation_description + "\">";
        std::ostringstream header_out;
        header_out << "##fileformat=VCFv4.2\n";
        std::sort(header_filter_lines.begin(), header_filter_lines.end());
        for (const auto& line : header_filter_lines) header_out << line << '\n';
        for (const auto& [id, line] : header_format_lines) header_out << line << '\n';
        for (const auto& [id, line] : header_info_lines) header_out << line << '\n';
        for (const auto& line : header_contig_lines) header_out << line << '\n';
        header_out << "##Funcotator Version=fastgatk-native | " << data_source_info << '\n';
        header_out << "##source=Funcotator\n";
        for (const auto& line : header_other_lines) header_out << line << '\n';
        header_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        if (bcf_hdr_nsamples(header) > 0) {
            header_out << "\tFORMAT";
            for (int index = 0; index < bcf_hdr_nsamples(header); ++index)
                header_out << '\t' << bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
        }
        header_out << '\n';

        if (options.output.size() > 3 &&
            options.output.substr(options.output.size() - 3) == ".gz")
            fail("compressed output is not supported; write an uncompressed .vcf");
        std::ofstream output(options.output);
        if (!output) fail("cannot open output VCF: " + options.output);
        {
            const std::string rendered = header_out.str();
            output.write(rendered.data(), static_cast<std::streamsize>(rendered.size()));
        }

        bcf1_t* record = bcf_init();
        while (bcf_read(input, header, record) == 0) {
            if (bcf_unpack(record, BCF_UN_ALL) != 0) continue;
            const char* contig_name = bcf_hdr_id2name(header, record->rid);
            const std::string contig = contig_name == nullptr ? "" : contig_name;
            const std::int64_t start = record->pos + 1;
            const std::string ref = record->d.allele[0] == nullptr
                ? "" : record->d.allele[0];

            // Per-ALT funcotation blocks.
            std::vector<std::string> per_alt;
            for (int alt_index = 1; alt_index < record->n_allele; ++alt_index) {
                const std::string alt = record->d.allele[alt_index] == nullptr
                    ? "" : record->d.allele[alt_index];
                if (alt == "*") continue;

                // Transcript selection: ALL keeps every transcript at the
                // locus; CANONICAL keeps the appris_principal one, else the
                // longest, else the first (GencodeFuncotationFactory's
                // canonical preference on the fixture surface).
                std::vector<const Transcript*> candidates;
                const std::int64_t end = start + static_cast<std::int64_t>(ref.size()) - 1;
                for (const auto& transcript : transcripts) {
                    const auto context = classify_position(transcript, start, end);
                    if (context.overlaps) candidates.push_back(&transcript);
                }
                std::stable_sort(candidates.begin(), candidates.end(),
                    [](const Transcript* left, const Transcript* right) {
                        if (left->appris_principal != right->appris_principal)
                            return left->appris_principal;
                        return (left->end - left->start) > (right->end - right->start);
                    });
                std::vector<const Transcript*> selected;
                if (options.transcript_selection_mode == "ALL") {
                    selected = candidates;
                } else if (!candidates.empty()) {
                    selected.push_back(candidates.front());
                }

                std::string context_ref;
                std::string context_alt;
                const std::string context = reference_context(
                    reference, contig, start, ref, alt, context_ref, context_alt);
                const double gc = gc_content(reference, contig, start, ref, alt);

                std::vector<FuncotationBlock> blocks;
                if (selected.empty()) {
                    FuncotationBlock block;
                    std::string trimmed_ref;
                    std::string trimmed_alt;
                    {
                        std::int64_t position = start;
                        const std::string change = trim_alt_for_change(ref, alt, position);
                        const auto separator = change.find('>');
                        trimmed_ref = change.substr(0, separator);
                        trimmed_alt = change.substr(separator + 1);
                    }
                    block.fields = {
                        {"hugoSymbol", "Unknown"},
                        {"ncbiBuild", data_sources.empty() ? "" : data_sources.front().ncbi_build},
                        {"chromosome", contig},
                        {"start", std::to_string(start)},
                        {"end", std::to_string(end)},
                        {"variantClassification", "IGR"},
                        {"secondaryVariantClassification", ""},
                        {"variantType", variant_type(ref, alt)},
                        {"refAllele", ref},
                        {"tumorSeqAllele1", ref},
                        {"tumorSeqAllele2", alt},
                        {"genomeChange", genome_change(contig, start, ref, alt)},
                        {"annotationTranscript", "no_transcript"},
                        {"transcriptStrand", ""},
                        {"transcriptExon", ""},
                        {"transcriptPos", ""},
                        {"cDnaChange", ""},
                        {"codonChange", ""},
                        {"proteinChange", ""},
                        {"gcContent", double_text(gc)},
                        {"referenceContext", context},
                        {"otherTranscripts", ""},
                    };
                    blocks.push_back(std::move(block));
                } else {
                    for (const Transcript* transcript : selected) {
                        const auto context_info = classify_position(*transcript, start, end);
                        FuncotationBlock block;
                        std::string classification;
                        std::string c_dna_change;
                        std::string codon_change;
                        std::string protein_change;
                        std::string exon_number;
                        std::string transcript_pos;
                        std::string strand(1, transcript->strand == "-" ? '-' : '+');
                        std::int64_t c_dna_position = 0;

                        // Transcript coordinates: sum exon lengths.
                        const bool forward = transcript->strand != "-";
                        std::int64_t position_in_transcript = 0;
                        std::int64_t coding_position = 0;
                        std::string transcript_bases;
                        for (const auto& exon : transcript->exons) {
                            const std::int64_t exon_length = exon.end - exon.start + 1;
                            if (start >= exon.start && start <= exon.end) {
                                const std::int64_t offset_in_exon = start - exon.start;
                                position_in_transcript += forward
                                    ? offset_in_exon + 1
                                    : exon_length - offset_in_exon;
                                exon_number = std::to_string(exon.number);
                            } else if (start > exon.end) {
                                position_in_transcript += exon_length;
                            }
                            if (!transcripts_fasta.has_value()) continue;
                        }
                        for (const auto& span : transcript->coding) {
                            const std::int64_t span_length = span.second - span.first + 1;
                            if (start >= span.first && start <= span.second) {
                                const std::int64_t offset_in_span = start - span.first;
                                coding_position += forward
                                    ? offset_in_span + 1
                                    : span_length - offset_in_span;
                            } else if (start > span.second) {
                                coding_position += span_length;
                            }
                        }
                        c_dna_position = position_in_transcript;
                        transcript_pos = std::to_string(position_in_transcript);
                        if (context_info.in_coding && coding_position > 0) {
                            // Codon and protein changes from the transcript
                            // FASTA (GencodeFuncotation's SequenceComparison).
                            std::string tx_sequence;
                            if (transcripts_fasta.has_value())
                                tx_sequence = transcripts_fasta->fetch(
                                    transcript->id, 1,
                                    transcript->end - transcript->start + 1);
                            if (!tx_sequence.empty()) {
                                if (!forward) tx_sequence = reverse_complement(tx_sequence);
                                const std::int64_t codon_start =
                                    ((coding_position - 1) / 3) * 3;
                                if (codon_start + 3 <=
                                    static_cast<std::int64_t>(tx_sequence.size())) {
                                    const std::string ref_codon =
                                        tx_sequence.substr(static_cast<std::size_t>(codon_start), 3);
                                    std::string alt_codon = ref_codon;
                                    const std::size_t base_in_codon =
                                        static_cast<std::size_t>((coding_position - 1) % 3);
                                    std::string event_ref = ref;
                                    std::string event_alt = alt;
                                    if (!forward) {
                                        event_ref = reverse_complement(ref);
                                        event_alt = reverse_complement(alt);
                                    }
                                    alt_codon.replace(base_in_codon, 1, event_alt.substr(0, 1));
                                    const std::string ref_amino = translate_codon(ref_codon);
                                    const std::string alt_amino = translate_codon(alt_codon);
                                    const std::int64_t codon_number = codon_start / 3 + 1;
                                    std::ostringstream codon_text;
                                    std::string ref_styled = to_lower(ref_codon);
                                    std::string alt_styled = to_lower(alt_codon);
                                    ref_styled[base_in_codon] =
                                        static_cast<char>(std::toupper(ref_styled[base_in_codon]));
                                    alt_styled[base_in_codon] =
                                        static_cast<char>(std::toupper(alt_styled[base_in_codon]));
                                    codon_text << "c.(" << (codon_start + 1) << '-'
                                               << (codon_start + 3) << ')' << ref_styled
                                               << '>' << alt_styled;
                                    codon_change = codon_text.str();
                                    std::ostringstream protein_text;
                                    if (alt_amino == "*") {
                                        protein_text << "p." << ref_amino << codon_number << "*";
                                    } else {
                                        protein_text << "p." << ref_amino << codon_number
                                                     << alt_amino;
                                    }
                                    protein_change = protein_text.str();
                                    if (ref_amino == alt_amino)
                                        classification = "SYNONYMOUS";
                                    else if (alt_amino == "*")
                                        classification = "NONSENSE";
                                    else
                                        classification = "MISSENSE";
                                }
                            }
                            std::int64_t change_position = start;
                            const std::string change =
                                trim_alt_for_change(ref, alt, change_position);
                            std::ostringstream cdna_text;
                            cdna_text << "c." << c_dna_position << change;
                            c_dna_change = cdna_text.str();
                        } else if (context_info.in_exon) {
                            // UTR side is decided below against the CDS span.
                            std::int64_t first_coding = 0;
                            std::int64_t last_coding = 0;
                            for (const auto& span : transcript->coding) {
                                if (first_coding == 0 || span.first < first_coding)
                                    first_coding = span.first;
                                if (span.second > last_coding) last_coding = span.second;
                            }
                            const bool before_coding = first_coding != 0 && end < first_coding;
                            const bool after_coding = last_coding != 0 && start > last_coding;
                            if (forward)
                                classification = before_coding ? "FIVE_PRIME_UTR"
                                    : (after_coding ? "THREE_PRIME_UTR" : "INTRON");
                            else
                                classification = before_coding ? "THREE_PRIME_UTR"
                                    : (after_coding ? "FIVE_PRIME_UTR" : "INTRON");
                        } else {
                            classification = "INTRON";
                        }
                        {
                            std::int64_t change_position = start;
                            const std::string change =
                                trim_alt_for_change(ref, alt, change_position);
                            std::ostringstream cdna_text;
                            cdna_text << "c." << c_dna_position << change;
                            c_dna_change = cdna_text.str();
                        }

                        // otherTranscripts: <gene>_<tx>_<classification>_<protein>
                        // for the non-selected transcripts at this locus.
                        std::string other;
                        for (const Transcript* candidate : candidates) {
                            if (candidate == transcript) continue;
                            if (!other.empty()) other += ",";
                            other += candidate->gene_name + "_" + candidate->id + "_" +
                                classification + "_" + protein_change;
                        }

                        block.fields = {
                            {"hugoSymbol", transcript->gene_name},
                            {"ncbiBuild", data_sources.empty() ? "" : data_sources.front().ncbi_build},
                            {"chromosome", contig},
                            {"start", std::to_string(start)},
                            {"end", std::to_string(end)},
                            {"variantClassification", classification},
                            {"secondaryVariantClassification", ""},
                            {"variantType", variant_type(ref, alt)},
                            {"refAllele", ref},
                            {"tumorSeqAllele1", ref},
                            {"tumorSeqAllele2", alt},
                            {"genomeChange", genome_change(contig, start, ref, alt)},
                            {"annotationTranscript", transcript->id},
                            {"transcriptStrand", strand},
                            {"transcriptExon", exon_number},
                            {"transcriptPos", transcript_pos},
                            {"cDnaChange", c_dna_change},
                            {"codonChange", codon_change},
                            {"proteinChange", protein_change},
                            {"gcContent", double_text(gc)},
                            {"referenceContext", context},
                            {"otherTranscripts", other},
                        };
                        blocks.push_back(std::move(block));
                    }
                }

                // XSV + default/override fields appended to every block.
                for (auto& block : blocks) {
                    for (const auto& source : data_sources) {
                        if (source.type != "simpleXSV") continue;
                        std::string key_value;
                        for (const auto& field : block.fields) {
                            if (field.first == "hugoSymbol") {
                                key_value = field.second;
                                break;
                            }
                        }
                        std::string row;
                        bool found = false;
                        for (const auto& [key, raw] : source.xsv_rows) {
                            if (key == key_value && !key_value.empty()) {
                                row = raw;
                                found = true;
                                break;
                            }
                        }
                        const auto columns = found
                            ? split(row, source.xsv_delimiter.front())
                            : std::vector<std::string>{};
                        for (std::size_t column = 0; column < source.xsv_header.size(); ++column) {
                            if (column == source.xsv_key_column) continue;
                            std::string name = source.xsv_header[column];
                            std::replace(name.begin(), name.end(), ' ', '_');
                            block.fields.emplace_back(
                                source.name + "_" + name,
                                found && column < columns.size() ? columns[column] : "");
                        }
                    }
                    for (const auto& [key, value] : defaults) {
                        (void)value;
                        bool present = false;
                        for (const auto& field : block.fields)
                            if (field.first == key) { present = true; break; }
                        if (!present) block.fields.emplace_back(key, defaults.at(key));
                    }
                    for (const auto& [key, value] : overrides) {
                        bool present = false;
                        for (auto& field : block.fields)
                            if (field.first == key) {
                                field.second = value;
                                present = true;
                                break;
                            }
                        if (!present) block.fields.emplace_back(key, value);
                    }
                }

                std::string alt_text;
                for (std::size_t block_index = 0; block_index < blocks.size(); ++block_index) {
                    if (block_index != 0) alt_text += '#';
                    alt_text += '[' + render_block(blocks[block_index]) + ']';
                }
                per_alt.push_back(alt_text);
            }

            // Serialize the record with FUNCOTATION merged into INFO.
            kstring_t line = {0, 0, nullptr};
            if (vcf_format(header, record, &line) < 0 || line.s == nullptr) {
                free(line.s);
                fail("cannot format VCF record");
            }
            std::string text(line.s, line.l);
            free(line.s);
            while (!text.empty() && (text.back() == '\n' || text.back() == '\r'))
                text.pop_back();
            std::vector<std::string> fields;
            std::size_t field_start = 0;
            for (std::size_t index = 0; index <= text.size(); ++index) {
                if (index == text.size() || text[index] == '\t') {
                    fields.push_back(text.substr(field_start, index - field_start));
                    field_start = index + 1;
                }
            }
            if (fields.size() < 8) fail("malformed VCF record");
            std::string info = fields[7];
            if (info == "." ) info.clear();
            if (!info.empty()) info += ';';
            std::string alt_values;
            for (std::size_t index = 0; index < per_alt.size(); ++index) {
                if (index != 0) alt_values += ',';
                alt_values += per_alt[index];
            }
            info += "FUNCOTATION=" + alt_values;
            fields[7] = info;
            std::ostringstream out;
            for (std::size_t index = 0; index < fields.size(); ++index) {
                if (index != 0) out << '\t';
                out << fields[index];
            }
            out << '\n';
            const std::string rendered = out.str();
            output.write(rendered.data(), static_cast<std::streamsize>(rendered.size()));
        }
        bcf_destroy(record);
        bcf_hdr_destroy(header);
        hts_close(input);
        output.close();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
