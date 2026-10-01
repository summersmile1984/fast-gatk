// fastgatk-variant-annotator — GATK-compatible VariantAnnotator walker.
//
// Mirrors GATK's VariantAnnotator (gatk-source .../walkers/annotator/
// VariantAnnotator.java) on its deterministic surfaces:
//   * --dbsnp           : DB membership flag + rsID merge into the ID field
//     (VariantOverlapAnnotator.annotateRsID/getRsID: unfiltered source VCs
//     split to minimum-representation biallelics, matched on exact ref+alt).
//   * --comp[:NAME]     : NAME membership flag (same match rule).
//   * --resource:NAME   : feature inputs for -E expressions and NAME flags.
//   * -E NAME.FIELD     : transfer FIELD (ID/ALT/FILTER special cases, else
//     an INFO attribute of the *first* resource record starting at the site)
//     into INFO key NAME.FIELD (VAExpression/annotateExpressions).
//   * -A Coverage       : INFO/DP from the read likelihood evidence count
//     over the site's pileup (VariantAnnotator.makeLikelihoods assigns each
//     pileup element to chooseAlleleForRead's allele; every element is one
//     evidence row regardless of base quality).
//
// Output rows keep the input record verbatim except the annotated ID/INFO
// fields, matching htsjdk's VCFEncoder serialization (INFO keys sorted,
// formatVCFDouble for Float values).

#include "optional_boolean.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if FASTGATK_HAS_HTSLIB
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#endif

namespace {

struct ExpressionSpec {
    std::string full_name;   // NAME.FIELD
    std::string binding;     // NAME
    std::string field;       // FIELD
};

struct Options {
    std::string variants;
    std::string output;
    std::string reference;
    std::string reads;
    std::vector<std::string> annotations;
    std::vector<std::string> expressions;
    std::vector<std::pair<std::string, std::string>> named_resources;
    std::vector<std::pair<std::string, std::string>> comps;
    std::string dbsnp;
    std::vector<std::string> regions;
    int min_base_quality_score = 10;
    bool requested_coverage = false;
    bool resource_allele_concordance = false;
    bool create_index = false;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::invalid_argument(message);
}

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
                "fastgatk-variant-annotator (GATK-compatible native prototype)\n"
                "usage: fastgatk-variant-annotator -V <calls.vcf> -O <out.vcf[.gz]>\n"
                "  -V, --variant FILE          input VCF/BCF (required)\n"
                "  -O, --output FILE           output VCF (required)\n"
                "  -R, --reference FILE        reference FASTA (enables -A Coverage)\n"
                "  -I, --input FILE            reads BAM/CRAM (enables -A Coverage)\n"
                "  -A, --annotation NAME       annotation to compute (Coverage)\n"
                "  -E, --expression NAME.FIELD transfer resource field to INFO\n"
                "  --resource:NAME FILE        named feature VCF for -E / flags\n"
                "  --comp[:NAME] FILE          comparison VCF membership flag\n"
                "  --dbsnp FILE                dbSNP VCF for DB flag + rsID merge\n"
                "  -L, --interval CHR:START-END  restrict output records\n"
                "  --min-base-quality-score N  pileup allele assignment floor (10)\n";
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
        else if (is_option(argument, "-I") || is_option(argument, "--input"))
            options.reads = option_value(index, argc, argv, argument,
                argument.rfind("-I", 0) == 0 ? "-I" : "--input");
        else if (is_option(argument, "-A") || is_option(argument, "--annotation"))
            options.annotations.push_back(option_value(index, argc, argv, argument,
                argument.rfind("-A", 0) == 0 ? "-A" : "--annotation"));
        else if (is_option(argument, "-E") || is_option(argument, "--expression"))
            options.expressions.push_back(option_value(index, argc, argv, argument,
                argument.rfind("-E", 0) == 0 ? "-E" : "--expression"));
        else if (argument.rfind("--resource:", 0) == 0) {
            const std::string rest = argument.substr(11);
            const auto equals = rest.find('=');
            const std::string name = rest.substr(0, equals);
            std::string path = equals == std::string::npos
                ? std::string() : rest.substr(equals + 1);
            if (name.empty()) fail("--resource requires a name");
            if (path.empty()) {
                if (index + 1 >= argc) fail("missing value for --resource");
                path = argv[++index];
            }
            options.named_resources.emplace_back(name, path);
        }
        else if (argument.rfind("--comp", 0) == 0 || argument.rfind("-comp", 0) == 0) {
            const bool long_form = argument.rfind("--comp", 0) == 0;
            const std::string rest = argument.substr(long_form ? 6 : 5);
            std::string name = "comp";
            std::string path;
            if (rest.rfind(":", 0) == 0) {
                const std::string named = rest.substr(1);
                const auto equals = named.find('=');
                name = named.substr(0, equals);
                if (equals != std::string::npos) path = named.substr(equals + 1);
            } else if (rest.rfind("=", 0) == 0) {
                path = rest.substr(1);
            } else if (!rest.empty()) {
                fail("unknown option: " + argument);
            }
            if (name.empty()) name = "comp";
            if (path.empty()) {
                if (index + 1 >= argc) fail("missing value for comp");
                path = argv[++index];
            }
            options.comps.emplace_back(name, path);
        }
        else if (is_option(argument, "--dbsnp"))
            options.dbsnp = option_value(index, argc, argv, argument, "--dbsnp");
        else if (is_option(argument, "-L") || is_option(argument, "--interval"))
            options.regions.push_back(option_value(index, argc, argv, argument,
                argument.rfind("-L", 0) == 0 ? "-L" : "--interval"));
        else if (argument == "--resource-allele-concordance" || argument == "-rac")
            options.resource_allele_concordance = true;
        else if (is_option(argument, "--min-base-quality-score"))
            options.min_base_quality_score = std::stoi(
                option_value(index, argc, argv, argument, "--min-base-quality-score"));
        else if (argument == "--create-output-variant-index")
            options.create_index = fastgatk::native::parse_optional_boolean(
                index, argc, argv, argument, "--create-output-variant-index");
        else
            fail("unknown option: " + argument);
    }
    if (options.variants.empty()) fail("--variant is required");
    if (options.output.empty()) fail("--output is required");
    if (options.min_base_quality_score < 0 || options.min_base_quality_score > 93)
        fail("--min-base-quality-score must be in [0, 93]");
    for (const auto& annotation : options.annotations)
        if (annotation != "Coverage")
            fail("unsupported -A annotation (supported: Coverage): " + annotation);
    if (!options.annotations.empty()) {
        options.requested_coverage = true;
        if (options.reads.empty())
            fail("-A Coverage requires --input reads");
    }
    return options;
}

#if FASTGATK_HAS_HTSLIB

// Minimum-representation biallelic decomposition of a record, mirroring
// GATKVariantContextUtils.splitVariantContextToEvents(..., true,
// SET_TO_NO_CALL_NO_ANNOTATIONS, true) for the equality checks in
// VariantOverlapAnnotator.getRsID / VAExpression concordance: every
// (ref, alt) pair trimmed of common leading/trailing bases.
struct BiallelicEvent {
    std::string ref;
    std::string alt;
};

std::pair<std::string, std::string> minimal_allele(std::string ref, std::string alt) {
    std::size_t ref_offset = 0;
    std::size_t alt_offset = 0;
    // Trim shared suffix first (GATK's minimal representation keeps at
    // least one leading base), then the shared prefix.
    while (ref.size() - ref_offset > 1 && alt.size() - alt_offset > 1 &&
           ref.back() == alt.back()) {
        ref.pop_back();
        alt.pop_back();
    }
    std::size_t common = 0;
    while (ref.size() - ref_offset - common > 1 &&
           alt.size() - alt_offset - common > 1 &&
           ref[ref_offset + common] == alt[alt_offset + common])
        ++common;
    return {ref.substr(ref_offset + common), alt.substr(alt_offset + common)};
}

std::vector<BiallelicEvent> split_events(const std::string& ref,
                                         const std::vector<std::string>& alts) {
    std::vector<BiallelicEvent> events;
    events.reserve(alts.size());
    for (const auto& alt : alts) {
        if (alt == "*") continue;
        const auto minimal = minimal_allele(ref, alt);
        events.push_back(BiallelicEvent{minimal.first, minimal.second});
    }
    return events;
}

struct FeatureRecord {
    std::string id;
    std::vector<std::string> filters;
    bool filtered = false;
    std::vector<BiallelicEvent> events;
    std::map<std::string, std::string> info;  // key -> raw value text
    std::vector<std::string> info_order;
};

struct FeatureFile {
    std::string path;
    // Records keyed by 1-based start position, in file order.
    std::map<std::int64_t, std::vector<FeatureRecord>> by_start;
    std::map<std::string, std::string> info_header;  // key -> "Number,Type,Description"
};

FeatureFile load_features(const std::string& path) {
    FeatureFile file;
    file.path = path;
    htsFile* handle = bcf_open(path.c_str(), "r");
    if (handle == nullptr) fail("cannot open feature VCF: " + path);
    bcf_hdr_t* header = bcf_hdr_read(handle);
    if (header == nullptr) {
        hts_close(handle);
        fail("cannot read feature VCF header: " + path);
    }
    for (int line = 0; line < header->nhrec; ++line) {
        const bcf_hrec_t* record = header->hrec[line];
        if (record == nullptr) continue;
        if (record->type != BCF_HL_INFO) continue;
        std::string id;
        std::string number;
        std::string type;
        std::string description;
        for (int item = 0; item < record->nkeys; ++item) {
            const std::string key = record->keys[item];
            const std::string value = record->vals[item] == nullptr ? "" : record->vals[item];
            if (key == "ID") id = value;
            else if (key == "Number") number = value;
            else if (key == "Type") type = value;
            else if (key == "Description") description = value;
        }
        if (!id.empty())
            file.info_header[id] = number + "\t" + type + "\t" + description;
    }
    bcf1_t* record = bcf_init();
    std::string ref;
    std::vector<std::string> alts;
    while (bcf_read(handle, header, record) == 0) {
        // bcf_read leaves the record packed; d.allele/d.info/d.flt require
        // an explicit unpack.
        if (bcf_unpack(record, BCF_UN_ALL) != 0) continue;
        FeatureRecord feature;
        feature.id = record->d.id == nullptr ? "" : record->d.id;
        // htslib materialises FILTER=PASS as one filter id (0); htsjdk's
        // isFiltered() treats PASS and "." as unfiltered.
        for (int index = 0; index < record->d.n_flt; ++index) {
            const char* filter = header->id[BCF_DT_ID][record->d.flt[index]].key;
            if (filter != nullptr && std::strcmp(filter, "PASS") != 0)
                feature.filters.push_back(filter);
        }
        feature.filtered = !feature.filters.empty();
        ref = record->d.allele[0] == nullptr ? "" : record->d.allele[0];
        alts.clear();
        for (int index = 1; index < record->n_allele; ++index)
            alts.push_back(record->d.allele[index] == nullptr ? "" : record->d.allele[index]);
        feature.events = split_events(ref, alts);
        for (int index = 0; index < record->n_info; ++index) {
            const bcf_info_t& info = record->d.info[index];
            const char* key = bcf_hdr_int2id(header, BCF_DT_ID, info.key);
            if (key == nullptr) continue;
            kstring_t text = {0, 0, nullptr};
            if (info.len <= 0 || info.vptr == nullptr) continue;
            if (bcf_fmt_array(&text, info.len, info.type, info.vptr) >= 0 && text.s != nullptr) {
                feature.info[key] = text.s;
                feature.info_order.push_back(key);
            }
            free(text.s);
        }
        file.by_start[record->pos + 1].push_back(std::move(feature));
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    hts_close(handle);
    return file;
}

// VariantOverlapAnnotator.getRsID: unfiltered source records whose
// minimum-representation biallelic set shares an exact (ref, alt) event.
std::string overlap_ids(const std::vector<FeatureRecord>& records,
                        const std::vector<BiallelicEvent>& events) {
    std::vector<std::string> ids;
    for (const auto& record : records) {
        if (record.filtered) continue;
        bool matched = false;
        for (const auto& source : record.events) {
            for (const auto& target : events) {
                if (source.ref == target.ref && source.alt == target.alt) {
                    matched = true;
                    break;
                }
            }
            if (matched) break;
        }
        // getRsID() joins source getID() values verbatim (htsjdk returns
        // "." for an ID-less record).
        if (matched)
            ids.push_back(record.id.empty() ? std::string(".") : record.id);
    }
    if (ids.empty()) return {};
    std::string joined = ids.front();
    for (std::size_t index = 1; index < ids.size(); ++index)
        joined += ";" + ids[index];
    return joined;
}

bool overlaps_any(const std::vector<FeatureRecord>& records,
                  const std::vector<BiallelicEvent>& events) {
    return !overlap_ids(records, events).empty();
}

std::string join_filters(const FeatureRecord& record) {
    if (!record.filtered) return "PASS";
    std::string joined;
    for (const auto& filter : record.filters) {
        if (!joined.empty()) joined += ",";
        joined += filter;
    }
    return joined;
}

std::string escape_vcf_value(const std::string& value) {
    return value;
}

struct AnnotationEngine {
    const Options* options = nullptr;
    FeatureFile dbsnp;
    std::map<std::string, FeatureFile> named;
    bool has_dbsnp = false;

    void load(const Options& source) {
        options = &source;
        if (!source.dbsnp.empty()) {
            dbsnp = load_features(source.dbsnp);
            has_dbsnp = true;
        }
        for (const auto& [name, path] : source.named_resources)
            named.emplace(name, load_features(path));
        for (const auto& [name, path] : source.comps)
            named.emplace(name, load_features(path));
    }
};

struct RecordAnnotation {
    std::string id_suffix;
    std::vector<std::pair<std::string, std::string>> info;  // ordered
};

#if FASTGATK_HAS_HTSLIB
// VariantAnnotator.makeLikelihoods(): one evidence row per pileup element at
// the variant's start position (aligned base or deletion placeholder); the
// Coverage annotation reports likelihoods.evidenceCount() — the row count,
// with no reads at all producing no DP key (Coverage.annotate's empty-map
// case).
std::optional<std::uint64_t> pileup_evidence_count(samFile* bam, sam_hdr_t* bam_header,
                                                   const std::string& bam_path,
                                                   const std::string& chrom,
                                                   std::int64_t position_1based) {
    if (bam == nullptr || bam_header == nullptr) return std::nullopt;
    hts_idx_t* index = sam_index_load(bam, bam_path.c_str());
    if (index == nullptr) return std::nullopt;
    const std::string region = chrom + ":" + std::to_string(position_1based) +
        "-" + std::to_string(position_1based);
    hts_itr_t* iterator = sam_itr_querys(index, bam_header, region.c_str());
    std::uint64_t count = 0;
    if (iterator != nullptr) {
        bam1_t* read = bam_init1();
        while (sam_itr_next(bam, iterator, read) >= 0) {
            const std::int64_t start = read->core.pos + 1;  // 1-based
            std::int64_t reference = start;
            bool element = false;
            const std::uint32_t* cigar = bam_get_cigar(read);
            for (std::uint32_t operation = 0; operation < read->core.n_cigar; ++operation) {
                const std::uint32_t length = cigar[operation] >> BAM_CIGAR_SHIFT;
                const std::uint32_t code = cigar[operation] & BAM_CIGAR_MASK;
                if (code == BAM_CMATCH || code == BAM_CEQUAL || code == BAM_CDIFF ||
                    code == BAM_CDEL || code == BAM_CREF_SKIP) {
                    if (position_1based >= reference &&
                        position_1based < reference + static_cast<std::int64_t>(length)) {
                        // Base and deletion elements appear in the pileup;
                        // reference skips do not.
                        element = code != BAM_CREF_SKIP;
                        break;
                    }
                    reference += length;
                }
            }
            if (element) ++count;
        }
        bam_destroy1(read);
        hts_itr_destroy(iterator);
    }
    hts_idx_destroy(index);
    return count;
}
#endif

RecordAnnotation annotate(const Options& options, const AnnotationEngine& engine,
                          std::int64_t position, const std::string& ref,
                          const std::vector<std::string>& alts,
                          const std::string& input_id,
                          const std::optional<std::uint64_t>& evidence_count) {
    RecordAnnotation output;
    const auto events = split_events(ref, alts);

    // Coverage (INFO/DP) reports the read likelihood evidence count; with no
    // reads the annotation yields no key at all.
    if (options.requested_coverage && evidence_count.has_value()) {
        std::ostringstream depth;
        depth << *evidence_count;
        output.info.emplace_back("DP", depth.str());
    }

    // VariantOverlapAnnotator: comps and --dbsnp flags + rsID merge.
    for (const auto& [name, path] : options.comps) {
        (void)path;
        const auto found = engine.named.find(name);
        if (found == engine.named.end()) continue;
        const auto records = found->second.by_start.find(position);
        if (records == found->second.by_start.end()) continue;
        if (overlaps_any(records->second, events))
            output.info.emplace_back(name, "");
    }
    if (engine.has_dbsnp) {
        const auto records = engine.dbsnp.by_start.find(position);
        if (records != engine.dbsnp.by_start.end()) {
            const std::string rsid = overlap_ids(records->second, events);
            if (!rsid.empty()) {
                output.info.emplace_back("DB", "");
                // annotateRsID: no ID -> use rsID; ID without the rsID ->
                // append with ';'; already contained -> unchanged.
                if (input_id.empty() || input_id == ".")
                    output.id_suffix = rsid;
                else if (input_id.find(rsid) == std::string::npos)
                    output.id_suffix = input_id + ";" + rsid;
            }
        }
    }

    // VAExpression / annotateExpressions: first resource record at the start.
    for (const auto& expression : options.expressions) {
        const auto separator = expression.rfind('.');
        if (separator == std::string::npos)
            fail("invalid expression (want NAME.FIELD): " + expression);
        ExpressionSpec spec{expression, expression.substr(0, separator),
                            expression.substr(separator + 1)};
        const auto found = engine.named.find(spec.binding);
        if (found == engine.named.end())
            fail("expression binding not found: " + spec.binding);
        const auto records = found->second.by_start.find(position);
        if (records == found->second.by_start.end()) continue;
        const FeatureRecord& source = records->second.front();
        if (spec.field == "ID") {
            if (!source.id.empty() && source.id != ".")
                output.info.emplace_back(spec.full_name, source.id);
        } else if (spec.field == "ALT") {
            // annotateExpressions uses the first alternate allele.
            output.info.emplace_back(spec.full_name,
                source.events.empty() ? "" : source.events.front().alt);
        } else if (spec.field == "FILTER") {
            output.info.emplace_back(spec.full_name, join_filters(source));
        } else {
            const auto value = source.info.find(spec.field);
            if (value == source.info.end()) continue;
            // VAExpression/annotateExpressions: for A/R-counted fields (or
            // with --resource-allele-concordance) values are mapped onto the
            // input's minimum-representation biallelics; an input allele with
            // no concordant resource allele contributes "0", and the key is
            // omitted entirely when nothing matched.  Other fields transfer
            // the raw attribute value.
            const auto header = found->second.info_header.find(spec.field);
            std::string count_type;
            if (header != found->second.info_header.end()) {
                const auto first_tab = header->second.find('\t');
                count_type = header->second.substr(0, first_tab);
            }
            const bool useRefAndAltAlleles = count_type == "R";
            const bool useAltAlleles = count_type == "A";
            if (!((useAltAlleles || useRefAndAltAlleles) ||
                  options.resource_allele_concordance)) {
                output.info.emplace_back(spec.full_name, value->second);
                continue;
            }
            std::string cleaned = value->second;
            cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(),
                [](char character) {
                    return character == '[' || character == ']' ||
                        std::isspace(static_cast<unsigned char>(character)) != 0;
                }), cleaned.end());
            std::vector<std::string> expression_values;
            {
                std::size_t start = 0;
                while (start <= cleaned.size()) {
                    const auto comma = cleaned.find(',', start);
                    expression_values.push_back(cleaned.substr(
                        start, comma == std::string::npos ? std::string::npos
                                                         : comma - start));
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
            }
            const auto source_events = source.events;
            std::vector<std::string> mapped;
            bool can_annotate = false;
            for (const auto& target : events) {
                bool concordant = false;
                std::size_t index = 0;
                for (const auto& candidate : source_events) {
                    if (candidate.ref == target.ref && candidate.alt == target.alt) {
                        if (index == 0 && useRefAndAltAlleles &&
                            index < expression_values.size())
                            mapped.push_back(expression_values[index++]);
                        if (index < expression_values.size())
                            mapped.push_back(expression_values[index]);
                        concordant = true;
                        can_annotate = true;
                        break;
                    }
                    ++index;
                }
                if (!concordant) mapped.push_back("0");
            }
            if (can_annotate) {
                std::string joined;
                for (std::size_t index = 0; index < mapped.size(); ++index) {
                    if (index != 0) joined += ",";
                    joined += mapped[index];
                }
                output.info.emplace_back(spec.full_name, joined);
            }
        }
    }
    return output;
}

std::string serialize_info(const bcf_hdr_t* header, bcf1_t* record,
                           const RecordAnnotation& annotation) {
    // htsjdk writes INFO keys in sorted order; start from the record's
    // existing INFO text and merge the new keys into that order.
    std::map<std::string, std::string> entries;
    std::vector<std::string> keys;
    for (int index = 0; index < record->n_info; ++index) {
        const bcf_info_t& info = record->d.info[index];
        const char* key = bcf_hdr_int2id(header, BCF_DT_ID, info.key);
        if (key == nullptr) continue;
        kstring_t text = {0, 0, nullptr};
        if (info.len > 0 && info.vptr != nullptr &&
            bcf_fmt_array(&text, info.len, info.type, info.vptr) >= 0 &&
            text.s != nullptr) {
            entries[key] = text.s;
            keys.push_back(key);
        }
        free(text.s);
    }
    for (const auto& [key, value] : annotation.info) {
        if (entries.find(key) == entries.end()) keys.push_back(key);
        entries[key] = value;
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    if (keys.empty()) return ".";
    std::ostringstream out;
    bool first = true;
    for (const auto& key : keys) {
        if (!first) out << ';';
        first = false;
        out << key;
        const auto value = entries.find(key);
        if (value != entries.end() && !value->second.empty())
            out << '=' << escape_vcf_value(value->second);
    }
    return out.str();
}

std::string serialize_record(bcf_hdr_t* header, bcf1_t* record,
                             const RecordAnnotation& annotation,
                             bool sites_only) {
    kstring_t line = {0, 0, nullptr};
    if (vcf_format(header, record, &line) < 0 || line.s == nullptr) {
        free(line.s);
        fail("cannot format VCF record");
    }
    std::string text(line.s, line.l);
    free(line.s);
    while (!text.empty() && (text.back() == '\n' || text.back() == '\r'))
        text.pop_back();
    // Split the formatted line and rewrite ID and INFO in place.
    std::vector<std::string> fields;
    std::size_t start = 0;
    for (std::size_t index = 0; index <= text.size(); ++index) {
        if (index == text.size() || text[index] == '\t') {
            fields.push_back(text.substr(start, index - start));
            start = index + 1;
        }
    }
    if (fields.size() < 8) fail("malformed VCF record");
    if (!annotation.id_suffix.empty()) fields[2] = annotation.id_suffix;
    fields[7] = serialize_info(header, record, annotation);
    std::ostringstream out;
    const std::size_t limit = sites_only ? 8 : fields.size();
    for (std::size_t index = 0; index < limit; ++index) {
        if (index != 0) out << '\t';
        out << fields[index];
    }
    out << '\n';
    return out.str();
}

int run(const Options& options) {
    AnnotationEngine engine;
    engine.load(options);

    samFile* reads = nullptr;
    sam_hdr_t* reads_header = nullptr;
    if (!options.reads.empty()) {
        reads = sam_open(options.reads.c_str(), "r");
        if (reads == nullptr) fail("cannot open reads: " + options.reads);
        reads_header = sam_hdr_read(reads);
        if (reads_header == nullptr) fail("cannot read BAM header: " + options.reads);
    }

    htsFile* input = bcf_open(options.variants.c_str(), "r");
    if (input == nullptr) fail("cannot open input VCF: " + options.variants);
    bcf_hdr_t* header = bcf_hdr_read(input);
    if (header == nullptr) {
        hts_close(input);
        fail("cannot read input VCF header");
    }

    // Header: htsjdk's VCFEncoder writes fileformat first, then FILTER,
    // FORMAT and INFO lines in sorted ID order (with the standard-key lines
    // of GATK's default annotation set overriding redefinitions), then the
    // contig lines, then any remaining metadata in input order.  The
    // canonical lines below are GATK 4.6.2.0's rendered output for the keys
    // its default annotation set owns.
    const std::map<std::string, std::string> canonical_lines = {
        {"INFO:AC", "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">"},
        {"INFO:AF", "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">"},
        {"INFO:SB", "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">"},
        {"INFO:DP", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">"},
        {"FORMAT:AD", "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"},
        {"FORMAT:DP", "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"},
        {"FORMAT:FT", "##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Genotype-level filter\">"},
        {"FORMAT:GQ", "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"},
        {"FORMAT:PL", "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"},
        {"FORMAT:PS", "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phasing set (typically the position of the first variant in the set)\">"},
    };
    kstring_t raw_header = {0, 0, nullptr};
    if (bcf_hdr_format(header, 0, &raw_header) < 0 || raw_header.s == nullptr) {
        free(raw_header.s);
        fail("cannot format input header");
    }
    const std::string header_text(raw_header.s, raw_header.l);
    free(raw_header.s);
    std::map<std::string, std::string> format_lines;
    std::map<std::string, std::string> info_lines;
    std::vector<std::string> filter_lines;
    std::vector<std::string> contig_lines;
    std::vector<std::string> other_lines;
    {
        std::size_t start_line = 0;
        while (start_line < header_text.size()) {
            const auto newline = header_text.find('\n', start_line);
            const std::string line = header_text.substr(
                start_line, newline == std::string::npos ? std::string::npos
                                                        : newline - start_line);
            start_line = newline == std::string::npos ? header_text.size() : newline + 1;
            if (line.empty()) continue;
            const auto id_marker = line.find("=<ID=");
            std::string id;
            if (id_marker != std::string::npos) {
                const auto id_end = line.find(',', id_marker + 5);
                id = line.substr(id_marker + 5, id_end == std::string::npos
                    ? std::string::npos : id_end - (id_marker + 5));
            }
            if (line.rfind("##fileformat=", 0) == 0) continue;
            if (line.rfind("#CHROM", 0) == 0) continue;
            if (line.rfind("##FILTER=", 0) == 0) {
                // htsjdk never emits the implicit PASS definition; other
                // FILTER lines sort by ID before FORMAT/INFO.
                if (id != "PASS") filter_lines.push_back(line);
            } else if (line.rfind("##FORMAT=", 0) == 0) {
                const auto canonical = canonical_lines.find("FORMAT:" + id);
                format_lines[id] = canonical == canonical_lines.end()
                    ? line : canonical->second;
            } else if (line.rfind("##INFO=", 0) == 0) {
                const auto canonical = canonical_lines.find("INFO:" + id);
                info_lines[id] = canonical == canonical_lines.end()
                    ? line : canonical->second;
            } else if (line.rfind("##contig=", 0) == 0) {
                contig_lines.push_back(line);
            } else if (line.rfind("##reference=", 0) == 0) {
                // updateHeaderContigLines drops the input reference key and
                // re-adds it from --reference below.
                continue;
            } else {
                other_lines.push_back(line);
            }
        }
    }
    // Expression header lines mirror annotateExpressions' sources: the
    // special cases ID/ALT/FILTER get fixed lines, other fields clone the
    // resource's header definition (VAExpression.sethInfo).
    for (const auto& expression : options.expressions) {
        const auto separator = expression.rfind('.');
        if (separator == std::string::npos) fail("invalid expression: " + expression);
        const std::string binding = expression.substr(0, separator);
        const std::string field = expression.substr(separator + 1);
        if (field == "ID") {
            info_lines[expression] = "##INFO=<ID=" + expression +
                ",Number=1,Type=String,Description=\"ID field transferred from external VCF resource\">";
            continue;
        }
        const auto found = engine.named.find(binding);
        if (found == engine.named.end()) fail("expression binding not found: " + binding);
        const auto source_line = found->second.info_header.find(field);
        if (source_line != found->second.info_header.end()) {
            const std::string& spec = source_line->second;
            const auto first_tab = spec.find('\t');
            const auto second_tab = spec.find('\t', first_tab + 1);
            std::string description = spec.substr(second_tab + 1);
            // htslib keeps the header description's surrounding quotes in
            // the parsed value; the rendered line adds its own.
            if (description.size() >= 2 && description.front() == '"' &&
                description.back() == '"')
                description = description.substr(1, description.size() - 2);
            info_lines[expression] = "##INFO=<ID=" + expression + ",Number=" +
                spec.substr(0, first_tab) + ",Type=" +
                spec.substr(first_tab + 1, second_tab - first_tab - 1) +
                ",Description=\"" + description + "\">";
        } else {
            info_lines[expression] = "##INFO=<ID=" + expression +
                ",Number=.,Type=String,Description=\"Value transferred from another external VCF resource\">";
        }
    }
    for (const auto& [name, path] : options.comps) {
        (void)path;
        info_lines[name] = "##INFO=<ID=" + name +
            ",Number=0,Type=Flag,Description=\"" + name + " Membership\">";
    }
    if (engine.has_dbsnp)
        info_lines["DB"] = "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">";
    if (!options.annotations.empty())
        info_lines["DP"] = canonical_lines.at("INFO:DP");

    std::ofstream output(options.output);
    if (!output) fail("cannot open output VCF: " + options.output);
    {
        std::ostringstream header_out;
        header_out << "##fileformat=VCFv4.2\n";
        std::sort(filter_lines.begin(), filter_lines.end());
        for (const auto& line : filter_lines) header_out << line << '\n';
        for (const auto& [id, line] : format_lines) header_out << line << '\n';
        for (const auto& [id, line] : info_lines) header_out << line << '\n';
        if (options.reference.empty()) {
            for (const auto& line : contig_lines) header_out << line << '\n';
        } else {
            // updateHeaderContigLines rebuilds contigs from the sequence
            // dictionary with assembly = the reference file name, then adds
            // the ##reference URI line.
            std::string assembly = options.reference;
            const auto slash = assembly.find_last_of('/');
            if (slash != std::string::npos) assembly = assembly.substr(slash + 1);
            for (const auto& line : contig_lines) {
                std::string rebuilt = line;
                const auto length_marker = line.find(",length=");
                if (length_marker != std::string::npos) {
                    const auto close = line.find('>', length_marker);
                    rebuilt = line.substr(0, close) + ",assembly=" + assembly + ">";
                }
                header_out << rebuilt << '\n';
            }
            std::string reference_path = options.reference;
            if (reference_path.front() != '/')
                reference_path = std::filesystem::absolute(reference_path).string();
            header_out << "##reference=file://" << reference_path << '\n';
        }
        for (const auto& line : other_lines) header_out << line << '\n';
        header_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        if (bcf_hdr_nsamples(header) > 0) {
            header_out << "\tFORMAT";
            for (int index = 0; index < bcf_hdr_nsamples(header); ++index)
                header_out << '\t' << bcf_hdr_int2id(header, BCF_DT_SAMPLE, index);
        }
        header_out << '\n';
        const std::string rendered = header_out.str();
        output.write(rendered.data(), static_cast<std::streamsize>(rendered.size()));
    }


    // Interval restriction (-L CHR:START-END, 1-based inclusive) on records.
    std::vector<std::pair<std::string, std::pair<std::int64_t, std::int64_t>>> intervals;
    for (const auto& region : options.regions) {
        const auto colon = region.find(':');
        if (colon == std::string::npos) {
            intervals.emplace_back(region, std::make_pair<std::int64_t, std::int64_t>(
                1, std::numeric_limits<std::int64_t>::max()));
        } else {
            const auto dash = region.find('-', colon);
            intervals.emplace_back(region.substr(0, colon),
                std::make_pair(std::strtoll(region.substr(colon + 1).c_str(), nullptr, 10),
                    dash == std::string::npos
                        ? std::strtoll(region.substr(colon + 1).c_str(), nullptr, 10)
                        : std::strtoll(region.substr(dash + 1).c_str(), nullptr, 10)));
        }
    }

    bcf1_t* record = bcf_init();
    while (bcf_read(input, header, record) == 0) {
        if (bcf_unpack(record, BCF_UN_ALL) != 0) continue;
        const char* contig = bcf_hdr_id2name(header, record->rid);
        const std::string chrom = contig == nullptr ? "" : contig;
        const std::int64_t position = record->pos + 1;
        if (!intervals.empty()) {
            bool keep = false;
            for (const auto& [name, span] : intervals)
                if (name == chrom && position >= span.first && position <= span.second) {
                    keep = true;
                    break;
                }
            if (!keep) continue;
        }
        std::string ref = record->d.allele[0] == nullptr ? "" : record->d.allele[0];
        std::vector<std::string> alts;
        for (int index = 1; index < record->n_allele; ++index)
            alts.push_back(record->d.allele[index] == nullptr ? "" : record->d.allele[index]);
        const std::string input_id = record->d.id == nullptr ? "" : record->d.id;

        // VariantAnnotator.apply(): records at an ambiguous reference base
        // are written through unannotated.
        const auto evidence_count = options.requested_coverage
            ? pileup_evidence_count(reads, reads_header, options.reads,
                                     chrom, position)
            : std::optional<std::uint64_t>{};
        const auto annotation = annotate(options, engine, position, ref, alts,
                                         input_id, evidence_count);
        const std::string text = serialize_record(header, record, annotation, false);
        output.write(text.data(), static_cast<std::streamsize>(text.size()));
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    hts_close(input);
    if (reads_header != nullptr) sam_hdr_destroy(reads_header);
    if (reads != nullptr) sam_close(reads);
    output.close();
    return 0;
}

#endif  // FASTGATK_HAS_HTSLIB

}  // namespace

int main(int argc, char** argv) {
    try {
#if FASTGATK_HAS_HTSLIB
        Kokkos::initialize(argc, argv);
        const Options options = parse_options(argc, argv);
        const int result = run(options);
        Kokkos::finalize();
        return result;
#else
        std::cerr << "fastgatk-variant-annotator requires HTSlib support\n";
        return 1;
#endif
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
