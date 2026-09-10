#include "fastgatk/io/hts_reader.hpp"
#include "fastgatk/io/flow_codec.hpp"

#include <cstdint>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

struct Counts {
    std::size_t records = 0;
    std::size_t bases = 0;
    std::size_t projected_bases = 0;
    std::size_t cigar_ops = 0;
    bool indexed = false;
    bool index_available = false;
};

Counts verify(const std::string& input, const std::string& region,
              const bool require_index_capability = false) {
    fastgatk::io::HtsReader reader(input, {}, region, 17);
    require(reader.batch_records() == 17, "initial HTS batch size was not retained");
    reader.set_batch_records(9);
    require(reader.batch_records() == 9, "HTS batch size update was not retained");
    fastgatk::io::ReadBatch batch;
    Counts counts;
    counts.indexed = reader.indexed();
    counts.index_available = reader.has_index();
    if (require_index_capability)
        require(counts.index_available,
                "unbounded HTS reader did not expose an available BAM/CRAM index");
    while (reader.next(batch)) {
        require(batch.offsets.size() == batch.records() + 1,
                "ReadBatch base offsets do not match records");
        require(batch.cigar_offsets.size() == batch.records() + 1,
                "HTSlib batch has no per-record CIGAR offsets");
        require(batch.name_offsets.size() == batch.records() + 1,
                "HTSlib batch has no per-record name offsets");
        require(batch.read_group_offsets.size() == batch.records() + 1,
                "HTSlib batch has no per-record RG offsets");
        require(batch.original_alignment_offsets.size() == batch.records() + 1 &&
                    batch.mate_contig_offsets.size() == batch.records() + 1,
                "HTSlib batch has no per-record OA/XM offsets");
        require(batch.original_alignment_present.size() == batch.records() &&
                    batch.mate_contig_present.size() == batch.records(),
                "HTSlib batch has no per-record OA/XM presence metadata");
        require(batch.insertion_quality_offsets.size() == batch.records() + 1 &&
                    batch.deletion_quality_offsets.size() == batch.records() + 1 &&
                    batch.indel_quality_layout_valid(),
                "HTSlib batch has no per-record BI/BD offsets");
        require(batch.flow_tp_offsets.size() == batch.records() + 1 &&
                    batch.flow_t0_offsets.size() == batch.records() + 1 &&
                    batch.flow_order_offsets.size() == batch.records() + 1 &&
                    batch.flow_max_hmer.size() == batch.records(),
                "HTSlib batch has no per-record flow metadata offsets");
        require(batch.flags.size() == batch.records() &&
                    batch.mate_tids.size() == batch.records() &&
                    batch.mate_positions.size() == batch.records() &&
                    batch.template_lengths.size() == batch.records(),
                "HTSlib batch core metadata does not match records");
        require(batch.cigar_layout_valid(), "malformed packed CIGAR layout");
        require(batch.bytes() >= batch.bases.size() + batch.qualities.size(),
                "ReadBatch byte accounting under-reported base payload");
        counts.records += batch.records();
        counts.bases += batch.bases_count();
        counts.cigar_ops += batch.cigar_ops.size();
        for (std::size_t record = 0; record < batch.records(); ++record) {
            const auto begin = batch.offsets[record];
            const auto end = batch.offsets[record + 1];
            const auto insertion_quality_begin = batch.insertion_quality_offsets[record];
            const auto insertion_quality_end = batch.insertion_quality_offsets[record + 1];
            const auto deletion_quality_begin = batch.deletion_quality_offsets[record];
            const auto deletion_quality_end = batch.deletion_quality_offsets[record + 1];
            if (insertion_quality_end != insertion_quality_begin ||
                deletion_quality_end != deletion_quality_begin) {
                require(insertion_quality_end - insertion_quality_begin == end - begin &&
                            deletion_quality_end - deletion_quality_begin == end - begin,
                        "BI/BD quality arrays do not match read bases");
                require(batch.insertion_qualities[insertion_quality_begin] == 0 &&
                            batch.insertion_qualities[insertion_quality_begin + 1] == 20 &&
                            batch.deletion_qualities[deletion_quality_begin] == 3 &&
                            batch.deletion_qualities[deletion_quality_begin + 1] == 4,
                        "BI/BD FASTQ qualities were not decoded to Phred bytes");
            }
            const auto oa_begin = batch.original_alignment_offsets[record];
            const auto oa_end = batch.original_alignment_offsets[record + 1];
            const auto xm_begin = batch.mate_contig_offsets[record];
            const auto xm_end = batch.mate_contig_offsets[record + 1];
            if (oa_end > oa_begin || xm_end > xm_begin) {
                require(std::string(batch.original_alignments.begin() + oa_begin,
                                    batch.original_alignments.begin() + oa_end) ==
                            "chr1,1,+,2M,60,0" &&
                        std::string(batch.mate_contigs.begin() + xm_begin,
                                    batch.mate_contigs.begin() + xm_end) == "chr1",
                        "OA/XM auxiliary tags were not decoded into ReadBatch");
            }
            const auto flow_begin = batch.flow_order_offsets[record];
            const auto flow_end = batch.flow_order_offsets[record + 1];
            if (flow_end > flow_begin) {
                require(batch.flow_tp_offsets[record + 1] - batch.flow_tp_offsets[record] == end - begin,
                        "flow tp array does not match read bases");
                require(batch.flow_max_hmer[record] > 0, "flow RG max hmer is missing");
                const auto tags = fastgatk::io::flow_tags_from_batch(batch, record, true);
                const auto decoded = fastgatk::io::decode_flow_read(tags);
                require(decoded.key.size() == decoded.flow_order.size() &&
                            decoded.probabilities.size() == decoded.key.size() * 256,
                        "flow tag decoder returned malformed flat matrix");
            }
            const auto reference_end = fastgatk::io::reference_end(batch, record);
            if (batch.positions[record] >= 0)
                require(reference_end >= batch.positions[record], "reference end precedes read start");
            for (std::size_t offset = begin; offset < end; ++offset) {
                fastgatk::io::ReadProjection projection;
                if (!fastgatk::io::project_read_offset(batch, record, offset - begin, projection)) continue;
                require(projection.reference_position >= batch.positions[record],
                        "projected coordinate precedes read start");
                require(projection.reference_position < reference_end,
                        "projected coordinate is outside reference span");
                require(fastgatk::io::CigarOp{1, projection.operation}.projects_base(),
                        "projected operation does not consume one reference base");
                ++counts.projected_bases;
            }
        }
    }
    require(counts.records != 0, "fixture produced no records");
    require(reader.records_read() >= counts.records,
            "HTS reader records_read telemetry under-counted decoded records");
    require(counts.projected_bases <= counts.bases, "projection emitted more bases than input");
    return counts;
}

std::vector<std::int32_t> read_flow_key_oracle(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open flow key oracle: " + path.string());
    std::vector<std::int32_t> key;
    std::int32_t value = 0;
    while (input >> value) key.push_back(value);
    return key;
}

void verify_flow_oracle(const std::string& input_path,
                        const std::filesystem::path& expected_dir,
                        bool use_t0,
                        const std::string& flow_order_override = {}) {
    fastgatk::io::HtsReader reader(input_path, std::string{}, std::string{}, 17);
    fastgatk::io::ReadBatch batch;
    std::size_t record = 0;
    std::size_t checked = 0;
    while (reader.next(batch)) {
        for (std::size_t local = 0; local < batch.records(); ++local, ++record) {
            // Match FlowBasedArgumentCollection defaults used by the GATK
            // fixture: boundary flows are spread across possible hmer calls.
            auto tags = fastgatk::io::flow_tags_from_batch(batch, local, use_t0, 0.001, false);
            // The historical GATK fixture intentionally supplies the flow
            // order as a test argument (its BAM RG header uses another
            // order). Production callers leave this empty and use @RG FO.
            if (!flow_order_override.empty()) tags.flow_order = flow_order_override;
            const auto decoded = fastgatk::io::decode_flow_read(tags);
            const auto stem = std::string(use_t0 ? "sample.t0." : "sample.") +
                std::to_string(record);
            const auto expected_key = read_flow_key_oracle(expected_dir / (stem + ".key.txt"));
            // GATK's checked-in flow fixture is emitted after phase
            // normalization: zero-length flows before the first called flow
            // are omitted from both key and matrix files. HtsReader retains
            // the literal RG cycle, so perform that representation-only
            // normalization at the oracle boundary.
            std::size_t leading_zero_flows = 0;
            while (leading_zero_flows < decoded.key.size() && decoded.key[leading_zero_flows] == 0)
                ++leading_zero_flows;
            const auto normalized_key_begin = decoded.key.begin() + leading_zero_flows;
            const auto normalized_key_end = decoded.key.end();
            if (static_cast<std::size_t>(normalized_key_end - normalized_key_begin) != expected_key.size() ||
                !std::equal(expected_key.begin(), expected_key.end(), normalized_key_begin)) {
                std::ostringstream detail;
                detail << "FlowBasedKeyCodec key mismatch at record " << record <<
                    " flow_order=" << tags.flow_order << " observed_len=" << decoded.key.size()
                    << " expected_len=" << expected_key.size() << " observed=";
                for (const auto value : decoded.key) detail << value << ',';
                detail << " expected=";
                for (const auto value : expected_key) detail << value << ',';
                throw std::runtime_error(detail.str());
            }
            std::ifstream matrix(expected_dir / (stem + ".matrix.txt"));
            require(static_cast<bool>(matrix), "cannot open flow matrix oracle: " + stem);
            std::string line;
            require(static_cast<bool>(std::getline(matrix, line)),
                    "empty flow matrix oracle: " + stem);
            std::size_t matrix_rows = 0;
            while (std::getline(matrix, line)) {
                if (line.empty()) continue;
                if (line == "C,R,F,B,Bi,Q,ti") continue;
                const auto separator = line.rfind(' ');
                if (separator == std::string::npos) throw std::runtime_error(
                    "malformed flow matrix oracle row: " + line);
                const auto comma = line.find(',');
                const auto second_comma = line.find(',', comma + 1);
                if (comma == std::string::npos || second_comma == std::string::npos)
                    throw std::runtime_error("malformed flow matrix coordinates: " + line);
                const auto flow = static_cast<std::size_t>(std::stoul(line.substr(0, comma)));
                const auto hmer = static_cast<std::size_t>(std::stoul(
                    line.substr(comma + 1, second_comma - comma - 1)));
                const auto expected = std::stod(line.substr(separator + 1));
                require(flow < decoded.key.size() && hmer < 256,
                        "flow matrix oracle coordinate exceeds decoded matrix");
                const auto observed_flow = flow + leading_zero_flows;
                require(observed_flow < decoded.key.size(),
                        "flow matrix oracle phase exceeds decoded matrix");
                const auto observed = decoded.probabilities[observed_flow * 256 + hmer];
                if (!(std::isfinite(observed) && std::abs(observed - expected) <= 6.0e-6)) {
                    std::ostringstream detail;
                    detail << "FlowBasedRead probability mismatch at record " << record
                           << " flow=" << flow << " hmer=" << hmer
                           << " observed=" << observed << " expected=" << expected;
                    throw std::runtime_error(detail.str());
                }
                ++matrix_rows;
            }
            require(matrix_rows == expected_key.size() * 13,
                    "flow matrix oracle row count does not match maxHmer=12");
            ++checked;
        }
    }
    require(checked > 0, "flow oracle fixture produced no records");
}

void verify_projection_edges() {
    using fastgatk::io::CigarOp;
    using fastgatk::io::CigarOpCode;
    fastgatk::io::ReadBatch batch;
    batch.offsets = {0, 12};
    batch.bases.assign(12, static_cast<std::uint8_t>('A'));
    batch.qualities.assign(12, 30);
    batch.positions = {100};
    batch.tids = {0};
    batch.mapq = {60};
    batch.cigar_offsets = {0, 5};
    batch.cigar_ops = {
        CigarOp{5, CigarOpCode::Match}.pack(),
        CigarOp{2, CigarOpCode::Insertion}.pack(),
        CigarOp{3, CigarOpCode::SequenceMatch}.pack(),
        CigarOp{1, CigarOpCode::Deletion}.pack(),
        CigarOp{2, CigarOpCode::SequenceMismatch}.pack(),
    };
    require(batch.cigar_layout_valid(), "synthetic CIGAR layout rejected");
    require(fastgatk::io::reference_end(batch, 0) == 111,
            "synthetic CIGAR reference end is incorrect");
    for (std::size_t offset = 0; offset < 12; ++offset) {
        fastgatk::io::ReadProjection projection;
        const bool projected = fastgatk::io::project_read_offset(batch, 0, offset, projection);
        if (offset == 5 || offset == 6) {
            require(!projected, "insertion was projected as a reference base");
        } else {
            require(projected, "reference-consuming CIGAR base was dropped");
            const std::int64_t expected = offset < 5 ? 100 + offset :
                (offset < 10 ? 100 + offset - 2 : 100 + offset - 1);
            require(projection.reference_position == expected,
                    "synthetic CIGAR projection coordinate is incorrect");
        }
    }

    // Existing callers may construct a ReadBatch without optional metadata;
    // the compatibility projection remains the old contiguous mapping.
    fastgatk::io::ReadBatch legacy;
    legacy.offsets = {0, 3};
    legacy.bases.assign(3, static_cast<std::uint8_t>('C'));
    legacy.qualities.assign(3, 30);
    legacy.positions = {7};
    legacy.tids = {0};
    legacy.mapq = {20};
    require(fastgatk::io::reference_end(legacy, 0) == 10,
            "legacy reference end changed");
    fastgatk::io::ReadProjection projection;
    require(fastgatk::io::project_read_offset(legacy, 0, 2, projection) &&
                projection.reference_position == 9,
            "legacy contiguous projection changed");

    // Per-record CIGAR predicates must not make a malformed sibling poison a
    // valid read's coordinate projection.  The aggregate layout remains
    // false (useful for strict batch-boundary validation), while the first
    // record is independently usable by read filters and walkers.
    fastgatk::io::ReadBatch mixed;
    mixed.offsets = {0, 3, 6};
    mixed.bases.assign(6, static_cast<std::uint8_t>('A'));
    mixed.qualities.assign(6, 30);
    mixed.positions = {100, 100};
    mixed.tids = {0, 0};
    mixed.mapq = {60, 60};
    mixed.flags = {0, 0};
    mixed.cigar_offsets = {0, 1, 2};
    mixed.cigar_ops = {
        fastgatk::io::CigarOp{3, fastgatk::io::CigarOpCode::Match}.pack(),
        fastgatk::io::CigarOp{3, static_cast<fastgatk::io::CigarOpCode>(15)}.pack(),
    };
    require(!mixed.cigar_layout_valid(), "aggregate malformed CIGAR unexpectedly accepted");
    require(mixed.cigar_record_layout_valid(0) && !mixed.cigar_record_layout_valid(1),
            "per-record CIGAR validation did not isolate malformed sibling");
    require(fastgatk::io::reference_end(mixed, 0) == 103,
            "per-record reference end was poisoned by malformed sibling");
    require(fastgatk::io::project_read_offset(mixed, 0, 2, projection) &&
                projection.reference_position == 102,
            "per-record projection was poisoned by malformed sibling");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 3 && std::string(argv[2]) == "--flow-oracle") {
            if (argc < 4 || argc > 6)
                throw std::invalid_argument("usage: hts_reader_projection_test BAM --flow-oracle EXPECTED_DIR [t0] [flow-order=FO]");
            bool use_t0 = false;
            std::string flow_order_override;
            for (int index = 4; index < argc; ++index) {
                const std::string option = argv[index];
                if (option == "t0") use_t0 = true;
                else if (option.rfind("flow-order=", 0) == 0)
                    flow_order_override = option.substr(std::string("flow-order=").size());
                else throw std::invalid_argument("flow oracle option must be t0 or flow-order=FO");
            }
            verify_flow_oracle(argv[1], argv[3], use_t0, flow_order_override);
            std::cout << "{\"status\":\"pass\",\"suite\":\"flow-codec-gatk-oracle\"}\n";
            return 0;
        }
        if (argc < 2 || argc > 3)
            throw std::invalid_argument(
                "usage: hts_reader_projection_test BAM [region|--require-index]");
        verify_projection_edges();
        const bool require_index_capability = argc == 3 &&
            std::string(argv[2]) == "--require-index";
        const auto counts = verify(argv[1], argc == 3 && !require_index_capability
                                             ? argv[2] : std::string{},
                                   require_index_capability);
        std::cout << "{\"status\":\"pass\",\"records\":" << counts.records
                  << ",\"bases\":" << counts.bases
                  << ",\"projected_bases\":" << counts.projected_bases
                  << ",\"cigar_ops\":" << counts.cigar_ops
                  << ",\"indexed\":" << (counts.indexed ? "true" : "false")
                  << ",\"index_available\":"
                  << (counts.index_available ? "true" : "false") << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
