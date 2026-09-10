#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <string>
#include <vector>

namespace fastgatk::io {

struct ReadBatch;

// Host-side representation of the production flow tags used by GATK's
// FlowBasedRead.  `tp` is the signed per-base offset from the called
// homopolymer length; `t0_phred` is optional and is only consumed when
// use_t0_tag is true.  The Kokkos layer never sees these tags or strings.
struct FlowReadTags {
    std::vector<std::uint8_t> bases;
    std::vector<std::uint8_t> qualities;
    std::vector<std::int8_t> tp;
    std::vector<std::uint8_t> t0_phred;
    std::string flow_order;
    std::size_t max_hmer = 12;
    double filling_value = 0.001;
    bool use_t0_tag = false;
    bool keep_boundary_flows = false;
};

// Decoded flow key and calibrated probability table.  `probabilities` is
// flattened as [flow][256] so it can be passed directly to
// FlowPairHmmRead::probabilities.  Values above max_hmer repeat the max_hmer
// row, matching FlowBasedRead.getProb()'s hmer clamp.
struct DecodedFlowRead {
    std::vector<std::int32_t> key;
    std::vector<std::uint8_t> flow_order;
    std::vector<double> probabilities;
    // Semantic decoder bound; the transport matrix remains 256 rows but
    // FlowBasedRead spreads only through maxHmer.
    std::size_t max_hmer = 12;
    double filling_value = 0.001;
};

// Host-side flow key encoder shared by read decoding and haplotype materialization.
// The returned order is one flow symbol per key element (not the repeated FO
// string), which is the representation consumed by the flow PairHMM kernel.
struct EncodedFlowKey {
    std::vector<std::int32_t> key;
    std::vector<std::uint8_t> flow_order;
    // Index of the first base represented by each flow.  These are the
    // FlowBasedRead/FlowBasedHaplotype `flow2base` arrays and are retained on
    // Host for clipping/alignment decisions; the numeric PairHMM kernel only
    // consumes key and flow_order.
    std::vector<std::int32_t> flow_to_base;
    std::vector<std::int32_t> reverse_key;
    std::vector<std::int32_t> reverse_flow_to_base;
};

EncodedFlowKey encode_flow_key(const std::vector<std::uint8_t>& bases,
                               const std::string& flow_order);

// Apply GATK LongHomopolymerHaplotypeCollapsingEngine.collapseBases().  The
// first homopolymer is preserved; subsequent runs are capped at
// hmer_size_threshold.  A threshold of zero is an explicit no-op.  This is a
// Host-side assembly transform performed before flow-key encoding.
std::vector<std::uint8_t> collapse_flow_homopolymers(
    const std::vector<std::uint8_t>& bases, std::size_t hmer_size_threshold);

// Equivalent to GATK LongHomopolymerHaplotypeCollapsingEngine.needsCollapsing:
// true only when a run is longer than the configured threshold.
bool flow_homopolymer_exceeds_threshold(
    const std::vector<std::uint8_t>& bases, std::size_t hmer_size_threshold);

// Return {number_of_flows_to_remove, bases_to_remove_from_first_remaining
// flow}, matching FlowBasedReadUtils.findLeftClipping/findRightClipping.
std::pair<std::size_t, std::size_t> find_left_flow_clipping(
    std::size_t base_clipping, const EncodedFlowKey& encoded);
std::pair<std::size_t, std::size_t> find_right_flow_clipping(
    std::size_t base_clipping, const EncodedFlowKey& encoded);

// Decode the production BAM representation described by GATK FlowBasedRead.
// This is deliberately a Host utility: it owns validation, base/quality
// conversion, the tp/t0 rules and boundary-flow policy before a flat table is
// copied into Kokkos Views.  Unsupported legacy kr/ti encodings fail closed at
// the caller instead of being silently interpreted as tp.
DecodedFlowRead decode_flow_read(const FlowReadTags& input);

// Apply GATK FlowBasedRead.applyBaseClipping semantics to an already decoded
// flow read. The operation is Host-side: base clipping is converted to
// boundary-flow removal, the boundary homopolymer count is shifted, and the
// probability columns are shifted in lockstep before the flat table enters
// the Kokkos flow PairHMM kernel.  Zero base clips are meaningful: GATK
// still runs applyClipping(0, 0, 0, 0), which strips terminal zero flows.
void clip_decoded_flow_read(DecodedFlowRead& decoded,
                            std::size_t left_bases,
                            std::size_t right_bases,
                            bool spread_boundary_probabilities = true);

// Convert one HTSlib-backed ReadBatch record into the raw Host tag view.  This
// is intentionally a copying adapter at the batch boundary; callers can then
// decode/calibrate once and hand only flat numeric arrays to Kokkos.
FlowReadTags flow_tags_from_batch(const ReadBatch& batch, std::size_t record,
                                  bool use_t0_tag = false,
                                  double filling_value = 0.001,
                                  bool keep_boundary_flows = false);

}  // namespace fastgatk::io
