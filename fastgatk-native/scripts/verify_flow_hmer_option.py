#!/usr/bin/env python3
"""Regression for the Host flow haplotype HMER-restoration option."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


# This has a unique-kmer reference flank around a 12-base HMER.  The former
# all-A synthetic reference was intentionally useful for the codec, but it is
# rejected by ReadThreadingAssembler's non-unique-reference-kmer retry and
# cannot prove the caller's Flow PairHMM ordering.
FLOW_REFERENCE = (
    "GCTAGTCGATCGTACGATCGGATCCTAGCGTACCTGATCGTACG"
    "AAAAAAAAAAAA"
    "CGTACGATGCTAGGCTTACGATCGTCCGATGCTAGCTACGATCG")
FLOW_READ_START = 21  # 1-based; the 41M read spans the long A HMER.
FLOW_READ_BASES = (
    FLOW_REFERENCE[FLOW_READ_START - 1:FLOW_READ_START - 1 + 5] + "A" +
    FLOW_REFERENCE[FLOW_READ_START - 1 + 6:FLOW_READ_START - 1 + 41])


def flow_sam(read_count: int) -> str:
    tp = ",".join("0" for _ in FLOW_READ_BASES)
    qualities = "I" * len(FLOW_READ_BASES)
    indel_qualities = "!" * len(FLOW_READ_BASES)
    t0 = "K" * len(FLOW_READ_BASES)
    records = "\n".join(
        f"flow-read-{index}\t0\tchr1\t{FLOW_READ_START}\t60\t41M\t*\t0\t0\t"
        f"{FLOW_READ_BASES}\t{qualities}\tRG:Z:flow\tBI:Z:{indel_qualities}\t"
        f"BD:Z:{indel_qualities}\ttp:B:c,{tp}\tt0:Z:{t0}"
        for index in range(read_count))
    return "\n".join((
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:chr1\tLN:{len(FLOW_REFERENCE)}",
        "@RG\tID:flow\tPL:ULTIMA\tFO:ACGT\tmc:12",
        records,
        ""))


def run(binary: Path, root: Path, work: Path, value: int,
        reference: Path | None = None, partial_mode: bool = False,
        input_path: Path | None = None, native_double_precision: bool = True) -> dict:
    manifest = work / f"manifest-{value}.json"
    output = work / f"output-{value}.json"
    command = [
        str(binary),
        "-I", str(input_path or root / "fastgatk-native/tests/flow_reader_fixture.sam"),
        "-R", str(reference or root / "fastgatk-native/tests/flow_reader_reference.fa"),
        "-L", f"chr1:1-{len(FLOW_REFERENCE)}",
        "-O", str(output),
        "--output-manifest", str(manifest),
        "--min-depth", "1",
        "--min-alt-support", "1",
        "--flow-assembly-collapse-hmer-size", str(value),
    ]
    if not native_double_precision:
        command.extend(["--native-pair-hmm-use-double-precision", "false"])
    if partial_mode:
        command.append("--flow-assembly-collapse-partial-mode")
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    if not binary.exists():
        raise SystemExit("missing HC binary")
    with tempfile.TemporaryDirectory(prefix="fastgatk-flow-hmer-") as directory:
        work = Path(directory)
        # The pinned reader fixture is a parser probe.  Build a valid,
        # reference-unique long-HMER/SNP fixture here so the regression reaches
        # the GATK-compatible assembly -> Flow PairHMM option boundary.
        flow_reference = work / "flow-reference.fa"
        flow_remap_input = work / "flow-remap.sam"
        flow_remap_input.write_text(flow_sam(4), encoding="ascii")
        flow_reference.write_text(">chr1\n" + FLOW_REFERENCE + "\n", encoding="ascii")
        # Default ReadThreadingAssembler pruning requires two supporting
        # observations; use four identical flow reads for every caller path.
        disabled = run(binary, root, work, 0, flow_reference, input_path=flow_remap_input)
        enabled = run(binary, root, work, 2, flow_reference, input_path=flow_remap_input)
        partial = run(binary, root, work, 2, flow_reference, partial_mode=True,
                      input_path=flow_remap_input)
        auto = run(binary, root, work, -1, flow_reference, input_path=flow_remap_input)
        float_requested = run(binary, root, work, 2, flow_reference,
                              input_path=flow_remap_input,
                              native_double_precision=False)
        short_reference = work / "short-reference.fa"
        short_reference.write_text(
            ">chr1\n" + ("ACGT" * 30) + "\n",
            encoding="ascii")
        short_reference_result = run(binary, root, work, 2, short_reference,
                                     input_path=flow_remap_input)
    if disabled["pairhmm_error_model"] != "flow":
        raise AssertionError(
            "disabled flow fixture did not reach Flow PairHMM: " +
            json.dumps(disabled, sort_keys=True))
    assert disabled["pairhmm_flow_haplotypes_collapsed"] == 0
    assert enabled["pairhmm_error_model"] == "flow"
    # GATK's FlowBasedAlignmentLikelihoodEngine never reads the native
    # base-space PairHMM precision argument.  The native Host adapter must
    # therefore preserve Flow Kokkos dispatch when the common CLI spelling
    # explicitly supplies false.
    assert float_requested["pairhmm_error_model"] == "flow"
    assert float_requested["pairhmm_skip_reason"] == "executed-flow"
    assert enabled["flow_assembly_collapse_hmer_size"] == 2
    # GATK scores FlowBasedHaplotype's assembled bases first, then restores
    # bounded HMERs before allele marginalization.  A pre-PairHMM collapse
    # would change the likelihood matrix, so the compatibility counter remains
    # zero.  The direct deleted-HMER restoration recurrence is covered by
    # fastgatk-calling-indel-smoke.
    assert enabled["pairhmm_flow_haplotypes_collapsed"] == 0
    for telemetry in (disabled, enabled, partial, auto, float_requested,
                      short_reference_result):
        assert telemetry["pairhmm_flow_haplotype_remaps"] >= 0
        assert telemetry["pairhmm_flow_identical_haplotype_groups"] >= 0
        assert telemetry["pairhmm_flow_haplotype_remaps"] >= telemetry[
            "pairhmm_flow_identical_haplotype_groups"]
    assert partial["flow_assembly_collapse_partial_mode"] is True
    assert auto["flow_assembly_collapse_hmer_size"] == -1
    assert auto["pairhmm_flow_haplotypes_collapsed"] == 0
    assert short_reference_result["pairhmm_flow_haplotypes_collapsed"] == 0
    print(json.dumps({
        "status": "pass",
        "disabled_collapsed": disabled["pairhmm_flow_haplotypes_collapsed"],
        "pre_score_collapsed": enabled["pairhmm_flow_haplotypes_collapsed"],
        "enabled_uncollapsed": enabled["pairhmm_flow_haplotypes_uncollapsed"],
        "enabled_remaps": enabled["pairhmm_flow_haplotype_remaps"],
        "enabled_identical_groups": enabled[
            "pairhmm_flow_identical_haplotype_groups"],
        "auto_uncollapsed": auto["pairhmm_flow_haplotypes_uncollapsed"],
        "short_reference_uncollapsed": short_reference_result[
            "pairhmm_flow_haplotypes_uncollapsed"],
        "flow_haplotypes": enabled["pairhmm_flow_haplotypes"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
