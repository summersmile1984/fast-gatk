#!/usr/bin/env python3
"""Verify opt-in HC contig streaming against the aggregate HostBatch path."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


SAM = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:40
@SQ\tSN:chr2\tLN:40
@RG\tID:rg1\tSM:streaming
q1\t0\tchr1\t1\t60\t21M\t*\t0\t0\tCGTATCGGCTAGCTTACGATC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
q2\t0\tchr1\t1\t60\t21M\t*\t0\t0\tCGTATCGGCTCGCTTACGATC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
q3\t0\tchr1\t1\t60\t21M\t*\t0\t0\tCGTATCGGCTCGCTTACGATC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
q4\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCGATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
q5\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCTATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
q6\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCTATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1
"""


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True, check=False)


def data_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-streaming-") as directory:
        work = Path(directory)
        sam = work / "two-contig.sam"
        reference = work / "reference.fa"
        sam.write_text(SAM, encoding="utf-8")
        reference.write_text(
            ">chr1\nCGTATCGGCTAGCTTACGATCGGTACCTGATCAGTCCGTA\n"
            ">chr2\nTGCACTGATCGATGCTAGTCCGATTCGAGCTTAGCATGCA\n",
            encoding="utf-8")

        common = ["-I", str(sam), "-R", str(reference), "--threads", "2",
                  "--min-depth", "1", "--min-alt-support", "1",
                  "--standard-min-confidence-threshold-for-calling", "0",
                  # Short, non-repetitive contigs ensure these records use
                  # the normal assembly EventMap path instead of a removed
                  # pileup/CIGAR fallback.
                  ]
        normal = work / "normal.vcf"
        streamed = work / "streamed.vcf"
        normal_manifest = work / "normal.manifest.json"
        streamed_manifest = work / "streamed.manifest.json"
        normal_result = run(binary, [*common, "-O", str(normal),
                                     "--output-manifest", str(normal_manifest)])
        assert normal_result.returncode == 0, normal_result.stderr
        streamed_result = run(binary, [*common, "-O", str(streamed), "--stream-by-contig",
                                       "--output-manifest", str(streamed_manifest)])
        assert streamed_result.returncode == 0, streamed_result.stderr
        assert normal.read_text(encoding="utf-8") == streamed.read_text(encoding="utf-8")
        records = data_records(streamed)
        assert len(records) == 2 and [row.split("\t")[0] for row in records] == ["chr1", "chr2"]

        normal_metadata = json.loads(normal_manifest.read_text(encoding="utf-8"))
        streamed_metadata = json.loads(streamed_manifest.read_text(encoding="utf-8"))
        normal_telemetry = normal_metadata["telemetry"]
        stream_telemetry = streamed_metadata["telemetry"]
        assert normal_telemetry["stream_by_contig"] is False
        assert streamed_metadata["compatibility"]["stream_by_contig"] is True
        assert stream_telemetry["streamed_contigs"] == 2
        assert stream_telemetry["streamed_peak_host_bytes"] > 0
        assert stream_telemetry["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert stream_telemetry["pipeline_decoded_items"] == stream_telemetry["pipeline_computed_items"]
        assert stream_telemetry["pipeline_computed_items"] == stream_telemetry["pipeline_encoded_items"]
        assert stream_telemetry["pipeline_decoded_items"] == 2
        assert stream_telemetry["pipeline_peak_decoded_bytes"] > 0
        assert stream_telemetry["pipeline_peak_computed_bytes"] > 0
        assert stream_telemetry["pipeline_peak_encoded_bytes"] > 0
        assert stream_telemetry["reads"] == normal_telemetry["reads"] == 6
        assert stream_telemetry["variant_calls"] == normal_telemetry["variant_calls"] == 2

        # Re-entry is rejected rather than silently duplicating a contig's
        # reference-confidence state.  The aggregate path remains available
        # for such non-coordinate-sorted inputs.
        unsorted = work / "unsorted.sam"
        unsorted.write_text(SAM.replace(
            "q5\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCTATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n"
            "q6\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCTATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1",
            "q5\t0\tchr1\t1\t60\t21M\t*\t0\t0\tCGTATCGGCTAGCTTACGATC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n"
            "q6\t0\tchr2\t1\t60\t21M\t*\t0\t0\tTGCACTGATCTATGCTAGTCC\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1"),
            encoding="utf-8")
        rejected = run(binary, ["-I", str(unsorted), "-R", str(reference),
                                "-O", str(work / "rejected.vcf"), "--stream-by-contig",
                                "--min-depth", "1", "--min-alt-support", "1"])
        assert rejected.returncode != 0
        assert "stream-by-contig" in rejected.stderr

    print(json.dumps({"status": "pass", "streamed_contigs": 2,
                      "records": 2, "aggregate_output_exact": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
