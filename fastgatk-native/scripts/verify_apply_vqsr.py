#!/usr/bin/env python3
"""Contract checks for ApplyVQSR's scalar and allele-specific recal-VCF paths."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build")) / "fastgatk-apply-vqsr"
HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def semantic_records(path: pathlib.Path) -> list[tuple[str, ...]]:
    """Normalize VCF records so numeric INFO formatting is not an oracle false positive."""
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, ...]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, object] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    if separator:
                        if key in {"VQSLOD", "AS_VQSLOD"}:
                            try:
                                info[key] = float(value)
                            except ValueError:
                                info[key] = value
                        else:
                            info[key] = value
                    else:
                        info[key] = True
            records.append((fields[0], fields[1], fields[2], fields[3], fields[4],
                            fields[5], fields[6], json.dumps(info, sort_keys=True)))
    return records


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-") as temporary:
        work = pathlib.Path(temporary)
        input_vcf = work / "input.vcf"
        recal_vcf = work / "recal.vcf"
        tranches = work / "tranches.csv"
        output = work / "output.vcf.gz"
        manifest = work / "output.json"
        write_vcf(input_vcf, [
            "chr1\t1\t.\tA\tG\t50\tPASS\t.",
            "chr1\t2\t.\tC\tT\t50\tPASS\t.",
            "chr1\t3\t.\tA\tAC\t50\tPASS\t.",
        ])
        write_vcf(recal_vcf, [
            "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.0",
            "chr1\t2\t.\tC\tT\t.\tPASS\tVQSLOD=-1.0",
            "chr1\t3\t.\tA\tAC\t.\tPASS\tVQSLOD=-2.0",
        ])
        tranches.write_text(
            "# Variant quality score tranches file\n# Version number 5\n"
            "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
            "90.00,1,1,2.0,1.0,1.0000,VQSRTrancheSNP90.00to99.00,SNP,100,90,0.9000\n"
            "99.00,1,1,2.0,1.0,0.0000,VQSRTrancheSNP99.00to100.00,SNP,100,99,0.9900\n",
            encoding="utf-8",
        )
        result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--tranches-file", str(tranches), "--truth-sensitivity-filter-level", "99",
             "--mode", "SNP", "--use-allele-specific-annotations=false",
             "-O", str(output), "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=True,
        )
        text = gzip.open(output, "rt", encoding="utf-8").read()
        rows = [line.split("\t") for line in text.splitlines()
                if line and not line.startswith("#")]
        assert rows[0][6] == "PASS" and "VQSLOD=2" in rows[0][7], rows
        assert "VQSRTrancheSNP99.00to100.00" in rows[1][6], rows
        assert rows[1][6].endswith("+"), rows
        assert rows[2][6] == "PASS", rows
        assert (output.with_suffix(output.suffix + ".tbi")).is_file()
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        assert payload["tool"] == "ApplyVQSR"
        assert payload["scored_records"] == 2
        assert payload["filtered_records"] == 1
        assert payload["recal_allele_specific"] is False
        assert payload["cutoff_source"] == "tranches"
        assert payload["tranche_count"] == 1
        assert payload["provenance"]["recal_file_size"] == recal_vcf.stat().st_size
        assert payload["provenance"]["tranches_file_size"] == tranches.stat().st_size
        assert json.loads(result.stdout)["filtered_records"] == 1

        # The shared plain-VCF output boundary publishes a Tribble `.idx`
        # and must be usable by a real GATK interval query.
        plain_output = work / "output-plain.vcf"
        plain_manifest = work / "output-plain.json"
        plain_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "--mode", "SNP", "-O", str(plain_output),
             "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.is_file() and pathlib.Path(f"{plain_output}.idx").is_file()
        assert not pathlib.Path(f"{plain_output}.tbi").exists()
        java = ROOT / "third_party" / "jdk17" / "bin" / "java"
        gatk_jar = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
        plain_query = work / "output-plain-query.vcf"
        query = subprocess.run(
            [str(java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
             "-L", "chr1:1-1", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert query.returncode == 0, query.stderr
        assert len(semantic_records(plain_query)) == 1
        plain_payload = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_payload["create_output_variant_index"] is True

        # GATK's index switch is an optional boolean.  A false value must
        # suppress the Tabix sidecar while preserving identical record
        # semantics and making the choice auditable in the manifest.
        no_index_output = work / "no-index.output.vcf.gz"
        no_index_manifest = work / "no-index.output.json"
        no_index_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "--mode", "SNP",
             "--create-output-variant-index=false", "-O", str(no_index_output),
             "--output-manifest", str(no_index_manifest)],
            text=True, capture_output=True, check=True,
        )
        assert not no_index_output.with_suffix(no_index_output.suffix + ".tbi").exists()
        no_index_payload = json.loads(no_index_manifest.read_text(encoding="utf-8"))
        assert no_index_payload["create_output_variant_index"] is False
        assert json.loads(no_index_result.stdout)["create_output_variant_index"] is False
        indexed_cutoff_output = work / "indexed-cutoff.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "--mode", "SNP", "-O", str(indexed_cutoff_output)],
            text=True, capture_output=True, check=True,
        )
        assert indexed_cutoff_output.with_suffix(indexed_cutoff_output.suffix + ".tbi").is_file()
        assert semantic_records(no_index_output) == semantic_records(indexed_cutoff_output)
        invalid_index = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "--create-output-variant-index=maybe",
             "-O", str(work / "invalid-index.vcf.gz")],
            text=True, capture_output=True,
        )
        assert invalid_index.returncode != 0
        assert "invalid boolean" in invalid_index.stderr

        # Repeated GATK interval selectors are evaluated before mode scoring;
        # INTERSECTION of chr1:1-3 and chr1:2-2 therefore keeps only POS=2.
        interval_output = work / "interval.output.vcf.gz"
        interval_manifest = work / "interval.output.json"
        interval_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--lod-score-cutoff", "0", "--mode", "SNP", "-L", "chr1:1-3",
             "-L", "chr1:2-2", "--interval-set-rule", "INTERSECTION",
             "-O", str(interval_output), "--output-manifest", str(interval_manifest)],
            text=True, capture_output=True, check=True,
        )
        interval_rows = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8") if line and not line.startswith("#")]
        assert len(interval_rows) == 1 and interval_rows[0][1] == "2", interval_rows
        interval_payload = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_payload["interval_set_rule"] == "INTERSECTION"
        assert interval_payload["interval_skipped_records"] == 2
        assert json.loads(interval_result.stdout)["interval_skipped_records"] == 2

        compressed_tranches = work / "tranches.csv.gz"
        with gzip.open(compressed_tranches, "wt", encoding="utf-8") as stream:
            stream.write(tranches.read_text(encoding="utf-8"))
        compressed_output = work / "compressed-tranches.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "--tranches-file", str(compressed_tranches), "--truth-sensitivity-filter-level", "99",
             "--mode", "SNP", "-O", str(compressed_output)],
            text=True, capture_output=True, check=True,
        )
        assert gzip.open(compressed_output, "rt", encoding="utf-8").read() == text

        # Scalar records use the same GATK mode classification as AS records:
        # MNPs are part of SNP mode and symbolic structural records are part
        # of INDEL mode.
        mnp_scalar_input = work / "mnp.scalar.input.vcf"
        mnp_scalar_recal = work / "mnp.scalar.recal.vcf"
        write_vcf(mnp_scalar_input, ["chr1\t10\t.\tAA\tGG\t50\tPASS\t."])
        write_vcf(mnp_scalar_recal, ["chr1\t10\t.\tAA\tGG\t.\tPASS\tVQSLOD=2.0"])
        mnp_scalar_output = work / "mnp.scalar.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(mnp_scalar_input), "--recal-file", str(mnp_scalar_recal),
             "--lod-score-cutoff", "0", "--mode", "SNP", "-O", str(mnp_scalar_output)],
            text=True, capture_output=True, check=True,
        )
        mnp_scalar_row = next(line.split("\t") for line in gzip.open(mnp_scalar_output, "rt", encoding="utf-8")
                              if line and not line.startswith("#"))
        assert "VQSLOD=2" in mnp_scalar_row[7]

        symbolic_scalar_input = work / "symbolic.scalar.input.vcf"
        symbolic_scalar_recal = work / "symbolic.scalar.recal.vcf"
        write_vcf(symbolic_scalar_input, ["chr1\t11\t.\tA\t<DEL>\t50\tPASS\t."])
        write_vcf(symbolic_scalar_recal, ["chr1\t11\t.\tA\t<DEL>\t.\tPASS\tVQSLOD=2.0"])
        symbolic_scalar_output = work / "symbolic.scalar.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(symbolic_scalar_input), "--recal-file", str(symbolic_scalar_recal),
             "--lod-score-cutoff", "0", "--mode", "INDEL", "-O", str(symbolic_scalar_output)],
            text=True, capture_output=True, check=True,
        )
        symbolic_scalar_row = next(line.split("\t") for line in gzip.open(symbolic_scalar_output, "rt", encoding="utf-8")
                                   if line and not line.startswith("#"))
        assert "VQSLOD=2" in symbolic_scalar_row[7]

        as_recal_vcf = work / "as.recal.vcf"
        as_recal_vcf.write_text(
            HEADER.replace(
                "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n",
                "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
                "##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=allele-score>\n"
                "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=allele-status>\n"
                "##INFO=<ID=AS_culprit,Number=A,Type=String,Description=allele-culprit>\n",
            ) + "\n".join([
                "chr1\t1\t.\tA\tG\t.\tPASS\tAS_VQSLOD=2.0;AS_culprit=QD",
                "chr1\t2\t.\tC\tT\t.\tPASS\tAS_VQSLOD=-1.0;AS_culprit=MQ",
                "chr1\t3\t.\tA\tAC\t.\tPASS\tAS_VQSLOD=-2.0;AS_culprit=FS",
            ]) + "\n",
            encoding="utf-8",
        )
        as_output = work / "as.output.vcf.gz"
        as_manifest = work / "as.output.json"
        as_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(as_recal_vcf),
             "--lod-score-cutoff", "0", "--AS", "--mode", "SNP", "-O", str(as_output),
             "--output-manifest", str(as_manifest)],
            text=True, capture_output=True, check=True,
        )
        as_text = gzip.open(as_output, "rt", encoding="utf-8").read()
        as_rows = [line.split("\t") for line in as_text.splitlines()
                   if line and not line.startswith("#")]
        assert "AS_VQSLOD=2" in as_rows[0][7] and "AS_FilterStatus=PASS" in as_rows[0][7]
        assert "AS_culprit=QD" in as_rows[0][7]
        assert "AS_FilterStatus=LOW_VQSLOD" in as_rows[1][7]
        as_payload = json.loads(as_manifest.read_text(encoding="utf-8"))
        assert as_payload["allele_specific"] is True
        assert as_payload["applicable_alleles"] == 2
        assert as_payload["mode_skipped_alleles"] == 1
        assert json.loads(as_result.stdout)["allele_specific"] is True

        # GATK's scalar recalibration VCF uses a dummy N/<VQSR> allele and
        # carries POSITIVE/NEGATIVE_TRAIN_SITE flags.  ApplyVQSR matches these
        # records by location (including END), propagates the labels, and
        # classifies scores using the selected tranche interval rather than a
        # single fixed cutoff filter.
        dummy_recal = work / "dummy.recal.vcf"
        dummy_recal.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=end>\n"
            "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
            "##INFO=<ID=CULPRIT,Number=1,Type=String,Description=culprit>\n"
            "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=positive>\n"
            "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=negative>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tN\t<VQSR>\t.\tPASS\tEND=1;VQSLOD=2.0;CULPRIT=QD;POSITIVE_TRAIN_SITE\n"
            "chr1\t2\t.\tN\t<VQSR>\t.\tPASS\tEND=2;VQSLOD=-1.0;CULPRIT=MQ;NEGATIVE_TRAIN_SITE\n"
            "chr1\t3\t.\tN\t<VQSR>\t.\tPASS\tEND=3;VQSLOD=-2.0;CULPRIT=FS\n",
            encoding="utf-8",
        )
        dummy_output = work / "dummy.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(dummy_recal),
             "--tranches-file", str(tranches), "--truth-sensitivity-filter-level", "99",
             "--mode", "SNP", "-O", str(dummy_output)],
            text=True, capture_output=True, check=True,
        )
        dummy_rows = [line.split("\t") for line in gzip.open(dummy_output, "rt", encoding="utf-8")
                      if line and not line.startswith("#")]
        assert dummy_rows[0][6] == "PASS" and "VQSLOD=2" in dummy_rows[0][7]
        assert "POSITIVE_TRAIN_SITE" in dummy_rows[0][7]
        assert dummy_rows[1][6].endswith("+") and "NEGATIVE_TRAIN_SITE" in dummy_rows[1][7]
        # The third input record is an INDEL and is untouched in SNP mode.
        assert "VQSLOD" not in dummy_rows[2][7] and dummy_rows[2][6] == "PASS"

        filtered_input = work / "filtered.input.vcf"
        filtered_input.write_text(
            HEADER.replace("#CHROM", "##FILTER=<ID=LowQual,Description=low quality>\n#CHROM")
            + "chr1\t1\t.\tA\tG\t50\tLowQual\t.\n", encoding="utf-8")
        filtered_recal = work / "filtered.recal.vcf"
        filtered_recal.write_text(
            HEADER + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.0\n", encoding="utf-8")
        untouched = work / "filtered.untouched.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(filtered_input), "--recal-file", str(filtered_recal),
             "--lod-score-cutoff", "0", "-O", str(untouched)],
            text=True, capture_output=True, check=True,
        )
        untouched_row = next(line.split("\t") for line in gzip.open(untouched, "rt", encoding="utf-8")
                             if line and not line.startswith("#"))
        assert untouched_row[6] == "LowQual" and "VQSLOD" not in untouched_row[7]
        ignored = work / "filtered.ignored.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(filtered_input), "--recal-file", str(filtered_recal),
             "--ignore-filter", "LowQual", "--lod-score-cutoff", "0", "-O", str(ignored)],
            text=True, capture_output=True, check=True,
        )
        ignored_row = next(line.split("\t") for line in gzip.open(ignored, "rt", encoding="utf-8")
                           if line and not line.startswith("#"))
        assert ignored_row[6] == "PASS" and "VQSLOD=2" in ignored_row[7]

        mixed_input = work / "mixed.input.vcf"
        mixed_input.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t4\t.\tA\tG,AC\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        mixed_recal = work / "mixed.recal.vcf"
        mixed_recal.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t4\t.\tA\tG\t.\tPASS\tAS_VQSLOD=2.0;AS_culprit=QD\n",
            encoding="utf-8",
        )
        mixed_output = work / "mixed.output.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(mixed_input), "--recal-file", str(mixed_recal),
             "--lod-score-cutoff", "0", "--AS", "--mode", "SNP", "-O", str(mixed_output)],
            text=True, capture_output=True, check=True,
        )
        mixed_row = next(line.split("\t") for line in gzip.open(mixed_output, "rt", encoding="utf-8")
                         if line and not line.startswith("#"))
        assert "AS_FilterStatus=PASS,NA" in mixed_row[7], mixed_row
        assert "AS_culprit=QD,NA" in mixed_row[7], mixed_row

        # Allele-based mode classification follows GATK: symbolic structural
        # alleles are INDELs, while a spanning deletion '*' is SNP-shaped
        # because REF and ALT have equal length.
        symbolic_input = work / "symbolic.input.vcf"
        symbolic_input.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t5\t.\tA\t<DEL>\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        symbolic_recal = work / "symbolic.recal.vcf"
        symbolic_recal.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t5\t.\tA\t<DEL>\t.\tPASS\tAS_VQSLOD=2.0;AS_culprit=SV\n",
            encoding="utf-8",
        )
        symbolic_output = work / "symbolic.output.vcf.gz"
        symbolic_result = subprocess.run(
            [str(BINARY), "-V", str(symbolic_input), "--recal-file", str(symbolic_recal),
             "--lod-score-cutoff", "0", "--AS", "--mode", "INDEL", "-O", str(symbolic_output)],
            text=True, capture_output=True, check=True,
        )
        symbolic_row = next(line.split("\t") for line in gzip.open(symbolic_output, "rt", encoding="utf-8")
                            if line and not line.startswith("#"))
        assert "AS_VQSLOD=2" in symbolic_row[7]
        assert json.loads(symbolic_result.stdout)["applicable_alleles"] == 1

        spanning_input = work / "spanning.input.vcf"
        spanning_input.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t6\t.\tA\t*\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        spanning_recal = work / "spanning.recal.vcf"
        spanning_recal.write_text(
            as_recal_vcf.read_text(encoding="utf-8").split("#CHROM")[0]
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t6\t.\tA\t*\t.\tPASS\tAS_VQSLOD=2.0;AS_culprit=SPAN\n",
            encoding="utf-8",
        )
        spanning_output = work / "spanning.output.vcf.gz"
        spanning_result = subprocess.run(
            [str(BINARY), "-V", str(spanning_input), "--recal-file", str(spanning_recal),
             "--lod-score-cutoff", "0", "--AS", "--mode", "SNP", "-O", str(spanning_output)],
            text=True, capture_output=True, check=True,
        )
        spanning_row = next(line.split("\t") for line in gzip.open(spanning_output, "rt", encoding="utf-8")
                            if line and not line.startswith("#"))
        assert "AS_VQSLOD=2" in spanning_row[7]
        assert json.loads(spanning_result.stdout)["applicable_alleles"] == 1

        unsupported = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
             "-O", str(work / "unsupported.vcf"), "--use-allele-specific-annotations"],
            text=True, capture_output=True,
        )
        assert unsupported.returncode != 0
        assert "BACKEND_UNAVAILABLE" in unsupported.stderr
        java = ROOT / "third_party" / "jdk17" / "bin" / "java"
        gatk = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
        java_oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
        java_as_oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
        if oracle_guard.oracle_ready('verify_apply_vqsr.py', java, gatk):
            # ApplyVQSR queries the recalibration VCF by interval, so make the
            # exact same small fixture random-accessible before invoking Java.
            subprocess.run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_vcf)],
                           text=True, capture_output=True, check=True)
            subprocess.run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(recal_vcf)],
                           text=True, capture_output=True, check=True)
            java_output = work / "java.output.vcf.gz"
            subprocess.run(
                [str(java), "-jar", str(gatk), "ApplyVQSR", "-V", str(input_vcf),
                 "--recal-file", str(recal_vcf), "--tranches-file", str(tranches),
                 "--truth-sensitivity-filter-level", "99", "--mode", "SNP", "-O", str(java_output),
                 "--create-output-variant-index", "false"],
                text=True, capture_output=True, check=True,
            )
            native_oracle_output = work / "native.oracle.output.vcf.gz"
            subprocess.run(
                [str(BINARY), "-V", str(input_vcf), "--recal-file", str(recal_vcf),
                 "--tranches-file", str(tranches), "--truth-sensitivity-filter-level", "99",
                 "--mode", "SNP", "-O", str(native_oracle_output)],
                text=True, capture_output=True, check=True,
            )
            assert semantic_records(java_output) == semantic_records(native_oracle_output)
            java_oracle = {
                "status": "pass",
                "records_compared": len(semantic_records(native_oracle_output)),
                "semantic_records_bit_identical": True,
            }

            # Keep the AS oracle deliberately single-ALT.  GATK's AS ApplyVQSR
            # requires the recalibration record's allele set to exactly match
            # the input record; the mixed-allele fixture above is intentionally
            # testing native NA padding and is therefore not a valid Java
            # oracle input.  A dot culprit is used because Java emits that
            # provenance value for this minimal hand-authored recal VCF while
            # still exercising AS score/status classification.
            as_java_input = work / "as.java.input.vcf"
            as_java_recal = work / "as.java.recal.vcf"
            as_java_header = (
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=1000>\n"
                "##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=allele-score>\n"
                "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=allele-status>\n"
                "##INFO=<ID=AS_culprit,Number=A,Type=String,Description=allele-culprit>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )
            as_java_input.write_text(
                as_java_header
                + "chr1\t1\t.\tA\tG\t50\tPASS\t.\n"
                + "chr1\t2\t.\tC\tT\t50\tPASS\t.\n",
                encoding="utf-8",
            )
            as_java_recal.write_text(
                as_java_header.replace(
                    "##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=allele-score>\n",
                    "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>\n"
                    "##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=allele-score>\n",
                )
                + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.0;AS_VQSLOD=2.0;AS_culprit=.\n"
                + "chr1\t2\t.\tC\tT\t.\tPASS\tVQSLOD=-1.0;AS_VQSLOD=-1.0;AS_culprit=.\n",
                encoding="utf-8",
            )
            subprocess.run(
                [str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(as_java_input)],
                text=True, capture_output=True, check=True,
            )
            subprocess.run(
                [str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(as_java_recal)],
                text=True, capture_output=True, check=True,
            )
            as_java_output = work / "as.java.output.vcf.gz"
            subprocess.run(
                [str(java), "-jar", str(gatk), "ApplyVQSR", "-V", str(as_java_input),
                 "--recal-file", str(as_java_recal),
                 "--use-allele-specific-annotations", "true", "--lod-score-cutoff", "0",
                 "--mode", "SNP", "-O", str(as_java_output),
                 "--create-output-variant-index", "false"],
                text=True, capture_output=True, check=True,
            )
            as_native_output = work / "as.native.output.vcf.gz"
            subprocess.run(
                [str(BINARY), "-V", str(as_java_input), "--recal-file", str(as_java_recal),
                 "--use-allele-specific-annotations", "true", "--lod-score-cutoff", "0", "--mode", "SNP",
                 "-O", str(as_native_output)],
                text=True, capture_output=True, check=True,
            )
            assert semantic_records(as_java_output) == semantic_records(as_native_output)
            java_as_oracle = {
                "status": "pass",
                "records_compared": len(semantic_records(as_native_output)),
                "semantic_records_bit_identical": True,
            }
        print(json.dumps({
            "status": "pass",
            "input_records": 3,
            "scored_records": 2,
            "filtered_records": 1,
            "indexed": True,
            "allele_specific_applied": True,
            "missing_allele_specific_rejected": True,
            "java_oracle": java_oracle,
            "java_as_oracle": java_as_oracle,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
