#!/usr/bin/env python3
"""Verify SelectVariants' GATK sites-only writer boundary.

GATK 4.6.2.0 keeps the complete FORMAT/sample payload available while
selecting records, then drops those columns only at the final VCF writer when
``--sites-only-vcf-output`` is enabled.  The option is an optional Boolean
whose bare and separated forms are accepted, while Barclay rejects an
embedded ``=`` spelling.  This oracle fixes all three behaviors against a
small FORMAT-rich VCF.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FILTER=<ID=q10,Description=Low quality>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##INFO=<ID=AS_INT,Number=A,Type=Integer,Description=Per-ALT score>
##INFO=<ID=ADINFO,Number=R,Type=Integer,Description=Per-allele score>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t1\trs1\tA\tG,T\t50\tPASS\tDP=20;AS_INT=11,22;ADINFO=100,11,22\tGT:AD:PL:GQ\t0/1:10,10,0:50,0,50,99,99,99:50\t0/2:8,0,8:80,99,80,99,99,0:1
chr1\t2\trs2\tC\tT\t10\tq10\tDP=5;AS_INT=7;ADINFO=50,3\tGT:AD:PL:GQ\t0/1:2,3:20,0,20:20\t0/0:5,0:0,20,80:80
"""


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=True)


def records(path: Path) -> list[tuple[str, ...]]:
    rows: list[tuple[str, ...]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                # HTSJDK and HTSlib may serialize INFO keys in different
                # orders; compare the site payload while preserving every
                # key/value and all fixed VCF columns.
                if len(fields) >= 8 and fields[7] not in ("", "."):
                    fields[7] = ";".join(sorted(fields[7].split(";")))
                rows.append(tuple(fields))
    return rows


def chrom_header(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")
    raise AssertionError(f"missing #CHROM header in {path}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SELECT_VARIANTS_BINARY",
        root / "fastgatk-native/build/fastgatk-select-variants"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"SelectVariants sites-only oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip", "reason": "GATK sites-only oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-select-variants-sites-only-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")

        gatk_outputs: dict[str, Path] = {}
        native_outputs: dict[str, Path] = {}
        forms = {
            "bare": ["--sites-only-vcf-output"],
            "true": ["--sites-only-vcf-output", "true"],
            "false": ["--sites-only-vcf-output", "false"],
        }
        for label, switch in forms.items():
            output = work / f"gatk-{label}.vcf"
            result = run([str(java), "-jar", str(jar), "SelectVariants", "-V", str(source),
                          *switch, "-O", str(output)])
            if result.returncode != 0:
                raise RuntimeError(f"GATK SelectVariants {label} failed: {result.stderr[-3000:]}")
            gatk_outputs[label] = output

            native_output = work / f"native-{label}.vcf"
            manifest = work / f"native-{label}.manifest.json"
            result = run([str(binary), "-V", str(source), *switch,
                          "--create-output-variant-index", "false", "-O", str(native_output),
                          "--output-manifest", str(manifest)])
            if result.returncode != 0:
                raise RuntimeError(f"native SelectVariants {label} failed: {result.stderr[-3000:]}")
            native_outputs[label] = native_output
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if metadata["compatibility"]["sites_only_vcf_output"] is not (label != "false"):
                raise AssertionError(f"native {label} manifest lost sites-only choice")
            if metadata["telemetry"]["sites_only_vcf_output"] is not (label != "false"):
                raise AssertionError(f"native {label} telemetry lost sites-only choice")

        for label in ("bare", "true"):
            if len(chrom_header(gatk_outputs[label])) != 8:
                raise AssertionError(f"GATK {label} did not emit an 8-column sites-only header")
            if len(chrom_header(native_outputs[label])) != 8:
                raise AssertionError(f"native {label} did not emit an 8-column sites-only header")
            if any(len(row) != 8 for row in records(gatk_outputs[label])):
                raise AssertionError(f"GATK {label} emitted FORMAT fields")
            if any(len(row) != 8 for row in records(native_outputs[label])):
                raise AssertionError(f"native {label} emitted FORMAT fields")
            if records(native_outputs[label]) != records(gatk_outputs[label]):
                raise AssertionError(f"native {label} site rows differ from GATK")
            if records(native_outputs[label]) != records(gatk_outputs["true"]):
                raise AssertionError(f"GATK bare/true sites-only rows differ")

        if len(chrom_header(gatk_outputs["false"])) != 11:
            raise AssertionError("GATK false unexpectedly removed sample columns")
        if len(chrom_header(native_outputs["false"])) != 11:
            raise AssertionError("native false unexpectedly removed sample columns")
        if records(native_outputs["false"]) != records(gatk_outputs["false"]):
            raise AssertionError("native false full records differ from GATK")

        # Barclay treats this argument name as a complete token and rejects an
        # embedded '=' rather than binding an inline Boolean value.
        embedded_output = work / "native-embedded.vcf"
        embedded = run([str(binary), "-V", str(source),
                        "--sites-only-vcf-output=true", "--create-output-variant-index", "false",
                        "-O", str(embedded_output)], check=False)
        if embedded.returncode == 0:
            raise AssertionError("native accepted embedded sites-only Boolean")
        gatk_embedded = run([str(java), "-jar", str(jar), "SelectVariants", "-V", str(source),
                             "--sites-only-vcf-output=true", "-O", str(work / "gatk-embedded.vcf")],
                            check=False)
        if gatk_embedded.returncode == 0:
            raise AssertionError("GATK unexpectedly accepted embedded sites-only Boolean")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "sites_only_forms_exact": ["bare", "true"],
            "false_full_shape_exact": True,
            "site_records_compared": len(records(native_outputs["true"])),
            "sites_only_columns": 8,
            "embedded_equals_fail_closed": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
