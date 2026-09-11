#!/usr/bin/env python3
"""Contract checks for the Kokkos VQSR Gaussian/GMM prototype."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_VARIANT_RECALIBRATOR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-variant-recalibrator")))
APPLY_BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_APPLY_VQSR_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-apply-vqsr")))
HEADER = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "##FILTER=<ID=LowQual,Description=low quality>\n"
    "##INFO=<ID=QD,Number=1,Type=Float,Description=QD>\n"
    "##INFO=<ID=MQ,Number=1,Type=Float,Description=MQ>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)
JITTER_HEADER = HEADER.replace(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    "##INFO=<ID=FS,Number=1,Type=Float,Description=FS>\n"
    "##INFO=<ID=SOR,Number=1,Type=Float,Description=SOR>\n"
    "##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=HaplotypeScore>\n"
    "##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=InbreedingCoeff>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
)


def write_vcf(path: pathlib.Path, rows: list[str]) -> None:
    path.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")


def write_vcf_with_header(path: pathlib.Path, header: str, rows: list[str]) -> None:
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    gatk_java = ROOT / "third_party" / "jdk17" / "bin" / "java"
    gatk_jar = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-recalibrator-") as temporary:
        work = pathlib.Path(temporary)
        input_vcf = work / "input.vcf"
        training_vcf = work / "training.vcf"
        known_vcf = work / "known.vcf"
        output = work / "recal.vcf.gz"
        tranches = work / "recal.tranches"
        rscript = work / "plot.R"
        manifest = work / "recal.json"
        model = work / "recal.model"
        rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ={mq}"
            for index, (qd, mq) in enumerate([(30, 60), (28, 58), (25, 55), (4, 25), (3, 20), (2, 15)], 1)
        ]
        write_vcf(input_vcf, rows)
        write_vcf(training_vcf, rows[:3])
        write_vcf(known_vcf, rows[:2])
        result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
             "--max-attempts", "20", "--k-means-iterations", "20", "-O", str(output),
             "--tranches-file", str(tranches), "--truth-sensitivity-tranche", "100",
             "--truth-sensitivity-tranche", "90", "--rscript-file", str(rscript),
             "--output-manifest", str(manifest), "--output-model", str(model), "--threads", "2"],
            text=True, capture_output=True, check=True,
        )
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line for line in text.splitlines() if line and not line.startswith("#")]
        assert len(records) == 6, records
        # GATK writes the NAME of the datum's worst-scoring annotation
        # dimension here, never a model/provenance string:
        # VariantRecalibratorEngine.calculateWorstPerformingAnnotation
        # (VariantRecalibratorEngine.java:80-93) minimises
        # goodModel.evaluateDatumInOneDimension - badModel.evaluateDatumInOneDimension
        # (the per-dimension log10 mixture likelihood,
        # GaussianMixtureModel.java:208-222) and VariantDataManager.java:485
        # writes annotationKeys.get(worstAnnotation).  The previous value,
        # "full-covariance-gmm", was a native-only provenance string (0
        # occurrences in the pinned gatk-package-4.6.2.0-local.jar); pinned GATK
        # 4.6.2.0 emits culprit=MQ for all six records of this fixture
        # (byte-identical gate:
        # verify_variant_recalibrator_culprit_gatk_oracle.py, case
        # tiny-truth-known).
        assert all("VQSLOD=" in line and "culprit=MQ" in line for line in records), records
        assert all("\tN\t<VQSR>\t.\t.\t" in line for line in records), records
        assert "POSITIVE_TRAIN_SITE" in records[0].split("\t")[7]
        assert "POSITIVE_TRAIN_SITE" in records[2].split("\t")[7]
        assert "POSITIVE_TRAIN_SITE" not in records[3].split("\t")[7]

        # Uncompressed recalibration VCFs use the shared Tribble LinearIndex
        # v3 writer and are directly queryable by GATK.
        plain_output = work / "recal-plain.vcf"
        plain_manifest = work / "recal-plain.json"
        plain_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
             "--max-attempts", "20", "--k-means-iterations", "20", "-O", str(plain_output),
             "--tranches-file", str(work / "recal-plain.tranches"),
             "--output-manifest", str(plain_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert plain_result.returncode == 0, plain_result.stderr
        assert plain_output.exists() and pathlib.Path(f"{plain_output}.idx").exists()
        assert not pathlib.Path(f"{plain_output}.tbi").exists()
        plain_query = work / "recal-plain-query.vcf"
        plain_query_result = subprocess.run(
            [str(gatk_java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
             "-L", "chr1:1-6", "-O", str(plain_query)],
            text=True, capture_output=True, check=False,
        )
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == 6
        plain_payload = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_payload["vcf_index_enabled"] is True
        assert plain_payload["vcf_index"] == str(plain_output) + ".idx"
        no_index_output = work / "recal-no-index.vcf"
        no_index_manifest = work / "recal-no-index.json"
        no_index_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "--max-gaussians", "1",
             "--max-attempts", "20", "--k-means-iterations", "20",
             "--create-output-variant-index=false", "-O", str(no_index_output),
             "--tranches-file", str(work / "recal-no-index.tranches"),
             "--output-manifest", str(no_index_manifest)],
            text=True, capture_output=True, check=False,
        )
        assert no_index_result.returncode == 0, no_index_result.stderr
        assert no_index_output.exists() and not pathlib.Path(f"{no_index_output}.idx").exists()
        no_index_payload = json.loads(no_index_manifest.read_text(encoding="utf-8"))
        assert no_index_payload["vcf_index_enabled"] is False and no_index_payload["vcf_index"] == ""
        invalid_index = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--create-output-variant-index=maybe",
             "-O", str(work / "invalid-index.vcf"), "--tranches-file", str(work / "invalid-index.tranches")],
            text=True, capture_output=True, check=False,
        )
        assert invalid_index.returncode != 0 and "invalid boolean" in invalid_index.stderr

        # Resource records that are already filtered are not valid GATK
        # training/known sites; they must not leak into POSITIVE_TRAIN_SITE.
        filtered_training = work / "filtered.training.vcf"
        filtered_training_rows = rows[:3] + [rows[3].replace("\tPASS\t", "\tLowQual\t")]
        write_vcf(filtered_training, filtered_training_rows)
        filtered_output = work / "filtered-resource.recal.vcf.gz"
        filtered_manifest = work / "filtered-resource.json"
        subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(filtered_training),
             "-an", "QD", "-an", "MQ", "-O", str(filtered_output),
             "--tranches-file", str(work / "filtered-resource.tranches"),
             "--output-manifest", str(filtered_manifest)],
            text=True, capture_output=True, check=True,
        )
        filtered_payload = json.loads(filtered_manifest.read_text(encoding="utf-8"))
        assert filtered_payload["training_records"] == 3, filtered_payload
        # With no advanced overrides, native defaults match GATK's
        # max-gaussians=8/max-iterations=150/k-means=100 surface (three
        # available training rows therefore yield three bounded components).
        assert filtered_payload["components"] == 3, filtered_payload
        assert filtered_payload["em_iterations"] == 150, filtered_payload
        assert filtered_payload["kmeans_iterations"] == 100, filtered_payload
        filtered_text = gzip.open(filtered_output, "rt", encoding="utf-8").read()
        filtered_records = [line for line in filtered_text.splitlines() if line and not line.startswith("#")]
        assert "POSITIVE_TRAIN_SITE" not in filtered_records[3].split("\t")[7]

        capped_output = work / "capped-training.recal.vcf.gz"
        capped_manifest = work / "capped-training.json"
        capped_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--maximum-training-variants", "2",
             "-O", str(capped_output), "--tranches-file", str(work / "capped.tranches"),
             "--output-manifest", str(capped_manifest)],
            text=True, capture_output=True, check=True,
        )
        capped_payload = json.loads(capped_manifest.read_text(encoding="utf-8"))
        assert capped_payload["maximum_training_variants"] == 2
        assert capped_payload["training_records"] == 2
        assert capped_payload["bad_cutoff_selected"] + capped_payload["bad_fallback_selected"] > 0
        assert json.loads(capped_result.stdout)["scored_records"] == 6

        # GATK's SNP mode includes MNPs (checkVariationClass), while a direct
        # type-string equality would incorrectly discard this record.
        mnp_input = work / "mnp.input.vcf"
        mnp_training = work / "mnp.training.vcf"
        # A one-row training set is now intentionally rejected by the same
        # zero-variance guard as GATK.  Keep this mode-classification fixture
        # statistically valid with three distinct MNP annotations so it tests
        # SNP-mode inclusion rather than the separate failure boundary.
        mnp_rows = [
            f"chr1\t{position}\t.\tAA\tGG\t50\tPASS\tQD={qd};MQ={mq}"
            for position, (qd, mq) in enumerate(((18, 45), (20, 50), (22, 55)), 100)
        ]
        write_vcf(mnp_input, mnp_rows)
        write_vcf(mnp_training, mnp_rows)
        mnp_output = work / "mnp.recal.vcf.gz"
        mnp_manifest = work / "mnp.json"
        mnp_result = subprocess.run(
            [str(BINARY), "-V", str(mnp_input),
             "--resource:truth,training=true,truth=true", str(mnp_training),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "-O", str(mnp_output),
             "--tranches-file", str(work / "mnp.tranches"), "--output-manifest", str(mnp_manifest)],
            text=True, capture_output=True, check=True,
        )
        assert json.loads(mnp_result.stdout)["scored_records"] == 3
        assert "VQSLOD=" in gzip.open(mnp_output, "rt", encoding="utf-8").read()
        assert (output.with_suffix(output.suffix + ".tbi")).is_file()
        tranche_rows = [line for line in tranches.read_text(encoding="utf-8").splitlines()
                        if line and not line.startswith("#") and not line.startswith("targetTruthSensitivity")]
        assert len(tranche_rows) == 2, tranche_rows

        # Scattered VariantRecalibrator runs use the VQSLOD tranche schema
        # (version 6) rather than truth-sensitivity targets.  Explicit slices
        # keep this contract test small while exercising both populated and
        # empty-above-max-LOD behavior and the deterministic row ordering.
        scatter_output = work / "scatter.recal.vcf.gz"
        scatter_tranches = work / "scatter.tranches"
        scatter_manifest = work / "scatter.json"
        scatter_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--output-tranches-for-scatter",
             "--vqslod-tranche", "20", "--vqslod-tranche=0", "--vqslod-tranche", "-5",
             "-O", str(scatter_output), "--tranches-file", str(scatter_tranches),
             "--output-manifest", str(scatter_manifest)],
            text=True, capture_output=True, check=True,
        )
        scatter_lines = scatter_tranches.read_text(encoding="utf-8").splitlines()
        assert "# Version number 6" in scatter_lines, scatter_lines[:3]
        scatter_header = "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity"
        assert scatter_header in scatter_lines, scatter_lines[:4]
        scatter_rows = [line for line in scatter_lines if line and not line.startswith("#") and line != scatter_header]
        assert len(scatter_rows) == 3, scatter_rows
        requested_vqslods = [float(row.split(",", 1)[0]) for row in scatter_rows]
        assert set(requested_vqslods) == {20.0, 0.0, -5.0}, requested_vqslods
        scatter_payload = json.loads(scatter_manifest.read_text(encoding="utf-8"))
        assert scatter_payload["scatter_tranches"] is True, scatter_payload
        assert scatter_payload["vqslod_tranches"] == 3, scatter_payload
        assert json.loads(scatter_result.stdout)["scatter_tranches"] is True
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        assert manifest_payload["tool"] == "VariantRecalibrator"
        assert manifest_payload["model"] == "full-covariance-gmm"
        assert manifest_payload["training_execution_policy"] == "VBEM-TeamPolicy-KMeans"
        assert manifest_payload["variational_normal_wishart"] is True
        assert manifest_payload["scoring_execution_policy"] == "RangePolicy"
        assert manifest_payload["training_execution_space"] == manifest_payload["execution_space"]
        assert manifest_payload["scoring_execution_space"] == manifest_payload["execution_space"]
        assert manifest_payload["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert manifest_payload["kernel_execution_policy"] == "TeamPolicy+RangePolicy"
        assert manifest_payload["kernel_batches"] > 0
        assert manifest_payload["kernel_observations"] > 0
        assert manifest_payload["kernel_prepare_seconds"] >= 0.0
        assert manifest_payload["kernel_execute_seconds"] >= 0.0
        assert manifest_payload["training_records"] == 3
        assert manifest_payload["scoreable_records"] == 6
        assert manifest_payload["normalization"]["means"] == [
            # GATK VariantDataManager.normalizeData() reorders annotation
            # dimensions by absolute training/non-training mean shift before
            # model initialization; the fixture therefore emits MQ,QD.
            sum((60, 58, 55)) / 3,
            sum((30, 28, 25)) / 3,
        ], manifest_payload["normalization"]
        assert all(value > 0.0 for value in manifest_payload["normalization"]["stdevs"])
        assert manifest_payload["resource_prior"]["log10_odds_max"] > manifest_payload["resource_prior"]["log10_odds_min"]
        assert json.loads(result.stdout)["scored_records"] == 6
        assert rscript.read_text(encoding="utf-8").startswith("# fastgatk VariantRecalibrator")
        model_text = model.read_text(encoding="utf-8")
        assert model_text.startswith("#:GATKReport.v1.1:8")
        assert "#:GATKTable:PositiveModelMeans:" in model_text
        assert "2.7666666666666668e+01" in model_text
        assert "5.7666666666666664e+01" in model_text

        # GATKReport model files are text but are routinely staged as gzip in
        # pipelines; native read/write must preserve the same table payload.
        compressed_model = work / "recal.model.gz"
        with gzip.open(compressed_model, "wt", encoding="utf-8") as stream:
            stream.write(model_text)
        compressed_output_model = work / "compressed-output.model.gz"
        compressed_output = work / "compressed-model.recal.vcf.gz"
        compressed_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--input-model", str(compressed_model),
             "-an", "QD", "-an", "MQ", "-O", str(compressed_output),
             "--tranches-file", str(work / "compressed-model.tranches"),
             "--output-model", str(compressed_output_model)],
            text=True, capture_output=True, check=True,
        )
        assert json.loads(compressed_result.stdout)["input_model"] is True
        assert gzip.open(compressed_output_model, "rt", encoding="utf-8").read() == model_text
        assert gzip.open(compressed_output, "rt", encoding="utf-8").read().count("VQSLOD=") == 6

        # Re-running the same training/scoring plan must preserve the model,
        # recal VCF payload and tranche ordering.  This guards the deterministic
        # VBEM/K-means M-step against backend-dependent accumulation order.
        repeat_output = work / "repeat.recal.vcf.gz"
        repeat_tranches = work / "repeat.recal.tranches"
        repeat_model = work / "repeat.recal.model"
        subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true,known=false,prior=15.0", str(training_vcf),
             "--resource:known,training=false,truth=false,known=true", str(known_vcf),
             "-an", "QD", "-an", "MQ", "--mode", "SNP", "-O", str(repeat_output),
             "--max-gaussians", "1", "--max-attempts", "20", "--k-means-iterations", "20",
             "--tranches-file", str(repeat_tranches), "--truth-sensitivity-tranche", "100",
             "--truth-sensitivity-tranche", "90", "--output-model", str(repeat_model),
             "--threads", "2"],
            text=True, capture_output=True, check=True,
        )
        assert gzip.open(repeat_output, "rt", encoding="utf-8").read() == text
        assert repeat_tranches.read_bytes() == tranches.read_bytes()
        assert repeat_model.read_bytes() == model.read_bytes()

        # GATK retains records with a missing annotation and marginalizes that
        # dimension during model scoring; native must not silently drop the
        # record or turn it into a hard input failure.
        missing_rows = list(rows)
        missing_rows[2] = missing_rows[2].replace("MQ=55", "MQ=.")
        missing_input = work / "missing.input.vcf"
        missing_training = work / "missing.training.vcf"
        write_vcf(missing_input, missing_rows)
        write_vcf(missing_training, missing_rows[:3])
        missing_output = work / "missing.vcf.gz"
        missing_manifest = work / "missing.json"
        missing_result = subprocess.run(
            [str(BINARY), "-V", str(missing_input),
             "--resource:truth,training=true,truth=true,prior=15.0", str(missing_training),
             "-an", "QD", "-an", "MQ", "-O", str(missing_output),
             "--tranches-file", str(work / "missing.tranches"),
             "--output-manifest", str(missing_manifest)],
            text=True, capture_output=True, check=True,
        )
        missing_payload = json.loads(missing_result.stdout)
        assert missing_payload["scoreable_records"] == 6, missing_payload
        assert missing_payload["missing_annotation_records"] == 1, missing_payload
        assert gzip.open(missing_output, "rt", encoding="utf-8").read().count("VQSLOD=") == 6
        missing_manifest_payload = json.loads(missing_manifest.read_text(encoding="utf-8"))
        assert missing_manifest_payload["scoreable_records"] == 6
        assert missing_manifest_payload["missing_annotation_records"] == 1

        # VariantDataManager.decodeAnnotation applies seeded jitter to the
        # zero/endpoint modes that otherwise collapse a Gaussian component
        # (FS/SOR/HaplotypeScore/InbreedingCoeff and MQ/AS_MQ). Exercise both
        # the MQ logit-cap path and zero-valued annotation paths, and retain
        # the controls in the manifest for deterministic replay.
        jitter_input = work / "jitter.input.vcf"
        jitter_training = work / "jitter.training.vcf"
        jitter_rows = [
            f"chr1\t{index}\t.\tA\tG\t50\tPASS\tQD={qd};MQ=60;FS=0;SOR=0.6931472;HaplotypeScore=0;InbreedingCoeff=0"
            for index, qd in enumerate((30, 28, 25, 4, 3, 2), 1)
        ]
        write_vcf_with_header(jitter_input, JITTER_HEADER, jitter_rows)
        write_vcf_with_header(jitter_training, JITTER_HEADER, jitter_rows[:3])
        jitter_output = work / "jitter.vcf.gz"
        jitter_manifest = work / "jitter.json"
        jitter_args = [
            str(BINARY), "-V", str(jitter_input),
            "--resource:truth,training=true,truth=true", str(jitter_training),
            "-an", "QD", "-an", "MQ", "-an", "FS", "-an", "SOR",
            "-an", "HaplotypeScore", "-an", "InbreedingCoeff", "--mq-cap", "60",
            "--mq-jitter", "0.05", "-O", str(jitter_output),
            "--tranches-file", str(work / "jitter.tranches"),
            "--output-manifest", str(jitter_manifest),
        ]
        jitter_result = subprocess.run(jitter_args, text=True, capture_output=True, check=True)
        jitter_payload = json.loads(jitter_manifest.read_text(encoding="utf-8"))
        assert jitter_payload["annotation_jitter"] is True, jitter_payload
        assert jitter_payload["gatk_random_seed"] == 47382911, jitter_payload
        assert jitter_payload["mq_cap"] == 60, jitter_payload
        assert jitter_payload["mq_jitter"] == 0.05, jitter_payload
        assert json.loads(jitter_result.stdout)["scored_records"] == 6
        jitter_repeat = work / "jitter-repeat.vcf.gz"
        jitter_args[jitter_args.index(str(jitter_output))] = str(jitter_repeat)
        jitter_args[jitter_args.index(str(jitter_manifest))] = str(work / "jitter-repeat.json")
        subprocess.run(jitter_args, text=True, capture_output=True, check=True)
        assert jitter_repeat.read_bytes() == jitter_output.read_bytes()

        # Missing-annotation marginalization consumes a GATK-compatible seeded
        # Java-Random sequence on Host.  Repeating the same run must therefore
        # be byte-stable across execution spaces rather than depending on
        # device scheduling or an unseeded libc RNG.
        missing_repeat = work / "missing-repeat.vcf.gz"
        subprocess.run(
            [str(BINARY), "-V", str(missing_input),
             "--resource:truth,training=true,truth=true,prior=15.0", str(missing_training),
             "-an", "QD", "-an", "MQ", "-O", str(missing_repeat),
             "--tranches-file", str(work / "missing-repeat.tranches")],
            text=True, capture_output=True, check=True,
        )
        assert missing_repeat.read_bytes() == missing_output.read_bytes()

        reused_output = work / "reused.vcf.gz"
        reused_tranches = work / "reused.tranches"
        reused_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "--input-model", str(model),
             "-an", "MQ", "-an", "QD", "--mode", "SNP", "-O", str(reused_output),
             "--tranches-file", str(reused_tranches)],
            text=True, capture_output=True, check=True,
        )
        reused_payload = json.loads(reused_result.stdout)
        assert reused_payload["input_model"] is True
        reused_text = gzip.open(reused_output, "rt", encoding="utf-8").read()
        assert reused_text.count("VQSLOD=") == 6

        gmm_output = work / "gmm.vcf.gz"
        gmm_tranches = work / "gmm.tranches"
        gmm_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--max-gaussians", "2", "--max-attempts", "8",
             "-O", str(gmm_output), "--tranches-file", str(gmm_tranches)],
            text=True, capture_output=True, check=True,
        )
        gmm_payload = json.loads(gmm_result.stdout)
        assert gmm_payload["model"] == "full-covariance-gmm"
        assert gmm_payload["components"] == 2
        gmm_text = gzip.open(gmm_output, "rt", encoding="utf-8").read()
        # Same GATK contract as above (annotation name, not the native model
        # provenance string): pinned GATK 4.6.2.0 writes culprit=MQ for this
        # fixture.
        assert "culprit=MQ" in gmm_text
        assert "culprit=full-covariance-gmm" not in gmm_text

        full_output = work / "full.vcf.gz"
        full_tranches = work / "full.tranches"
        full_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--max-gaussians", "2",
             "--full-covariance", "--max-attempts", "8", "-O", str(full_output),
             "--tranches-file", str(full_tranches)],
            text=True, capture_output=True, check=True,
        )
        full_payload = json.loads(full_result.stdout)
        assert full_payload["model"] == "full-covariance-gmm"
        assert full_payload["full_covariance"] is True
        # GATK contract again: the culprit column carries the annotation name
        # GATK's calculateWorstPerformingAnnotation selects (measured
        # culprit=MQ on this fixture with pinned GATK 4.6.2.0), not the native
        # model provenance string asserted before this fix.
        full_text = gzip.open(full_output, "rt", encoding="utf-8").read()
        assert "culprit=MQ" in full_text
        assert "culprit=full-covariance-gmm" not in full_text

        as_input = work / "as.input.vcf"
        as_training = work / "as.training.vcf"
        as_header = HEADER.replace("Number=1", "Number=A")
        as_input.write_text(as_header + "\n".join(rows) + "\n", encoding="utf-8")
        as_training.write_text(as_header + "\n".join(rows[:3]) + "\n", encoding="utf-8")
        as_output = work / "as.vcf.gz"
        as_tranches = work / "as.tranches"
        as_result = subprocess.run(
            [str(BINARY), "-V", str(as_input),
             "--resource:truth,training=true,truth=true", str(as_training),
             "-an", "QD", "-an", "MQ", "--AS", "--max-gaussians", "2",
             "-O", str(as_output), "--tranches-file", str(as_tranches)],
            text=True, capture_output=True, check=True,
        )
        as_payload = json.loads(as_result.stdout)
        assert as_payload["allele_specific"] is True
        as_text = gzip.open(as_output, "rt", encoding="utf-8").read()
        # GATK's per-allele culprit is the same annotation NAME per (ref,alt)
        # datum (VariantDataManager.java:485 collects one entry per allele into
        # AS_culprit); the previous assertion pinned the native-only model
        # provenance string "full-covariance-gmm".  The per-record values below
        # are native's deterministic output for exactly these arguments; GATK
        # itself cannot build a negative model for them (measured against the
        # pinned jar: exit 2, UserException "No data found", because no datum
        # passes the default --bad-lod-score-cutoff while the failing-STD rows
        # are excluded from selectWorstVariants), so the AS column is asserted
        # to the GATK *contract* -- requested annotation names only -- plus the
        # exact native regression values.  Byte-identical GATK gating for the
        # scalar column lives in
        # verify_variant_recalibrator_culprit_gatk_oracle.py.
        assert "AS_VQSLOD=" in as_text
        assert "AS_culprit=full-covariance-gmm" not in as_text
        as_culprits = [line.split("\t")[7].split("AS_culprit=")[1].split(";")[0]
                       for line in as_text.splitlines() if line and not line.startswith("#")]
        assert as_culprits == ["QD", "MQ", "MQ", "QD", "QD", "QD"], as_culprits
        assert all(value in {"QD", "MQ"} for value in as_culprits), as_culprits

        # The Java spelling is an optional boolean, so an explicit false must
        # select scalar scoring and an invalid literal must fail closed.
        as_false_output = work / "as-false.vcf.gz"
        as_false_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--use-allele-specific-annotations", "false",
             "-O", str(as_false_output), "--tranches-file", str(work / "as-false.tranches")],
            text=True, capture_output=True, check=True,
        )
        assert json.loads(as_false_result.stdout)["allele_specific"] is False
        as_false_text = gzip.open(as_false_output, "rt", encoding="utf-8").read()
        assert "AS_VQSLOD=" not in as_false_text
        as_invalid_result = subprocess.run(
            [str(BINARY), "-V", str(input_vcf),
             "--resource:truth,training=true,truth=true", str(training_vcf),
             "-an", "QD", "-an", "MQ", "--AS=maybe",
             "-O", str(work / "as-invalid.vcf.gz"),
             "--tranches-file", str(work / "as-invalid.tranches")],
            text=True, capture_output=True, check=False,
        )
        assert as_invalid_result.returncode != 0
        assert "invalid boolean" in as_invalid_result.stderr

        applied = work / "applied.vcf.gz"
        apply_result = subprocess.run(
            [str(APPLY_BINARY), "-V", str(input_vcf), "--recal-file", str(output),
             "--tranches-file", str(tranches), "--truth-sensitivity-filter-level", "90",
             "--mode", "SNP", "-O", str(applied)],
            text=True, capture_output=True, check=True,
        )
        applied_text = gzip.open(applied, "rt", encoding="utf-8").read()
        applied_rows = [line for line in applied_text.splitlines()
                        if line and not line.startswith("#")]
        assert len(applied_rows) == 6
        assert json.loads(apply_result.stdout)["scored_records"] == 6

        unsupported = subprocess.run(
            [str(BINARY), "-V", str(input_vcf), "-O", str(work / "unsupported.vcf"),
             "--tranches-file", str(work / "unsupported.tranches"), "--AS"],
            text=True, capture_output=True,
        )
        assert unsupported.returncode != 0
        assert "BACKEND_UNAVAILABLE" in unsupported.stderr
        print(json.dumps({
            "status": "pass",
            "input_records": 6,
            "scoreable_records": 6,
            "scored_records": 6,
            "training_records": 3,
            "tranches": 2,
            "indexed": True,
            "allele_specific_applied": True,
            "allele_specific_missing_annotation_rejected": True,
            "gmm_components": gmm_payload["components"],
            "full_covariance": full_payload["full_covariance"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
