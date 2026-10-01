#!/bin/bash
# Full germline short-variant pipeline: fast-gatk native vs GATK 4.6.2.0
# Input: public_pipeline/NA12878.chr20.10M-11M.bam (1000 Genomes, hs37d5 slice)
# Fixes vs pipeline_v3:
#   1. native mkdup output is indexed (pysam) so GATK interval traversal works
#      (v3 GATK HC exited 2 "not indexed", leaving a header-only g.vcf.gz)
#   2. mkdir -p before every step (native apply-bqsr does not create parent dirs;
#      v3 failed with hts_open "No such file or directory")
#   3. native HC runs on BQSR-applied BAM, matching the GATK side
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

OUT=public_pipeline/pipeline_v4
RAW=public_pipeline/NA12878.chr20.10M-11M.bam
REF=testdata/downloads/reference/hs37d5.fa.gz
KS1=testdata/real/1000g_phase3/1000G.phase3.broad.withGenotypes.chr20.10100000.vcf.gz
KS2=testdata/real/mills_indels/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.gz
REGION=20:10000001-11000000
REGION_GATK=20:10000001-11000000

JAVA=third_party/jdk17/bin/java
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
NB=fastgatk-native/build

gatk() { "$JAVA" -Xmx4g -jar "$JAR" "$@"; }

run_step() { # run_step <logfile> <cmd...>
  local log="$1"; shift
  echo "== $(date +%H:%M:%S) $*" | tee -a "$OUT/TIMINGS.txt"
  local t0=$SECONDS
  "$@" 2>"$log"
  echo "   wall: $((SECONDS - t0))s" | tee -a "$OUT/TIMINGS.txt"
}

rm -rf "$OUT"
mkdir -p "$OUT"/{native/{mkdup,bqsr,applied,hc,gt},gatk,diag}
: > "$OUT/TIMINGS.txt"

######################## NATIVE ########################
run_step $OUT/native/mkdup/err.log \
  $NB/fastgatk-mark-duplicates -I "$RAW" -O $OUT/native/mkdup/out.bam \
    --metrics-file $OUT/native/mkdup/metrics.txt

# FIX v3#1: index native BAM so GATK (and any index-requiring consumer) works
python3 -c "import pysam; pysam.index('$OUT/native/mkdup/out.bam')"

run_step $OUT/native/bqsr/err.log \
  $NB/fastgatk-bqsr -I $OUT/native/mkdup/out.bam -R "$REF" \
    --known-sites "$KS1" --known-sites "$KS2" \
    -O $OUT/native/bqsr/recal.txt

# FIX v3#2: parent dir created above (mkdir -p .../applied)
run_step $OUT/native/applied/err.log \
  $NB/fastgatk-apply-bqsr -I $OUT/native/mkdup/out.bam \
    --bqsr-recal-file $OUT/native/bqsr/recal.txt \
    -O $OUT/native/applied/out.bam

# FIX v3#3: HC on BQSR-applied BAM (parity with GATK side)
run_step $OUT/native/hc/err.log \
  $NB/fastgatk-hc-call -I $OUT/native/applied/out.bam -R "$REF" \
    -L "$REGION" --gvcf -O $OUT/native/hc/out.g.vcf.gz

run_step $OUT/native/gt/err.log \
  $NB/fastgatk-genotype-gvcf -V $OUT/native/hc/out.g.vcf.gz -R "$REF" \
    -O $OUT/native/gt/out.vcf.gz

######################## GATK 4.6.2.0 ########################
run_step $OUT/gatk/mkdup.err.log \
  gatk MarkDuplicates -I "$RAW" -O $OUT/gatk/mkdup.bam -M $OUT/gatk/mkdup_metrics.txt

run_step $OUT/gatk/bqsr.err.log \
  gatk BaseRecalibrator -I $OUT/gatk/mkdup.bam -R "$REF" \
    --known-sites "$KS1" --known-sites "$KS2" -O $OUT/gatk/recal.txt

run_step $OUT/gatk/apply.err.log \
  gatk ApplyBQSR -I $OUT/gatk/mkdup.bam --bqsr-recal-file $OUT/gatk/recal.txt \
    -O $OUT/gatk/applied.bam

run_step $OUT/gatk/hc.err.log \
  gatk HaplotypeCaller -I $OUT/gatk/applied.bam -R "$REF" \
    -L "$REGION_GATK" -ERC GVCF -O $OUT/gatk/full.g.vcf.gz

run_step $OUT/gatk/gt.err.log \
  gatk GenotypeGVCFs -R "$REF" -V $OUT/gatk/full.g.vcf.gz -O $OUT/gatk/full.vcf.gz

######################## DIAGNOSTIC ########################
# v3 root-cause check: GATK HC on native mkdup output (now indexed) must emit records
run_step $OUT/diag/gatk_on_native_mkdup.err.log \
  gatk HaplotypeCaller -I $OUT/native/mkdup/out.bam -R "$REF" \
    -L "$REGION_GATK" -ERC GVCF -O $OUT/diag/gatk_on_native_mkdup.g.vcf.gz

echo "ALL DONE" | tee -a "$OUT/TIMINGS.txt"
