#!/bin/bash
# Full-chr20 germline pipeline: fast-gatk native vs GATK 4.6.2.0
# Input: complete 1000 Genomes NA12878 chr20 BAM (3,227,627 mapped reads)
# Known sites: full-chr20 1000G phase3 v5b (1.8M sites) + Mills b37 chr20 (27.6k sites)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

OUT=public_pipeline/full_chr20
RAW=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
REF=testdata/downloads/reference/hs37d5.fa.gz
KS1=public_pipeline/1000G.phase3.chr20.full.vcf.gz
KS2=testdata/real/mills_indels/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.gz
REGION=20:1-63025520

JAVA=third_party/jdk17/bin/java
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
NB=fastgatk-native/build

gatk() { "$JAVA" -Xmx8g -jar "$JAR" "$@"; }

run_step() { # run_step <logfile> <cmd...>
  local log="$1"; shift
  echo "== $(date +%H:%M:%S) $*" | tee -a "$OUT/TIMINGS.txt"
  local t0=$SECONDS
  "$@" 2>"$log"
  echo "   wall: $((SECONDS - t0))s" | tee -a "$OUT/TIMINGS.txt"
}

mkdir -p "$OUT"/{native/{mkdup,bqsr,applied,hc,gt},gatk}
: > "$OUT/TIMINGS.txt"

######################## NATIVE ########################
run_step $OUT/native/mkdup/err.log \
  $NB/fastgatk-mark-duplicates -I "$RAW" -O $OUT/native/mkdup/out.bam \
    --metrics-file $OUT/native/mkdup/metrics.txt
python3 -c "import pysam; pysam.index('$OUT/native/mkdup/out.bam')"

run_step $OUT/native/bqsr/err.log \
  $NB/fastgatk-bqsr -I $OUT/native/mkdup/out.bam -R "$REF" \
    --known-sites "$KS1" --known-sites "$KS2" \
    -O $OUT/native/bqsr/recal.txt

run_step $OUT/native/applied/err.log \
  $NB/fastgatk-apply-bqsr -I $OUT/native/mkdup/out.bam \
    --bqsr-recal-file $OUT/native/bqsr/recal.txt \
    -O $OUT/native/applied/out.bam

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
    -L "$REGION" -ERC GVCF -O $OUT/gatk/full.g.vcf.gz

run_step $OUT/gatk/gt.err.log \
  gatk GenotypeGVCFs -R "$REF" -V $OUT/gatk/full.g.vcf.gz -O $OUT/gatk/full.vcf.gz

echo "ALL DONE" | tee -a "$OUT/TIMINGS.txt"
