#!/bin/bash
# Fair native-vs-GATK comparison: BOTH sides pinned to exactly 1 CPU core,
# full chr20 data (311.5MB BAM, 3.23M reads; -L 20:1-63025520 covers 100% of
# the reads, no subset slicing).
# Native: core 0 (OMP_NUM_THREADS=1, hc --threads 1)
# GATK:   core 2 (-XX:ActiveProcessorCount=1 so JVM pools size to 1 CPU)
set -uo pipefail
cd /home/turing-agents/Documents/fast-gatk
JAVA=third_party/jdk17/bin/java
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
REF=testdata/downloads/reference/hs37d5.fa.gz
KS1=public_pipeline/1000G.phase3.chr20.full.vcf.gz
KS2=testdata/real/mills_indels/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.gz
RAW=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
REGION=20:1-63025520
NDIR=public_pipeline/full_chr20/fair/native
GDIR=public_pipeline/full_chr20/fair/gatk
mkdir -p $NDIR/mkdup $NDIR/bqsr $NDIR/applied $NDIR/hc $NDIR/gt $GDIR

native_chain() {
  set -e
  local T=$NDIR/TIMINGS.txt; : > $T
  local t0=$SECONDS
  OMP_NUM_THREADS=1 taskset -c 0 ./fastgatk-native/build/fastgatk-mark-duplicates -I $RAW -O $NDIR/mkdup/out.bam --metrics-file $NDIR/mkdup/metrics.txt 2>$NDIR/mkdup/err.log
  python3 -c "import pysam; pysam.index('$NDIR/mkdup/out.bam')"
  echo "mkdup wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  OMP_NUM_THREADS=1 taskset -c 0 ./fastgatk-native/build/fastgatk-bqsr -I $NDIR/mkdup/out.bam -R $REF --known-sites $KS1 --known-sites $KS2 -O $NDIR/bqsr/recal.txt 2>$NDIR/bqsr/err.log
  echo "bqsr wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  OMP_NUM_THREADS=1 taskset -c 0 ./fastgatk-native/build/fastgatk-apply-bqsr -I $NDIR/mkdup/out.bam --bqsr-recal-file $NDIR/bqsr/recal.txt -O $NDIR/applied/out.bam 2>$NDIR/applied/err.log
  echo "apply wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  OMP_NUM_THREADS=1 taskset -c 0 ./fastgatk-native/build/fastgatk-hc-call -I $NDIR/applied/out.bam -R $REF -L $REGION --gvcf --stream-by-region 1000000 --threads 1 -O $NDIR/hc/out.g.vcf.gz --telemetry $NDIR/hc/telemetry.json 2>$NDIR/hc/err.log
  echo "hc wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  OMP_NUM_THREADS=1 taskset -c 0 ./fastgatk-native/build/fastgatk-genotype-gvcf -V $NDIR/hc/out.g.vcf.gz -R $REF -O $NDIR/gt/out.vcf.gz 2>$NDIR/gt/err.log
  echo "gt wall: $((SECONDS-t0))s" >> $T
  echo NATIVE_DONE >> $T
}

gatk_chain() {
  set -e
  local T=$GDIR/TIMINGS.txt; : > $T
  local t0=$SECONDS
  taskset -c 2 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR MarkDuplicates -I $RAW -O $GDIR/mkdup.bam -M $GDIR/mkdup_metrics.txt 2>$GDIR/mkdup.err.log
  echo "mkdup wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  taskset -c 2 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR BaseRecalibrator -I $GDIR/mkdup.bam -R $REF --known-sites $KS1 --known-sites $KS2 -O $GDIR/recal.txt 2>$GDIR/bqsr.err.log
  echo "bqsr wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  taskset -c 2 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR ApplyBQSR -I $GDIR/mkdup.bam --bqsr-recal-file $GDIR/recal.txt -O $GDIR/applied.bam 2>$GDIR/apply.err.log
  echo "apply wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  taskset -c 2 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR HaplotypeCaller -I $GDIR/applied.bam -R $REF -L $REGION -ERC GVCF -O $GDIR/full.g.vcf.gz 2>$GDIR/hc.err.log
  echo "hc wall: $((SECONDS-t0))s" >> $T; t0=$SECONDS
  taskset -c 2 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR GenotypeGVCFs -R $REF -V $GDIR/full.g.vcf.gz -O $GDIR/full.vcf.gz 2>$GDIR/gt.err.log
  echo "gt wall: $((SECONDS-t0))s" >> $T
  echo GATK_DONE >> $T
}

native_chain & NP=$!
gatk_chain & GP=$!
wait $NP; echo "native exit=$?"
wait $GP; echo "gatk exit=$?"
