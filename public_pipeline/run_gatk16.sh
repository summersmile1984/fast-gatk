#!/bin/bash
# GATK 16-process scatter-gather, full chr20 (cores 16-31). Rerun of chain A
# with BuildBamIndex after MarkDuplicates (scatter -L requires the index).
set -euo pipefail
cd /home/turing-agents/Documents/fast-gatk
JAVA=third_party/jdk17/bin/java
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
REF=testdata/downloads/reference/hs37d5.fa.gz
KS1=public_pipeline/1000G.phase3.chr20.full.vcf.gz
KS2=testdata/real/mills_indels/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.gz
RAW=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
GDIR=public_pipeline/full_chr20/fair16/gatk
mapfile -t CHUNKS < public_pipeline/full_chr20/fair16_intervals.txt

run_step() {
  local T=$1 name=$2; shift 2
  local t0=$SECONDS
  "$@"
  local rc=$?
  echo "$name wall: $((SECONDS-t0))s rc=$rc" >> $T
  return $rc
}

rm -rf $GDIR; mkdir -p $GDIR/shards
T=$GDIR/TIMINGS.txt; : > $T
run_step $T mkdup taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR MarkDuplicates -I $RAW -O $GDIR/mkdup.bam -M $GDIR/mkdup_metrics.txt 2>$GDIR/mkdup.err.log
run_step $T mkdup-index taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR BuildBamIndex -I $GDIR/mkdup.bam -O $GDIR/mkdup.bai 2>$GDIR/mkdup-index.err.log
t0=$SECONDS
for i in "${!CHUNKS[@]}"; do
  taskset -c $((16+i)) $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR BaseRecalibrator -I $GDIR/mkdup.bam -R $REF --known-sites $KS1 --known-sites $KS2 -L "${CHUNKS[$i]}" -O $GDIR/shards/recal.$i.txt 2>$GDIR/shards/bqsr.$i.err.log &
done
wait
echo "bqsr-scatter wall: $((SECONDS-t0))s" >> $T
run_step $T bqsr-gather taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR GatherBQSRReports $(for i in "${!CHUNKS[@]}"; do echo -n "-I $GDIR/shards/recal.$i.txt "; done) -O $GDIR/recal.txt 2>$GDIR/bqsr-gather.err.log
run_step $T apply taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR ApplyBQSR -I $GDIR/mkdup.bam --bqsr-recal-file $GDIR/recal.txt -O $GDIR/applied.bam 2>$GDIR/apply.err.log
t0=$SECONDS
for i in "${!CHUNKS[@]}"; do
  taskset -c $((16+i)) $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR HaplotypeCaller -I $GDIR/applied.bam -R $REF -L "${CHUNKS[$i]}" -ERC GVCF -O $GDIR/shards/hc.$i.g.vcf.gz 2>$GDIR/shards/hc.$i.err.log &
done
wait
echo "hc-scatter wall: $((SECONDS-t0))s" >> $T
run_step $T hc-gather taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR GatherVcfs $(for i in "${!CHUNKS[@]}"; do echo -n "-I $GDIR/shards/hc.$i.g.vcf.gz "; done) -O $GDIR/full.g.vcf.gz 2>$GDIR/hc-gather.err.log
run_step $T gt taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR GenotypeGVCFs -R $REF -V $GDIR/full.g.vcf.gz -O $GDIR/full.vcf.gz 2>$GDIR/gt.err.log
echo GATK16_DONE >> $T
