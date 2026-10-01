#!/bin/bash
# 16-core comparison, full chr20, no slicing of the INPUT BAM (scatter by
# interval over 20:1-63025520, 16 shards):
#   A) GATK 4.6.2.0, 16-process scatter-gather (cores 16-31)  [its only
#      multi-core mode; mkdup/apply/gt stay serial as in production WDL]
#   B) native, single process, 16 OpenMP threads (cores 0-15)
#   C) native, 16-process scatter + native GatherVcfs (cores 0-15, after B)
# Each GATK shard: -XX:ActiveProcessorCount=1 -Xmx2g, pinned to one core.
set -uo pipefail
cd /home/turing-agents/Documents/fast-gatk
JAVA=third_party/jdk17/bin/java
JAR=third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar
NB=./fastgatk-native/build
REF=testdata/downloads/reference/hs37d5.fa.gz
KS1=public_pipeline/1000G.phase3.chr20.full.vcf.gz
KS2=testdata/real/mills_indels/Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.gz
RAW=NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
BASE=public_pipeline/full_chr20/fair16
NSHARD=16

python3 - <<'EOF'
end=63025520; n=16; step=(end+n-1)//n
with open('public_pipeline/full_chr20/fair16_intervals.txt','w') as f:
    for i in range(n):
        s=i*step+1; e=min((i+1)*step,end)
        f.write(f"20:{s}-{e}\n")
EOF
mapfile -t CHUNKS < public_pipeline/full_chr20/fair16_intervals.txt

run_step() { # run_step <timings> <name> <cmd...>
  local T=$1 name=$2; shift 2
  local t0=$SECONDS
  "$@"
  local rc=$?
  echo "$name wall: $((SECONDS-t0))s rc=$rc" >> $T
  return $rc
}

######################## A) GATK 16-proc scatter (cores 16-31) ########################
gatk_chain() {
  set -e
  local GDIR=$BASE/gatk
  mkdir -p $GDIR/shards
  local T=$GDIR/TIMINGS.txt; : > $T
  run_step $T mkdup taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR MarkDuplicates -I $RAW -O $GDIR/mkdup.bam -M $GDIR/mkdup_metrics.txt 2>$GDIR/mkdup.err.log
  run_step $T mkdup-index taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR BuildBamIndex -I $GDIR/mkdup.bam -O $GDIR/mkdup.bai 2>$GDIR/mkdup-index.err.log
  # BQSR scatter
  local t0=$SECONDS i core
  for i in "${!CHUNKS[@]}"; do
    core=$((16+i))
    taskset -c $core $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR BaseRecalibrator -I $GDIR/mkdup.bam -R $REF --known-sites $KS1 --known-sites $KS2 -L "${CHUNKS[$i]}" -O $GDIR/shards/recal.$i.txt 2>$GDIR/shards/bqsr.$i.err.log &
  done
  wait
  echo "bqsr-scatter wall: $((SECONDS-t0))s" >> $T
  run_step $T bqsr-gather taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR GatherBQSRReports $(for i in "${!CHUNKS[@]}"; do echo -n "-I $GDIR/shards/recal.$i.txt "; done) -O $GDIR/recal.txt 2>$GDIR/bqsr-gather.err.log
  run_step $T apply taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR ApplyBQSR -I $GDIR/mkdup.bam --bqsr-recal-file $GDIR/recal.txt -O $GDIR/applied.bam 2>$GDIR/apply.err.log
  # HC scatter
  t0=$SECONDS
  for i in "${!CHUNKS[@]}"; do
    core=$((16+i))
    taskset -c $core $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR HaplotypeCaller -I $GDIR/applied.bam -R $REF -L "${CHUNKS[$i]}" -ERC GVCF -O $GDIR/shards/hc.$i.g.vcf.gz 2>$GDIR/shards/hc.$i.err.log &
  done
  wait
  echo "hc-scatter wall: $((SECONDS-t0))s" >> $T
  run_step $T hc-gather taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR GatherVcfs $(for i in "${!CHUNKS[@]}"; do echo -n "-I $GDIR/shards/hc.$i.g.vcf.gz "; done) -O $GDIR/full.g.vcf.gz 2>$GDIR/hc-gather.err.log
  run_step $T gvcf-index taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx2g -jar $JAR IndexFeatureFile -I $GDIR/full.g.vcf.gz 2>$GDIR/idx.err.log
  run_step $T gt taskset -c 16 $JAVA -XX:ActiveProcessorCount=1 -Xmx8g -jar $JAR GenotypeGVCFs -R $REF -V $GDIR/full.g.vcf.gz -O $GDIR/full.vcf.gz 2>$GDIR/gt.err.log
  echo GATK16_DONE >> $T
}

################ B) native single-process 16 threads (cores 0-15) ################
native_threaded() {
  set -e
  local NDIR=$BASE/native16t
  mkdir -p $NDIR/mkdup $NDIR/bqsr $NDIR/applied $NDIR/hc $NDIR/gt
  local T=$NDIR/TIMINGS.txt; : > $T
  run_step $T mkdup env OMP_NUM_THREADS=16 taskset -c 0-15 $NB/fastgatk-mark-duplicates -I $RAW -O $NDIR/mkdup/out.bam --metrics-file $NDIR/mkdup/metrics.txt 2>$NDIR/mkdup/err.log
  python3 -c "import pysam; pysam.index('$NDIR/mkdup/out.bam')"
  run_step $T bqsr env OMP_NUM_THREADS=16 taskset -c 0-15 $NB/fastgatk-bqsr -I $NDIR/mkdup/out.bam -R $REF --known-sites $KS1 --known-sites $KS2 -O $NDIR/bqsr/recal.txt 2>$NDIR/bqsr/err.log
  run_step $T apply env OMP_NUM_THREADS=16 taskset -c 0-15 $NB/fastgatk-apply-bqsr -I $NDIR/mkdup/out.bam --bqsr-recal-file $NDIR/bqsr/recal.txt -O $NDIR/applied/out.bam 2>$NDIR/applied/err.log
  run_step $T hc env OMP_NUM_THREADS=16 taskset -c 0-15 $NB/fastgatk-hc-call -I $NDIR/applied/out.bam -R $REF -L 20:1-63025520 --gvcf --stream-by-region 1000000 --threads 16 -O $NDIR/hc/out.g.vcf.gz --telemetry $NDIR/hc/telemetry.json 2>$NDIR/hc/err.log
  run_step $T gt env OMP_NUM_THREADS=16 taskset -c 0-15 $NB/fastgatk-genotype-gvcf -V $NDIR/hc/out.g.vcf.gz -R $REF -O $NDIR/gt/out.vcf.gz 2>$NDIR/gt/err.log
  echo NATIVE16T_DONE >> $T
}

################ C) native 16-proc scatter (cores 0-15) ################
native_scatter() {
  set -e
  local NDIR=$BASE/native16p
  mkdir -p $NDIR/mkdup $NDIR/shards $NDIR/bqsr $NDIR/applied $NDIR/hc $NDIR/gt
  local T=$NDIR/TIMINGS.txt; : > $T
  run_step $T mkdup taskset -c 0 $NB/fastgatk-mark-duplicates -I $RAW -O $NDIR/mkdup/out.bam --metrics-file $NDIR/mkdup/metrics.txt 2>$NDIR/mkdup/err.log
  python3 -c "import pysam; pysam.index('$NDIR/mkdup/out.bam')"
  local t0=$SECONDS i
  for i in "${!CHUNKS[@]}"; do
    taskset -c $i $NB/fastgatk-bqsr -I $NDIR/mkdup/out.bam -R $REF --known-sites $KS1 --known-sites $KS2 -L "${CHUNKS[$i]}" -O $NDIR/shards/recal.$i.txt 2>$NDIR/shards/bqsr.$i.err.log &
  done
  wait
  echo "bqsr-scatter wall: $((SECONDS-t0))s" >> $T
  run_step $T bqsr-gather taskset -c 0 $NB/fastgatk-gather-bqsr-reports $(for i in "${!CHUNKS[@]}"; do echo -n "-I $NDIR/shards/recal.$i.txt "; done) -O $NDIR/bqsr/recal.txt 2>$NDIR/shards/bqsr-gather.err.log
  run_step $T apply taskset -c 0 $NB/fastgatk-apply-bqsr -I $NDIR/mkdup/out.bam --bqsr-recal-file $NDIR/bqsr/recal.txt -O $NDIR/applied/out.bam 2>$NDIR/applied/err.log
  python3 -c "import pysam; pysam.index('$NDIR/applied/out.bam')"
  t0=$SECONDS
  for i in "${!CHUNKS[@]}"; do
    taskset -c $i $NB/fastgatk-hc-call -I $NDIR/applied/out.bam -R $REF -L "${CHUNKS[$i]}" --gvcf --stream-by-region 1000000 --threads 1 -O $NDIR/shards/hc.$i.g.vcf.gz 2>$NDIR/shards/hc.$i.err.log &
  done
  wait
  echo "hc-scatter wall: $((SECONDS-t0))s" >> $T
  run_step $T hc-gather taskset -c 0 $NB/fastgatk-gather-vcfs --allow-overlaps $(for i in "${!CHUNKS[@]}"; do echo -n "-I $NDIR/shards/hc.$i.g.vcf.gz "; done) -O $NDIR/hc/out.g.vcf.gz 2>$NDIR/hc/gather.err.log
  t0=$SECONDS
  for i in "${!CHUNKS[@]}"; do
    taskset -c $i $NB/fastgatk-genotype-gvcf -V $NDIR/shards/hc.$i.g.vcf.gz -R $REF -O $NDIR/shards/gt.$i.vcf.gz 2>$NDIR/shards/gt.$i.err.log &
  done
  wait
  echo "gt-scatter wall: $((SECONDS-t0))s" >> $T
  run_step $T gt-gather taskset -c 0 $NB/fastgatk-gather-vcfs --allow-overlaps $(for i in "${!CHUNKS[@]}"; do echo -n "-I $NDIR/shards/gt.$i.vcf.gz "; done) -O $NDIR/gt/out.vcf.gz 2>$NDIR/gt/gather.err.log
  echo NATIVE16P_DONE >> $T
}

mkdir -p $BASE
gatk_chain & GP=$!
native_threaded & NP=$!
wait $GP; echo "gatk16 exit=$?"
wait $NP; echo "native16t exit=$?"
native_scatter; echo "native16p exit=$?"
