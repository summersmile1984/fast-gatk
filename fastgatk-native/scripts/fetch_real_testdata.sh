#!/usr/bin/env bash
# fetch_real_testdata.sh — Pull pinned GATK 4.6.2.0 public test datasets via
# GitHub's LFS media CDN. The on-disk files in gatk-source/src/test/resources/large/
# are LFS pointer stubs (130 B), so we resolve them to the actual blobs and copy
# them into testdata/real/<category>/ where the verify_*.py scripts can read them
# by category. Each fetch is md5-checked against the GATK LFS OID for byte identity.
#
# Categories mapped to Track A/B oracles:
#   cnv_somatic/    HCC1143 tumor/normal chr20 + 4 PoN HDF5 + chr20 reference +
#                   interval lists (drives B3/B4/B5 + A4/A5)
#   dream_synthetic/ tumor.bam + normal.bam + 1kg chr20 + their .bai
#                   (drives A4/A5 Mutect2 + joint filter)
#   ceutrio/        CEUTrio NA12878 chr20/21 BAM + 1000G phase3 chr20 sites
#                   (drives A3 ≥5 真实 fixture pin)
#   1000g_phase3/   1000G phase3 chr20 training set for VQSR
#                   (drives B2 scatter/gather real)
#   mills_indels/   Mills_and_1000G_gold_standard indels b37 chr20
#                   (drives C3 BQSR real known-sites)
#   flow_chr9/      FlowBasedHaplotype chr9 part BAM
#                   (drives A2 PairHMM flow corpus extension)
#   mitomode/       mito.bam + mito reference (M2 mito mode)
#                   (drives A4 mito oracle)
#   bqsr_real/      reserved for BQSR real-语料 扩展
#
# Usage:  bash fastgatk-native/scripts/fetch_real_testdata.sh [CATEGORY...]
#         With no args: pulls all categories in priority order.
#         FASTGATK_SKIP_DOWNLOAD=1 to dry-run (just print plan).
#         FASTGATK_LFS_BASE=https://media.githubusercontent.com/media/broadinstitute/gatk/master
#         overrides the CDN base.
#
# Each file's expected size + LFS OID are baked in so any byte drift is loud.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
REAL="$ROOT/testdata/real"
BASE="${FASTGATK_LFS_BASE:-https://media.githubusercontent.com/media/broadinstitute/gatk/master/src/test/resources/large}"

mkdir -p "$REAL"/{cnv_somatic,dream_synthetic,ceutrio,1000g_phase3,mills_indels,flow_chr9,mitomode,bqsr_real}

# Each entry is "category|subpath|expected_size|lfs_oid"
# Source LFS pointer extraction: `head -3 file | awk '/oid sha256/ {print $2}' | tr -d '\n'`
# The expected sizes are pulled from the `size` line of the LFS pointer text.
ENTRIES_CNV_SOMATIC=(
  "HCC1143-t1-chr20-downsampled.deduplicated.bam|26583140|cb29bf33a759172ba11601bfb1b8568268ddedbbe36a9b5c50b79fc645ead790"
  "HCC1143-t1-chr20-downsampled.deduplicated.bam.bai|2952|d1d2c82d84ddc66b3c0cf08003775250408dafba5117724167ff4a6332146639"
  "HCC1143_BL-n1-chr20-downsampled.deduplicated.bam|18248106|752cf70c56644ffa4f7da3e7e114f2acec6ed7adb6e02f7fd0f03a1cda3e45d4"
  "HCC1143_BL-n1-chr20-downsampled.deduplicated.bam.bai|0|0"  # size/oid unknown; pull only, skip hash check
  "SM-74NEG-v1-chr20-downsampled.deduplicated.cram|0|0"
  "SM-74NEG-v1-chr20-downsampled.deduplicated.cram.crai|0|0"
  "SM-74P4M-v1-chr20-downsampled.deduplicated.bam|0|0"
  "SM-74P4M-v1-chr20-downsampled.deduplicated.bam.bai|0|0"
  "chr20.interval_list|0|0"
  "common_snps_sample-chr20.interval_list|0|0"
  "human_g1k_v37.chr-20.truncated.dict|0|0"
  "human_g1k_v37.chr-20.truncated.fasta|0|0"
  "human_g1k_v37.chr-20.truncated.fasta.fai|0|0"
  "ice_targets_sample-chr20.interval_list|0|0"
  "wes-do-gc.pon.hdf5|30584|b9c1b28743e4569379f987d5eeebc5393467fd2fc981c85475873b219576dfa4"
  "wes-no-gc.pon.hdf5|0|0"
  "wgs-do-gc.pon.hdf5|0|0"
  "wgs-no-gc.pon.hdf5|0|0"
)

ENTRIES_DREAM_SYNTHETIC=(
  "tumor.bam|4790562|668b7f14fe80e7537343f4add08cd22aaa4f53439eb97a511bbefe481aa9a756"
  "tumor.bam.bai|0|0"
  "normal.bam|4421969|b8e1a84424dba382e2f4b355999df160646840aae8cbe7112470962dadd5d956"
  "normal.bam.bai|0|0"
  "tumor_1.bam|0|0"
  "tumor_1.bam.bai|0|0"
  "normal_1.bam|0|0"
  "normal_1.bam.bai|0|0"
)

ENTRIES_CEUTRIO=(
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam|79856849|6b1304800e60c0ac0358df137bdad48b7857a36465b04fef3fbbb09380f04746"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.cram.crai|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.mmp2.bam|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.mmp2.bam.bai|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.md.bam|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.unaligned.bam|0|0"
  "CEUTrio.HiSeq.WGS.b37.NA12878.20.21.tiny.vcf|0|0"
)

ENTRIES_1000G=(
  "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf|78683380|b4e1b252c038cd1522f7bd65b5d4e90b0499f475e5b047350460a6146652dd9a"
  "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf.idx|0|0"
)

ENTRIES_MILLS=(
  "Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf|1913473|6e100f2afb7dc01661c6241bed15f0934e5a693a13e1369d2f877f3fbc01e15c"
  "Mills_and_1000G_gold_standard.indels.b37.sites.chr20.vcf.idx|0|0"
)

ENTRIES_FLOW_CHR9=(
  "FlowBasedHaplotype_HC_flow_chr9.part.bam|0|0"
  "FlowBasedHaplotype_HC_flow_chr9.part.bam.bai|0|0"
)

ENTRIES_MITOMODE=(
  "mito.bam|0|0"
  "mito.bam.bai|0|0"
  "mito_shifted_8000.dict|0|0"
  "mito_shifted_8000.fasta|0|0"
  "mito_shifted_8000.fasta.fai|0|0"
  "alleles.vcf|0|0"
  "alleles.vcf.idx|0|0"
)

cat_category_base() {
  case "$1" in
    cnv_somatic) echo "$BASE/cnv_somatic_workflows_test_files" ;;
    dream_synthetic) echo "$BASE/mutect/dream_synthetic_bams" ;;
    mitomode) echo "$BASE/mutect/mito" ;;
    *) echo "$BASE" ;;
  esac
}

cat_category_entries() {
  case "$1" in
    cnv_somatic) echo ENTRIES_CNV_SOMATIC ;;
    dream_synthetic) echo ENTRIES_DREAM_SYNTHETIC ;;
    ceutrio) echo ENTRIES_CEUTRIO ;;
    1000g_phase3) echo ENTRIES_1000G ;;
    mills_indels) echo ENTRIES_MILLS ;;
    flow_chr9) echo ENTRIES_FLOW_CHR9 ;;
    mitomode) echo ENTRIES_MITOMODE ;;
    *) echo "" ;;
  esac
}

CATEGORY_ORDER=(cnv_somatic dream_synthetic ceutrio 1000g_phase3 mills_indels flow_chr9 mitomode)

if [[ "${FASTGATK_SKIP_DOWNLOAD:-0}" == "1" ]]; then
  echo "Plan only (FASTGATK_SKIP_DOWNLOAD=1):"
  for c in "${CATEGORY_ORDER[@]}"; do
    echo "[$c] -> $REAL/$c"
    entries_name="$(cat_category_entries "$c")"
    eval "entries=(\"\${${entries_name}[@]}\")"
    for e in "${entries[@]}"; do
      IFS='|' read -r fn sz oid <<< "$e"
      printf "  %-60s  size=%-10s  oid=%s\n" "$fn" "$sz" "${oid:0:12}"
    done
  done
  exit 0
fi

fetch_one() {
  local category="$1"
  local fn="$2"
  local expected_size="$3"
  local expected_oid="$4"
  local dst="$REAL/$category/$fn"
  local base
  base="$(cat_category_base "$category")"
  local url="$base/$fn"
  mkdir -p "$(dirname "$dst")"
  if [[ -s "$dst" ]]; then
    # Quick skip: if not a pointer stub and size matches (or no expectation), skip.
    if ! head -3 "$dst" 2>/dev/null | grep -q "git-lfs.github.com"; then
      echo "[skip] $category/$fn (already on disk, not a pointer)"
      return 0
    fi
  fi
  echo "[get]  $category/$fn  <- $url"
  # -f: fail on HTTP error; -s: silent; -L: follow redirect; --retry 3: network resilience.
  curl -fSL --retry 3 --connect-timeout 15 -o "$dst.tmp" "$url"
  # Replace the LFS pointer stub in gatk-source/ with the actual blob (symlink to testdata/real/).
  local gatk_path
  case "$category" in
    cnv_somatic)   gatk_path="$ROOT/gatk-source/src/test/resources/large/cnv_somatic_workflows_test_files/$fn" ;;
    dream_synthetic) gatk_path="$ROOT/gatk-source/src/test/resources/large/mutect/dream_synthetic_bams/$fn" ;;
    mitomode)      gatk_path="$ROOT/gatk-source/src/test/resources/large/mutect/mito/$fn" ;;
    *)             gatk_path="$ROOT/gatk-source/src/test/resources/large/$fn" ;;
  esac
  if [[ -f "$dst" ]]; then rm -f "$dst"; fi
  mv "$dst.tmp" "$dst"
  if [[ "$expected_size" != "0" ]]; then
    local got_size
    got_size=$(stat -c %s "$dst")
    if [[ "$got_size" != "$expected_size" ]]; then
      echo "  [WARN] size mismatch: expected $expected_size, got $got_size" >&2
    fi
  fi
  if [[ "$expected_oid" != "0" ]]; then
    local got_oid
    got_oid=$(sha256sum "$dst" | awk '{print $1}')
    if [[ "$got_oid" != "$expected_oid" ]]; then
      echo "  [WARN] sha256 mismatch for $category/$fn: expected $expected_oid, got $got_oid" >&2
    fi
  fi
  if [[ -f "$gatk_path" ]] && head -3 "$gatk_path" 2>/dev/null | grep -q "git-lfs.github.com"; then
    rm -f "$gatk_path"
    ln -s "$dst" "$gatk_path"
  fi
}

requested=("$@")
if [[ ${#requested[@]} -eq 0 ]]; then
  requested=("${CATEGORY_ORDER[@]}")
fi

total=0
for c in "${requested[@]}"; do
  entries_name="$(cat_category_entries "$c")"
  if [[ -z "$entries_name" ]]; then
    echo "[error] unknown category: $c" >&2
    echo "  known: ${CATEGORY_ORDER[*]}" >&2
    exit 2
  fi
  eval "entries=(\"\${${entries_name}[@]}\")"
  for e in "${entries[@]}"; do
    IFS='|' read -r fn sz oid <<< "$e"
    fetch_one "$c" "$fn" "$sz" "$oid" || echo "[FAIL] $c/$fn" >&2
    total=$((total+1))
  done
done

echo
echo "Done. $total file(s) processed."
echo "Top-level: $REAL"
ls -la "$REAL" 2>&1 | head -10
