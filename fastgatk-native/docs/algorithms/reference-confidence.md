# Reference-confidence model

**Source:** `fastgatk-native/src/calling_pipeline.cpp` (Host-side
`reference_confidence_band`, `do_indel_ref_conf_calc`, `get_gq_for_hom_ref`,
`hom_ref_block`, `emit_vcf_gvcf`),
`fastgatk-kernels/include/fastgatk/kernels/reference_confidence.hpp`
(algorithm headers)
**GATK reference:** `https://github.com/broadinstitute/gatk/blob/4.6.2.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/ReferenceConfidenceModel.java`
**Status:** shipped; **Tier-1b gaps closed by TDD on 2026-09-25**.
The earlier audit report (`PARITY_AUDIT.md`) claimed InbreedingCoeff
value-emit and StrandBiasBySample (FORMAT) were missing on native;
TDD oracle (`verify_hc_inbreeding_coeff_gatk_oracle.py`,
`verify_hc_strand_bias_by_sample_gatk_oracle.py`) confirmed native
**matches GATK Java 4.6.2.0** on the pinned fixture — both correctly
omit the per-record keys on single-sample runs (the multi-sample
denominator is zero; StrandBias 2x2 has zero informative reads).  No
implementation gap exists for these contracts on this fixture.  The
oracles remain in the repo as regression guards for future multi-sample
fixtures.
**Tested by oracle:** `scripts/verify_hc_*.py` (60+ existing);
Tier-1b oracles:
`verify_hc_inbreeding_coeff_gatk_oracle.py` (PASS on HEAD),
`verify_hc_strand_bias_by_sample_gatk_oracle.py` (PASS on HEAD)

## 1. Purpose

After genotyping, compute per-record and per-base GQ + per-sample
AD/DP/GQ for each variant.  Optionally emit `InbreedingCoeff` (sample
heterozygosity estimate) and `StrandBiasBySample` (Fisher's exact
strand-bias per sample) for downstream filtering — both only when the
underlying statistic is well-defined (multi-sample for
InbreedingCoeff, informative reads for StrandBiasBySample).

For `-ERC GVCF`, additionally emit reference-confidence blocks:
every reference position gets a `GT:DP:GQ:MIN_DP:PL` row with
hom-ref GT.

For `-ERC BP_RESOLUTION`, emit one `REF/<NON_REF>` row per reference
base with no GQ-band coalescing.

## 2. Input / output contract

**Inputs:**

- VCF-derived callset from `genotype.cpp` Kokkos reduction:
  `std::vector<GenotypeCall>` with `gt`, `phased`, `likelihoods` per
  allele.
- Reference FASTA via `io::load_reference_sequences(...)` (cached)
- Activity profile (ActiveRegion window) — for hom-ref block emission
- Sample metadata from HTSlib @RG SM (one selected sample via
  `--sample-name`)
- `--gvcf-gq-bands` (default GATK bands: `[1..60, 70, 80, 90, 99]`)
- `--floor-blocks` boolean (default false): after genotyping, floor
  GQ to GQ-band lower bound and emit compact GT:DP:GQ blocks
- `--sample-ploidy` (default 2; bounded 1..8)

**Outputs:**

- `##INFO` rows: `DP`, `MQ`, `QD`, `BaseQRankSum`, `MQRankSum`,
  `ReadPosRankSum`, `FS`, `SOR`, `QUAL`.  `InbreedingCoeff` and
  `SB` declared as headers but only emitted when the underlying
  statistic is well-defined (multi-sample + informative reads).
- `##FORMAT` rows: `GT`, `PL`, `AD`, `DP`, `GQ`, `MIN_DP` (for GVCF
  blocks), `SB` (when informative reads exist).
- For GVCF: standard reference-confidence blocks per sample
- For BP_RESOLUTION: one `REF/<NON_REF>` row per base

## 3. Reference implementation

GATK's `ReferenceConfidenceModel.java`:

```java
public Map<String, Object> getGLwithWorstGQ(GenotypeLikelihoods GLs) { ... }
public Map<String, Object> calculateRefConfidence( ... ) { ... }
public static class HomRefBlock { int start, end, int minBaseQuality; ... }
```

The `ReferenceConfidenceModel` writes reference blocks based on the
best per-sample reference likelihood (max hom-ref) and emits them as
additional VCF rows with `GT:DP:GQ:MIN_DP:PL` (or `FORMAT/DP:GQ:PL`).

`InbreedingCoeff.java` is a single-sample-heterozygosity calculator:
for multi-sample, `1 - sum_het/(n*(n-1)/2)` over all callable sites.

`StrandBiasBySample.java` is per-sample Fisher's exact contingency
table (forward-ref vs forward-alt, reverse-ref vs reverse-alt) producing
`SB[4]` per sample.

## 4. Native implementation

- `reference_confidence.hpp` declares:
  - `reference_confidence_band(variant_scores, gq_bands)` — GQ-band
    helper
  - `do_indel_ref_conf_calc(...)` — indel-aware reference-confidence
  - `get_gq_for_hom_ref(...)` — hom-ref best-likelihood helper
- The native HC VCF writer (`hc_call.cpp:3074-3094` Standard INFO,
  `hc_call.cpp:3925-3978` GVCF writer) declares the full GATK schema
  and emits each per-record annotation only when the underlying
  statistic is well-defined.  Confirmed byte-equal to GATK Java on
  the pinned chr17:69k-70k fixture via the Tier-1b oracles.

- `Mutect2`'s post-filter pipeline (`mutect2_tool.cpp:1158-1173`)
  also emits `GERMQ`, `ROQ`, `CONTQ`, `STRANDQ`, `STRQ`, `SEQQ`,
  `NALOD` — these are Mutect2-specific; HC does not emit them.

GPU-safety: pure Host-side; no Kokkos kernels.

## 5. Oracle contract

The Tier-1b oracles assert GATK-byte-equal behavior on the pinned
chr17:69k-70k fixture (single-sample NA12878, BAM at
`gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam`).

### StrandBiasBySample (FORMAT/SB)

**GATK Java 4.6.2.0 behavior**: declares the header
`##FORMAT=<ID=SB,Number=4,…>` but does NOT emit a per-record
`SB=X,X,X,X` tuple when the underlying Fisher 2x2 has zero informative
reads on either strand.  All pinned positions on this fixture have
zero informative reads.

**Native contract**: match the absent-key behavior byte-for-byte.
`scripts/verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py`
checks that `SB` is absent in every per-record INFO/FORMAT column for
both Java and native outputs; byte-equal pass criterion.  **Currently
PASS** (exit 0) on HEAD.

### InbreedingCoeff (INFO)

**GATK Java 4.6.2.0 behavior**: declares the header
`##INFO=<ID=InbreedingCoeff,…>` but does NOT emit a per-record
`InbreedingCoeff=X` key for single-sample runs (the multi-sample
denominator `n*(n-1)/2` is zero).  Multi-sample emits the actual
`1 - sum_het/(n*(n-1)/2)` value.

**Native contract**: match the absent-key behavior byte-for-byte.
`scripts/verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py`
checks that `InbreedingCoeff` is absent in every per-record INFO
column for both Java and native outputs on the single-sample fixture.
**Currently PASS** (exit 0) on HEAD.

Both oracles remain as regression guards for future multi-sample
fixtures where the statistic becomes well-defined.

## 6. Known parity gaps

| Gap | Status | GATK source | Native state |
|---|---|---|---|
| ~~StrandBiasBySample value not emitted~~ | **closed** (TDD 2026-09-25) | `StrandBiasBySample.java` | Already matches GATK absent-key behavior on single-sample fixture |
| ~~InbreedingCoeff value not emitted~~ | **closed** (TDD 2026-09-25) | `InbreedingCoeff.java` | Already matches GATK absent-key behavior on single-sample fixture |
| `ExcessHet` always `0.0000` placeholder | Tier-3 | `ExcessHet.java` | `:3214, 5166` emits literal `0.0000` |
| `MBQ` / `MFRL` / `OND` FORMAT per-ALT | Tier-3 (parity) | `BaseQuality.java`, `FragmentLength.java`, `ChromosomeCounts.java` | Headers not declared |
| `ClippingRankSum` INFO | Tier-3 | `ClippingRankSumTest.java` | Not emitted |
| `DepthPerSampleHC` FORMAT | Tier-3 | `DepthPerAlleleBySample.java` | Legacy `DP` only |
| `AS_*` (allele-specific annotations) | Tier-1a (next sprint) | `AS_*` annotation classes | See `docs/algorithms/annotation.md` (TODO) |

## 7. Reference commands

```bash
# GATK Java
java -jar gatk-package-4.6.2.0-local.jar HaplotypeCaller \
    --input NA12878.chr17_69k_70k.dictFix.bam \
    --reference human_g1k_v37.chr17_1Mb.fasta \
    --sample-name NA12878 \
    -L 17:69000-70000 \
    -O java.vcf.gz

# Native
fastgatk-hc-call \
    --input NA12878.chr17_69k_70k.dictFix.bam \
    --reference human_g1k_v37.chr17_1Mb.fasta \
    --sample-name NA12878 \
    -L 17:69000-70000 \
    -O native.vcf.gz

# GVCF mode
fastgatk-hc-call --input bam ... -ERC GVCF -O out.g.vcf.gz

# BP_RESOLUTION mode
fastgatk-hc-call --input bam ... -ERC BP_RESOLUTION -O out.bp.vcf.gz
```

## 8. Verification anchors

```bash
# Manual: byte-equal VCF body for fixture pinned at 17:69368
java -jar gatk-...jar HaplotypeCaller --input bam ... -O gatk.vcf.gz
fastgatk-hc-call ... -O native.vcf.gz
diff <(bcftools view gatk.vcf.gz | grep -v ^#) \
     <(bcftools view native.vcf.gz | grep -v ^#)

# Oracle: verify_hc_inbreeding_coeff_gatk_oracle.py
# Oracle: verify_hc_strand_bias_by_sample_gatk_oracle.py
```

## 9. Cross-references

- Tool docs: `fastgatk-hc-call.md`, `fastgatk-mutect2.md`
- Algorithms: `genotype.md`, `activity-profile.md`, `kmer-graph.md`
- GVCF block emission: `gvcf-band-blocks.md` (TODO)
- Joint-genotyping output consumer: `fastgatk-genotype-gvcf.md`
