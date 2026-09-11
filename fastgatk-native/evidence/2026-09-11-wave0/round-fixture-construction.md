# Round: fixture construction for the seven "source-level divergence, no reproducing fixture" sites

Date: 2026-09-11 (local). Repo `/home/turing-agents/Documents/fast-gatk`.
**No production code was modified** (`git status --short` shows only new untracked files under
`fastgatk-native/scripts/` and `.diag/`), **nothing was rebuilt**, `fastgatk-native/CMakeLists.txt`
was not touched, no existing gate was edited or weakened, no root `*.md` was touched, no `git`
commit/branch. Every fixture is built inside `tempfile.TemporaryDirectory`; the printed logs are
under `.diag/`.

Binaries: `fastgatk-native/build/fastgatk-hc-call` (OpenMP/Kokkos), plus
`fastgatk-native/build/fastgatk-genotype-gvcf` for the GenotypeGVCFs probe.
Oracle: pinned GATK 4.6.2.0
(`third_party/jdk17/bin/java -Xmx1g -jar third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`).

Baseline report: `.diag/round-length-prior-sweep.md` (sections A / C.2 / F).

---

## 0. Bottom line

| # | Site | Status | Observable result |
|---|---|---|---|
| 1 | `hc_call.cpp:3358-3360` ordinary-VCF arbitrary-ploidy `*` prior | **REACHED** (divergence) | `chr1:698` QUAL **368.90 vs 369.37**, QD 8.20 vs 8.21; ploidy 4: 391.84 vs 391.94 |
| 2 | `hc_call.cpp:4485-4488` polyploid/hidden-span gVCF `*` prior | **REACHED** (divergence) | `chr1:698` QUAL **368.90 vs 369.37** (del 3), **411.25 vs 416.91** (del 2), 391.84 vs 391.94 (ploidy 4) |
| 3 | `hc_call.cpp:682-694` `genotype_priors_for_group('*')` | **REACHED** (divergence) | `chr1:698` **GQ 29 vs 33**, 61 vs 65, 67 vs 63, 27 vs 31 — the predicted 4.26-unit prior gap |
| 4 | `hc_call.cpp:3521` ordinary-VCF multi-ALT `annotations_from_owner` | **REACHED, no divergence** | multi-ALT ordinary record built; guard provably satisfied; +1 annotation pass measured; all annotations identical to GATK |
| 5 | `calling_pipeline.cpp:405-424` `prior_allele_type('*')` | **NOT REACHED** | no input found that puts a `*`/`<…>` allele into `result.candidates`; `--alleles` with `ALT=*` changes nothing |
| 6 | `calling_pipeline.cpp:6059-6074` `gvcf_candidate_is_non_monomorphic` `<NON_REF>` max() | **NOT REACHED** | retention-threshold boundary never crossed in the inputs tried |
| 7 | `hc_call.cpp:2915-2926` `estimate_mle_allele_counts` `include_non_ref` max() | **NOT REACHED** | fallback path only; the two options that would force it are rejected by native (see §D.3) |
| — | `genotype_gvcf_tool.cpp:1585-1615` `'*'`→Other (out of HC scope) | **NOT REACHED** (probe failed to launch) | `GenotypeGVCFs` invocation mismatch; see §C.4 |

The three REACHED-with-divergence sites share one regime, which is the whole point of this round:

> **The `*` pseudocount/prior is only observable when `*` is inside the AF/PL matrix but NOT inside
> the called genotype.** Once a genotype containing `*` is the called one, the allele's effective
> count is data-dominated and an 8x prior change moves the posterior by ~1e-4 — invisible at the
> printed precision. Every previous "could not reproduce" attempt was in the data-dominated regime.

---

## A. REACHED with a measured GATK-vs-native difference

### A.1 Site #1 — `hc_call.cpp:3358-3360`, arbitrary-ploidy ordinary-VCF `*` pseudocount

Deliverable script:
`fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py`

Fixture (built by the script; identical synthetic reference to the two previous rounds —
1500 bp `chr1`, `random.seed(11)`, 10×`A` at 600, 8×`CAG` at 900, base `chr1:698 = C`):

* 30 × 300M reference reads from 440, step 3
* 12 × 300M reads carrying `chr1:698 C>A` from 440, step 3
* 3 × reads from 560, step 4 with a **real `138M1D161M`** CIGAR (one-base deletion of `chr1:700`)
* `--sample-ploidy 3`, ordinary VCF (no `-ERC`), `-L chr1:500-780 --min-pruning 1`

```
timeout 1800 python3 fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py --expect-divergence   # exit 0
timeout 1800 python3 fastgatk-native/scripts/verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py                        # exit 1 (red gate)
```

Literal rows (`.diag/fixture_p5_expect.log`):

```
GATK   chr1 698 . C A 368.90 . AC=1;AF=0.333;AN=3;BaseQRankSum=0.000;DP=45;FS=0.000;MLEAC=1;MLEAF=0.333;MQ=60.00;MQRankSum=0.000;QD=8.20;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244
NATIVE chr1 698 . C A 369.37 . AC=1;AF=0.333;AN=3;BaseQRankSum=0.000;DP=45;FS=0.000;MLEAC=1;MLEAF=0.333;MQ=60.00;MQRankSum=0.000;QD=8.21;ReadPosRankSum=-2.995;SOR=0.162  GT:AD:DP:GQ:PL  0/0/1:33,12:45:58:424,0,58,1244
```

Every field except QUAL (and QUAL-derived QD) is byte-identical; the record's ALT list is `A` on
both sides, i.e. the `*` is **internal only**. Ploidy 4 with 2 deletion reads gives the same shape:
GATK `391.84` / native `391.94` (QD 8.91 on both sides).

Evidence that the divergent line executes (not merely "probably"):

1. Native's ordinary-VCF QUAL can only depend on a deletion that contributes no output ALT through
   the internal `REF/concrete/*` AF matrix, whose `*` pseudocount is assigned at exactly
   `hc_call.cpp:3358-3360`.
2. Cross-writer identity: the triploid **gVCF** record for the same pileup — where `*` is printed in
   the ALT list and therefore `include_spanning_deletion` is provably true — carries the *same*
   QUAL pair `368.90 / 369.37` (site #2 below). Two different writers producing the same
   QUAL to the cent is only possible with the same `REF/A/*` AF matrix in both.
3. Controls (all byte-identical, and they are): `ord-ploidy3-no-deletion` (no `*` at all;
   421.02/421.02), `ord-ploidy2-cigar-del3` (diploid path, where `has_spanning_deletion` is
   hard-coded false at `hc_call.cpp:3224`; 384.86/384.86), `ord-ploidy1-snp30-del6`
   (haploid; 409.04/409.04).

Expected fix location: `fastgatk-native/src/hc_call.cpp:3358-3360`
(the sweep's §C.2 item 1 one-liner).

### A.2 Site #2 — `hc_call.cpp:4485-4488`, polyploid / hidden-spanning gVCF `*` pseudocount

Deliverable script:
`fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py`

Fixture: same reference and pileup family, `-ERC BP_RESOLUTION --sample-ploidy 3`
(deletion counts 2/3, two encodings), plus ploidy 4 with 2 deletion reads.

```
timeout 1800 python3 fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py --expect-divergence   # exit 0
timeout 1800 python3 fastgatk-native/scripts/verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py                        # exit 1 (red gate)
```

Literal rows — `bp-ploidy3-cigar-del3` (`.diag/fixture_p6_expect.log`):

```
GATK   chr1 698 . C *,A,<NON_REF> 368.90 . BaseQRankSum=0.000;DP=45;MLEAC=0,1,0;MLEAF=0.00,0.333,0.00;MQRankSum=0.000;RAW_MQandDP=162000,45;ReadPosRankSum=-1.505  GT:AD:DP:GQ:PL:SB  0/1/2:30,3,12,0:45:6:430,366,431,1340,6,0,963,64,1005,1250,474,450,1385,99,1019,1205,558,1444,1197,1617:30,0,15,0
NATIVE chr1 698 . C A,*,<NON_REF> 369.37 . BaseQRankSum=0.000;DP=45;MLEAC=1,0,0;MLEAF=0.333,0.00,0.00;MQRankSum=0.000;RAW_MQandDP=162000,45;ReadPosRankSum=-1.505  GT:AD:DP:GQ:PL:SB  0/1/2:30,12,3,0:45:6:430,6,64,1250,366,0,1005,431,963,1340,474,99,1205,450,1019,1385,558,1197,1444,1617:30,0,15,0
```

After the documented canonicalization (concrete ALTs sorted, `<NON_REF>` pinned last, every
allele-indexed field permuted with them — including Number=G for arbitrary ploidy) the **only**
difference is QUAL. The `MLEAC=0,1,0` vs `1,0,0` pair is the same ALT-order effect and vanishes
under the permutation.

Two-deletion-read probe (`bp-ploidy3-cigar-del2`, reported as a report-only probe because native
does not retain the flanking `chr1:697` site and therefore prints no annotations there):

```
GATK   chr1 698 . C *,A,<NON_REF> 411.25 ... GT 0/0/2 ...
NATIVE chr1 698 . C A,*,<NON_REF> 416.91 ... GT 0/0/1 ...      # delta +5.66
```

Reachability evidence: `*` appears in **native's** emitted ALT list. That string is written by
`if (include_spanning_deletion) out << ",*";` in the same writer, so the emitted
`NATIVE ... C A,*,<NON_REF>` record is by itself proof that `include_spanning_deletion == true`,
hence that the guard at `hc_call.cpp:4485` was entered and `prior_pseudocounts[spanning]` received
the unconditional indel pseudocount before the very block hands that vector to
`calculate_allele_frequency_kokkos` (whose `qual` is what the writer prints). All four gated cases
emit `*` at `chr1:698` (the script prints `native emits '*' at chr1:698: True`).

Controls (byte-identical): `bp-ploidy2-cigar-del3` (diploid → the *already-fixed* gVCF writer at
`hc_call.cpp:4820-4836`; only the phasing triple is reported as a note on the reordered record) and
`bp-ploidy3-no-deletion` (no `*`; 421.02/421.02).

Expected fix location: `fastgatk-native/src/hc_call.cpp:4485-4488`
(same ternary shape as the landed sibling fix; sweep §C.2 item 2).

### A.3 Site #3 — `hc_call.cpp:682-694`, `genotype_priors_for_group` gives `*` the indel prior

Deliverable script:
`fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py`

This is the site the sweep called "not observable: `GT`/`GQ` byte-identical in all six gated
fixtures". The reason is now clear and mechanical: `result.genotype_priors_used =
options.use_genotype_priors && options.use_posterior_genotype_assignment`
(`calling_pipeline.cpp:14046`), and HaplotypeCaller's default is `USE_PLS_TO_ASSIGN`, in which case
the callers use `derive_genotype_gt_gq_kokkos` and **never read the priors**. Adding
`--genotype-assignment-method USE_POSTERIOR_PROBABILITIES` makes the priors live, and then the
predicted `0.426` log10 gap between GATK's
`GenotypePriorCalculator` SNP value `log10(1e-3) - log10(3) = -3.4771` for the called,
non-symbolic, 1 bp `*` allele and native's `log10(1/8000) = -3.9031` shows up directly in GQ:

```
timeout 1800 python3 fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py --expect-divergence   # exit 0
timeout 1800 python3 fastgatk-native/scripts/verify_hc_spanning_prior_genotype_gq_fixture_oracle.py                        # exit 1 (red gate)
```

Literal rows — `priors-p3-cigar-del3` (`.diag/fixture_gq_expect.log`):

```
GATK   chr1 698 . C *,A,<NON_REF> 368.90 . ... GT:AD:DP:GP:GQ:PG:PL:SB  0/0/2:30,3,12,0:45:389.23,360,455,1394,0,28.77,1021.77,88,1063.77,1304,463.23,474,1439,123,1077.77,1259,577.23,1498,1251,1666.23:29:0,34.77,...
NATIVE chr1 698 . C A,*,<NON_REF> 369.37 . ... GT:AD:DP:GQ:PL:SB        0/0/1:30,12,3,0:45:33:430,6,64,1250,366,0,1005,431,963,1340,474,99,1205,450,1019,1385,558,1197,1444,1617:30,0,15,0
```

After canonicalization the called genotype is the same (`REF/REF/A`), and **GQ is the only
differing gated field**: GATK `29` vs native `33` (= 4 phred ≈ the predicted 4.26).

| case | `*` present | priors used | GATK GQ | native GQ | verdict |
|---|---|---|---|---|---|
| `priors-p3-cigar-del3` | yes | yes | 29 | 33 | **divergence +4** |
| `priors-p3-cigar-del2` | yes | yes | 61 | 65 | **divergence +4** |
| `priors-p3-cigar-del6` | yes (in the called GT) | yes | 67 | 63 | **divergence −4** (opposite sign, same cause) |
| `priors-p4-cigar-del2` | yes | yes | 27 | 31 | **divergence +4** |
| `priors-no-spanning-deletion` | no | yes | 84 | 84 | control, identical |
| `priors-disabled-span-del` | yes | no | 6 | 6 | control, identical |

The two controls bracket the divergence exactly to the guarded line: it needs the spanning deletion
**and** posterior genotype assignment. Unsupported by anything else — the del6 case shows the
divergence is not one-signed, which is what a prior (not a likelihood) defect looks like.

Expected fix location: `fastgatk-native/src/hc_call.cpp:682-687` (`spanning_heterozygosity` must be
`log10(snpHet) - log10(3)` for a 1 bp REF, i.e. GATK's SNP class), with the sibling
`calling_pipeline.cpp:405-424` (see §C.1).

---

## B. REACHED but harmless: `hc_call.cpp:3521` (ordinary-VCF multi-ALT `annotations_from_owner`)

Deliverable script:
`fastgatk-native/scripts/verify_hc_multialt_owner_annotation_fixture_oracle.py` (green, exit 0).

The owner-pattern census could not produce a multi-ALT ordinary-VCF record from reads alone. This
round produces one **deterministically with `--alleles`**: a feature VCF whose single record carries
three concrete ALT alleles at `chr1:698` (`A,G,T`). GATK splits it into three EventMap events and
emits ONE multi-allelic ordinary record; native reproduces it byte for byte.

Guard satisfaction is provable from the output alone:

* the emitted record is multi-ALT **and** carries `BaseQRankSum`/`MQRankSum`/`ReadPosRankSum` ⇒
  `annotation_reads != nullptr`;
* the emitted ALT list is **shorter** than the forced allele list (3 forced → 2 emitted with
  `--max-alternate-alleles 2`) ⇒ `max_allele_subset_changed == true`, so the guard
  `(annotation_reads != nullptr && (max_allele_subset_changed || joint_af_allele_subset_changed))`
  is satisfied.

Independent instrumented evidence, using the existing Host debug hook
(`FASTGATK_DEBUG_ANNOTATION_POSITION=697`, i.e. 0-based POS of the multi-ALT record), which prints
one `[FASTGATK_ANNOTATION_EVIDENCE]` line per read per annotation-boundary invocation:

| case | forced ALTs | `--max-alternate-alleles` | emitted ALT | evidence lines | annotation passes |
|---|---|---|---|---|---|
| `multialt-subset-max2` | A,G,T | 2 | `A,G` | 210 (42 reads) | **5** |
| `multialt-subset-max2-noreads` | A,G,T | 2 | `A,G` | 150 (30 reads) | **5** |
| `multialt-subset-max2-ploidy3` | A,G,T | 2 | `A,G` | 210 | **5** |
| `multialt-subset-max1` | A,G,T | 1 | `A` | 210 | **5** |
| `multialt-nosubset-max3` (control) | A,G,T | 3 | `A,G,T` | 168 (42 reads) | 4 |
| `biallelic-forced-control` (control) | A,T | 6 | `A,T` | 126 (42 reads) | 3 |

Exactly **one extra annotation pass** appears iff the ALT set was reduced — the extra call is the
`hc_call.cpp:3521` block. The gate asserts that delta, so it fails if the guard stops executing.

Result: every gated field of every record is identical to GATK (GT/AD/PL/INFO annotations), i.e. the
owner-based annotation is not observably different from the flattened one in any configuration
tried (ploidy 2 and 3, with and without read support, wide and narrow `-L`, 1/2/3 emitted ALTs).
Only the separately tracked `AF`/`MLEAF` zero formatting (§D.1) differs, and it is reported as a
note rather than gated.

---

## C. NOT REACHED (inputs tried, and why the negative result is trustworthy)

### C.1 `calling_pipeline.cpp:405-424` `prior_allele_type` — `"*"`/`"<…>"` branch

`prior_allele_type` is applied only to `result.candidates[i]`
(`calling_pipeline.cpp:16603`, one call site), i.e. to the concrete AssemblyCandidates. Candidate
alternate alleles are built at three places: the SNP/indel normalizer
(`calling_pipeline.cpp:4837`), the MNP/indel builder (`...:4990`) and the **forced** allele replay
(`...:16241`, from `--alleles`). The `value == "*" || value.front() == '<'` branch can therefore only
fire if a `*`/symbolic string becomes a candidate.

Inputs tried (all with `-ERC BP_RESOLUTION --sample-ploidy 3` unless noted):

| input | result |
|---|---|
| `--alleles` feature `chr1 698 . C * . PASS .` (bgzip+tabix'd) | accepted by both callers; output identical to the no-`--alleles` run — no new candidate, no field change |
| same + `--genotype-assignment-method USE_POSTERIOR_PROBABILITIES` | identical to the no-`--alleles` posterior run (GQ 29/33) — the forced `*` neither creates a candidate nor changes the priors |
| `--alleles` feature `C *,A` | same as above |
| `--alleles` feature `C <DEL>` (symbolic deletion) | not attempted (out of the HC scope of this round; a symbolic ALT in `--alleles` is the mutect2/recovery path) |

Conclusion: no fixture found in which `prior_allele_type` receives a symbolic candidate. The
divergence is **not proven unreachable** — the `*` branch of `prior_allele_type` is dead for every
input tried, but a `--alleles` symbolic-allele path could still exist outside HC's ordinary flow
(that path is separately tracked as the symbolic-ALT emission gap). The *observable* part of this
root class is covered instead by site #3, which is the same rule stated for the spanning-deletion
allele inside `genotype_priors_for_group`.

### C.2 `calling_pipeline.cpp:6059-6074` `gvcf_candidate_is_non_monomorphic` `<NON_REF>` max()

The flag only flips when the `<NON_REF>` absent posterior crosses `-0.1*stand_call_conf = -3`.
Inputs tried: the del 1/2/3/4/6 × ploidy 1/2/3/4 × `-ERC GVCF`/`BP_RESOLUTION` matrix of §E.1, plus
shallow reference-only pileups (`ref` counts 5–30 with 0–12 variant reads), where the site is either
clearly retained or clearly dropped in **both** callers — no crossing was observed. The two callers
never disagreed on a retention decision in any of those runs. No fixture produced a
`[FASTGATK_GVCF_EMISSION]` trace pair that differed. Reported as NOT REACHED rather than unreachable:
locating the boundary needs a targeted low-depth scan that this round did not have time for.

### C.3 `hc_call.cpp:2915-2926` `estimate_mle_allele_counts(..., include_non_ref)` max()

Only reachable through the `else if (!include_spanning_deletion)` fallback at
`hc_call.cpp:4519` (polyploid gVCF), `...:3551` (ordinary) and `...:4866` (diploid gVCF), i.e. only
when the full symbolic AF matrix is unusable. Inputs tried to force that:

| input | GATK | native |
|---|---|---|
| `--heterozygosity-stdev 0` | exit 0 (produces `<NON_REF>` blocks) | **exit 2**: `error: call-confidence prior/threshold options are invalid` |
| `--indel-heterozygosity 0` | exit 0 | **exit 2**: `error: calling thresholds must be positive` |
| `--heterozygosity 0` | exit 0 | **exit 2**: `error: calling thresholds must be positive` |
| `--indel-heterozygosity 0` + `USE_POSTERIOR_PROBABILITIES` | exit 3 (GATK exception) | exit 2 (option validation) |

`genotype_heterozygosity_stdev <= 0` (the `hc_call.cpp:4497` gate that disables the AF block) is
therefore not reachable through the CLI in native, and no other input was found that makes the
kernel return `samples_with_likelihoods == 0` or a size-mismatched `integer_allele_counts` while the
record is still written. The sweep's own recommendation ("localize with a unit-level call rather
than a whole-pipeline fixture") stands.

### C.4 `genotype_gvcf_tool.cpp:1585-1615` (out of HC scope) — probe not completed

`fastgatk-native/build/fastgatk-genotype-gvcf` exists, and the obvious fixture is a GATK
`-ERC BP_RESOLUTION` gVCF carrying a `*` record (produced by the §A.2 fixture) fed to both
GenotypeGVCFs implementations. The one probe attempted failed on **argument mismatch**, not on the
site: GATK's run errored on the unindexed input VCF path and native rejected
`--add-output-vcf-command-line`
(`fastgatk-genotype-gvcf: unknown option: --add-output-vcf-command-line`). Re-running this needs the
input gVCF indexed (`IndexFeatureFile`) and the native tool's own option surface; not completed
inside this round's budget.

---

## D. Adjacent divergences found while building these fixtures (not this class, not fixed, not gated)

**(1) `INFO AF`/`MLEAF` zero formatting on multi-ALT records.** Found by the §B fixture (and it is
the only field difference there). Whenever a record has ≥2 ALTs and one allele's value needs three
decimals while another is exactly zero:

```
GATK   ... AF=0.500,0.00;... MLEAF=0.500,0.00 ...
NATIVE ... AF=0.500,0.000;... MLEAF=0.500,0.000 ...
```

GATK trims each value independently to a minimum of two decimals; native appears to apply one shared
precision to the whole list (with all-zero lists it prints `0.00,0.00` like GATK). Reproduced in
`multialt-subset-max2` (2 ALTs), `multialt-noset-max3` (3 ALTs) and `multialt-subset-max2-ploidy3`,
and in the ploidy-3/4 span-del gVCF records (`MLEAF=0.333,0.00,0.00` vs `0.333,0.000,0.000`).
Reported as a note by all four new oracles.

**(2) Missing `GP`/`PG` FORMAT fields under `--genotype-assignment-method
USE_POSTERIOR_PROBABILITIES`.** GATK emits `GT:AD:DP:GP:GQ:PG:PL[:SB]`; native emits
`GT:AD:DP:GQ:PL[:SB]`. Both `GP` and `PG` are Number=G posterior arrays. Reproduced in every
posterior-mode case (both gVCF and ordinary). Reported as a note by the §A.3 oracle.

**(3) Native rejects numeric options that GATK accepts.** `--heterozygosity-stdev 0`,
`--heterozygosity 0` and `--indel-heterozygosity 0` each make native exit 2
(`call-confidence prior/threshold options are invalid` / `calling thresholds must be positive`)
while GATK runs and produces output. Found while probing C.3.

**(4) Forced-allele selection differs when the feature record contains an insertion.**
`--alleles` `A,G,T,CA` with `--max-alternate-alleles 3` and 12 reads supporting `A`:
GATK emits `C A,G,T` and native emits `C A,CA,G` — different ALT *sets*, so the comparison is a hard
violation. Not this round's class; recorded so a future round does not rediscover it.

**(5) Zero-QUAL vs missing-QUAL on a flanking deletion record.** With 1–2 deletion reads the
`chr1:697` record carries QUAL `0` (GATK) vs `.` (native) plus missing
`BaseQRankSum`/`MQRankSum`/`ReadPosRankSum`/`RAW_MQandDP`/`SB` in native, i.e. native does not retain
that site. Visible in the `del2` probes of §A.2/§A.3; they are therefore **report-only probes**
(their QUAL divergence is still asserted, their other differences are reported as notes).

---

## E. Commands, exit statuses, and the inputs that failed

### E.1 Deliverable scripts (all reproducible, all in `fastgatk-native/scripts/`)

| script | mode | exit | log |
|---|---|---|---|
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py --expect-divergence` | 5 cases, 2 gated | **0** | `.diag/fixture_p5_expect.log` |
| same, `--case ord-ploidy3-cigar-del3` (strict, no flag) | strict parity | **1** | red gate |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py --expect-divergence` | 6 cases, 4 gated + 2 probes | **0** | `.diag/fixture_p6_expect.log` |
| same, `--case bp-ploidy3-cigar-del3` (strict) | strict parity | **1** | `.diag/fixture_p6_strict.log` |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py --expect-divergence` | 6 cases, 4 gated | **0** | `.diag/fixture_gq_expect.log` |
| same, `--case priors-p3-cigar-del3` (strict) | strict parity | **1** | red gate |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | 6 cases, green by design | **0** | `.diag/fixture_p7.log` |

Shared helpers: `fastgatk-native/scripts/hc_symbolic_prior_fixture_lib.py` (fixture builders +
GATK-conformant arbitrary-ploidy VCF canonicalizer).

### E.2 Sweeps that produced the negative results (raw logs kept)

| log | content |
|---|---|
| `.diag/fixture_sweep3.log`, `.diag/fixture_sweep3b.log` | ploidy 3/4 × deletion 0–6 × SNP 12–36 × ordinary/gVCF × `--max-alternate-alleles` 1–3 × `--alleles` 2–4 ALTs |
| `.diag/fixture_sweepB.log`, `.diag/fixture_sweepB2.log` | `--genotype-assignment-method USE_POSTERIOR_PROBABILITIES` across ploidy 2/3/4 and deletion counts 0–6 |
| `.diag/fixture_sweepC.log` | the control matrix of §A.3 (no-deletion / no-flag controls) |

Inputs tried and failed to reach a site (the list that makes the negatives trustworthy):

* ordinary VCF **without** `--alleles`: no multi-ALT record at any ploidy/window tried (reproduces
  the owner-pattern census's observation) — the ALT sets are always biallelic from assembly alone;
  two substitution alleles at one base from reads alone (`C>A` + `C>T` reads) also produced **zero
  records** (the site was not called), and `--alleles A,T --max-alternate-alleles 1` collapsed to a
  single ALT;
* `--alleles` with an insertion ALT (`A,G,T,CA`) → different ALT sets in native (§D.4);
* haploid `-ERC BP_RESOLUTION` with 6 deletion reads → GATK emits `<NON_REF>` blocks where native
  emits a concrete `chr1:697 TC>T` record, so the haploid hidden-spanning path cannot be compared
  cleanly at that site (expected divergence §F(e) of the previous report);
* ploidy 3/4 with **6 deletion reads** (the previous round's probe): the called genotype contains
  `*` on both sides → the prior divergence is real but invisible (QUAL identical to the cent);
* `--heterozygosity-stdev 0` / `--indel-heterozygosity 0` / `--heterozygosity 0`: native rejects the
  options, so the `estimate_mle_allele_counts` fallback cannot be entered that way;
* `--alleles` with `ALT=*` (alone and with posterior assignment): no observable change.

---

## F. What remains unproven

* **The exact magnitude of the fix's effect is proven only to the printed precision** of the
  compared fields (QUAL/QD 2–3 decimals, GQ integer). The internal `double` posteriors are not
  proven bit-identical, and the GQ divergence is measured as an integer difference (4 phred), not
  as a bitwise posterior.
* **Site #1's line execution is inferred, not instrumented.** The binary has no DWARF line
  information, and the ordinary-path helper was constant-propagated so its `include_spanning_deletion`
  argument is no longer visible in a register; the evidence is behavioural plus the cross-writer
  QUAL identity argument of §A.1. Sites #2 and #3 have stronger, direct evidence (`*` in the
  emitted ALT list; a prior-dependent GQ that requires the flag **and** the `*` allele).
* **Sites #5 (prior_allele_type), #6 (`gvcf_candidate_is_non_monomorphic`) and #7
  (`estimate_mle_allele_counts`)** are still **source-level divergences without a fixture**. #5 may
  be genuinely dead for HaplotypeCaller; #6 and #7 are "not reached", not "proven unreachable".
  Nothing in this round demonstrates that they are harmless.
* **The `*` prior fix is proven only for the AF/PL matrices reached by these fixtures** (ploidy 1 is
  not covered for #1/#2 beyond the haploid control, and ploidy 4 only with 2 deletion reads).
* **The `hc_call.cpp:3521` result is a negative**: the owner-based annotation is identical to the
  flattened one on 6 configurations. That does not prove the branch is harmless for partitioned
  scatter (`--scatter`-style multi-region ownership) or for the max-ALT-subset paths of the gVCF
  writers, which were not exercised.
* **`genotype_gvcf_tool.cpp:1585-1615`** was not exercised at all (§C.4).
* Every new oracle is **not registered in CTest** (editing `fastgatk-native/CMakeLists.txt` is
  forbidden this round), so `ctest` alone does not protect these rules.
