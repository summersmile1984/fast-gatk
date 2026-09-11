# Round — oracle integrity: a missing GATK oracle can no longer pass for a parity check

Scope of this round: **measurement integrity of the test suite only.** No expected value was
re-valued, no assertion was weakened or deleted, no production code (`fastgatk-native/src`,
`include`, `kernels`) was touched, nothing was rebuilt, `fastgatk-native/CMakeLists.txt` and the
root `.md` files are untouched, and nothing was committed or branched — the change is left in the
working tree.

## 0. Headline measurements

| measurement | value |
| --- | --- |
| `fastgatk-native/scripts/verify_*.py` scanned | 268 |
| scripts whose GATK comparison sits behind a presence guard that **never asserts** the oracle exists | **183** — 176 registered, 7 unregistered |
| registered CTest tests (of 294) whose command is one of those scripts | **187** (63.6%) |
| of those 176 registered scripts, how many had *no* `FASTGATK_REQUIRE_GATK_ORACLE` escape at all before this round | **51** (the other 125 mentioned the variable, mostly as a bare `== "1"` early exit; their skip lines go to stdout, which `--output-on-failure` discards, so a pass stayed indistinguishable either way) |
| of those, scripts where the guard **wraps the comparison itself** (the shape named in the audit) | **25**, of which **23** had no escape at all |
| of those, scripts whose only oracle guards are **skips** (shape B) | 151 |
| **registered tests currently running with no real GATK comparison while appearing to make one** | **0** — see §1.3. The exposure was latent, not active |
| scripts changed | 176 modified + `fastgatk-native/scripts/oracle_guard.py` added |

The audit's premise is confirmed, and the exposure is larger than "11 or more": it is **176
registered scripts / 187 registered tests**. What is *not* true today is that they are actually
running without the oracle — the vendored JDK and the pinned jar are both present and every one of
those scripts resolves to them (§1.3). The defect is that nothing in the suite would have noticed
if they had not.

## 1. TASK 1 — measuring the exposure

### 1.1 Method

I did not grep for one string. The scan parses each of the 268 `verify_*.py` files with `ast` and
looks for every `if` whose condition tests the existence of a path that **resolves to the vendored
JDK or the pinned GATK jar**. Resolution is done by evaluating the module's own `name = expr`
assignments (`root = Path(__file__).resolve().parents[2]` and friends) with `JAVA` and `GATK_JAR`
cleared, then checking the resolved value's string form for the markers `third_party/jdk17` and
`gatk-package`. The classification is therefore semantic rather than textual, which matters because:

* the corpus spells the jar a dozen ways — `jar`, `JAR`, `gatk`, `gatk_jar`, `GATK_JAR`, or the
  full `root / "third_party/gatk-package/gatk-4.6.2.0/..."` path written inline;
* the corpus reaches the oracle through containers — `required = (JAVA, GATK, NATIVE, ...)` then
  `if not all(path.is_file() for path in required):` — where the condition text never mentions
  java or gatk at all. A text-only scan misses 92 of them;
* the corpus sometimes computes a boolean first — `java_oracle = JAVA.exists() and GATK_JAR.exists()`
  then `if java_oracle:` — which no `if java.exists()` pattern matches.

Scanner and patcher are the same file (`.diag/oracle_guard_patch.py`), so the exposure numbers and
the applied patch cannot disagree. Raw artefacts: `.diag/oracle_table.json` (enriched scan),
`.diag/oracle_scan.json` (first, text-oriented pass), `.diag/oracle_patch.applied.diff`,
`.diag/oracle_patch_apply.log`.

### 1.2 The two guard shapes found

**Shape A — the guard wraps the GATK comparison (the audit's pattern).** When the oracle is absent
the comparison is silently not executed and the test carries on against its recorded native
expectations:

```python
java = root / "third_party/jdk17/bin/java"
jar  = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
java_checked = False
if java.exists() and jar.exists():        # never asserted
    ...run GATK, compare against native output, set java_checked = True...
print(json.dumps({"status": "pass", "java_oracle": java_checked}))   # -> "pass" either way
```

The script even reports `"java_oracle": false` — but the *exit code is 0*, so CTest prints
`Passed` and nobody reads the JSON. This is precisely the `verify_indel.py:173` failure mode made
invisible.

**Shape B — the guard skips the whole test.** `if not JAVA.exists() or not GATK.exists():` → print
a skip line → `return 0`. Green, no GATK, and in most of these scripts the skip line is the only
trace, and it is discarded by `--output-on-failure`.

Both shapes fall under the requested policy, so both were fixed. Only shape A silently degrades
into "native compared against native expectations"; shape B silently reports nothing at all.

### 1.3 Does the jar path each script computes actually exist? — yes

Resolved for every guarded script with a clean environment: every oracle expression lands on
`third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar` (425 767 494 bytes, present)
and `third_party/jdk17/bin/java` (executable, present). **No registered script produced a missing
oracle path** — `paths_missing` is empty for every row of the table below. Therefore:

> **registered tests currently running with no real GATK comparison while appearing to make one: 0.**

Their guards evaluate true, the GATK comparisons do run, and the scripts' own JSON says
`"java_oracle": true`. The exposure is what happens the day the jar or the JDK is not vendored,
moves, is overridden by a stale `JAVA`/`GATK_JAR`, or is pinned to a different GATK release —
measured directly in §3.2, and confirmed suite-wide in §3.1.

### 1.4 The scripts the audit named — confirmed

| audit entry | measured |
| --- | --- |
| `verify_filter_intervals.py:177` | confirmed; there is a second site at `:242` — both shape A |
| `verify_preprocess_intervals.py:137` | confirmed, shape A |
| `verify_annotate_intervals.py:89,134,182,208` | confirmed, all four shape A |
| `verify_fasta_reference_tools.py:105-236` | confirmed: 11 `if java_oracle:` sites (105,116,130,150,156,178,184,203,208,231,236), shape A via `java_oracle = JAVA.exists() and GATK_JAR.exists()` at `:89` |
| `verify_count_bases_in_reference.py:121` | confirmed, shape A |
| `verify_compare_references.py:70` | confirmed, shape A |
| `verify_check_reference_compatibility.py:84` | confirmed, shape A |
| `verify_left_align.py:147,271` | confirmed, plus a third at `:333` — all shape A |
| `verify_validate_variants.py:115,191` | confirmed, both shape A |
| `verify_variant_recalibrator.py` ("etc.") | **not** a guard case: it invokes the jar only as a *reader* (`SelectVariants`, `:109`), unconditionally. With the jar missing the subprocess fails loudly. Its issue is the already-catalogued "GATK as reader, not comparator" one, out of scope here |
| `verify_indel.py:173` (the incident) | **not** a guard case either: `verify_indel.py` never invokes the jar, so it has no oracle guard to fail closed. Its `:173` incident is the assertion-level one the audit already describes — the *consequence* of an unverified expectation, not a guard |

### 1.5 Full exposure table

`guard` column: line number and shape (`positive` = shape A, `negative` = shape B).
`behaviour when the oracle is missing` describes the **pre-fix** behaviour.

| script | guard line(s) | registered | jar/JDK path exists | behaviour when the oracle is missing |
| --- | --- | --- | --- | --- |
| `verify_analyze_covariates.py` | :32 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_analyze_covariates_bqsr_alias_gatk_oracle.py` | :44 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_annotate_intervals.py` | :89 (positive), :134 (positive), :182 (positive), :208 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_apply_bqsr_alias_gatk_oracle.py` | :40 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_apply_vqsr.py` | :415 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_apply_vqsr_default_cutoff_gatk_oracle.py` | :63 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_apply_vqsr_exclude_intervals_gatk_oracle.py` | :99 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_apply_vqsr_filter_booleans_gatk_oracle.py` | :77 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_apply_vqsr_gatk_oracle.py` | :70 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_apply_vqsr_sites_only_gatk_oracle.py` | :72 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr.py` | :433 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_bqsr_context_size_gatk_oracle.py` | :71 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_cram_gatk_oracle.py` | :67 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_gatk_oracle.py` | :63 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_indel_gatk_oracle.py` | :52 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_long_read_gatk_oracle.py` | :78 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_preserve_gatk_oracle.py` | :41 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_read_filter_gatk_oracle.py` | :57 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_bqsr_report_roundtrip_gatk_oracle.py` | :39 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_calculate_contamination_gatk_oracle.py` | :56 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py` | :71 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_call_copy_ratio_segments_gatk_oracle.py` | :102 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py` | :29 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py` | :48 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_check_reference_compatibility.py` | :84 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_collect_allelic_counts_gatk_oracle.py` | :121 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_collect_f1r2_counts.py` | :231 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_collect_read_counts_gatk_oracle.py` | :39 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_combine_gvcfs_gatk_oracle.py` | :82 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_combine_gvcfs_interval_refblock_gatk_oracle.py` | :67 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_combine_gvcfs_plless_oracle.py` | :54 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_compare_references.py` | :70 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_count_bases_in_reference.py` | :121 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_count_reads.py` | :421 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_create_read_count_panel_of_normals.py` | :102 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py` | :48 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_depth_of_coverage.py` | :111 (positive), :266 (positive), :384 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py` | :64 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_depth_of_coverage_read_filter_gatk_oracle.py` | :32 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_fasta_alternate_iupac_hom_gatk_oracle.py` | :48 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_fasta_reference_tools.py` | :105 (positive), :116 (positive), :130 (positive), :150 (positive), :156 (positive), :178 (positive), :184 (positive), :203 (positive), :208 (positive), :231 (positive), :236 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_filter_intervals.py` | :177 (positive), :242 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_filter_mutect_calls_contamination_oracle.py` | :41 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_calls_germline_oracle.py` | :45 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_contamination_joint_oracle.py` | :60 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_dream_synthetic_joint_oracle.py` | :186 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_hcc1143_joint_oracle.py` | :148 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_normal_artifact_oracle.py` | :45 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_orientation_gatk_oracle.py` | :104 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_orientation_joint_oracle.py` | :51 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_filter_mutect_variant_index_alias_gatk_oracle.py` | :40 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_flag_stat.py` | :453 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_gather_bqsr_gatk_oracle.py` | :61 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gather_pileup_gatk_oracle.py` | :33 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gather_tranches.py` | :99 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gather_vcfs_cli_boundary_gatk_oracle.py` | :38 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_genotype_gvcf.py` | :62 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_genotype_gvcf_inbreeding.py` | :65 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_genotype_gvcf_legacy_qual.py` | :55 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_genotype_gvcf_multiallelic.py` | :73 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_genotype_gvcf_multisample.py` | :55 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gatk_oracle.py` | :150 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_bridge.py` | :50 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_import_gatk_oracle.py` | :31 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_import_native_interval_gatk_oracle.py` | :66 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_import_sample_map_gatk_oracle.py` | :47 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_import_update_workspace_gatk_oracle.py` | :56 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genomicsdb_native_storage_boundary.py` | :32 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_assignment_gatk_oracle.py` | :56 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_exclude_intervals_gatk_oracle.py` | :43 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_gp_input_gatk_oracle.py` | :48 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_max_alternate_alleles_gatk_oracle.py` | :35 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_multisample_reference_confidence_oracle.py` | :35 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_spanning_deletion_gatk_oracle.py` | :71 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_genotype_gvcf_starts_in_intervals_gatk_oracle.py` | :65 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_get_pileup_gatk_oracle.py` | :25 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_gvcf_stream_overlapping_indels_gatk_oracle.py` | :141 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_af_zero_format_gatk_oracle.py` | :379 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_alleles_gatk_oracle.py` | :116 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_alleles_overlap_gate_oracle.py` | :273 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_arbitrary_ploidy_span_del_prior_fixture_oracle.py` | :232 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_assembly_region_boundary_gatk_oracle.py` | :37 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_bp_resolution_gatk_oracle.py` | :40 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_broad_gatk_oracle.py` | :84 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr17_69k_70k_gatk_oracle.py` | :141 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr20_100k_nocall_gatk_oracle.py` | :38 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr20_max_mnp_gvcf_gatk_oracle.py` | :47 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr20_max_mnp_polyploid_gvcf_gatk_oracle.py` | :130 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py` | :47 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_chr20_real_contract.py` | :191 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_cigar_indel_activity_gatk_oracle.py` | :71 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_complex_multiallelic_oracle.py` | :78 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_dense_gvcf_genotype_gatk_oracle.py` | :72 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_genotype_priors.py` | :28 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_gq_bands_gatk_oracle.py` | :36 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_gvcf_indel_end_gatk_oracle.py` | :246 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_gvcf_symbolic_prior_gatk_oracle.py` | :560 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_indel_zero_gatk_oracle.py` | :36 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_informative_overlap_margin_gatk_oracle.py` | :37 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_issue3845_gatk_oracle.py` | :41 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_kmer_list_gatk_oracle.py` | :39 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_likelihood_filter.py` | :42 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_min_base_quality_boundary_gatk_oracle.py` | :38 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_min_base_quality_gatk_oracle.py` | :42 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_min_pruning_gatk_oracle.py` | :44 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_multialt_owner_annotation_fixture_oracle.py` | :276 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_polyploid_gatk_oracle.py` | :51 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_polyploid_gvcf_span_del_prior_fixture_oracle.py` | :257 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_rcm_pl_range_gatk_oracle.py` | :46 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_region_streaming.py` | :174 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_hc_sites_only_gatk_oracle.py` | :46 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_softclip_contig_start_gatk_oracle.py` | :40 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_softclip_gatk_oracle.py` | :60 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_softclip_low_quality_gatk_oracle.py` | :38 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_span_del_qual_oracle.py` | :393 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_spanning_prior_genotype_gq_fixture_oracle.py` | :277 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_window_invariance_gatk_oracle.py` | :336 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_index_feature_file.py` | :79 (positive), :117 (positive), :163 (positive), :198 (positive), :219 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_learn_read_orientation_model.py` | :139 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_learn_read_orientation_model_multisample_gatk_oracle.py` | :83 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_left_align.py` | :147 (positive), :271 (positive), :333 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_left_align_cli_boundary_gatk_oracle.py` | :45 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_left_align_gatk_oracle.py` | :118 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_left_align_sites_only_gatk_oracle.py` | :39 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mark_duplicates.py` | :286 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_mark_duplicates_gatk_oracle.py` | :262 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mark_duplicates_pair_key_gatk_oracle.py` | :53 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mark_duplicates_tagging_policy_gatk_oracle.py` | :38 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2.py` | :933 (positive), :1090 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_mutect2_bp_resolution_gatk_contract.py` | :39 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_dream_low_bq_gatk_oracle.py` | :153 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_dream_synthetic_oracle.py` | :175 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_feature_resource_gatk_oracle.py` | :131 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_force_active_gatk_oracle.py` | :56 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_gatk_oracle.py` | :172 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_gvcf_eventmap_gatk_contract.py` | :90 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_gvcf_reference_blocks_gatk_contract.py` | :59 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_hcc1143_chr20_oracle.py` | :172 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_issue3845_gatk_oracle.py` | :42 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_kmer_list_gatk_oracle.py` | :55 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_mismapping_rate_boundary.py` | :71 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_mito_interval_halo_gatk_oracle.py` | :40 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_normal_lod_gatk_oracle.py` | :94 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_reference_confidence_eventmap_full_gatk_contract.py` | :45 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_mutect2_zero_lod_indel_activity_gatk_oracle.py` | :96 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_pairhmm_default_indel_quality_gatk_oracle.py` | :47 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_preprocess_intervals.py` | :137 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_reblock_gatk_oracle.py` | :53 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_reblock_gatk_shards.py` | :85 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_reblock_overlap_gatk_oracle.py` | :23 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_select_variants_filtered_nocall_oracle.py` | :121 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_select_variants_filtered_oracle.py` | :62 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_select_variants_gatk_oracle.py` | :74 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_select_variants_sites_only_gatk_oracle.py` | :72 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_shift_fasta.py` | :74 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_sort_sam.py` | :153 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_sort_sam_cli_boundary_gatk_oracle.py` | :83 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_sort_sam_duplicate_gatk_oracle.py` | :70 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_sort_sam_gatk_oracle.py` | :152 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_validate_variants.py` | :115 (positive), :191 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_validate_variants_gatk_oracle.py` | :49 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_validate_variants_symbolic_gatk_oracle.py` | :43 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_eval_gatk_oracle.py` | :67 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_eval_keep_ac0_gatk_oracle.py` | :80 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_eval_validation_report_gatk_oracle.py` | :137 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_filtration_gatk_oracle.py` | :156 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_filtration_missing_boolean_gatk_oracle.py` | :57 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_filtration_set_nocall_gatk_oracle.py` | :44 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_annotation_order_gatk_oracle.py` | :75 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_attempt_iteration_gatk_oracle.py` | :47 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_gatk_model_oracle.py` | :93 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_sample_every_gatk_oracle.py` | :69 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_vbem_gatk_oracle.py` | :110 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variant_recalibrator_zero_variance_gatk_oracle.py` | :37 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_variants_to_table.py` | :201 (positive), :288 (positive), :323 (positive), :335 (positive), :383 (positive), :412 (positive) | yes | yes | positive guard not taken -> GATK comparison silently skipped, assertions fall back to recorded native expectations, exit 0 |
| `verify_variants_to_table_gatk_oracle.py` | :55 (negative) | yes | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_alleles_deep_boundary.py` | :210 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_alleles_deep_limits.py` | :106 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_forced_alleles_emission_gate_oracle.py` | :296 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py` | :279 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_hc_multiallelic_gatk_oracle.py` | :24 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_pairhmm_results_oracle.py` | :36 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |
| `verify_reblock_gatk_multisample.py` | :84 (negative) | **no** | yes | skip branch -> whole test prints a skip line and exits 0 |

## 2. TASK 2 — the mechanism

### 2.1 The shared helper

`fastgatk-native/scripts/oracle_guard.py` (new, 195 lines) implements the requested policy in two
entry points. It deliberately does **not** re-implement path resolution: scripts keep computing
`root / "third_party/..."` exactly as before, and `oracle_toolchain.py` — which already resolves
the vendored JDK and jar while honouring `JAVA`/`GATK_JAR` — is untouched and still the resolver
for the scripts that use it.

* `oracle_ready(test, *paths) -> bool` — drop-in replacement for `java.exists() and jar.exists()`.
  * all paths present → `True`; the caller behaves **exactly** as before;
  * missing and `FASTGATK_REQUIRE_GATK_ORACLE` set → raises `GATKOracleMissing` (uncaught ⇒ exit 1)
    naming **every missing path**;
  * missing and the variable unset → prints the NOT-VERIFIED banner on stderr and returns `False`,
    so the script keeps its historical exit-0 path.
* `oracle_not_verified(test, *paths) -> None` — annotates an already-taken skip branch. It is a
  no-op unless one of the paths handed to it is really absent, so a branch entered because the
  *native binary* or a fixture is missing keeps its old behaviour and is never mislabelled as a
  missing oracle.

`FASTGATK_REQUIRE_GATK_ORACLE` counts as set for any value other than empty, `0`, `false`, `no`,
`off` — a superset of the pre-existing `== "1"` convention. Scripts that already carried their own
`if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1": raise SystemExit(...)` line keep it; the
helper now raises first, with the missing path named, which is strictly better evidence.

### 2.2 The two diff shapes

Shape A — one condition rewritten, nothing else changes:

```diff
-        if java.exists() and jar.exists():
+        if oracle_guard.oracle_ready('verify_filter_intervals.py', java, jar):
```

and for the boolean-variable variant (`verify_fasta_reference_tools.py`, `verify_shift_fasta.py`,
`verify_depth_of_coverage.py`):

```diff
-    java_oracle = JAVA.exists() and GATK_JAR.exists()
+    java_oracle = oracle_guard.oracle_ready('verify_fasta_reference_tools.py', JAVA, GATK_JAR)
```

Shape B — one line inserted at the top of the existing skip branch:

```diff
     if not jar.exists() or not java.exists() or not known_sites.exists():
+        oracle_guard.oracle_not_verified('verify_analyze_covariates.py', jar, java)
         if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
             raise SystemExit("missing GATK/JDK/known-sites oracle fixture")
```

plus `import oracle_guard` inserted after the last top-level import of the file. That is the entire
per-file change: **+369 / −42 lines across 176 files** — 176 `import oracle_guard` lines, 42
`oracle_ready` call sites (39 rewritten conditions + 3 boolean assignments) and 151
`oracle_not_verified` insertions; about 2 lines per file.

The generator refuses to guess. It skips any guard whose branch is not a conjunction of existence
calls (never `or`, never `any(...)`, never `all(...)`, never an `assert`), any branch that is not a
skip/exit branch, and any guard whose expressions do not resolve to a real JDK/jar path — 150 such
non-oracle guards (native binary, fixture, report-file and index checks) were correctly left
untouched. Every patched file is re-parsed with `ast.parse` before it is written, every changed
file passes `python3 -m py_compile`, and every `oracle_guard` reference was checked to lie after
its import and inside a function (no module-level use, no use-before-import).

### 2.3 Files changed

`git status --porcelain` outside the changed scripts contains exactly one entry — ['?? fastgatk-native/scripts/oracle_guard.py'] —
the new helper, which is expected to be untracked. No production file, no `CMakeLists.txt`, no root
`.md` and no branch/commit is touched.

```
176 files changed, 369 insertions(+), 42 deletions(-)
```

The complete diff of all 176 files is kept at `.diag/oracle_patch.applied.diff`; the per-guard edit
log is `.diag/oracle_patch_apply.log`.

### 2.4 What was deliberately **not** changed

* **No expected value, literal or assertion.** The patch only adds a presence check; it never edits
  an `assert`, a comparison, a fixture or a recorded value. `git diff` contains exactly 369 added
  lines (176 `import oracle_guard` + 193 helper calls) and 42 removed lines, every one of them a
  pure `A.exists()/is_file() and B.exists()/is_file()` guard expression.
* **No production code, no rebuild, no `CMakeLists.txt`, no root `.md`.**
* **Unregistered scripts were left alone** — the 7 with a real oracle guard are
  `verify_hc_alleles_deep_boundary.py`, `verify_hc_alleles_deep_limits.py`,
  `verify_hc_forced_alleles_emission_gate_oracle.py`,
  `verify_hc_gvcf_forced_allele_refblock_gatk_oracle.py`, `verify_hc_multiallelic_gatk_oracle.py`,
  `verify_pairhmm_results_oracle.py`, `verify_reblock_gatk_multisample.py`. None appears in
  `fastgatk-native/CMakeLists.txt`, so CTest never runs them and they cannot produce a green test
  on their own — there is nothing for them to falsely certify. Two of them
  (`verify_hc_alleles_deep_boundary.py`, `verify_hc_alleles_deep_limits.py`) are imported by
  registered scripts, but only as function libraries: their guards live inside `main()`, so the
  import executes nothing. Applying the same one-line edit there would be harmless; I left them
  alone to keep the blast radius to exactly what CTest executes.
* **Scripts that never guard the oracle** were untouched, including `verify_variant_recalibrator.py`
  (jar used as a reader) and the ~85 scripts with no oracle path at all.

## 3. TASK 3 — verification

### 3.1 Suite result on both backends

Command, exactly as requested:

```bash
fastgatk-native/scripts/run_regression.sh --label oracle-integrity
```

**Run 1 — the required green run**, exactly as requested:

```
# 双后端回归证据

- 时间：2026-09-11 13:42:15 CST
- 标签：oracle-integrity
- git：`243f973`（未提交变更 177 项）
- 过滤：`<全量>`　并行度：8
- ctest：`/home/turing-agents/Documents/fast-gatk/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest`
- pinned GATK 4.6.2.0：`third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`　JDK：`third_party/jdk17`

| 后端 | 构建目录 | 结果 | 通过/总数 | 耗时 |
| --- | --- | --- | --- | --- |
| omp | `OpenMP (fastgatk-native/build)` | 通过 | 294/294 | 1195.67 sec |
| serial | `Serial (fastgatk-native/build-serial)` | 通过 | 294/294 | 1198.17 sec |

日志：`20260911-132217/{omp,serial}.log`

> 说明：本块由 `fastgatk-native/scripts/run_regression.sh` 自动生成，可直接粘贴进
> `IMPLEMENTATION_STATUS.md` 的「最近验证状态」。结论只对本块记录的 git 版本与二进制有效。
```

Both backends: **294 / 294 passed**, 0 failed, ~20 minutes per backend, with the pinned GATK
4.6.2.0 jar present. **No test went red.** So no guard I changed had been masking a real GATK
comparison failure: the guards were only ever hiding the *absence of the oracle*, and the oracle is
present here. Had any test flipped red it would have meant that script had a comparison that no
longer matched GATK and had been getting away with it while green.

**Run 2 — the same suite with `FASTGATK_REQUIRE_GATK_ORACLE=1`** (the second and last of the two
full-suite runs this task allows). This is the decisive runtime measurement of §1.3: with the
oracle *required*, any registered script that cannot actually find its JDK or jar must go red.

```
# 双后端回归证据

- 时间：2026-09-11 14:02:44 CST
- 标签：oracle-integrity-require-oracle
- git：`243f973`（未提交变更 177 项）
- 过滤：`<全量>`　并行度：8
- ctest：`/home/turing-agents/Documents/fast-gatk/third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/ctest`
- pinned GATK 4.6.2.0：`third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar`　JDK：`third_party/jdk17`

| 后端 | 构建目录 | 结果 | 通过/总数 | 耗时 |
| --- | --- | --- | --- | --- |
| omp | `OpenMP (fastgatk-native/build)` | 通过 | 294/294 | 1223.09 sec |
| serial | `Serial (fastgatk-native/build-serial)` | 通过 | 294/294 | 1223.14 sec |

日志：`20260911-134221/{omp,serial}.log`

> 说明：本块由 `fastgatk-native/scripts/run_regression.sh` 自动生成，可直接粘贴进
> `IMPLEMENTATION_STATUS.md` 的「最近验证状态」。结论只对本块记录的 git 版本与二进制有效。
```

Also **294 / 294 on both backends, 0 failures**. Two conclusions, both measured rather than
inferred:

1. **No registered test is running without a real GATK comparison.** If any of the 187 affected
   tests had been unable to resolve its oracle path — a wrong `parents[N]`, a stale `JAVA`, a
   renamed jar — it would have failed loudly on this run. None did. The "0" of §1.3 is confirmed
   at runtime under fail-closed conditions, not just statically.
2. **`FASTGATK_REQUIRE_GATK_ORACLE=1` does not break the suite**, so it is usable today as a CI
   gate for this property. That is the cheapest way to make this round's guarantee permanent;
   I did not add it to the runner because that would mean editing files outside the allowed scope.


### 3.2 Proof that the new behaviour bites (one script, both directions — and the same check through CTest)

Target: **`verify_filter_intervals.py`** — registered, shape A, two guarded comparison sites.
The oracle is hidden by running the script from a `tempfile.TemporaryDirectory` that mirrors the
repo root: `fastgatk-native/build` and `third_party/jdk17` are symlinked in and
`third_party/gatk-package/gatk-4.6.2.0/` is created **empty**, so the script's own
`Path(__file__).resolve().parents[2] / "third_party/gatk-package/..."` genuinely does not exist.
Nothing in the working tree is renamed, moved or written; the temporary tree is removed by the
context manager. The only other manipulation is a per-invocation `JAVA=` and
`FASTGATK_REQUIRE_GATK_ORACLE=` environment variable, and no environment change outlives the
command. Cases (v)/(vi) additionally use `verify_analyze_covariates.py` (registered, shape B, the
script honours `JAVA`) with `JAVA=/nonexistent/jdk17/bin/java` — the "bogus JAVA" route from the
task.

Additionally, the same mechanism was driven **through the real CTest harness** (which also proves
that `run_regression.sh` propagates the environment into the tests), saved verbatim at
`.diag/oracle_ctest_demo.txt`:

```
### A) through the real CTest harness: JAVA=/nonexistent/jdk17/bin/java FASTGATK_REQUIRE_GATK_ORACLE=1 -> RED
$ JAVA=/nonexistent/jdk17/bin/java FASTGATK_REQUIRE_GATK_ORACLE=1 ctest --test-dir fastgatk-native/build -R '^fastgatk-analyze-covariates-contract$' --output-on-failure
[ORACLE REQUIRED] verify_analyze_covariates.py: FASTGATK_REQUIRE_GATK_ORACLE='1'
The GATK oracle is required for this test but was not found:
    missing: /nonexistent/jdk17/bin/java

The GATK comparison is the whole point of this test, so it fails loudly
instead of passing without it.  Unset the variable to fall back to the old
skip-and-pass behaviour (which prints a NOT-VERIFIED notice).
==============================================================================


0% tests passed, 1 tests failed out of 1
	 48 - fastgatk-analyze-covariates-contract (Failed)
(ctest exit code: 8)

### B) JAVA=/nonexistent/jdk17/bin/java, FASTGATK_REQUIRE_GATK_ORACLE unset, --output-on-failure (the run_regression.sh mode)
$ JAVA=/nonexistent/jdk17/bin/java ctest --test-dir fastgatk-native/build -R '^fastgatk-analyze-covariates-contract$' --output-on-failure
Test project /home/turing-agents/Documents/fast-gatk/fastgatk-native/build
    Start 48: fastgatk-analyze-covariates-contract
1/1 Test #48: fastgatk-analyze-covariates-contract ...   Passed    0.03 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.03 sec

### C) same as B but ctest -V (the output --output-on-failure discards)
$ JAVA=/nonexistent/jdk17/bin/java ctest --test-dir fastgatk-native/build -R '^fastgatk-analyze-covariates-contract$' -V
48: ==============================================================================
48: [NOT VERIFIED AGAINST GATK] verify_analyze_covariates.py
48: NOT verified against GATK 4.6.2.0.  A PASS from this script is NOT evidence
48: of GATK parity -- it only means native output matched native expectations
48: recorded earlier.
48:     missing: /nonexistent/jdk17/bin/java
48: Set FASTGATK_REQUIRE_GATK_ORACLE=1 to turn this into a hard failure.
48: ==============================================================================
48: {"status": "skipped", "reason": "GATK oracle unavailable"}
1/1 Test #48: fastgatk-analyze-covariates-contract ...   Passed    0.02 sec
```

Cases A/B/C drive that whole chain through the harness: **(A)** with the variable set the CTest
test is genuinely `Failed` (ctest exit code 8) and prints the `[ORACLE REQUIRED]` banner naming
the path; because A and B run the *same* test with the *same* bogus `JAVA` and differ only in
`FASTGATK_REQUIRE_GATK_ORACLE`, A also proves that `run_regression.sh` propagates the environment
into the tests — which is what makes run 2 in §3.1 meaningful. **(B)** without the variable the
same test `Passed` in exactly the mode `run_regression.sh` uses. **(C)** the NOT-VERIFIED banner
is there the moment you ask for it with `-V`. Case B is the residual gap of §4: the warning exists
but `--output-on-failure` discards it, so the durable fix is the env-var gate, not the banner.

Literal output of the six direct cases, verbatim (also saved at `.diag/oracle_bite_demo.txt`):

```
==============================================================================
### (i) POSITIVE-GUARD SCRIPT, oracle hidden, FASTGATK_REQUIRE_GATK_ORACLE=1 -> loud failure
$ FASTGATK_REQUIRE_GATK_ORACLE=1 /usr/bin/python3 /tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/verify_filter_intervals.py
--- exit code: 1
--- stdout:
(empty)
--- stderr:
Traceback (most recent call last):
  File "/tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/verify_filter_intervals.py", line 307, in <module>
    raise SystemExit(main())
                     ~~~~^^
  File "/tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/verify_filter_intervals.py", line 178, in main
    if oracle_guard.oracle_ready('verify_filter_intervals.py', java, jar):
       ~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/oracle_guard.py", line 168, in oracle_ready
    _raise(test, reported)
    ~~~~~~^^^^^^^^^^^^^^^^
  File "/tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/oracle_guard.py", line 111, in _raise
    raise GATKOracleMissing(
    ...<13 lines>...
    )
oracle_guard.GATKOracleMissing: ==============================================================================
[ORACLE REQUIRED] verify_filter_intervals.py: FASTGATK_REQUIRE_GATK_ORACLE='1'
The GATK oracle is required for this test but was not found:
    missing: /tmp/fastgatk-oracle-bite-lgtmq2ud/repo/third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar

The GATK comparison is the whole point of this test, so it fails loudly
instead of passing without it.  Unset the variable to fall back to the old
skip-and-pass behaviour (which prints a NOT-VERIFIED notice).
==============================================================================

==============================================================================
### (ii) POSITIVE-GUARD SCRIPT, oracle hidden, env var unset -> NOT-VERIFIED notice, exit 0
$ FASTGATK_REQUIRE_GATK_ORACLE=<unset> /usr/bin/python3 /tmp/fastgatk-oracle-bite-lgtmq2ud/repo/fastgatk-native/scripts/verify_filter_intervals.py
--- exit code: 0
--- stdout:
{"count_input_intervals": 6, "count_java_oracle": false, "count_output_intervals": 5, "dictionary_order_union": true, "hdf5_output_intervals": 2, "input_intervals": 5, "java_oracle": false, "merge_validation": true, "output_intervals": 2, "padding_rejection": true, "set_rule_validation": true, "short_aliases": true, "status": "pass"}
--- stderr:
==============================================================================
[NOT VERIFIED AGAINST GATK] verify_filter_intervals.py
The GATK oracle is absent, so the expectations asserted by this test were
NOT verified against GATK 4.6.2.0.  A PASS from this script is NOT evidence
of GATK parity -- it only means native output matched native expectations
recorded earlier.

    missing: /tmp/fastgatk-oracle-bite-lgtmq2ud/repo/third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar

Set FASTGATK_REQUIRE_GATK_ORACLE=1 to turn this into a hard failure.
==============================================================================

==============================================================================
### (iii) POSITIVE-GUARD SCRIPT, oracle present -> unchanged behaviour, silent, exit 0
$ FASTGATK_REQUIRE_GATK_ORACLE=<unset> /usr/bin/python3 /home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_filter_intervals.py
--- exit code: 0
--- stdout:
{"count_input_intervals": 6, "count_java_oracle": true, "count_output_intervals": 5, "dictionary_order_union": true, "hdf5_output_intervals": 2, "input_intervals": 5, "java_oracle": true, "merge_validation": true, "output_intervals": 2, "padding_rejection": true, "set_rule_validation": true, "short_aliases": true, "status": "pass"}
--- stderr:
(empty)

==============================================================================
### (iv) POSITIVE-GUARD SCRIPT, oracle present + FASTGATK_REQUIRE_GATK_ORACLE=1 -> still green
$ FASTGATK_REQUIRE_GATK_ORACLE=1 /usr/bin/python3 /home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_filter_intervals.py
--- exit code: 0
--- stdout:
{"count_input_intervals": 6, "count_java_oracle": true, "count_output_intervals": 5, "dictionary_order_union": true, "hdf5_output_intervals": 2, "input_intervals": 5, "java_oracle": true, "merge_validation": true, "output_intervals": 2, "padding_rejection": true, "set_rule_validation": true, "short_aliases": true, "status": "pass"}
--- stderr:
(empty)

==============================================================================
### (v) SKIP-GUARD SCRIPT (verify_analyze_covariates.py), JAVA=/nonexistent/jdk17/bin/java, env var unset -> notice, exit 0
$ FASTGATK_REQUIRE_GATK_ORACLE=<unset> /usr/bin/python3 /home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_analyze_covariates.py
--- exit code: 0
--- stdout:
{"status": "skipped", "reason": "GATK oracle unavailable"}
--- stderr:
==============================================================================
[NOT VERIFIED AGAINST GATK] verify_analyze_covariates.py
The GATK oracle is absent, so the expectations asserted by this test were
NOT verified against GATK 4.6.2.0.  A PASS from this script is NOT evidence
of GATK parity -- it only means native output matched native expectations
recorded earlier.

    missing: /nonexistent/jdk17/bin/java

Set FASTGATK_REQUIRE_GATK_ORACLE=1 to turn this into a hard failure.
==============================================================================

==============================================================================
### (vi) SKIP-GUARD SCRIPT, JAVA=/nonexistent/jdk17/bin/java, FASTGATK_REQUIRE_GATK_ORACLE=1 -> loud failure
$ FASTGATK_REQUIRE_GATK_ORACLE=1 /usr/bin/python3 /home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_analyze_covariates.py
--- exit code: 1
--- stdout:
(empty)
--- stderr:
Traceback (most recent call last):
  File "/home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_analyze_covariates.py", line 139, in <module>
    raise SystemExit(main())
                     ~~~~^^
  File "/home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/verify_analyze_covariates.py", line 34, in main
    oracle_guard.oracle_not_verified('verify_analyze_covariates.py', jar, java)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/oracle_guard.py", line 184, in oracle_not_verified
    _raise(test, reported)
    ~~~~~~^^^^^^^^^^^^^^^^
  File "/home/turing-agents/Documents/fast-gatk/fastgatk-native/scripts/oracle_guard.py", line 111, in _raise
    raise GATKOracleMissing(
    ...<13 lines>...
    )
oracle_guard.GATKOracleMissing: ==============================================================================
[ORACLE REQUIRED] verify_analyze_covariates.py: FASTGATK_REQUIRE_GATK_ORACLE='1'
The GATK oracle is required for this test but was not found:
    missing: /nonexistent/jdk17/bin/java

The GATK comparison is the whole point of this test, so it fails loudly
instead of passing without it.  Unset the variable to fall back to the old
skip-and-pass behaviour (which prints a NOT-VERIFIED notice).
==============================================================================
```

Reading of the six cases:

* **(i)** oracle hidden + `FASTGATK_REQUIRE_GATK_ORACLE=1` → **exit 1**, loud traceback ending in
  `GATKOracleMissing` and naming the missing jar path. The GATK comparison is the point of the test,
  so it refuses to pass.
* **(ii)** oracle hidden, variable unset → **exit 0** with the unmistakable
  `[NOT VERIFIED AGAINST GATK]` banner on stderr. Note the stdout line still says
  `"status": "pass", "java_oracle": false` — which is exactly the misleading green this round
  removes; the banner is now the counterweight.
* **(iii)** and **(iv)** oracle present (real repo) → **exit 0**, `"java_oracle": true`, empty
  stderr, byte-identical output with and without the env var. The change is inert when the oracle
  is there.
* **(v)** shape B with a bogus `JAVA`, variable unset → **exit 0**, banner names
  `/nonexistent/jdk17/bin/java`, alongside the script's own
  `{"status": "skipped", "reason": "GATK oracle unavailable"}`.
* **(vi)** same, variable set → **exit 1**, loud failure naming the same path.

## 4. What remains unproven

* **Not proven: that the guarded comparisons are exhaustive.** This round proves a missing oracle
  is now loud. It says nothing about whether the comparisons inside those guards cover all the
  literals their scripts assert — the audit's §4.4 findings (native-only literals, positional
  FORMAT reads, `SelectVariants` used as a reader) are untouched.
* **Not proven: perfect oracle coverage of every comparison site.** The `FASTGATK_REQUIRE_GATK_ORACLE=1`
  run proves that every registered script *that reaches a patched guard* found its oracle. A script
  whose oracle guard sits behind an earlier `return` (e.g. a missing native binary) never reaches
  one. That is the intended contract of `oracle_not_verified` — it will not libel a missing native
  binary as a missing oracle — but it means the env-var run is not proof that every comparison in
  the corpus executed.
* **Not proven: that no other silent-skip mechanism exists.** I scanned for existence guards on
  java/jar. A script could still swallow a GATK failure another way (a bare `except`, a swallowed
  non-zero return code, `if result.returncode != 0: pass`). I did not audit for that.
* **Not proven: the notice is visible in the regression log.** `run_regression.sh` drives CTest
  with `--output-on-failure`, which discards the output of *passing* tests. On the pass path the
  NOT-VERIFIED banner therefore lands in a captured-then-dropped stream: it is visible when the
  script is run by hand or with `ctest -V -R <name>`, not in `omp.log`/`serial.log`. Fixing that
  needs a change to `CMakeLists.txt` or the runner, both out of scope for this round. **This is the
  one real residual gap in "no reader can mistake it for a verified parity pass"**, and the
  cheapest follow-up is a CTest `PASS_REGULAR_EXPRESSION`/`FAIL_REGULAR_EXPRESSION` or a
  `FASTGATK_REQUIRE_GATK_ORACLE=1` CI gate.
* **Speculation (labelled as such):** I believe the 125 scripts that already carried a bare
  `if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1": raise SystemExit(...)` were added
  reactively, one audit round at a time, rather than from a single policy — the wording and the
  `== "1"` comparison are inconsistent across files, and a few use a bare truthiness test instead.
  I have no commit-history evidence for that.
* **Speculation:** the 7 unregistered guarded scripts are probably exercised by hand during
  investigation rounds; nothing in the repository enforces that.

## 5. Artefacts produced

| path | content |
| --- | --- |
| `fastgatk-native/scripts/oracle_guard.py` | the new shared helper (the only new tracked-to-be file) |
| `.diag/round-oracle-integrity.md` | this report |
| `.diag/oracle_patch.applied.diff` | complete diff of all 176 changed scripts |
| `.diag/oracle_patch_apply.log` | per-guard edit log (every rewritten condition and inserted line) |
| `.diag/oracle_patch.diff`, `.diag/oracle_patch.log` | the dry run that was reviewed before applying |
| `.diag/oracle_table.json` | enriched per-script scan (guard lines, shapes, resolved paths) |
| `.diag/oracle_scan.json`, `.diag/oracle_scan.err` | first, text-oriented scan pass |
| `.diag/oracle_exposure_table.md` | the table embedded in §1.5 |
| `.diag/oracle_stats.json` | the counts quoted in §0 |
| `.diag/oracle_bite_demo.txt` | verbatim output of the six direct bite-demo cases |
| `.diag/oracle_ctest_demo.txt` | verbatim output of the three cases driven through the real CTest harness |
| `.diag/suite_run1_evidence.md`, `.diag/suite_run2_evidence.md` | the two `run_regression.sh` evidence blocks |
| `.diag/oracle_guard_scan.py`, `.diag/oracle_guard_patch.py`, `.diag/build_report.py` | the scan/patcher/report tooling |
| `.diag/regression/20260911-132217/` | run 1 (green) logs and evidence block |
| `.diag/regression/20260911-134221/` | run 2 (`FASTGATK_REQUIRE_GATK_ORACLE=1`) logs and evidence block |
| `.diag/ctest_list.json` | `ctest -N --show-only=json-v1`, source of the "187 of 294" figure |

Nothing under `/tmp` was left behind: every experiment ran inside a `tempfile.TemporaryDirectory`.
