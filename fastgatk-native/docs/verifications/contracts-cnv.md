# CNV verification contracts index

Each row points at oracle scripts (`scripts/verify_*.py`) that pin a
specific contract of the copy-number chain against GATK 4.6.2.0, in the
style of `contracts-hc.md` / `contracts-mutect.md`.  Status reflects the
2026-09-28 full `ctest` run (370/371 green; the single failure is
`fastgatk-fixture-presence-check`, an unmaterialized-corpus environment
gate unrelated to these tools).  Rerun evidence lives under
`evidence/<date>-rerun/<tool>/`.

Covered tools: **CollectReadCounts, DenoiseReadCounts,
CreateReadCountPanelOfNormals, ModelSegments, CallCopyRatioSegments**.
The germline gCNV trio (DetermineGermlineContigPloidy / GermlineCNVCaller /
PostprocessGermlineCNVCalls) is a separate, not-yet-ported family whose
contracts and fixture-pin strategy are specified in `local://gcnv-plan.md`.

## CollectReadCounts

| Oracle | Contract | Status |
|---|---|---|
| `verify_collect_read_counts.py` | Native TSV/interval output contract (columns, interval identity, determinism) | pass |
| `verify_collect_read_counts_gatk_oracle.py` | Native TSV boundary byte-compared against pinned GATK output | pass |

## DenoiseReadCounts

| Oracle | Contract | Status |
|---|---|---|
| `verify_denoise_read_counts.py` | Native TSV denoised-copy-ratio path contract | pass |
| `verify_denoise_read_counts_hdf5_metadata_gatk_oracle.py` | `SimpleCountCollection` stores its SAM sequence dictionary in HDF5 metadata; metadata propagation into the denoised output must match GATK | pass |
| `verify_denoise_read_counts_integer_input_gatk_oracle.py` | COUNT decoding uses `DataLine.getInt`; integer-count input boundaries match GATK exactly | pass |
| `verify_denoise_read_counts_interval_identity_gatk_oracle.py` | GATK 4.6.2.0 requires the case interval list to be identical to the PoN's; identity violation behavior matches | pass |

## CreateReadCountPanelOfNormals

| Oracle | Contract | Status |
|---|---|---|
| `verify_create_read_count_panel_of_normals.py` | Native GATK-readable HDF5 PoN writer contract (v7 layout, Java-reader round-trip) | pass |
| `verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py` | Degenerate-SVD boundary: GATK refuses to write a multi-sample PoN when SVD is requested but degenerate; native matches the refusal | pass |
| `verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py` | v7 PoN carries both original and filtered panel sample filenames; metadata bytes match GATK | pass |

## ModelSegments (bounded parity)

ModelSegments' stochastic surfaces cannot be byte-compared against GATK's
random draws; the oracles follow the bounded-parity convention
(Pinned-vs-Recorded: deterministic slices byte-equal, random slices
tolerance-pinned with a determinism self-check).

| Oracle | Contract | Status |
|---|---|---|
| `verify_model_segments.py` | Native TSV segmentation boundary contract | pass |
| `verify_model_segments_mcmc_oracle.py` | MCMC determinism (two runs byte-equal) + convergence pins at the documented tolerance | pass |
| `verify_model_segments_default_kernel_gatk_oracle.py` | Default `KernelSegmenter` dispatch when `--segments` omitted matches GATK's dispatch | pass |
| `verify_model_segments_input_segments_gatk_oracle.py` | `--segments` partition input: GATK skips kernel segmentation when a Picard interval list is supplied; native matches the partitioned path | pass |
| `verify_model_segments_copy_ratio_conditionals_gatk_oracle.py` | `CopyRatioModeller` conditional-model slice pinned (native does not claim Java's release-specific MCMC trajectory) | pass |
| `verify_model_segments_allele_fraction_initialization_gatk_oracle.py` | `AlleleFractionInitializer` alt-minor/ref-minor integration pinned | pass |
| `verify_model_segments_allele_fraction_likelihood_gatk_oracle.py` | Count-backed allele-fraction likelihood pinned (bounded oracle) | pass |
| `verify_model_segments_first_alt_fraction_gatk_oracle.py` | Combined allelic segmentation uses the first alt-allele fraction of each site, matching `MultisampleMultidimensionalKernelSegmenter` | pass |
| `verify_model_segments_multisample_gatk_oracle.py` | Multisample vector-kernel segmentation behavior pinned | pass |
| `verify_model_segments_smoothing_gatk_oracle.py` | `MultidimensionalModeller` credible-interval segment smoothing/merging pinned | pass |

## CallCopyRatioSegments

| Oracle | Contract | Status |
|---|---|---|
| `verify_call_copy_ratio_segments.py` | Native arithmetic path contract (calling caps, rounding) | pass |
| `verify_call_copy_ratio_segments_gatk_oracle.py` | Byte comparison against pinned GATK segment calls | pass |
| `verify_call_copy_ratio_segments_compensated_sum_gatk_oracle.py` | Commons-Math compensated summation parity in the segment arithmetic | pass |
| `verify_call_copy_ratio_segments_interval_validation_gatk_oracle.py` | Interval validation boundary matches GATK | pass |
| `verify_call_copy_ratio_segments_nonfinite_gatk_oracle.py` | Non-finite input handling boundary matches GATK | pass |

## Convention notes

- **Pinned-vs-Recorded** (see `verify_model_segments_mcmc_oracle.py`
  docstring): deterministic outputs must be byte-identical across runs;
  stochastic outputs are recorded once and compared at a documented
  tolerance, with a determinism self-check on the same invocation.
- **Fixture-pin** (see `local://gcnv-plan.md`): comparisons against the
  checked-in gold corpora
  (`gatk-source/src/test/resources/.../copynumber/gcnv-postprocess/`,
  `gcnv-sim-data/`, `gcnv-numerical-accuracy/`) do not require a pinned
  gcnvkernel Python environment; the Java-run `*_gatk_oracle.py` route is
  a separate, explicitly out-of-scope project.
