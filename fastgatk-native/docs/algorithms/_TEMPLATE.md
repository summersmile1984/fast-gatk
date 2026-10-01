# <algorithm-name> — algorithm doc template

**Copy this file to `docs/algorithms/<name>.md` and fill in every section
before opening the code editor for the corresponding feature gap.**

---

# <algorithm name>

**Source:** `fastgatk-kernels/<relative-path>.cpp` (or `fastgatk-native/src/...`)
**GATK reference:** `https://github.com/broadinstitute/gatk/blob/4.6.2.0/<gatk-path>.java`
**Status:** shipped | partial | not-shipped
**Tested by oracle:** `scripts/verify_<...>.py` (link in `docs/verifications/contracts-*.md`)

## 1. Purpose

What the algorithm computes. The GATK Java class name. The math
formula in TeX if non-trivial. One-paragraph summary.

## 2. Input / output contract

Inputs:

- type
- shape
- NaN / null conventions
- edge cases (empty input, single element, etc.)

Outputs:

- type
- shape
- index / key conventions

## 3. Reference implementation

Link to GATK source. Excerpt of the Java method (or the math
reference, e.g. Wikipedia for Smith-Waterman). Document the **exact
numerical policy**: log10 vs natural, FTZ, exact-zero substitution,
denormal handling, etc. Anything the oracle must verify byte-equal.

## 4. Native implementation

Architecture summary: which Kokkos kernels, where on Host/device,
whether a chunker exists, whether slab release happens. Algorithmic
differences from GATK (must be explicit, with file:line references).
GPU-safety analysis: any host-only intrinsics, any pure-host
codepaths, any divergence between OpenMP / CUDA / HIP / SYCL.

## 5. Oracle contract

The byte-level (or semantic-equivalent) fields the oracle verifies.
Cross-reference: `docs/verifications/contracts-*.md`.

## 6. Known parity gaps

For each open gap:

- gap description
- GATK source path
- native state (shipped / partial / not-shipped)
- blocking oracle (the next step)
- priority tier (1 / 2 / 3 / 4 / 5)

## 7. Reference commands

CLI invocations for both GATK Java and native that exercise this
algorithm. Output artifact paths. Command-line flags relevant to the
algorithm.

Example:

```bash
# GATK Java
java -jar gatk-package-4.6.2.0-local.jar <Tool> \
    --input fixture.bam --reference ref.fasta -O out.vcf
# Native
fastgatk-native/build/fastgatk-<tool> \
    --input fixture.bam --reference ref.fasta -O out.vcf
```

## 8. Verification anchors

The shell + python one-liners that prove the contract. Oracle
script path. Manual reproduction snippet for offline debugging.

## 9. Cross-references

Tool docs that depend on this algorithm. Other algorithm docs that
this one depends on.
