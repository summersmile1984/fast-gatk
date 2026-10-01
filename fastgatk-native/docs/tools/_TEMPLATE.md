# fastgatk-<tool-name> — tool doc template

**Copy this file to `docs/tools/<binary>.md` and fill in every section
before opening the code editor for the corresponding feature gap.**

---

# fastgatk-<tool-name>

**GATK equivalent:** `org.broadinstitute.hellbender.tools.<Tool>.java`
(URL)
**Source:** `fastgatk-native/src/<tool>_tool.cpp` (or wrapped file)
**Status:** shipped | partial | not-shipped
**Oracle(s):** `scripts/verify_<...>.py`

## 1. Purpose

One sentence: what this tool does and when a user invokes it.

## 2. CLI parity

Table or link to `docs/cli-alignment.md` section. Highlight only the
flags that differ from GATK 4.6.2.0. For each difference: native
behavior + fallback policy from `dispatcher/tool_registry.json`.

## 3. Algorithm(s) used

Bullet list of `docs/algorithms/<name>.md` linked.

## 4. Output contract

What the tool writes: file format, compression, headers, sample
ordering, etc. Pin to GATK 4.6.2.0 schema.

## 5. Oracle(s) and expected parity

For each oracle script that targets this tool:

- script path
- pinned fixture(s) (with SHA1 if pinned)
- what fields are byte-identical
- what fields are semantic-equivalent (with documented threshold)
- last-verified timestamp from `evidence/<date>-rerun/`

## 6. Known gaps

Same shape as algorithm-doc §6. Cross-reference to
`docs/algorithms/` and `work/mutect2-priority/PARITY_AUDIT.md`.

## 7. Build / install

- CMake target name (`fastgatk-<tool-name>`)
- Source dependencies (other `fastgatk-*` libs)
- Kokkos backend requirement (default OpenMP; CUDA/HIP/SYCL
  optional but not in this sprint's CI gate)
- Runtime requirements (HTSlib, etc.)

## 8. CLI quick reference

Example invocation with minimal flags. Link to `cli-alignment.md`
for full option surface.

## 9. Cross-references

Pipeline-position: which upstream tools feed it, which downstream
tools consume its output.
