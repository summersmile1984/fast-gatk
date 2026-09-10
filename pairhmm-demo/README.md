# Kokkos PairHMM prototype

This directory now contains two implementations:

- `pairhmm-kokkos`: the primary one-source implementation. Execution,
  allocation/copies, views, and CPU SIMD all use Kokkos 5.2.0 APIs.
- `pairhmm-demo`: the earlier raw-intrinsics experiment, retained as a legacy
  comparison and not part of the portable design.

The Kokkos kernel processes independent read/haplotype pairs in
`Kokkos::Experimental::simd<double>` lanes. It recognizes contiguous and
broadcast input groups through the same portable SIMD API, updates M/I/D in
place using a three-row workspace, and launches groups with
`Kokkos::RangePolicy<DefaultExecutionSpace>`. There is no `_mm*` intrinsic or
Rust/C ABI in the Kokkos path.

## Pinned build

The repository pins Kokkos 5.2.0 and CMake 4.3.4 under `third_party/`. Build the
same source for generic scalar, AVX2, and AVX-512:

```bash
bash pairhmm-demo/scripts/build_kokkos_variants.sh
```

The architecture selection is a build property, not an `#ifdef` in the
algorithm:

| Build | Kokkos architecture | `simd<double>` width on this host |
|---|---|---:|
| `build-kokkos-scalar` | generic | 1 |
| `build-kokkos-avx2` | `ZEN3` | 4 |
| `build-kokkos-avx512` | `ZEN4` | 8 |

Example:

```bash
pairhmm-demo/build-kokkos-avx512/pairhmm-kokkos \
  --workload=matrix8 --pairs=512 --read-len=150 --hap-len=160 \
  --threads=16 --iterations=10
```

`matrix8` creates the same 8-read × 8-haplotype blocked workload used by the
GKL benchmark. Input sequences remain unique and pairs refer to them by ID, so
the matrix does not duplicate bases or qualities.

## Bit-identical verification

```bash
python3 pairhmm-demo/scripts/verify_kokkos.py
```

The default test compares all 4096 matrix values from scalar, AVX2, and
AVX-512 with the real GATK 4.6.2.0 Java `LoglessPairHMM`; all currently report
`bit_different=0`. The strict oracle disables HotSpot's optional x86
`Math.log10` intrinsic so Java uses its specified StrictMath/fdlibm path:

```text
-XX:+UnlockDiagnosticVMOptions -XX:DisableIntrinsic=_dlog10
```

This distinction is necessary: the default HotSpot x86 intrinsic is permitted
to select an adjacent double. On this 4096-value corpus, DP scaled sums are all
bit-identical, while default HotSpot changes three final likelihoods by one ULP.
The Kokkos kernel embeds the portable fdlibm operation order for a stable
cross-libm contract and is compiled with `-ffp-contract=off -fno-fast-math`.

## Benchmark

```bash
python3 pairhmm-demo/scripts/benchmark_kokkos.py \
  --independent-iterations=30 --matrix-iterations=10 --repeats=3 \
  --cores=1,2,4,8,16 --variants=avx2,avx512 \
  --output=pairhmm-demo/results/benchmark-kokkos-5.2.0-range-final.json
```

For a compact same-core Java comparison that records the actual portable API
and Kokkos execution space, use:

```bash
python3 pairhmm-demo/scripts/benchmark_same_cores.py \
  --native-backend=auto --cores=1,2,4,8,16 --pairs=4096 \
  --iterations=3 --repeats=3
```

`auto` selects the highest ISA build that the host advertises and that exists
on disk (AVX-512, then AVX2, then scalar); it will never start an unsafe ISA
binary. The command invokes `pairhmm-kokkos` and reports
`native_api=Kokkos::Experimental::simd<double>`; it no longer measures the
legacy raw-intrinsics executable. An explicit `--native-backend=avx512` or
`avx2` fails closed when the host does not advertise that ISA.

The script pins both sides to identical physical CPUs, alternates run order,
uses a complete untimed warmup, and reports medians. Independent pairs compare
with pure-Java GATK. Matrix throughput compares with GKL's actual `n*n` work,
not its misleading diagonal-record count. Kokkos kernel-only and amortized
input-preparation rates are reported separately.

On the current AMD Ryzen 9 7945HX, the recorded AVX-512 independent-pair
speedup is workload-, core-count-, and thermal-state-dependent (the checked-in
artifact is the source of truth; this short 1/2/4-core run is roughly 5–9×).
Matrix8 and GKL comparisons must be reported separately because they include a
different pair cardinality. This is a PairHMM kernel result, not yet a
whole-HaplotypeCaller or whole-GATK result.
The reproducible same-core artifact for the current host is
[`benchmark-kokkos-same-core-20260825.json`](/home/turing-agents/Documents/fast-gatk/pairhmm-demo/results/benchmark-kokkos-same-core-20260825.json); it records the selected host ISA, execution space, SIMD width and preparation cost.

## Reusable kernel target

The implementation is exported by the top-level `fastgatk-kernels` CMake target.
`pairhmm-kokkos` links that target, while the legacy headers and Makefile remain
source-compatible forwarders. The target also exports the
`fastgatk::kernels::smith_waterman_score_reference` and
`smith_waterman_align_reference` score/traceback/CIGAR APIs. Exact GATK/GKL
overhang and tie-breaking remain an oracle task. Enable the smoke test with
`-DFASTGATK_KERNELS_BUILD_TESTS=ON` and run the resulting
`fastgatk-kernels-api-smoke` executable.

A small reproducible P0 benchmark artifact is stored at
`results/benchmark-p0-kernels-smoke.json`; it records scalar, AVX2 and AVX-512
Kokkos SIMD widths plus kernel-only and amortized throughput.
