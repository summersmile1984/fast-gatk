# Resource parity workflow

Second verification axis alongside the GATK output-parity oracles: given the
same input, the same functional output and the **same resource envelope**,
how fast is each implementation and what does it consume (CPU, memory, disk
I/O)?  The goal is a *more efficient* implementation under functional
equality — so functional equality is asserted before any ratio is reported.

## Protocol (equal-resource-v1)

* Same fixture and functionally identical CLI per case.  The only permitted
  differences are the implementation entry points and the recorded fixed JVM
  heap limits (`-Xmx1g`/`-Xmx2g`, per case).
* Both sides of a lane receive the **same** resource envelope, applied
  outside the CLI (no per-tool flags leak in):

  | Lane | Envelope | Answers |
  | --- | --- | --- |
  | `1-core` | `taskset -c 0`, `OMP_NUM_THREADS=1` | per-core algorithmic efficiency |
  | `full-machine` | all cores, `OMP_NUM_THREADS=nproc` | as-shipped throughput in the same envelope |

  GATK's inability to exploit extra cores is its own limitation, not an
  unequal envelope.
* Runs are strictly serialized (one workload at a time), warm cache,
  `--repeat N` (default 2), p50/min/max reported.
* Paired, binary-identified runs: SHA-256 of each native binary exercised and
  of the GATK jar, java version, git HEAD, CPU model, nproc and loadavg
  before/after — per the evidence doctrine in
  [regression-evidence.md](regression-evidence.md).
* **Correctness first**: each case runs the byte/record equivalence checks of
  `run_headline_comparison.py`; speedup ratios appear only on output-
  equivalent rows.  Unlike workloads never justify a speedup claim.

## Measurement

GNU `time -v` root statistics (wall, user+sys CPU, %CPU, max RSS, filesystem
block I/O) plus a 25 ms process-tree sampler (walks
`/proc/<pid>/task/<tid>/children`, aggregates every descendant, native and
java alike):

| Field | Meaning |
| --- | --- |
| `wall_p50_s` | throughput |
| `cpu_p50_s` | user+sys CPU (wait4 tree rollup) — the equal-resource efficiency measure |
| `tree_rss_peak_kb` | peak of the summed per-process RSS over the tree |
| `tree_io_read/write_bytes` | physical disk bytes (`read_bytes`/`write_bytes`); near zero on warm cache |
| `tree_rchar/wchar_bytes` | bytes through read/write syscalls (cache-inclusive; mmap reads are not counted) |
| `fs_in/out_blocks` | GNU time `%I/%O`, kept for continuity with historical records |

Sampler limits are the same as the validated
`work/mutect2-priority/sample_tree` harness: 25 ms cadence can miss very
short-lived helpers; counters are per-process cumulative maxima, not cgroup
accounting.

## Reproducing

```bash
# full suite (9 cases × 2 lanes; background-friendly)
python3 fastgatk-native/scripts/resource_parity_report.py --repeat 2

# focused run
python3 fastgatk-native/scripts/resource_parity_report.py --cases hc bqsr --lanes 1-core
```

Output lands in `report/resource_parity_<ts>/` with
`RESOURCE_PARITY.md` (human tables, one per lane), `summary.json`
(per-run raw data + ratios), `work/<case>/` (raw artifacts), and the
`report/resource_parity_latest` symlink.

Cases are the nine definitions of `run_headline_comparison.py`: HC (chr20
mnp real), HC aligned-lineage variant, BQSR (CEUTrio chr20), Mutect2
(1000G low-coverage chr17 window), Mutect2 DREAM chr20 tumor+normal,
SortSam, MarkDuplicates, MarkDuplicates CEUTrio chr20, GenotypeGVCFs.

## OpenMP spin budget (measured frontier)

libgomp's default spin count trades CPU for barrier latency.  On the DREAM
somatic workload (barrier-heavy, ~2.9 cores effective):

| Spin budget | Wall (s) | CPU (s) | Note |
| --- | --- | --- | --- |
| default (spin ~3e5) | 11.6 | 33.1 | wall-optimal |
| `GOMP_SPINCOUNT=100000` | 12.3 | 29.0 | |
| `GOMP_SPINCOUNT=5000` | 12.9 | 18.9 | launcher default: −43% CPU for +10% wall |
| `OMP_WAIT_POLICY=passive` | 13.7 | 14.7 | CPU-minimal |
| `GOMP_SPINCOUNT=0` | 15.0 | 15.5 | dominated by passive |

Thread count (8 vs 32) barely moves the wall: the workload is barrier-latency
bound, so the spin policy is the effective knob.  Output is byte-identical
across every configuration.

libgomp parses `GOMP_SPINCOUNT` at startup (`setenv` after `main` begins is
NOT equivalent — verified interleaved: in-binary setenv left 6-36× more CPU),
so the budget must be applied before exec.  The dispatcher
(`dispatcher/fastgatk.py`) therefore sets `GOMP_SPINCOUNT=5000` when the user
has not set `GOMP_SPINCOUNT`/`OMP_WAIT_POLICY`; an explicit setting always
wins.  The `full-machine-tuned` lane in `resource_parity_report.py` measures
this shipped configuration.

## Reading the tables

* Wall speedup = GATK p50 / native p50 (>1 ⇒ native faster).
* CPU ratio is the equal-resource efficiency measure: same work, same
  envelope.  Prefer it for "is the implementation more efficient"; prefer
  wall for practical throughput under a given envelope.
* `rchar` differences can reflect I/O strategy (htslib `read()` vs JVM mmap),
  not just volume; physical disk I/O is the `read_bytes`/`write_bytes` pair.
