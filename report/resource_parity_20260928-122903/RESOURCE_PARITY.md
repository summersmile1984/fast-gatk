# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [8.8115234375, 5.8916015625, 4.33203125] -> [5.609375, 6.07177734375, 4.7001953125]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.60 | 6.22 | 10.33× | 0.51 | 3.62 | 7.11× | 146.12 | 388.86 | 109.47 | 51.46 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 3.19 | 23.97 | n/a | 3.17 | 20.62 | n/a | 142.00 | 540.29 | 89.99 | 76.45 | 28.79 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.30 | 5.68 | 19.14× | 0.28 | 3.01 | 10.75× | 32.57 | 392.75 | 1.24 | 50.90 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 6.25 | 619.29× | 0.00 | 2.90 | n/a | 1.79 | 392.80 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.52 | 4.23 | 8.07× | 0.82 | 5.06 | 6.13× | 158.04 | 339.48 | 109.51 | 51.45 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.81 | 17.28 | n/a | 3.77 | 23.38 | n/a | 142.18 | 757.51 | 89.99 | 76.46 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.25 | 4.04 | 16.15× | 1.12 | 4.07 | 3.62× | 24.89 | 273.54 | 1.23 | 50.92 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.23 | 4.72 | 20.83× | 5.32 | 4.23 | 0.80× | 10.40 | 335.32 | 0.02 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
