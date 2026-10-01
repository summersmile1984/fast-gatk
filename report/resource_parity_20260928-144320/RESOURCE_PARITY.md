# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [6.02685546875, 6.44189453125, 5.68017578125] -> [5.16650390625, 6.55419921875, 6.05810546875]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.55 | 6.18 | 11.25× | 0.53 | 3.39 | 6.41× | 129.34 | 389.32 | 109.47 | 51.46 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.64 | 6.16 | 9.60× | 0.30 | 3.38 | 11.47× | 158.18 | 365.55 | 109.50 | 51.44 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 3.18 | 23.84 | n/a | 3.15 | 20.16 | n/a | 141.56 | 524.44 | 89.99 | 76.45 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.34 | 5.92 | 17.46× | 0.32 | 3.08 | 9.64× | 29.57 | 365.74 | 1.24 | 50.93 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 14.32 | 22.39 | 1.56× | 14.27 | 19.60 | 1.37× | 717.64 | 656.28 | 37.41 | 123.18 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.03 | 5.78 | n/a | 0.00 | 2.40 | n/a | 7.62 | 325.67 | 0.13 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 5.74 | 428.98× | 0.00 | 3.15 | n/a | 1.29 | 681.08 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.89 | 11.13 | n/a | 2.87 | 8.33 | n/a | 101.11 | 1156.77 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 6.80 | 592.89× | 0.00 | 3.04 | n/a | 0.87 | 391.66 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.53 | 4.59 | 8.63× | 0.79 | 5.38 | 6.77× | 145.85 | 356.43 | 109.51 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.32 | 4.18 | 13.11× | 0.56 | 4.75 | 8.48× | 159.63 | 300.57 | 109.50 | 51.45 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.77 | 17.17 | n/a | 3.70 | 22.73 | n/a | 142.59 | 749.90 | 89.99 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.26 | 4.15 | 16.08× | 1.06 | 3.88 | 3.64× | 30.29 | 270.03 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 12.07 | 19.60 | 1.62× | 35.38 | 36.36 | 1.03× | 717.90 | 944.16 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.02 | 4.47 | n/a | 0.00 | 3.73 | n/a | 0.61 | 273.19 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 4.02 | 365.89× | 0.00 | 3.96 | n/a | 0.73 | 542.92 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.84 | 6.67 | n/a | 2.82 | 15.11 | n/a | 101.25 | 1367.06 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.20 | 5.17 | 25.29× | 4.93 | 4.48 | 0.91× | 10.52 | 331.75 | 0.02 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
