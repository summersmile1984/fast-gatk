# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [2.5830078125, 2.8583984375, 2.81396484375] -> [3.04833984375, 2.81494140625, 2.83837890625]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.46 | 5.45 | 11.85× | 0.44 | 2.79 | 6.34× | 154.57 | 388.37 | 109.47 | 51.43 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.36 | 5.41 | 15.10× | 0.33 | 2.83 | 8.43× | 138.72 | 366.37 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 43.45 | 21.02 | 0.48× | 43.38 | 17.77 | 0.41× | 398.70 | 523.39 | 95.26 | 76.47 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.28 | 8.29 | 29.68× | 0.25 | 2.46 | 9.86× | 30.11 | 390.11 | 1.24 | 50.89 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 13.32 | 19.60 | 1.47× | 13.29 | 17.05 | 1.28× | 717.62 | 653.89 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.01 | 5.83 | 399.37× | 0.00 | 2.01 | n/a | 1.05 | 323.15 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.02 | 5.16 | 279.05× | 0.00 | 2.56 | n/a | 1.60 | 685.85 | 0.00 | 48.19 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.76 | 9.34 | 3.38× | 2.73 | 6.68 | 2.44× | 97.07 | 1154.43 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 5.82 | 605.08× | 0.00 | 2.36 | n/a | 1.19 | 391.83 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.46 | 4.10 | 8.86× | 0.58 | 4.51 | 7.78× | 144.96 | 381.99 | 109.51 | 51.43 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.29 | 3.93 | 13.63× | 0.49 | 3.93 | 7.94× | 132.50 | 298.86 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 43.00 | 15.95 | 0.37× | 49.22 | 20.77 | 0.42× | 399.48 | 759.26 | 95.26 | 76.46 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.22 | 5.06 | 22.51× | 0.89 | 3.62 | 4.10× | 24.87 | 304.15 | 1.23 | 50.89 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 10.80 | 14.76 | 1.37× | 31.95 | 28.27 | 0.89× | 717.85 | 960.20 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.01 | 3.90 | 269.98× | 0.00 | 2.70 | n/a | 0.93 | 247.46 | 0.00 | 48.04 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 3.80 | 341.90× | 0.00 | 3.34 | n/a | 1.11 | 524.13 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.77 | 5.90 | 2.13× | 2.74 | 12.92 | 4.71× | 97.21 | 1393.76 | 51.68 | 99.88 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.18 | 5.45 | 29.86× | 5.13 | 3.49 | 0.68× | 10.39 | 309.40 | 0.02 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine-tuned`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.48 | 4.06 | 8.45× | 0.47 | 4.53 | 9.65× | 144.54 | 325.47 | 109.47 | 51.43 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.28 | 4.02 | 14.59× | 0.27 | 4.24 | 15.98× | 156.64 | 340.57 | 109.47 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 41.06 | 15.61 | 0.38× | 46.86 | 20.20 | 0.43× | 400.38 | 731.84 | 95.26 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.22 | 4.61 | 20.63× | 0.38 | 3.54 | 9.33× | 23.12 | 310.74 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 11.26 | 18.02 | 1.60× | 18.05 | 28.86 | 1.60× | 717.75 | 932.30 | 37.41 | 123.15 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.01 | 3.68 | 245.90× | 0.00 | 2.99 | n/a | 0.93 | 257.54 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 3.72 | 372.38× | 0.00 | 3.19 | n/a | 0.55 | 521.07 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.80 | 5.93 | 2.12× | 2.78 | 12.54 | 4.51× | 97.10 | 1404.87 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.02 | 4.57 | 238.92× | 0.14 | 3.58 | 26.48× | 0.93 | 331.16 | 0.00 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
