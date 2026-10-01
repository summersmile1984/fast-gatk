# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [5.416015625, 5.689453125, 5.76171875] -> [6.93505859375, 8.48974609375, 7.46826171875]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.53 | 8.03 | 15.06× | 0.50 | 3.84 | 7.68× | 158.55 | 388.89 | 109.48 | 51.46 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.33 | 6.40 | 19.18× | 0.32 | 3.58 | 11.38× | 153.59 | 363.63 | 109.48 | 51.44 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 3.39 | 25.21 | n/a | 3.37 | 21.74 | n/a | 142.21 | 529.97 | 89.99 | 76.47 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.29 | 5.76 | 19.64× | 0.28 | 3.17 | 11.51× | 32.55 | 367.95 | 1.24 | 50.93 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 14.26 | 22.16 | 1.55× | 14.23 | 19.59 | 1.38× | 717.58 | 652.71 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.02 | 5.44 | 282.95× | 0.00 | 2.72 | n/a | 0.93 | 323.66 | 0.00 | 48.06 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 6.58 | 533.04× | 0.00 | 3.92 | n/a | 1.60 | 682.46 | 0.00 | 48.21 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.97 | 11.63 | 3.91× | 2.94 | 8.40 | 2.86× | 97.07 | 1156.52 | 51.68 | 99.87 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 6.81 | 512.11× | 0.00 | 3.06 | n/a | 1.05 | 389.75 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.54 | 4.85 | 9.03× | 0.77 | 6.07 | 7.83× | 151.74 | 320.71 | 109.51 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.34 | 6.84 | 20.16× | 0.64 | 5.58 | 8.72× | 150.01 | 372.14 | 109.47 | 51.45 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.84 | 17.38 | n/a | 3.78 | 23.80 | n/a | 142.53 | 773.23 | 89.99 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.25 | 4.09 | 16.55× | 1.11 | 4.35 | 3.91× | 25.80 | 318.94 | 1.23 | 50.90 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 12.18 | 16.95 | 1.39× | 35.70 | 34.33 | 0.96× | 717.43 | 956.84 | 37.41 | 123.15 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.07 | 4.08 | 55.12× | 0.00 | 3.72 | n/a | 7.68 | 270.16 | 0.19 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.02 | 4.01 | 207.49× | 0.00 | 4.24 | n/a | 0.68 | 547.06 | 0.00 | 48.21 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.92 | 6.99 | 2.40× | 2.91 | 15.61 | 5.37× | 97.19 | 1526.72 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.21 | 9.46 | 44.91× | 4.88 | 4.06 | 0.83× | 10.50 | 297.00 | 0.02 | 49.88 | 0.00 | 0.00 |

## Lane `full-machine-tuned`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.56 | 4.38 | 7.88× | 0.54 | 5.88 | 10.88× | 160.84 | 317.64 | 109.47 | 51.43 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.33 | 4.34 | 13.33× | 0.31 | 5.08 | 16.40× | 146.41 | 294.96 | 109.50 | 51.45 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.93 | 17.78 | n/a | 3.41 | 24.80 | n/a | 142.79 | 751.36 | 89.99 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.26 | 4.21 | 16.11× | 0.41 | 4.47 | 10.78× | 28.62 | 332.24 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 12.76 | 17.08 | 1.34× | 20.16 | 35.47 | 1.76× | 717.87 | 948.14 | 37.41 | 123.18 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.02 | 4.00 | 245.33× | 0.00 | 3.81 | n/a | 1.40 | 273.67 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 4.00 | 295.20× | 0.00 | 3.90 | n/a | 0.80 | 535.37 | 0.00 | 48.21 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.95 | 6.57 | 2.23× | 2.92 | 16.61 | 5.68× | 97.21 | 1397.55 | 51.68 | 99.88 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.03 | 6.01 | 221.51× | 0.14 | 4.58 | 32.68× | 10.32 | 346.50 | 0.02 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
