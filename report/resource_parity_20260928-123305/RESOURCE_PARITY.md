# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [4.64892578125, 5.8125, 4.708984375] -> [8.56689453125, 7.77001953125, 6.1484375]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.58 | 6.30 | 10.96× | 0.53 | 3.62 | 6.89× | 148.95 | 383.65 | 109.51 | 51.44 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.30 | 5.90 | 19.62× | 0.29 | 3.29 | 11.54× | 144.99 | 368.38 | 109.47 | 51.46 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 3.18 | 24.38 | n/a | 3.17 | 20.67 | n/a | 142.18 | 529.03 | 89.99 | 76.45 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.30 | 5.48 | 18.43× | 0.28 | 2.80 | 10.20× | 32.56 | 365.57 | 1.24 | 50.93 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 39.62 | 21.88 | 0.55× | 39.58 | 19.39 | 0.49× | 717.56 | 653.11 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.02 | 5.06 | n/a | 0.00 | 2.36 | n/a | 1.06 | 327.83 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 5.37 | 445.20× | 0.00 | 2.78 | n/a | 1.79 | 686.25 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.87 | 11.38 | n/a | 2.85 | 8.52 | n/a | 101.18 | 1155.61 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 6.09 | 514.52× | 0.00 | 2.85 | n/a | 1.06 | 393.43 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.58 | 4.24 | 7.25× | 0.82 | 5.23 | 6.34× | 144.62 | 390.08 | 109.50 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.34 | 4.39 | 13.01× | 0.61 | 5.12 | 8.40× | 147.09 | 320.38 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.77 | 17.08 | n/a | 3.71 | 22.50 | n/a | 141.98 | 749.10 | 89.99 | 76.46 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.25 | 4.56 | 18.31× | 1.07 | 3.92 | 3.67× | 30.29 | 324.24 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 37.52 | 16.90 | 0.45× | 62.45 | 32.97 | 0.53× | 717.73 | 948.98 | 37.41 | 123.18 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.02 | 4.01 | n/a | 0.00 | 3.51 | n/a | 0.93 | 268.89 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 3.98 | 319.84× | 0.00 | 3.99 | n/a | 0.99 | 535.61 | 0.00 | 48.18 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.84 | 6.62 | n/a | 2.81 | 14.86 | n/a | 101.19 | 1486.97 | 51.68 | 99.89 | 0.00 | 0.03 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.22 | 4.83 | 22.46× | 5.01 | 3.95 | 0.79× | 10.35 | 337.85 | 0.02 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
