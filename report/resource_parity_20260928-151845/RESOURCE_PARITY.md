# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [6.88525390625, 13.19921875, 13.2255859375] -> [7.248046875, 9.57421875, 11.32666015625]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.57 | 6.79 | 11.93× | 0.53 | 4.00 | 7.55× | 158.80 | 384.75 | 109.47 | 51.46 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.33 | 6.03 | 18.01× | 0.31 | 3.40 | 10.98× | 151.29 | 365.35 | 109.48 | 51.46 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 3.33 | 24.67 | n/a | 3.32 | 20.78 | n/a | 140.98 | 526.92 | 89.99 | 76.45 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.32 | 5.88 | 18.55× | 0.29 | 3.08 | 10.44× | 32.59 | 391.42 | 1.24 | 50.93 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 14.92 | 23.40 | 1.57× | 14.88 | 20.67 | 1.39× | 717.63 | 651.75 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.02 | 5.51 | n/a | 0.01 | 2.55 | n/a | 0.68 | 325.39 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.02 | 5.82 | 373.98× | 0.00 | 2.98 | n/a | 0.93 | 682.82 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.92 | 11.51 | n/a | 2.91 | 8.56 | n/a | 101.25 | 1156.68 | 51.68 | 99.88 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 6.37 | 468.83× | 0.00 | 3.17 | n/a | 0.74 | 392.86 | 0.00 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.53 | 4.83 | 9.14× | 0.82 | 5.53 | 6.70× | 144.82 | 404.23 | 109.51 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.32 | 4.49 | 14.01× | 0.64 | 4.82 | 7.54× | 146.82 | 375.61 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.78 | 19.65 | n/a | 3.66 | 23.51 | n/a | 142.39 | 801.07 | 89.99 | 76.46 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.25 | 4.50 | 18.36× | 1.09 | 4.79 | 4.38× | 30.20 | 332.24 | 1.23 | 50.90 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 12.35 | 17.77 | 1.44× | 36.16 | 38.52 | 1.06× | 717.72 | 928.80 | 37.41 | 123.15 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.08 | 3.97 | n/a | 0.00 | 3.75 | n/a | 7.57 | 268.46 | 0.13 | 48.04 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.02 | 4.28 | 178.63× | 0.00 | 4.40 | n/a | 7.65 | 554.31 | 0.13 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.89 | 6.89 | n/a | 2.87 | 16.54 | n/a | 101.25 | 1470.65 | 51.68 | 99.88 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.23 | 4.85 | 20.79× | 5.17 | 4.43 | 0.86× | 10.54 | 329.00 | 0.02 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine-tuned`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.52 | 4.54 | 8.71× | 0.50 | 5.38 | 10.76× | 129.56 | 347.84 | 109.47 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.37 | 4.56 | 12.32× | 0.30 | 5.27 | 17.57× | 147.11 | 332.07 | 109.47 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | no | 2.80 | 17.75 | n/a | 3.28 | 24.23 | n/a | 142.61 | 805.25 | 89.99 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.29 | 9.65 | 33.40× | 0.46 | 4.54 | 9.75× | 30.18 | 271.71 | 1.23 | 50.93 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 12.61 | 16.91 | 1.34× | 19.95 | 34.34 | 1.72× | 717.94 | 962.36 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | no | 0.02 | 3.91 | n/a | 0.00 | 3.78 | n/a | 0.53 | 279.97 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.07 | 4.39 | 63.30× | 0.00 | 4.31 | n/a | 7.69 | 528.67 | 0.13 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | no | 2.89 | 6.60 | n/a | 2.87 | 15.35 | n/a | 101.25 | 1470.94 | 51.68 | 99.88 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.02 | 4.83 | 238.29× | 0.13 | 4.38 | 33.73× | 0.54 | 293.73 | 0.00 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
