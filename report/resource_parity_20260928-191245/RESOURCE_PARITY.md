# Resource parity: native fastgatk vs GATK 4.6.2.0

Report date: 2026-09-28  |  protocol: equal-resource-v1  |  repetitions: 2

Equal-resource protocol (v1): identical fixture and functional CLI per
case; both sides of a lane get the SAME envelope (taskset affinity +
OMP_NUM_THREADS, applied outside the CLI); runs serialized; warm cache;
p50 reported; output-equivalence asserted per case before any ratio is
reported. Lane `1-core` = one core each (algorithmic efficiency); lane
`full-machine` = the whole machine for both sides (as-shipped
throughput).

Machine: AMD Ryzen 9 7945HX with Radeon Graphics | 32 cpus | loadavg [1.12744140625, 1.67919921875, 2.29833984375] -> [2.29541015625, 2.220703125, 2.25927734375]

git HEAD: `2e0da9215b730fce6eea2b0d4fe4255903951ee2`  |  GATK jar sha256: `40b494b1e356931c…`  |  java: openjdk version "17.0.20" 2026-07-21

## Lane `1-core`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.45 | 5.38 | 11.86× | 0.43 | 2.77 | 6.38× | 158.47 | 390.82 | 109.47 | 51.46 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.28 | 5.24 | 19.03× | 0.25 | 2.61 | 10.44× | 153.56 | 366.69 | 109.47 | 51.46 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 19.17 | 20.31 | 1.06× | 19.09 | 16.95 | 0.89× | 399.48 | 523.28 | 95.26 | 76.45 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.27 | 4.90 | 18.26× | 0.25 | 2.37 | 9.46× | 24.80 | 394.68 | 1.24 | 50.90 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 13.10 | 19.80 | 1.51× | 13.07 | 16.84 | 1.29× | 717.74 | 653.62 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.01 | 4.74 | 323.67× | 0.00 | 2.04 | n/a | 0.87 | 326.14 | 0.00 | 48.05 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 4.89 | 472.42× | 0.00 | 2.31 | n/a | 0.55 | 685.12 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.75 | 9.48 | 3.45× | 2.72 | 6.59 | 2.42× | 97.18 | 1155.36 | 51.68 | 99.87 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.01 | 8.22 | 893.37× | 0.00 | 2.45 | n/a | 0.93 | 390.10 | 0.00 | 49.91 | 0.00 | 0.00 |

## Lane `full-machine`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.47 | 4.02 | 8.47× | 0.78 | 4.38 | 5.62× | 131.47 | 374.93 | 109.51 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.28 | 3.90 | 13.85× | 0.44 | 4.09 | 9.39× | 149.55 | 373.90 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 17.30 | 15.74 | 0.91× | 23.58 | 20.39 | 0.86× | 400.73 | 734.58 | 95.26 | 76.46 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.22 | 3.78 | 17.48× | 0.89 | 3.33 | 3.75× | 28.79 | 278.78 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 10.88 | 14.66 | 1.35× | 32.17 | 27.75 | 0.86× | 717.61 | 932.77 | 37.41 | 123.15 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.02 | 10.30 | 516.60× | 0.00 | 2.90 | n/a | 0.54 | 259.70 | 0.00 | 48.04 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.01 | 3.91 | 390.74× | 0.00 | 3.32 | n/a | 0.37 | 541.25 | 0.00 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.81 | 7.40 | 2.63× | 2.79 | 13.25 | 4.76× | 97.21 | 1495.68 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.13 | 4.74 | 35.48× | 3.48 | 3.69 | 1.06× | 10.55 | 353.26 | 0.02 | 49.90 | 0.00 | 0.00 |

## Lane `full-machine-tuned`

| Case | Output equal | Wall native (s) | Wall GATK (s) | Wall speedup | CPU native (s) | CPU GATK (s) | CPU ratio | TreeRSS native (MiB) | TreeRSS GATK (MiB) | rchar native (MiB) | rchar GATK (MiB) | disk-read native (MiB) | disk-read GATK (MiB) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HC | yes | 0.45 | 4.07 | 8.95× | 0.44 | 4.33 | 9.85× | 153.06 | 305.45 | 109.47 | 51.45 | 0.00 | 0.00 |
| HC (A-line regression window) | yes | 0.32 | 4.02 | 12.76× | 0.27 | 4.02 | 14.91× | 145.94 | 362.26 | 109.50 | 51.43 | 0.00 | 0.00 |
| BQSR (BaseRecalibrator) | yes | 16.76 | 15.55 | 0.93× | 22.57 | 20.39 | 0.90× | 400.76 | 742.97 | 95.26 | 76.44 | 0.00 | 0.00 |
| Mutect2 (tumor+normal) | yes | 0.21 | 3.77 | 17.57× | 0.36 | 3.17 | 8.82× | 24.96 | 312.21 | 1.23 | 50.92 | 0.00 | 0.00 |
| Mutect2 (DREAM somatic, hs37d5 ref) | yes | 10.95 | 14.73 | 1.34× | 17.56 | 27.68 | 1.58× | 717.82 | 954.72 | 37.41 | 123.16 | 0.00 | 0.00 |
| SortSam (coordinate) | yes | 0.01 | 4.46 | 303.75× | 0.00 | 2.91 | n/a | 1.40 | 256.35 | 0.00 | 48.04 | 0.00 | 0.00 |
| MarkDuplicates (small fixture) | yes | 0.02 | 3.79 | 164.99× | 0.00 | 3.18 | n/a | 7.70 | 530.83 | 0.13 | 48.20 | 0.00 | 0.00 |
| MarkDuplicates (CEUTrio chr20, 222k reads) | yes | 2.82 | 6.03 | 2.14× | 2.80 | 12.77 | 4.56× | 97.19 | 1407.83 | 51.68 | 99.89 | 0.00 | 0.00 |
| GenotypeGVCFs (HC + CombineGVCFs + GenotypeGVCFs) | yes | 0.02 | 4.55 | 246.15× | 0.11 | 3.29 | 31.33× | 1.29 | 304.19 | 0.00 | 49.90 | 0.00 | 0.00 |

Notes:
* Wall speedup = GATK p50 / native p50 (>1 = native faster); CPU ratio
  is the equal-resource efficiency measure (same work, same envelope).
* TreeRSS = peak of the 25 ms-sampled sum of per-process RSS over the
  whole workload process tree (native/java alike).
* rchar = bytes requested via read syscalls (cache-inclusive volume);
  disk-read = physical `read_bytes` (near zero on warm cache).
* CPU = GNU time user+sys (wait4 tree rollup). Full per-run data and
  FS-block counters in `summary.json`.
