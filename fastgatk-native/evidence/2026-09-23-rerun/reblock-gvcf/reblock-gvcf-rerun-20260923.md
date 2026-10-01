# ReblockGVCF rerun report — 2026-09-23

## Summary
- Total scripts: 6
- Passed: 6    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 172.245s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_reblock_gatk_oracle.py` | 0 | 82.498s | 6568.82s | 7.6 MiB | 454.5 MiB | 82.41s | ✅ |
| `verify_reblock_gatk_shards.py` | 0 | 35.745s | 6523.24s | 6.9 MiB | 304.9 MiB | 35.68s | ✅ |
| `verify_reblock_gvcf.py` | 0 | 9.286s | 6481.08s | 7.0 MiB | 376.6 MiB | 8.94s | ✅ |
| `verify_reblock_gvcf_droplowqual_gatk_oracle.py` | 0 | 22.263s | 6504.27s | 5.5 MiB | 314.2 MiB | 22.18s | ✅ |
| `verify_reblock_gvcf_triploid_gatk_oracle.py` | 0 | 13.391s | 6491.07s | 7.1 MiB | 276.8 MiB | 13.33s | ✅ |
| `verify_reblock_overlap_gatk_oracle.py` | 0 | 9.062s | 6470.33s | 5.5 MiB | 417.7 MiB | 8.99s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool reblock-gvcf
```

## JSON sidecar

See `reblock-gvcf-rerun-20260923.json` for full stdout/stderr tails.
