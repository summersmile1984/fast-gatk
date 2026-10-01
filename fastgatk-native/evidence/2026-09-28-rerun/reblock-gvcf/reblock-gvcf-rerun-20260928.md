# ReblockGVCF rerun report — 2026-09-28

## Summary
- Total scripts: 6
- Passed: 5    Failed: 1    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 124.169s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_reblock_gatk_oracle.py` | 0 | 63.437s | 99.15s | 7.6 MiB | 447.8 MiB | 63.34s | ✅ |
| `verify_reblock_gatk_shards.py` | 0 | 20.435s | 55.80s | 5.2 MiB | 322.9 MiB | 15.90s | ✅ |
| `verify_reblock_gvcf.py` | 0 | 10.697s | 24.99s | 7.2 MiB | 381.8 MiB | 10.48s | ✅ |
| `verify_reblock_gvcf_droplowqual_gatk_oracle.py` | 0 | 14.858s | 37.91s | 5.7 MiB | 286.6 MiB | 14.79s | ✅ |
| `verify_reblock_gvcf_triploid_gatk_oracle.py` | 1 | 9.593s | 16.08s | 7.0 MiB | 315.8 MiB | 9.52s | ❌ |
| `verify_reblock_overlap_gatk_oracle.py` | 0 | 5.149s | 6.05s | 0.0 MiB | 410.0 MiB | 5.11s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool reblock-gvcf
```

## JSON sidecar

See `reblock-gvcf-rerun-20260928.json` for full stdout/stderr tails.
