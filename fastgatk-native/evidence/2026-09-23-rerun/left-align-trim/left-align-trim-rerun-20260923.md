# LeftAlignAndTrimVariants rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 199.623s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_left_align.py` | 0 | 55.200s | 4367.23s | 7.3 MiB | 377.4 MiB | 54.99s | ✅ |
| `verify_left_align_cli_boundary_gatk_oracle.py` | 0 | 67.792s | 4394.35s | 6.7 MiB | 308.1 MiB | 67.69s | ✅ |
| `verify_left_align_gatk_oracle.py` | 0 | 53.933s | 4340.73s | 7.4 MiB | 311.7 MiB | 53.76s | ✅ |
| `verify_left_align_sites_only_gatk_oracle.py` | 0 | 22.698s | 4314.60s | 5.6 MiB | 318.4 MiB | 22.57s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool left-align-trim
```

## JSON sidecar

See `left-align-trim-rerun-20260923.json` for full stdout/stderr tails.
