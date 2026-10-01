# AnalyzeCovariates rerun report — 2026-09-25

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 35.037s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_analyze_covariates.py` | 0 | 19.731s | 32.64s | 6.3 MiB | 328.4 MiB | 19.24s | ✅ |
| `verify_analyze_covariates_bqsr_alias_gatk_oracle.py` | 0 | 15.306s | 17.31s | 14.8 MiB | 532.6 MiB | 12.91s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool analyze-covariates
```

## JSON sidecar

See `analyze-covariates-rerun-20260925.json` for full stdout/stderr tails.
