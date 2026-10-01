# AnalyzeCovariates rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 35.147s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_analyze_covariates.py` | 0 | 21.623s | 23.79s | 6.4 MiB | 337.1 MiB | 21.57s | ✅ |
| `verify_analyze_covariates_bqsr_alias_gatk_oracle.py` | 0 | 13.524s | 12.85s | 14.6 MiB | 500.9 MiB | 13.31s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool analyze-covariates
```

## JSON sidecar

See `analyze-covariates-rerun-20260928.json` for full stdout/stderr tails.
