# AnalyzeCovariates rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 30.294s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_analyze_covariates.py` | 0 | 18.681s | 24.29s | 6.4 MiB | 392.3 MiB | 18.62s | ✅ |
| `verify_analyze_covariates_bqsr_alias_gatk_oracle.py` | 0 | 11.613s | 13.46s | 15.2 MiB | 510.9 MiB | 11.40s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool analyze-covariates
```

## JSON sidecar

See `analyze-covariates-rerun-20260923.json` for full stdout/stderr tails.
