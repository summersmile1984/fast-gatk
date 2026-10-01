# VariantAnnotator rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 22.795s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_variant_annotator_coverage_gatk_oracle.py` | 0 | 9.656s | 7801.22s | 6.5 MiB | 273.0 MiB | 9.62s | ✅ |
| `verify_variant_annotator_resource_expression_gatk_oracle.py` | 0 | 13.139s | 7807.54s | 2.2 MiB | 318.0 MiB | 13.10s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variant-annotator
```

## JSON sidecar

See `variant-annotator-rerun-20260928.json` for full stdout/stderr tails.
