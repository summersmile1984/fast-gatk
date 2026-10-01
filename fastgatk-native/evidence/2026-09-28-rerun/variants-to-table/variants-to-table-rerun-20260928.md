# VariantsToTable rerun report — 2026-09-28

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 84.464s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_malformed_input_fail_loud_oracle.py` | 0 | 0.884s | 16.16s | 10.1 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_variants_to_table.py` | 0 | 64.293s | 66.21s | 5.2 MiB | 321.9 MiB | 64.21s | ✅ |
| `verify_variants_to_table_gatk_oracle.py` | 0 | 19.287s | 29.45s | 2.9 MiB | 313.4 MiB | 19.24s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variants-to-table
```

## JSON sidecar

See `variants-to-table-rerun-20260928.json` for full stdout/stderr tails.
