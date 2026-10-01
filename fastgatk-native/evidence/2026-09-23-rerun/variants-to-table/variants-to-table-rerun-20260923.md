# VariantsToTable rerun report — 2026-09-23

## Summary
- Total scripts: 3
- Passed: 3    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 125.516s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_malformed_input_fail_loud_oracle.py` | 0 | 2.228s | 7426.55s | 9.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_variants_to_table.py` | 0 | 91.429s | 7478.81s | 5.4 MiB | 310.4 MiB | 91.35s | ✅ |
| `verify_variants_to_table_gatk_oracle.py` | 0 | 31.859s | 7440.35s | 0.0 MiB | 313.1 MiB | 31.82s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool variants-to-table
```

## JSON sidecar

See `variants-to-table-rerun-20260923.json` for full stdout/stderr tails.
