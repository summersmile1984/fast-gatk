# CountReads rerun report — 2026-09-23

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 153.796s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_cli_parity_gatk_oracle.py` | 0 | 10.218s | 13.92s | 6.5 MiB | 365.6 MiB | 1.96s | ✅ |
| `verify_count_reads.py` | 0 | 143.578s | 103.06s | 6.7 MiB | 371.6 MiB | 138.40s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool count-reads
```

## JSON sidecar

See `count-reads-rerun-20260923.json` for full stdout/stderr tails.
