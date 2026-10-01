# CountReads rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 110.991s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_cli_parity_gatk_oracle.py` | 0 | 10.245s | 13.63s | 5.8 MiB | 358.2 MiB | 1.88s | ✅ |
| `verify_count_reads.py` | 0 | 100.746s | 97.09s | 6.6 MiB | 367.0 MiB | 96.89s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool count-reads
```

## JSON sidecar

See `count-reads-rerun-20260928.json` for full stdout/stderr tails.
