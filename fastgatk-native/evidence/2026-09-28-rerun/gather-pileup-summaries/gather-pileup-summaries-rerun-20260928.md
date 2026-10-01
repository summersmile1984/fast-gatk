# GatherPileupSummaries rerun report — 2026-09-28

## Summary
- Total scripts: 2
- Passed: 2    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 11.147s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_gather_pileup_gatk_oracle.py` | 0 | 11.107s | 2129.99s | 4.4 MiB | 335.3 MiB | 11.07s | ✅ |
| `verify_gather_pileup_summaries.py` | 0 | 0.040s | 2121.11s | 2.1 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool gather-pileup-summaries
```

## JSON sidecar

See `gather-pileup-summaries-rerun-20260928.json` for full stdout/stderr tails.
