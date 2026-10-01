# LeftAlignAndTrimVariants rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 164.203s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_left_align.py` | 0 | 43.231s | 5575.26s | 7.4 MiB | 355.2 MiB | 43.03s | ✅ |
| `verify_left_align_cli_boundary_gatk_oracle.py` | 0 | 53.825s | 5601.28s | 7.2 MiB | 325.1 MiB | 53.76s | ✅ |
| `verify_left_align_gatk_oracle.py` | 0 | 38.596s | 5549.88s | 7.3 MiB | 324.9 MiB | 38.48s | ✅ |
| `verify_left_align_sites_only_gatk_oracle.py` | 0 | 28.551s | 5525.64s | 7.2 MiB | 308.5 MiB | 28.47s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool left-align-trim
```

## JSON sidecar

See `left-align-trim-rerun-20260928.json` for full stdout/stderr tails.
