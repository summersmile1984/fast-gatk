# CreateReadCountPanelOfNormals rerun report — 2026-09-28

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 29.326s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_create_read_count_panel_of_normals.py` | 0 | 18.487s | 922.18s | 16.6 MiB | 613.4 MiB | 18.37s | ✅ |
| `verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py` | 0 | 5.912s | 902.65s | 11.9 MiB | 569.9 MiB | 5.85s | ✅ |
| `verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py` | 0 | 4.757s | 890.72s | 14.3 MiB | 422.8 MiB | 4.68s | ✅ |
| `verify_pon_gatk_bundled_byte_structure_oracle.py` | 0 | 0.170s | 884.42s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool create-read-count-panel-of-normals
```

## JSON sidecar

See `create-read-count-panel-of-normals-rerun-20260928.json` for full stdout/stderr tails.
