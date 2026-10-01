# CreateReadCountPanelOfNormals rerun report — 2026-09-23

## Summary
- Total scripts: 4
- Passed: 4    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 63.448s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_create_read_count_panel_of_normals.py` | 0 | 28.584s | 855.44s | 16.1 MiB | 560.3 MiB | 28.47s | ✅ |
| `verify_create_read_count_panel_of_normals_degenerate_gatk_oracle.py` | 0 | 18.302s | 836.17s | 12.7 MiB | 547.9 MiB | 18.24s | ✅ |
| `verify_create_read_count_panel_of_normals_sample_metadata_gatk_oracle.py` | 0 | 16.404s | 823.92s | 12.7 MiB | 391.3 MiB | 16.34s | ✅ |
| `verify_pon_gatk_bundled_byte_structure_oracle.py` | 0 | 0.158s | 818.08s | 0.0 MiB | 0.0 MiB | 0.00s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool create-read-count-panel-of-normals
```

## JSON sidecar

See `create-read-count-panel-of-normals-rerun-20260923.json` for full stdout/stderr tails.
