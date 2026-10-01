# GenomicsDBImport rerun report — 2026-09-23

## Summary
- Total scripts: 7
- Passed: 7    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 128.648s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_genomicsdb_bridge.py` | 0 | 14.597s | 1405.25s | 317.2 MiB | 628.8 MiB | 14.09s | ✅ |
| `verify_genomicsdb_import.py` | 0 | 1.231s | 1368.93s | 8.5 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_genomicsdb_import_gatk_oracle.py` | 0 | 12.893s | 1384.08s | 4.8 MiB | 616.7 MiB | 12.84s | ✅ |
| `verify_genomicsdb_import_native_interval_gatk_oracle.py` | 0 | 14.418s | 1394.82s | 314.4 MiB | 398.0 MiB | 14.07s | ✅ |
| `verify_genomicsdb_import_sample_map_gatk_oracle.py` | 0 | 18.083s | 1415.61s | 8.9 MiB | 635.2 MiB | 18.05s | ✅ |
| `verify_genomicsdb_import_update_workspace_gatk_oracle.py` | 0 | 61.966s | 1441.09s | 9.6 MiB | 639.8 MiB | 61.88s | ✅ |
| `verify_genomicsdb_native_storage_boundary.py` | 0 | 5.460s | 1373.95s | 6.6 MiB | 275.6 MiB | 5.42s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool genomicsdb-import
```

## JSON sidecar

See `genomicsdb-import-rerun-20260923.json` for full stdout/stderr tails.
