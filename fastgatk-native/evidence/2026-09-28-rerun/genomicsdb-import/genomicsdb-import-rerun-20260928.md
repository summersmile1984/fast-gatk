# GenomicsDBImport rerun report — 2026-09-28

## Summary
- Total scripts: 7
- Passed: 7    Failed: 0    Skipped (exit 77 / oracle-guard skip): 0
- Total elapsed: 84.031s

## Per-script results

| Script | Exit | Elapsed | CPU | Native peak RSS | Java peak RSS | Java wall-clock | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `verify_genomicsdb_bridge.py` | 0 | 10.358s | 2201.29s | 320.3 MiB | 651.3 MiB | 9.84s | ✅ |
| `verify_genomicsdb_import.py` | 0 | 1.248s | 2156.52s | 8.3 MiB | 0.0 MiB | 0.00s | ✅ |
| `verify_genomicsdb_import_gatk_oracle.py` | 0 | 9.926s | 2188.32s | 4.7 MiB | 605.4 MiB | 9.87s | ✅ |
| `verify_genomicsdb_import_native_interval_gatk_oracle.py` | 0 | 15.599s | 2215.22s | 319.4 MiB | 427.2 MiB | 15.25s | ✅ |
| `verify_genomicsdb_import_sample_map_gatk_oracle.py` | 0 | 9.820s | 2175.11s | 8.8 MiB | 656.7 MiB | 9.76s | ✅ |
| `verify_genomicsdb_import_update_workspace_gatk_oracle.py` | 0 | 32.748s | 2243.17s | 9.3 MiB | 613.3 MiB | 32.66s | ✅ |
| `verify_genomicsdb_native_storage_boundary.py` | 0 | 4.332s | 2161.97s | 5.3 MiB | 273.7 MiB | 4.30s | ✅ |

## Re-run command

```bash
python3 fastgatk-native/scripts/rerun_all_verify.py --tool genomicsdb-import
```

## JSON sidecar

See `genomicsdb-import-rerun-20260928.json` for full stdout/stderr tails.
