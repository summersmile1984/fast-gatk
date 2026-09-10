# GenomicsDB native storage audit

The `--fastgatk-native-workspace` format is deliberately **not** a TileDB or
GenomicsDB workspace.  It is `fastgatk-portable-sparse-v1`: materialized local
VCF/GVCF shards, an ordered input index, and a coarse contig/span index used by
the native GenotypeGVCFs reader to skip non-overlapping shards.

## Dependency result

The pinned GATK 4.6.2.0 package contains `libtiledbgenomicsdb.so.1`.  Its stable
native surface available here is the legacy-ABI `GenomicsDB::generate_vcf`
query path used by `fastgatk-genomicsdb-export`.  Import is exposed through JNI
entry points and Java protobuf configuration; the repository has no matching
GenomicsDB importer or TileDB schema headers/SDK.  Reconstructing the array
schema, cell encoding, fragment metadata, and consolidation ABI from exported
symbols would be release-specific and unsafe.

Therefore native TileDB storage is a no-go with the present dependencies.  A
real workspace continues to use the resource-aware external GATK adapter.  The
isolated bridge continues to query those real workspaces.  Neither fallback is
changed by the portable sparse implementation.

## Executable gates

- `fastgatk-genomicsdb-import-gatk-oracle` creates a real workspace through the
  adapter and reopens it with GATK 4.6.2.0.
- `fastgatk-genomicsdb-bridge-gatk-oracle` queries a real workspace through the
  isolated legacy-ABI bridge.
- `fastgatk-genomicsdb-native-storage-boundary` proves that portable sparse
  metadata explicitly reports no TileDB/query compatibility and that both the
  bridge and GATK 4.6.2.0 reject it.
