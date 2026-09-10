# CNV HDF5 SimpleCountCollection metadata slice

This slice closes the file boundary shared by `CollectReadCounts`,
`DenoiseReadCounts`, and the count-input side of
`CreateReadCountPanelOfNormals`:

- `/sample_metadata/sample_name` is a one-element variable-length string;
- `/locatable_metadata/sequence_dictionary` is the HTSJDK
  `SAMTextHeaderCodec` dictionary (`@HD\tVN:1.6` plus every normalized `@SQ`
  record and its optional `AS`, `M5`, `UR`, and `SP` tags);
- `/intervals/indexed_contig_names` and
  `/intervals/transposed_index_start_end` retain ordered interval metadata;
- `/counts/values` is a `1 x N` finite, non-negative count matrix.

The native writer now preserves the complete `@SQ` metadata from the HTSlib
header, removes BAM-only `SO`/`RG` records from the dictionary field, and
validates dictionary membership, interval uniqueness, matrix dimensions, and
finite values.  The reader validates the same contract before exposing values
to the native denoising/PoN paths.

`verify_hdf5_simple_count_collection.py` is the focused oracle.  It compares a
native and GATK 4.6.2.0 `CollectReadCounts --format HDF5` dictionary, then
decodes both files through Java `DenoiseReadCounts` and requires identical
standardized/denoised tables.  The gate is registered as
`fastgatk-hdf5-simple-count-collection-metadata` and is included in
`verify_all.sh`; the existing contract remains separate for read-start/filter
semantics.

This does not claim complete CNV equivalence.  TSV remains an interchange
approximation, while HDF5 SVD panel/model state, GC correction parity, cloud
I/O, and release-specific HDF5 byte identity remain outside this slice and are
covered (where available) by the dedicated Denoise/PoN contracts.
