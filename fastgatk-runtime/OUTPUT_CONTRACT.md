# Runtime output contract

This specification defines the smallest reusable output boundary for native
tools. It is intentionally independent of Kokkos, HTSlib, JSON libraries and
tool-specific VCF/CRAM validation.

## Bundle

An `OutputBundle` contains:

- `primary`: a regular, non-empty primary output;
- `index`: an optional sidecar, or a required sidecar when `require_index=true`;
- `manifest`: a regular, non-empty JSON-object `OutputManifest`;
- `reject_temporary_residue`: whether conventional incomplete siblings are an
  error for the destination bundle.

The manifest is required to contain a numeric `schema_version`, a non-empty
`outputs` array, and `"complete":true` on every output entry. Missing or
`"complete":false` entries fail closed.
This envelope check does not replace tool-level validation of VCF headers,
coordinate order, CRC, tabix/BAI contents, or manifest-specific telemetry.

## Publication

Writers should create all artifacts under staging names on the destination
filesystem, then call:

```cpp
fastgatk::runtime::publish_output_bundle(staging, destination);
```

The call validates every staged artifact before the first rename, rejects an
already-existing destination, renames primary and index first, and renames the
manifest last. The manifest therefore acts as the completion marker consumed
by gather/resume paths. If a rename or post-publish validation fails, already
published files are moved back to their staging names and
`OUTPUT_CONTRACT_FAILURE` is thrown. A failed validation before publication
does not modify either staging or destination files.

The API is deliberately opt-in: adding it to the runtime library does not
change existing native tool behavior until a tool adopts the boundary.
