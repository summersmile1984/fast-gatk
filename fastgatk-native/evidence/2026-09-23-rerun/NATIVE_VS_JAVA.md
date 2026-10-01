# Native and Java resource report — historical data withdrawn

## Evidence status (2026-09-24)

The previous resource table is **unvalidated and invalid as performance
 evidence**. It contained sums of per-script RSS peaks mislabeled as peak
memory and ratios of unlike native/Java invocation totals. The historical
records also include manually estimated Mutect2 values, cumulative child CPU
accounting and incomplete sampler finalization. Those defects cannot be
repaired by recomputing a ratio or taking a maximum of the same source values.

The prior “native wins” and anomaly conclusions are withdrawn. No native
optimization, speedup, memory reduction or memory leak has been established
by this report. The old resource table is removed rather than presenting
unsupported numbers as verified measurements. JSON sidecars are retained
unchanged for traceability, not endorsed as valid resource evidence.

## Retained context

The historical report covered 48 tool groups. `fastgatk-genomicsdb-export`
is a native-only bridge without a verify script; that is a tool-registry fact,
not an inference from sampler counts.

The rerun driver grouped descendant processes by `comm`: `java`,
`fastgatk-*`, or other helpers. Recorded pid counts and lifetimes can miss
short-lived processes and processes still tracked at sampler shutdown.
A zero class count does not establish that a script ran only native code.

Tool groups can contain different numbers of native and Java invocations,
different parameters, and cross-tool pipelines. Their aggregate process
lifetimes are not paired benchmark elapsed times. Script exit-status counts
are distinct from output parity: comparisons must be assessed from the
specific oracle assertions and retained results, with skips reported separately.

## Corrected aggregation semantics

`compare_native_vs_java.py --repo <repository>` reads sidecars under that
repository and emits descriptive sampled values only:

* RSS is the maximum of available per-script class peaks, not their sum.
  It is an observed individual-process maximum, not simultaneous tree memory.
* Missing/null fields, zero RSS and unsampled records are unavailable, not
  imputed zero measurements. Partial aggregates show contributing script counts.
* Recorded lifetime sums and counts are descriptive and subject to the source
  sampler's limitations. No speedup ratios or anomaly rankings are produced.
* Re-aggregation does not authenticate manually patched or incomplete records.

## Fresh evidence pending

Fresh isolated Mutect2 measurements and independent output-parity checks are
pending in [work/mutect2-priority/](../../../work/mutect2-priority/).
No completed comparison or optimization result is claimed here. See the
[regression evidence workflow](../../docs/regression-evidence.md) for the
historical-data caveat and comparison requirements.
