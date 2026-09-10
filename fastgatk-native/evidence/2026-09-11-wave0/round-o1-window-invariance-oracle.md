# O1 — Acceptance oracle for the D2 `-L`-window-dependent annotation (write-only round)

Status: **delivered**. The oracle is in place, runs the pinned GATK 4.6.2.0 jar and the native
binary per window, and **fails (exit 1) on the current binaries exactly at the documented D2
divergence**; `--expect-divergence` reports the same numbers and exits 0.

Discipline: no production source edited, nothing rebuilt, no `git commit`, no repo-root `*.md`
touched, `CMakeLists.txt` untouched (registration snippet only, in §5). The only files created are

* `fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py` (the oracle)
* this report, `.diag/round-o1-window-invariance-oracle.md`

Every experiment ran inside a single `bash` invocation with `timeout`, with all scratch in
`tempfile.TemporaryDirectory`.

---

## 1. Deliverable and invocation

```bash
# strict assertion (default; this is the gate) — non-zero on mismatch
python3 fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py

# diagnostic report of the current, unfixed state — always exit 0
python3 fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py --expect-divergence

# optional extras used below
python3 fastgatk-native/scripts/verify_hc_window_invariance_gatk_oracle.py \
    --serial-native fastgatk-native/build-serial/fastgatk-hc-call \
    --gate-structural-windows --windows 10019901,10020381,10020421
```

Observed exit codes on the current binaries:

| mode | exit code | JSON `status` | violations |
| --- | ---: | --- | ---: |
| strict (default) | **1** | `fail` | 8 |
| `--expect-divergence` | **0** | `diagnostic` | 8 (reported, not fatal) |
| `--gate-structural-windows --serial-native ...` (3 windows) | **1** | `fail` | 6 |

Flags: `--native` (default `$FASTGATK_HC_BINARY`, else `fastgatk-native/build/fastgatk-hc-call`),
`--serial-native` (optional second backend), `--expect-divergence`,
`--gate-structural-windows`, `--windows` (start list / range override), `--threads` (default 1).
House conventions are followed: `main() -> int`, `records()` data-row reader, `tempfile.TemporaryDirectory`,
`FASTGATK_REQUIRE_GATK_ORACLE` to turn an asset-absence skip into a hard error, and a final JSON status line.

## 2. What the oracle pins

* **GATK**: `third_party/jdk17/bin/java -Xmx1g -jar
  third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar HaplotypeCaller`, with
  `--native-pair-hmm-threads 1 --create-output-variant-index false
  --add-output-vcf-command-line false` (the same baseline flags as the sibling oracles).
* **Native**: the OpenMP `fastgatk-hc-call` (Serial opt-in via `--serial-native`).
* **Fixture**: `fixtures/chr20/ref20mnp.fasta`, `fixtures/chr20/mnp.bam`.
* **Per-window command** (identical on both sides apart from the backend's own threading flag):
  `-R <ref> -I <bam> -L 20:<start>-10020710 --emit-ref-confidence GVCF --max-mnp-distance 1
  --threads 1 --add-output-vcf-command-line false`. No streaming flag anywhere.
* **Comparison unit**: data rows only (`#`-comment/provenance lines are skipped by construction), full
  tab-split field equality; `--add-output-vcf-command-line false` removes the one provenance
  difference the harness itself could introduce.
* **Assertions** (all strict-mode):
  1. for every **gating** window, native data rows are byte-identical to GATK data rows;
  2. for **every** window, the `POS 10020680` row (ALT `AT,*,<NON_REF>`) is byte-identical, i.e.
     `RAW_MQandDP` and `SB` match GATK;
  3. the live GATK values at `POS 10020680` equal the pinned per-window expectation table
     (`GATK_D2_EXPECTATION`), so an accidental GATK-invocation drift cannot make a broken native build
     "agree".
  A failing window prints both sides' `RAW_MQandDP`/`SB`, the ROW counts, the divergent POS list and the
  field-level diff of the first divergent row (`INFO RAW_MQandDP '28800,8' -> '97200,27'; FORMAT SB '0,0,3,3' -> '0,0,0,0'`).

**Window set.** 10 gating windows — `10019901` (window A, the 48-row reference baseline), `10020201`,
`10020301`, `10020351`, `10020381` (window B, the documented divergence), `10020391`, `10020401`,
`10020411`, `10020441`, `10020501` — plus 2 reported-but-non-gating windows `10020421`, `10020431`
(see §3 for why). The pinned `POS 10020680` expectation is `RAW_MQandDP=28800,8` / `SB=0,0,3,3` for
every gating window; the two structural windows are pinned at `97200,27` / `1,2,6,16`.

## 3. The extra windows `10020421` / `10020431`: the prior audit's "third wrong state" is **refuted**

I verified the audit's claim directly against the pinned GATK before including those windows. The
audit (`track-d-streaming-rootcause.md:61,327`) states `-L 20:10020421-…` / `-L 20:10020431-…`
produce "a third, distinct wrong state (`SB=1,2,6,16`)". That is **not a native divergence**: the
pinned GATK itself emits that exact state in exactly those two windows, and native matches it
field-for-field.

| `-L 20:<start>-10020710` | GATK rows | native rows | `POS 10020680` GATK | `POS 10020680` native | rows differ |
| ---: | ---: | ---: | --- | --- | --- |
| 10019901 | 48 | 48 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020201 | 44 | 44 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020301 | 13 | 13 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020351 | 13 | 13 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| **10020381** | 13 | 13 | `28800,8` / `0,0,3,3` | **`97200,27` / `0,0,0,0`** | 1 (`10020680`) |
| **10020391** | 13 | 13 | `28800,8` / `0,0,3,3` | **`97200,27` / `0,0,0,0`** | 1 (`10020680`) |
| **10020401** | 13 | 13 | `28800,8` / `0,0,3,3` | **`97200,27` / `0,0,0,0`** | 1 (`10020680`) |
| **10020411** | 13 | 13 | `28800,8` / `0,0,3,3` | **`97200,27` / `0,0,0,0`** | 1 (`10020680`) |
| 10020421 | 53 | 12 | **`97200,27` / `1,2,6,16`** | `97200,27` / `1,2,6,16` ✔ | 50 (block granularity) |
| 10020431 | 51 | 9 | **`97200,27` / `1,2,6,16`** | `97200,27` / `1,2,6,16` ✔ | 50 (block granularity) |
| 10020441 | 4 | 4 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020451 | 4 | 4 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020501 | 4 | 4 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020601 | 4 | 4 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |
| 10020661 | 4 | 4 | `28800,8` / `0,0,3,3` | `28800,8` / `0,0,3,3` | 0 |

GATK's own `10020421` output at `POS 10020680` (exactly reproduced by native):

```text
20  10020680  .  CA  AT,*,<NON_REF>  ... DP=27;ExcessHet=0.0000;MLEAC=1,1,0;MLEAF=0.500,0.500,0.00;RAW_MQandDP=97200,27 ...
   1|2:...:434,507,937:1,2,6,16
```

The **real** divergence in those two windows is not the D2 annotation row but reference-block
granularity: native collapses GATK's fine-grained `<NON_REF>` blocks (e.g. GATK
`10020432 … END=10020432 0/0:12:36:12:0,36,491`, `10020433 … END=10020433 0/0:13:39:13:0,39,532`,
…, `10020675 … END=10020678 /0:26:63:25:0,63,945`) into one degenerate block
`10020439 … END=10020678 0/0:26:0:20:0,0,0` (`GQ=0 / PL=0,0,0`). That is the D3 family (tile-local /
interval-window reference-confidence degeneracy, localized in `.diag/round-r2-d3-localization.md` as
"a pure function of the interval window handed to `calling::run`"), a different defect from the D2
annotation-evidence defect this oracle gates.

Consequences for the oracle design (deliberate, documented):

* both windows are **run and fully reported** in every mode, and their `POS 10020680` row **is gated**
  (it must stay byte-identical to GATK — it currently is);
* their *full row set* is **not** in the default gate, because gating it would tie the D2 fix to the
  unrelated D3 defect and the oracle could not go green when D2 is fixed;
* `--gate-structural-windows` promotes them into the full-row gate, so they can be folded in as soon
  as D3 is also addressed.

## 4. Observed output on the current binaries

### 4.1 strict mode (default) — `EXIT=1`, JSON `status=fail`, 8 violations

Report header and per-window block (verbatim stdout; the per-window numbers are identical in both
modes — only the trailing JSON status and the exit code differ, see §4.2):

```text
# verify_hc_window_invariance_gatk_oracle.py: pinned GATK 4.6.2.0 vs native, fixture fixtures/chr20/mnp.bam, windows '-L 20:<start>-10020710'
# mode: strict assertion; gate-structural-windows=False; threads=1
# backends: OpenMP
[window 20:10019901-10020710] gatk_rows=48 native_rows=48 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  pinned GATK expectation at POS 10020680: {'RAW_MQandDP': '28800,8', 'SB': '0,0,3,3'}
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=48
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=48
[window 20:10020201-10020710] gatk_rows=44 native_rows=44 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=44
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=44
[window 20:10020301-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
[window 20:10020351-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
[window 20:10020381-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=False d2_row_identical=False gating=True
  [OpenMP] divergent positions: ['10020680']
  pinned GATK expectation at POS 10020680: {'RAW_MQandDP': '28800,8', 'SB': '0,0,3,3'}
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=97200,27 SB=0,0,0,0  rows=13
    first divergent row field diff (GATK -> NATIVE): INFO RAW_MQandDP '28800,8' -> '97200,27'; FORMAT SB '0,0,3,3' -> '0,0,0,0'
[window 20:10020391-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=False d2_row_identical=False gating=True
  [OpenMP] divergent positions: ['10020680']
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=97200,27 SB=0,0,0,0  rows=13
    first divergent row field diff (GATK -> NATIVE): INFO RAW_MQandDP '28800,8' -> '97200,27'; FORMAT SB '0,0,3,3' -> '0,0,0,0'
[window 20:10020401-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=False d2_row_identical=False gating=True
  [OpenMP] divergent positions: ['10020680']
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=97200,27 SB=0,0,0,0  rows=13
    first divergent row field diff (GATK -> NATIVE): INFO RAW_MQandDP '28800,8' -> '97200,27'; FORMAT SB '0,0,3,3' -> '0,0,0,0'
[window 20:10020411-10020710] gatk_rows=13 native_rows=13 data_rows_byte_identical=False d2_row_identical=False gating=True
  [OpenMP] divergent positions: ['10020680']
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=13
    NATIVE RAW_MQandDP=97200,27 SB=0,0,0,0  rows=13
    first divergent row field diff (GATK -> NATIVE): INFO RAW_MQandDP '28800,8' -> '97200,27'; FORMAT SB '0,0,3,3' -> '0,0,0,0'
[window 20:10020441-10020710] gatk_rows=4 native_rows=4 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=4
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=4
[window 20:10020501-10020710] gatk_rows=4 native_rows=4 data_rows_byte_identical=True d2_row_identical=True gating=True
  [OpenMP] divergent positions: []
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=28800,8 SB=0,0,3,3  rows=4
    NATIVE RAW_MQandDP=28800,8 SB=0,0,3,3  rows=4
[window 20:10020421-10020710] gatk_rows=53 native_rows=12 data_rows_byte_identical=False d2_row_identical=True gating=False
  [OpenMP] divergent positions: ['10020430', '10020432', '10020433/10020434', '10020434/10020435', '10020435/10020438', '10020438/10020439', '10020439/10020679', '10020445/10020680', '10020447/10020682', 'missing-in-native@12:POS=10020449'] (+40 more)
  pinned GATK expectation at POS 10020680: {'RAW_MQandDP': '97200,27', 'SB': '1,2,6,16'}
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=97200,27 SB=1,2,6,16  rows=53
    NATIVE RAW_MQandDP=97200,27 SB=1,2,6,16  rows=12
    first divergent row field diff (GATK -> NATIVE): FORMAT DP '12' -> '31'; FORMAT GQ '36' -> '0'; FORMAT MIN_DP '12' -> '31'; FORMAT PL '0,36,498' -> '0,0,0'
[window 20:10020431-10020710] gatk_rows=51 native_rows=9 data_rows_byte_identical=False d2_row_identical=True gating=False
  [OpenMP] divergent positions: ['10020432', '10020433/10020434', '10020434/10020435', '10020435/10020438', '10020438/10020439', '10020439/10020679', '10020445/10020680', '10020447/10020682', 'missing-in-native@9:POS=10020449', 'missing-in-native@10:POS=10020451'] (+40 more)
  pinned GATK expectation at POS 10020680: {'RAW_MQandDP': '97200,27', 'SB': '1,2,6,16'}
  [OpenMP] POS 10020680 ALT=AT,*,<NON_REF>
    GATK   RAW_MQandDP=97200,27 SB=1,2,6,16  rows=51
    NATIVE RAW_MQandDP=97200,27 SB=1,2,6,16  rows=9
    first divergent row field diff (GATK -> NATIVE): INFO END '10020432' -> '10020433'; FORMAT DP '11' -> '31'; FORMAT GQ '33' -> '0'; FORMAT MIN_DP '11' -> '31'; FORMAT PL '0,33,449' -> '0,0,0'
# 8 violation(s):
#   - 20:10020381-10020710 [OpenMP]: POS 10020680 annotation differs (GATK RAW_MQandDP=28800,8 SB=0,0,3,3 vs native RAW_MQandDP=97200,27 SB=0,0,0,0)
#   - 20:10020381-10020710 [OpenMP]: data rows differ at ['10020680'] (gatk_rows=13 native_rows=13)
#   - 20:10020391-10020710 [OpenMP]: POS 10020680 annotation differs (GATK RAW_MQandDP=28800,8 SB=0,0,3,3 vs native RAW_MQandDP=97200,27 SB=0,0,0,0)
#   - 20:10020391-10020710 [OpenMP]: data rows differ at ['10020680'] (gatk_rows=13 native_rows=13)
#   - 20:10020401-10020710 [OpenMP]: POS 10020680 annotation differs (GATK RAW_MQandDP=28800,8 SB=0,0,3,3 vs native RAW_MQandDP=97200,27 SB=0,0,0,0)
#   - 20:10020401-10020710 [OpenMP]: data rows differ at ['10020680'] (gatk_rows=13 native_rows=13)
#   - 20:10020411-10020710 [OpenMP]: POS 10020680 annotation differs (GATK RAW_MQandDP=28800,8 SB=0,0,3,3 vs native RAW_MQandDP=97200,27 SB=0,0,0,0)
#   - 20:10020411-10020710 [OpenMP]: data rows differ at ['10020680'] (gatk_rows=13 native_rows=13)
```

JSON status line (the `results` array was additionally emitted in full; the whole 819-line block
**parses as JSON** — verified):

```json
{
  "acceptance_criterion": "POS 10020680 must be byte-identical to pinned GATK in every window; data rows must be byte-identical for every gating window",
  "status": "fail",
  "strict_mode": true,
  "violations": [ ... 8 entries as listed above ... ],
  "gating_windows": ["20:10019901-10020710", "20:10020201-10020710", "20:10020301-10020710",
    "20:10020351-10020710", "20:10020381-10020710", "20:10020391-10020710",
    "20:10020401-10020710", "20:10020411-10020710", "20:10020441-10020710",
    "20:10020501-10020710"],
  "windows": [ ...the 10 above + "20:10020421-10020710", "20:10020431-10020710" ]
}
```

The `POS 10020680` requirement of the task is therefore enforced in the two documented windows:

| window | GATK | native | verdict |
| --- | --- | --- | --- |
| `20:10019901-10020710` (A) | `RAW_MQandDP=28800,8`, `SB=0,0,3,3` (48 rows) | identical (48 rows) | pass |
| `20:10020381-10020710` (B) | `RAW_MQandDP=28800,8`, `SB=0,0,3,3` (13 rows) | `RAW_MQandDP=97200,27`, `SB=0,0,0,0` (13 rows) | **fail (D2)** |

### 4.2 `--expect-divergence` — `EXIT=0`, JSON `status=diagnostic`

```text
# mode: expect-divergence (diagnostic); gate-structural-windows=False; threads=1
# backends: OpenMP
```

Per-window output is **byte-identical to §4.1** (same windows, same numbers; verified by parsing both
runs). Only two things change: the header line above, and the trailing JSON status/exit code:

```json
{
  "acceptance_criterion": "POS 10020680 must be byte-identical to pinned GATK in every window; data rows must be byte-identical for every gating window",
  "status": "diagnostic",
  "strict_mode": false,
  "violations": [ ...same 8 entries, reported but non-fatal... ]
}
```

Condensed per-window observation (identical in both modes), from the JSON `per_window` block:

| window | d2 row identical | data rows identical | GATK (RAW_MQandDP \| SB) | native (RAW_MQandDP \| SB) |
| --- | --- | --- | --- | --- |
| 20:10019901-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |
| 20:10020201-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |
| 20:10020301-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |
| 20:10020351-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |
| 20:10020381-10020710 | **False** | **False** | `28800,8` \| `0,0,3,3` | **`97200,27` \| `0,0,0,0`** |
| 20:10020391-10020710 | **False** | **False** | `28800,8` \| `0,0,3,3` | **`97200,27` \| `0,0,0,0`** |
| 20:10020401-10020710 | **False** | **False** | `28800,8` \| `0,0,3,3` | **`97200,27` \| `0,0,0,0`** |
| 20:10020411-10020710 | **False** | **False** | `28800,8` \| `0,0,3,3` | **`97200,27` \| `0,0,0,0`** |
| 20:10020421-10020710 | True | False (D3) | `97200,27` \| `1,2,6,16` | `97200,27` \| `1,2,6,16` |
| 20:10020431-10020710 | True | False (D3) | `97200,27` \| `1,2,6,16` | `97200,27` \| `1,2,6,16` |
| 20:10020441-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |
| 20:10020501-10020710 | True | True | `28800,8` \| `0,0,3,3` | `28800,8` \| `0,0,3,3` |

### 4.3 Extra checks: Serial backend and the promoted structural windows

`--gate-structural-windows --serial-native fastgatk-native/build-serial/fastgatk-hc-call
--windows 10019901,10020381,10020421` → `EXIT=1`, JSON `status=fail`, 6 violations, and:

* Serial and OpenMP produce identical results in every window (the defect is Host-side, as the
  evidence doc already recorded);
* window `10020421` promoted to gating produces full-row violations on **both** backends
  (`gatk_rows=53 native_rows=12`), i.e. the D3 block-granularity divergence, while its `POS 10020680`
  row stays identical to GATK;
* the check that the live GATK baseline equals `GATK_D2_EXPECTATION` produced **no** violations in any
  run, so the pinned expectation table is correct for all 12 windows (including `97200,27` / `1,2,6,16`
  for the two structural windows).

## 5. CMake registration snippet (proposed, **not applied**)

Belongs inside the existing `if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../fixtures/chr20/ref20mnp.fasta"
AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../fixtures/chr20/mnp.bam")` block
(`fastgatk-native/CMakeLists.txt:1452`), next to the other chr20 max-MNP oracles:

```cmake
        # D2 acceptance oracle: the POS 10020680 (ALT AT,*,<NON_REF>) annotation
        # evidence (RAW_MQandDP / SB) must not depend on the -L window start, and
        # every gating window must stay byte-identical to pinned GATK 4.6.2.0 with
        # no streaming at all.  Red until the interval-dependent annotation defect
        # is fixed; --expect-divergence turns it into a non-fatal diagnostic.
        add_test(NAME fastgatk-hc-window-invariance-gatk-oracle
            COMMAND ${Python3_EXECUTABLE}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_hc_window_invariance_gatk_oracle.py")
        set_tests_properties(fastgatk-hc-window-invariance-gatk-oracle PROPERTIES
            ENVIRONMENT "FASTGATK_HC_BINARY=$<TARGET_FILE:fastgatk-hc-call>;FASTGATK_REQUIRE_GATK_ORACLE=1")
```

Register it **only after the fix lands** — in strict mode it is red on the current binaries by design.
If a temporary non-fatal registration is ever wanted before that, append `--expect-divergence` to the
`COMMAND` (JSON `status` becomes `diagnostic`, exit 0) and drop `FASTGATK_REQUIRE_GATK_ORACLE=1`;
the strict form above is the one that should be the permanent gate.

## 6. Notes, limitations, uncertainty

* **Verified, not assumed**: every number in §3 and §4 was produced by running the pinned jar and the
  binaries during this round; nothing is quoted from earlier reports except as an explicit citation.
* **GATK is not globally window-invariant at `POS 10020680`.** It is invariant for every gating window
  (`28800,8` / `0,0,3,3`), but it *legitimately* changes to `97200,27` / `1,2,6,16` at window starts
  `10020421` / `10020431`, where native agrees. Any later statement that "GATK is always
  `28800,8`/`0,0,3,3` for this locus" would be wrong; the oracle encodes the per-window truth instead.
* **Prior-audit correction**: the claim that `10020421`/`10020431` show a third *native* wrong state
  is refuted (§3). The `SB=1,2,6,16` state is GATK's own value there.
* **Deliberate non-gating treatment of the two structural windows** is the one judgement call in this
  oracle. The alternative (gating them by default) would make the gate depend on the *unrelated* D3
  block-granularity defect and would leave the suite red after a correct D2 fix; the chosen default
  still gates their `POS 10020680` row, and `--gate-structural-windows` makes the stricter behaviour
  available in one flag. If the orchestrator prefers the strict-by-default variant, flipping
  `STRUCTURAL_WINDOWS` into `GATING_WINDOWS` (or defaulting `--gate-structural-windows` to true) is a
  two-line change in the oracle.
* **Runtime**: ~12 windows × (≈4.2 s GATK + native) ≈ 90 s per full run in this environment; the
  Serial extra roughly doubles the native time only.
* **No TIMEOUT property** is set in the snippet because this `CMakeLists.txt` does not use `TIMEOUT`
  properties anywhere.
* **Not verified here**: whether the planned D2 fix changes the structural windows at all, and whether
  the D3 defect (`10020421`/`10020431` block granularity) is in scope for the same fix. The oracle is
  built so that either outcome is observable and reported.
