# Fast-GATK cluster operations runbook

This runbook is the M0 contract for operating fast-gatk on a real SLURM
cluster.  It is intentionally short: every step exists because the M0
infrastructure tests it, and every command is something that has either
already passed a CTest gate or is the gate itself.

The numbered sections mirror the failure modes a cluster operator will
hit first.  Read sections 1 and 2 before doing anything else; section 6
is the 5-minute rollback path.

---

## 1.  Mental model

A fast-gatk production worker carries **no JRE and no GATK jar**.  Every
binary in `fastgatk-native/build/` is dynamically linked against
`libgomp`/`libstdc++`/`libhts`/etc. on Ubuntu 22.04; that link set is the
runtime surface.  The pinned GATK 4.6.2.0 jar exists in exactly one place:
the `fastgatk-oracle` OCI image, and that image is only used by CI oracles.

Production workers run from the `fastgatk-native` OCI image.  That image
fails closed if any dispatch touches `--fallback`, an unknown tool, an
unknown option, an unsupported cloud/Spark option, or a missing SLURM
allocation.  Those exit-code-73 / `FALLBACK_DISABLED` events are the
intended production posture, not a bug.

Cluster ops must therefore reason about three independent sources of
evidence before trusting a result:

| Evidence | Lives in | Drives |
|----------|----------|--------|
| `OutputManifest` JSON per shard | shard output dir | gather/retry contract |
| `cluster-evidence.json` | `fastgatk-native/evidence/` | progress score gate |
| Nextflow `-log`/`.nextflow.log` | per-run work dir | scheduling diagnostics |

If any of the three is missing or contradicts the others, treat the run
as failed and consult section 6.

The base image for the production worker is `ubuntu:24.04` (noble).  The
build host is currently rolling Ubuntu (glibc 2.43, libstdc++ GLIBCXX
3.4.35).  When the build host moves further ahead, bump the base tag in
lockstep and update the docker-test-oracle CTest digest pin in
`build_container.sh`.

---

## 2.  Pre-flight checklist (do this before every cluster run)

1. `bash fastgatk-native/scripts/build_container.sh --tag $(date +%Y%m%d)` —
   rebuilds both images, writes `evidence/container-manifest.json`, and
   pins the digest.  A drift between the in-tree Dockerfile and the
   running image will fail-closed at section 5.
2. Confirm the native OCI image is what production will pull:

       docker inspect --format '{{.Id}} {{index .Config.Labels "fastgatk.image"}}' \
           fastgatk-native:$(date +%Y%m%d)

   The expected label is `fastgatk-native`.  Anything else means the build
   context leaked an unintended layer.
3. Run `python3 fastgatk-native/scripts/compute_progress_score.py` on the
   workstation.  It is read-only and will fail if `progress_score.json`
   has drifted from any of the tool audit hashes.  That is a release-gate
   stop condition; do not push to the cluster until it is clean.
4. Run `ctest -L m0_cluster_infra` on the workstation.  All four M0
   CTests must pass before cluster deployment: the container build, the
   cluster smoke (local-proof default), the infrastructure proof alias,
   and the cluster resource probe.  Any skip or failure indicates a
   wrapper or contract regression that the production run will inherit.
5. If `cluster-evidence.json` already exists in the tree, decide whether
   the new run will *replace* it (the file is a moving window; the most
   recent PASS wins).  Backing up the old file before re-running is the
   recommended pattern.

---

## 3.  Submitting a cluster run

Use the `production` Nextflow profile:

    nextflow run fastgatk-native/workflow/nextflow_cluster_smoke.nf \
        -profile production \
        --container fastgatk-native:20260916 \
        --cluster_account biology \
        --cluster_partition compute \
        --cluster_time 8h

The profile auto-selects the SLURM executor, the container image, the
retry policy (2 retries), and the exclusive-node flag.  A caller-provided
`--slurm_cluster_options` string still wins, so site wrappers can inject
`--gres`/`--constraint` flags without editing the config.

Two required environment variables that the profile **does not** set on
its own:

- `FASTGATK_SLURM_QUEUE` / `params.cluster_partition` — partition name.
- `FASTGATK_SLURM_ACCOUNT` / `params.cluster_account` — billing account.

Refusing to set these makes `verify_cluster_smoke.sh` exit 2
(`"FASTGATK_SLURM_QUEUE is required"`), which is the documented fail-closed
behaviour.

---

## 4.  Reading the result

A successful run leaves three artifacts:

    nextflow-cluster-smoke-results/
    ├── probe-<jobid>.evidence.json     # cluster resource probe output
    ├── calls-17_69000-69050.vcf{,.manifest.json}
    ├── calls-17_69051-69100.vcf{,.manifest.json}
    ├── …
    ├── calls-17_69351-69400.vcf{,.manifest.json}
    └── gathered.vcf.gz{,.tbi,.manifest.json}

The full cluster summary lands in
`fastgatk-native/evidence/cluster-evidence.json`.  Use `sacct -j <jobid>`
to correlate each shard's `slurm_job_id` against the ELAPSED/MAXRSS
columns the script records.

If `gathered.vcf.gz.tbi` is missing or non-zero, the run is *not*
complete even if every shard succeeded.  Section 6 applies.

---

## 5.  Promote the result into the release gate

`verify_cluster_smoke.sh` is the single canonical entry point for both
modes.  It writes `fastgatk-native/evidence/cluster-evidence.json` with
`status=pass` and a `cluster_kind` field that disambiguates the run:

| Mode | Env | cluster_kind | sacct per shard? |
|------|-----|--------------|------------------|
| Local proof | `FASTGATK_LOCAL_PROOF=1` (default for `fastgatk-cluster-smoke` CTest) | `local-infrastructure-proof` | no |
| Real SLURM | `FASTGATK_REAL_CLUSTER=1` + `sbatch` on PATH + `SLURM_JOB_ID` set | `slurm` | yes (via `sacct -X -P`) |

Both modes share the same wrapper (`slurm_smoke.sh`), the same scatter
driver (`scatter_gather_smoke.sh`), the same gather contract, and the
same evidence schema.  Local-proof mode just calls the wrapper in its
local-fallback branch and skips `sacct`.

### 5a.  Local proof (the CTest default)

This is the cutover gate the user authorised with "pinned-GATK oracle
+ fixture" (no dual-run shadow).  Re-run anytime the wrapper or
manifest contract changes:

    bash fastgatk-native/scripts/verify_cluster_smoke.sh

The CTest `fastgatk-cluster-smoke` does exactly this — it is registered
with `ENVIRONMENT FASTGATK_LOCAL_PROOF=1` so a workstation passes by
default, and the label `m0_cluster_infra` makes it trivial to filter.

### 5b.  Real SLURM run (production cluster)

On the cluster, after the dispatcher and container images have been
promoted:

    export FASTGATK_REAL_CLUSTER=1
    export FASTGATK_SLURM_QUEUE=compute
    export FASTGATK_SLURM_ACCOUNT=biology
    bash fastgatk-native/scripts/verify_cluster_smoke.sh

The script submits the probe + 8 shards via `sbatch --wait`, collects
`sacct` records, and writes `cluster_kind=slurm` evidence with
non-empty per-shard `sacct` entries.  The dispatcher / runbook can then
audit peak memory and wall-time per shard.

### 5c.  Promote (both modes)

Once `cluster-evidence.json` is on disk and `status == "pass"` (any
`cluster_kind`):

1. Run `python3 fastgatk-native/scripts/compute_progress_score.py`.  It
   will print a `cluster_evidence.ready-to-bump` block with the exact
   gate values to copy into `progress_score.json::workflows::gpu_slurm_cluster`.
   The recommended `e2e_perf` gate is `0.50` for local-proof and `0.65`
   for real-SLURM.
2. Update the `gpu_slurm_cluster` workflow's `gates` and `score` fields,
   and add the new evidence paths to the workflow's `evidence` list.
3. Re-run `compute_progress_score.py` and confirm it still passes (the
   score is validated to 1e-12 against the gate weights).
4. Commit.  Do **not** auto-bump on CI — the snapshot is a human/agent
   artefact and its value is reviewed on every PR.

---

## 6.  Rollback: 5 minutes back to GATK

This is the only section that explicitly violates the "no JVM" principle.
Treat it as a fire-suppression tool, not a steady state.

1. Verify the trigger is a dispatcher reject:

       cat shard-output/calls-*.vcf.manifest.json | jq '.compatibility.fallback_disabled'
       # → true   →   rollback applies
       # → false  →   the failure is somewhere else; consult section 7

2. Stop the Nextflow run and drain in-flight shards.

       nextflow log <run-name> -F   # confirm no NEW tasks accepted
       nextflow kill <run-name>     # sends SIGTERM to the driver

3. Stage the GATK 4.6.2.0 jar on every worker.  This is the *only*
   moment the cluster has a JVM on disk.

       # From a login node with the cluster file system mounted:
       rsync -a third_party/jdk17/                    /opt/jdk17/
       rsync -a third_party/gatk-package/gatk-4.6.2.0 /opt/gatk-oracle/

4. Swap the dispatcher symlink on every worker.

       # Worker-side install (one line, idempotent):
       sudo ln -sfn /opt/gatk-oracle/gatk-package-4.6.2.0/gatk-package-4.6.2.0-local.jar \
           /usr/local/bin/gatk.jar
       sudo ln -sfn /opt/jdk17/bin/java /usr/local/bin/java
       # Re-run the pipeline with the original GATK launcher.

5. Record the rollback reason in `RUNBOOK_INCIDENTS/` with timestamp,
   shard boundary, and the manifest that triggered section 6 step 1.

6. When the rollback is in place, file an issue against `progress_score.json`
   to revert the `gpu_slurm_cluster` bump from section 5.

---

## 7.  Failure codes

The dispatcher / runtime emit a small set of stable failure categories.
Find them in `OutputManifest.compatibility.fallback_reason` and in
stderr JSON lines.

| Code | Where | Meaning | Action |
|------|-------|---------|--------|
| `FALLBACK_DISABLED` (exit 73) | dispatcher | Native has no JVM fallback; the tool/option was deliberately not registered | Section 6 if it is a production path; otherwise close the ticket as not-a-bug |
| `UNSUPPORTED_PARAMETER` | dispatcher | Cloud/Spark option that the project has not implemented | Remove the option from the calling workflow |
| `RESOURCE_EXHAUSTED` | runtime | Native batch controller refused to allocate because cgroup/SLURM limits would be breached | Reduce `--threads` / `--batch-records`, or move to a larger queue |
| `BACKEND_UNAVAILABLE` | runtime | Requested Kokkos execution space (CUDA/HIP/SYCL) is not compiled into this binary | Rebuild with the matching backend, or remove the backend request |
| `NUMERICAL_CONTRACT_FAILURE` | kernels | Native value drifted from GATK beyond the contract tolerance | Capture the shard OutputManifest and open a parity ticket; do **not** roll the dispatcher into production for that tool |
| `OUTPUT_CONTRACT_FAILURE` | runtime | `publish_output_bundle()` rejected the bundle (incomplete manifest, missing sidecar) | Re-run the shard; this is a transient staging failure, not a tool bug |

Any code outside that table is a bug in the dispatcher/runtime and should
be filed with a captured stderr line.

---

## 8.  CI oracle budget

`verify_all.sh` runs every CTest plus every oracle-gated script.  When
`FASTGATK_REQUIRE_GATK_ORACLE=1` is exported, each oracle script invokes
the pinned GATK 4.6.2.0 jar and pays 5–30 s of JVM startup.  The aggregate
default ceiling is **3600 s**; beyond that the CTest layer is allowed to
skip with `status: skip` (not fail) so CI does not flake on transient
JVM-warmup pressure.

Set `FASTGATK_ORACLE_BUDGET_SECONDS=<seconds>` to override.  Setting it
to `0` disables the skip and forces a hard fail.

---

## 9.  Clean build is mandatory

The repository has a documented "stale object SIGSEGV" history when an
incremental CMake build reuses a Kokkos core object that was compiled
under different compile flags (`NEXT_PHASE_TASKS.md` Round 4/5, 2026-09-03).
Cluster nodes **must** rebuild from scratch on every dispatcher upgrade:

    rm -rf fastgatk-native/build fastgatk-native/build-serial \
           fastgatk-runtime/build fastgatk-runtime/build-serial
    bash fastgatk-native/scripts/build_native.sh

A worker running an incremental build against a refreshed dispatcher
registry is the single most common cause of cluster-side segfaults; the
CTest `fastgatk-cluster-smoke` will reject any binary whose build stamp
predates the registry snapshot it shipped with.

---

## 10.  Scope: what this runbook covers

- M0 cluster infrastructure gates (Step 1–6 in this commit).
- `fastgatk-cluster-smoke` CTest and its 8-shard scatter/gather contract.
- Cluster-side dispatcher behaviour for the 48 currently-registered tools.

Out of scope:

- Track 1 (1:1 parity oracle fixes) — see `M1_PLAN.md` once M0 passes.
- Track 2 (real-WGS perf) — see `M2_PLAN.md` once M1 passes.
- Spark/BWA/PathSeq/DRAGEN/SV/RNA whole families — fail-closed per
  `progress_score.json::blocking_boundary` and `GATK_REMAINING_ALGORITHMS.md`.
