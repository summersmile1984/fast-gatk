#!/usr/bin/env bash
# ============================================================================
# Thin wrapper around verify_cluster_smoke.sh for the local-infrastructure
# proof.  Kept as a separate CTest target so a caller that wants the longer
# name can still drive the same code path.
#
# Exit codes match verify_cluster_smoke.sh exactly (0 = pass, 1 = fail,
# 2 = fail-closed).
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

export FASTGATK_LOCAL_PROOF=1
export FASTGATK_NATIVE_BUILD="${FASTGATK_NATIVE_BUILD:-$ROOT/fastgatk-native/build}"

bash "$ROOT/fastgatk-native/scripts/verify_cluster_smoke.sh"