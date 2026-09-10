#!/usr/bin/env bash
set -euo pipefail

# Compatibility launcher for existing `gatk HaplotypeCaller ...` processes.
# Keep this stable path for existing Nextflow/SLURM processes, but delegate all
# parsing, registry decisions, @args expansion, and explicit fallback behavior
# to the real dispatcher. There must not be a second, subtly different parser
# here.
project_root=$(cd "$(dirname "$0")/../.." && pwd)
dispatcher="$project_root/fastgatk-native/dispatcher/fastgatk"
[[ -x "$dispatcher" ]] || {
    echo "fastgatk dispatcher unavailable: $dispatcher" >&2
    exit 69
}
exec "$dispatcher" "$@"
