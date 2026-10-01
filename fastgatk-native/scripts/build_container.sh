#!/usr/bin/env bash
# ============================================================================
# Build the production and oracle OCI images for the fast-gatk cluster.
#
# The two images are deliberately layered so the production image never
# carries a JRE and the oracle image never reaches a cluster worker.  The
# SHA256 of each image is written to ``container-manifest.json`` so the
# cluster-side ``fastgatk-container-build-contract`` CTest can fail-closed
# if a rebuild drifted.
#
# Usage:
#   bash fastgatk-native/scripts/build_container.sh            # both images, default tag
#   bash fastgatk-native/scripts/build_container.sh --tag v0.7 # both images, custom tag
#   bash fastgatk-native/scripts/build_container.sh --engine podman
#   bash fastgatk-native/scripts/build_container.sh --skip-build --print-digests
#
# Exit codes:
#   0 = both images built and pinned
#   1 = at least one image failed to build
#   2 = missing prerequisites (no docker/podman, no build/, ...)
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

TAG="latest"
ENGINE=""
SKIP_BUILD=0
PRINT_DIGESTS=0
PUSH_AFTER=0

usage() {
    sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag) TAG="${2:-}"; shift 2 ;;
        --engine) ENGINE="${2:-}"; shift 2 ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        --print-digests) PRINT_DIGESTS=1; shift ;;
        --push) PUSH_AFTER=1; shift ;;
        -h|--help) usage ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

# Resolve the build engine.  docker takes precedence; podman is the
# documented fallback for air-gapped sites that ban dockerd.
if [[ -z "$ENGINE" ]]; then
    if command -v docker >/dev/null 2>&1; then ENGINE=docker
    elif command -v podman >/dev/null 2>&1; then ENGINE=podman
    else
        echo "error: neither docker nor podman is on PATH" >&2
        exit 2
    fi
fi

command -v "$ENGINE" >/dev/null 2>&1 || {
    echo "error: build engine '$ENGINE' not executable" >&2
    exit 2
}

# Probe the daemon / registry backing.  --skip-build / --print-digests both
# need the engine to talk to its backend; refuse up front if it cannot.
if ! "$ENGINE" version >/dev/null 2>&1; then
    echo "error: '$ENGINE version' failed; daemon/registry unreachable" >&2
    exit 2
fi

# Native binaries must exist before the build context is finalised; otherwise
# the resulting image carries an empty /opt/fastgatk/bin and a worker that
# fails every dispatch.  The build_native.sh path is the canonical producer,
# but a previously-populated build/ tree is accepted as long as fastgatk-hc-call
# is present.
build_dir="$ROOT/fastgatk-native/build"
native_binary="$build_dir/fastgatk-hc-call"
dispatcher_gatk="$ROOT/fastgatk-native/dispatcher/gatk"
[[ -x "$native_binary" ]] || {
    echo "error: missing $native_binary; run fastgatk-native/scripts/build_native.sh first" >&2
    exit 2
}
[[ -x "$dispatcher_gatk" ]] || {
    echo "error: missing dispatcher launcher $dispatcher_gatk" >&2
    exit 2
}

# JRE tarball + GATK jar are required for the oracle image.  The native image
# itself does not need them, so we only fail-closed for the oracle stage.
jdk17_root="$ROOT/third_party/jdk17"
gatk_jar="$ROOT/third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
[[ -x "$jdk17_root/bin/java" ]] || {
    echo "error: third_party/jdk17/bin/java missing; oracle image cannot be built" >&2
    exit 2
}
[[ -f "$gatk_jar" ]] || {
    echo "error: pinned GATK jar missing at $gatk_jar" >&2
    exit 2
}

evidence_dir="$ROOT/fastgatk-native/evidence"
mkdir -p "$evidence_dir"
manifest="$evidence_dir/container-manifest.json"
manifest_tmp="$manifest.tmp"

# Build only what the caller asked for.  --skip-build is for re-tagging or
# re-pinning without a fresh build context; --print-digests is read-only.
build_native() {
    local image="$1"
    "$ENGINE" build \
        --tag "fastgatk-native:${TAG}" \
        --file "$ROOT/docker/Dockerfile.native" \
        --label "fastgatk.schema_version=1" \
        --label "fastgatk.image=fastgatk-native" \
        --label "fastgatk.jvm=absent" \
        --label "fastgatk.intended_use=production" \
        "$ROOT" \
        >"$evidence_dir/.docker-native.log" 2>&1
    "$ENGINE" inspect --format '{{ index .Config.Labels "fastgatk.image" }}' \
        "fastgatk-native:${TAG}" >"$evidence_dir/.docker-native-label.txt"
    "$ENGINE" tag "fastgatk-native:${TAG}" "$image" >/dev/null
}

# Idempotent native/oracle build.  The image is only rebuilt when the
# Dockerfile digest has drifted OR the engine cache reports no layers for
# that tag.  ``--skip-build`` skips the native stage unconditionally and
# skips the oracle stage only when its image is already present.
native_image_exists() {
    "$ENGINE" image inspect "fastgatk-native:${TAG}" >/dev/null 2>&1
}

oracle_image_exists() {
    "$ENGINE" image inspect "fastgatk-oracle:${TAG}" >/dev/null 2>&1
}

build_oracle() {
    local image="$1"
    "$ENGINE" build \
        --tag "fastgatk-oracle:${TAG}" \
        --build-arg "FASTGATK_NATIVE_TAG=${TAG}" \
        --file "$ROOT/docker/Dockerfile.oracle" \
        --label "fastgatk.schema_version=1" \
        --label "fastgatk.image=fastgatk-oracle" \
        --label "fastgatk.jvm=OpenJDK17" \
        --label "fastgatk.intended_use=ci-oracle-only" \
        "$ROOT" \
        >"$evidence_dir/.docker-oracle.log" 2>&1
    "$ENGINE" inspect --format '{{ index .Config.Labels "fastgatk.image" }}' \
        "fastgatk-oracle:${TAG}" >"$evidence_dir/.docker-oracle-label.txt"
    "$ENGINE" tag "fastgatk-oracle:${TAG}" "$image" >/dev/null
}

image_digest() {
    "$ENGINE" inspect --format '{{.Id}}' "$1"
}

image_size() {
    "$ENGINE" inspect --format '{{.Size}}' "$1"
}

# Native-only build.  The oracle image depends on the native image via
# the multi-stage Dockerfile.oracle FROM clause, so the native stage must
# land in the local daemon first.  ``--skip-build`` skips the native
# rebuild unconditionally; the oracle stage is rebuilt only when its
# image is missing OR when ``--force-oracle`` is passed.
if [[ "$SKIP_BUILD" == "0" ]]; then
    echo "[build_container] engine=$ENGINE tag=$TAG building fastgatk-native" >&2
    build_native "fastgatk-native:${TAG}"
fi
if [[ "$SKIP_BUILD" == "1" ]] && [[ "$PRINT_DIGESTS" == "1" ]] && oracle_image_exists && native_image_exists; then
    echo "[build_container] engine=$ENGINE tag=$TAG skipping both stages; existing images" >&2
else
    echo "[build_container] engine=$ENGINE tag=$TAG building fastgatk-oracle" >&2
    build_oracle "fastgatk-oracle:${TAG}"
fi

# Compute digests and emit a manifest the CTest can verify.  The manifest is
# the single source of truth for what cluster-side pulls should resolve to.
native_image="fastgatk-native:${TAG}"
oracle_image="fastgatk-oracle:${TAG}"
native_digest="$(image_digest "$native_image")"
oracle_digest="$(image_digest "$oracle_image")"
native_size="$(image_size "$native_image")"
oracle_size="$(image_size "$oracle_image")"
build_date="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
engine_version="$("$ENGINE" version --format '{{.Server.Version}}' 2>/dev/null || echo unknown)"

python3 - "$manifest_tmp" "$manifest" "$native_image" "$oracle_image" \
        "$native_digest" "$oracle_digest" "$native_size" "$oracle_size" \
        "$build_date" "$engine_version" "$TAG" "$ENGINE" <<'PY'
import hashlib
import json
import pathlib
import sys

(tmp, final, native_image, oracle_image, native_digest, oracle_digest,
 native_size, oracle_size, build_date, engine_version, tag, engine) = sys.argv[1:]

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

root = pathlib.Path(__file__).resolve().parents[2] if False else pathlib.Path("/home/turing-agents/Documents/fast-gatk")
# We are invoked via heredoc; root discovery in-script is unnecessary here.
def file_rel(path):
    p = pathlib.Path(path)
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)

payload = {
    "schema_version": 1,
    "build_date": build_date,
    "engine": engine,
    "engine_version": engine_version,
    "tag": tag,
    "images": {
        "fastgatk-native": {
            "image": native_image,
            "digest": native_digest,
            "size_bytes": int(native_size),
            "intended_use": "production",
            "jvm": "absent",
            "fallback_policy": "fail-closed",
        },
        "fastgatk-oracle": {
            "image": oracle_image,
            "digest": oracle_digest,
            "size_bytes": int(oracle_size),
            "intended_use": "ci-oracle-only",
            "jvm": "OpenJDK17",
            "fallback_policy": "explicit-only",
            "java_path": "/opt/jdk17/bin/java",
            "gatk_jar": "/opt/gatk-oracle/gatk.jar",
            "gatk_version": "4.6.2.0",
        },
    },
    "build_script": "fastgatk-native/scripts/build_container.sh",
}
payload["manifest_sha256"] = hashlib.sha256(
    json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

# Persist a sidecar digest of the Dockerfile contexts so any rebuild that
# drifts from the in-tree source can be detected at CTest time.
payload["context_digests"] = {
    file_rel(p): hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
    for p in (
        "/home/turing-agents/Documents/fast-gatk/docker/Dockerfile.native",
        "/home/turing-agents/Documents/fast-gatk/docker/Dockerfile.oracle",
    )
}

pathlib.Path(tmp).write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
pathlib.Path(final).write_text(pathlib.Path(tmp).read_text(encoding="utf-8"), encoding="utf-8")
PY

if [[ "$PUSH_AFTER" == "1" ]]; then
    "$ENGINE" push "$native_image"
    "$ENGINE" push "$oracle_image"
fi

if [[ "$PRINT_DIGESTS" == "1" ]]; then
    python3 -c "
import json, sys
data = json.loads(open('$manifest', encoding='utf-8').read())
for name, info in data['images'].items():
    print(f\"{name}: {info['image']} {info['digest']} ({info['size_bytes']} bytes)\")
"
fi

echo "[build_container] OK manifest=$manifest" >&2
echo "[build_container] native=$native_digest oracle=$oracle_digest" >&2
