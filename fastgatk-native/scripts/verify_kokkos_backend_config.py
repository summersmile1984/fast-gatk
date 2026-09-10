#!/usr/bin/env python3
"""Validate a configured Kokkos backend without requiring a GPU.

This is intentionally a configure-time check.  It reads the CMake cache and
the generated Kokkos package metadata, so CUDA/HIP/SYCL builds can be checked
in a CPU-only CI runner for selection, compiler provenance, and architecture
consistency.  It never creates a device context and never treats a missing
GPU as a test failure; an explicitly requested backend that was not actually
configured is a hard failure instead.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


BACKENDS = ("SERIAL", "OPENMP", "CUDA", "HIP", "SYCL")
DEVICE_BACKENDS = ("CUDA", "HIP", "SYCL")
HOST_ARCHITECTURES = {
    "NATIVE", "AMDAVX", "ARMV80", "ARMV81", "ARMV84", "ARMV84_SVE",
    "ARMV8_THUNDERX", "ARMV8_THUNDERX2", "ARMV9_GRACE", "A64FX", "SNB",
    "HSW", "BDW", "ICL", "ICX", "SKL", "SKX", "KNC", "KNL", "SPR",
    "POWER8", "POWER9", "ZEN", "ZEN2", "ZEN3", "ZEN4", "ZEN5",
    "RISCV_SG2042", "RISCV_RVA22V", "RISCV_U74MC",
}


def parse_cache(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("//") or ":" not in line or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        key, _, _type = key_type.partition(":")
        if key and _type:
            values[key] = value
    return values


def is_on(value: str | None) -> bool:
    return (value or "").strip().upper() in {"ON", "TRUE", "YES", "1"}


def executable_probe(value: str | None) -> dict[str, object]:
    if not value:
        return {"configured": False, "path": "", "executable": False}
    command = value.strip().split()[0]
    resolved = Path(command)
    if not resolved.is_absolute():
        located = shutil.which(command)
        resolved = Path(located) if located else resolved
    return {"configured": True, "path": str(resolved),
            "executable": resolved.is_file() and os.access(resolved, os.X_OK)}


def compiler_version(command: str | None) -> str:
    if not command:
        return ""
    try:
        result = subprocess.run([command, "--version"], text=True,
                                capture_output=True, check=False, timeout=3)
    except (OSError, subprocess.SubprocessError):
        return ""
    return (result.stdout or result.stderr).splitlines()[0][:200] if result.returncode == 0 else ""


def package_metadata(path: Path) -> tuple[list[str], str, str]:
    if not path.is_file():
        return [], "", ""
    text = path.read_text(encoding="utf-8")
    devices_match = re.search(r'set\(Kokkos_DEVICES\s+([^\)]*)\)', text)
    id_match = re.search(r'set\(Kokkos_CXX_COMPILER_ID\s+"([^"]*)"\)', text)
    version_match = re.search(r'set\(Kokkos_CXX_COMPILER_VERSION\s+"([^"]*)"\)', text)
    devices_text = devices_match.group(1).strip().strip('"') if devices_match else ""
    devices = [item for item in re.split(r"[;\s]+", devices_text) if item]
    return devices, id_match.group(1) if id_match else "", version_match.group(1) if version_match else ""


def fail(report: dict[str, object], *reasons: str) -> int:
    report["status"] = "fail"
    report["reasons"] = list(reasons)
    print(json.dumps(report, sort_keys=True))
    return 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--kokkos-config", type=Path, required=True)
    args = parser.parse_args()
    if not args.cache.is_file():
        raise SystemExit(f"CMake cache not found: {args.cache}")

    cache = parse_cache(args.cache)
    configured_backend = cache.get("FASTGATK_KOKKOS_BACKEND", "AUTO").upper()
    enabled = {backend: is_on(cache.get(f"Kokkos_ENABLE_{backend}"))
               for backend in BACKENDS}
    enabled_backends = [backend for backend, value in enabled.items() if value]
    enabled_devices = [backend for backend in DEVICE_BACKENDS if enabled[backend]]
    arch_values = {key.removeprefix("Kokkos_ARCH_"): is_on(value)
                   for key, value in cache.items() if key.startswith("Kokkos_ARCH_")}
    enabled_arches = sorted(key for key, value in arch_values.items() if value)
    host_arches = sorted(key for key in enabled_arches if key in HOST_ARCHITECTURES)
    gpu_arches = sorted(key for key in enabled_arches if key not in HOST_ARCHITECTURES)
    devices_in_package, compiler_id, compiler_version_value = package_metadata(args.kokkos_config)

    report: dict[str, object] = {
        "status": "pass",
        "configured_backend": configured_backend,
        "enabled_backends": enabled_backends,
        "device_backends": enabled_devices,
        "expected_execution_space": (
            "Cuda" if enabled["CUDA"] else
            "HIP" if enabled["HIP"] else
            "SYCL" if enabled["SYCL"] else
            "OpenMP" if enabled["OPENMP"] else
            "Serial" if enabled["SERIAL"] else "Host"),
        "architectures": {"enabled": enabled_arches, "host": host_arches,
                           "gpu": gpu_arches},
        "compiler": {
            "cxx": executable_probe(cache.get("CMAKE_CXX_COMPILER")),
            "cuda": executable_probe(cache.get("CMAKE_CUDA_COMPILER")),
            "hip": executable_probe(cache.get("CMAKE_HIP_COMPILER")),
            "version": compiler_version(cache.get("CMAKE_CXX_COMPILER")),
            "kokkos_id": compiler_id,
            "kokkos_version": compiler_version_value,
        },
        "kokkos_package_devices": devices_in_package,
        "device_runtime_probe": "not-run (configuration-only gate)",
    }

    if not enabled_backends:
        return fail(report, "no Kokkos execution space is enabled")
    if len(host_arches) > 1:
        return fail(report, "more than one host architecture is enabled")
    if configured_backend != "AUTO":
        valid_aliases = set(BACKENDS) | {f"{device}_OPENMP" for device in DEVICE_BACKENDS}
        if configured_backend not in valid_aliases:
            return fail(report, f"unknown FASTGATK_KOKKOS_BACKEND={configured_backend}")
        required = set()
        if configured_backend.endswith("_OPENMP"):
            required.update({"OPENMP", configured_backend.removesuffix("_OPENMP")})
        else:
            required.add(configured_backend)
        missing = sorted(backend for backend in required if not enabled[backend])
        if missing:
            return fail(report, "explicit backend was not enabled: " + ",".join(missing))
        if configured_backend in DEVICE_BACKENDS and not enabled["SERIAL"] and not enabled["OPENMP"]:
            return fail(report, "device-only build has no Host execution space")
        unexpected_devices = sorted(set(enabled_devices) -
                                    ({configured_backend.removesuffix("_OPENMP")} if configured_backend.endswith("_OPENMP") else
                                     ({configured_backend} if configured_backend in DEVICE_BACKENDS else set())))
        if unexpected_devices:
            return fail(report, "explicit backend also enabled unexpected devices: " + ",".join(unexpected_devices))
    if args.kokkos_config.is_file():
        package_upper = {item.upper() for item in devices_in_package}
        missing_package = sorted(set(enabled_backends) - package_upper)
        if missing_package:
            return fail(report, "generated Kokkos package omits enabled backends: " + ",".join(missing_package))
    elif configured_backend != "AUTO":
        return fail(report, "generated Kokkos package metadata is missing")
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
