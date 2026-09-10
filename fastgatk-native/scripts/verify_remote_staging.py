#!/usr/bin/env python3
"""Contract test for dispatcher URI input staging and fail-closed outputs."""

from __future__ import annotations

import http.server
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import threading
import time
from urllib.parse import unquote, urlsplit


def write_reference(path: pathlib.Path) -> None:
    sequence = "ACGTACGTGGCCNNNNATAT"
    header = f">chr1\n".encode()
    body = (sequence + "\n").encode()
    with path.open("wb") as output:
        output.write(header)
        output.write(body)
    path.with_suffix(path.suffix + ".fai").write_text(
        f"chr1\t{len(sequence)}\t{len(header)}\t{len(sequence)}\t{len(sequence) + 1}\n",
        encoding="utf-8",
    )
    path.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\n",
        encoding="utf-8",
    )


_active_requests = 0
_max_active_requests = 0
_request_lock = threading.Lock()


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, _format: str, *_args: object) -> None:
        return

    def do_GET(self) -> None:  # noqa: N802 - stdlib handler hook
        global _active_requests, _max_active_requests
        with _request_lock:
            _active_requests += 1
            _max_active_requests = max(_max_active_requests, _active_requests)
        try:
            # Make overlap observable even on a very fast local filesystem.
            time.sleep(0.05)
            super().do_GET()
        finally:
            with _request_lock:
                _active_requests -= 1

    def do_PUT(self) -> None:  # noqa: N802 - stdlib handler hook
        length = int(self.headers.get("Content-Length", "0"))
        payload = self.rfile.read(length)
        target = pathlib.Path(self.directory) / unquote(urlsplit(self.path).path.lstrip("/"))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        self.send_response(200)
        self.end_headers()


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    dispatcher = root / "fastgatk-native/dispatcher/fastgatk"
    annotate = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-annotate-intervals"
    with tempfile.TemporaryDirectory(prefix="fastgatk-remote-staging-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        write_reference(reference)
        intervals = work / "targets.interval_list"
        intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:20\nchr1\t1\t8\t+\tbin1\n",
            encoding="utf-8",
        )
        handler = lambda *args, **kwargs: QuietHandler(
            *args, directory=str(work), **kwargs
        )
        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            base = f"http://127.0.0.1:{server.server_port}"
            sys.path.insert(0, str(pathlib.Path(dispatcher).parent))
            import fastgatk as dispatcher_module

            resource_args = [
                "--resource:fixture,training=true,truth=true",
                f"{base}/targets.interval_list?sig=fixture",
            ]
            staged_resource_args, resource_sources = dispatcher_module.stage_remote_inputs(resource_args)
            if staged_resource_args[1] == resource_args[1] or resource_sources != [resource_args[1]]:
                raise AssertionError({"resource_staging": staged_resource_args,
                                      "sources": resource_sources})
            output = work / "annotated.tsv"
            environment = os.environ.copy()
            environment.update({
                "FASTGATK_REMOTE_CACHE": str(work / "cache"),
                "FASTGATK_ANNOTATEINTERVALS_BINARY": str(annotate),
                "FASTGATK_REMOTE_DOWNLOAD_THREADS": "2",
            })
            invalid_environment = environment.copy()
            invalid_environment["FASTGATK_REMOTE_DOWNLOAD_THREADS"] = "0"
            invalid = subprocess.run([
                str(dispatcher), "AnnotateIntervals", "-R", f"{base}/reference.fasta?sig=fixture",
                "-L", f"{base}/targets.interval_list?sig=fixture", "-O", str(output),
                "--interval-merging-rule", "OVERLAPPING_ONLY",
            ], check=False, text=True, capture_output=True, env=invalid_environment)
            if invalid.returncode != 2:
                raise AssertionError({"invalid_worker_returncode": invalid.returncode,
                                      "stderr": invalid.stderr})
            invalid_error = json.loads(invalid.stderr.splitlines()[-1])
            if invalid_error.get("category") != "BAD_ARGUMENT":
                raise AssertionError({"invalid_worker_error": invalid_error})
            result = subprocess.run([
                str(dispatcher), "AnnotateIntervals", "-R", f"{base}/reference.fasta?sig=fixture",
                "-L", f"{base}/targets.interval_list?sig=fixture", "-O", str(output),
                "--interval-merging-rule", "OVERLAPPING_ONLY",
            ], check=True, text=True, capture_output=True, env=environment)
            summary = json.loads(result.stdout.splitlines()[-1])
            rows = [line for line in output.read_text(encoding="utf-8").splitlines()
                    if line and not line.startswith("@") and not line.startswith("CONTIG")]
            if summary["status"] != "prototype" or rows != ["chr1\t1\t8\t0.500000"]:
                raise AssertionError({"summary": summary, "rows": rows})
            cached = list((work / "cache").glob("*"))
            if len(cached) < 3:  # primary FASTA + required .fai + canonical .dict
                raise AssertionError({"cache_entries": [str(path) for path in cached]})

            remote_output = subprocess.run([
                str(dispatcher), "AnnotateIntervals", "-R", str(reference),
                "-L", str(intervals), "-O", f"{base}/remote-output.tsv",
            ], text=True, capture_output=True, env=environment)
            if remote_output.returncode != 69:
                raise AssertionError({"returncode": remote_output.returncode,
                                      "stderr": remote_output.stderr})
            error = json.loads(remote_output.stderr.splitlines()[-1])
            if error["category"] != "BACKEND_UNAVAILABLE":
                raise AssertionError(error)

            environment["FASTGATK_REMOTE_OUTPUT_MODE"] = "upload"
            uploaded = subprocess.run([
                str(dispatcher), "AnnotateIntervals", "-R", str(reference),
                "-L", str(intervals), "-O", f"{base}/uploaded-output.tsv",
                "--output-manifest", f"{base}/uploaded-output.manifest.json",
                "--interval-merging-rule", "OVERLAPPING_ONLY",
            ], text=True, capture_output=True, env=environment)
            if uploaded.returncode != 0:
                raise AssertionError({"returncode": uploaded.returncode,
                                      "stderr": uploaded.stderr})
            uploaded_rows = [line for line in
                             (work / "uploaded-output.tsv").read_text(encoding="utf-8").splitlines()
                             if line and not line.startswith("@") and not line.startswith("CONTIG")]
            if uploaded_rows != ["chr1\t1\t8\t0.500000"]:
                raise AssertionError({"uploaded_rows": uploaded_rows})
            manifest = (work / "uploaded-output.manifest.json").read_text(encoding="utf-8")
            if f"{base}/uploaded-output.tsv" not in manifest:
                raise AssertionError({"manifest": manifest})
            if _max_active_requests < 2:
                raise AssertionError({"max_concurrent_remote_gets": _max_active_requests})
            print(json.dumps({"status": "pass", "cached_entries": len(cached),
                              "max_concurrent_remote_gets": _max_active_requests,
                              "remote_output_fail_closed": True,
                              "remote_output_upload": True}, sort_keys=True))
        finally:
            server.shutdown()
            server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
