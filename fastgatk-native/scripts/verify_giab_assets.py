#!/usr/bin/env python3
"""Exercise real HTTP resume/publication and archive containment boundaries."""
from __future__ import annotations

import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
import json
from pathlib import Path
import tarfile
import tempfile
import threading
import unittest

import giab_assets as assets


class DownloadContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='giab-download-')
        self.root = Path(self.temporary.name)
        self.payload = b'locked-genomic-object\n' * 100
        self.mode = 'interrupt'
        self.requests = []
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self):
                byte_range = self.headers.get('Range')
                owner.requests.append((byte_range, self.headers.get('If-Match')))
                start = int(byte_range[6:-1]) if byte_range else 0
                body = owner.payload[start:]
                status = 206 if byte_range else 200
                if owner.mode == 'ignore_range':
                    status, start, body = 200, 0, owner.payload
                self.send_response(status)
                self.send_header('ETag', '"changed"' if owner.mode == 'changed_etag' else '"locked"')
                self.send_header('Content-Length', str(len(body)))
                if status == 206:
                    self.send_header('Content-Range', f'bytes {start}-{len(owner.payload)-1}/{len(owner.payload)}')
                self.end_headers()
                if owner.mode == 'interrupt':
                    self.wfile.write(body[:37])
                    self.close_connection = True
                else:
                    self.wfile.write(body)

            def log_message(self, *args):
                pass

        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.worker = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.worker.start()
        self.path = self.root / 'asset.bam'
        self.asset = {'url': f'http://127.0.0.1:{self.server.server_port}/asset',
                      'md5': hashlib.md5(self.payload).hexdigest(), 'size_bytes': len(self.payload), 'large': False}

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.worker.join()
        self.temporary.cleanup()

    def transfer(self):
        return assets._download_asset(self.path, self.asset, manifest=False, small_limit=100000)

    def interrupted(self):
        with self.assertRaises(OSError):
            self.transfer()
        self.assertFalse(self.path.exists(), 'partial input must never be a published asset')
        self.assertEqual(self.path.with_name(self.path.name + '.part').read_bytes(), self.payload[:37])

    def test_interrupted_body_resumes_only_matching_object(self):
        self.interrupted()
        self.mode = 'complete'
        self.transfer()
        self.assertEqual(self.requests[-1], ('bytes=37-', '"locked"'))
        self.assertEqual(self.path.read_bytes(), self.payload)
        self.assertFalse(self.path.with_name(self.path.name + '.part').exists())
        receipt = json.loads(self.path.with_name(self.path.name + '.receipt.json').read_text())
        self.assertEqual(receipt['sha256'], hashlib.sha256(self.payload).hexdigest())

    def test_server_ignoring_range_cannot_publish_or_reuse_partial(self):
        self.interrupted()
        self.mode = 'ignore_range'
        with self.assertRaises(assets.AssetError) as caught:
            self.transfer()
        self.assertEqual(caught.exception.code, 'invalid_range')
        self.assertFalse(self.path.exists())
        self.assertFalse(self.path.with_name(self.path.name + '.part').exists())

    def test_changed_object_cannot_append_to_existing_partial(self):
        self.interrupted()
        self.mode = 'changed_etag'
        with self.assertRaises(assets.AssetError) as caught:
            self.transfer()
        self.assertEqual(caught.exception.code, 'object_changed')
        self.assertFalse(self.path.exists())
        self.assertFalse(self.path.with_name(self.path.name + '.part').exists())

    def test_wrong_checksum_never_publishes_complete_body(self):
        self.mode = 'complete'
        self.asset['md5'] = hashlib.md5(b'wrong-object').hexdigest()
        with self.assertRaises(assets.AssetError) as caught:
            self.transfer()
        self.assertEqual(caught.exception.code, 'checksum_mismatch')
        self.assertFalse(self.path.exists())
        self.assertFalse(self.path.with_name(self.path.name + '.part').exists())


class ArchiveContract(unittest.TestCase):
    def test_path_traversal_and_symlink_never_escape_staging(self):
        for kind in ('traversal', 'symlink'):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(prefix='giab-archive-') as temporary:
                root = Path(temporary)
                archive = root / 'input.tar.gz'
                with tarfile.open(archive, 'w:gz') as tar:
                    entry = tarfile.TarInfo('../escaped' if kind == 'traversal' else 'GRCh38@all/link')
                    if kind == 'symlink':
                        entry.type, entry.linkname = tarfile.SYMTYPE, str(root / 'escaped')
                    else:
                        entry.size = 4
                    tar.addfile(entry, io.BytesIO(b'bad!') if kind == 'traversal' else None)
                entrypoint = 'GRCh38@all/layers.tsv'
                digest = hashlib.md5(b'layer\ttrack.bed\n').hexdigest()
                lock = {'stratification': {'archive_root': 'GRCh38@all', 'entrypoint': entrypoint,
                                           'entrypoint_md5': digest}}
                manifests = {'stratification_manifest': {entrypoint: digest},
                             'stratification_member_manifest': {'layers.tsv': digest}}
                recorded = {key: {'path': str(archive), 'sha256': 'evidence'} for key in
                            ('stratifications', 'stratification_manifest', 'stratification_member_manifest')}
                with self.assertRaises(assets.AssetError):
                    assets._prepare_stratifications(root, lock, manifests, recorded, extract=True)
                self.assertFalse((root / 'escaped').exists())
                self.assertFalse((root / 'prepared/stratifications').exists())


if __name__ == '__main__':
    unittest.main()
