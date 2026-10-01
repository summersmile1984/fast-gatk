#!/usr/bin/env python3
"""Small real pinned-tool checks of input conservation and strict reference gates."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from giab_assets import load_lock
from giab_environment import PinnedTools, ROOT
from giab_preflight import _audit_sam, _benchmark_bed, _truth_scope, _validate_vcf


class PreflightContract(unittest.TestCase):
    def setUp(self):
        (ROOT / 'work').mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix='giab-preflight-', dir=ROOT / 'work')
        self.root = Path(self.temporary.name)
        self.prepared = self.root / 'prepared'
        self.prepared.mkdir()
        self.tools = PinnedTools(self.root, load_lock()['images'], list(range(16)), 64 * 1024**3)
        self.reference = {'chr20': {'length': 100}}
        self.header = ('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:100\n'
                       '@RG\tID:rg\tSM:HG002\tLB:lib\tPL:ILLUMINA\n')

    def tearDown(self):
        self.temporary.cleanup()

    def test_duplicate_reset_preserves_qualities_mapping_and_other_tags(self):
        source = self.root / 'source.sam'
        source.write_text(self.header + 'read:1\t1024\tchr20\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg\tNM:i:1\tDT:Z:SQ\tDI:i:3\tDS:i:2\n')
        common = self.root / 'common.bam'
        self.tools('samtools', ['view', '-b', '--remove-flags', '0x400', '-x', 'DT', '-x', 'DI', '-x', 'DS',
                               '-o', str(common), str(source)], 'reset')
        decoded = self.root / 'common.sam'
        self.tools('samtools', ['view', '-h', str(common)], 'decode', stdout=decoded)
        before = _audit_sam(source, self.root, 'source')
        after = _audit_sam(decoded, self.root, 'common', cleared=True)
        self.assertEqual(after['status'], 'eligible')
        self.assertEqual(before['canonical_ordered_sha256'], after['canonical_ordered_sha256'])
        altered = self.root / 'altered.sam'
        altered.write_text(decoded.read_text().replace('\tIIII\t', '\tHIII\t'))
        changed = _audit_sam(altered, self.root, 'altered', cleared=True)
        self.assertNotEqual(before['canonical_ordered_sha256'], changed['canonical_ordered_sha256'])

    def test_bed_overlap_is_union_not_double_counted(self):
        bed = self.root / 'benchmark.bed'
        bed.write_text('chr20\t0\t10\nchr20\t5\t15\nchr20\t15\t20\n')
        result = _benchmark_bed(bed, self.reference, ('chr20',), self.prepared, self.root)['chr20']
        self.assertEqual(result['confident_bases'], 20)
        self.assertEqual(Path(result['evaluation_bed']).read_text(), 'chr20\t0\t20\n')
        self.assertEqual(Path(result['scope_bed']).read_text(), 'chr20\t0\t100\n')

    def test_truth_ref_mismatch_outside_confident_bed_still_blocks_scope(self):
        reference = self.prepared / 'reference.fa'
        reference.write_text('>chr20\n' + 'A' * 100 + '\n')
        self.tools('samtools', ['faidx', str(reference)], 'fai')
        bed = self.root / 'benchmark.bed'
        bed.write_text('chr20\t0\t10\n')
        _benchmark_bed(bed, self.reference, ('chr20',), self.prepared, self.root)
        plain = self.root / 'truth.vcf'
        plain.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr20,length=100>\n'
                         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
                         'chr20\t80\t.\tC\tG\t60\tPASS\t.\tGT\t0/1\n')
        raw = self.root / 'truth.vcf.gz'
        self.tools('bcftools', ['view', '-Oz', '-o', str(raw), str(plain)], 'compress')
        self.tools('bcftools', ['index', '-t', str(raw)], 'index')
        full = _validate_vcf(raw, self.reference, sample='HG002')
        with self.assertRaises(RuntimeError):
            _truth_scope(raw, full, 'chr20', self.reference, self.prepared, self.root, self.tools)
        result = json.loads((self.root / 'preparation-logs/chr20-truth-ref-check.command.json').read_text())
        self.assertNotEqual(result['returncode'], 0)
        # Alleles are preserved for diagnosis, not repaired or BED-filtered away.
        retained = _validate_vcf(self.prepared / 'chr20/truth.vcf.gz', self.reference, sample='HG002')
        self.assertEqual(retained['records_sha256'], full['records_sha256'])


if __name__ == '__main__':
    unittest.main()
