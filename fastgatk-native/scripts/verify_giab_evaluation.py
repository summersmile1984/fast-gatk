#!/usr/bin/env python3
"""Real pinned vcfeval regressions; synthetic results are never GIAB accuracy."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import unittest
import uuid

from giab_assets import load_lock
from giab_environment import PinnedTools, ROOT, save_json
from giab_evaluation import score_callset


class ComparatorSemantics(unittest.TestCase):
    data_root: Path
    output: Path

    @classmethod
    def setUpClass(cls):
        cls.output.mkdir(parents=True, exist_ok=False)
        cls.inputs = cls.output / 'inputs'
        cls.inputs.mkdir()
        cls.environment = json.loads((cls.data_root / 'environment-manifest.json').read_text())
        cls.tools = PinnedTools(cls.data_root, load_lock()['images'], cls.environment['cpus'],
                                cls.environment['memory_bytes'], mounts=[cls.output])
        cls.tools.logs = cls.output / 'preparation-logs'
        rng = random.Random(501)
        sequence = list(''.join(rng.choice('ACGT') for _ in range(600)))
        sequence[200:205] = 'ACGTA'
        cls.sequence = ''.join(sequence)
        cls.reference = cls.inputs / 'reference.fa'
        cls.reference.write_text('>chr20\n' + cls.sequence + '\n')
        cls.tools('samtools', ['faidx', str(cls.reference)], 'reference-index')
        cls.tools('samtools', ['dict', '-o', str(cls.reference.with_suffix('.dict')), str(cls.reference)], 'reference-dictionary')
        cls.scope = cls.inputs / 'scope.bed'
        cls.scope.write_text('chr20\t0\t600\n')
        cls.bed = cls.inputs / 'evaluation.bed'
        cls.bed.write_text('chr20\t0\t350\n')
        left, right = cls.inputs / 'left.bed', cls.inputs / 'right.bed'
        left.write_text('chr20\t0\t310\n')
        right.write_text('chr20\t290\t350\n')
        cls.strata = cls.inputs / 'strata.tsv'
        cls.strata.write_text(f'left\t{left}\nright\t{right}\n')
        cls.header = ('##fileformat=VCFv4.2\n##contig=<ID=chr20,length=600>\n'
                      '##FILTER=<ID=q10,Description="Synthetic low quality">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')

    def variant(self, position, ref=None, alt=None, gt='0/1', filter_value='PASS'):
        ref = ref or self.sequence[position - 1]
        alt = alt or {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[ref]
        self.assertEqual(self.sequence[position - 1:position - 1 + len(ref)], ref)
        return f'chr20\t{position}\t.\t{ref}\t{alt}\t60\t{filter_value}\t.\tGT\t{gt}\n'

    def score(self, name, truth, query):
        paths = {}
        for side, rows in [('truth', truth), ('query', query)]:
            plain = self.inputs / (name + '-' + side + '.vcf')
            plain.write_text(self.header + ''.join(rows))
            path = Path(str(plain) + '.gz')
            self.tools('bcftools', ['view', '-Oz', '-o', str(path), str(plain)], name + '-' + side + '-compress')
            self.tools('bcftools', ['index', '-t', str(path)], name + '-' + side + '-index')
            paths[side] = path
        result = score_callset(**paths, reference=self.reference, scope_bed=self.scope,
                               evaluation_bed=self.bed, data_root=self.data_root, output=self.output / name,
                               environment=self.environment, experiment_kind='synthetic_smoke',
                               stratification_tsv=self.strata, mounts=[self.output])
        self.assertEqual(result['accuracy_status'], 'measured', result.get('errors'))
        return {(row['variant_type'], row['filter'], row['subset']): row
                for row in result['metrics'] if row['subtype'] == '*' and row['genotype'] == '*'}

    def test_complex_equivalent_deletions_use_query_side_precision(self):
        truth = [self.variant(201, 'ACGTA', 'AGA', gt='1/1')]
        extra = self.sequence[249:251]
        query = [self.variant(201, 'AC', 'A', gt='1/1'), self.variant(203, 'GT', 'G', gt='1/1'),
                 self.variant(250, extra, extra[0], gt='1/1')]
        row = self.score('complex', truth, query)['INDEL', 'PASS', '*']
        self.assertEqual(row['truth_fn'], 0)
        self.assertEqual(row['query_fp'], 1)
        self.assertEqual(row['recall'], 1)
        self.assertEqual((row['truth_tp'], row['query_tp']), (1, 2))
        self.assertAlmostEqual(row['precision'], 2 / 3)

    def test_wrong_genotype_is_not_allele_match_success(self):
        row = self.score('genotype', [self.variant(301)], [self.variant(301, gt='1/1')])['SNP', 'PASS', '*']
        self.assertEqual((row['truth_tp'], row['truth_fn'], row['query_tp'], row['query_fp']), (0, 1, 0, 1))
        self.assertEqual((row['fp_gt'], row['fp_al']), (1, 0))
        self.assertEqual((row['precision'], row['recall']), (0, 0))
        self.assertIsNone(row['f1'])

    def test_outside_confident_bed_and_overlapping_strata(self):
        rows = self.score('outside-and-filter', [self.variant(301)],
                          [self.variant(301), self.variant(325, filter_value='q10'), self.variant(401)])
        passed, all_calls = rows['SNP', 'PASS', '*'], rows['SNP', 'ALL', '*']
        self.assertEqual((passed['query_fp'], passed['query_unk'], passed['query_tp']), (0, 1, 1))
        self.assertEqual((all_calls['query_fp'], all_calls['query_unk']), (1, 1))
        self.assertEqual(passed['confident_bases'], 350)
        self.assertEqual(rows['SNP', 'PASS', 'left']['truth_tp'], 1)
        self.assertEqual(rows['SNP', 'PASS', 'right']['truth_tp'], 1)
        self.assertEqual(passed['truth_tp'], 1, 'overlapping strata must not inflate the overall count')

    def test_missing_truth_variant_remains_false_negative(self):
        row = self.score('missing', [self.variant(301), self.variant(330)], [self.variant(301)])['SNP', 'PASS', '*']
        self.assertEqual((row['truth_tp'], row['truth_fn']), (1, 1))
        self.assertEqual(row['recall'], .5)
        self.assertEqual(row['precision'], 1)

    def test_empty_query_has_false_negatives_and_undefined_precision(self):
        rows = self.score('empty', [self.variant(201, 'AC', 'A'), self.variant(301)], [])
        for kind in ['SNP', 'INDEL']:
            for filter_value in ['ALL', 'PASS']:
                row = rows[kind, filter_value, '*']
                self.assertEqual((row['truth_tp'], row['truth_fn'], row['query_tp'], row['query_fp']), (0, 1, 0, 0))
                self.assertEqual(row['recall'], 0)
                self.assertIsNone(row['precision'])
                self.assertIsNone(row['f1'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=ROOT / 'testdata/giab-hg002-v5')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    ComparatorSemantics.data_root = args.data_root.resolve()
    ComparatorSemantics.output = (args.output or args.data_root / 'smoke' / ('comparator-semantics-' + uuid.uuid4().hex)).resolve()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ComparatorSemantics)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    save_json(ComparatorSemantics.output / 'verification.json', {
        'experiment_kind': 'synthetic_smoke', 'tests_run': result.testsRun,
        'failures': [(str(test), detail) for test, detail in result.failures],
        'errors': [(str(test), detail) for test, detail in result.errors],
        'successful': result.wasSuccessful(), 'accuracy_claim': 'No HG002/GIAB accuracy inference'})
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
