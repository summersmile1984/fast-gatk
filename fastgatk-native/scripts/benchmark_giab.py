#!/usr/bin/env python3
"""HG002/GRCh38 GIAB validation, with fail-closed acquisition and eligibility gates."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from giab_assets import fetch_assets, load_lock, verify_assets
from giab_environment import PinnedTools, preflight_environment, save_json
from giab_preflight import prepare_inputs, prepare_reference_inputs

REFERENCE_ASSETS = {'reference', 'truth', 'truth_index', 'benchmark_bed',
                    'dbsnp', 'dbsnp_index', 'mills', 'mills_index'}


def failure_report(root: Path, gate: str, details: dict) -> dict:
    """No completed execution or biological scores are manufactured for blocked gates."""
    output = root / 'reports/preflight'
    save_json(output / 'failure-details.json', details)
    report = {'schema_version': 1, 'sample': 'HG002', 'truth_release': 'v5.0q',
              'truth_release_status': 'NIST recommended replacement for v4.2.1; draft/work-in-progress',
              'reference': 'GRCh38_no_alt', 'blocked_gate': gate,
              'failure_details': str(output / 'failure-details.json'),
              'input_status': 'blocked', 'execution_status': None, 'run_attempted': False,
              'contract_status': 'not_evaluated', 'accuracy_status': 'blocked',
              'resource_status': 'not_measured', 'scores': None,
              'scopes_not_run': ['chr20', 'autosomes'],
              'downstream_execution': 'not_started_due_to_preflight',
              'full_bam_audit': details.get('gates', {}).get('source_bam', {}).get('status', 'not_evaluated'),
              'quality_provenance': details.get('gates', {}).get('provenance', {}).get('quality_eligibility', 'unknown'),
              'failing_gates': {name: value for name, value in details.get('gates', {}).items()
                               if value.get('status') not in ('eligible', 'completed')},
              'errors': details.get('errors', []),
              'scope_prerequisites': {
                  name: {key: value.get(key) for key in
                         ('input_status', 'reference_eligibility', 'scope_bases', 'confident_bases', 'blocking_gates')}
                  for name, value in details.get('scopes', {}).items()},
              'reproduce': [sys.executable, str(Path(__file__).resolve()), 'fetch',
                            '--data-root', str(root)]}
    environment = root / 'environment-manifest.json'
    if environment.is_file():
        report['environment_manifest'] = str(environment)
        report['resource_status'] = json.loads(environment.read_text()).get('resource_status', 'not_measured')
        report['resource_evidence_scope'] = 'preflight only; no five-stage pipeline resource measurement'
    report['artifacts'] = {name: str(root / name) for name in (
        'asset-manifest.json', 'environment-manifest.json', 'build-manifest.json',
        'reference-input-manifest.json', 'input-manifest.json',
        'preflight/ref-check-correction.json', 'preflight/reference-mismatch-diagnostic.json')
        if (root / name).is_file()}
    diagnostic = root / 'preflight/reference-mismatch-diagnostic.json'
    if diagnostic.is_file():
        report['first_reference_mismatch'] = json.loads(diagnostic.read_text())
    stratification_manifest = root / 'prepared/stratifications-manifest.json'
    if stratification_manifest.is_file():
        stratifications = json.loads(stratification_manifest.read_text())
        report['stratifications'] = {
            'manifest': str(stratification_manifest),
            'members': len(stratifications['members']), 'layers': len(stratifications['layers']),
            'source_sha256': stratifications['source_sha256'],
            'validation': stratifications['bed_validation'],
            'reference_eligibility': 'pending', 'used_for_scoring': False}
    verification = root / 'preflight/verification-summary.json'
    if verification.is_file():
        report['verification'] = json.loads(verification.read_text())
    save_json(output / 'giab_summary.json', report)
    with (output / 'giab_summary.tsv').open('w') as table:
        table.write('sample\tscope\tinput_status\trun_attempted\tcontract_status\taccuracy_status'
                    '\treference_eligibility\tconfident_bases\tblocked_gate\n')
        for scope in ('chr20', 'autosomes'):
            scoped = report['scope_prerequisites'].get(scope, {})
            table.write(f"HG002\t{scope}\t{scoped.get('input_status', 'blocked')}\tfalse\tnot_evaluated\tblocked\t"
                        f"{scoped.get('reference_eligibility', 'pending')}\t"
                        f"{scoped.get('confident_bases', 'NA')}\t{gate}\n")
    print(json.dumps({'status': 'blocked', 'gate': gate, 'report': str(output / 'giab_summary.json')}))
    return report


def fetch(root: Path, *, scope: str | None = None) -> int:
    small = fetch_assets(root, large=False)
    if small.get('errors'):
        failure_report(root, 'small_assets', small)
        return 2
    environment = preflight_environment(root, load_lock()['images'])
    if environment['status'] != 'eligible':
        failure_report(root, 'environment', environment)
        return 2
    selected = fetch_assets(root, large=True, only=REFERENCE_ASSETS)
    if selected.get('errors'):
        failure_report(root, 'reference_assets', selected)
        return 2
    tools = PinnedTools(root, load_lock()['images'], environment['cpus'], environment['memory_bytes'])
    reference = prepare_reference_inputs(root, tools, scopes=(scope,) if scope else ('chr20', 'autosomes'))
    if reference.get('reference_inputs_status') != 'eligible':
        failure_report(root, 'reference_eligibility', reference)
        return 2
    assets = fetch_assets(root, large=True)
    if assets.get('errors'):
        failure_report(root, 'full_assets', assets)
        return 2
    print(json.dumps({'status': assets['status'], 'manifest': str(root / 'asset-manifest.json')}))
    return 0


def preflight(root: Path, cpus: int, memory_gib: int, *, scope: str | None = None) -> int:
    assets = verify_assets(root)
    if assets.get('errors'):
        failure_report(root, 'assets', assets)
        return 2
    environment = preflight_environment(root, load_lock()['images'], cpus=cpus, memory_gib=memory_gib)
    if environment['status'] != 'eligible':
        failure_report(root, 'environment', environment)
        return 2
    tools = PinnedTools(root, load_lock()['images'], environment['cpus'], environment['memory_bytes'])
    inputs = prepare_inputs(root, tools, scopes=(scope,) if scope else ('chr20', 'autosomes'))
    if inputs.get('input_status') != 'eligible':
        failure_report(root, 'input_eligibility', inputs)
        return 2
    print(json.dumps({'status': 'eligible', 'manifest': str(root / 'input-manifest.json')}))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('fetch', 'preflight'):
        sub = commands.add_parser(name)
        sub.add_argument('--data-root', type=Path, required=True)
        sub.add_argument('--scope', choices=('chr20', 'autosomes'),
                         help='Audit this declared scope only; omission audits both, without automatic downgrade')
        if name == 'preflight':
            sub.add_argument('--cpus', type=int, default=16)
            sub.add_argument('--memory-gib', type=int, default=64)
    run = commands.add_parser('run')
    run.add_argument('--data-root', type=Path, required=True)
    run.add_argument('--scope', choices=('chr20', 'autosomes'), required=True)
    run.add_argument('--profile', default='region1m-locus-float32')
    run.add_argument('--cpus', type=int, default=16)
    run.add_argument('--memory-gib', type=int, default=64)
    run.add_argument('--output', type=Path, required=True)
    evaluate = commands.add_parser('evaluate')
    evaluate.add_argument('--run', type=Path, required=True)
    evaluate.add_argument('--output', type=Path)
    report = commands.add_parser('report')
    report.add_argument('--run', type=Path, required=True)
    report.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        if args.command in ('fetch', 'preflight'):
            root = args.data_root.resolve()
            root.mkdir(parents=True, exist_ok=True)
            return (fetch(root, scope=args.scope) if args.command == 'fetch' else
                    preflight(root, args.cpus, args.memory_gib, scope=args.scope))
        from giab_reporting import report_run
        if args.command == 'run':
            from giab_pipeline import run_pipeline
            result = run_pipeline(args.data_root, args.output, scope=args.scope, profile=args.profile,
                                  cpus=args.cpus, memory_gib=args.memory_gib)
            report_run(args.output)
            print(json.dumps({**{key: result[key] for key in
                                 ('input_status', 'execution_status', 'contract_status', 'resource_status', 'accuracy_status')},
                              'run': str(args.output.resolve() / 'run.json'),
                              'report': str(args.output.resolve() / 'giab_summary.json')}))
            return 0 if result['execution_status'] == 'completed' else 2
        if args.command == 'evaluate':
            from giab_evaluation import evaluate_run
            result = evaluate_run(args.run, args.output)
            report_run(args.run)
            print(json.dumps({'accuracy_status': result['accuracy_status'], 'output': result['output'],
                              'evaluations': [{key: item.get(key) for key in ('id', 'accuracy_status', 'errors')}
                                              for item in result['evaluations']]}))
            return 0 if result['accuracy_status'] == 'measured' else 2
        source = args.run.resolve()
        if source.is_dir() and not (source / 'run.json').is_file():
            # Preflight reports have no attempted run and remain independently readable.
            content = json.loads((source / 'giab_summary.json').read_text())
            if args.output is not None:
                raise ValueError('--output requires a run manifest, not a preflight report')
            print(json.dumps(content, indent=2, sort_keys=True))
        else:
            content = report_run(source, args.output)
            print(json.dumps({key: content.get(key) for key in
                              ('input_status', 'execution_status', 'contract_status', 'resource_status',
                               'accuracy_status', 'artifacts', 'report_errors')}, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        if args.command in ('fetch', 'preflight'):
            failure_report(args.data_root.resolve(), args.command, {'error': str(error)})
        else:
            print(json.dumps({'command': args.command, 'error': str(error),
                              'policy': 'Previous run and evaluation evidence is not overwritten'}), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
