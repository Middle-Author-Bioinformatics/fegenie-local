"""Run legacy FeGenie in its conda environment with explicit resource paths.

Only the worker's disposable bins/ copies are prepared. Uploaded input/ files
are retained verbatim, and missing protein translations are never inferred.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys


def resource_directory(explicit, inherited, candidate, required, label):
    choices = [Path(explicit)] if explicit else [candidate]
    if not explicit and inherited:
        choices.append(Path(inherited))
    for directory in choices:
        directory = directory.expanduser().resolve()
        if directory.is_dir() and all((directory / name).is_file() for name in required):
            return directory
    raise RuntimeError(
        f"Cannot locate FeGenie's {label}. Check the source installation or set "
        f"--{'hmm-dir' if label == 'HMM library' else 'rscripts-dir'} explicitly."
    )


def execution_environment(fegenie_bin, hmm_dir=None, rscripts_dir=None):
    found = shutil.which(fegenie_bin)
    script = Path(found or fegenie_bin).expanduser().resolve()
    if not script.is_file():
        raise RuntimeError("Cannot locate the configured FeGenie.py script.")
    root = script.parent
    hmms = resource_directory(hmm_dir, os.environ.get('iron_hmms'), root / 'hmms/iron',
                              ['HMM-bitcutoffs.txt', 'FeGenie-map.txt'], 'HMM library')
    rscripts = resource_directory(rscripts_dir, os.environ.get('rscripts'), root / 'rscripts',
                                  ['dendro-heatmap.R'], 'R scripts')
    env = os.environ.copy()
    # Set these after conda activation, including for older FeGenie versions.
    env['iron_hmms'] = str(hmms)
    env['rscripts'] = str(rscripts)
    return script, env


def check_genbank_dependencies():
    if not shutil.which('gtt-genbank-to-AA-seqs'):
        raise RuntimeError('Missing gtt-genbank-to-AA-seqs in the FeGenie execution environment.')
    try:
        from Bio import SeqIO  # noqa: F401
    except ImportError as error:
        raise RuntimeError('Missing Biopython in the FeGenie execution environment.') from error


def prepare_genbank(bins_dir, extension, output_dir):
    """Exclude CDS without /translation from run copies; document every exclusion.

    Other features and supplied translations are preserved. Upstream GToTree's
    remaining protein filters (partial calls, pseudogenes, etc.) still apply.
    """
    from Bio import SeqIO

    fields = ['genome', 'contig', 'locus_tag', 'protein_id', 'location', 'reason']
    exclusions = []
    files = sorted(Path(bins_dir).glob(f'*.{extension}'))
    if not files:
        raise RuntimeError('No GenBank files found for the configured extension.')
    for path in files:
        with path.open() as handle:
            records = list(SeqIO.parse(handle, 'genbank'))
        if not records:
            raise RuntimeError(f'No GenBank records found in {path.name}.')
        missing = 0
        translated = 0
        for record in records:
            kept = []
            for feature in record.features:
                translation = feature.qualifiers.get('translation', [])
                if feature.type == 'CDS' and not (translation and str(translation[0]).strip()):
                    exclusions.append({
                        'genome': path.name, 'contig': record.name,
                        'locus_tag': feature.qualifiers.get('locus_tag', [''])[0],
                        'protein_id': feature.qualifiers.get('protein_id', [''])[0],
                        'location': str(feature.location), 'reason': 'missing_translation',
                    })
                    missing += 1
                    continue
                if feature.type == 'CDS':
                    translated += 1
                kept.append(feature)
            record.features = kept
        print(f'{path.name}: {translated} CDS with supplied translations; '
              f'{missing} CDS excluded for missing translation.', flush=True)
        if not translated:
            raise RuntimeError(f'{path.name} has no CDS with supplied protein translations.')
        if missing:
            temporary = path.with_name(path.name + '.tmp')
            try:
                SeqIO.write(records, temporary, 'genbank')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / 'genbank-exclusions.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(exclusions)
    return exclusions


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fegenie-bin', required=True)
    parser.add_argument('--hmm-dir')
    parser.add_argument('--rscripts-dir')
    parser.add_argument('--check-dependencies', action='store_true')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    try:
        script, env = execution_environment(args.fegenie_bin, args.hmm_dir, args.rscripts_dir)
        if args.check_dependencies or '--gbk' in command:
            try:
                check_genbank_dependencies()
            except RuntimeError as error:
                print(str(error), file=sys.stderr)
                return 3
        if args.check_dependencies:
            print('FeGenie resources and GenBank dependencies are available.', flush=True)
            return 0
        inputs = argparse.ArgumentParser(add_help=False)
        inputs.add_argument('-bin_dir', required=True)
        inputs.add_argument('-bin_ext', required=True)
        inputs.add_argument('-out', required=True)
        paths, _ = inputs.parse_known_args(command)
        # The child runs in the job directory, so resolve worker-relative paths
        # before changing cwd. Replace only these named CLI arguments.
        paths.bin_dir = str(Path(paths.bin_dir).resolve())
        paths.out = str(Path(paths.out).resolve())
        for option, value in [('-bin_dir', paths.bin_dir), ('-out', paths.out)]:
            command[command.index(option) + 1] = value
        if '--gbk' in command:
            prepare_genbank(paths.bin_dir, paths.bin_ext, paths.out)
        # Legacy FeGenie creates temporary HMMlib.txt/rscripts.txt in its cwd.
        # Keep those in this job, rather than the caller's directory.
        return subprocess.run([sys.executable, str(script), *command], env=env,
                              cwd=Path(paths.bin_dir).resolve().parent).returncode
    except (RuntimeError, ValueError, OSError) as error:
        print(f'FeGenie preparation failed: {error}', file=sys.stderr, flush=True)
        return 4


if __name__ == '__main__':
    sys.exit(main())
