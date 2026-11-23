"""Run legacy FeGenie in its conda environment with explicit resource paths.

Only the worker's disposable bins/ copies are prepared. Uploaded input/ files
are retained verbatim. Missing translations can be recovered from validated,
annotated CDS; gene prediction is never rerun.
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


def qualifier(feature, name, default=''):
    values = feature.qualifiers.get(name, [])
    return values[0] if values else default


def recover_translation(record, feature):
    """Translate a complete, ordinary annotated CDS, or explain its exclusion."""
    from Bio.Data.CodonTable import TranslationError
    from Bio.SeqFeature import ExactPosition

    if 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers:
        return None, 'pseudogene'
    if 'transl_except' in feature.qualifiers or 'exception' in feature.qualifiers:
        return None, 'translation_exception'
    note = str(qualifier(feature, 'note')).lower()
    if any(term in note for term in ['frameshifted', 'internal stop', 'incomplete']):
        return None, 'noted_incomplete_or_disrupted'
    location = feature.location
    if location is None:
        return None, 'missing_location'
    # The existing GToTree converter excludes joins and fuzzy/partial calls.
    if len(location.parts) != 1 or location.ref is not None:
        return None, 'compound_or_remote_location'
    if not isinstance(location.start, ExactPosition) or not isinstance(location.end, ExactPosition):
        return None, 'partial_or_uncertain_location'
    if location.strand not in (1, -1):
        return None, 'unknown_strand'
    if not 0 <= int(location.start) < int(location.end) <= len(record.seq):
        return None, 'location_outside_sequence'
    try:
        nucleotides = str(feature.extract(record.seq)).upper()
        if not nucleotides or set(nucleotides) - set('ACGT'):
            return None, 'missing_or_ambiguous_nucleotides'
        # SeqFeature.translate respects strand, codon_start and transl_table.
        # GenBank's documented default when transl_table is absent is table 1.
        protein = str(feature.translate(record.seq, table=1, cds=True))
    except (ValueError, TranslationError, TypeError, KeyError) as error:
        return None, 'invalid_cds: ' + str(error)
    if not protein:
        return None, 'empty_translation'
    return protein, ''


def prepare_genbank(bins_dir, extension, output_dir):
    """Recover valid missing CDS translations in run copies and report changes.

    Other features and supplied translations are preserved. Upstream GToTree's
    remaining protein filters (partial calls, pseudogenes, etc.) still apply.
    """
    from Bio import SeqIO

    fields = ['genome', 'contig', 'locus_tag', 'protein_id', 'location', 'reason']
    exclusions = []
    recoveries = []
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def write_reports():
        for name, columns, rows in [
            ('genbank-exclusions.csv', fields, exclusions),
            ('genbank-translations.csv', fields[:-1] + [
                'translation_table', 'codon_start', 'protein_length', 'source'], recoveries),
        ]:
            with (output_dir / name).open('w', newline='', encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                writer.writerows(rows)
    files = sorted(Path(bins_dir).glob(f'*.{extension}'))
    if not files:
        raise RuntimeError('No GenBank files found for the configured extension.')
    for path in files:
        with path.open() as handle:
            records = list(SeqIO.parse(handle, 'genbank'))
        if not records:
            raise RuntimeError(f'No GenBank records found in {path.name}.')
        excluded = 0
        supplied = 0
        recovered = 0
        for record in records:
            kept = []
            for feature in record.features:
                translation = feature.qualifiers.get('translation', [])
                if feature.type == 'CDS' and not (translation and str(translation[0]).strip()):
                    row = {
                        'genome': path.name, 'contig': record.name,
                        'locus_tag': qualifier(feature, 'locus_tag'),
                        'protein_id': qualifier(feature, 'protein_id'),
                        'location': str(feature.location),
                    }
                    protein, reason = recover_translation(record, feature)
                    if protein is None:
                        exclusions.append(dict(row, reason=reason))
                        excluded += 1
                        continue
                    feature.qualifiers['translation'] = [protein]
                    recoveries.append(dict(row,
                        translation_table=qualifier(feature, 'transl_table', '1'),
                        codon_start=qualifier(feature, 'codon_start', '1'),
                        protein_length=len(protein), source='annotated_CDS_translation'))
                    recovered += 1
                elif feature.type == 'CDS':
                    supplied += 1
                kept.append(feature)
            record.features = kept
        print(f'{path.name}: {supplied} supplied CDS translations; '
              f'{recovered} recovered from annotated DNA; {excluded} CDS excluded.', flush=True)
        # Keep QC reports available even when a genome has no usable proteins.
        write_reports()
        if not supplied + recovered:
            raise RuntimeError(f'{path.name} has no usable CDS translations. '
                               'See genbank-exclusions.csv for reasons.')
        if excluded or recovered:
            temporary = path.with_name(path.name + '.tmp')
            try:
                SeqIO.write(records, temporary, 'genbank')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
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

