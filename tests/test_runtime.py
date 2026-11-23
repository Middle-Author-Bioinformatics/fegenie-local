import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import fegenie_runtime as runtime
import fegenie_worker as worker


def write_genbank(path, translated=True):
    record = SeqRecord(Seq('ATG' * 30), id='contig.1', name='contig')
    record.annotations['molecule_type'] = 'DNA'
    record.features = [
        SeqFeature(SimpleLocation(0, 9), type='gene', qualifiers={'locus_tag': ['good']}),
        SeqFeature(SimpleLocation(0, 9), type='CDS', qualifiers={
            'locus_tag': ['good'], **({'translation': ['MMA']} if translated else {})}),
        SeqFeature(SimpleLocation(15, 24, strand=-1), type='gene', qualifiers={'locus_tag': ['missing']}),
        SeqFeature(SimpleLocation(15, 24, strand=-1), type='CDS', qualifiers={'locus_tag': ['missing']}),
        SeqFeature(SimpleLocation(30, 39), type='CDS', qualifiers={'locus_tag': ['pseudo'], 'pseudo': ['']}),
        SeqFeature(SimpleLocation(45, 54), type='tRNA', qualifiers={'product': ['tRNA-Met']}),
    ]
    SeqIO.write([record], path, 'genbank')


def source_tree(root):
    source = root / 'FeGenie'; source.mkdir()
    hmms = source / 'hmms/iron'; hmms.mkdir(parents=True)
    (hmms / 'HMM-bitcutoffs.txt').write_text('test cutoff')
    (hmms / 'FeGenie-map.txt').write_text('test map')
    rscripts = source / 'rscripts'; rscripts.mkdir()
    (rscripts / 'dendro-heatmap.R').write_text('# fixture')
    script = source / 'FeGenie.py'; script.write_text('# fixture')
    return script, hmms, rscripts


class RuntimeTests(unittest.TestCase):
    def test_missing_translations_excluded_only_from_run_copy_and_documented(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root / 'input'; inputs.mkdir()
            original = inputs / 'sample.gb'; write_genbank(original)
            original_bytes = original.read_bytes()
            extension, paths = worker.prepare_bins(inputs, root / 'bins')
            exclusions = runtime.prepare_genbank(root / 'bins', extension, root / 'out')
            self.assertEqual(original.read_bytes(), original_bytes)
            record = SeqIO.read(paths[0], 'genbank')
            cds = [feature for feature in record.features if feature.type == 'CDS']
            self.assertEqual(len(cds), 1)
            self.assertEqual(cds[0].qualifiers['translation'], ['MMA'])
            self.assertEqual(len([f for f in record.features if f.type == 'gene']), 2)
            self.assertTrue(any(f.type == 'tRNA' for f in record.features))
            self.assertEqual([row['locus_tag'] for row in exclusions], ['missing', 'pseudo'])
            with (root / 'out/genbank-exclusions.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertTrue(rows[0]['reason'].startswith('invalid_cds:'))
            self.assertEqual(rows[1]['reason'], 'pseudogene')
            self.assertNotIn(str(root), json.dumps(rows))

    def test_no_translations_fails_without_inventing_proteins_or_rewriting_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); bins = root / 'bins'; bins.mkdir()
            genome = bins / 'sample.gbk'; write_genbank(genome, translated=False)
            original = genome.read_bytes()
            with self.assertRaisesRegex(RuntimeError, 'no usable CDS translations'):
                runtime.prepare_genbank(bins, 'gbk', root / 'out')
            self.assertEqual(genome.read_bytes(), original)

    def test_malformed_or_empty_genbank_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root / 'empty.gbk').write_text('')
            with self.assertRaisesRegex(RuntimeError, 'No GenBank records'):
                runtime.prepare_genbank(root, 'gbk', root / 'out')

    def test_resource_paths_found_beside_script_over_stale_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            script, hmms, rscripts = source_tree(Path(tmp))
            with patch.dict(os.environ, {'iron_hmms': '/stale/hmms', 'rscripts': '/stale/rscripts'}):
                resolved, env = runtime.execution_environment(str(script))
            self.assertEqual(resolved, script)
            self.assertEqual(env['iron_hmms'], str(hmms))
            self.assertEqual(env['rscripts'], str(rscripts))

    def test_explicit_missing_resources_fail_instead_of_using_other_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            script, _, _ = source_tree(Path(tmp))
            with self.assertRaisesRegex(RuntimeError, '--hmm-dir'):
                runtime.execution_environment(str(script), str(Path(tmp) / 'absent'))

    def test_resource_overrides_and_installed_environment_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); script, hmms, rscripts = source_tree(root)
            installed = root / 'installed.py'; installed.write_text('# fixture')
            with patch.dict(os.environ, {'iron_hmms': str(hmms), 'rscripts': str(rscripts)}):
                _, env = runtime.execution_environment(str(installed))
            self.assertEqual(env['iron_hmms'], str(hmms))
            _, env = runtime.execution_environment(str(installed), str(hmms), str(rscripts))
            self.assertEqual(env['rscripts'], str(rscripts))

    def test_launcher_prepares_genbank_and_sets_resources_in_child_process(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); script, hmms, rscripts = source_tree(root)
            job = root / 'job'; bins = job / 'bins'; bins.mkdir(parents=True)
            write_genbank(bins / 'sample.gbk')
            helper = root / 'gtt-genbank-to-AA-seqs'
            # Reproduce the upstream helper's direct translation lookup.
            helper.write_text('#!' + sys.executable + '\n' + '''
import sys
from Bio import SeqIO
from pathlib import Path
args = sys.argv
proteins = [f.qualifiers['translation'][0] for rec in SeqIO.parse(args[args.index('-i')+1], 'genbank') for f in rec.features if f.type == 'CDS']
Path(args[args.index('-o')+1]).write_text('\\n'.join(proteins))
''')
            helper.chmod(0o755)
            script.write_text('''
import json, os, subprocess, sys
from pathlib import Path
args = sys.argv
out = Path(args[args.index('-out')+1])
bins = Path(args[args.index('-bin_dir')+1])
assert (Path(os.environ['iron_hmms'])/'HMM-bitcutoffs.txt').read_text() == 'test cutoff'
assert (Path(os.environ['rscripts'])/'dendro-heatmap.R').exists()
subprocess.run(['gtt-genbank-to-AA-seqs', '-i', str(bins/'sample.gbk'), '-o', str(out/'proteins.faa')], check=True)
(out/'observed.json').write_text(json.dumps({'cwd': os.getcwd(), 'python': sys.executable}))
''')
            env = dict(os.environ, PATH=str(root) + os.pathsep + os.environ['PATH'])
            result = subprocess.run([sys.executable, str(Path(runtime.__file__)),
                '--fegenie-bin', str(script), '--', '-bin_dir', str(bins), '-bin_ext', 'gbk',
                '-out', str(job / 'out'), '--gbk'], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((job / 'out/proteins.faa').read_text(), 'MMA')
            self.assertEqual(json.loads((job / 'out/observed.json').read_text())['cwd'], str(job))
            self.assertTrue((job / 'out/genbank-exclusions.csv').exists())

    def test_fasta_does_not_need_genbank_helper_or_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); script, _, _ = source_tree(root)
            script.write_text('import sys; sys.exit(7)')
            bins = root / 'job/bins'; bins.mkdir(parents=True)
            with patch.object(runtime, 'check_genbank_dependencies') as check, patch.object(runtime, 'prepare_genbank') as prepare:
                result = runtime.main(['--fegenie-bin', str(script), '--', '-bin_dir', str(bins),
                                       '-bin_ext', 'fa', '-out', str(root / 'out')])
            self.assertEqual(result, 7)
            check.assert_not_called(); prepare.assert_not_called()

    def test_worker_places_runtime_after_conda_prefix_and_preserves_overrides(self):
        args = worker.parse_args(['--command-prefix', '"/path with space/conda" run -n fegenie',
                                  '--hmm-dir', '/data/iron', '--rscripts-dir', '/data/rscripts'])
        cmd = worker.build_fegenie_command(args, worker.JobManifest('Example123'), Path('/job/bins'), Path('/job/out'), 'gbk')
        self.assertEqual(cmd[:5], ['/path with space/conda', 'run', '-n', 'fegenie', 'python'])
        self.assertIn('fegenie_runtime.py', cmd[5])
        self.assertEqual(cmd[cmd.index('--hmm-dir')+1], '/data/iron')
        self.assertEqual(cmd[cmd.index('--rscripts-dir')+1], '/data/rscripts')
        self.assertIn('--gbk', cmd[cmd.index('--')+1:])

    def test_missing_cutoffs_stops_before_preparing_inputs_or_running_fegenie(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); script, hmms, _ = source_tree(root)
            (hmms / 'HMM-bitcutoffs.txt').unlink()
            with patch.dict(os.environ, {'iron_hmms': ''}), patch.object(runtime, 'prepare_genbank') as prepare, patch.object(runtime.subprocess, 'run') as run:
                result = runtime.main(['--fegenie-bin', str(script), '--', '-bin_dir', str(root),
                                       '-bin_ext', 'gbk', '-out', str(root / 'out'), '--gbk'])
            self.assertEqual(result, 4)
            prepare.assert_not_called(); run.assert_not_called()


if __name__ == '__main__':
    unittest.main()

