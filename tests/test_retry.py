import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from botocore.exceptions import ClientError
import fegenie_worker as worker


class S3:
    def __init__(self, status=None):
        self.status = status
        self.writes = []
        self.uploads = []
    def head_object(self, **kwargs):
        if self.status is None:
            raise ClientError({'Error': {'Code': '404'}}, 'HeadObject')
        return {}
    def get_object(self, **kwargs):
        if self.status is None:
            raise ClientError({'Error': {'Code': 'NoSuchKey'}}, 'GetObject')
        return {'Body': io.BytesIO(json.dumps(self.status).encode())}
    def put_object(self, **kwargs):
        self.status = json.loads(kwargs['Body'])
        self.writes.append(self.status.copy())
    def upload_file(self, source, bucket, key, **kwargs):
        self.uploads.append(key)


class RetryTests(unittest.TestCase):
    def args(self, root, *extra):
        return worker.parse_args(['--work-root', str(root), '--ses-sender', '', *extra])

    def download(self, s3, bucket, prefix, dest):
        dest.mkdir(parents=True)
        slug = prefix.removeprefix('fegenie-').rstrip('/')
        (dest / f'manifest-{slug}.txt').write_text(f'job_slug: {slug}\n')
        (dest / 'sample.gb').write_text('LOCUS example\n')

    def run_success(self, cmd, log):
        if '-out' not in cmd:
            Path(cmd[cmd.index('-o') + 1]).write_text('<html>scientific report</html>')
            return
        out = Path(cmd[cmd.index('-out') + 1]);out.mkdir()
        (out / 'FeGenie-heatmap-data.csv').write_text('X,sample\niron,0\n')
        (out / 'FeGenie-geneSummary.csv').write_text('category,genome/assembly,orf,HMM\n')

    def test_retry_and_force_require_a_selected_safe_job(self):
        for flags in [['--retry-failed'], ['--force'], ['--job', '../bad'], ['--job', 'A1', '--retry-failed', '--force']]:
            with self.subTest(flags=flags), self.assertRaises(SystemExit):worker.parse_args(flags)

    def test_job_selection_does_not_scan_queue(self):
        args = self.args('.', '--job', 'Example123', '--retry-failed')
        s3 = Mock()
        self.assertTrue(args.once)
        self.assertEqual(worker.selected_prefixes(args, s3), ['fegenie-Example123/'])
        s3.get_paginator.assert_not_called()

    def test_retry_refuses_unfailed_or_missing_status_before_mutations(self):
        for state in ['complete', 'running', 'queued', None]:
            with self.subTest(state=state), tempfile.TemporaryDirectory() as tmp:
                args = self.args(tmp, '--job', 'Example123', '--retry-failed')
                s3 = S3({'state': state} if state else None)
                with self.assertRaises(RuntimeError):worker.process_job(args, s3, 'fegenie-Example123/')
                self.assertEqual(s3.writes, [])

    def test_unrelated_job_and_completed_jobs_are_not_reprocessed(self):
        args = self.args('.', '--job', 'Example123', '--retry-failed')
        s3 = S3({'state':'failed'})
        worker.process_job(args, s3, 'fegenie-Other123/')
        self.assertEqual(s3.writes, [])
        worker.process_job(self.args('.'), S3({'state':'complete'}), 'fegenie-Example123/')

    def test_successful_retry_preserves_old_inputs_and_publishes_complete_last(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.args(tmp, '--job', 'Example123', '--retry-failed')
            previous = Path(tmp) / 'fegenie-Example123';previous.mkdir()
            (previous / 'run.log').write_text('original failure')
            (previous / 'input').mkdir();(previous / 'input/original.gb').write_text('original genome')
            s3 = S3({'state':'failed'})
            with patch.object(worker, 'find_manifest_key', return_value='fegenie-Example123/manifest-Example123.txt'), patch.object(worker, 'download_prefix', side_effect=self.download), patch.object(worker, 'check_genbank_helper') as check, patch.object(worker, 'run_command', side_effect=self.run_success):
                worker.process_job(args, s3, 'fegenie-Example123/')
            check.assert_called_once()
            archived = list((Path(tmp)/'previous_attempts').iterdir())
            self.assertEqual(len(archived),1)
            self.assertEqual((archived[0]/'run.log').read_text(),'original failure')
            self.assertEqual((archived[0]/'input/original.gb').read_text(),'original genome')
            self.assertEqual([s['state'] for s in s3.writes],['running','complete'])
            self.assertTrue(any(k.endswith('fegenie-report.html') for k in s3.uploads))

    def test_notification_failure_does_not_change_completed_analysis_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.args(tmp, '--job', 'Example123', '--retry-failed', '--ses-sender', 'sender@example.org')
            s3 = S3({'state':'failed'})
            with patch.object(worker, 'find_manifest_key', return_value='fegenie-Example123/manifest-Example123.txt'), patch.object(worker, 'download_prefix', side_effect=self.download), patch.object(worker, 'parse_manifest', return_value=worker.JobManifest('Example123', submitter_email='recipient@example.org')), patch.object(worker, 'check_genbank_helper'), patch.object(worker, 'run_command', side_effect=self.run_success), patch.object(worker, 'send_ses_email', side_effect=RuntimeError('simulated notification failure')):
                worker.process_job(args, s3, 'fegenie-Example123/')
            self.assertEqual(s3.status['state'],'complete')

    def test_missing_helper_fails_before_analysis_keeps_debug_files_and_publishes_safe_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.args(tmp, '--job', 'Example123', '--retry-failed', '--clean', '--continue-on-error')
            s3 = S3({'state':'failed'})
            with patch.object(worker, 'find_manifest_key', return_value='fegenie-Example123/manifest-Example123.txt'), patch.object(worker, 'download_prefix', side_effect=self.download), patch.object(worker, 'check_genbank_helper', side_effect=worker.MissingGenBankHelper('/private/path/helper missing')), patch.object(worker, 'run_command') as run:
                with self.assertRaises(worker.MissingGenBankHelper):worker.process_job(args,s3,'fegenie-Example123/')
            run.assert_not_called()
            self.assertEqual(s3.status['error_code'],'missing_genbank_helper')
            self.assertNotIn('/private',json.dumps(s3.status))
            self.assertTrue((Path(tmp)/'fegenie-Example123/run.log').exists())
            self.assertEqual(s3.uploads,[])

    def test_probe_runs_in_selected_environment_and_preserves_quoted_arguments(self):
        args = self.args('.', '--command-prefix', '"/path with space/conda" run -n fegenie')
        with patch.object(worker.subprocess, 'run', return_value=Mock(returncode=0)) as run:
            worker.check_genbank_helper(args)
        self.assertEqual(run.call_args.args[0][:4],['/path with space/conda','run','-n','fegenie'])
        with patch.object(worker.subprocess,'run',return_value=Mock(returncode=3)):
            with self.assertRaises(worker.MissingGenBankHelper):worker.check_genbank_helper(args)

    def test_dependency_check_needs_no_aws_or_lock(self):
        args = self.args('.', '--check-dependencies')
        with patch.object(worker, 'parse_args', return_value=args), patch.object(worker,'configure_logging'), patch.object(worker,'check_genbank_helper') as check, patch.object(worker,'make_s3') as aws, patch.object(worker,'acquire_lock') as lock:
            self.assertEqual(worker.main(),0)
        check.assert_called_once();aws.assert_not_called();lock.assert_not_called()

    def test_os_lock_accepts_stale_pid_excludes_second_process_and_releases(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'worker.lock';path.write_text('99999999')
            handle=worker.acquire_lock(path)
            try:
                code='from pathlib import Path; import fegenie_worker as w; w.acquire_lock(Path('+repr(str(path))+'))'
                result=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
            finally:worker.release_lock(handle)
            worker.release_lock(worker.acquire_lock(path))
            self.assertTrue(path.exists())

    def test_no_manifest_does_not_move_local_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'fegenie-Example123';path.mkdir();(path/'run.log').write_text('keep')
            with patch.object(worker,'find_manifest_key',return_value=None):
                with self.assertRaises(RuntimeError):worker.process_job(self.args(tmp,'--job','Example123','--retry-failed'),S3({'state':'failed'}),'fegenie-Example123/')
            self.assertEqual((path/'run.log').read_text(),'keep')


if __name__ == '__main__':unittest.main()
