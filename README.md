# FeGenie worker

This worker reads uploaded jobs from S3, runs FeGenie on the worker host, and publishes the report and scientific result files. It requires the existing FeGenie environment, AWS credentials and access to the configured input/results buckets.

## GenBank dependency

GenBank jobs use FeGenie's `--gbk` mode. That mode requires `gtt-genbank-to-AA-seqs` on PATH **inside the environment that runs FeGenie**. Installing GToTree into a different conda environment does not make the helper available in `fegenie`.

On The Belly, install the helper through GToTree:

```bash
/home/ark/miniconda3/bin/conda install -n fegenie -c astrobiomike -c conda-forge -c bioconda gtotree
```

Upstream installation documentation: https://github.com/AstrobioMike/GToTree/wiki/installation

The worker checks this dependency before starting a GenBank analysis. FASTA jobs do not require the GenBank helper. You can run the check alone, without accessing AWS or submitting an analysis:

```bash
/home/ark/miniconda3/envs/fegenie/bin/python fegenie_worker.py \
  --command-prefix "/home/ark/miniconda3/bin/conda run -n fegenie" \
  --check-dependencies
```

## Recover Rebecca's failed job

After merging this backend PR, pull the update and install the dependency. The existing cron entry does not need changing. The updated wrapper accepts worker arguments:

```bash
cd /home/ark/MAB/bin/fegenie-local
git pull --ff-only
/home/ark/miniconda3/bin/conda install -n fegenie -c astrobiomike -c conda-forge -c bioconda gtotree
bash run_fegenie_worker.sh --job rydYvnC9xG --retry-failed
tail -n 40 /home/ark/MAB/bin/fegenie_worker_cron.log
```

`--job` processes only that exact result code and exits. `--retry-failed` requires its current S3 status to be `failed`; it refuses completed/running jobs. The worker downloads the selected job's retained S3 inputs and uses the same result code. It does not process other jobs or require a new upload while those inputs remain in S3.

The previous local attempt is moved into `/home/ark/MAB/fegenie/previous_attempts/` before starting the new attempt. This preserves its input files, logs and partial outputs. A missing S3 manifest stops the retry before moving the local folder. Failed attempts remain available even with `--clean`.

The wrapper and retry share the worker's lock. If another worker is active, a targeted invocation exits nonzero without processing anything; retry it when that worker finishes. The updated lock also respects the PID file of an older worker during rollout and tolerates stale PID files. Do not delete an active worker's lock file.

The log should end with `Completed job rydYvnC9xG`. The published `status.json` then has `"state": "complete"`. A failed targeted retry exits nonzero and records a safe failure message; detailed errors stay in its local log. If inputs have expired from S3, restore the saved local `input/` files and manifest to the original input prefix before retrying.

Rebecca's result page: https://main.d3csq8rzw3uxyt.amplifyapp.com/#fegenie/results-rydYvnC9xG

## Tests

Run in a Python environment with boto3 installed (the existing `fegenie` environment is suitable):

```bash
python -m unittest discover -s tests -v
```

Tests use local fixtures and mocked S3/analysis commands; they do not run analyses or send notifications.
