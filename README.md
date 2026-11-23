# FeGenie worker

This worker reads uploaded jobs from S3, runs FeGenie on the worker host, and publishes the report and scientific result files. It requires the existing FeGenie environment, AWS credentials and access to the configured input/results buckets.

## GenBank dependency

GenBank jobs use FeGenie's `--gbk` mode. That mode requires `gtt-genbank-to-AA-seqs` on PATH **inside the environment that runs FeGenie**. Installing GToTree into a different conda environment does not make the helper available in `fegenie`.

On The Belly, install the helper through GToTree:

```bash
/home/ark/miniconda3/bin/conda install -n fegenie -c astrobiomike -c conda-forge -c bioconda gtotree
```

Upstream installation documentation: https://github.com/AstrobioMike/GToTree/wiki/installation

The worker checks this dependency, Biopython and FeGenie's resource files before starting a GenBank analysis. FASTA jobs do not require the GenBank helper. You can run the check alone, without accessing AWS or submitting an analysis:

```bash
/home/ark/miniconda3/envs/fegenie/bin/python fegenie_worker.py \
  --command-prefix "/home/ark/miniconda3/bin/conda run -n fegenie" \
  --check-dependencies
```

## Resource paths and GenBank preparation

The worker runs `fegenie_runtime.py` **inside** the selected conda environment. It finds `hmms/iron/` and `rscripts/` beside the configured `FeGenie.py`, validates the cutoff/map and R-script files, then sets `iron_hmms` and `rscripts` explicitly for the child process. Existing valid environment paths are a fallback for packaged installations. This supports older FeGenie source checkouts that otherwise look for `/HMM-bitcutoffs.txt` when their activation setup is missing.

For a different resource layout, pass `--hmm-dir /path/to/hmms/iron --rscripts-dir /path/to/rscripts` to the worker/wrapper, or set `FEGENIE_HMM_DIR` and `FEGENIE_RSCRIPT_DIR`. Explicit invalid paths fail instead of silently selecting another library. FeGenie runs in the private job directory, isolating its temporary path-discovery files.

Some GenBank exports have CDS coordinates and nucleotide sequences but omit `/translation`, sometimes for every CDS in a genome. The upstream GToTree converter raises `KeyError: 'translation'` for these. Before calling FeGenie, the runtime recovers missing translations from complete annotated CDS in the disposable `bins/` copies. It respects strand, `/codon_start`, and `/transl_table`, using GenBank's standard-code default (table 1) if no table is supplied. Strict CDS validation requires a valid start, a terminal stop, a whole number of codons and no internal stops. No gene prediction is rerun, existing supplied protein translations are preserved, and original `input/` uploads are unchanged.

Missing translations are not recovered for pseudogenes, translation exceptions, partial/uncertain or compound/remote locations, ambiguous/missing nucleotide sequences, or invalid CDS. These features are excluded from the run copy and recorded with a reason in `genbank-exclusions.csv`. A genome with no supplied or recoverable proteins still fails explicitly; its exclusion report is retained locally. The converter's other existing filters continue to apply.

Recovered translations are recorded in `genbank-translations.csv` with genome, contig, locus tag, protein ID, location, genetic code, codon offset and protein length. Both CSV files are included in the downloadable `raw-results.tar.gz`. The local run log gives supplied/recovered/excluded counts per genome. A successfully processed genome can therefore have fewer analyzed CDS than annotated CDS.

Translation rules follow the [Biopython SeqFeature API](https://biopython.org/docs/latest/api/Bio.SeqFeature.html) and [NCBI genetic-code conventions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/).

## Recover Rebecca's failed job

After merging this backend PR, pull the update and retry the selected job. GToTree is already installed after the previous update; no additional installation or FeGenie source changes are needed for this fix. The existing cron entry does not need changing:

```bash
cd /home/ark/MAB/bin/fegenie-local
git pull --ff-only
bash run_fegenie_worker.sh --job rydYvnC9xG --retry-failed
tail -n 40 /home/ark/MAB/bin/fegenie_worker_cron.log
```

`--job` processes only that exact result code and exits. `--retry-failed` requires its current S3 status to be `failed`; it refuses completed/running jobs. The worker downloads the selected job's retained S3 inputs and uses the same result code. It does not process other jobs or require a new upload while those inputs remain in S3.

The previous local attempt is moved into `/home/ark/MAB/fegenie/previous_attempts/` before starting the new attempt. This preserves its input files, logs and partial outputs. A missing S3 manifest stops the retry before moving the local folder. Failed attempts remain available even with `--clean`.

The wrapper and retry share the worker's lock. If another worker is active, a targeted invocation exits nonzero without processing anything; retry it when that worker finishes. The updated lock also respects the PID file of an older worker during rollout and tolerates stale PID files. Do not delete an active worker's lock file.

The log should end with `Completed job rydYvnC9xG`. The published `status.json` then has `"state": "complete"`. A failed targeted retry exits nonzero and records a safe failure message; detailed errors stay in its local log. If inputs have expired from S3, restore the saved local `input/` files and manifest to the original input prefix before retrying.

Rebecca's result page: https://main.d3csq8rzw3uxyt.amplifyapp.com/#fegenie/results-rydYvnC9xG

## Tests

Run in a Python environment with boto3 and Biopython installed (the existing `fegenie` environment is suitable):

```bash
python -m unittest discover -s tests -v
```

Tests use local GenBank/resource fixtures, a subprocess extraction fixture and mocked S3/analysis commands; they do not run HMM analyses or send notifications.


