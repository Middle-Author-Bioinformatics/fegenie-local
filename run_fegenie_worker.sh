#!/usr/bin/env bash
set -u

export PATH="/home/ark/miniconda3/bin:/home/ark/miniconda3/envs/fegenie/bin:/usr/local/bin:/usr/bin:/bin"
export AWS_REGION="us-east-2"

LOG="/home/ark/MAB/bin/fegenie_worker_cron.log"

echo "===== $(date) starting FeGenie worker =====" >> "$LOG"

# The worker uses an OS lock; do not remove its lock file.

/home/ark/miniconda3/envs/fegenie/bin/python \
  /home/ark/MAB/bin/fegenie-local/fegenie_worker.py \
  --input-bucket midauthorbio-fegenie-input \
  --results-bucket midauthorbio-fegenie-results \
  --region us-east-2 \
  --work-root /home/ark/MAB/fegenie \
  --fegenie-bin /home/ark/MAB/bin/FeGenie/FeGenie.py \
  --command-prefix "/home/ark/miniconda3/bin/conda run -n fegenie" \
  --once \
  --continue-on-error \
  "$@" \
  >> "$LOG" 2>&1

status=$?
echo "===== $(date) finished FeGenie worker (exit $status) =====" >> "$LOG"
exit "$status"