#!/usr/bin/env bash
set -u

export PATH="/home/ark/miniconda3/bin:/home/ark/miniconda3/envs/fegenie/bin:/usr/local/bin:/usr/bin:/bin"
export AWS_REGION="us-east-2"

LOG="/home/ark/MAB/bin/fegenie_worker_cron.log"

echo "===== $(date) starting FeGenie worker =====" >> "$LOG"

# Remove stale lock if no matching worker is running
if [ -f /tmp/fegenie-worker.lock ]; then
    if ! pgrep -f "fegenie_worker.py" >/dev/null 2>&1; then
        echo "$(date) removing stale lock file" >> "$LOG"
        rm -f /tmp/fegenie-worker.lock
    fi
fi

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
  >> "$LOG" 2>&1

echo "===== $(date) finished FeGenie worker =====" >> "$LOG"