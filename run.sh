#!/bin/bash

set -o pipefail  # if any command in pipeline fails, get a failure status
set -e  # stop script when a command fails i.e.s the "errexit" option

MAX_LOGS=10
STAMP=$(date +%Y-%m-%d_%H-%M-%S)
LOG_TMP="log/run_${STAMP}.log"

echo -e "PULLING LATEST CODE\n"
git restore outputs/  # discard any uncommitted output from a previous run
git pull --rebase origin main  # get the latest repository

mkdir -p log/complete log/failed  # don’t complain if folder already exists

echo -e "\nRUNNING PROGRAM"

if python -u main.py 2>&1 | tee "$LOG_TMP"; then
    STATUS="complete"
    echo "PROGRAM EXECUTION COMPLETE"
else
    STATUS="failed"
    echo "PROGRAM EXECUTION FAILED"
fi

mv "$LOG_TMP" "log/${STATUS}/"

# Keep only the newest MAX_LOGS logs in this folder.
# The timestamps sort alphabetically = chronologically, so drop all but the last MAX_LOGS
find "log/${STATUS}" -maxdepth 1 -name 'run_*.log' | sort | head -n -"$MAX_LOGS" | xargs -r rm --


echo "PUSHING PROGRAM"
git add outputs/ log/
git diff --cached --quiet || git commit -m "Results from Nectar (${STATUS}, ${STAMP})"
git pull --rebase
git push
