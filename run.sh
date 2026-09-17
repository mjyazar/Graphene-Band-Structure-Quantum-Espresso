#!/bin/bash

set -o pipefail  # if any command in pipeline fails, get a failure status
set -e  # enable the "errexit" option

echo "PULLING LATEST CODE"
git pull --rebase origin main

mkdir -p log  # don’t complain if folder already exists

echo "RUNNING PROGRAM"
python main.py 2>&1 | tee log/run.log  # stderror goes where stdout goes, print in both terminal and save in run.log
echo "PROGRAM EXECUTION COMPLETE"

echo "PUSHING PROGRAM"
git pull --rebase origin main
git add outputs/
git commit -m "Results from Nectar"
git push
