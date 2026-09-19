#!/bin/bash

set -o pipefail  # if any command in pipeline fails, get a failure status
set -e  # enable the "errexit" option

echo -e "PULLING LATEST CODE\n"
git restore outputs/  # throw away any uncommitted output from a previous run
git pull --rebase origin main  # get the latest repository

mkdir -p log  # don’t complain if folder already exists

echo -e "\nRUNNING PROGRAM"
python main.py 2>&1 | tee log/run.log  # stderror goes where stdout goes, print in both terminal and save in run.log
echo "PROGRAM EXECUTION COMPLETE"

echo "PUSHING PROGRAM"
git add outputs/
git add log/
git commit -m "Results from Nectar"
git pull --rebase
git push
