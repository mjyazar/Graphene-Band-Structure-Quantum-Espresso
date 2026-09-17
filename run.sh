$!/bin/bash

echo "RUNNING PROGRAM"

python main.py

echo "PROGRAM EXECUTION COMPLETE"

git pull --rebase origin main
git add outputs/
git commit -m "Results from Nectar"
git push
