# Worker brief: survey TESS sectors for the ad-astra search

You are one of several cloud sessions running the same, already-validated survey pipeline on different
TESS sectors. A coordinating session wrote and tested the pipeline and assigned you a list of sectors.
Your job is operational: run it, keep it running, and report. Do not edit pipeline code. If something
breaks, record what happened (see Failures) and continue with the next sector.

## What it does (context)
`pipeline/worker.sh S1 S2 ...` sets up a Python venv, builds a local extract of the TESS Input Catalog
from the public MAST S3 bucket (in the background), then for each sector runs `pipeline/run_sector.py`:
it reads every QLP full-frame light curve of that sector for the ~2 million stars in
`results/samples/survey_sample.parquet` straight from S3, flattens them, catalogues every significant
dip, flags spacecraft systematics/blends/artefacts, and writes small candidate tables for three science
tracks to `results/survey/sNNNN/`, which it commits and pushes to branch
`claude/gifted-albattani-szay6x` (rebasing on other workers' pushes). Only S3 access is needed.

## Steps
1. Make sure you are on the branch and up to date:
   `git fetch origin claude/gifted-albattani-szay6x && git checkout claude/gifted-albattani-szay6x && git pull --rebase origin claude/gifted-albattani-szay6x`
2. Launch the worker in the background (it runs for hours):
   `nohup bash pipeline/worker.sh <your sectors, space separated> > /home/user/worker_stdout.log 2>&1 &`
3. Watch progress with a Monitor on `/home/user/work/worker.log` (each line is an event: start/done/pushed/FAILED
   per sector, then `ALL DONE`). Never use foreground sleep loops. Re-arm the monitor when it expires.
   Typical sector: 20-60 minutes.
4. If the worker process dies before `ALL DONE` (check with `pgrep -f worker.sh`), relaunch the same
   command; finished sectors are skipped automatically.

## Rules
- Only files under `results/survey/` are committed (the script does this). Never edit, revert or delete
  other files, never force-push, never open a pull request, never change branches.
- If a push keeps failing, leave the results committed locally and say so in your final report.

## Failures
The script retries a failed sector once, then commits `results/survey/sNNNN/FAILED.txt` with the log tail.
If many sectors fail the same way, stop and report the error instead of burning through the list.

## Final report
When `ALL DONE` appears, reply with a table: sector, light curves, n_trackb, n_trackc, n_trackd, runtime
(from each `results/survey/sNNNN/summary.json`), and any failures.
