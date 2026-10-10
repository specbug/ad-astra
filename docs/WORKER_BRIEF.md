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
3. Wait cheaply. Do NOT use Monitor or repeated status checks: every event costs a full model turn and all
   worker sessions share one usage limit (it was exhausted once by monitor traffic). Instead start ONE
   background Bash command (`run_in_background: true`, `timeout: 7200000`) that exits when the job ends:
   `until grep -q "ALL DONE" /home/user/work/worker.log || ! pgrep -f "pipeline/[w]orker.sh" >/dev/null; do sleep 120; done; tail -3 /home/user/work/worker.log`
   (the `[w]` keeps pgrep from matching this command itself; move an old `worker.log` aside before launching,
   so a stale `ALL DONE` cannot end the wait).
   Then end your turn and do nothing until it returns. Typical sector: 20-60 minutes.
4. When the wait returns: if the log shows `ALL DONE`, write the final report. If the worker died without
   `ALL DONE`, relaunch step 2 (finished sectors are skipped automatically) and wait again. If the wait timed
   out while the worker is still running, start the same wait again.

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

## Track A workers
Same rules, different script: `nohup bash pipeline/worker_a.sh START STOP > /home/user/worker_a_stdout.log 2>&1 &`
searches rows [START, STOP) of `results/track_a/targets_search.csv` (stars already in
`results/track_a/stars/` are skipped) and pushes every 30 minutes. Wait with the same single background
command, using `/home/user/work/worker_a.log` and `pipeline/[w]orker_a.sh` in place of the survey names. Final
report: rows done, stars searched, any errors.
