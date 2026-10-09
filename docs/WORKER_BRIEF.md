# Worker brief: survey sectors for the ad-astra TESS search

You are one of several cloud sessions running the same survey pipeline on different TESS sectors.
A coordinating session wrote this pipeline, validated it, and assigned you a list of sectors. Your job
is mechanical: run the pipeline on your sectors and push the compact results. Do not change the
pipeline code. If something is broken, record it (see "Failures") and move on; the coordinator will
read your notes.

## What the pipeline does (context only)
`pipeline/run_sector.py S` reads every QLP full-frame light curve of sector S, for stars in
`results/samples/survey_sample.parquet`, straight from the public MAST S3 bucket (stpubdata). It
flattens each light curve, catalogues every significant dip, flags dips shared by many stars
(spacecraft systematics), blends and artefacts, then writes candidate tables for three science tracks
to `results/survey/sNNNN/`. About 10% of light curves are also searched upside-down to measure the
false-alarm rate. Only S3 is needed; the other astronomy archives are blocked in this environment.

## Setup (once)
```bash
cd <repo root>                 # the checked-out ad-astra repository
git checkout claude/gifted-albattani-szay6x && git pull origin claude/gifted-albattani-szay6x
python3 -m venv /home/user/venv
/home/user/venv/bin/pip install -q -r requirements.txt
export ADASTRA_WORK=/home/user/work
cd pipeline
nohup /home/user/venv/bin/python build_tic.py 8 > /home/user/work_tic.log 2>&1 &   # ~15-25 min, needed by the reduce step
```
`run_sector.py` waits for the TIC table by itself only if you start it after `build_tic.py` finished.
So: start `build_tic.py`, wait for its log to print `DONE` (use a Monitor or a background
`until grep -q DONE ...; do sleep 20; done` command — never foreground sleep loops), then run
`/home/user/venv/bin/python build_manifest.py` is NOT needed; only `tic_master.parquet` is. Build it with:
```bash
/home/user/venv/bin/python -c "import build_manifest as b; b.master()"
```

## Run your sectors
For each sector S in your list, in order:
```bash
cd <repo root>/pipeline
ADASTRA_WORK=/home/user/work nohup /home/user/venv/bin/python run_sector.py S 4 > /home/user/work/run_S.log 2>&1
```
Run it in the background and wait for it with a Monitor/until-loop on the log (look for the final JSON
summary line or a Traceback). One sector takes roughly 0.5-2 hours. Then commit and push its results:
```bash
cd <repo root>
git add results/survey/sNNNN
git commit -m "Survey results for sector S" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
for i in 1 2 3 4 5; do git pull --rebase -q origin claude/gifted-albattani-szay6x && git push -q origin claude/gifted-albattani-szay6x && break; sleep $((i*5)); done
```
Only ever add files under `results/survey/`. Never modify, revert or delete other files; never force-push;
never open a pull request.

## Failures
If a sector fails, read the log, retry it once. If it fails again, write the last 40 lines of the log to
`results/survey/sNNNN/FAILED.txt`, commit and push that file, and continue with the next sector.
Disk: delete `/home/user/work/survey/sNNNN` if a failed run left it behind.

## When done
Reply with a short table: sector, number of light curves, n_trackb, n_trackc, n_trackd, runtime, status.
