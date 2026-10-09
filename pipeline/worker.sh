#!/usr/bin/env bash
# Survey worker: set up, then run run_sector.py for each sector given and push results after each.
#   bash pipeline/worker.sh 12 13 14 ...
# Logs: $ADASTRA_WORK/worker.log (progress), $ADASTRA_WORK/run_<S>.log (per sector).
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ADASTRA_WORK="${ADASTRA_WORK:-/home/user/work}"
BRANCH="claude/gifted-albattani-szay6x"
PY=/home/user/venv/bin/python
LOG="$ADASTRA_WORK/worker.log"
mkdir -p "$ADASTRA_WORK"
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }

cd "$REPO"
git checkout -q "$BRANCH" 2>/dev/null || git checkout -q -b "$BRANCH" "origin/$BRANCH"
git pull -q --rebase origin "$BRANCH" || true

if [ ! -x "$PY" ]; then
  log "creating venv"
  python3 -m venv /home/user/venv && /home/user/venv/bin/pip install -q -r requirements.txt >> "$LOG" 2>&1
fi

cd "$REPO/pipeline"
if [ ! -f "$ADASTRA_WORK/catalogs/tic_v2/_SUCCESS" ]; then
  log "building TIC extract from S3 in the background (run_sector waits for it before reducing)"
  ( $PY build_tic.py 8 >> "$ADASTRA_WORK/tic.log" 2>&1 || $PY build_tic.py 4 >> "$ADASTRA_WORK/tic.log" 2>&1
    echo "$(date -u +%FT%TZ) TIC extract finished: $(tail -1 "$ADASTRA_WORK/tic.log")" >> "$LOG" ) &
fi

push_results() {
  local s4="$1"
  cd "$REPO"
  git add "results/survey/$s4" >> "$LOG" 2>&1
  git commit -q -m "Survey results for sector ${s4#s}" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" >> "$LOG" 2>&1
  for i in 1 2 3 4 5 6; do
    if git pull -q --rebase origin "$BRANCH" >> "$LOG" 2>&1 && git push -q origin "$BRANCH" >> "$LOG" 2>&1; then
      log "pushed $s4"; break
    fi
    sleep $((i * 7))
  done
  cd "$REPO/pipeline"
}

for S in "$@"; do
  s4=$(printf "s%04d" "$S")
  if [ -f "$REPO/results/survey/$s4/summary.json" ]; then log "$s4 already done"; continue; fi
  for attempt in 1 2; do
    log "start $s4 (attempt $attempt)"
    if $PY run_sector.py "$S" 4 > "$ADASTRA_WORK/run_$S.log" 2>&1; then
      log "done $s4: $(tail -1 "$ADASTRA_WORK/run_$S.log" | cut -c1-300)"
      push_results "$s4"
      break
    else
      log "FAILED $s4 attempt $attempt"
      rm -rf "$ADASTRA_WORK/survey/$s4"
      if [ "$attempt" = 2 ]; then
        mkdir -p "$REPO/results/survey/$s4"
        tail -40 "$ADASTRA_WORK/run_$S.log" > "$REPO/results/survey/$s4/FAILED.txt"
        push_results "$s4"
      fi
    fi
  done
done
log "ALL DONE: $*"
