#!/usr/bin/env bash
# Track A worker: multi-sector BLS search of rows [START, STOP) of results/track_a/targets_search.csv.
#   bash pipeline/worker_a.sh START STOP
# Star results go to results/track_a/stars/<tic>.json and are pushed every ~30 minutes.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export ADASTRA_WORK="${ADASTRA_WORK:-/home/user/work}"
BRANCH="claude/gifted-albattani-szay6x"
PY=/home/user/venv/bin/python
LOG="$ADASTRA_WORK/worker_a.log"
mkdir -p "$ADASTRA_WORK"
log() { echo "$(date -u +%FT%TZ) $*" | tee -a "$LOG"; }
cd "$REPO"
git checkout -q "$BRANCH" 2>/dev/null || git checkout -q -b "$BRANCH" "origin/$BRANCH"
git pull -q --rebase origin "$BRANCH" || true
if [ ! -x "$PY" ]; then
  log "creating venv"
  python3 -m venv /home/user/venv && /home/user/venv/bin/pip install -q -r requirements.txt >> "$LOG" 2>&1
fi
push() {
  cd "$REPO"
  git add results/track_a/stars >> "$LOG" 2>&1
  if ! git diff --cached --quiet; then
    git commit -q -m "Track A search results ($1)" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01HPQdHQVehmk4raGYZB42gv" >> "$LOG" 2>&1
    for i in 1 2 3 4 5 6; do
      if git pull -q --rebase origin "$BRANCH" >> "$LOG" 2>&1 && git push -q origin "$BRANCH" >> "$LOG" 2>&1; then
        log "pushed ($1)"; break
      fi
      sleep $((i * 7))
    done
  fi
}
cd "$REPO/pipeline"
log "start Track A rows $1-$2"
$PY run_track_a.py "$REPO/results/track_a/targets_search.csv" "$REPO/results/track_a" 4 "$1" "$2" > "$ADASTRA_WORK/track_a_$1.log" 2>&1 &
PID=$!
while kill -0 $PID 2>/dev/null; do
  sleep 1800
  log "progress: $(tail -1 "$ADASTRA_WORK/track_a_$1.log" | cut -c1-120)"
  push "rows $1-$2, in progress"
done
wait $PID; rc=$?
push "rows $1-$2, finished rc=$rc"
log "ALL DONE rows $1-$2 rc=$rc"
