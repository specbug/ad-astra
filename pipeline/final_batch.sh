#!/usr/bin/env bash
# End-of-survey processing, in order. Each step logs to $LOG; a failed step is reported and the rest still run.
#   bash pipeline/final_batch.sh
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
PY=/home/user/venv/bin/python
LOG=/home/user/work/logs/final_batch.log
step() { echo "== $(date -u +%T) $1" | tee -a "$LOG"; shift; "$@" >> "$LOG" 2>&1 || echo "   FAILED: $*" | tee -a "$LOG"; }

cd .. && git pull -q --rebase origin claude/gifted-albattani-szay6x; cd pipeline
step "Track B stage 2"   $PY -c "import aggregate; aggregate.track_b()"
step "Track C"           $PY -c "import aggregate; aggregate.track_c()"
step "Track D stage 2"   $PY -c "import aggregate; aggregate.track_d()"
step "finalize B"        $PY finalize_b.py 6
step "finalize D"        $PY finalize_d.py 6
step "known TOIs"        $PY -c "
import env, known, pandas as pd
from env import RESULTS
C = RESULTS / 'candidates'
for n in ('trackb_candidates.csv', 'trackb_final.csv', 'trackc_candidates.csv', 'trackd_final.csv'):
    d = pd.read_csv(C / n)
    if 'ra' not in d:
        continue
    d = known.flag_known(d.drop(columns=[c for c in ('toi_exofop', 'toi_hf', 'is_toi') if c in d]))
    d.to_csv(C / n, index=False, float_format='%.6g')
    print(n, len(d), 'TOIs', int(d.is_toi.sum()))
"
step "single-transit periods" $PY vet_singles.py 4
step "deep vetting (T<10.5)"  $PY vet_deep.py --shortlist 10.5
step "completeness"      $PY completeness.py
step "Track A triage"    $PY triage_a.py ../results/track_a ../results/track_a/targets_search.csv
step "report extras"     $PY report_extras.py
step "report"            $PY report.py
echo "== $(date -u +%T) final batch done" | tee -a "$LOG"
