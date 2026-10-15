"""Paths and process settings shared by every pipeline script.

Bulk data (light curves, catalogs, caches) lives outside the repository in WORK; only compact
results, figures and docs are committed.
"""
import os
import warnings
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORK = Path(os.environ.get("ADASTRA_WORK", REPO.parent / "work"))
CATALOGS = WORK / "catalogs"
CACHE = WORK / "cache"
LCCACHE = WORK / "lc"
SURVEY = WORK / "survey"
RESULTS = REPO / "results"
FIGURES = REPO / "results" / "figures"

for d in (CATALOGS, CACHE, LCCACHE, SURVEY, RESULTS, FIGURES, CACHE / "mpl"):
    d.mkdir(parents=True, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(CACHE / "numba"))
# Parallelism comes from process pools; keep numerical libraries single-threaded inside workers.
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(var, "1")
warnings.filterwarnings("ignore", message=".*oktopus.*")
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*(Mean of empty|invalid value|divide by zero).*")
