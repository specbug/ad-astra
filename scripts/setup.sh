#!/usr/bin/env bash
set -euo pipefail
cd /workspace/ad-astra
if [[ ! -x .venv/bin/python ]]; then
  python -m venv .venv
fi
.venv/bin/python -m pip install --disable-pip-version-check --no-cache-dir -r requirements-validated.txt
.venv/bin/python -m pip install --disable-pip-version-check --no-deps -e .
if [[ "${1:-}" != "--install-only" ]]; then
  .venv/bin/python -m pytest -q
fi

