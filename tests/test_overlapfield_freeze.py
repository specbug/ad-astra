"""The offline freeze must fail closed on altered raw metadata or selection order."""
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
SPEC = importlib.util.spec_from_file_location('freeze_overlap', ROOT / 'scripts/freeze_qlp_overlapfield.py')
freeze = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(freeze)


@pytest.fixture
def archived_selection(tmp_path, monkeypatch):
    source = ROOT / 'data/qlp-overlapfield'
    if not source.exists():
        pytest.skip('Archived metadata experiment unavailable')
    monkeypatch.chdir(ROOT)
    destination = tmp_path / 'selection'
    shutil.copytree(source, destination)
    return destination


def test_reordered_pool_rejected(archived_selection):
    path = archived_selection / 'shuffled-pool.json'
    pool = json.loads(path.read_text())
    pool[0], pool[1] = pool[1], pool[0]
    path.write_text(json.dumps(pool))
    with pytest.raises(ValueError, match='deterministic shuffle'):
        freeze.review(archived_selection)


def test_changed_raw_metadata_rejected(archived_selection):
    state = json.loads((archived_selection / 'selection-ledger.json').read_text())[0]
    raw = next(relative for relative in state['metadata_files_sha256']
               if relative.startswith('observations/') and relative.endswith('.json')
               and not relative.endswith('.provenance.json'))
    path = archived_selection / raw
    path.write_text(path.read_text() + ' ')
    with pytest.raises(ValueError, match='Recorded hash mismatch'):
        freeze.review(archived_selection)


def test_promoted_target_state_rejected(archived_selection):
    path = archived_selection / 'selection-ledger.json'
    states = json.loads(path.read_text())
    rejected = next(state for state in states if not state['eligible'])
    rejected['eligible'] = True
    rejected['selected'] = True
    path.write_text(json.dumps(states))
    with pytest.raises(ValueError, match='persisted target state'):
        freeze.review(archived_selection)
