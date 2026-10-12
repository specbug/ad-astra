#!/usr/bin/env python3
"""Offline review and freeze of the separately declared S1/S2 overlap field."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import random

from adastra.archive import sha256_file, utc_now, write_json
from qlp_survey import FROZEN_STATUS, source_hashes, validate_frozen
import select_qlp_pilot as base


def read(path):
    return json.loads(Path(path).read_text())


def checked(path, digest):
    if sha256_file(Path(path)) != digest:
        raise ValueError(f"Recorded hash mismatch: {path}")


def review(directory):
    directory = Path(directory)
    provisional = directory / 'provisional-manifest.json'
    manifest = read(provisional)
    declaration_path = directory / 'selection-declaration.json'
    declaration = read(declaration_path)
    declaration_sha = sha256_file(declaration_path)
    checked(declaration_path, manifest['selection_declaration_sha256'])
    checked(provisional, read(directory / 'provisional-manifest.provenance.json')['sha256'])
    for path, digest in declaration['sources_sha256'].items():
        checked(path, digest)
    checked('scripts/select_qlp_overlapfield.py', manifest['selection_script_sha256'])
    for record in declaration['catalog_snapshot_sources']:
        checked(record['path'], record['sha256'])
    checked(declaration['literature_exclusions_path'], declaration['literature_exclusions_sha256'])
    for key in ('known_science_manifest', 'previous_coverage_pilot'):
        checked(declaration[key]['path'], declaration[key]['sha256'])
    if (manifest['selected_count'] != 24 or manifest['requested_count'] != 24
            or manifest['flux_accessed'] is not False or manifest['selection_complete'] is not True
            or declaration['science_flux_accessed'] is not False
            or declaration['discovery_sectors'] != [1, 2]):
        raise ValueError('Expected complete, metadata-only, literal S1/S2 selection')
    catalog = read(directory / 'catalog-ledger.json')
    pool = sorted((row['tic'] for row in catalog if not row['excluded_reasons']))
    random.Random(declaration['seed']).shuffle(pool)
    if pool != read(directory / 'shuffled-pool.json') or len(pool) != manifest['post_exclusion_pool']:
        raise ValueError('Pool or deterministic shuffle changed')
    ledger = read(directory / 'selection-ledger.json')
    if len(ledger) != manifest['seed_ranks_examined']:
        raise ValueError('Examined count differs from ledger')
    eligible = []
    for rank, state in enumerate(ledger, 1):
        if state['rank'] != rank or state['tic'] != pool[rank - 1] or state['declaration_sha256'] != declaration_sha:
            raise ValueError('Selection ledger rank, identity, or declaration mismatch')
        saved = read(directory / 'target-state' / f"{state['tic']}.json")
        if saved != {key: value for key, value in state.items() if key != 'selected'}:
            raise ValueError('Selection ledger conflicts with persisted target state')
        for relative, digest in state['metadata_files_sha256'].items():
            checked(directory / relative, digest)
        if state['eligible']:
            if state['excluded_reasons']:
                raise ValueError('Eligible target has exclusion reasons')
            eligible.append(state['target'])
        if state['selected'] != (state['eligible'] and len(eligible) <= 24):
            raise ValueError('Selected flags violate first24 eligible rule')
    if eligible[:24] != manifest['targets']:
        raise ValueError('Manifest is not the first24 eligible targets in seed order')
    for target in manifest['targets']:
        tic = target['tic']
        gaia = read(target['gaia_audit'])
        if (gaia['query_status'] != 'complete' or gaia['radius_arcsec'] != 90
                or gaia['variable_neighbors'] or gaia['excluded_reasons']):
            raise ValueError('Gaia preselection filter is not complete and passing')
        observations = []
        snapshots = []
        for path in sorted((directory / 'observations' / str(tic)).glob('response-[0-9][0-9][0-9].json')):
            provenance = read(path.with_name(path.name + '.provenance.json'))
            checked(path, provenance['sha256'])
            response = read(path)
            if response['status'] != 'COMPLETE':
                raise ValueError('Incomplete observation metadata')
            observations.extend(response['data'])
            snapshots.append((path, provenance, response['data']))
        usable, reason = base.eligible_observations(observations)
        if reason or not usable:
            raise ValueError('Selected target fails recomputed observation eligibility')
        for group in ('discovery', 'holdout'):
            bare = [{key: value for key, value in product.items() if key != 'provenance_name'} for product in target[group]]
            if bare != usable[group]:
                raise ValueError('Products differ from raw-metadata season selection')
            for product in target[group]:
                matches = [(path, prov, row) for path, prov, rows in snapshots for row in rows
                           if row['obsid'] == product['obsid'] and row['dataURL'] == product['uri']]
                if not matches:
                    raise ValueError('No exact raw observation identity')
                for _, _, row in matches:
                    if (row['provenance_name'] != 'QLP' or int(row['target_name']) != tic
                            or int(row['sequence_number']) != product['sector']
                            or row['t_min'] != product['t_min'] or row['t_max'] != product['t_max']):
                        raise ValueError('Raw observation identity conflict')
                path, prov, row = matches[0]
                product.update(metadata_response=str(path), metadata_response_sha256=prov['sha256'],
                               metadata_retrieved_utc=prov['retrieved_utc'])
    manifest['selection_review_artifacts_sha256'] = {
        str(path): sha256_file(path) for path in directory.rglob('*') if path.is_file()}
    return manifest


def main():
    output = Path('data/qlp-overlapfield-frozen')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Refusing to overwrite an existing freeze')
    manifest = review('data/qlp-overlapfield')
    previous_config = Path('data/qlp-frozen/config.json')
    config = deepcopy(read(previous_config))
    if config['source_sha256'] != source_hashes():
        raise ValueError('Reviewed screen source changed')
    checked(config['control_results'], config['control_results_sha256'])
    controls = read(config['control_results'])
    for path, digest in controls['declaration']['sources_sha256'].items():
        checked(path, digest)
    if controls['summary']['cases'] != 18 or controls['summary']['exact_period_and_phase_recoveries'] != 18:
        raise ValueError('Injection control outcome changed')
    manifest.update(status=FROZEN_STATUS, review_status='approved', frozen_utc=utc_now(),
                    provisional_manifest_sha256=sha256_file(Path('data/qlp-overlapfield/provisional-manifest.json')),
                    freeze_helper_sha256=sha256_file(Path(__file__)),
                    review_meaning='Approved for exploratory screening; exact metadata and ordered selection verified; no novelty claim')
    config.update(frozen_utc=utc_now(), reused_reviewed_config=str(previous_config),
                  reused_reviewed_config_sha256=sha256_file(previous_config),
                  novelty_policy='Scoped NASA/CTOI, SIMBAD and literature exclusions; Gaia DR3 VARIABLE neighbors within90 arcsec excluded. Further literature review remains necessary; no exhaustive novelty clearance.',
                  selection_limitations='Gaia NOT_AVAILABLE/unknown flags and multiple or absent close Gaia sources remain unchecked. Epoch2016 positions are not propagated. Three metadata HTTP503 failures excluded without replacement or retry, creating availability bias. Conservative neighborhood filtering also biases the sample.',
                  protected_holdout_policy='Sectors68/69 remain metadata-only; science flux only literal sectors1/2.')
    serialized = json.dumps(manifest, indent=2, allow_nan=False) + '\n'
    config['manifest_sha256'] = sha256(serialized.encode()).hexdigest()
    validate_frozen(manifest, config, config['manifest_sha256'], source_hashes())
    write_json(output / 'manifest.json', manifest)
    checked(output / 'manifest.json', config['manifest_sha256'])
    write_json(output / 'config.json', config)
    print(f'Frozen24 targets offline: {output}; protected sectors68/69 unopened')


if __name__ == '__main__':
    main()
