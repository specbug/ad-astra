#!/usr/bin/env python3
"""Build the scoped pre-flux exclusions from preserved public catalog responses."""
import json
from pathlib import Path
from adastra.archive import sha256_file, utc_now, write_json


def tsv_table(path, name):
    rows, header, skip, selected = [], None, 0, False
    for line in path.read_text().splitlines():
        if line.startswith('#Name: '):
            selected = line[7:].strip() == name
            header = None
            continue
        if not selected or not line or line.startswith('#'):
            continue
        values = [x.strip() for x in line.split('\t')]
        if header is None:
            header, skip = values, 2
        elif skip:
            skip -= 1
        else:
            if len(values) != len(header):
                raise ValueError(f'Malformed table row in {path}')
            rows.append(dict(zip(header, values)))
    if not rows:
        raise ValueError(f'No rows in {name}')
    return rows


def main():
    root = Path('data/qlp-prior-review')
    output = root / 'literature-exclusions.json'
    if output.exists():
        raise ValueError('Preserving existing exclusion declaration')
    avallone = Path('data/qlp-pilot-v2/literature/2022ApJ-930-7-table1.tsv')
    chespa = root / 'chespa-tables.tsv'
    review = root / 'catalog-review-v1.json'
    sources = [avallone, chespa, review, root/'avallone-ReadMe.txt', root/'chespa-ReadMe.txt',
               root/'macho-31380486.tsv', root/'macho-41028187.tsv', root/'macho-ReadMe.txt']
    for path in sources:
        sidecar = path.with_name(path.name + '.provenance.json')
        if sidecar.exists() and json.loads(sidecar.read_text())['sha256'] != sha256_file(path):
            raise ValueError(f'Raw response changed: {path}')
    arows = tsv_table(avallone, 'J/ApJ/930/7/table1')
    crows = tsv_table(chespa, 'J/ApJS/240/17/cand')
    tess = [r for r in arows if r['Mission'] == 'TESS']
    if len(arows) != 4216 or len(tess) != 2115 or len(crows) != 222:
        raise ValueError(f'Unexpected public table sizes: {len(arows)}, {len(tess)}, {len(crows)}')
    exclusions = {}
    def add(tic, reason):
        exclusions.setdefault(str(int(tic)), []).append(reason)
    for row in tess:
        add(row['ID'], 'Avallone2022_TESS_class:' + row['Class'])
    for row in crows:
        if row['TIC']:
            add(row['TIC'], 'CHESPA2019_reported_signal:' + row['Label'] + ':P=' + row['Per'] + 'd')
    add(41028187, 'Fraser2008_MACHO73.12999.9_prior_periods_18.01_and_7.49d_unclassified')
    add(31380486, 'Fraser2008_MACHO28.11765.344_one_year_artifact_flag_not_confirmed_variability')
    matches = []
    for target in json.loads(review.read_text())['variable_catalogs']:
        if target['status'] != 'complete':
            continue
        for table in target['tables']:
            for row in table['rows']:
                # _r is explicitly arcsec in the preserved VizieR XML.
                if row.get('_r') is not None and row['_r'] <= 5:
                    add(target['tic'], 'initial_cohort_catalog_match:' + table['catalog'] + ':' + str(row.get('Type', '')))
                    matches.append({'tic': target['tic'], 'catalog': table['catalog'], 'row': row})
    write_json(output, {
        'created_utc': utc_now(), 'science_flux_accessed': False,
        'builder_sha256': sha256_file(__file__),
        'rules': [
            'Exclude all exact TIC identities in Avallone2022 TESS table1; preserve original classes. This conservative study-membership filter does not imply every entry has a measured period or confirmed binary.',
            'Exclude exact TIC identities in all 222 returned CHESPA2019 table rows, including uncertain and likely-binary labels; the paper title refers to 116 candidates, not the full returned table size.',
            'Exclude two exact SIMBAD MACHO aliases verified against Fraser2008 tables: TIC41028187 prior period report; TIC31380486 prior one-year artifact flag.',
            'Exclude <=5 arcsec variable-catalog matches in the completed initial 24-target review. This additional review is scoped to that cohort; replacement catalog queries occur separately.',
        ],
        'sources': [{'path': str(p), 'sha256': sha256_file(p)} for p in sources],
        'catalog_counts': {'Avallone_table1': len(arows), 'Avallone_TESS': len(tess), 'CHESPA': len(crows)},
        'initial_cohort_variable_matches': matches,
        'exclusions': {k: sorted(set(v)) for k, v in sorted(exclusions.items(), key=lambda kv: int(kv[0]))},
        'limitations': 'Scoped pre-flux selection with catalog-dependent bias; not exhaustive literature clearance; coordinates are not propagated across epochs. RUWE>1.2 is a potential-binary flag, not proof of a physical binary.'})
    print(f'Saved {len(exclusions)} catalog identities to {output}')


if __name__ == '__main__':
    main()
