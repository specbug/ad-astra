#!/usr/bin/env python3
"""Postselection S1/S2 fixed-frequency diagnostics; never retrieves new flux."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import json
import argparse
import numpy as np
from adastra.archive import validated_fits, sha256_file, utc_now, write_json
from qlp_diagnostics import load_cached_sector, ROOT

LEADS = (149392671, 150104580)
EPOCH = 1340.0


def design(time, period):
    cuts = np.r_[0, np.flatnonzero(np.diff(time) > .5) + 1, len(time)]
    columns = [np.sin(2*np.pi*(time-EPOCH)/period), np.cos(2*np.pi*(time-EPOCH)/period)]
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        indicator = np.zeros(len(time)); indicator[lo:hi] = 1
        slope = np.zeros(len(time)); slope[lo:hi] = time[lo:hi] - np.mean(time[lo:hi])
        columns.append(indicator)
        # An isolated sample has an identifiable offset, but no independent
        # slope. Retain the sample without adding an all-zero design column.
        if hi - lo > 1:
            columns.append(slope)
    return np.array(columns).T


def fit(time, values, period, error=None, robust=False):
    x = design(time, period)
    base_weight = np.ones(len(time)) if error is None else 1 / error
    rw = np.ones(len(time)); converged = not robust
    for iteration in range(30 if robust else 1):
        w = base_weight*np.sqrt(rw)
        beta, _, rank, _ = np.linalg.lstsq(x*w[:, None], values*w, rcond=None)
        if rank != x.shape[1]:
            raise ValueError('Rank-deficient fixed-frequency design')
        resid = (values-x@beta)*base_weight
        if not robust:
            break
        center = np.median(resid); scale = 1.4826*np.median(np.abs(resid-center))
        if scale == 0:
            converged = True; break
        new = np.minimum(1, 2.5*scale/np.maximum(np.abs(resid-center), np.finfo(float).eps))
        if np.max(np.abs(new-rw)) < 1e-6:
            converged = True; break
        rw = new
    sine, cosine = beta[:2]
    residual = values-x@beta
    nuisance = x[:, 2:]@beta[2:]
    amplitude = float(np.hypot(sine, cosine))
    return {'amplitude': amplitude, 'sine': float(sine), 'cosine': float(cosine),
            'phase_cycles_at_btjd1340': float(np.arctan2(cosine, sine)/(2*np.pi)),
            'residual_rms': float(np.sqrt(np.mean(residual**2))), 'samples': len(time),
            'rank': int(rank), 'robust_huber': robust, 'converged': converged,
            'iterations': iteration+1, 'model': 'sin+cos at fixed S1 frequency, independent offset and linear slope for each >0.5d gap segment'}, nuisance


def raw_columns(row, tic):
    curve = load_cached_sector(row, tic, ROOT)
    with validated_fits(row['input_path']) as hdus:
        table = hdus[1].data
        columns = {name: np.array(table[name], dtype=float) for name in table.names}
        camera, ccd = hdus[0].header.get('CAMERA'), hdus[0].header.get('CCD')
    # QLP reader has already validated TDB/BJD reference. The table is TIME in
    # BTJD for these exact products; bind to its parsed times independently.
    from adastra.qlp import read_qlp
    lc = read_qlp(row['input_path'])
    columns['TIME'] = lc.time_btjd
    return curve, columns, camera, ccd


def column_arrays(columns, name):
    t, f, e, q = (columns[k] for k in ('TIME', name, 'KSPSAP_FLUX_ERR', 'QUALITY'))
    good = (q == 0) & np.isfinite(t) & np.isfinite(f)
    photometry = name in ('KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG')
    if photometry:
        good &= f > 0
    if name == 'KSPSAP_FLUX':
        good &= np.isfinite(e) & (e > 0)
    t, f, e = t[good], f[good], e[good]
    order = np.argsort(t, kind='stable'); t, f, e = t[order], f[order], e[order]
    keep = np.r_[True, np.diff(t) > 0]; t, f, e = t[keep], f[keep], e[keep]
    if len(t) < 100:
        raise ValueError('Too few usable samples')
    if photometry:
        norm = np.median(f); f = (f/norm - 1)*1e6; e = e/norm*1e6
    else:
        f = f-np.median(f)
    return t, f, e if name == 'KSPSAP_FLUX' else None, 'ppm' if photometry else ('pixel' if name in ('SAP_X', 'SAP_Y') else 'archive_background_units')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('runs/qlp-variability-followup'))
    args = parser.parse_args()
    path = Path('runs/qlp-pilot/results.json'); pilot = json.loads(path.read_text())
    if pilot['status'] != 'complete' or not pilot['original_sources_unchanged']:
        raise ValueError('Need complete original pilot')
    out = args.output
    if out.exists() and any(out.iterdir()):
        raise ValueError('Preserving existing followup')
    targets = {t['tic']: t for t in pilot['results']}
    periods = {tic: next(s for s in targets[tic]['sectors'] if s['sector'] == 1)['analysis']['lomb_scargle']['period_days'] for tic in LEADS}
    declaration = {'created_utc': utc_now(), 'stage': 'Postselection exploratory diagnostics; both sectors already screened',
        'source_sha256': sha256_file(__file__), 'pilot_sha256': sha256_file(path),
        'lead_periods_fixed_from_sector1': periods, 'reference_epoch_btjd': EPOCH,
        'photometric_columns': ['KSPSAP_FLUX', 'SAP_FLUX', 'KSPSAP_FLUX_SML', 'KSPSAP_FLUX_LAG'],
        'ancillary_columns': ['SAP_X', 'SAP_Y', 'SAP_BKG'],
        'methods': 'Fixed-frequency linear regression with per-gap segment offsets and linear slopes (offset only for isolated one-sample segments); weighted only for explicit KSPSAP errors. Parallel Huber2.5 IRLS; no sample clipping. Leave-one-one-day-bin-out ordinary fits for leads. All24curves evaluated at both lead frequencies.',
        'limitations': 'No FAP, independent validation, source localization or intrinsic variability classification. Alternative QLP apertures are already processed. SAP/centroid/background have no borrowed KSPSAP errors. Comparison stars were selected in this same exploratory sample.',
        'allowed_sectors': [1, 2], 'protected_holdouts_accessed': False}
    write_json(out/'declaration.json', declaration)
    results, comparisons = [], []
    for tic, target in sorted(targets.items()):
        for row in target['sectors']:
            curve, columns, camera, ccd = raw_columns(row, tic)
            time, flux, error, units = column_arrays(columns, 'KSPSAP_FLUX')
            for lead, period in periods.items():
                ordinary, _ = fit(time, flux, period, error)
                robust, _ = fit(time, flux, period, error, robust=True)
                comparisons.append({'tic': tic, 'sector': row['sector'], 'camera': camera, 'ccd': ccd,
                                    'frequency_selected_for_tic': lead, 'period_days': period,
                                    'ordinary': ordinary, 'robust': robust})
            if tic not in LEADS:
                continue
            fits = {}
            for name in declaration['photometric_columns'] + declaration['ancillary_columns']:
                try:
                    t, y, e, units = column_arrays(columns, name)
                    ordinary, _ = fit(t, y, periods[tic], e)
                    robust, _ = fit(t, y, periods[tic], e, robust=True)
                    diagnostic = {'units': units, 'ordinary': ordinary, 'robust': robust}
                    if name in ('KSPSAP_FLUX', 'SAP_FLUX'):
                        days = np.floor(t-t[0]).astype(int); daily = []
                        for day in np.unique(days):
                            keep = days != day
                            a, _ = fit(t[keep], y[keep], periods[tic], None if e is None else e[keep])
                            daily.append({'omitted_day_index': int(day), 'fit': a})
                        diagnostic['leave_one_day_out'] = daily
                    fits[name] = diagnostic
                except Exception as exc:
                    fits[name] = {'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}
            results.append({'tic': tic, 'sector': row['sector'], 'camera': camera, 'ccd': ccd,
                            'input_sha256': row['input_sha256'], 'period_days': periods[tic], 'columns': fits})
            write_json(out/'partial-results.json', {'status': 'partial', 'leads': results, 'comparisons': comparisons})
            print(f'Fixed-frequency diagnostics TIC {tic} S{row["sector"]}', flush=True)
    write_json(out/'results.json', {'declaration': declaration, 'completed_utc': utc_now(), 'leads': results, 'comparisons': comparisons})


if __name__ == '__main__':
    main()
