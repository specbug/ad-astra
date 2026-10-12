#!/usr/bin/env python3
"""Explore fixed published neighbor periods and exact daily aliases in S1/S2."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adastra.archive import write_json, utc_now, sha256_file
from qlp_variability_followup import raw_columns, column_arrays, design, fit, EPOCH

PERIODS = {149392671: .47711359612012993, 150104580: .5442901668164593}


def joint_fit(t, y, e, rr, alias):
    columns = [design(t, rr/h)[:, :2] for h in (1, 2, 3)]
    columns += [design(t, alias)[:, :2], design(t, rr)[:, 2:]]
    x = np.column_stack(columns); w = 1/e
    beta, _, rank, _ = np.linalg.lstsq(x*w[:, None], y*w, rcond=None)
    if rank != x.shape[1]:
        raise ValueError('Joint design rank deficient')
    return {'rr_harmonic_amplitudes_ppm': [float(np.hypot(*beta[i:i+2])) for i in (0, 2, 4)],
            'alias_amplitude_ppm': float(np.hypot(*beta[6:8])),
            'alias_phase_cycles': float(np.arctan2(beta[7], beta[6])/(2*np.pi)),
            'rr_fundamental_phase_cycles': float(np.arctan2(beta[1], beta[0])/(2*np.pi)),
            'residual_rms_ppm': float(np.sqrt(np.mean((y-x@beta)**2)))}


def leakage(t, e, injected_period, fitted_period):
    """Worst-phase linear transfer from unit sinusoid at one fixed frequency."""
    x = design(t, fitted_period); z = design(t, injected_period)[:, :2]; w=1/e
    coefficients = np.linalg.lstsq(x*w[:, None], z*w[:, None], rcond=None)[0][:2]
    return float(np.linalg.svd(coefficients, compute_uv=False)[0])


def fold(ax, t, y, period, title):
    phase = ((t-EPOCH)/period+.5)%1-.5
    ax.scatter(phase, y, s=2, alpha=.2)
    edges=np.linspace(-.5,.5,41); b=np.digitize(phase, edges)-1
    points=[((edges[i]+edges[i+1])/2, np.median(y[b==i])) for i in range(40) if np.count_nonzero(b==i)>=3]
    if points: ax.plot(*np.array(points).T, 'o', ms=3, color='darkorange')
    ax.set(title=title, xlabel='Phase at BTJD 1340', ylabel='ppm after segment baselines')


def main():
    out=Path('runs/qlp-catalog-period-tests')
    if out.exists() and any(out.iterdir()): raise ValueError('Preserve previous output')
    path=Path('runs/qlp-pilot/results.json');data=json.loads(path.read_text())
    write_json(out/'declaration.json', {'created_utc':utc_now(), 'source_sha256':sha256_file(__file__),
        'pilot_sha256':sha256_file(path), 'stage':'Postselection fixed published-neighbor-period diagnostics',
        'rr_periods_days':PERIODS, 'alias_formula':'1/abs(1/P_rr - 1/day)',
        'joint_model':'RR fundamental plus second/third harmonics, exact daily alias, per-gap segment offset/slope',
        'window_test':'Fit unit sin/cos at RR frequency onto exact daily alias using actual sample times and weights; maximum singular value is worst-phase amplitude leakage',
        'limitation':'Sampling/linear-baseline leakage only; does not replay nonlinear upstream QLP corrections; no FAP or source localization',
        'allowed_sectors':[1,2], 'protected_holdouts_accessed':False})
    result=[]
    for tic, rr in PERIODS.items():
        target=next(r for r in data['results'] if r['tic']==tic)
        fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
        for i,row in enumerate(target['sectors']):
            _,columns,cam,ccd=raw_columns(row,tic);t,y,e,_=column_arrays(columns,'KSPSAP_FLUX')
            alias=1/abs(1/rr-1)
            rr_fit, rr_nuisance=fit(t,y,rr,e);alias_fit,alias_nuisance=fit(t,y,alias,e)
            lead=row['analysis']['lomb_scargle']['period_days']
            joint=joint_fit(t,y,e,rr,alias)
            row_out={'tic':tic,'sector':row['sector'],'rr_period':rr,'daily_alias_period':alias,
                     'screen_period':lead,'rr_single_fit':rr_fit,'alias_single_fit':alias_fit,'joint_fit':joint,
                     'rr_harmonic_to_alias_window_leakage':[leakage(t,e,rr/h,alias) for h in (1,2,3)]}
            result.append(row_out)
            fold(axes[i,0],t,y-rr_nuisance,rr,f'S{row["sector"]}: published RR period {rr:.8f} d')
            fold(axes[i,1],t,y-alias_nuisance,alias,f'S{row["sector"]}: daily alias {alias:.8f} d')
        fig.suptitle(f'TIC {tic}: aperture photometry; periods of known nearby stars tested after selection\nNo source localization or independent validation; orange = 40-bin medians')
        fig.savefig(out/f'tic-{tic}.png',dpi=150);plt.close(fig)
    write_json(out/'results.json',{'completed_utc':utc_now(),'results':result})
    for row in result:
        print(row['tic'],row['sector'],'RR',round(row['rr_single_fit']['amplitude'],2),'alias',round(row['alias_single_fit']['amplitude'],2),
              'leakage',row['rr_harmonic_to_alias_window_leakage'],'joint',row['joint_fit'],flush=True)


if __name__=='__main__':main()
