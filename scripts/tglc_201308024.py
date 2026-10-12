#!/usr/bin/env python3
"""Declared native TGLC fixed-frequency reduction comparison on S1/S2."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adastra.archive import write_json, sha256_file, utc_now, ArchiveClient
from adastra.tglc import read_tglc, series, FLUX_COLUMNS
from qlp_variability_followup import fit

DATA=Path('data/tglc-201308024'); OUT=Path('runs/tglc-201308024')
PERIODS=(2.03270208338902,4.194743460281676)
def main():
    if OUT.exists():raise ValueError('Preserving prior declared run')
    declaration=json.loads((DATA/'declaration.json').read_text())
    if declaration['periods_days']!=list(PERIODS) or declaration['allowed_sectors']!=[1,2]:raise ValueError('Declaration mismatch')
    sources=[Path(__file__),Path('src/adastra/tglc.py'),Path('scripts/qlp_variability_followup.py'),Path('src/adastra/archive.py')]
    record={'created_utc':utc_now(),'source_sha256':{str(p):sha256_file(p) for p in sources},'download_declaration_sha256':sha256_file(DATA/'declaration.json'),'method':'Unweighted fixed-frequency sinusoid and independent >0.5day-gap offset/slope; ordinary and Huber2.5; no clipping, both flag fields zero. Calibrated curves already filtered by producer; no new search.','periods_days':list(PERIODS),'protected_holdouts_accessed':False}
    write_json(OUT/'declaration.json',record); rows=[]
    for product in declaration['products']:
        sector=product['sector']; path=ArchiveClient(DATA/'cache').download(product['product']['dataURI'])
        if sha256_file(product['metadata_path'])!=product['metadata_sha256']:raise ValueError('Metadata changed')
        columns, metadata=read_tglc(path,201308024,4906011953900670976,sector)
        fig,axs=plt.subplots(5,3,figsize=(16,14))
        row={'sector':sector,'input_path':str(path),'input_sha256':sha256_file(path),'metadata':metadata,'total_cadences':len(columns['time']),'both_flags_zero':int(np.sum((columns['TESS_flags']==0)&(columns['TGLC_flags']==0))),'columns':{}}
        for i,name in enumerate(FLUX_COLUMNS+('background',)):
            t,y,unit=series(columns,name); values=[]
            axs[i,0].plot(t,y,'.',ms=2);axs[i,0].set_ylabel(name+' '+unit)
            for j,period in enumerate(PERIODS):
                ordinary,nuisance=fit(t,y,period);robust,_=fit(t,y,period,robust=True)
                values.append({'period_days':period,'ordinary':ordinary,'robust':robust})
                phase=((t-1340)/period)%1; clean=y-nuisance
                axs[i,j+1].plot(phase,clean,'.',ms=2);axs[i,j+1].set_title(f'P={period:.6f} d; A={ordinary["amplitude"]:.1f} {unit}')
            row['columns'][name]={'unit':unit,'fixed_period_fits':values}
        fig.suptitle(f'TIC201308024 TGLC S{sector}: same TESS data, Gaia-informed reduction');fig.tight_layout();fig.savefig(OUT/f'sector-{sector}.png',dpi=140);plt.close(fig)
        rows.append(row);write_json(OUT/'partial-results.json',rows)
    write_json(OUT/'results.json',{'completed_utc':utc_now(),'declaration':record,'results':rows,'sources_unchanged':all(sha256_file(p)==record['source_sha256'][str(p)] for p in sources),'claim_limits':'No detection FAP, physical interpretation, novel object, or independent-instrument validation. The raw and calibrated reductions have different filtering transfer functions.'})
if __name__=='__main__':main()
