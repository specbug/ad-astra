"""Postselection fixed-frequency comparison of already opened QLP SAP curves."""
import json
from pathlib import Path
from adastra.archive import utc_now,write_json,sha256_file
from qlp_variability_followup import raw_columns,column_arrays,fit

def main():
    out=Path('runs/sap-field-frequency-201308024')
    if out.exists():raise ValueError('Preserve previous comparison')
    pilot=Path('runs/qlp-overlapfield-pilot/results.json')
    period=4.194743460281676
    write_json(out/'declaration.json',{'created_utc':utc_now(),
        'stage':'Postselection same-frequency field check, not independent validation or a null population',
        'period_days':period,'pilot_sha256':sha256_file(pilot),
        'sources_sha256':{__file__:sha256_file(__file__),'scripts/qlp_variability_followup.py':sha256_file('scripts/qlp_variability_followup.py')},
        'allowed_sectors':[1,2],'protected_holdouts_accessed':False,
        'method':'All 24 QLP SAP curves, per-gap offsets/slopes, ordinary and Huber2.5, no frequency optimization; compare detector groups and phase without calibrated FAP'})
    rows=[]
    for target in json.loads(pilot.read_text())['results']:
        for s in target['sectors']:
            _,c,camera,ccd=raw_columns(s,target['tic'])
            t,y,e,_=column_arrays(c,'SAP_FLUX')
            a,_=fit(t,y,period,e);b,_=fit(t,y,period,e,robust=True)
            rows.append({'tic':target['tic'],'sector':s['sector'],'camera':camera,'ccd':ccd,
                         'input_sha256':s['input_sha256'],'ordinary':a,'robust':b})
    write_json(out/'results.json',{'completed_utc':utc_now(),'period_days':period,'results':rows})
    for sector in [1,2]:
        vals=sorted([r for r in rows if r['sector']==sector],key=lambda r:r['robust']['amplitude'],reverse=True)
        print('S',sector,'top',[(r['tic'],round(r['robust']['amplitude'],1),r['camera'],r['ccd']) for r in vals[:5]],
              'target rank',next(i+1 for i,r in enumerate(vals) if r['tic']==201308024))

if __name__=='__main__':main()
