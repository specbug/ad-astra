import sys,json,io,urllib.parse
from pathlib import Path
from astropy.io import fits
sys.path.insert(0,'scripts')
from review_qlp_catalogs import retrieve
out=Path('data/rotation-prior-review'); rows=[]
for target in json.loads((out/'tars-observation-summary.json').read_text()):
 for obs in target.get('observations',[]):
  assert obs['sequence_number'] in [1,2]
  uri=obs['dataURL'];fn=out/'tars-sector-products'/Path(uri).name
  try:
   b=retrieve('https://mast.stsci.edu/api/v0.1/Download/file?'+urllib.parse.urlencode({'uri':uri}),fn)
   with fits.open(io.BytesIO(b)) as h:
    primary=dict(h[0].header);keys=['PER'+str(i) for i in range(1,6)]+['POW'+str(i) for i in range(1,6)]+['AMPFIT','SNR1','TMIN1','PKWID1','NAXIS2']
    assert int(primary['TICID'])==target['tic'] and int(primary['SECTOR'])==obs['sequence_number']
    r={'tic':target['tic'],'sector':obs['sequence_number'],'uri':uri,'gaia_dr3':primary.get('GAIADR3'),'period_header':{k:h[1].header.get(k) for k in keys},'flux_arrays_inspected':False};rows.append(r);print(r)
  except Exception as e:rows.append({'tic':target['tic'],'sector':obs['sequence_number'],'error':str(e)})
(out/'tars-header-periods.json').write_text(json.dumps(rows,indent=2))
