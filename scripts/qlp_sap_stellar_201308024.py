import sys,json,urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,'scripts')
import qlp_stellar_metadata as b
import qlp_neighbor_metadata as n
import review_qlp_catalogs as p
from adastra.archive import write_json,utc_now,sha256_file
out=Path('data/qlp-sap-stellar-201308024');out.mkdir(exist_ok=True)
ids=[4906011953900670976,4906017863775668992,4906017863775668864]
write_json(out/'declaration.json',{'created_utc':utc_now(),'purpose':'Post-inspection exact-source Gaia and spectroscopy metadata; no flux/holdouts','source_ids':ids,'tic':201308024,'source_sha256':sha256_file(Path(__file__)),'coordinates_icrs':[3.80292450856,-59.705286388],'periods_exploratory_days':[2.032702083,4.194743460]})
q='SELECT * FROM gaiadr3.gaia_source WHERE source_id IN ('+','.join(map(str,ids))+')'
def gaia():return b.tap(b.GAIA,q,out/'gaia-exact.json')
def galah():
 url='https://vizier.cds.unistra.fr/viz-bin/votable?'+urllib.parse.urlencode({'-source':'J/MNRAS/506/150/stars,J/other/PASA/42.51/allstar','-c':'3.80292450856 -59.705286388','-c.rs':5,'-out.all':'1','-out.max':100})
 data=n.parse_spectroscopy(p.retrieve(url,out/'galah.xml'))
 data['missing_requested_tables']=sorted({'J/MNRAS/506/150/stars','J/other/PASA/42.51/allstar'}-{t['catalog'] for t in data['tables']})
 return data
def simbad():
 q="SELECT b.main_id,b.otype,b.ra,b.dec,ids.ids,alltypes.otypes FROM basic AS b LEFT JOIN ids ON b.oid=ids.oidref LEFT JOIN alltypes ON b.oid=alltypes.oidref WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',3.80292450856,-59.705286388,0.03))"
 return b.tap(b.SIMBAD,q,out/'simbad.json')
with ThreadPoolExecutor(max_workers=3) as e:
 fs={k:e.submit(fn) for k,fn in [('gaia',gaia),('galah',galah),('simbad',simbad)]}
 result={}
 for k,f in fs.items():
  try:result[k]=f.result()
  except Exception as ex:result[k]={'failure':str(ex)}
write_json(out/'summary.json',result)
print(json.dumps(result,indent=2))
