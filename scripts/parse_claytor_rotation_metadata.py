import io,json,warnings
from pathlib import Path
from astropy.io.votable import parse
p=Path('data/rotation-prior-review');results=[]
for f in sorted(p.glob('2024-table*.raw')):
 with warnings.catch_warnings(record=True) as ww:
  warnings.simplefilter('always');d=parse(io.BytesIO(f.read_bytes()),verify='warn')
 reviewed=[]
 for w in ww:
  s=str(w.message)
  if "W03:" in s or "W50: Invalid unit string '(10-6)'" in s or "W50: Invalid unit string 'log(Sun)'" in s:reviewed.append(s)
  else:raise ValueError(s)
 infos=[{'name':i.name,'value':i.value,'content':i.content} for i in d.infos]
 if any(i['value'] in ['ERROR','OVERFLOW'] for i in infos):raise ValueError(infos)
 ts=[]
 for t in d.iter_tables():
  a=t.to_table();rows=[]
  for row in a:
   rr={}
   for c in a.colnames:
    v=row[c];rr[c]=None if getattr(v,'mask',False) else v.item() if hasattr(v,'item') else v
   rows.append(rr)
  ts.append({'table':t.ID,'rows':rows,'fields':a.colnames})
 results.append({'file':f.name,'status':'parsed_metadata_warning_reviewed','warnings':reviewed,'tables':ts})
(p/'claytor-results.json').write_text(json.dumps(results,indent=2,default=str))
for r in results:print(r['file'],[(t['table'],t['rows']) for t in r['tables']])
