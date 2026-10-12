"""Scoped public binary-catalogue checks; no light curves or orbital fits."""
import io
import json
from pathlib import Path
import urllib.parse
import requests
from adastra.archive import utc_now, write_json, sha256_file
from review_qlp_catalogs import parse_vizier, primitive

OUT = Path('data/binary-prior-201308024')

def main():
    if (OUT / 'query-declaration.json').exists():
        raise ValueError('Preserve completed attempts')
    queries = [
        ('galah-spectrum', {'-source':'J/A+A/638/A145/catalog', 'spectID':'170909003201311'}),
        ('galah-gaia', {'-source':'J/A+A/638/A145/catalog', 'GaiaDR2':'4906011953900670976'}),
        ('sb9-position', {'-source':'B/sb9/main', '-c':'3.80292450856055 -59.7052863880729', '-c.rs':'10'}),
    ]
    write_json(OUT/'query-declaration.json', {'created_utc':utc_now(),
        'source_sha256':sha256_file(__file__), 'queries':queries,
        'metadata_basis':['galah-sb2-ReadMe.txt','sb9-ReadMe.txt'],
        'scope':'Exact GALAH spectrum/Gaia DR2 identity and a 10 arcsec SB9 positional query; nonmatches do not clear binary novelty'})
    rows=[]
    for name,params in queries:
        params.update({'-out.all':'1','-out.max':'100'})
        url='https://vizier.cds.unistra.fr/viz-bin/votable?'+urllib.parse.urlencode(params)
        row={'name':name,'url':url}
        try:
            r=requests.get(url,timeout=60); r.raise_for_status()
            path=OUT/(name+'.xml');path.write_bytes(r.content)
            write_json(path.with_suffix('.provenance.json'), {'retrieved_utc':utc_now(),
                'url':url,'http_status':r.status_code,'sha256':sha256_file(path),'bytes':len(r.content)})
            doc,warnings=parse_vizier(r.content)
            tables=[]
            for t in doc.iter_tables():
                if len(t.array)>=100:raise ValueError('Reached query row cap')
                tables.append({'catalogue':t.name,'rows':[
                    {n:primitive(rec[n]) for n in t.array.dtype.names} for rec in t.array]})
            if params['-source'] not in [t['catalogue'] for t in tables]:
                raise ValueError('Requested table missing; not a zero result')
            row.update(status='complete',tables=tables,metadata_warnings=warnings)
        except Exception as e:
            row.update(status='failed',error=str(e))
        rows.append(row);print(name,row['status'],sum(len(t['rows']) for t in row.get('tables',[])))
    write_json(OUT/'results.json',{'completed_utc':utc_now(),'queries':rows})

if __name__=='__main__':main()
