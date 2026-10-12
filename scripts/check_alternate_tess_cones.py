#!/usr/bin/env python3
"""Supplement exact native-TESS queries with all-provider IDs and narrow cones."""
from pathlib import Path
from select_qlp_pilot import mast_pages, write_json, sha, utc
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data/alternate-tess-201308024'
coords={201308024:(3.80292450856055,-59.7052863880729),201308023:(3.84844604196237,-59.6952750760412)}
def main():
    write_json(OUT/'supplement-declaration.json',{'created_utc':utc(),'script_sha256':sha(Path(__file__)),'sectors':[1,2],'cone_radius_arcsec':5,'flux_accessed':False,'coordinates':coords})
    results=[]
    for tic,(ra,dec) in coords.items():
        for mode in ('all-provider-id','cone'):
            filters=[{'paramName':'sequence_number','values':['1','2']}]
            params={'columns':'*','filters':filters}
            service='Mast.Caom.Filtered'
            if mode=='all-provider-id': filters.append({'paramName':'target_name','values':[str(tic),'TIC '+str(tic)]})
            else:
                service='Mast.Caom.Filtered.Position'; params.update(ra=ra,dec=dec,radius=5/3600)
            try:
                rows=mast_pages({'service':service,'params':params},OUT/str(tic)/mode,page_size=1000)
                for r in rows:
                    if int(r.get('sequence_number') or -1) not in (1,2): raise ValueError('Sector guard failed')
                products=[]
                for r in rows:
                    products.extend(mast_pages({'service':'Mast.Caom.Products','params':{'obsid':str(r['obsid'])}},OUT/str(tic)/'all-products'/str(r['obsid']),page_size=1000))
                results.append({'tic':tic,'query':mode,'status':'complete','observations':rows,'products':products})
            except Exception as e: results.append({'tic':tic,'query':mode,'status':'failed','error':str(e)})
            write_json(OUT/'supplement-results.json',{'created_utc':utc(),'flux_accessed':False,'queries':results})
            print(tic,mode,results[-1]['status'],len(results[-1].get('observations',[])),flush=True)
if __name__=='__main__': main()
