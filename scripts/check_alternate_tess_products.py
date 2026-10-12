#!/usr/bin/env python3
"""Metadata-only, sector-1/2 alternate TESS product availability audit."""
import json
from pathlib import Path
from select_qlp_pilot import mast_pages, write_json, sha, utc

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/alternate-tess-201308024'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    declaration={'created_utc':utc(),'script_sha256':sha(Path(__file__)),'targets':[201308024,201308023],'sectors':[1,2],'flux_accessed':False,'purpose':'Alternate official TESS reductions and target-pixel metadata only; later sectors excluded from requests.'}
    write_json(OUT/'declaration.json',declaration)
    results=[]
    common=[{'paramName':'obs_collection','values':['TESS']},{'paramName':'sequence_number','values':['1','2']}]
    for tic in declaration['targets']:
        query={'service':'Mast.Caom.Filtered','params':{'columns':'*','filters':common+[{'paramName':'target_name','values':[str(tic),'TIC '+str(tic)]}]}}
        try:
            rows=mast_pages(query,OUT/str(tic)/'observations',page_size=1000)
            write_json(OUT/str(tic)/'observation-summary.json',rows)
            products=[]
            for row in rows:
                if int(row.get('sequence_number') or -1) not in (1,2):
                    raise ValueError('Unexpected sector in filtered metadata')
                q={'service':'Mast.Caom.Products','params':{'obsid':str(row['obsid'])}}
                p=mast_pages(q,OUT/str(tic)/'products'/str(row['obsid']),page_size=1000)
                products.extend(p)
            write_json(OUT/str(tic)/'products-summary.json',products)
            results.append({'tic':tic,'status':'complete','observation_count':len(rows),'product_count':len(products),'observations':rows,'products':products})
        except Exception as exc:
            results.append({'tic':tic,'status':'failed','error':f'{type(exc).__name__}: {exc}'})
        write_json(OUT/'results.json',{'created_utc':utc(),'declaration_sha256':sha(OUT/'declaration.json'),'flux_accessed':False,'targets':results})
        print(tic,results[-1]['status'],results[-1].get('observation_count'),results[-1].get('product_count'),flush=True)
if __name__=='__main__': main()
