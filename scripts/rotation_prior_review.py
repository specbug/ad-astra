#!/usr/bin/env python3
"""Scoped published rotation metadata queries, without opening withheld sectors."""
from pathlib import Path
import json
from concurrent.futures import ThreadPoolExecutor
from select_qlp_pilot import mast_pages
OUT=Path('data/rotation-prior-review')
IDS=[201308024,355767420,281593180,149392671,150104580]
def query(tic):
 q={'service':'Mast.Caom.Filtered','params':{'columns':'obsid,obs_id,target_name,sequence_number,t_min,t_max,provenance_name,dataURL','filters':[{'paramName':'provenance_name','values':['tars']},{'paramName':'target_name','values':[str(tic)]},{'paramName':'sequence_number','values':['1','2']}]}}
 try:return {'tic':tic,'status':'complete','observations':mast_pages(q,OUT/'tars-observations'/str(tic),page_size=1000)}
 except Exception as e:return {'tic':tic,'status':'failed','error':str(e)}
if __name__=='__main__':
 with ThreadPoolExecutor(4) as pool:r=list(pool.map(query,IDS))
 (OUT/'tars-observation-summary.json').write_text(json.dumps(r,indent=2))
 for x in r:print(x)
