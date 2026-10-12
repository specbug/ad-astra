#!/usr/bin/env python3
"""Declared native SPOC neighbor audit at two previously selected periods."""
import json, sys, re
from pathlib import Path
import numpy as np
from astropy.wcs import WCS
from adastra.archive import ArchiveClient, validated_fits, write_json, sha256_file, utc_now
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/native-spoc-201308023'; OUT=ROOT/'runs/native-spoc-201308023-v2'
PERIODS=[2.03270208338902,4.194743460281676]
PRIMARY=(3.80292450856055,-59.7052863880729)

def native_time(hdus,tic,sector):
    """Fail closed on identity, seconds/days and TDB-reference ambiguity."""
    h=hdus[0].header; ext=hdus[1].header
    if int(h.get('TICID',-1))!=tic or int(h.get('SECTOR',-1))!=sector: raise ValueError('Native TIC/sector identity mismatch')
    if ext.get('TIMESYS')!='TDB' or ext.get('TIMEUNIT')!='d': raise ValueError('Native TIME must be TDB days')
    if 'BJDREFI' not in ext or 'BJDREFF' not in ext: raise ValueError('Missing native BJD reference')
    for k in ('BJDREFI','BJDREFF','TIMESYS','TIMEUNIT'):
        if k in h and h[k]!=ext[k]: raise ValueError('Conflicting native time metadata')
    idx=list(hdus[1].columns.names).index('TIME')
    unit=hdus[1].columns[idx].unit
    if unit not in ('d','day','days'):
        match=re.fullmatch(r'BJD - (\d+(?:\.\d+)?), days',str(unit))
        if not match or float(match[1])!=float(ext['BJDREFI'])+float(ext['BJDREFF']): raise ValueError('Native TIME column unit/reference conflict')
    return np.asarray(hdus[1].data['TIME'],float)+float(ext['BJDREFI'])+float(ext['BJDREFF'])

def main():
    DATA.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
    metadata=ROOT/'data/alternate-tess-201308024/201308023/products-summary.json'
    products=json.loads(metadata.read_text())
    declaration={'created_utc':utc_now(),'source_sha256':sha256_file(Path(__file__)),'metadata_sha256':sha256_file(metadata),'tic':201308023,'sectors':[1,2],'periods_days':PERIODS,'primary_coordinates_icrs':PRIMARY,'method':'Only public native SPOC LC/TP products; original FITS and publisher checksums; native TICID, SECTOR, TIMESYS TDB, day units and BJDREF independently validated. Fixed two-period sinusoid fits to quality-zero SAP and PDCSAP with respective own errors; gap>.5d segment offset and slope. Optional own position/background covariates, no clipping/search/FAP. Target-pixel WCS primary inclusion checked before any primary extraction. Formal regression errors omit correlated noise.','limitations':'Same TESS instrument; independent reduction/cadence rather than independent validation. WCS inclusion is geometric, not source attribution.'}
    if (OUT/'declaration.json').exists(): raise ValueError('Existing run; refusing overwrite')
    write_json(OUT/'declaration.json',declaration)
    result={'created_utc':utc_now(),'products':[],'failures':[]}
    client=ArchiveClient(DATA/'cache')
    for p in products:
        try:
            uri=p['dataURI']; m=re.search(r'-s(\d{4})-0*(\d+)-',uri)
            if not m or int(m[1]) not in (1,2) or int(m[2])!=201308023 or p['project']!='SPOC' or p['dataRights']!='PUBLIC': raise ValueError('Product metadata identity guard')
            sector=int(m[1]); path=client.download(uri)
            with validated_fits(path) as hdus:
                time=native_time(hdus,201308023,sector)
                record={'uri':uri,'sector':sector,'kind':p['productSubGroupDescription'],'path':str(path.relative_to(ROOT)),'sha256':sha256_file(path),'rows':len(time),'columns':list(hdus[1].columns.names),'primary_header':dict(hdus[0].header),'time_header':{k:hdus[1].header.get(k) for k in ('TIMESYS','TIMEUNIT','BJDREFI','BJDREFF','TIMEDEL','TSTART','TSTOP')}}
                if record['kind']=='TP':
                    wcs=WCS(hdus[2].header); x,y=wcs.world_to_pixel_values(*PRIMARY)
                    shape=hdus[1].data['FLUX'].shape[1:]; ny,nx=shape
                    record.update(pixel_shape=list(shape),primary_pixel_xy=[float(x),float(y)],primary_center_inside=bool(-.5<=x<nx-.5 and -.5<=y<ny-.5),aperture_header=dict(hdus[2].header))
                result['products'].append(record)
                print(sector,record['kind'],record.get('primary_pixel_xy'),record.get('primary_center_inside'),flush=True)
        except Exception as exc:
            result['failures'].append({'uri':p.get('dataURI'),'error':f'{type(exc).__name__}: {exc}'})
            print('FAILED',result['failures'][-1],flush=True)
        write_json(OUT/'retrieval.json',result)
if __name__=='__main__': main()
