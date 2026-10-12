#!/usr/bin/env python3
"""Declared 11x11 source/background sensitivity; no holdouts or period search."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json,zipfile
from urllib.parse import urlparse
import requests
import numpy as np
from astropy.wcs import WCS
from astropy.time import Time
from adastra.archive import sha256_file,utc_now,validated_fits,write_json
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/qlp-sap-pixels-201308024-11x11'
OUT=ROOT/'runs/qlp-sap-pixels-201308024-11x11'
API='https://mast.stsci.edu/tesscut/api/v0.1/astrocut'
RA,DEC=3.80292450856055,-59.7052863880729
MAX_BYTES=128*1024*1024
PERIODS=[2.032702083,4.194743460]
TARGET=4906011953900670976
NEIGHBOR=4906017863775668992
def fetch(sector):
    if type(sector) is not int or sector not in (1, 2):
        raise ValueError("Literal sector 1/2 guard protects the old holdouts")
    directory = DATA / f"s{sector:04d}"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "tesscut.zip"
    sidecar = directory / "provenance.json"
    parameters = dict(ra=RA, dec=DEC, x=11, y=11, units="px", sector=sector)
    if path.exists():
        saved = json.loads(sidecar.read_text())
        if saved["parameters"] != parameters or saved["zip_sha256"] != sha256_file(path):
            raise ValueError("Cached archive mismatch")
    else:
        part = path.with_suffix(".zip.part")
        if part.exists():
            raise ValueError("Preserve existing partial for explicit recovery")
        initial = requests.get(API, params=parameters, stream=True, timeout=(30, 150), allow_redirects=False)
        initial_url = initial.url
        delivery_url = initial.headers.get("Location")
        regional_delivery = False
        if initial.is_redirect:
            location = urlparse(delivery_url)
            if (location.scheme != "https" or location.hostname != "stpubdata.s3.amazonaws.com"
                    or location.query or not location.path.startswith("/tess/public/staged_cutouts/")):
                initial.close()
                raise ValueError("Unexpected or signed TESSCut delivery redirect")
            # The same public bucket's official regional endpoint is already in
            # this environment's NASA-data allowlist. Preserve both exact URLs.
            delivery_url = "https://stpubdata.s3.us-east-1.amazonaws.com" + location.path
            regional_delivery = True
            initial.close()
            response = requests.get(delivery_url, stream=True, timeout=(30, 120))
        else:
            response = initial
        with response:
            record = dict(requested_utc=utc_now(), parameters=parameters,
                          initial_url=initial_url, delivery_url=delivery_url,
                          regional_delivery=regional_delivery,
                          resolved_url=response.url, http_status=response.status_code)
            if not response.ok:
                record["error_body"] = response.text[:4000]
                write_json(directory / "error.json", record)
                response.raise_for_status()
            size = 0
            with part.open("xb") as stream:
                for block in response.iter_content(1024 * 1024):
                    size += len(block)
                    if size > MAX_BYTES:
                        raise ValueError("TESSCut archive exceeds explicit 128 MiB cap")
                    stream.write(block)
            if not zipfile.is_zipfile(part):
                write_json(directory / "error.json", {**record, "error": "Non-ZIP body",
                           "body": part.read_bytes()[:4000].decode(errors="replace")})
                raise ValueError("TESSCut did not return a ZIP archive")
            record.update(completed_utc=utc_now(), zip_bytes=size, zip_sha256=sha256_file(part))
            part.replace(path)
            write_json(sidecar, record)
    products = []
    with zipfile.ZipFile(path) as archive:
        members = [x for x in archive.infolist() if x.filename.lower().endswith(".fits")]
        if not members or sum(x.file_size for x in members) > MAX_BYTES:
            raise ValueError("Unexpected or oversized FITS archive")
        for member in members:
            # Flat basename only: never extract user/server supplied path components.
            target = directory / Path(member.filename).name
            payload = archive.read(member)
            if target.exists() and target.read_bytes() != payload:
                raise ValueError("Existing extracted artifact differs from preserved ZIP")
            if not target.exists():
                target.write_bytes(payload)
            with validated_fits(target) as hdus:
                sectors = {int(h.header["SECTOR"]) for h in hdus if "SECTOR" in h.header}
                if sectors != {sector}:
                    raise ValueError(f"Downloaded sector mismatch: {sectors}")
                integrity = [{"hdu": i, "checksum_supplied": "CHECKSUM" in h.header,
                              "datasum_supplied": "DATASUM" in h.header}
                             for i, h in enumerate(hdus)]
            products.append(dict(path=str(target.relative_to(ROOT)), sha256=sha256_file(target),
                                 original_zip_member=member.filename, integrity=integrity))
    write_json(directory / "products.json", products)
    return products


def fit(t,v,e,p,nuisances):
 label=np.cumsum(np.r_[0,np.diff(t)>.5]);cols=[]
 for k in np.unique(label):
  m=label==k;cols.extend([m.astype(float),m*(t-np.mean(t[m]))/max(np.ptp(t[m]),1)])
 for n in nuisances:
  if np.std(n)>0:cols.append((n-np.median(n))/np.std(n))
 phase=2*np.pi*(t-t[0])/p;cols.extend(f(k*phase) for k in (1,2,3) for f in (np.cos,np.sin));a=np.array(cols).T
 good=np.isfinite(v)&np.isfinite(e)&(e>0)&np.all(np.isfinite(a),axis=1);a,v,e=a[good],v[good],e[good];aw=a/e[:,None];vw=v/e;weight=np.ones(len(v))
 for _ in range(12):
  c,_,rank,singular=np.linalg.lstsq(aw*np.sqrt(weight[:,None]),vw*np.sqrt(weight),rcond=None)
  res=vw-aw@c;scale=max(1.4826*np.median(abs(res-np.median(res))),1e-10);weight=np.minimum(1,3*scale/np.maximum(abs(res),1e-10))
 return dict(amplitude_electrons_per_s=float(np.linalg.norm(c[-6:-4])),harmonics=c[-6:].reshape(3,2).tolist(),rank=int(rank),columns=a.shape[1],condition_number=float(singular[0]/singular[-1]),downweighted_samples=int(np.sum(weight<1)),coefficient_nuisances=c[2*len(np.unique(label)):-6].tolist())
def analyse(product,s,sources):
 import matplotlib
 matplotlib.use('Agg')
 from matplotlib import pyplot as plt
 path=ROOT/product['path']
 if sha256_file(path)!=product['sha256']:raise ValueError('Changed file')
 with validated_fits(path) as h:
  tab=h[1];data=tab.data
  if tab.header.get('TIMESYS')!='TDB' or tab.header.get('BJDREFI')!=2457000 or tab.header.get('BJDREFF',0)!=0 or tab.header.get('TIMEUNIT') not in ('d','days','day'):raise ValueError('Timing mismatch')
  if h[0].header['RA_OBJ']!=RA or h[0].header['DEC_OBJ']!=DEC:raise ValueError('Coordinate mismatch')
  t=np.array(data['TIME'],float);good=np.isfinite(t)&(data['QUALITY']==0);t=t[good];cube=np.array(data['FLUX'][good],float);e=np.array(data['FLUX_ERR'][good],float);wcs=WCS(h[2].header)
  if tab.columns['FLUX'].unit!='e-/s' or tab.columns['FLUX_ERR'].unit!='e-/s':raise ValueError('Flux unit mismatch')
 order=np.argsort(t);t,cube,e=t[order],cube[order],e[order]
 if np.any(np.diff(t)<=0) or cube.shape[1:]!=(11,11):raise ValueError('Dimensions mismatch')
 positions=[];year=float(Time(2457000+np.median(t),format='jd',scale='tdb').jyear)
 for g in sources:
  dt=year-g['ref_epoch'];ra=g['ra']+(g.get('pmra') or 0)*dt/3.6e6/np.cos(np.deg2rad(g['dec']));dec=g['dec']+(g.get('pmdec') or 0)*dt/3.6e6;x,y=wcs.world_to_pixel_values(ra,dec);positions.append(dict(**g,x=float(x),y=float(y)))
 target=next(g for g in positions if g['source_id']==TARGET);neighbor=next(g for g in positions if g['source_id']==NEIGHBOR)
 yy,xx=np.indices((11,11));border=np.zeros((11,11),bool);border[[0,-1],:]=True;border[:,[0,-1]]=True
 clean=border.copy()
 for g in positions:
  if g['phot_g_mean_mag'] is not None and g['phot_g_mean_mag']<17:clean &= np.hypot(xx-g['x'],yy-g['y'])>2.5
 if clean.sum()<5:raise ValueError('Insufficient clean border pixels')
 bgs={'all_border':np.nanmedian(cube[:,border],axis=1),'bright_masked_border':np.nanmedian(cube[:,clean],axis=1)}
 masks={};photometry={}
 for label,g in [('target',target),('neighbor',neighbor)]:
  for n in (3,5):
   mask=(abs(xx-round(g['x']))<=n//2)&(abs(yy-round(g['y']))<=n//2)
   other=neighbor if label=='target' else target
   # Exclude other star's inner 1.5 pixels; recordactualmask andcount.
   mask &= np.hypot(xx-other['x'],yy-other['y'])>1.5
   masks[f'{label}_{n}']=mask;photometry[f'{label}_{n}']=(np.sum(cube[:,mask],axis=1),np.sqrt(np.sum(e[:,mask]**2,axis=1)))
 bg=bgs['bright_masked_border'];positive=np.maximum(cube-bg[:,None,None],0)
 centroid_x=np.sum(positive*xx,axis=(1,2))/np.sum(positive,axis=(1,2));centroid_y=np.sum(positive*yy,axis=(1,2))/np.sum(positive,axis=(1,2))
 apertures=[];pixels=[];joint=[]
 for p in PERIODS:
  for label,(v,sigma) in photometry.items():
   for bglabel,bg in bgs.items():
    for mode,value,nuisances in [('subtracted',v-masks[label].sum()*bg,[]),('nuisance',v,[bg])]:
     f=fit(t,value,sigma,p,nuisances);f['amplitude_ppm']=f['amplitude_electrons_per_s']/np.median(value)*1e6;apertures.append(dict(period_days=p,aperture=label,actual_pixel_count=int(masks[label].sum()),background=bglabel,mode=mode,**f))
  for bglabel,bg in bgs.items():
   for y in range(11):
    for x in range(11):pixels.append(dict(period_days=p,background=bglabel,x=x,y=y,**fit(t,cube[:,y,x],e[:,y,x],p,[bg])))
  for n in (3,5):
   v,sigma=photometry[f'target_{n}'];nv,_=photometry[f'neighbor_{n}'];bg=bgs['bright_masked_border'];neighbor_sub=nv-masks[f'neighbor_{n}'].sum()*bg
   for mode,nuisances in [('background_neighbor',[bg,neighbor_sub]),('background_neighbor_centroids',[bg,neighbor_sub,centroid_x,centroid_y])]:
    f=fit(t,v,sigma,p,nuisances);f['amplitude_ppm']=f['amplitude_electrons_per_s']/np.median(v)*1e6;joint.append(dict(period_days=p,aperture=n,mode=mode,target_neighbor_pearson=float(np.corrcoef(v,neighbor_sub)[0,1]),**f))
 fig,axes=plt.subplots(2,3,figsize=(14,9));axes[0,0].imshow(np.median(cube,axis=0),origin='lower');axes[0,0].set_title('Median flux')
 for ax,(p,bglabel) in zip(axes.flat[1:5],[(p,b) for p in PERIODS for b in bgs]):
  image=np.full((11,11),np.nan)
  for pix in pixels:
   if pix['period_days']==p and pix['background']==bglabel:image[pix['y'],pix['x']]=pix['amplitude_electrons_per_s']
  im=ax.imshow(image,origin='lower');fig.colorbar(im,ax=ax);ax.set_title(f'{p:.5f}d {bglabel}')
  for g in positions:
   if -.5<g['x']<10.5 and -.5<g['y']<10.5:ax.plot(g['x'],g['y'],'r+' if g['source_id']==TARGET else ('co' if g['source_id']==NEIGHBOR else 'wx'))
 for k,bg in bgs.items():axes[1,2].plot(t,bg,'.',ms=1,label=k)
 axes[1,2].legend();fig.tight_layout();fig.savefig(OUT/f's{s:04d}-maps.png',dpi=130);plt.close(fig)
 return dict(sector=s,status='ok',samples=len(t),product=product,gaia_positions=positions,clean_border_pixels=int(clean.sum()),masks={k:v.astype(int).tolist() for k,v in masks.items()},apertures=apertures,pixels=pixels,joint=joint,limitations='No PRF deblending; companionnearbrightneighborunresolved; time/interpixelcovarianceomitted;neighborregressionmayremove real correlatedsignal; finitecutout masks truncatesneighbor5x5 inS2; conditionsreported notconfidence.')
def main():
 DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
 if not (OUT/'declaration.json').exists():write_json(OUT/'declaration.json',dict(utc=utc_now(),source_sha256=sha256_file(__file__),sectors=[1,2],periods=PERIODS,plan='11x11 cutouts; literalS1/S2 only. Quality0 noclip; Huber3 residualMAD12iterations. Fixedperiod3harmonics per-gapoffsetslopes. Target/neighbor3x3/5x5 masks excludingotherinner1.5px. Allborder vsborderexcludeGaiaG<17within2.5px median. Bothsubtracted+nuisance. Jointtargetneighborbackground optionalglobalpositivefluxcentroids withsinesincluded. Rank/condition/correlations. Pixelmapsbothbackgrounds. No PRF orperiodoptimization.',snapshot_inventory_sha256=sha256_file(OUT/'source-snapshot/inventory.json'),gaia_snapshot='data/qlp-sap-pixels-201308024-sensitivity/gaia180.json',gaia_sha256=sha256_file(ROOT/'data/qlp-sap-pixels-201308024-sensitivity/gaia180.json')))
 raw=json.loads((ROOT/'data/qlp-sap-pixels-201308024-sensitivity/gaia180.json').read_text());names=[g['name'] for g in raw['metadata']];sources=[dict(zip(names,row)) for row in raw['data']];rows=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  jobs={s:pool.submit(fetch,s) for s in (1,2)}
  for s,job in jobs.items():
   try:
    for product in job.result():rows.append(analyse(product,s,sources))
   except Exception as e:rows.append(dict(sector=s,status='failed',error_type=type(e).__name__,error=str(e)))
   write_json(OUT/'results.json',dict(utc=utc_now(),source_sha256=sha256_file(__file__),results=rows));print(json.dumps(dict(sector=s,status=rows[-1]['status'])),flush=True)
if __name__=='__main__':main()
