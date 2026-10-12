#!/usr/bin/env python3
"""Declared robust sensitivity after inspecting first pixel maps; no period search."""
from pathlib import Path
import json
import numpy as np
from astropy.wcs import WCS
from astropy.time import Time
from adastra.archive import write_json,utc_now,sha256_file,validated_fits
from select_qlp_overlapfield import gaia_query
import qlp_sap_pixels_201308024 as base
ROOT=base.ROOT;OUT=ROOT/'runs/qlp-sap-pixels-201308024-sensitivity';DATA=ROOT/'data/qlp-sap-pixels-201308024-sensitivity'
def fit(t,v,e,p,bg=None):
 label=np.cumsum(np.r_[0,np.diff(t)>.5]);cols=[]
 for k in np.unique(label):
  member=label==k;cols.extend([member.astype(float),member*(t-np.mean(t[member]))/max(np.ptp(t[member]),1)])
 if bg is not None:cols.append((bg-np.median(bg))/np.std(bg))
 phase=2*np.pi*(t-t[0])/p;cols.extend(f(k*phase) for k in (1,2,3) for f in (np.cos,np.sin));a=np.array(cols).T
 good=np.isfinite(v)&np.isfinite(e)&(e>0);a,v,e=a[good],v[good],e[good];aw=a/e[:,None];vw=v/e;weight=np.ones(len(v))
 for _ in range(12):
  c=np.linalg.lstsq(aw*np.sqrt(weight[:,None]),vw*np.sqrt(weight),rcond=None)[0]
  res=vw-aw@c;scale=max(1.4826*np.median(abs(res-np.median(res))),1e-10);weight=np.minimum(1,3*scale/np.maximum(abs(res),1e-10))
 return dict(amplitude_electrons_per_s=float(np.linalg.norm(c[-6:-4])),harmonics=c[-6:].reshape(3,2).tolist(),downweighted_samples=int(np.sum(weight<1)),iterations=12,huber_threshold_mad_sigma=3)
def main():
 OUT.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
 if not (OUT/'declaration.json').exists():write_json(OUT/'declaration.json',dict(utc=utc_now(),script_sha256=sha256_file(__file__),parent_sha256=sha256_file(base.OUT/'results.json'),plan='Post-inspection sensitivity: unchanged two fixed periods, Huber3 12 iterations on error-normalized residuals, same per-gap slopes and3harmonics. All49pixels plus actual3/5/7square masks both border-subtraction/nuisance. Expanded180arcsec Gaia PM cone+independent count to identify bright source wings seen at edge. No new flux or period optimization.',reason='S1 majoroutliers near1348 andupperleftbrightstar beyondoriginal90arcsec possible. S2 targetpeaked4.1947d modulation merits robustness.'))
 condition=f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{base.RA},{base.DEC},0.05))"
 sources=gaia_query('SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,phot_g_mean_mag,phot_variable_flag FROM gaiadr3.gaia_source WHERE '+condition,DATA/'gaia180.json')
 count=gaia_query('SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE '+condition,DATA/'gaia180-count.json')[0]['n']
 if len(sources)!=count:raise ValueError('Incomplete Gaia cone')
 rows=[]
 for original in json.loads((base.OUT/'results.json').read_text())['results']:
  sector=original['sector'];path=ROOT/original['product']['path']
  with validated_fits(path) as h:
   tab=h[1].data;good=np.isfinite(tab['TIME'])&(tab['QUALITY']==0);t=np.array(tab['TIME'][good],float);cube=np.array(tab['FLUX'][good],float);e=np.array(tab['FLUX_ERR'][good],float);wcs=WCS(h[2].header)
  order=np.argsort(t);t,cube,e=t[order],cube[order],e[order];year=float(Time(2457000+np.median(t),format='jd',scale='tdb').jyear);positions=[]
  for s in sources:
   dt=year-s['ref_epoch'];ra=s['ra']+(s.get('pmra') or 0)*dt/3.6e6/np.cos(np.deg2rad(s['dec']));dec=s['dec']+(s.get('pmdec') or 0)*dt/3.6e6;x,y=wcs.world_to_pixel_values(ra,dec)
   positions.append(dict(**s,x=float(x),y=float(y)))
  target=next(s for s in positions if s['source_id']==base.SOURCE);cx,cy=round(target['x']),round(target['y']);y,x=np.indices((7,7));border=np.zeros((7,7),bool);border[[0,-1],:]=True;border[:,[0,-1]]=True;bg=np.nanmedian(cube[:,border],axis=1);apertures=[];pixels=[]
  for p in base.PERIODS:
   for mode,background in [('baseline_only',None),('border_nuisance',bg)]:
    for yy in range(7):
     for xx in range(7):pixels.append(dict(period_days=p,mode=mode,x=xx,y=yy,**fit(t,cube[:,yy,xx],e[:,yy,xx],p,background)))
   for n in (3,5,7):
    mask=(abs(x-cx)<=n//2)&(abs(y-cy)<=n//2);v=np.sum(cube[:,mask],axis=1);sigma=np.sqrt(np.sum(e[:,mask]**2,axis=1));number=int(mask.sum())
    for mode,value,background in [('raw',v,None),('border_subtracted',v-number*bg,None),('border_nuisance',v,bg)]:
     result=fit(t,value,sigma,p,background);result['amplitude_ppm']=result['amplitude_electrons_per_s']/np.median(value)*1e6;apertures.append(dict(period_days=p,aperture=n,pixel_count=number,mode=mode,**result))
  rows.append(dict(sector=sector,gaia_positions=positions,pixels=pixels,apertures=apertures))
  write_json(OUT/'results.json',dict(utc=utc_now(),source_sha256=sha256_file(__file__),gaia_count=count,results=rows))
if __name__=='__main__':main()
