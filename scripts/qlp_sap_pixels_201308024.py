#!/usr/bin/env python3
"""Declared alternative FFI extraction; literal S1/S2, fixed SAP periods."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import numpy as np
from astropy.wcs import WCS
from astropy.time import Time
from adastra.archive import write_json, utc_now, sha256_file, validated_fits
import qlp_rrlyrae_pixels as helper
from select_qlp_overlapfield import gaia_query
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/qlp-sap-pixels-201308024'
OUT=ROOT/'runs/qlp-sap-pixels-201308024'
RA,DEC=3.80292450856055,-59.7052863880729
PERIODS=[2.032702083,4.194743460]
SOURCE=4906011953900670976
helper.DATA,helper.OUT,helper.RA,helper.DEC=DATA,OUT,RA,DEC

def analyse(product,sector,sources):
 import matplotlib
 matplotlib.use('Agg')
 from matplotlib import pyplot as plt
 path=ROOT/product['path']
 if sha256_file(path)!=product['sha256']: raise ValueError('Changed FITS')
 with validated_fits(path) as h:
  tab=h[1]
  if tab.header.get('TIMESYS')!='TDB' or tab.header.get('BJDREFI')!=2457000 or tab.header.get('BJDREFF',0)!=0 or tab.header.get('TIMEUNIT') not in ('d','day','days'): raise ValueError('Timing declaration mismatch')
  time=np.array(tab.data['TIME'],float);cube=np.array(tab.data['FLUX'],float);err=np.array(tab.data['FLUX_ERR'],float)
  good=np.isfinite(time)&(tab.data['QUALITY']==0)
  time,cube,err=time[good],cube[good],err[good]
  order=np.argsort(time);time,cube,err=time[order],cube[order],err[order]
  if np.any(np.diff(time)<=0) or cube.shape[1:]!=(7,7): raise ValueError('Unexpected dimensions')
  wcs=WCS(h[2].header)
 positions=[];year=float(Time(2457000+np.median(time),format='jd',scale='tdb').jyear)
 for s in sources:
  dt=year-s['ref_epoch'];ra=s['ra']+(s.get('pmra') or 0)*dt/3.6e6/np.cos(np.deg2rad(s['dec']));dec=s['dec']+(s.get('pmdec') or 0)*dt/3.6e6
  x,y=wcs.world_to_pixel_values(ra,dec)
  positions.append(dict(source_id=s['source_id'],x=float(x),y=float(y),gmag=s['phot_g_mean_mag'],epoch_jyear=year,proper_motion_available=s.get('pmra') is not None and s.get('pmdec') is not None))
 target=next(s for s in positions if s['source_id']==SOURCE)
 border=np.zeros((7,7),bool);border[[0,-1],:]=True;border[:,[0,-1]]=True
 bg=np.nanmedian(cube[:,border],axis=1)
 # Fixed square apertures centered on nearest proper-motion-adjusted target pixel.
 cx,cy=int(round(target['x'])),int(round(target['y']))
 y,x=np.indices((7,7));apertures={str(n):(abs(x-cx)<=n//2)&(abs(y-cy)<=n//2) for n in (3,5,7)}
 models=[];aperture_results=[]
 for period in PERIODS:
  for mode,background in [('baseline_only',None),('border_background_nuisance',bg)]:
   result=helper.fit_pixels(time,cube,err,background,period)
   pixels=[p for p in result['pixels'] if p['status']=='ok']
   peak=max(pixels,key=lambda p:p['fundamental_semiamplitude_electrons_per_s'])
   result.update(period_days=period,mode=mode,maximum_pixel=peak,nearest_target_pixel=min(pixels,key=lambda p:(p['x']-target['x'])**2+(p['y']-target['y'])**2))
   models.append(result)
  for label,mask in apertures.items():
   raw=np.sum(cube[:,mask],axis=1);sigma=np.sqrt(np.sum(err[:,mask]**2,axis=1));count=int(mask.sum())
   for mode,values,background in [('raw',raw,None),('border_subtracted',raw-count*bg,None),('border_nuisance',raw,bg)]:
    result=helper.fit_pixels(time,values[:,None,None],sigma[:,None,None],background,period)['pixels'][0]
    median=float(np.median(values));result.update(aperture=label,pixel_count=count,mode=mode,period_days=period,median_flux_electrons_per_s=median,fundamental_semiamplitude_ppm=result['fundamental_semiamplitude_electrons_per_s']/median*1e6)
    aperture_results.append(result)
 # Border, flux centroids and total flux temporal diagnostics at both periods.
 positive=np.maximum(cube-bg[:,None,None],0);total=np.sum(positive,axis=(1,2))
 centroid_x=np.sum(positive*x,axis=(1,2))/total;centroid_y=np.sum(positive*y,axis=(1,2))/total
 ancillary=[]
 for label,values in [('border_median',bg),('positive_flux_centroid_x',centroid_x),('positive_flux_centroid_y',centroid_y)]:
  for period in PERIODS:
   fit=helper.fit_pixels(time,values[:,None,None],np.ones((len(time),1,1)),None,period)['pixels'][0]
   ancillary.append(dict(quantity=label,period_days=period,fit=fit))
 fig,axes=plt.subplots(2,3,figsize=(13,8))
 axes[0,0].imshow(np.nanmedian(cube,axis=0),origin='lower');axes[0,0].set_title('Median FFI flux')
 for ax,model in zip(axes.flat[1:5],models):
  image=np.full((7,7),np.nan)
  for p in model['pixels']:
   if p['status']=='ok':image[p['y'],p['x']]=p['fundamental_semiamplitude_electrons_per_s']
  im=ax.imshow(image,origin='lower');fig.colorbar(im,ax=ax);ax.set_title(f"{model['period_days']:.5f} d {model['mode']}")
  for s in positions:
   if -.5<s['x']<6.5 and -.5<s['y']<6.5:ax.plot(s['x'],s['y'],'r+' if s['source_id']==SOURCE else 'wx')
 axes[1,2].plot(time,bg,'.',ms=1);axes[1,2].set_title('Outer-border median')
 fig.tight_layout();fig.savefig(OUT/f's{sector:04d}-maps.png',dpi=140);plt.close(fig)
 fig,axes=plt.subplots(3,1,figsize=(12,8),sharex=True)
 for ax,(label,mask) in zip(axes,apertures.items()):
  values=np.sum(cube[:,mask],axis=1)-int(mask.sum())*bg
  ax.plot(time,(values/np.median(values)-1)*1e6,'.',ms=1);ax.set_ylabel(label+'x'+label+' ppm')
 axes[-1].set_xlabel('BTJD');fig.tight_layout();fig.savefig(OUT/f's{sector:04d}-apertures.png',dpi=140);plt.close(fig)
 return dict(sector=sector,status='ok',samples=len(time),time_range_btjd=[float(time[0]),float(time[-1])],product=product,gaia_positions=positions,models=models,apertures=aperture_results,ancillary=ancillary)

def main():
 DATA.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
 decl=OUT/'declaration.json'
 if not decl.exists():write_json(decl,dict(declared_utc=utc_now(),tic=201308024,sectors=[1,2],periods_days=PERIODS,script_sha256=sha256_file(__file__),helper_sha256=sha256_file(helper.__file__),plan='Quality zero, no clipping. Fixed SAP-selected periods in BOTH sectors. Three harmonic pixel fits, per-gap offset+slope, baseline-only and outer border nuisance. Fixed target-centered 3/5/7 squares raw/border-subtracted/nuisance. Proper-motion Gaia90 arcsec cone. Border/positive flux centroids. No period search or holdout data.',limitations='Same FFI input as QLP; alternative processing, not independent instrument. Border stellar wings, unresolved blends, covariance and correlated noise prevent white-error localization confidence. Fixed periods chosen post-selection.'))
 try:
  query=f"SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.025))"
  sources=gaia_query(query,DATA/'gaia-pm.json')
  rows=[]
  with ThreadPoolExecutor(max_workers=2) as pool:
   jobs={s:pool.submit(helper.fetch,s) for s in (1,2)}
   for s,job in jobs.items():
    try:
     for product in job.result():rows.append(analyse(product,s,sources))
    except Exception as e:rows.append(dict(sector=s,status='failed',error_type=type(e).__name__,error=str(e)))
    write_json(OUT/'results.json',dict(completed_utc=utc_now(),source_sha256=sha256_file(__file__),gaia_sha256=sha256_file(DATA/'gaia-pm.json'),results=rows))
    print(json.dumps(dict(sector=s,status=rows[-1]['status'])),flush=True)
 except Exception as e:
  write_json(OUT/'failure.json',dict(utc=utc_now(),error_type=type(e).__name__,error=str(e)));raise
if __name__=='__main__':main()
