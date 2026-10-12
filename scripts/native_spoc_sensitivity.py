#!/usr/bin/env python3
"""Declared postselection Huber and day-block sensitivity, no new searches."""
import json
from pathlib import Path
import numpy as np
from astropy.wcs import WCS
from adastra.archive import validated_fits,sha256_file,write_json,utc_now
from native_spoc_201308023 import native_time,ROOT
OUT=ROOT/'runs/native-spoc-201308023-sensitivity';P=4.194743460281676;T0=2458325.

def regression(t,y,e,controls,robust=False):
    scale=np.median(y);y=y/scale-1;e=e/scale
    group=np.cumsum(np.r_[False,np.diff(t)>.5]);columns=[]
    for g in np.unique(group):
        m=group==g;columns += [m.astype(float),np.where(m,t-np.mean(t[m]),0)]
    for c in controls:
        sd=np.std(c)
        if sd>0:columns.append((c-np.mean(c))/sd)
    angle=2*np.pi*(t-T0)/P;columns += [np.sin(angle),np.cos(angle)]
    X=np.column_stack(columns);weights=1/e**2;rw=np.ones(len(t))
    converged=not robust
    for iteration in range(100 if robust else 1):
        a=X*np.sqrt(weights*rw)[:,None];b=y*np.sqrt(weights*rw)
        coef,_,rank,_=np.linalg.lstsq(a,b,rcond=None)
        if rank<X.shape[1]:raise ValueError('Rank deficient sensitivity model')
        residual=y-X@coef
        if not robust:break
        # Huber residuals standardized by own photon errors; robust scale MAD
        standardized=residual/e;median=np.median(standardized)
        scatter=1.4826*np.median(np.abs(standardized-median))
        if not np.isfinite(scatter) or scatter<=0:raise ValueError('Degenerate robust residual scale')
        z=np.abs((standardized-median)/scatter);new=np.minimum(1,2.5/np.maximum(z,1e-300))
        if np.max(np.abs(new-rw))<1e-6:converged=True;rw=new;break
        rw=new
    return {'n':len(t),'semiamplitude_ppm':float(np.hypot(*coef[-2:])*1e6),'phase_cycles':float(np.arctan2(coef[-1],coef[-2])/(2*np.pi)),'iterations':iteration+1,'converged':converged,'downweighted_fraction':float(np.mean(rw<1))}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'declaration.json').exists():raise ValueError('Existing sensitivity run')
    retrieval=ROOT/'runs/native-spoc-201308023-v2/retrieval.json';rec=json.loads(retrieval.read_text())['products']
    inputs={r['path']:sha256_file(ROOT/r['path']) for r in rec if r['sector']==2}
    declaration={'created_utc':utc_now(),'source_sha256':sha256_file(Path(__file__)),'reader_source_sha256':sha256_file(ROOT/'scripts/native_spoc_201308023.py'),'retrieval_sha256':sha256_file(retrieval),'input_sha256':inputs,'sector':2,'period_days':P,'epoch_bjd_tdb':T0,'method':'Primary2x2/3x3 joint neighbor3x3,border median,POS_CORR1/2 plus gap>.5d offsets/slopes. Native neighbor SAP/PDCSAP jointSAP_BKG,MOM_CENTR1/2,POS_CORR1/2. Quality0 finite ownpositiveerrors no clipping. WLS and Huber2.5 IRLS standardized ownerrors with residual MADscale,100iterations,weight convergence1e-6. Fullfit plus deleting each floor(BTJD)1-day bin and refitting all nuisance coefficients and gap groups. No periodoptimization; postselection sensitivity, noFAP or independentvalidation.'}
    write_json(OUT/'declaration.json',declaration);results=[];failures=[]
    for r in rec:
      if r['sector']!=2:continue
      try:
       with validated_fits(ROOT/r['path']) as h:
        t=native_time(h,201308023,2);d=h[1].data;q=np.asarray(d['QUALITY'])
        if r['kind']=='LC':
            controls=[np.asarray(d[k],float) for k in ('SAP_BKG','MOM_CENTR1','MOM_CENTR2','POS_CORR1','POS_CORR2')]
            series=[(name,np.asarray(d[name],float),np.asarray(d[name+'_ERR'],float),controls) for name in ('SAP_FLUX','PDCSAP_FLUX')]
        else:
            f=np.asarray(d['FLUX'],float);e=np.asarray(d['FLUX_ERR'],float);ny,nx=f.shape[1:];yy,xx=np.mgrid[:ny,:nx]
            px,py=r['primary_pixel_xy'];nxp,nyp=WCS(h[2].header).world_to_pixel_values(3.84844604196237,-59.6952750760412)
            p3=(abs(xx-round(px))<=1)&(abs(yy-round(py))<=1);p2=(abs(xx-px)<1)&(abs(yy-py)<1)
            neighbor=(abs(xx-round(float(nxp)))<=1)&(abs(yy-round(float(nyp)))<=1)
            border=((xx==0)|(xx==nx-1)|(yy==0)|(yy==ny-1))&~(p3|neighbor)
            controls=[np.sum(f[:,neighbor],axis=1),np.median(f[:,border],axis=1),np.asarray(d['POS_CORR1'],float),np.asarray(d['POS_CORR2'],float)]
            series=[(name,np.sum(f[:,mask],axis=1),np.sqrt(np.sum(e[:,mask]**2,axis=1)),controls) for name,mask in [('primary2x2',p2),('primary3x3',p3)]]
        for label,flux,error,controls in series:
            good=(q==0)&np.isfinite(t)&np.isfinite(flux)&np.isfinite(error)&(error>0)
            for c in controls:good &= np.isfinite(c)
            tt,yy0,ee=t[good],flux[good],error[good];cc=[c[good] for c in controls];order=np.argsort(tt)
            tt,yy0,ee=tt[order],yy0[order],ee[order];cc=[c[order] for c in cc]
            bins=np.floor(tt-2457000).astype(int);full={};blocks=[]
            for robust in (False,True):
                mode='huber2.5' if robust else 'wls';full[mode]=regression(tt,yy0,ee,cc,robust)
                for day in np.unique(bins):
                    m=bins!=day;value=regression(tt[m],yy0[m],ee[m],[c[m] for c in cc],robust)
                    value.update(omitted_btjd_day=int(day),omitted_n=int(np.sum(~m)),mode=mode)
                    # unwrap relative full-fit phase for interpretable ranges
                    value['phase_difference_cycles']=float((value['phase_cycles']-full[mode]['phase_cycles']+.5)%1-.5)
                    blocks.append(value)
            entry={'series':label,'full':full,'leave_one_day':blocks,'day_bins':len(np.unique(bins))}
            entry['summary']={mode:{'amplitude_min_ppm':min(x['semiamplitude_ppm'] for x in blocks if x['mode']==mode),'amplitude_max_ppm':max(x['semiamplitude_ppm'] for x in blocks if x['mode']==mode),'phase_difference_min_cycles':min(x['phase_difference_cycles'] for x in blocks if x['mode']==mode),'phase_difference_max_cycles':max(x['phase_difference_cycles'] for x in blocks if x['mode']==mode),'all_converged':all(x['converged'] for x in blocks if x['mode']==mode)} for mode in full}
            results.append(entry);print(label,entry['full'],entry['summary'],flush=True)
      except Exception as exc:failures.append({'uri':r['uri'],'error':f'{type(exc).__name__}: {exc}'})
      write_json(OUT/'results.json',{'created_utc':utc_now(),'results':results,'failures':failures,'inputs_unchanged':all(sha256_file(ROOT/p)==s for p,s in inputs.items())})
if __name__=='__main__':main()
