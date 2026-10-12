#!/usr/bin/env python3
"""Fixed-frequency native SPOC LC/TPF descriptive source-attribution checks."""
import json
from pathlib import Path
import numpy as np
from astropy.wcs import WCS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adastra.archive import validated_fits, write_json, sha256_file, utc_now
from native_spoc_201308023 import native_time, PERIODS, ROOT
OUT=ROOT/'runs/native-spoc-201308023-analysis'
T0=2458325.0

def fit(time,flux,error,quality,period,controls=None):
    controls=[] if controls is None else controls
    good=(quality==0)&np.isfinite(time)&np.isfinite(flux)&np.isfinite(error)&(error>0)
    for c in controls: good &= np.isfinite(c)
    t,y,e=time[good],flux[good],error[good]
    order=np.argsort(t); t,y,e=t[order],y[order],e[order]
    scale=np.median(y)
    if len(t)<100 or not np.isfinite(scale) or scale<=0: raise ValueError('Insufficient positive finite native flux')
    y=y/scale-1;e=e/scale
    group=np.cumsum(np.r_[False,np.diff(t)>.5]);columns=[]
    for g in np.unique(group):
        mask=group==g;columns += [mask.astype(float),np.where(mask,t-np.mean(t[mask]),0)]
    for c in controls:
        v=np.asarray(c)[good][order]; sd=np.std(v)
        if sd>0: columns.append((v-np.mean(v))/sd)
    phase=2*np.pi*(t-T0)/period
    columns += [np.sin(phase),np.cos(phase)]
    X=np.column_stack(columns); A=X/e[:,None]; b=y/e
    beta,_,rank,_=np.linalg.lstsq(A,b,rcond=None)
    if rank<X.shape[1]:raise ValueError('Rank deficient descriptive fit')
    residual=y-X@beta; dof=len(y)-len(beta)
    covariance=np.linalg.inv(A.T@A)*np.sum((residual/e)**2)/dof
    amp=np.hypot(*beta[-2:]); gradient=beta[-2:]/amp if amp else np.zeros(2)
    amp_error=float(np.sqrt(gradient@covariance[-2:,-2:]@gradient))
    return {'n':len(t),'period_days':period,'epoch_bjd_tdb':T0,'semiamplitude_ppm':float(amp*1e6),'formal_residual_scaled_white_error_ppm':amp_error*1e6,'phase_cycles':float(np.arctan2(beta[-1],beta[-2])/(2*np.pi)),'median_flux_e_per_s':float(scale),'residual_rms_ppm':float(np.std(residual)*1e6),'control_count':len(controls),'segments':int(group.max()+1)},(t,y,X[:,:-2]@beta[:-2],beta[-2:])

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    retrieval=ROOT/'runs/native-spoc-201308023-v2/retrieval.json'; records=json.loads(retrieval.read_text())['products']
    declaration={'created_utc':utc_now(),'source_sha256':sha256_file(Path(__file__)),'native_reader_sha256':sha256_file(ROOT/'scripts/native_spoc_201308023.py'),'retrieval_sha256':sha256_file(retrieval),'periods_days':PERIODS,'method':'quality=0 finite positive own errors, no clipping; weighted separate fixed-period sin/cos, gap>.5d offsets/slopes. LC SAP/PDCSAP own errors; optional SAP_BKG/MOM_CENTR/POS_CORR controls. TPF flux already pipeline-background-subtracted; 2x2 nearest or3x3 rounded-center primary aperture, quadrature native pixel errors, border median excludes3x3 primary/neighbor masks. Base and border+position controls and joint bright-neighbor-aperture controls. No optimization/FAP. Errors residual-scaled white-noise formal only. No aperture correction, PRF deblend, covariance model or intrinsic-source claim. Time is native bright-neighbor barycentric TDB; negligible differential direction correction is not recomputed.'}
    if (OUT/'declaration.json').exists():raise ValueError('Existing analysis run')
    write_json(OUT/'declaration.json',declaration)
    results=[];failures=[]
    for rec in records:
      try:
       with validated_fits(ROOT/rec['path']) as h:
        t=native_time(h,201308023,rec['sector']); d=h[1].data;q=np.asarray(d['QUALITY'])
        if rec['kind']=='LC':
            controls=[np.asarray(d[k],float) for k in ('SAP_BKG','MOM_CENTR1','MOM_CENTR2','POS_CORR1','POS_CORR2') if k in d.names]
            series=[(field,np.asarray(d[field],float),np.asarray(d[field+'_ERR'],float),{'base':[],'instrument_controls':controls}) for field in ('SAP_FLUX','PDCSAP_FLUX')]
        else:
            f=np.asarray(d['FLUX'],float);e=np.asarray(d['FLUX_ERR'],float);ny,nx=f.shape[1:];yy,xx=np.mgrid[:ny,:nx]
            px,py=rec['primary_pixel_xy'];nxp,nyp=WCS(h[2].header).world_to_pixel_values(3.84844604196237,-59.6952750760412)
            primary3=(abs(xx-round(px))<=1)&(abs(yy-round(py))<=1)
            primary2=(abs(xx-px)<1)&(abs(yy-py)<1)
            neighbor=(abs(xx-round(float(nxp)))<=1)&(abs(yy-round(float(nyp)))<=1)
            border=((xx==0)|(xx==nx-1)|(yy==0)|(yy==ny-1))&~(primary3|neighbor)
            bg=np.median(f[:,border],axis=1); nf=np.sum(f[:,neighbor],axis=1); ne=np.sqrt(np.sum(e[:,neighbor]**2,axis=1))
            pointing=[np.asarray(d[k],float) for k in ('POS_CORR1','POS_CORR2')]
            series=[('neighbor3x3',nf,ne,{'base':[],'background_pointing':[bg,*pointing]})]
            for label,mask in [('primary2x2',primary2),('primary3x3',primary3)]:
                ff=np.sum(f[:,mask],axis=1);ee=np.sqrt(np.sum(e[:,mask]**2,axis=1))
                if int(mask.sum())!=(4 if label.endswith('2x2') else 9):raise ValueError('Primary aperture clipped by stamp edge')
                series.append((label,ff,ee,{'base':[],'background_pointing':[bg,*pointing],'joint_neighbor_background_pointing':[nf,bg,*pointing]}))
            write_json(OUT/f"s{rec['sector']}-apertures.json",{'primary2x2_indices_yx':np.argwhere(primary2).tolist(),'primary3x3_indices_yx':np.argwhere(primary3).tolist(),'neighbor_indices_yx':np.argwhere(neighbor).tolist(),'border_indices_yx':np.argwhere(border).tolist(),'primary_neighbor_overlap_pixels':int(np.sum(primary3&neighbor))})
        fig,axes=plt.subplots(len(series),len(PERIODS),figsize=(12,3*len(series)),squeeze=False)
        for i,(label,flux,err,modes) in enumerate(series):
            for j,period in enumerate(PERIODS):
                for mode,controls in modes.items():
                    value,plot=fit(t,flux,err,q,period,controls)
                    value.update(sector=rec['sector'],product=rec['kind'],series=label,mode=mode,product_sha256=rec['sha256'])
                    results.append(value)
                    if mode=='base':
                        tt,yy0,baseline,coef=plot;phase=((tt-T0)/period)%1;ax=axes[i,j]
                        ax.scatter(phase,(yy0-baseline)*1e6,s=1,alpha=.15)
                        grid=np.linspace(0,1,200);ax.plot(grid,(coef[0]*np.sin(2*np.pi*grid)+coef[1]*np.cos(2*np.pi*grid))*1e6,color='red')
                        ax.set_title(f'{label}: P={period:.6f}d, A={value["semiamplitude_ppm"]:.0f}ppm');ax.set_xlabel('Fixed-period phase');ax.set_ylabel('Baseline-removed flux (ppm)')
        fig.tight_layout();fig.savefig(OUT/f"s{rec['sector']}-{rec['kind']}-fixed-periods.png",dpi=160);plt.close(fig)
      except Exception as exc:failures.append({'uri':rec['uri'],'error':f'{type(exc).__name__}: {exc}'})
      write_json(OUT/'results.json',{'created_utc':utc_now(),'fits':results,'failures':failures})
      print(rec['sector'],rec['kind'],'fits',len(results),'failures',len(failures),flush=True)
if __name__=='__main__':main()
