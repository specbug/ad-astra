#!/usr/bin/env python3
"""Exploratory SAP diagnostics on already opened overlap-field S1/S2 only."""
from pathlib import Path
import json
import numpy as np
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adastra.archive import write_json,sha256_file,utc_now
from qlp_variability_followup import raw_columns,column_arrays,fit,design
LEADS=(201308024,355767420,281593180)
OUT=Path('runs/qlp-overlapfield-sap-followup-v3')
PILOT=Path('runs/qlp-overlapfield-pilot/results.json')
def common_mode(t, controls, min_controls=2):
    """Nearest cadence within 0.005 d (7.2 minutes); no interpolation."""
    keys=np.round(t,6); rows=[]
    for ct,cy in controls:
        idx=np.searchsorted(ct,t); hi=np.clip(idx,0,len(ct)-1); lo=np.clip(idx-1,0,len(ct)-1)
        nearest=np.where(np.abs(ct[lo]-t)<=np.abs(ct[hi]-t),lo,hi)
        delta=np.abs(ct[nearest]-t); values=cy[nearest].copy(); values[delta>.005]=np.nan
        rows.append(values)
    if not rows:return np.full(len(t),np.nan),np.zeros(len(t),int)
    matrix=np.asarray(rows); count=np.sum(np.isfinite(matrix),axis=0)
    with np.errstate(all='ignore'): mode=np.nanmedian(matrix,axis=0)
    mode[count<min_controls]=np.nan
    return mode,count

def residualize(t,y,p,mode,robust=False):
    good=np.isfinite(mode); x=design(t[good],p)
    # Joint signal+nuisance fit preserves a target sinusoid unless collinear.
    scale=np.std(mode[good]); z=mode[good]/scale if scale else mode[good]
    xx=np.column_stack([x,z]); weights=np.ones(len(z))
    for iteration in range(30 if robust else 1):
        beta,_,rank,_=np.linalg.lstsq(xx*np.sqrt(weights)[:,None],y[good]*np.sqrt(weights),rcond=None)
        residual=y[good]-xx@beta; center=np.median(residual); scatter=1.4826*np.median(abs(residual-center))
        if scatter==0: break
        new=np.minimum(1,2.5*scatter/np.maximum(abs(residual-center),np.finfo(float).eps))
        if np.max(abs(weights-new))<1e-6: break
        weights=new
    if rank<xx.shape[1]: raise ValueError('Common mode is rank deficient with sinusoid and segment terms')
    clean=y[good]-z*beta[-1]
    return good,clean,{'common_mode_coefficient':float(beta[-1]/scale) if scale else 0.,'rank':int(rank),'design_columns':xx.shape[1],'condition_number':float(np.linalg.cond(xx)),'joint_huber':robust,'iterations':iteration+1}

def main():
    if OUT.exists():raise ValueError('Output exists; preserve previous attempt')
    pilot=json.loads(PILOT.read_text()); assert pilot['status']=='complete'
    declaration={'created_utc':utc_now(),'stage':'Postselection: visual SAP leads; both sectors already viewed','pilot_sha256':sha256_file(PILOT),'source_sha256':sha256_file(__file__),'leads':LEADS,'allowed_sectors':[1,2],'protected_holdouts_accessed':False,'methods':'SAP unweighted LS 0.5 to min(13, baseline/2) days; fixed S1 frequency for S2. Same camera/CCD leave-target-out median SAP, nearest cadence matches within 0.005d to accommodate target-dependent barycentric time; no interpolation, at least 2 controls. Joint sinusoid + per-gap offsets/slopes + fixed common mode ordinary and Huber2.5 joint regression; robust Huber fit of residual sensitivity. No FAP or localization. All24 SAP evaluated at lead S1 frequencies.'}
    write_json(OUT/'declaration.json',declaration)
    data={}
    for target in pilot['results']:
        for row in target['sectors']:
            assert row['sector'] in (1,2)
            _,cols,cam,ccd=raw_columns(row,target['tic']); t,y,_,_=column_arrays(cols,'SAP_FLUX')
            data[(target['tic'],row['sector'])]=(t,y,cols,cam,ccd,row)
    periods={}; search={}
    for tic in LEADS:
        t,y,*_=data[tic,1]; ls=LombScargle(t,y); f,power=ls.autopower(minimum_frequency=1/min(13,(t[-1]-t[0])/2),maximum_frequency=2,samples_per_peak=10)
        k=np.argmax(power); periods[tic]=float(1/f[k]); search[tic]={'period_days':periods[tic],'power':float(power[k]),'period_at_search_upper_bound':bool(k==0)}
    write_json(OUT/'s1-selected-periods.json',search)
    results=[]; comparisons=[]
    for tic in LEADS:
        p=periods[tic]; fig,axs=plt.subplots(2,3,figsize=(15,7))
        for sector in (1,2):
            t,y,cols,cam,ccd,row=data[tic,sector]
            peers=[(pt,py) for (other,ss),(pt,py,pc,pca,pcc,pr) in data.items() if other!=tic and ss==sector and pca==cam and pcc==ccd]
            mode,count=common_mode(t,peers); offsets=[]
            for ct,cy in peers:
                ix=np.searchsorted(ct,t); hi=np.clip(ix,0,len(ct)-1); lo=np.clip(ix-1,0,len(ct)-1); dd=np.minimum(abs(ct[hi]-t),abs(ct[lo]-t)); offsets.extend(dd[dd<=.005].tolist())
            result={'tic':tic,'sector':sector,'camera':cam,'ccd':ccd,'period_days':p,'input_sha256':row['input_sha256'],'same_detector_controls':len(peers),'pairing_max_offset_days':max(offsets) if offsets else None,'pairing_median_offset_days':float(np.median(offsets)) if offsets else None,'minimum_required_controls':2,'cadence_coverage_fraction':float(np.mean(np.isfinite(mode))),'control_count_min':int(count.min()),'control_count_max':int(count.max()),'columns':{}}
            for name in ('SAP_FLUX','KSPSAP_FLUX','KSPSAP_FLUX_SML','KSPSAP_FLUX_LAG','SAP_BKG','SAP_X','SAP_Y'):
                try:
                    tt,yy,ee,units=column_arrays(cols,name); a,_=fit(tt,yy,p,ee); b,_=fit(tt,yy,p,ee,robust=True); result['columns'][name]={'units':units,'ordinary':a,'robust':b}
                except Exception as exc:result['columns'][name]={'error':str(exc)}
            ax=axs[sector-1]; ax[0].plot(t,y,'.',ms=2,label='SAP'); kt,ky,_,_=column_arrays(cols,'KSPSAP_FLUX'); ax[0].plot(kt,ky,'.',ms=2,label='KSP'); ax[0].legend(); ax[0].set_title(f'TIC {tic} S{sector} C{cam}/{ccd}')
            if np.sum(np.isfinite(mode))>=100:
                try:
                    good,clean,reg=residualize(t,y,p,mode); a,_=fit(t[good],clean,p); b,_=fit(t[good],clean,p,robust=True); rg,rc,rr=residualize(t,y,p,mode,robust=True); ra,_=fit(t[rg],rc,p,robust=True)
                    result['common_mode_fit']={'joint_robust_regression':rr,'joint_robust_residual_fit':ra,'regression':reg,'ordinary':a,'robust':b,'samples':int(good.sum())}
                    ax[1].plot(t,mode,'.',ms=2,label='detector median'); ax[1].plot(t[good],clean,'.',ms=2,label='SAP minus fitted common mode'); ax[1].legend()
                    ax[2].plot(((t[good]-1340)/p)%1,clean,'.',ms=2); ax[2].set_title('Residual at fixed S1 period')
                except Exception as exc:result['common_mode_fit']={'error':str(exc)}
            else:result['common_mode_fit']={'status':'insufficient exact-cadence controls'}
            results.append(result); write_json(OUT/'partial-results.json',{'leads':results})
        fig.tight_layout();fig.savefig(OUT/f'tic-{tic}.png',dpi=140);plt.close(fig)
        for (other,sector),(t,y,cols,cam,ccd,row) in data.items():
            a,_=fit(t,y,p); b,_=fit(t,y,p,robust=True);comparisons.append({'lead_tic':tic,'tic':other,'sector':sector,'camera':cam,'ccd':ccd,'period_days':p,'ordinary':a,'robust':b})
    write_json(OUT/'results.json',{'declaration':declaration,'completed_utc':utc_now(),'s1_periods':search,'leads':results,'all24_comparisons':comparisons})
if __name__=='__main__':main()
