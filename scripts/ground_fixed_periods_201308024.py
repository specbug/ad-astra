#!/usr/bin/env python3
"""Descriptive fixed-TESS-period tests of public legacy ASAS-SN V data.
No period optimization, significance calibration, or TESS phase extrapolation.
"""
import csv, hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA=Path('data/ground-validation-201308024')
OUT=Path('runs/ground-validation-201308024-v2')
PERIODS=[2.03270208338902,4.194743460281676]
PHASES=[0,.25,.5,.75]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(label):
 rows=list(csv.DictReader(open(DATA/(label+'.csv'))))
 t=np.array([float(x['hjd']) for x in rows]);y=np.array([float(x['mag']) for x in rows]);e=np.array([float(x['mag_err']) for x in rows]);c=np.array([x['camera'] for x in rows]);ok=np.isfinite(t)&np.isfinite(y)&np.isfinite(e)&(e>0)
 t,y,e,c=t[ok],y[ok],e[ok],c[ok];order=np.argsort(t,kind='stable');t,y,e,c=t[order],y[order],e[order],c[order];season=np.cumsum(np.r_[0,np.diff(t)>90]);g=np.array([f'{s}/{cam}' for s,cam in zip(season,c)]);base=np.column_stack([g==a for a in sorted(set(g))]).astype(float)
 return t,y,e,c,base,{'rows_raw':len(rows),'rows_used':len(t),'rows_rejected':len(rows)-len(t),'baseline_days':float(np.ptp(t)),'unique_times':len(set(t)),'groups':{str(a):int(sum(g==a)) for a in sorted(set(g))},'median_quoted_error_mag':float(np.median(e)),'raw_rms_mag':float(np.std(y))}
def fit(t,y,e,base,p):
 z=2*np.pi*(t-t[0])/p;X=np.column_stack([base,np.sin(z),np.cos(z)]);a=X/e[:,None];b=y/e;coef,_,rank,_=np.linalg.lstsq(a,b,rcond=None)
 if rank!=X.shape[1]:raise ValueError('Rank-deficient fixed-period model')
 residual=y-X@coef;cov=np.linalg.inv(a.T@a);s,co=coef[-2:];amp=float(np.hypot(s,co));v=np.array([s,co])/amp if amp else np.array([1.,0.]);unc=float(np.sqrt(v@cov[-2:,-2:]@v));red=float(sum((residual/e)**2)/(len(t)-rank))
 return {'period_days':p,'sin_mag':float(s),'cos_mag':float(co),'semiamplitude_mag':amp,'approx_fractional_flux_ppm':amp*np.log(10)/2.5*1e6,'formal_white_noise_amp_error_mag':unc,'residual_scaled_amp_error_mag':unc*np.sqrt(red),'reduced_chi2':red,'residual_rms_mag':float(np.std(residual)),'degrees_of_freedom':int(len(t)-rank),'interpretation':'Conditional descriptive fit; quoted errors/floor, correlated ground noise, offsets and fixed-frequency choice prevent a calibrated false-alarm claim.'},coef,X@coef

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'declaration.json').exists():raise FileExistsError('Preserve existing experiment; no overwrite')
 labels=['target','bright-neighbor','field-control-AP37233363']
 save(OUT/'declaration.json',{'utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(Path(__file__)),'input_sha256':{n:sha(DATA/(n+'.csv')) for n in labels},'periods_days':PERIODS,'analysis':'Weighted magnitude sin/cos at exactly two TESS-selected periods, additive offsets per camera and season split at gaps >90d; no clipping, no optimization; finite HJD/mag and positive finite own mag error required.','time':'HJD field-center, timescale unspecified; use arbitrary first-ground-observation epoch, no TESS phase extrapolation.','injection':'Multiply inferred relative flux by 1+0.001754*sin(2*pi*(hjd-hjd0)/4.194743460281676+2*pi*phase); exact magnitude transform, original errors/times/offset design. Phases [0,.25,.5,.75]; paired coefficient increments assess linear recovery, not a detection-completeness estimate.','phases_cycles':PHASES,'flux_semiamplitude_injection':.001754,'no_calibrated_fap':True})
 results={'sources':{},'injections':[]}
 for label in labels:
  t,y,e,c,base,meta=load(label);nullcoef=np.linalg.lstsq(base/e[:,None],y/e,rcond=None)[0];res=y-base@nullcoef;meta['offset_only_rms_mag']=float(np.std(res));meta['offset_only_reduced_chi2']=float(sum((res/e)**2)/(len(t)-base.shape[1]));records=[]
  for p in PERIODS:
   rec,coef,model=fit(t,y,e,base,p);records.append(rec)
  results['sources'][label]={'metadata':meta,'fits':records}
  fig,axes=plt.subplots(1,3,figsize=(13,3.5));axes[0].errorbar(t-t[0],res,yerr=e,fmt='.',alpha=.5);axes[0].set(xlabel='Days from first HJD',ylabel='V mag minus camera/season offset',title=label)
  for ax,p,rec in zip(axes[1:],PERIODS,records):
   ph=((t-t[0])/p)%1;ax.errorbar(ph,res,yerr=e,fmt='.',alpha=.45);xx=np.linspace(0,1,400);ax.plot(xx,rec['sin_mag']*np.sin(2*np.pi*xx)+rec['cos_mag']*np.cos(2*np.pi*xx));ax.set(xlabel='Ground-only phase',title=f'Fixed P={p:.6f} d')
  fig.tight_layout();fig.savefig(OUT/(label+'.png'),dpi=150);plt.close(fig)
  if label=='target':
   p=PERIODS[1];baseline,_,_=fit(t,y,e,base,p)
   for ph in PHASES:
    delta=-2.5*np.log10(1+.001754*np.sin(2*np.pi*(t-t[0])/p+2*np.pi*ph));rec,_,_=fit(t,y+delta,e,base,p);ds=rec['sin_mag']-baseline['sin_mag'];dc=rec['cos_mag']-baseline['cos_mag'];results['injections'].append({'phase_cycles':ph,'injected_flux_semiamplitude':.001754,'injected_mag_approx_semiamplitude':.001754*2.5/np.log(10),'recovered_paired_coefficient_semiamplitude_mag':float(np.hypot(ds,dc)),'fit_after_injection':rec,'note':'Paired recovery tests fitting response only; whether a signal stands out in real noise remains uncalibrated.'})
 results['window']={'target_fixed_frequency_spectral_window':{str(f):float(abs(np.mean(np.exp(2j*np.pi*f*load('target')[0])))) for f in [1,2,1/365.25,1/PERIODS[0],1/PERIODS[1],1-1/PERIODS[1],1+1/PERIODS[1]]},'note':'Unweighted window modulus, qualitative daily/yearly alias warning; not calibrated noise.'}
 save(OUT/'results.json',results)
if __name__=='__main__':main()
