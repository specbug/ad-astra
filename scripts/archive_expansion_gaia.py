#!/usr/bin/env python3
"""Validate published RR-neighbor morphology in independent Gaia DR3 epochs; no TESS access."""
from pathlib import Path
import json,warnings,xml.etree.ElementTree as ET
import numpy as np
from astropy.io.votable import parse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adastra.archive import write_json,sha256_file,utc_now
ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/'data/archive-expansion';OUT=ROOT/'runs/archive-expansion';OUT.mkdir(parents=True,exist_ok=True)
META=ROOT/'data/qlp-stellar-metadata/completion/summary.json'
reports=json.loads(META.read_text())['reports']
write_json(OUT/'declaration.json',{'created_utc':utc_now(),'script_sha256':sha256_file(Path(__file__)),'period_metadata_sha256':sha256_file(META),'hypothesis':'Known published RRab periods and morphology; descriptive 3-harmonic fixed-period fits. No novel source search and no cross-mission phase comparison.'})
result=[]
for report in reports:
 rr=report['published_rrlyrae_neighbors'][0];sid=rr['source_id'];path=IN/f'gaia-epoch-{sid}.xml';p=rr['pf']
 raw=ET.fromstring(path.read_bytes());ns={'v':'http://www.ivoa.net/xml/VOTable/v1.3'}
 ts=raw.find('.//v:TIMESYS',ns)
 if ts is None or ts.attrib!={'ID':'time_frame','refposition':'BARYCENTER','timeorigin':'2455197.5','timescale':'TCB'}:raise ValueError('Unexpected Gaia time definition')
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always');doc=parse(str(path),verify='warn')
 table=doc.get_first_table();params={x.name:x.value for x in table.params}
 if int(params['source_id'])!=sid:raise ValueError('Identity mismatch')
 a=table.array
 summary={'tic':report['tic'],'rr_gaia_source_id':sid,'published_period_days':p,'input_sha256':sha256_file(path),'time_definition':ts.attrib,'metadata_warnings':[str(x.message) for x in caught],'bands':{}}
 fig,axes=plt.subplots(1,3,figsize=(13,3.7))
 for ax,band in zip(axes,['g','bp','rp']):
  prefix='g_transit' if band=='g' else band
  def values(name):return np.ma.filled(a[name],np.nan).astype(float)
  t=values(prefix+'_time' if band=='g' else band+'_obs_time');flux=values(prefix+'_flux');err=values(prefix+'_flux_error')
  reject=np.ma.filled(a['variability_flag_'+band+'_reject'],True)|np.ma.filled(a['rejected_by_photometry'],True)
  mask=np.isfinite(t)&np.isfinite(flux)&np.isfinite(err)&(flux>0)&(err>0)&~reject
  t,flux,err=t[mask],flux[mask],err[mask]
  if len(t)<8:summary['bands'][band]={'retained':len(t),'status':'insufficient'};continue
  phase=((t-t.min())/p)%1.;norm=np.median(flux);y=flux/norm;e=err/norm
  design=np.column_stack([np.ones(len(t))]+[z for k in range(1,4) for z in [np.cos(2*np.pi*k*phase),np.sin(2*np.pi*k*phase)]])
  beta=np.linalg.lstsq(design/e[:,None],y/e,rcond=None)[0];pred=design@beta
  grid=np.linspace(0,1,1000);dgrid=np.column_stack([np.ones(len(grid))]+[z for k in range(1,4) for z in [np.cos(2*np.pi*k*grid),np.sin(2*np.pi*k*grid)]])
  model=dgrid@beta
  summary['bands'][band]={'raw_rows':len(a),'retained':len(t),'rejected_or_nonfinite':len(a)-len(t),'baseline_days':float(np.ptp(t)),'median_single_epoch_fractional_error':float(np.median(e)),'harmonic_flux_semiamplitudes_fractional':[float(np.hypot(beta[2*k-1],beta[2*k])) for k in range(1,4)],'fitted_peak_to_peak_magnitude':float(2.5*np.log10(model.max()/model.min())) if model.min()>0 else None,'residual_rms_fractional':float(np.std(y-pred)),'reduced_chi2_white_noise_only':float(np.sum(((y-pred)/e)**2)/(len(t)-7))}
  ax.errorbar(phase,y,yerr=e,fmt='.',ms=4,alpha=.7);ax.plot(grid,model,color='black');ax.set(xlabel='Phase at published period (own band epoch)',ylabel='Flux / median',title=f'{band.upper()}: {len(t)} retained epochs')
 fig.suptitle(f'Known RR neighbor Gaia {sid}; TIC {report["tic"]}');fig.tight_layout();fig.savefig(OUT/f'gaia-rr-{sid}.png',dpi=160);plt.close(fig);result.append(summary)
write_json(OUT/'gaia-rr-validation.json',{'completed_utc':utc_now(),'reports':result,'limitations':['Fits are descriptive, at periods already published from Gaia; not an independent discovery or new period measurement.','Gaia epochs are independently observed relative to TESS, but Gaia catalog periods and epoch photometry are from the same mission.','G, BP and RP sampling, masks and passbands differ; raw processing flags are preserved, not guessed.','Time is barycentric TCB JD minus 2455197.5; no subtraction of a TESS TDB epoch or cross-mission phase alignment.','This confirms known contaminants can vary strongly; it does not localize or quantitatively assign TESS flux contamination.']})
print(json.dumps(result,indent=2))
