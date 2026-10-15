# Figure for the blog post: the single-transit box search on the public Sector 39 TESS-SPOC light curve.
# Expects the light curve at dl/s39.fits (hlsp_tess-spoc_tess_phot_0000000402710280-s0039_tess_v1_lc.fits from MAST).
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
d=fits.getdata('dl/s39.fits',1)
t=np.array(d['TIME'],float); f=np.array(d['PDCSAP_FLUX'],float); q=np.array(d['QUALITY'])
ok=np.isfinite(t)&np.isfinite(f)&(q==0); t,f=t[ok],f[ok]
x=t-t.min()
trend=np.array([np.median(f[np.abs(x-xi)<1.0]) for xi in x]); ppm=(f/trend-1)*1e6
sig=np.median(np.abs(np.diff(ppm)))*1.4826/np.sqrt(2)
D=5.27/24
def score_at(c):
    inn=np.abs(x-c)<D/2; out=(np.abs(x-c)>=D/2)&(np.abs(x-c)<1.5*D)
    if inn.sum()<8 or out.sum()<8: return np.nan,np.nan
    dep=ppm[out].mean()-ppm[inn].mean(); return dep, dep/(sig*np.sqrt(1/inn.sum()+1/out.sum()))
grid=np.arange(x.min()+0.2, x.max()-0.2, 0.005)
res=np.array([score_at(c) for c in grid]); depth,score=res[:,0],res[:,1]
k=np.nanargmax(score); tc=grid[k]
far=np.abs(grid-tc)>1
print('cadence %.0f s, sig %.0f ppm, peak %.1f sigma depth %.0f at %.3f, max elsewhere %.1f, min elsewhere %.1f'%(np.median(np.diff(t))*86400,sig,score[k],depth[k],tc,np.nanmax(score[far]),np.nanmin(score[far])))
def binned(xx,yy,w):
    e=np.arange(xx.min(),xx.max()+w,w); i=np.digitize(xx,e); u=[j for j in np.unique(i) if (i==j).sum()>=3]
    return np.array([xx[i==j].mean() for j in u]), np.array([yy[i==j].mean() for j in u])
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTE='#cfcec8'; BLUE='#2a78d6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.edgecolor':MUTE,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(11,6.2),facecolor=SURF)
gs=fig.add_gridspec(2,2,width_ratios=[2.1,1],height_ratios=[1,1],hspace=0.18,wspace=0.16)
aA=fig.add_subplot(gs[0,0]); aC=fig.add_subplot(gs[1,0],sharex=aA); aB=fig.add_subplot(gs[:,1])
for a in (aA,aB,aC): a.set_facecolor(SURF)
# A: the whole visit
aA.plot(x,ppm,'.',ms=1.8,color=MUTE)
bx,by=binned(x,ppm,1/24); aA.plot(bx,by,'.',ms=3,color=INK2)
aA.set_ylim(-3500,1800); aA.set_ylabel('brightness change (ppm)')
aA.set_title('1  The light curve: 27 days of one star (TESS Sector 39, 2021)',loc='left',fontsize=10.5,color=INK)
aA.annotate('the transit',(tc,-2300),xytext=(tc-6,-3100),color=INK,fontsize=9,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
plt.setp(aA.get_xticklabels(),visible=False)
# C: the score
aC.plot(grid,score,color=BLUE,lw=1.4)
aC.axhline(6,color=INK2,lw=0.9,ls=':'); aC.text(x.max(),7.5,'threshold: 6σ',ha='right',color=INK2,fontsize=9)
aC.annotate('%.0fσ'%score[k],(tc,score[k]),xytext=(tc+1.2,score[k]-4),color=INK,fontsize=9,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
aC.set_ylim(-8,score[k]+6); aC.set_ylabel('dip strength (σ)'); aC.set_xlabel('days since the start of the sector')
aC.set_title('3  The score at every box position',loc='left',fontsize=10.5,color=INK)
# B: zoom, sliding box
w=1.1; m=np.abs(x-tc)<w
aB.plot(x[m]-tc,ppm[m],'.',ms=3,color=MUTE)
bx,by=binned(x[m]-tc,ppm[m],30/1440); aB.plot(bx,by,'o',ms=3.6,color=INK2)
def box(c,ls,lw,bold):
    dep,sc=score_at(c); o=c-tc
    aB.axvspan(o-D/2,o+D/2,color=BLUE,alpha=0.10 if not bold else 0.16,lw=0)
    xs=np.array([o-1.5*D,o-D/2,o-D/2,o+D/2,o+D/2,o+1.5*D]); ys=[0,0,-dep,-dep,0,0]
    aB.plot(xs,ys,color=BLUE,lw=lw,ls=ls,solid_capstyle='round',zorder=4)
    lab='%dσ'%round(abs(sc)) if round(sc)==0 else '%.0fσ'%sc
    aB.text(o,1350,lab,ha='center',color=INK if bold else INK2,fontsize=10 if bold else 9,weight='bold' if bold else 'normal')
box(tc-0.8,'--',1.5,False); box(tc-0.4,'--',1.5,False); box(tc,'-',2.2,True)
aB.annotate('',xy=(-0.15,1650),xytext=(-0.85,1650),arrowprops=dict(arrowstyle='->',color=INK2,lw=1))
aB.text(-0.5,1720,'slide the box',ha='center',va='bottom',color=INK2,fontsize=9)
aB.set_ylim(-3500,2100); aB.set_xlim(-w,w)
aB.set_xlabel('hours from the transit'); aB.set_xticks([-1,-0.5,0,0.5,1]); aB.set_xticklabels(['-24','-12','0','12','24'])
aB.set_title('2  A 5-hour box slid across the data',loc='left',fontsize=10.5,color=INK)
fig.savefig('/home/user/ad-astra/docs/blog/box-search.png',dpi=200,bbox_inches='tight',facecolor=SURF)
