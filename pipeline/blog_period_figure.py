# Figure for the blog post: why only one period survives (TIC 402710280.01).
import csv, numpy as np, matplotlib, datetime as dt
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.dates as mdates
from astropy.time import Time
REL='/home/user/ad-astra/release/tic402710280_zenodo/data/ffi_photometry/tic402710280_s%d_ffi.csv'
P=36.33986; T0=2460838.2008
def date(jd): return mdates.date2num(Time(jd,format='jd').to_datetime())
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTE='#cfcec8'; BAR='#e4e3de'; BLUE='#2a78d6'; ORANGE='#eb6834'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.edgecolor':MUTE,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
fig,(a1,a2)=plt.subplots(2,1,figsize=(10.5,6.6),gridspec_kw={'height_ratios':[0.75,1.25],'hspace':0.45},facecolor=SURF)
for a in (a1,a2): a.set_facecolor(SURF)
# 1: timeline
sectors={'S12':('2019-05-24','2019-06-18'),'S39':('2021-05-27','2021-06-24'),'S66':('2023-06-02','2023-06-27'),'S93':('2025-06-04','2025-06-28'),'S100–103':('2026-02-03','2026-05-16')}
for name,(s,e) in sectors.items():
    s,e=mdates.date2num(dt.datetime.fromisoformat(s)),mdates.date2num(dt.datetime.fromisoformat(e))
    a1.axvspan(s,e,ymin=0.45,ymax=0.75,color=BAR,lw=0)
    a1.text(s+(e-s)/2,0.93,name,ha='center',va='bottom',color=INK2,fontsize=9,transform=a1.get_xaxis_transform())
obs={'S39':2459384.6094,'S93':2460838.1847,'S100':2461092.5845,'S102':2461128.9446}
for k,v in obs.items(): a1.plot(date(v),0.6,'o',ms=8,color=BLUE,transform=a1.get_xaxis_transform(),zorder=3)
a1.plot(date(T0-20*P),0.6,'o',ms=8,mfc=SURF,mec=BLUE,mew=1.6,transform=a1.get_xaxis_transform(),zorder=3)
a1.annotate('2023: lost to\nscattered light',(date(T0-20*P),0.5),xytext=(date(T0-20*P),0.36),ha='center',va='top',color=INK2,fontsize=8.5,xycoords=a1.get_xaxis_transform(),textcoords=a1.get_xaxis_transform())
a1.annotate('',xy=(date(obs['S93']),0.12),xytext=(date(obs['S39']),0.12),xycoords=a1.get_xaxis_transform(),arrowprops=dict(arrowstyle='<->',color=INK2,lw=0.9))
a1.text(date((obs['S39']+obs['S93'])/2),0.0,'1,453.58 days between the two public transits',ha='center',va='bottom',color=INK,fontsize=9,transform=a1.get_xaxis_transform())
a1.set_ylim(0,1); a1.set_yticks([]); a1.spines['left'].set_visible(False)
a1.set_xlim(mdates.date2num(dt.datetime(2019,1,1)),mdates.date2num(dt.datetime(2026,9,1)))
a1.xaxis.set_major_locator(mdates.YearLocator()); a1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
a1.set_title('1  Every time TESS has watched this star, and the transits it caught',loc='left',fontsize=10.5,color=INK)
# 2: four close-ups in 2026
a2.remove()
sub=gs2=fig.add_gridspec(2,1,height_ratios=[0.75,1.25],hspace=0.75)[1].subgridspec(1,4,wspace=0.08)
T=[];F=[]
for s_ in (100,101,102,103):
    r=list(csv.DictReader(open(REL%s_))); t=np.array([float(x['bjd_tdb']) for x in r]); f=np.array([float(x['flux']) for x in r])
    ok=np.isfinite(f); T.append(t[ok]); F.append(f[ok])
t=np.concatenate(T); f=np.concatenate(F)
cases=[(obs['S100'],'21 Feb: transit',True),(obs['S100']+P/2,'11 Mar: half-period\nprediction',False),(obs['S102'],'29 Mar: transit',True),(obs['S102']+P/2,'16 Apr: half-period\nprediction',False)]
W=1.0; D=5.27/24; axs=[]
for j,(tc,lab,real) in enumerate(cases):
    ax=fig.add_subplot(sub[0,j]); ax.set_facecolor(SURF)
    m=np.abs(t-tc)<W; tt=t[m]-tc; ff=f[m]
    oot=np.abs(tt)>0.3
    c=np.polyfit(tt[oot],ff[oot],2); y=(ff/np.polyval(c,tt)-1)*1e6
    ax.plot(tt*24,y,'.',ms=1.6,color=MUTE)
    e=np.arange(-W,W+1e-9,0.5/24); k=np.digitize(tt,e)
    bx=[tt[k==q].mean()*24 for q in np.unique(k) if (k==q).sum()>=5]; by=[y[k==q].mean() for q in np.unique(k) if (k==q).sum()>=5]
    ax.plot(bx,by,'o',ms=3.2,color=INK2)
    ax.axvspan(-D/2*24,D/2*24,color=BLUE if real else ORANGE,alpha=0.12,lw=0)
    inn=np.abs(tt)<D/2*0.8; d=-np.mean(y[inn]); sg=np.median(np.abs(np.diff(y)))*1.4826/np.sqrt(2)/np.sqrt(inn.sum())
    if not real: ax.text(0,-2250,'a transit would dip\n~2,000 ppm here\nmeasured: %s ± %d ppm'%(('%d'%round(d,-1)).replace('-','\u2212'),round(sg,-1)),ha='center',va='top',fontsize=8.5,color=INK2)
    ax.set_title(lab,fontsize=9.5,color=INK if real else ORANGE,loc='center')
    ax.set_ylim(-3600,1300); ax.set_xlim(-W*24,W*24); ax.set_xticks([-12,0,12])
    if j==0: ax.set_ylabel('brightness change (ppm)')
    else: ax.set_yticklabels([]); ax.tick_params(axis='y',length=0)
    if j==1: ax.set_xlabel('hours from the predicted time',x=1.04)
    print(lab.replace(chr(10),' '), round(d), round(sg)); axs.append(ax)
pos=axs[0].get_position(); fig.text(pos.x0-0.06,pos.y1+0.085,'2  My photometry of the 2026 images: two transits 36.36 days apart, and nothing at the half-period times',fontsize=10.5,color=INK)
fig.savefig('/home/user/ad-astra/docs/blog/period.png',dpi=200,bbox_inches='tight',facecolor=SURF)
