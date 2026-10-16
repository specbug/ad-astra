# Figure for the blog post: stacked transit of TIC 402710280.01 against the best planet model and the best
# eclipse of the 3.5" background neighbour (k = 1, diluted to the observed depth), from blend_shape_402710280.json.
import csv, json, math, numpy as np, matplotlib, batman
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R='/home/user/ad-astra/'
res=json.load(open(R+'results/candidates/blend_shape_402710280.json'))
pl=res['planet']; eb=[e for e in res['eb'] if e['k']==1.0 and e['u']==[0.4,0.25]][0]
T=[];F=[]
for s in (39,93,100,102):
    for r in csv.DictReader(open(R+f'release/tic402710280_zenodo/data/transit_windows/tic402710280_s{s:03d}_window.csv')):
        if r['status']!='used': continue
        T.append((float(r['bjd_tdb'])-float(r['centre_bjd_tdb']))*24); F.append((float(r['flux_normalised'])-1)*1e6)
T=np.array(T); F=np.array(F)
P=36.3399
def model(k,b,aR,u,Fd,th):
    p=batman.TransitParams(); p.t0=0; p.per=P; p.rp=k; p.a=aR; p.inc=math.degrees(math.acos(min(b/aR,1))); p.ecc=0; p.w=90; p.u=list(u); p.limb_dark='quadratic'
    m=batman.TransitModel(p,th/24,supersample_factor=7,exp_time=200/86400).light_curve(p)
    return (-Fd*(1-m))*1e6
th=np.linspace(-6,6,1500)
mp=model(pl['k'],pl['b'],pl['a_rs'],(0.39,0.21),1.0,th)
me=model(1.0,eb['b'],eb['a_rs'],eb['u'],eb['F'],th)
e=np.arange(-6,6.001,0.25); idx=np.digitize(T,e)
bx=[];by=[];be=[]
for j in np.unique(idx):
    m=idx==j
    if m.sum()>5: bx.append(T[m].mean()); by.append(F[m].mean()); be.append(F[m].std()/np.sqrt(m.sum()))
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTE='#cfcec8'; BLUE='#2a78d6'; ORANGE='#eb6834'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.edgecolor':MUTE,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(9,4.6),facecolor=SURF); ax.set_facecolor(SURF)
ax.errorbar(bx,by,yerr=be,fmt='o',ms=4,color=INK2,elinewidth=1,capsize=0,zorder=2,label='all four transits, 15-minute averages')
ax.plot(th,mp,color=BLUE,lw=2.2,zorder=3,label='planet crossing TIC 402710280')
ax.plot(th,me,color=ORANGE,lw=2.2,ls='--',zorder=3,label='best eclipse of the faint neighbour')
ax.set_xlim(-5,5); ax.set_ylim(-2700,700)
ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('brightness change (ppm)')
ax.annotate('sharp edges: about\n14 minutes to ingress',(-2.55,-1100),xytext=(-4.85,-2350),fontsize=9,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
ax.annotate('a 74% eclipse of a star\nfades in and out slowly',(1.9,-1000),xytext=(2.6,-2350),fontsize=9,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.8))
ax.legend(frameon=False,loc='upper center',ncol=1,fontsize=9,bbox_to_anchor=(0.5,1.02))
ax.text(4.95,550,'Δχ² = %.0f in favour of the planet'%eb['dchi2'],ha='right',fontsize=9,color=INK)
fig.savefig(R+'docs/blog/blend-shape.png',dpi=200,bbox_inches='tight',facecolor=SURF)
print('planet',pl,'eb',eb['b'],eb['a_rs'],eb['F'],eb['dchi2'])
