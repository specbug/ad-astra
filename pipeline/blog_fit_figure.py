# Figure for the blog post: the linear ephemeris of TIC 402710280.01 (paper values: Tables 2 and 3).
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
T=np.array([2459384.6094,2460838.1847,2461092.5845,2461128.9446]); s=np.array([0.0013,0.0012,0.0023,0.0017])
E=np.array([0,40,47,48]); OC=np.array([4.4,-23.2,6.8,35.9])
P,sP=36.33986,0.00051; T0,sT0=2460838.2008,0.0119     # reference epoch 40
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; MUTE='#cfcec8'; BLUE='#2a78d6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.edgecolor':MUTE,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.spines.top':False,'axes.spines.right':False})
fig,(a1,a2)=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={'wspace':0.3},facecolor=SURF)
for a in (a1,a2): a.set_facecolor(SURF)
ee=np.linspace(-2,50,100)
a1.plot(ee,T0+(ee-40)*P-T[0],color=BLUE,lw=2,zorder=1)
a1.plot(E,T-T[0],'o',ms=8,color=INK,mec=SURF,mew=2,zorder=3)
a1.annotate('Sector 39 (2021)',(0,0),xytext=(3,-60),fontsize=8.5,color=INK2)
a1.annotate('Sector 93 (2025)',(40,T[1]-T[0]),xytext=(37.5,T[1]-T[0]),ha='right',va='center',fontsize=8.5,color=INK2)
a1.annotate('Sectors 100 & 102 (2026)',(47.5,T[2]-T[0]+20),xytext=(45,T[3]-T[0]+120),ha='right',va='center',fontsize=8.5,color=INK2)
a1.text(4,1350,'slope = period\n$P$ = 36.33986 ± 0.00051 days',color=INK,fontsize=10)
a1.set_xlabel('epoch $E$ (orbits since June 2021)'); a1.set_ylabel('mid-time (days since the first transit)')
a1.set_xlim(-3,51); a1.set_ylim(-150,2050)
a1.set_title('1  Mid-time against epoch',loc='left',fontsize=10.5,color=INK)
band=np.sqrt(sT0**2+((ee-40)*sP)**2)*1440
a2.fill_between(ee,-band,band,color=BLUE,alpha=0.12,lw=0)
a2.axhline(0,color=BLUE,lw=2)
a2.errorbar(E,OC,yerr=s*1440,fmt='o',ms=8,color=INK,mec=SURF,mew=2,elinewidth=1.4,capsize=0,zorder=3)
for e,r in zip(E,OC):
    t=('%+.0f min'%r).replace('-','−')
    if e==0: xy=(0,10); ha='center'
    elif e==40: xy=(0,-16); ha='center'
    elif e==47: xy=(-10,6); ha='right'
    else: xy=(-10,-3); ha='right'
    a2.annotate(t,(e,r),xytext=xy,textcoords='offset points',ha=ha,fontsize=9,color=INK)
a2.set_xlim(-3,51); a2.set_ylim(-45,55)
a2.set_xlabel('epoch $E$'); a2.set_ylabel('observed − calculated (minutes)')
a2.set_title('2  Zoomed in: how far each transit sits off the line',loc='left',fontsize=10.5,color=INK)
a2.text(1,52,"error bars: each transit's timing uncertainty (2–3 min,\nmostly hidden behind the dots)\nshaded: uncertainty of the fitted line",fontsize=8.5,color=INK2,va='top')
fig.savefig('/home/user/ad-astra/docs/blog/period-fit.png',dpi=200,bbox_inches='tight',facecolor=SURF)
