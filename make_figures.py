"""Regenerates figure1_decision_curve_mice and figure2_subgroup_calibration_mice from the CSVs in this folder."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
dca=pd.read_csv('decision_curve_mice.csv'); x=dca.threshold*100
fig,ax=plt.subplots(figsize=(6.5,4.2),dpi=300)
ax.fill_between(x,dca.nb_model_lo*1000,dca.nb_model_hi*1000,color='#1b5e8a',alpha=0.15,linewidth=0,label='Model, 95% CI')
ax.plot(x,dca.nb_model*1000,color='#1b5e8a',lw=2,label='Model-directed testing')
ax.plot(x,dca.nb_all*1000,color='#c0392b',lw=1.6,label='Test all')
ax.plot(x,dca.nb_c2*1000,color='#e29b1a',lw=1.6,label='Rule C2 (C1, age ≥60 y, or CVD)')
ax.axhline(0,color='grey',lw=1.2,ls='--',label='Test none / Rule C1'); ax.axvline(5,color='grey',lw=0.8,ls=':')
ax.text(5.3,58,'primary threshold (5%)',color='grey',fontsize=8,va='top')
ax.set_xlabel('Threshold probability of undiagnosed CKD (%)'); ax.set_ylabel('Net benefit (cases per 1,000 adults)')
ax.set_xlim(2,20); ax.set_ylim(-60,70); ax.set_xticks([2,5,10,15,20]); ax.legend(frameon=False,fontsize=8,loc='lower left')
fig.tight_layout(); fig.savefig('figure1_decision_curve_mice.png'); fig.savefig('figure1_decision_curve_mice.pdf')
sub=pd.read_csv('subgroup_oe_mice_ci.csv')
panels=[('race','A. Race and ethnicity'),('income','B. Family income-to-poverty ratio'),('insurance','C. Health insurance coverage')]
fig,axes=plt.subplots(3,1,figsize=(6.5,6.6),dpi=300,gridspec_kw={'height_ratios':[6,4,2.2]},sharex=True)
for ax,(dim,title) in zip(axes,panels):
    s=sub[sub.dimension==dim].copy(); s=s.sort_values('oe') if dim=='race' else s.iloc[::-1]; yy=np.arange(len(s))
    ax.axvline(1,color='grey',ls='--',lw=1)
    ax.errorbar(s.oe,yy,xerr=[s.oe-s.oe_lo,s.oe_hi-s.oe],fmt='o',color='#1b5e8a',ecolor='#1b5e8a',elinewidth=1.2,capsize=3,ms=5)
    ax.set_yticks(yy); ax.set_yticklabels([f"{l}  (n={n:,}; {e} events)" for l,n,e in zip(s.level.str.replace('>=','≥').str.replace('<1.0 (below FPL)','<1.0'),s.n,s.events)],fontsize=8)
    ax.set_ylim(-0.6,len(s)-0.4); ax.set_xlim(0,2.1); ax.set_title(title,fontsize=9,loc='left',fontweight='bold'); ax.tick_params(axis='x',labelbottom=True)
axes[-1].set_xlabel('Observed / expected ratio (>1 = model under-predicts risk)')
fig.tight_layout(); fig.savefig('figure2_subgroup_calibration_mice.png'); fig.savefig('figure2_subgroup_calibration_mice.pdf')
