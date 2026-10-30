from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent
load=lambda n:json.loads((P/n).read_text())
b,a=load('cache-native-before.json'),load('cache-fixed2.json')
ob,oa=load('angle-native-before.json'),load('angle-fixed.json')
plt.rcParams.update({'font.size':11,'axes.labelsize':11,'axes.titlesize':12})
fig,ax=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
for rows,label,color,style,marker in [(b,'Before','#D55E00','-','o'),(a,'Fixed','#0072B2','--','s')]:
 ax[0].plot([x['phi0_rad'] for x in rows],[x['field_error_T'] for x in rows],label=label,color=color,linestyle=style,marker=marker,markersize=3,markevery=3)
ax[0].axhline(0,color='black',linestyle=':',label='Exact constant field')
ax[0].set(xlabel='Initial cumulative toroidal angle [rad]',ylabel='Cartesian field error [T]',title='Field cached after a collision')
ax[0].ticklabel_format(axis='both',style='sci',scilimits=(0,0));ax[0].legend(loc='upper left')
for rows,label,color,style in [(ob,'Before','#D55E00','-'),(oa,'Fixed','#0072B2','--')]:
 hn=np.array([x['hnext_s'] for x in rows]);ax[1].plot([x['turns'] for x in rows],hn/hn[0]-1,color=color,linestyle=style,label=label)
ax[1].axhline(0,color='black',linestyle=':',label='Exact axial symmetry')
ax[1].set(xlabel='Added toroidal turns [1]',ylabel='Change in suggested step / unshifted step [1]',title='Embedded orbit-error estimate')
ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[1].legend(loc='upper right')
for x in ax:x.grid(alpha=.2)
fig.savefig(P/'consistency.png',dpi=180)
print('cache peak T',max(x['field_error_T'] for x in b),max(x['field_error_T'] for x in a))
print('symmetry max rel',max(abs(x['hnext_s']/ob[0]['hnext_s']-1) for x in ob),max(abs(x['hnext_s']/oa[0]['hnext_s']-1) for x in oa))
# Exact prescribed orbit updates must be retained at fixed steps.
for i,(x,y) in enumerate(zip(ob,oa)):
 for k in ['r_m','z_m','ppar_kgms','mu_JT']:assert x[k]==y[k],(i,k,x[k],y[k])
print('201 prescribed synthetic updates bitwise unchanged in R, Z, p_parallel and mu')
