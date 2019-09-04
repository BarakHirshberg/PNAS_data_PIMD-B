"""
Created on Fri May 10 19:47:09 2019

@author: hirshb
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from auxfunctions_new import analyze_col_from_COLVAR
from matplotlib.ticker import FormatStrFormatter

#f, ((ax1,ax2)) = plt.subplots(2, 1,figsize=(2*3.375,0.66*2*3.375),dpi=600,  sharex=True, )
f, ((ax1,ax2)) = plt.subplots(2, 1,figsize=(2*3.375,0.75*2*3.375),dpi=600,  sharex=True, )

path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

#Now for the density
g_values=[0,1,3]
colors=['k','C0','C1']
count=0
#N=32
N=32
for g in g_values:
    fname = "N"+str(N)+"_g"+str(g)+"_dens_data.dat"
    (bins, dens, err) = pd.read_pickle(path + fname)
    ax1.errorbar(bins, dens, 3*err,  fmt='o', color=colors[count], capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
    count += 1
ax1.set_ylabel(r'$\rho(r)$',fontsize=12)
ax1.tick_params(axis="y", labelsize=12)
    

#N=32
N=32
count=0
for g in g_values:
    fname = "N"+str(N)+"_g"+str(g)+"_pair_corr_data.dat"
    (bins, dens, err) = pd.read_pickle(path + fname)
    ax2.errorbar(bins, dens, 3*err,  fmt='o', color=colors[count], capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
    count += 1
ax2.set_xlabel(r'$r\sqrt{m\omega _0 / \hbar}$',fontsize=12)
ax2.set_ylabel(r'$c(r)$',fontsize=12)
ax2.set_yticks([0.0,0.04,0.08,0.12,0.16])
#ax2.set_yticks([0.0,0.1,0.2])
ax2.tick_params(axis="y", labelsize=12)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    
plt.tight_layout()
#f.subplots_adjust(hspace=0.001)
plt.savefig(path + 'Figure4_NEW.pdf', format='pdf', dpi=600)
