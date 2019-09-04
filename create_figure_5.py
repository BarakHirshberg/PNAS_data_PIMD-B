#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 19:47:09 2019

@author: hirshb
"""

import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

#f, ((ax1, ax2)) = plt.subplots(2, 1, sharex='col',figsize=(2*3.375,2*3.375),dpi=600)
f, ax1 = plt.subplots(1, 1, sharex='col',figsize=(2*3.375,3.375),dpi=600)

fname = "N2PEoo_PEO_g0_data.dat"
path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

(bhw_val, mPEoo, BSE_PEoo, mPEO, BSE_PEO) = pd.read_pickle(path + fname)

fname = "N2PEoo_PEO_g8_data.dat"
path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

(bhw_val_g8, mPEoo_g8, BSE_PEoo_g8, mPEO_g8, BSE_PEO_g8) = pd.read_pickle(path + fname)

#Plot for N=2
Natoms = 2
#ax1.plot(bhw_val,math.factorial(Natoms)*np.array(mPEoo),'oC1', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.plot(bhw_val,Natoms*np.array(mPEO),'-oC0', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
ax1.errorbar(bhw_val, Natoms*np.array(mPEO), Natoms*np.array(BSE_PEO),  fmt='-oC0', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.errorbar(bhw_val, math.factorial(Natoms)*np.array(mPEoo), math.factorial(Natoms)*np.array(BSE_PEoo),  fmt='-oC0', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)

#ax1.plot(bhw_val,math.factorial(Natoms)*np.array(mPEoo_g8),'sC1', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.plot(bhw_val,Natoms*np.array(mPEO_g8),'-sC1', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
ax1.errorbar(bhw_val, Natoms*np.array(mPEO_g8), Natoms*np.array(BSE_PEO_g8),  fmt='-sC1', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.errorbar(bhw_val, math.factorial(Natoms)*np.array(mPEoo_g8), math.factorial(Natoms)*np.array(BSE_PEoo_g8),  fmt='-sC1', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)


#ax1.axhline(y=1/math.factorial(Natoms),linewidth=2, color='k')
#ax1.axhline(y=1/Natoms,linewidth=2, color='k')

ax1.set_ylabel(r'$P_{c} / P_{\infty}$',fontsize=12)
#plt.ylabel(r'$\langle E \rangle $')
#plt.xlabel(r'$ \beta \hbar \omega _0$',fontsize=24,fontweight='bold')
##ax1.set_xlabel(r'$ \beta \hbar \omega_0 $',fontsize=12)
##plt.show()
##ax1.tick_params(axis="x",labelsize=12)
ax1.tick_params(axis="y",labelsize=12)

fname = "N32PEoo_PEO_g0_data.dat"
path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

(bhw_val, mPEoo, mPEO, BSE_PEO) = pd.read_pickle(path + fname)

fname = "N32PEoo_PEO_g3_data.dat"
path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

(bhw_val_g3, mPEoo_g3, mPEO_g3, BSE_PEO_g3) = pd.read_pickle(path + fname)

#Plot for N=32
Natoms = 32
#ax2.plot(bhw_val,mPEoo,'oC0', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.plot(bhw_val,Natoms*np.array(mPEO),'-oC2', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
ax1.errorbar(bhw_val, Natoms*np.array(mPEO), Natoms*np.array(BSE_PEO),  fmt='-oC2', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax2.plot(bhw_val,mPEoo,'oC0', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
#ax1.plot(bhw_val_g3,Natoms*np.array(mPEO_g3),'-sC3', markersize=6, mfc='w', linewidth=2.0, mew=2.0)
ax1.errorbar(bhw_val_g3, Natoms*np.array(mPEO_g3), Natoms*np.array(BSE_PEO_g3),  fmt='-sC3', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)

#ax2.plot(bhw_val_g3,mPEoo_g3,'sC0', markersize=6, mfc='w', linewidth=2.0, mew=2.0)

#ax2.axhline(y=1/math.factorial(Natoms),linewidth=2, color='k')
#ax2.axhline(y=1/Natoms,linewidth=2, color='k')

ax1.set_ylabel(r'$P_{c}/P_{\infty}$',fontsize=12)
#plt.ylabel(r'$\langle E \rangle $')
#plt.xlabel(r'$ \beta \hbar \omega _0$',fontsize=24,fontweight='bold')
ax1.set_xlabel(r'$ \beta \hbar \omega_0 $',fontsize=12)
##plt.show()
#ax2.set_yticks([0.0, 0.01,0.02, 0.03])
ax1.tick_params(axis="x",labelsize=12)
ax1.tick_params(axis="y",labelsize=12)

f.tight_layout()
##f.subplots_adjust(hspace = 0.001)
f.savefig('/mnt/storage3/hirshb/NEW_LAMMPS_DATA/'+ 'Figure5.pdf', format='pdf', dpi=600)
