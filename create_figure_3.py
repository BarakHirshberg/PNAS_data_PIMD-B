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

#f, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2, 2,figsize=(2*3.375,1.33*2*3.375),dpi=600)
f, ((ax1,ax2,ax3,ax4)) = plt.subplots(4, 1,figsize=(2*3.375,2*2*3.375),dpi=600)
#f1 = plt.figure(figsize=(2*3.375,0.33*2*3.375),dpi=600)
#ax1 = f1.add_subplot(1, 1, 1)
#
#f2 = plt.figure(figsize=(2*3.375,2*3.375),dpi=600)
#ax2 = f2.add_subplot(3, 1, 1)
#ax3 = f2.add_subplot(3, 1, 2, sharex = ax2, sharey=ax2)
#ax4 = f2.add_subplot(3, 1, 3, sharex = ax2, sharey=ax2)

#f2.subplots_adjust(hspace=0.001)

#ax1 = f1.add_subplot(4, 1, 4)

#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax3.get_xticklabels(), visible=False)

path = "/mnt/storage3/hirshb/NEW_LAMMPS_DATA/"

#Plot for N=2
fname = "N2_E_vs_g_data.dat"
(g_values, E, err) = pd.read_pickle(path + fname)

ax1.errorbar(g_values, np.divide(E,1.0), 3*np.divide(err,1.0),  fmt='oC0',label=r'$\beta \hbar \omega_0 = 3$', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)

path_Mujal = "/mnt/storage3/hirshb/P_Mujal_Data/Fig1"

fig_E = 'Fig1'
panel='b'

#import and plot
data_file = path_Mujal + '/' +fig_E + '_panel_'+panel+'.txt'
data_col = [0,1]
tmp = pd.read_csv(data_file, sep='\s+', comment='#', usecols=[0,1], header=None)
tmp.columns=['0.5g','Egs']    
ax1.plot(2*tmp['0.5g'].iloc[:-1],tmp['Egs'].iloc[:-1],'-C0',markersize=10, mfc='w', linewidth=4.0, mew=2.0)

#Plot for N=3
fname = "N3_E_vs_g_data.dat"
(g_values, E, err) = pd.read_pickle(path + fname)

ax1.errorbar(g_values, np.divide(E,1.0), 3*np.divide(err,1.0),  fmt='oC1',label=r'$\beta \hbar \omega_0 = 3$', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)

path_Mujal = "/mnt/storage3/hirshb/P_Mujal_Data/Fig2_and_Fig3"

fig_E = 'Fig2'
panel='a'

#import and plot
data_file = path_Mujal + '/' +fig_E + '_panel_'+panel+'.txt'
data_col = [0,1]
tmp = pd.read_csv(data_file, sep='\s+', comment='#', usecols=[0,1], header=None)
tmp.columns=['0.5g','Egs']    
ax1.plot(2*tmp['0.5g'],tmp['Egs'],'-C1',markersize=10, mfc='w', linewidth=4.0, mew=2.0)

#Plot for N=4
fname = "N4_E_vs_g_data.dat"
(g_values, E, err) = pd.read_pickle(path + fname)

ax1.errorbar(g_values, np.divide(E,1.0), 3*np.divide(err,1.0),  fmt='oC2',label=r'$\beta \hbar \omega_0 = 3$', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0)
ax1.set_ylabel(r'$\langle E \rangle / \hbar \omega _0  $',fontsize=12)
ax1.set_xlabel(r'$ g$',fontsize=12)
ax1.set_xticks([0.0,4.0,8.0,12.0,16.0])
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax1.tick_params(axis="x", labelsize=12)
ax1.tick_params(axis="y", labelsize=12)

path_Mujal = "/mnt/storage3/hirshb/P_Mujal_Data/Fig2_and_Fig3"

fig_E = 'Fig2'
panel='b'

#import and plot
data_file = path_Mujal + '/' +fig_E + '_panel_'+panel+'.txt'
data_col = [0,1]
tmp = pd.read_csv(data_file, sep='\s+', comment='#', usecols=[0,1], header=None)
tmp.columns=['0.5g','Egs']    
ax1.plot(2*tmp['0.5g'],tmp['Egs'],'-C2',markersize=10, mfc='w', linewidth=4.0, mew=2.0)

##ax1.set_xticks([-2,0,2])
##ax1.set_xlim([-2.0,2.0])
##ax1.tick_params(axis="x", labelsize=12)

#f1.tight_layout()
#f1.savefig(path + 'Figure3a.pdf', format='pdf', dpi=600)

#Now for the density
g_values=[1,3,8,16]

#N=2
#colors={1: 'C0', 3: 'C1', 8:'C2', 16:'C3'}
colors={1: 'k', 3: 'k', 8:'k', 16:'k'}

N=2
path_Mujal = "/mnt/storage3/hirshb/P_Mujal_Data/"
panel='a'
mean_err_N2=[]
for g in g_values:
    fname = "N"+str(N)+"_g"+str(g)+"_dens_data.dat"
    (bins, dens, err) = pd.read_pickle(path + fname)
    ax2.errorbar(bins, dens, 3*err,  fmt='o', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0,zorder=0,alpha=1.0)

    data_file = 'Fig5' + '_g' + str(g) + '_N' + str(N) + '_panel_'+panel+'.txt'
    data_col = [0,1]    
    col_names = ['X','rho']
    tmp2 = analyze_col_from_COLVAR(path_Mujal + 'Fig5' + '/', data_file, data_col, col_names, start=0, end=101)

    ax2.plot(tmp2['X'],tmp2['rho'],'-', color=colors[g], markersize=6, mfc='w', linewidth=2.0, mew=2.0, zorder=1,alpha=1.0)
    ax2.set_ylabel(r'$\rho(r)$',fontsize=12)
    ax2.tick_params(axis="y", labelsize=12)
    ax2.tick_params(axis="x", labelsize=12)
    ax2.set_ylim([-0.01,0.31])
    ax2.set_yticks([0.0,0.1,0.2,0.3])
    
    mean_err_N2.append(np.mean(100*np.divide(np.abs(dens-tmp2['rho']),tmp2['rho'])))

#N=3
N=3
panel='b'
mean_err_N3=[]
for g in g_values:
    fname = "N"+str(N)+"_g"+str(g)+"_dens_data.dat"
    (bins, dens, err) = pd.read_pickle(path + fname)
    ax3.errorbar(bins, dens, 3*err,  fmt='o', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0,zorder=0,alpha=1.0)

    data_file = 'Fig5' + '_g' + str(g) + '_N' + str(N) + '_panel_'+panel+'.txt'
    data_col = [0,1]    
    col_names = ['X','rho']
    tmp2 = analyze_col_from_COLVAR(path_Mujal + 'Fig5' + '/', data_file, data_col, col_names, start=0, end=101)

    ax3.plot(tmp2['X'],tmp2['rho'],'-',color=colors[g], markersize=6, mfc='w', linewidth=2.0, mew=2.0,zorder=1,alpha=1.0)
    ax3.set_ylabel(r'$\rho(r)$',fontsize=12)
    ax3.tick_params(axis="y", labelsize=12)
    ax3.tick_params(axis="x", labelsize=12)
    ax3.set_ylim([-0.01,0.31])
    ax3.set_yticks([0.0,0.1,0.2,0.3])
    
    mean_err_N3.append(np.mean(100*np.divide(np.abs(dens-tmp2['rho']),tmp2['rho'])))
       
#N=4
N=4
panel='c'
mean_err_N4=[]
for g in g_values:
    fname = "N"+str(N)+"_g"+str(g)+"_dens_data.dat"
    (bins, dens, err) = pd.read_pickle(path + fname)
    ax4.errorbar(bins, dens, 3*err,  fmt='o', capsize=2, markersize=6, mfc='w', linewidth=2.0, mew=2.0,zorder=0,alpha=1.0)

    data_file = 'Fig5' + '_g' + str(g) + '_N' + str(N) + '_panel_'+panel+'.txt'
    data_col = [0,1]    
    col_names = ['X','rho']
    tmp2 = analyze_col_from_COLVAR(path_Mujal + 'Fig5' + '/', data_file, data_col, col_names, start=0, end=101)

    ax4.plot(tmp2['X'],tmp2['rho'],'-',color=colors[g], markersize=6, mfc='w', linewidth=2.0, mew=2.0, zorder=1,alpha=1.0)
    ax4.set_ylabel(r'$\rho(r)$',fontsize=12)
    ax4.set_xlabel(r'$r\sqrt{m\omega _0 / \hbar}$',fontsize=12)
    ax4.tick_params(axis="x", labelsize=12)
    ax4.tick_params(axis="y", labelsize=12)
    ax4.set_ylim([-0.01,0.31])
    ax4.set_yticks([0.0,0.1,0.2,0.3])
    
    mean_err_N4.append(np.mean(100*np.divide(np.abs(dens-tmp2['rho']),tmp2['rho'])))

#f2.tight_layout()
#f2.savefig(path + 'Figure3b.pdf', format='pdf', dpi=600)

f.tight_layout()
f.savefig(path + 'Figure3_Rossi.pdf', format='pdf', dpi=600)