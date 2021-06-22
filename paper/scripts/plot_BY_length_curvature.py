import pdb
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
plt.rcParams.update({'font.size': 16})

# Data file path
pth = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Test_fourier3/by_curvature.csv'

my_data = genfromtxt(pth, delimiter=',')

times = my_data[:,0::3]
lens = my_data[:,1::3]*0.05
curvature = np.abs(my_data[:,2::3]/0.05)
curvature = 1/curvature

nCells = times.shape[1]
if not os.path.exists('/Users/saadjansari/Desktop/by_plots'):
    os.mkdir('/Users/saadjansari/Desktop/by_plots')

for jc in range(nCells):

    fig,ax = plt.subplots(figsize=(8,3))
    ax.plot(times[:,jc],lens[:,jc],'r.', label = 'Length')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r'Length ($\mu m$)', color="red")
    # plt.legend()
    
    ax2 = ax.twinx()
    ax2.plot(times[:,jc],curvature[:,jc],'g.', label = 'Curvature')
    ax2.set(ylim=[0,10])
    ax2.set_ylabel(r'Radius of Curvature ($\mu m$)', color="green")

    # plt.legend()
    plt.tight_layout()
    plt.savefig('/Users/saadjansari/Desktop/by_plots/curvature_len_cell_{}.pdf'.format(jc))

all_curv = []
all_len = []
for jc in range(nCells):
    for lc,cc in zip( lens[:,jc], curvature[:,jc]):
        if not np.isnan(lc):
            all_len.append(lc)
            all_curv.append(cc)

fig,ax = plt.subplots(figsize=(8,3))
plt.scatter(lens[:,0], curvature[:,0], s=2, color='green')
ax.set_xlabel(r'Length ($\mu m$)')
ax.set_ylabel(r'$R_c$')
ax.set(xlim=[0,4], ylim=[0,20])

# plt.legend()
plt.tight_layout()
plt.savefig('/Users/saadjansari/Desktop/by_plots/len_vs_curvature.pdf')
    
# for jc in range(nCells):

    # dlen = np.abs(lens[1:,jc]-lens[:-1,jc])
    # dK = np.abs(curvature[1:,jc]-curvature[:-1,jc])
    # fig,ax = plt.subplots(figsize=(8,3))
    # ax.plot(times[1:,jc],dlen,'r.', label = 'Length')
    # ax.set_xlabel('Time (s)')
    # ax.set_ylabel(r'$\Delta$Length ($\mu m$)', color="red")
    # # plt.legend()
    
    # ax2 = ax.twinx()
    # ax2.plot(times[1:,jc],dK,'g.', label = 'Curvature')
    # ax2.set(ylim=[0,10])
    # ax2.set_ylabel(r'Changes in Radius of Curvature ($\mu m$)', color="green")

    # # plt.legend()
    # plt.tight_layout()
    # plt.savefig('/Users/saadjansari/Desktop/by_plots/changes_curvature_len_cell_{}.pdf'.format(jc))
