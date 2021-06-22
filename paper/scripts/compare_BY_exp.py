import pdb
from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
plt.rcParams.update({'font.size': 16})

# Data file path
pth = '/Users/saadjansari/Documents/Projects/ImageAnalysis/SingleCell/Results/Test_fourier3/by_length_comparison.csv'

my_data = genfromtxt(pth, delimiter=',')

times = my_data[:,0::3]
exp = my_data[:,1::3]*0.05
fit = my_data[:,2::3]*0.05

nCells = times.shape[1]
if not os.path.exists('/Users/saadjansari/Desktop/by_plots'):
    os.mkdir('/Users/saadjansari/Desktop/by_plots')

for jc in range(nCells):

    fig,ax = plt.subplots()
    ax.plot(times[:,jc],exp[:,jc],'r.', label = 'Hand-measured', markersize=10, alpha=0.7)
    ax.plot(times[:,jc],fit[:,jc],'g.', label = 'Automated', markersize=10, alpha=0.7)
    ax.set(xlabel='Time (s)', ylabel = r'Length ($\mu m$)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('/Users/saadjansari/Desktop/by_plots/cell_{}.pdf'.format(jc))

cmap = plt.cm.get_cmap('jet')
cols = cmap( np.linspace(0,255,nCells, dtype=int) )

fig,ax = plt.subplots()
for jc in range(nCells):

    ax.plot(times[:,jc],exp[:,jc],'--', alpha=0.7, color=cols[jc,:])
    ax.plot(times[:,jc],fit[:,jc],'.', alpha=0.7, color=cols[jc,:])
    ax.set(xlabel='Time (s)', ylabel = '2D Length (pixels)')

# plt.legend()
plt.tight_layout()
plt.savefig('/Users/saadjansari/Desktop/by_plots/all_cells.pdf')
