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
curvature = np.abs(my_data[:,2::3]/0.05)
# curvature = 1/curvature

nCells = times.shape[1]
if not os.path.exists('/Users/saadjansari/Desktop/by_plots'):
    os.mkdir('/Users/saadjansari/Desktop/by_plots')

for jc in range(nCells):

    fig,ax = plt.subplots()
    ax.plot(times[:,jc],curvature[:,jc],'g.', markersize=10, alpha=0.8)
    ax.set(ylim=[-0.1,3.1], yticks=[0,1,2,3])
    ax.set_ylabel(r'Curvature ($\mu m^{-1}$)', wrap=True)
    ax.set_xlabel('Time (s)')

    # plt.legend()
    plt.tight_layout()
    plt.savefig('/Users/saadjansari/Desktop/by_plots/curvature_cell_{}.pdf'.format(jc))

