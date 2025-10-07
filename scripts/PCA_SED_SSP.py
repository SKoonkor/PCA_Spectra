import numpy as np
import matplotlib.pyplot as plt
import os, sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')
from visualisations import plot_SED

# Define essential directories
data_dir = '../data/'
output_dir = '../outputs/'


# Reading the normalised SEDs and PCs
wave_norm = np.genfromtxt(data_dir + 'FSPS_wave_masked.csv')
PCs = np.genfromtxt(output_dir + 'PCA_PCs.csv')

n_PCA = 5


# Plot the PCs
fig, ax = plt.subplots(n_PCA, 1, figsize = (7, n_PCA))

for i in range(n_PCA):
    plot_SED(ax[i], wave_norm, PCs[i], yscale = 'linear')

    ax[i].set_ylim(-0.1, 0.1)
    ax[i].set_ylabel('PC {}'.format(i+1))


    if i == n_PCA - 1:
        ax[i].set_xlabel('wavelength [A]')
    else:
        ax[i].set_xticks([])
     
plt.show()


