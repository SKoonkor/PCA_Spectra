import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ipywidgets import interact
import ipywidgets as widgets
import sys, os

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from visualisations import plot_SED
from PCA_cal import denormalise_SED

####
# Define filenames for inputs
data_dir = '../data/'
FSPS_param_grid_filename = 'FSPS_param_grid_LHS.csv'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_norm_filename = 'FSPS_SED_norm_std.csv'
SED_mean_filename = 'FSPS_SED_norm_mean_std.csv'
SED_factor_std_filename = 'FSPS_SED_norm_factor_std.csv'

SED_masked_filename = 'FSPS_SED_masked.csv'

# Read the FSPS template parameters
FSPS_param_grid = np.genfromtxt(data_dir + FSPS_param_grid_filename)

print (FSPS_param_grid[:5])
# Read the input
wave = np.genfromtxt(data_dir + wave_filename)
SED_norm = np.genfromtxt(data_dir + SED_norm_filename)
SED_factor = np.genfromtxt(data_dir + SED_factor_std_filename)
SED_mean = np.genfromtxt(data_dir + SED_mean_filename)

SED_masked = np.genfromtxt(data_dir + SED_masked_filename)





# Define the colors
colors = ['red', 'orange', 'yellow', 'green', 'blue']*4



##############################################################################################
# Denormalise SED

SED_denormalised = denormalise_SED(SED_norm, SED_factor, SED_mean)





##############################################################################################
##### Plotting the SEDs

fig, ax = plt.subplots(4, 1, figsize = (12, 8), gridspec_kw = {'hspace': 0.01} )

##### Plotting the original SEDs
# This is a set of masked SEDs
plot_SED(ax[0], wave, SED_mean, color = 'black', lw = 2, label = 'mean SED')
for i in range(20):
    plot_SED(ax[0], wave, SED_masked[i], color = colors[i], label = '[{:.5f}, {:.3f}]'.format(np.log10(FSPS_param_grid[i][0]),
                                                                                      FSPS_param_grid[i][1]))
ax[0].legend(loc = 'upper right', fontsize = 5)
ax[0].set_ylabel('Original SEDs')

##### Plotting the normalised SEDs
for i in range(20): 
    plot_SED(ax[1], wave, SED_norm[i], color = colors[i])
ax[1].set_yscale('linear')
ax[1].set_ylabel('Normalised SEDs')

##### Plotting the reconstructed SEDs
for i in range(20):
    plot_SED(ax[2], wave, SED_denormalised[i], color = colors[i])
ax[2].set_ylabel('Reconstructed SEDs')
##### Plotting the reconstruction residuals

plot_SED(ax[3], wave, SED_denormalised[i]/SED_masked[i], color = 'black', lw = 2, label = 'mean SED')
ax[3].set_xlabel(r'Wavelengths $[AA]$')
ax[3].set_ylim(0.999, 1.001)
ax[3].set_ylabel('Residuals')
ax[3].set_yscale('linear')
for j in range(3):
    ax[j].set_xticks([])

plt.savefig('Denormalised_SEDs.png')


