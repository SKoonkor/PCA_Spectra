import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from visualisations import plot_SED

####
# Define filenames for inputs
data_dir = '../data/'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_filename = 'FSPS_SED_norm_std.csv'
SED_mean_filename = 'FSPS_SED_norm_mean_std.csv'
SED_mean_norm_factor_std = 'FSPS_SED_norm_factor_std.csv'

wave_original_filename = 'FSPS_wave.csv'
SED_original_filename = 'FSPS_SED_templates.csv'
SED_masked_filename = 'FSPS_SED_masked.csv'

# Read the input
wave = np.genfromtxt(data_dir + wave_filename)
SED = np.genfromtxt(data_dir + SED_filename)
SED_mean = np.genfromtxt(data_dir + SED_mean_filename)

wave_original = np.genfromtxt(data_dir + wave_original_filename)
SED_original = np.genfromtxt(data_dir + SED_original_filename)
SED_masked = np.genfromtxt(data_dir + SED_masked_filename)

# Define the colors
colors = ['red', 'orange', 'yellow', 'green', 'blue']


##### plotting the normalised SED
fig, ax  = plt.subplots(figsize = (12, 8))
for i in range(len(SED)):
    plot_SED(ax, wave, SED[i], yscale = 'linear', color = colors[i])

# plot_SED(ax, wave, SED_mean, color = 'black', lw = 5, yscale = 'linear'
ax.set_title('The normalised SEDs')
ax.set_xlabel('wavelength [AA]')
ax.set_ylabel('$F_{norm}$')
plt.savefig('Normalised_SEDs.png')

##### Plotting the original SEDs
fig2, ax2 = plt.subplots(figsize = (12, 8))

for i in range(len(SED_original)):
    plot_SED(ax2, wave, SED_masked[i], color = colors[i], lw = 3)

plot_SED(ax2, wave, SED_mean, color = 'black', lw = 3)

ax2.set_ylabel('F')
ax2.set_xlabel('wavelengths [AA]')
ax2.set_title('The original SEDs')
ax2.set_ylim(1e-6, 1e2)
plt.savefig('Original_SEDs.png')


##### Plotting the masked SEDs

