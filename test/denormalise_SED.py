import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from PCA_cal import denormalise_SED
from visualisations import plot_SED

###########################################################################################

# Define the data directory
data_dir = '../data/'

# Define the original data filenames
wave_filename = 'FSPS_wave_norm_std.csv'
SED_original_filename = 'FSPS_SED_masked.csv'

# Define the normalised data filenames
SED_norm_filename = 'FSPS_SED_norm_std.csv'
SED_norm_factor_filename = 'FSPS_SED_norm_factor_std.csv'
SED_norm_mean_filename = 'FSPS_SED_norm_mean_std.csv'


#######

# Read the files
wave = np.genfromtxt(data_dir + wave_filename)
SED_original = np.genfromtxt(data_dir + SED_original_filename)

SED_norm = np.genfromtxt(data_dir + SED_norm_filename)
SED_norm_factor = np.genfromtxt(data_dir + SED_norm_factor_filename)
SED_norm_mean = np.genfromtxt(data_dir + SED_norm_mean_filename)

############################################################################################

# Define the colors for plotting
colors = ['red', 'orange', 'yellow', 'green', 'blue']

############################################################################################
print ('Original SED dimensions')
print ('Wavelegnths')
print (wave.shape)
print ('SED ori')
print (SED_original.shape)

print ('\nNorm SED dimensions')
print ('SED norm')
print (SED_norm.shape)
print ('SED scale')
print (SED_norm_factor.shape)
print ('SED mean')
print (SED_norm_mean.shape)




#############################################################################################

SED_recon = denormalise_SED(SED_norm, SED_norm_factor, SED_norm_mean)

print ('\nRecon SED dimensions')
print ('SED recon')
print (SED_recon.shape)



#############################################################################################

fig, ax = plt.subplots(4, 1, figsize  =(12, 8), gridspec_kw = {'hspace': 0.01})

for i in range(len(SED_recon)):
    plot_SED(ax[0] ,wave, SED_original[i], color = colors[i])
    plot_SED(ax[1], wave, SED_recon[i], color = colors[i], ls = ':')
    plot_SED(ax[2], wave, SED_original[i], color = colors[i], lw = 5)
    plot_SED(ax[2], wave, SED_recon[i], color = 'grey')
    plot_SED(ax[3], wave, SED_recon[i]/SED_original[i], color = colors[i])

for j in range(2): ax[j].set_xticks([])
ax[3].set_yscale('linear')
ax[3].set_ylim(0.985, 1.015)
plt.savefig('SED_recon_comparison.png')
plt.show()
