import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys, os

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from visualisations import plot_SED
from PCA_cal import reconstruct_PCA_SED, denormalise_SED

############################################################################################
 
# Define input file filenames
data_dir = '../data/'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_norm_filename = 'FSPS_SED_norm_std.csv'
SED_norm_mean_filename = 'FSPS_SED_norm_mean_std.csv'
SED_norm_factor_filename = 'FSPS_SED_norm_factor_std.csv'

SED_masked_filename = 'FSPS_SED_masked.csv'

# Read the input files
wave = np.genfromtxt(data_dir + wave_filename)
SED_norm = np.genfromtxt(data_dir + SED_norm_filename)
SED_norm_mean = np.genfromtxt(data_dir + SED_norm_mean_filename)
SED_norm_factor = np.genfromtxt(data_dir + SED_norm_factor_filename)

SED_masked = np.genfromtxt(data_dir + SED_masked_filename)
###########################################################################################
# Define the PCA output directory
PCA_output_dir = '../outputs/'
PCA_PCs_filename = 'PCA_PCs.csv'
PCA_coeffs_filename = 'PCA_coeffs.csv'

###
# Read the PCA outputs
PCA_PCs = np.genfromtxt(PCA_output_dir + PCA_PCs_filename)
PCA_coeffs = np.genfromtxt(PCA_output_dir + PCA_coeffs_filename)

##########################################################################################


n_PC_recon = 50 # the number of PCs used for reconstructing the SED

SED_reconstructed = []

for i in range(len(PCA_coeffs)): # To reconstruct all the SEDs in the template
    SED_reconstructed.append(reconstruct_PCA_SED(PCA_coeffs[i], PCA_PCs, n_components = n_PC_recon))

SED_reconstructed = np.asarray(SED_reconstructed)


#########################################################################################

SED_recon_denorm = denormalise_SED(SED_reconstructed, SED_norm_factor, SED_norm_mean)





##########################################################################################
n_SED_plot = 5*5
colors = ['red', 'orange', 'yellow', 'green', 'blue']*int(n_SED_plot/5)

#### Plotting SEDs

fig, ax = plt.subplots(4, 1, figsize = (12, 8), gridspec_kw = {'hspace': 0.01})

#### Plotting norm SEDs
for i in range(n_SED_plot):
    plot_SED(ax[0], wave, SED_norm[i], color = colors[i])
ax[0].set_yscale('linear')
ax[0].set_ylabel('norm SEDs')

#### Plotted the reconstructed SEDs

for i in range(n_SED_plot):
    plot_SED(ax[1], wave, SED_reconstructed[i], color = colors[i])
ax[1].set_yscale('linear')
ax[1].set_ylabel('reconstructed SEDs')

for i in range(n_SED_plot):
    plot_SED(ax[2], wave, 
             SED_reconstructed[i]/SED_norm[i], color = colors[i])
ax[2].set_yscale('linear')
ax[2].set_ylabel('Frac. errors')

for i in range(n_SED_plot):
    plot_SED(ax[3], wave,  SED_recon_denorm[i]/SED_masked[i], color = colors[i])
ax[3].set_yscale('linear')
ax[3].set_ylabel('Frac. errors')
ax[3].axvline(1216)

for j in range(3):
    ax[j].set_xticks([])
ax[3].set_xlabel(r'Wavelengths [$\AA$]')

plt.show()
