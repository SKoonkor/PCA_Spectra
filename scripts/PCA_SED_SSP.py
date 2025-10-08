import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import os, sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')
from visualisations import plot_SED
from PCA_cal import denormalise_SED, reconstruct_PCA_SED

# Define essential directories
data_dir = '../data/'
output_dir = '../outputs/'

# Define data/output filenames
wave_original_filename = data_dir + 'FSPS_wave.csv'
SED_original_filename = data_dir + 'FSPS_SED_templates.csv'

SPS_grid_param_filename = data_dir + 'FSPS_param_grid_LHS.csv'
wave_norm_filename = data_dir + 'FSPS_wave_masked.csv'
SED_norm_filename = data_dir + 'FSPS_SED_norm_std.csv'
SED_norm_factor_filename = data_dir + 'FSPS_SED_norm_factor_std.csv'
SED_norm_mean_filename = data_dir + 'FSPS_SED_norm_mean_std.csv'
PC_filename = output_dir + 'PCA_PCs.csv'
coeff_filename = output_dir + 'PCA_coeffs.csv'
# Reading the normalised SEDs and PCs
wave_original = np.genfromtxt(wave_original_filename)
wave_norm = np.genfromtxt(wave_norm_filename)

PCs = np.genfromtxt(PC_filename)
coeffs = np.genfromtxt(coeff_filename)

SPS_grid_param = np.genfromtxt(SPS_grid_param_filename)

# Reading the SED templates as chunks
'''
Due to the large file size, it's more fleasible to read
the spectra in chunks instead of all together
'''

chunk_size = 10
####
# The original template
SED_original = pd.read_csv(SED_original_filename, 
                              header=None,
                              sep='\s+',
                              nrows = chunk_size).to_numpy()
####
# The normalised spectra
SED_norm = pd.read_csv(SED_norm_filename, 
                       header=None,
                       sep='\s+',
                       nrows = chunk_size).to_numpy()
SED_norm_factor = np.genfromtxt(SED_norm_factor_filename)
SED_norm_mean = np.genfromtxt(SED_norm_mean_filename)

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
ax[0].set_title('The first {} Principal Components of the SSP Spectra'.format(n_PCA))
plt.savefig(output_dir + '{}_SSP_PCs.png'.format(n_PCA))
plt.show()

# Plot the de-normalised SEDs
fig, ax = plt.subplots(chunk_size, 1, figsize = (7, chunk_size))
denorm_SED = denormalise_SED(SED_norm, SED_norm_factor, SED_norm_mean)
for i in range(chunk_size):
    plot_SED(ax[i], wave_norm, denorm_SED[i])
    plot_SED(ax[i], wave_original, SED_original[i], color = 'grey', ls = ':')

    ax[i].set_ylabel(r'$F_\lambda$')
    ax[i].set_ylim(1e-10, 1e1)
    
    if i == chunk_size - 1:
        ax[i].set_xlabel('wavelength [A]')
    else:
        ax[i].set_xticks([])
ax[0].set_title('The first {} SSP Spectra'.format(chunk_size))
plt.show()



####
# PCA SED re-contruction

n_PC_recon = 50
n_spec_recon = 10
print (coeffs.shape)
print (coeffs[0])

SED_recon_norm = []

for i in range(n_spec_recon):
    SED_recon_norm.append(reconstruct_PCA_SED(coeffs[0:n_spec_recon][i], PCs, n_PC_recon))
SED_recon_norm = np.array(SED_recon_norm)

SED_recon = denormalise_SED(SED_recon_norm, SED_norm_factor, SED_norm_mean)

fig, ax = plt.subplots(n_spec_recon, 1, figsize = (7, n_spec_recon))

for i in range(n_spec_recon):
    plot_SED(ax[i], wave_original, SED_original[i], color = 'grey')
    plot_SED(ax[i], wave_norm, SED_recon[i])
    

    ax[i].set_ylabel(r'$F_\lambda$')
    ax[i].set_ylim(1e-10, 1e1)
    ax[i].text(1e4, 1e-2, '{:.3f}'.format(SPS_grid_param[i]), fontsize = 12)

    if i == n_spec_recon - 1:
        ax[i].set_xlabel('Wavelength [A]')
    else:
        ax[i].set_xticks([])
ax[0].set_title('The SED reconstructed using the PCA')

plt.show()

    



