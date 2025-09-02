import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fsps
import sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_calculations import SEDInterpolator, FSPS_initializer
from visualisations import plot_SED
from PCA_cal import denormalise_SED, reconstruct_PCA_SED


###########################################################################################

# Define the data directory
data_dir = '../data/'
param_grid_filename = 'FSPS_param_grid_LHS.csv'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_norm_factor_std_filename = 'FSPS_SED_norm_factor_std.csv'
SED_norm_mean_std_filename = 'FSPS_SED_norm_mean_std.csv'

# Define the PCA output directory
PCA_output_dir = 'outputs/'
PCA_PC_filename = 'PCA_PCs.csv'
PCA_coeff_filename = 'PCA_coeffs.csv'

###########################################################################################

# Read the SED template files
wave_masked = np.genfromtxt(data_dir + wave_filename)
SED_norm_factor = np.genfromtxt(data_dir + SED_norm_factor_std_filename)
SED_norm_mean = np.genfromtxt(data_dir + SED_norm_mean_std_filename)

# Read the interpolation grid files
param_grid = np.genfromtxt(data_dir + param_grid_filename)
tage_grid, logzsol_grid = param_grid[:, 0], param_grid[:, 1]
PCA_coeffs = np.genfromtxt(PCA_output_dir + PCA_coeff_filename)
PCA_PCs = np.genfromtxt(PCA_output_dir + PCA_PC_filename)

###########################################################################################

print ('tage grid shape: ', tage_grid.shape)
print ('coeffs shape: ', PCA_coeffs.shape)

###########################################################################################

print ('Define the new parameters for interpolating')
tage_new, logzsol_new = 1.2, 0.01

print ('tage_new: ', tage_new)
print ('logzsol_new: ', logzsol_new)


interp_SED = SEDInterpolator(tage_grid, logzsol_grid, PCA_coeffs)

coeff_new = interp_SED.interpolate(tage_new, logzsol_new)

print ('\n INTERPOLATION DONE')

###########################################################################################
# Reconstruct the interpolated coeffs back into SED form
print ('\n Reconstructing the PCA SED')
SED_recon = reconstruct_PCA_SED(coeff_new, PCA_PCs, n_components = 150)

SED_recon_denorm = denormalise_SED(SED_recon, SED_norm_factor, SED_norm_mean)


###########################################################################################

sp = FSPS_initializer()
sp.params['logzsol'] = logzsol_new

wave_FSPS, SED_FSPS = sp.get_spectrum(tage = tage_new, peraa = True)
masked_index = np.where((wave_FSPS >= min(wave_masked)) & (wave_FSPS <= max(wave_masked)))[0]
SED_FSPS_masked = SED_FSPS[masked_index]


print (len(SED_norm_factor))
print (len(SED_FSPS_masked))

fig, ax = plt.subplots(2, 1, figsize = (24, 12), gridspec_kw = {'hspace': 0.01})
# Plot the reconstructed SED
plot_SED(ax[0], wave_masked, SED_recon_denorm) 
plot_SED(ax[0], wave_FSPS, SED_FSPS, color = 'grey', lw = 3)

plot_SED(ax[1], wave_masked, SED_recon_denorm/SED_FSPS_masked)

ax[1].axhline(1.1, color = 'grey')
ax[1].axhline(0.9, color = 'grey')
ax[1].axhline(1.05, color = 'green')
ax[1].axhline(0.95, color = 'green')
ax[1].axhline(1, color = 'black')

ax[0].set_title('MY BABY!')
ax[0].set_ylabel('F_lambda')
ax[0].set_xlabel(r'wavelength [$\AA$]')
for i in range(2): ax[i].set_xlim(min(wave_FSPS), max(wave_FSPS))
ax[1].set_ylabel('Fractional F')
ax[1].set_yscale('linear')
ax[1].set_ylim(0.5, 1.5)
plt.savefig('interp_SED.png')
plt.show()
