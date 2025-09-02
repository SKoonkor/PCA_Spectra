import numpy as np
import pandas as pd
import os, sys

# Local python libraries
sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')
from PCA_cal import mask_SED, normalise_SED 


############################################################################################
#### Script starts here
data_dir = '../data/'
os.makedirs(data_dir, exist_ok = True)
wave_filename = 'FSPS_wave.csv'
SED_filename = 'FSPS_SED_templates.csv'

# Define the names of the input files
wave = np.genfromtxt(data_dir + wave_filename)
SED = np.genfromtxt(data_dir + SED_filename)

# Mask the SED by applying the wavelength range
wave_min = 1e3  
wave_max = 5e4

# Define the masked data filenames
wave_masked_filename = 'FSPS_wave_masked.csv'
SED_masked_filename = 'FSPS_SED_masked.csv'

# Export the masked data
wave_masked, SED_masked = mask_SED(wave, SED, wave_min = wave_min, wave_max = wave_max)
np.savetxt(data_dir + wave_masked_filename, wave_masked)
print ('\nWavelengths (masked) stored at: ')
print (data_dir + wave_masked_filename)
np.savetxt(data_dir + SED_masked_filename, SED_masked)
print ('\nSEDs (masked) stored at: ')
print (data_dir + SED_masked_filename)

# normalise the SEDs
norm_type = 'std'
SED_norm, SED_norm_factor, SED_norm_mean = normalise_SED(wave_masked, SED_masked, norm_type = norm_type)


######
# Output the normalised SEDs and wavelengths information

#Define the output filenames for normalised data
wave_norm_output_filename = 'FSPS_wave_norm_{}.csv'.format(norm_type)
SED_norm_output_filename = 'FSPS_SED_norm_{}.csv'.format(norm_type)
SED_norm_factor_output_filename = 'FSPS_SED_norm_factor_{}.csv'.format(norm_type)
SED_norm_mean_output_filename = 'FSPS_SED_norm_mean_{}.csv'.format(norm_type)

# Store the files
np.savetxt(data_dir + wave_norm_output_filename, wave_masked)
print ('\nWavelengths (normalised) stored at:')
print (data_dir + wave_norm_output_filename)

np.savetxt(data_dir + SED_norm_output_filename, SED_norm)
print ('\nSEDs (normalised) stored at:')
print (data_dir + SED_norm_output_filename)

np.savetxt(data_dir + SED_norm_factor_output_filename, SED_norm_factor)
print ('\nSED norm factor stored at:')
print (data_dir + SED_norm_factor_output_filename)

np.savetxt(data_dir + SED_norm_mean_output_filename, SED_norm_mean)
print ('\nSED mean (normalised) stored at:')
print (data_dir + SED_norm_mean_output_filename)

print ('\nSED NORMALISATION COMPLETE')



# print (wave.shape)
# print (SED.shape)
# print (' ')
# print (wave_masked.shape)
# print (SED_masked.shape)
