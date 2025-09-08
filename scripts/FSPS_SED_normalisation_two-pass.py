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

# Mask the SED by applying the wavelength range
wave_min = 1e2  
wave_max = 1e7

# Define the masked data filenames
wave_masked_filename = 'FSPS_wave_masked.csv'
SED_masked_filename = 'FSPS_SED_masked.csv'

# Export the masked data


mask_SED(data_dir + wave_filename, data_dir + SED_filename, wave_min = wave_min, wave_max = wave_max, 
                  output_filename = data_dir + SED_masked_filename, chunk_size = 1000)


wave_masked = np.genfromtxt(data_dir + wave_masked_filename)








# normalise the SEDs
norm_type = 'std'

######

# Define the normalised SEDs using the two-pass algorithm
wave_norm_output_filename = 'FSPS_wave_norm_{}.csv'.format(norm_type)
SED_norm_output_filename = 'FSPS_SED_norm_{}.csv'.format(norm_type)
SED_norm_mean_output_filename = 'FSPS_SED_norm_mean_{}.csv'.format(norm_type)
SED_norm_factor_output_filename = 'FSPS_SED_norm_factor_{}.csv'.format(norm_type)

SED_factor, SED_mean = normalise_SED(wave_masked, data_dir + SED_masked_filename, 
                                                             norm_type = norm_type, chunk_size = 1000, output_filename =
                                                             data_dir + SED_norm_output_filename)




# Store the files
np.savetxt(data_dir + SED_norm_factor_output_filename, SED_factor)
np.savetxt(data_dir + SED_norm_mean_output_filename, SED_mean)
np.savetxt(data_dir + wave_norm_output_filename, wave_masked)
print ('\nWavelengths (normalised) stored at:')
print (data_dir + wave_norm_output_filename)

print ('\nSEDs (normalised) stored at:')
print (data_dir + SED_norm_output_filename)

print ('\nSED norm factor stored at:')
print (data_dir + SED_norm_factor_output_filename)

print ('\nSED mean (normalised) stored at:')
print (data_dir + SED_norm_mean_output_filename)









print ('\nSED NORMALISATION COMPLETE')



# print (wave.shape)
# print (SED.shape)
# print (' ')
# print (wave_masked.shape)
# print (SED_masked.shape)
