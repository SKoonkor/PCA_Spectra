import numpy as np
import pandas as pd
import os, sys

# Local python libraries
sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator



# Define the grid parameters

FSPS_param_names = ['tage',
                    'logzsol']

FSPS_param_ranges = np.array([[-4, np.log10(13.7)],
                              [-.5, 2.5]])

FSPS_param_spacings = ['log',
                       'linear']


FSPS_param_values = LHS_generator(FSPS_param_names, 
                                  FSPS_param_ranges,
                                  FSPS_param_spacings,
                                  n_samples = 1100)

output_dir = '/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/data/'


# Save the parameter space to the data directory
param_grid_output_filename = 'FSPS_param_grid_LHS.csv'
np.savetxt(output_dir + param_grid_output_filename, FSPS_param_values)
print ('The parameter space has been saved to: ' + output_dir + param_grid_output_filename)


# Output the parameter names to the data directory
param_name_filename = 'FSPS_param_grid_names.txt'
with open(output_dir + param_name_filename, 'w') as f_name:
    for name in FSPS_param_names:
        f_name.write('%s\n' %name)
f_name.close()
print ('The parameter names have been stored at: ' + output_dir + param_name_filename)


# Output the parameter spacings to the data directory
param_spacing_filename = 'FSPS_param_grid_spacings.txt'
with open(output_dir + param_spacing_filename, 'w') as f_spacing:
    for spacing in FSPS_param_spacings:
        f_spacing.write('%s\n' %spacing)
f_spacing.close()
print ('The parameter spacings have been stored at: ' + output_dir + param_spacing_filename)
