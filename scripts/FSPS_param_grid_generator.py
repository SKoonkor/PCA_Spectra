import numpy as np
import pandas as pd
import os, sys

# Local python libraries
sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator



# Define the grid parameters

FSPS_param_names = ['tage',
                    'logzsol']

FSPS_param_ranges = np.array([[-5, np.log10(13.7)],
                              [-2.5, 0.5]])

FSPS_param_spacings = ['log',
                       'linear']

N_samples = 20000


# Show the LHS setup on the screen
print ('\n####################################################')
print ('Generating parameter grid for SPS code')
print ('Using parameter space as follows')
print ('Param:     , Min.:       , Max.:       , Spacing:')
for i in range(len(FSPS_param_names)):
    print (f"'{FSPS_param_names[i]:<{10}}, {FSPS_param_ranges[i][0]:3.9f}, {FSPS_param_ranges[i][1]:3.10f}, {FSPS_param_spacings[i]:<{10}}'")
print ('####################################################')



FSPS_param_values = LHS_generator(FSPS_param_names, 
                                  FSPS_param_ranges,
                                  FSPS_param_spacings,
                                  n_samples = N_samples)

output_dir = '/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/data/'
os.makedirs(output_dir, exist_ok = True) # Create a data/ directory 

# Save the parameter space to the data directory
param_grid_output_filename = 'FSPS_param_grid_LHS.csv'
np.savetxt(output_dir + param_grid_output_filename, FSPS_param_values)
print ('\nParameter values stored at: \n' + output_dir + param_grid_output_filename)


# Output the parameter names to the data directory
param_name_filename = 'FSPS_param_grid_names.txt'
with open(output_dir + param_name_filename, 'w') as f_name:
    for name in FSPS_param_names:
        f_name.write('%s\n' %name)
f_name.close()
print ('\nParameter names stored at: \n' + output_dir + param_name_filename)


# Output the parameter spacings to the data directory
param_spacing_filename = 'FSPS_param_grid_spacings.txt'
with open(output_dir + param_spacing_filename, 'w') as f_spacing:
    for spacing in FSPS_param_spacings:
        f_spacing.write('%s\n' %spacing)
f_spacing.close()
print ('\nParameter spacings stored at: \n' + output_dir + param_spacing_filename + '\n')
