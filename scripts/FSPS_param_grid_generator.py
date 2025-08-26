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

np.savetxt(output_dir+'FSPS_param_grid_LHS.csv', FSPS_param_values)


with open(output_dir+'FSPS_param_grid_names.txt', 'w') as f:
    for name, space in zip(FSPS_param_names, FSPS_param_spacings):
        f.write('%s, %s \n' % (name, space))
f.close()
