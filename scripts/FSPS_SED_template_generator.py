import numpy as np
import pandas as pd
import os, sys

# Local python libraries
sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator
from FSPS_calculations import FSPS_initializer, FSPS_SED_generator


# Script stars here

# Define filenames and directory for storing SED template information
data_dir = '/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/data/'
os.makedirs(data_dir, exist_ok = True)

# Input files
param_name_filename = 'FSPS_param_grid_names.txt'
param_grid_filename = 'FSPS_param_grid_LHS.csv'

# Output files
output_wave_filename = 'FSPS_wave.csv'
output_SED_template_filename = 'FSPS_SED_templates.csv'

###########################################################################################
# Main calculations

# Read the param names into a list
param_grid_names = []
with open(data_dir + param_name_filename, 'r') as f_param_names:
    for name in f_param_names:
        param_grid_names.append(name)
f_param_names.close()

# Read the param grid values into an array
param_grid_values = np.genfromtxt(data_dir + param_grid_filename)

# Generate the SED templates
FSPS_SED_generator(param_names = param_grid_names,
                   param_grid = param_grid_values,
                   FSPS_initializer = FSPS_initializer(),
                   data_dir = data_dir,
                   output_wave_filename = output_wave_filename,
                   output_SED_template_filename = output_SED_template_filename)

print ('SED TEMPLATE GENERATION COMPLETE')
