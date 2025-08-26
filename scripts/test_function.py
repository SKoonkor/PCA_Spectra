import sys, os
import numpy as np



sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator





param_names = ['test_param']

param_ranges = np.array([
    [2, 5]
    ])

param_spacings = ['liner']


print (LHS_generator(param_names = param_names, param_ranges = param_ranges, param_spacings = param_spacings, n_samples
                     = 20))
