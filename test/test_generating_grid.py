import numpy as np
import pandas as pd
import os, sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator


FSPS_param_names = ['tage', 
                    'logzsol']
FSPS_param_ranges = [[1e-4, 13.7],
                     [-.5, 2.5]]
FSPS_param_spacings = ['linear',
                        'linear']



for name, param_range, param_spacing in zip(FSPS_param_names, FSPS_param_ranges, FSPS_param_spacings):
    print (name, param_range, param_spacing)



LHS_generator(FSPS_param_names, FSPS_param_ranges, FSPS_param_spacings)
