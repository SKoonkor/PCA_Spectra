import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys


sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from FSPS_paramgrid import LHS_generator
from visualisations import plot_LatinHypercube


FSPS_param_names = ['tage', 
                    'logzsol']
FSPS_param_ranges = np.array([[-4, np.log10(13.7)],
                     [-.5, 2.5]])
FSPS_param_spacings = ['log',
                        'linear']


for name, param_range, param_spacing in zip(FSPS_param_names, FSPS_param_ranges, FSPS_param_spacings):
    print (name, param_range, param_spacing)



FSPS_param_values = LHS_generator(FSPS_param_names, FSPS_param_ranges, FSPS_param_spacings, n_samples = 1100)
 


fig, ax = plot_LatinHypercube(FSPS_param_names, FSPS_param_values, FSPS_param_spacings, index_x = 0, index_y = 1)
plt.show()
