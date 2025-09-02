import numpy as np
import pandas as pd
import os, sys

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from PCA_cal import PCA_SED

############################################################################################

#### Define the input filenames
data_dir  = '../data/'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_norm_filename = 'FSPS_SED_norm_std.csv'


#### Read the input files
wave = np.genfromtxt(data_dir + wave_filename)
SED = np.genfromtxt(data_dir + SED_norm_filename)

###########################################################################################

#### Define the output directory
output_dir = '../outputs/'

#### RUN PCA on SEDs

PCA_coeffs, PCA_PCs = PCA_SED(SED, output_dir = output_dir, n_components = 200)

print ('coeff shape: ', PCA_coeffs.shape)
print ('coeff len: ', len(PCA_coeffs))
print ('PCs shape: ', PCA_PCs.shape)




print ('\n PCA SCRIPT COMPLETE')
