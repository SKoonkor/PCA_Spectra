import numpy as np
import pandas as pd
import sys, os

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')


from PCA_cal import PCA_SED, reconstruct_PCA_SED


############################################################################################

##### Define the input filenames
data_dir = '../data/'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_norm_filename = 'FSPS_SED_norm_std.csv'
SED_mean_filename = 'FSPS_SED_norm_mean_std.csv'
SED_factor_filename = 'FSPS_SED_norm_factor_std.csv'
SED_original_filename = 'FSPS_SED_masked.csv'

##### Read the input files
wave = np.genfromtxt(data_dir + wave_filename)
SED_norm = np.genfromtxt(data_dir + SED_norm_filename)
SED_mean = np.genfromtxt(data_dir + SED_mean_filename)
SED_factor = np.genfromtxt(data_dir + SED_factor_filename)
SED_original = np.genfromtxt(data_dir + SED_original_filename)

###########################################################################################

PCA_coeffs, PCA_PCs = PCA_SED(SED_norm, n_components = 200)


print ('Coeffs shape: ', PCA_coeffs.shape)
print ('PCs shape: ', PCA_PCs.shape)
print ('Length PCs: ', len(PCA_PCs))
print ('Length Coeffs: ', len(PCA_coeffs))
print ('Length PCs[i]: ', len(PCA_PCs[0]))
print ('Length Coeffs[i]: ', len(PCA_coeffs[0]))
