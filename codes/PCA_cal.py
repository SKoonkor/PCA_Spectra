import numpy as np
import pandas as pd
import os, sys

from sklearn import preprocessing
from sklearn.decomposition import PCA as skPCA

def mask_SED(wave, SED, wave_min, wave_max):
    # This applies a wavelength masking on the SEDs
    # Only information between min_wave < wavelength < max_wave is selected

    # Make use of numpy for masking which is more computationally efficient. 
    wave_mask = np.where((wave >= wave_min)& (wave <= wave_max))[0]


    wave_output = wave[wave_mask]
    SED_output = SED[:, wave_mask]
    print ('\n*************************************************')
    print ('****************  Masking Process ***************')
    print ('*************************************************')
    print ('\nOriginal wavelengths')
    print ('shape: {}'.format(wave.shape))
    print ('Min: {:.2e} AA,  Max: {:.2e} AA\n\n'.format(min(wave), max(wave)))
    print ('Output wavelengths')
    print ('shape: {}'.format(wave_output.shape))
    print ('Min: {:.2e} AA, MAX: {:.2e} AA\n'.format(min(wave_output), max(wave_output)))  
    print ('*************************************************')
    return wave[wave_mask], SED[:, wave_mask]

def normalise_SED(wave, SED, norm_type = 'std', eps = 0):
    # This function is for normalising SEDs with 3 normalisation options
    #   1. The L1 norm
    #   2. The L2 norm (Dafault)
    #   3. The flux at 5500A.

    print ('\n#################################################')
    print ('################ Normalise SEDs #################')
    print ('#################################################')
    
    SED = np.asarray(SED, dtype = float)

    # If using 'L1' or 'L2' for normalisation (which require the preprocessing package from scikit-learn)
    if norm_type in ['l1', 'l2']:
        SED_norm,  SED_norm_factor = preprocessing.normalize(SED, norm = norm_type, return_norm = True)
        print ('\nNormalisation method: ', norm_type)
        
    # If using the 5500A flux
    elif norm_type == '5500A':
        # Fnid the index for 5500A (or the closest wavelength)
        index_5500A = np.where(wave == min(abs(wave - 5500)) + 5500)[0][0]

        SED_norm_factor = SED.transpose()[index_5500A] # Take the flux at 5500A of each spectrum
        SED_norm = SED/SED_norm_factor.reshape(-1, 1)
        SED_norm_mean = np.mean(SED_norm, axis = 0)
        print ('\nNormalisation method: ', norm_type)

    elif norm_type == 'std':
        # This methed centres the SEDs to its mean (so the normalised SEDs have zero mean)
        # Then it scales the centered SEDs to have std = 1.
        SED_norm_mean = np.nanmean(SED, axis = 0)
        sigma = np.nanstd(SED, axis = 0 , ddof = 1)
        SED_norm = (SED - SED_norm_mean)/(sigma + eps)
        SED_norm_factor = sigma
        print ('\nNormalisation method: ', norm_type) 

    else:
        raise Exception('The normalisation should be \'l1\', \'l2\', \'std\' or \'5500A\' only')

    print ('\n#################################################')
    return SED_norm, SED_norm_factor, SED_norm_mean




def PCA_SED(n_components = 20):
    # The SED are dimentionally reduced to n_components
    # Instead of using the wavelengths as the vector for the spectra,
    # The spectra are projected onto the principal components (PCs)
    # The SED instead of represented as F_ij = f_ij * I_j,
    # They are described as S = c_ik* PC
    # notation: spectrum i-th, wavelentgh j-th.

    return 

def denormalise_SED(SED_norm, SED_std, SED_mean):
    # This de-normalise the normalised SEDs back to the original
    # This method can be applied at the very end as a postprocessing step
    # It is define (as simple as):
    #
    # SED_recon = SED_norm * SED_std + SED_mean.
    #

    print ('\n#################################################')
    print ('################ denormise SED ##################')
    print ('#################################################')
    return SED_norm*SED_std + SED_mean

def PCA_recon():
    # This reconstruct the PCA SEDs back to the 'original SEDs'
    

    return 



###### TEST ######
