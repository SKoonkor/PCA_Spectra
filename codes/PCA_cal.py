import numpy as np
import pandas as pd
import os, sys


def nornalise_SEDs(wave, SED, norm_type = 'l2'):
    # This function is for normalising SEDs with 3 normalisation options
    #   1. The L1 norm
    #   2. The L2 norm (Dafault)
    #   3. The flux at 5500A.

    # If using 'L1' or 'L2' for normalisation (which require the preprocessing package from scikit-learn)
    if norm_type in ['l1', 'l2']:
        SED_norm = SED_norm_factor = preprocessing.normalize(SED, norm = norm_type, return_norm = True)

    # If using the 5500A flux
    elif norm_type == '5500A':
        # Fnid the index for 5500A (or the closest wavelength)
        index_5500A = np.where(wave == min(abs(wave - 5500)) + 5500)[0][0]

        SED_norm_factor = SED.transpose()[index_5500A] # Take the flux at 5500A of each spectrum
        SED_norm = SED/SED_norm_factor.reshape(-1, 1)

    else:
        raise Exception('The normalisation should be \'l1\', \'l2\', or \'5500A\' only')

    SED_norm_mean = np.mean(SED_norm, axis = 0)

    return SED_norm, SED_norm_factor, SED_norm_mean
