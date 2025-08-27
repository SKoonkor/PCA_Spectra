import numpy as np
import pandas as pd
import os

import fsps



def FSPS_initializer():
    # This function initialize the SPS code using the FSPS.StellarPopulation.
    # This SPS has 
    #           NO DUST
    #           NO IGM
    #           NO NEBULAE
    sp = fsps.StellarPopulation(compute_vega_mags = False,
                                zcontinuous = 1,
                                sfh = 0,
                                logzsol = 0.0,
                                add_dust_emission = False,
                                add_neb_emission = -1,)

    return sp




def FSPS_SED_generator(FSPS_initializer = FSPS_initializer()):
    # This generates galaxy SEDs based on the Stellar Population code
    # It calls a FSPS_initializer code, depending on the input
    #           1. NO DUST?
    #           2. NO IGM?
    #           3. NO NEBULAE?

    sp = FSPS_initializer

    return sp





## TESTING ZONE

sp = FSPS_SED_generator()

print (sp)
