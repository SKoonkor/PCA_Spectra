import numpy as np
import pandas as pd
import os
from scipy.spatial import Delaunay

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




def FSPS_SED_generator(param_names,
                       param_grid,
                       FSPS_initializer,
                       data_dir,
                       output_wave_filename = 'FSPS_wave.csv',
                       output_SED_template_filename = 'FSPS_SED_templates.csv'):
    # This generates galaxy SEDs based on the Stellar Population code
    # It calls a FSPS_initializer code, depending on the input
    #           1. NO DUST?
    #           2. NO IGM?
    #           3. NO NEBULAE?
    print ('\n')


    # Initialize the stellar population synthesis code
    sp = FSPS_initializer

    # Split 'tage' and 'logzsol' from the param_grid
    tage_grid, logzsol_grid = param_grid[:, 0], param_grid[:, 1]
    N_samples = len(tage_grid)

    # Generate a SED based on the value of 'tage' and 'logzsol'
    # And store it in an array SED_templates
    SED_templates = []
    for i, (age, logzsol) in enumerate(zip(tage_grid, logzsol_grid)):
        # Update the metalicity of the stellar population
        sp.params['logzsol'] = logzsol

        # Calculate the spectrum of the population at the age 'tage'
        # 'peraa' kw is True for L_sun/AA, False for L_sun/Hz
        wave, SED = sp.get_spectrum(tage = age, peraa = True)

        # Add it to the template array
        SED_templates.append(SED)


        # Print the progress on screen
        Progress = int(i*100/N_samples)
        if Progress % 2 == 0:
            print (f'\rGenerating SED {Progress}% complete', end = '', flush = True)
    SED_templates = np.array(SED_templates)


    # Output the SED templates and its corresponding resolution (wavelength information)
    np.savetxt(data_dir + output_SED_template_filename, SED_templates)
    print ('\n')
    print ('SED templates stored at: ' + data_dir + output_SED_template_filename + '\n')
    np.savetxt(data_dir + output_wave_filename, wave)
    print ('Wavelengths stored at: ' + data_dir + output_wave_filename + '\n')


class SEDInterpolator:

    def __init__(self, tage, logzsol, coeffs):
        """
        Parameters
        __________
        tage : array (N,)
            Stellar ages for each template
        logzsol : array (N,)
            Log metallicities for each template
        coeffs : array (N, k)
            PCA coefficients for each template (k = number of PCs)
        """
        self.points = np.column_stack([tage, logzsol])
        self.coeffs = coeffs
        self.tri = Delaunay(self.points)


    def interpolate(self, tage_new, logzsol_new):
        """
        Interpolate PCA coefficients for new parameters.

        Parameters
        __________
        tage_new : float
        logzsol_new : float

        Returns
        _______
        coeff_interp : array (k,)
            Interpolated PCA coefficients

        """
        p = np.array([tage_new, logzsol_new]) # point_interp where we want to interpolate
        simplex = self.tri.find_simplex(p)
        
        if simplex == -1:
            raise ValueError('Point outside convex hull of training grid.')

        vertices = self.tri.simplices[simplex]
        X = self.tri.transform[simplex, :2]
        Y = p - self.tri.transform[simplex, 2]

        bary = np.dot(X, Y)
        bary = np.append(bary, 1 - bary.sum())

        coeff_interp = np.dot(bary, self.coeffs[vertices])

        return coeff_interp
        







## TESTING ZONE



