import numpy as np
import pandas as pd
import os
from scipy.spatial import Delaunay
from scipy.interpolate import RegularGridInterpolator
from scipy.constants import c

L_sun_erg_s = 3.828e33 # erg/s
pc_cm = 3.085677581e18 # cm
c_ang_per_s = c * 1e10   # speed of light in Angstrom/s

import fsps



def FSPS_initializer(zcontinuous = 1, sfh = 0, logzsol = 0.0, dust_emission = False, neb_emission = -1):
    # This function initialize the SPS code using the FSPS.StellarPopulation.
    # This SPS has 
    #           NO DUST
    #           NO IGM
    #           NO NEBULAE
    sp = fsps.StellarPopulation(compute_vega_mags = False,
                                zcontinuous = zcontinuous,
                                sfh = sfh,
                                logzsol = logzsol,
                                add_dust_emission = dust_emission,
                                add_neb_emission = neb_emission,)
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

def FSPS_SED_generator_SSP(param_names,
                       param_grid,
                       FSPS_initializer,
                       data_dir,
                       output_wave_filename = 'FSPS_wave_SSP.csv',
                       output_SED_template_filename = 'FSPS_SED_templates_SSP.csv'):
    # This generates galaxy SEDs based on the Stellar Population code
    # It calls a FSPS_initializer code, depending on the input
    #           1. NO DUST?
    #           2. NO IGM?
    #           3. NO NEBULAE?
    print ('\n')


    # Initialize the stellar population synthesis code
    sp = FSPS_initializer

    # Split 'tage' and 'logzsol' from the param_grid
    tage_grid = param_grid
    N_samples = len(tage_grid)

    # Generate a SED based on the value of 'tage' and 'logzsol'
    # And store it in an array SED_templates
    SED_templates = []
    for i, age in enumerate(tage_grid):
        # Update the metalicity of the stellar population

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

def reshape_coeff_grid(params, coeffs):
    """
    Reshape flattened parameter grid into 2D rectangular arrays.
    
    Parameters
    ----------
    params : (Nt*Nz, 2) array
        First column = tage, second = logz
    coeffs : (Nt*Nz, k) array
        PCA coefficients
        
    Returns
    -------
    tage_grid : (Nt,) unique sorted tage values
    logz_grid : (Nz,) unique sorted logz values
    coeff_grid : (Nt, Nz, k)
    """
    tage_vals = np.unique(params[:,0])
    logz_vals = np.unique(params[:,1])
    Nt, Nz = len(tage_vals), len(logz_vals)
    k = coeffs.shape[1]

    coeff_grid = np.zeros((Nt, Nz, k))
    for (t, z), c in zip(params, coeffs):
        i = np.where(tage_vals == t)[0][0]
        j = np.where(logz_vals == z)[0][0]
        coeff_grid[i,j,:] = c
    
    return tage_vals, logz_vals, coeff_grid
        
def SEDInterpolator_rect(tage_grid, logz_grid, coeff_grid, tage_new, logz_new):
    """
    Bilinear interpolation of PCA coefficients on a rectangular grid.

    Parameters
    ----------
    tage_grid : (Nt,) array of tage values (must be ascending)
    logz_grid : (Nz,) array of logZ values (must be ascending)
    coeff_grid : (Nt, Nz, k) PCA coefficients
    tage_new : float
    logz_new : float

    Returns
    -------
    coeff_new : (k,) interpolated coefficients
    """

    # find bracketing indices
    i1 = np.searchsorted(tage_grid, tage_new)
    j1 = np.searchsorted(logz_grid, logz_new)
    i0 = max(i1 - 1, 0)
    j0 = max(j1 - 1, 0)
    i1 = min(i1, len(tage_grid) - 1)
    j1 = min(j1, len(logz_grid) - 1)

    # interpolation weights
    if tage_grid[i1] == tage_grid[i0]:
        wt0, wt1 = 1.0, 0.0
    else:
        wt1 = (tage_new - tage_grid[i0]) / (tage_grid[i1] - tage_grid[i0])
        wt0 = 1.0 - wt1

    if logz_grid[j1] == logz_grid[j0]:
        wz0, wz1 = 1.0, 0.0
    else:
        wz1 = (logz_new - logz_grid[j0]) / (logz_grid[j1] - logz_grid[j0])
        wz0 = 1.0 - wz1

    # bilinear combination of 4 neighbors
    coeff_new = (
        wt0 * wz0 * coeff_grid[i0, j0, :]
        + wt0 * wz1 * coeff_grid[i0, j1, :]
        + wt1 * wz0 * coeff_grid[i1, j0, :]
        + wt1 * wz1 * coeff_grid[i1, j1, :]
    )

    return coeff_new

def spec_LsunA_to_f_lambda(spec_Lsun_per_A, distance_cm = 10*pc_cm):
    """
    Convert L (L_sun / Å) -> f_lambda (erg / s / cm^2 / Å) at given distance.
    distance_cm: distance in cm (e.g. 10*pc_cm for absolute mags; d_L for apparent mags).
    """
    L_lambda_erg_per_s_A = spec_Lsun_per_A * L_sun_erg_s           # erg / s / Å
    f_lambda = L_lambda_erg_per_s_A / (4.0 * np.pi * distance_cm**2) # erg / s / cm^2 / Å
    return f_lambda


def calc_ab_mag(wave, Lsun_per_A, filt_wave, filt_trans):
    """
    Compute AB magnitude for an arbitrary spectrum.

    wave, flux_lambda : arrays for spectrum [Angstrom], f_lambda in erg/s/cm^2/Å (or consistent units)
    filt_wave, filt_trans : arrays for filter [Angstrom], dimensionless transmission curve
    """

    flux_lambda = spec_LsunA_to_f_lambda(Lsun_per_A)
    # Interpolate spectrum onto filter grid
    spec_interp = np.interp(filt_wave, wave, flux_lambda, left=0, right=0)

    # Numerator and denominator of AB flux definition
    num = np.trapz(spec_interp * filt_trans * filt_wave, filt_wave)
    denom = np.trapz(filt_trans * c_ang_per_s / filt_wave, filt_wave)

    fnu = num / denom   # in erg/s/cm^2/Hz
    mag_ab = -2.5 * np.log10(fnu) - 48.6
    return mag_ab




## TESTING ZONE



