import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import os, sys

from astropy.cosmology import FlatLambdaCDM, z_at_value
from scipy.interpolate import interp1d


# Set up the COSMOLOGY
cosmo = FlatLambdaCDM(H0 = 67.8, Om0 = 0.307, Ob0 = 0.0483)


# Functions for cosmological calculations

def z_to_a(z):
    """
    Redshift to scale factor converter
    """
    return 1/(1+z)

def a_to_z(a):
    """
    Scale factor to redshift converter
    """
    return (1-a)/a

def age_to_z(age_Gyr):
    """
    Age to redshift converter
    Age (Gyr) -> redshift
    """

    if isinstance(age_Gyr, float) or isinstance(age_Gyr, int):
        return z_at_value(cosmo.age, age_Gyr * u.Gyr).value
    else:
        z = np.array([z_at_value(cosmo.age, a) for a in age_Gyr * u.Gyr])
        return z.values

def z_to_age(z):
    """
    Redshift to age converter
    Redshift -> age of the Universe (Gyr)
    """
    return cosmo.age(z).to(u.Gyr).value

def Universe_grid_generator(z_init = 127, z_final = 0, n_snapshots = 271):
    """
    To reproduce the snapshot grid for a simulation used
    Note: this expects the simulation to be uniform in scale factor
    Hence, need to work out the redshift of each snapshot and the Universe age herein.
    """
    a_init, a_final = z_to_a(z_init), z_to_a(z_final)

    a_grid = np.linspace(a_init, a_final, n_snapshots, endpoint = True)
    z_grid = a_to_z(a_grid)

    return z_to_age(z_grid)



def SFR_interp(Universe_age, SFH_age_grid_table, SFH_str_grid_table):

    """
    This function interpolates the star formation rate
    from the input star formation history using a linear
    interpolation
    """
    f = interp1d(SFH_age_grid_table, SFH_str_grid_table,
                 bounds_error = False,
                 fill_value = 0.0)

    return f(Universe_age)

def mass_formed(SFH_time_grid, SFH_sfr_grid, t_universe, bin_position = "middle"):
    """
    Calculate the stellar mass formed for each time step from the input star formation history.
    Parameters
    There are three time bin approaches
        -Front
    ........|<------|<------|<------|<------|<------|<------|--------
        -Middle
    ........|...!<--|---!<--|---!<--|---!<--|---!<--|---!<--|--------
        -Back
    ........|.......|<------|<------|<------|<------|<------|<-------
    
    tBB-----ti------ti+1----ti+2----|-------|-------tf-1----tf-------Tuniverse

        
    
    Time grid since big bang. 
    Each    | means the time at the snapshot output,
            < means when the SSP age is considered (Front, middle, or back of each bin),
            ! means the middle point between two time steps.
    The schematic pattern shows that the 
    
    
    ----------
    SFH_time: array_like
        Cosmic times (increasing) where the SFR is known. 
        Technically this should be the age of the universe at each snapshot of the GALFORM outputs
    SFH_str: array_like
        The star formation values corresponding to the SFH_time
        
        
    Returns
    -------
    mass_formed: array_like [same length as len(SFH_time) - 1]
        The amount of stellar mass formed between each time step
    """
    # Turn input time grid and sfr grid into arrays
    SFH_time_grid = np.asarray(SFH_time_grid)
    SFH_sfr_grid = np.asarray(SFH_sfr_grid)

    # Make sure that the SFH grid is ascending cosmic time
    order = np.argsort(SFH_time_grid)
    SFH_time_grid = SFH_time_grid[order]
    SFH_sfr_grid = SFH_sfr_grid[order]

    # Only select the time that are young than the input t_universe
    time_cut = np.where(SFH_time_grid < t_universe)[0]

    # Interpolate the last time step
    SFR_last = SFR_interp(t_universe, SFH_time_grid, SFH_sfr_grid)
    SFH_time_output = np.append(SFH_time_grid[time_cut], t_universe)
    SFH_sfr_output = np.append(SFH_sfr_grid[time_cut], SFR_last)

    # Calculate the mass formed between each bin
    mass_formed = []
    for i in range(len(SFH_time_output)-1):
        t = SFH_time_output[i+1] - SFH_time_output[i]
        sfr_avg = (SFH_sfr_output[i+1]+SFH_sfr_output[i])/2
        mass_formed.append(t*sfr_avg)
    mass_formed = np.array(mass_formed)

    if bin_position == 'front':
        mass_formed_time_grid = SFH_time_output[:-1]
    elif bin_position == 'middle':
        mass_formed_time_grid = (SFH_time_output[1:] + SFH_time_output[:-1])/2
    elif bin_position == 'back':
        mass_formed_time_grid = SFH_time_output[1:]
    else:
        raise ValueError("Incorrect bin_position value ['front', 'middle', 'back']")

    return SFH_time_output, SFH_sfr_output, t_universe - mass_formed_time_grid, mass_formed


def interpolate_SSP_SED_linear(T_age, ages, seds):
    """
    Linearly interpolate an SSP SED at age T_age.

    Parameters
    ----------
    T_age : float
        Target SSP age (same units as `ages`)
    ages : ndarray, shape (n_spectra,)
        Sorted SSP age grid
    seds : ndarray, shape (n_spectra, n_wavelength)
        SSP SED template array

    Returns
    -------
    ssp_T : ndarray, shape (n_wavelength,)
        Interpolated SSP spectrum at T_age
    """
    
    
    if T_age <= ages[0]:
        return seds[0].copy()
    elif T_age >= ages[-1]:
        return seds[-1].copy()
    else:
    
        i = np.searchsorted(ages, T_age) - 1
        j = i + 1

        age_i = ages[i]
        age_j = ages[j]

        # linear interpolation weight
        w = (T_age - age_i) / (age_j - age_i)

        # interpolate full spectrum
        ssp_T = (1.0 - w) * seds[i] + w*seds[j]


        return ssp_T



# Calculation
# SFR
SFH_input = np.genfromtxt("./../data/GALFORM_SFHs/gal_1_SFH_test.csv")
SFH_order = np.argsort(SFH_input[0])
SFH_time, SFH_sfr = SFH_input[0][SFH_order], 10**SFH_input[1][SFH_order]

Universe_time_grid = Universe_grid_generator()
# print (Universe_time_grid)

SFH_sfr_GALFORM = SFR_interp(Universe_time_grid, SFH_time, SFH_sfr)

# Mass formed
Galaxy_lb_time = 5.0 #Gyr
Universe_age_at_gal_lb = z_to_age(0) - Galaxy_lb_time


gal_mass_formed_f = mass_formed(SFH_time, SFH_sfr, Universe_age_at_gal_lb, bin_position = 'front')
gal_mass_formed_m = mass_formed(SFH_time, SFH_sfr, Universe_age_at_gal_lb, bin_position = 'middle')
gal_mass_formed_b = mass_formed(SFH_time, SFH_sfr, Universe_age_at_gal_lb, bin_position = 'back')

# The wavelength information
wave = np.genfromtxt('./../data/FSPS_wave.csv')
# The SED templates
SED_templates_file_path = './../data/FSPS_SED_templates.npy'
FSPS_SED_templates = np.load(SED_templates_file_path, mmap_mode = 'r')

SSP_age_grid_file_path = './../data/FSPS_param_grid_LHS.csv'
SSP_age_grid = np.genfromtxt(SSP_age_grid_file_path)


# Calculate the CSP SED from SFH and mass formed
gal_SED_f = np.zeros_like(wave)
for i, T_SSP in enumerate(gal_mass_formed_f[2]):
    # print (T_SSP, gal_mass_formed_f[3][i]/1e9)
    gal_SED_f += gal_mass_formed_f[3][i]*interpolate_SSP_SED_linear(T_SSP, SSP_age_grid, FSPS_SED_templates)

gal_SED_m = np.zeros_like(wave)
for i, T_SSP in enumerate(gal_mass_formed_m[2]):
    # print (T_SSP, gal_mass_formed_m[3][i]/1e9)
    gal_SED_m += gal_mass_formed_m[3][i]*interpolate_SSP_SED_linear(T_SSP, SSP_age_grid, FSPS_SED_templates)

gal_SED_b = np.zeros_like(wave)
for i, T_SSP in enumerate(gal_mass_formed_b[2]):
    # print (T_SSP, gal_mass_formed_b[3][i]/1e9)
    gal_SED_b += gal_mass_formed_b[3][i]*interpolate_SSP_SED_linear(T_SSP, SSP_age_grid, FSPS_SED_templates)

print (Universe_age_at_gal_lb)



plt.figure(figsize = (12, 4))
plt.plot(wave, gal_SED_f*wave, lw = 3, label = 'front')
plt.plot(wave, gal_SED_m*wave, lw = 2, label = 'middle')
plt.plot(wave, gal_SED_b*wave, lw = 1, label = 'back')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 'upper right', fontsize=12)
plt.ylim(1e0, 1e13)
plt.show()











