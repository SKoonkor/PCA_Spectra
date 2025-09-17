import numpy as np
import pandas as pd
import os, sys

import matplotlib.pyplot as plt



#########################################################################################################
def SFH_tau(t, T0, tau):
    '''
    This function returns the exponential declining (tau) SFH based on the time table provided.
    
       SFR(t) = exp(-(t-T0)/tau) if t > T0; else 0.

    Pamaters
    ________
    t : array (N_t,)
        A table containing the age bin
    T0 : float 
        The time where the initial star formation occurs
    tau : float
        The e-folding time-scale for the exponential decline.

    Returns
    _______
    SFR : array (N_t,)
        The star formation rate at time t
    '''

    SFR = np.zeros_like(t)

    for i, t_i in enumerate(t):
        if t_i > T0:
            SFR[i] = np.exp(-(t_i - T0)/tau)
        else:
            SFR[i] = 0
    return SFR

def SFH_delayed_tau(t, T0, tau):
    '''
    The delayed exponential declining (delayed tau) SFH

    Parameters
    __________
    t : array (N_t,)
        A table containing the age bin
    T0 : float 
        The time where the initial star formation occurs
    tau : float
        The e-folding time-scale for the exponential decline.

    Returns
    _______
    SFR : array (N_t,)
        The star formation rate at every time step in t
    '''

    SFR = np.zeros_like(t)
    
    SFR = (t-T0)*SFH_tau(t, T0, tau)

    return SFR

def SFH_lognormal(t, T0, tau):
    '''
    '''

    return np.exp(-(np.log(t) - T0)**2/(2*tau**2))/t

def SFH_DPL(t, tau, alpha, beta):

    return ((t/tau)**alpha + (t/tau)**(-beta))**-1

def SFH_burst(t, tburst, sigma = 0.05):
    '''
    Gaussian starburst centred at T_burst.

    Parameters
    ----------
    t : array
        Age grid (Gyr).
    tburst : float
        Burst time (Gyr).
    sigma : float
        Width of the burst (Gyr).

    Returns
    -------
    SFR : array
        Burst contribution to SFR(t)
    '''
    return np.exp(-0.5*((t - tburst)/sigma)**2)

def build_composite_SFH(t, components):
    '''
    Build composite SFH from multiple components.

    Parameters
    ----------
    t : array 
        Age grid (Gyr).
    components : list of dict
        Each dict specifies a component. Example:
        {'type': 'exp', 'T0': 0.5, 'tau': 2.0, 'strength': 0.5}
        {'type': 'burst', 'tburst': 5.0, 'strength': 0.5, 'sigma': 0.1}
    Returns
    ------- 
    stf_norm : array
        Normalised SFR(t) such that int sfr(t)*dt = 1.
    '''
    sfr_total = np.zeros_like(t)

    for comp in components:
        if comp['type'] == 'tau':
            sfr = SFH_tau(t, comp['T0'], comp['tau'])
            sfr_total += comp.get('strength', 1.0) * sfr
        elif comp['type'] == 'delayed_tau':
            sfr = SFH_delayed_tau(t, comp['T0'], comp['tau'])
            sfr_total += comp.get('strength', 1.0) * sfr
        elif comp['type'] == 'lognormal':
            sfr = SFH_lognormal(t, comp['T0'], comp['tau'])
            sfr_total += comp.get('strength', 1.0) * sfr
        elif comp['type'] == 'DPL':
            sfr = SFH_DPL(t, comp['tau'], comp['alpha'], comp['beta'])
            sfr_total += comp.get('strength', 1.0) * sfr
        elif comp['type'] == 'burst':
            sfr = SFH_burst(t, comp['tburst'], comp.get('sigma', 0.05))
            sfr_total += comp.get('strength', 1.0) * sfr
        else:
            raise ValueError(f"Unknown component type: {comp['type']}")

    # Normalise to 1 Msun formed
    dt = np.gradient(t) # non-uniform spacing safe
    mass_formed = np.sum(sfr_total * dt)
    sfr_norm = sfr_total/mass_formed
    
    return sfr_norm

############################################



