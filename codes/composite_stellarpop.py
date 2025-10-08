import numpy as np
import pandas as pd
import os, sys

import matplotlib.pyplot as plt



######################################################################################################

def t_ssp_log(t_SFH_grid):
    '''
    This function provides the time steps for the SSP age grid.
    The output is in log scale, to focus on the young stellar population that just recently form
    '''
    return max(t_SFH_grid) - min(t_SFH_grid) - np.append(0, np.logspace(-2, 
                                                      np.log10(max(t_SFH_grid) - min(t_SFH_grid)), 
                                                      30, 
                                                      endpoint = True))
def t_ssp_linear(t_SFH_grid):
    '''
    This function provides the time steps for the SSP age grid.
    The output is in linear scale, to focus on the average SFR at different time steps.
    '''
    return max(t_SFH_grid) - min(t_SFH_grid) - np.linspace(0, 
                                                           max(t_SFH_grid) - min(t_SFH_grid), 
                                                           30, 
                                                           endpoint = True)


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

def mass_formed_grid(sfr_ts, sfr_vals, age_edges, *, verbose = False):
    '''
    Compute stellar mass formed between each age bin for a CSP.

    Parameters
    ----------
    sfr_ts : array_like
        Cosmic times (increasing) where the SFR is known.
    sfr_vals : array_like
        SFR values corresponding to sfr_ts.
    age_edges : array_like
        Age grid edges [t_age_gal, ..., 0], *descending* order.

    Returns
    -------
    mass_formed : ndarray
        Mass formed in each age bin (same length as len(age_edge) - 1).
    '''
    sfr_ts = np.asarray(sfr_ts)
    sfr_vals = np.asarray(sfr_vals)
    age_edges = np.asarray(age_edges)
    
    # ensure SFR grid sorted in ascending cosmic time
    order = np.argsort(sfr_ts)
    sfr_ts = sfr_ts[order]
    sfr_vals = sfr_vals[order]



    # Sanity check
    if np.any(np.diff(sfr_ts) < 0):
        raise ValueError('sfr_ts must be strictly monotonic.')


    # Get t_sfr_start and t_sfr_finish (the earliest and latest time in the SFH)
    t_sfr_start, t_sfr_finish = sfr_ts[0], sfr_ts[-1]
    
    
    # if age_edges is descending, make an ascending copy for iteration
    descending = np.all(np.diff(age_edges) < 0)

    if descending:
        age_edges_sorted = age_edges[::-1] # youngest -> oldest
    else:
        age_edges_sorted = age_edges

    mass_bins_temp = np.zeros(len(age_edges_sorted) - 1)

    for i in range(len(age_edges_sorted) - 1):
        t_a = age_edges_sorted[i]
        t_b = age_edges_sorted[i + 1]

        # Skip if completely outside SFH range
        if (t_b < t_sfr_start) or (t_a > t_sfr_finish):
            continue

        # Clip bin edges to SFH limits
        t_a_clip = np.clip(t_a, t_sfr_start, t_sfr_finish)
        t_b_clip = np.clip(t_b, t_sfr_start, t_sfr_finish)

        # Time grid inside bin
        mask = (sfr_ts >= t_a_clip) & (sfr_ts <= t_b_clip)
        t_local = np.concatenate(([t_a_clip], sfr_ts[mask], [t_b_clip]))
        sfr_local = np.interp(t_local, sfr_ts, sfr_vals)

        # Integrate SFR(t) dt for this bin
        mass_bins_temp[i] = np.trapz(sfr_local, t_local)

        if verbose:
            print (f"Bin {i}: {t_a_clip:.3f} - {t_b_clip:.3f}, mass = {mass_bins_temp[i]:.4e}")

    mass_bins = mass_bins_temp[::-1] if descending else mass_bins_temp
    return mass_bins



############################################



