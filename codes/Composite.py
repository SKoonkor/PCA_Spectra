import numpy as np
import pandas as pd
import os, sys

import matplotlib.pyplot as plt



#########################################################################################################
def SFH_exp_declining(t, T0, tau):
    '''
    This function returns the exponential declining SFH based on the time table provided.
    
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

def SFH_delayed_exp_declining(t, T0, tau):
    '''
    The delayed exponential declining star formation rate

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
    
    SFR = (t-T0)*SFH_exp_declining(t, T0, tau)

    return SFR

def SFH_lognormal(t, T0, tau):
    '''
    '''

    return np.exp(-(np.log(t) - T0)**2/(2*tau**2))/t

def SFH_DPL(t, tau, alpha, beta):

    return ((t/tau)**alpha + (t/tau)**(-beta))**-1



############################################
# Test zone



t_grid = np.logspace(-5, 1, 100)
T0 = 1e-4
tau = 0.3

alpha = 0.5
beta = 0.3


SFR_exp_declining = SFH_exp_declining(t_grid, T0, tau)
SFR_delayed_exp_declining = SFH_delayed_exp_declining(t_grid, T0, tau)
SFR_lognormal = SFH_lognormal(t_grid, T0, tau)
SFR_DPL = SFH_DPL(t_grid, tau, alpha, beta)

plt.figure(figsize = (12, 8))
plt.plot(t_grid, SFR_exp_declining, label = 'Exp. decl.')
plt.plot(t_grid, SFR_delayed_exp_declining, label = 'Delay exp. decl.')
plt.plot(t_grid, SFR_lognormal, label = 'Lognormal')
plt.plot(t_grid, SFR_DPL, label = 'DPL')
plt.xscale('log')
plt.legend(loc = 'upper left')
plt.show()
