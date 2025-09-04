import numpy as np
from scipy.stats import qmc


def LHS_generator(param_names, param_ranges, param_spacings,
                  n_samples = 200):
    # This function generates a parameter grid using the Latin Hypercube Sampling technique
    # It take 3 lists containing
    #     1. Names of FSPS parameters
    #     2. Parameter ranges
    #     3. 'log' or 'linear' spacing
    
    # Convert param_ranges to Numpy array
    param_ranges = np.array(param_ranges) 
    
    # Create a Latin Hypercube sampler
    # This makes sure that the dimension of the LatinHypercube matches the dimension of parameter input
    sampler = qmc.LatinHypercube(d = len(param_names))
    
    # Generate samples in the unit cube [0, 1], that will be used for scaling later
    sample = sampler.random(n = n_samples)

    # Scale the unit cube to match the param_range defined in the input
    scaled_sample = qmc.scale(sample, param_ranges[:, 0], param_ranges[:, 1])


    # Apply log scaling if requested
    for i, spacing in enumerate(param_spacings):
        if spacing == 'log':
            scaled_sample[:, i] = 10**scaled_sample[:, i] # log10 scaling
    

    return scaled_sample

def normal_grid_generator(param_names, param_grid_dim, param_ranges, param_spacings):
    '''
    Parameters
    __________
    param_names : list (N_p, )
        List containing FSPS parameter names
    param_grid_dim: array (N_p, )
        Array containing the dimention of the grid
    param_ranges : array (N_p, 2)
        List containing the lower and upper limit of the ranges
    param_spacings : list (N_p, )
        List containing how each parameter is spaced ('linear' or 'log')

    Returns
    _______
    samples : array (N_s, N_p)
        Array containing the FSPS parameters for N_s SEDs
    '''

    param_ranges = np.asarray(param_ranges)

    param_value_set = []

    for i in range(len(param_names)):
        if param_spacings[i] == 'log':
            values = np.linspace(param_ranges[i][0], param_ranges[i][1], param_grid_dim[i])
            param_value_set.append(10**values)
        else:
            param_value_set.append(np.linspace(param_ranges[i][0], param_ranges[i][1], param_grid_dim[i]))
    param_value_set = np.array(param_value_set)

    tage_set = []
    logzsol_set = []
    for age in param_value_set[0]:
        for logzsol in param_value_set[1]:
            tage_set.append(age)
            logzsol_set.append(logzsol)
    
    param_values = np.array([np.array(tage_set), np.array(logzsol_set)]).transpose()
 
    return param_values



