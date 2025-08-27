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



