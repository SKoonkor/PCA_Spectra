import os, sys
import numpy as np
import matplotlib.pyplot as plt


def plot_LatinHypercube(param_names, param_values, param_spacings, 
                        index_x = 0,
                        index_y = 1,
                        figsize = (12, 8),):

    fig, ax = plt.subplots(figsize = figsize)

    ax.scatter(param_values[:, index_x], param_values[:, index_y])
    
    ax.set_xlabel(param_names[index_x])
    ax.set_ylabel(param_names[index_y])
    ax.set_xscale(param_spacings[index_x])
    ax.set_yscale(param_spacings[index_y])

    return fig, ax

def plot_SED(ax, wave, SED,
             ls = 'solid',
             lw = 1,
             label = None,
             xlabel = None,
             ylabel = None,
             xscale = 'log',
             yscale = 'log',
             color = None,):

    # Check if the color for the plot is defined, if not pick a random RGB color
    if color == None:
        color = np.random.rand(3,)

    # Plot the SED
    ax.plot(wave, SED, color = color, label = label, ls = ls, lw = lw)
    
    # Plot configurations
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
