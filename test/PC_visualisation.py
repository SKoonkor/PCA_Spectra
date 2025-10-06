import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys, os

sys.path.append('/Users/suteepornz/Documents/Suttikoon/Research_Projects/PCA_Spectra/codes')

from visualisations import plot_SED
from PCA_cal import denormalise_SED


###########################################################################################
# Define filenames containing PCs and Coeffs
PCA_outputs = '../outputs/'
PC_filename = 'PCA_PCs.csv'
coeff_filename = 'PCA_coeffs.csv'


# Define the data directory constainging SED scaling and mean
data_dir = '../data/'
param_name_filename = 'FSPS_param_grid_names.txt'
wave_filename = 'FSPS_wave_norm_std.csv'
SED_mean_filename = 'FSPS_SED_norm_mean_std.csv'
SED_factor_filename = 'FSPS_SED_norm_factor_std.csv'
SED_param_filename = 'FSPS_param_grid_LHS.csv'

###########################################################################################

param_names = np.genfromtxt(data_dir + param_name_filename)
wave = np.genfromtxt(data_dir + wave_filename)
PCs = np.genfromtxt(PCA_outputs + PC_filename)
coeffs = np.genfromtxt(PCA_outputs + coeff_filename)

params = np.genfromtxt(data_dir + SED_param_filename)
print (param_names)

if len(param_names) == 2:

    tage, logzsol = params[:, 0], params[:, 1]


    print (wave.shape)
    print (params[0].shape)
    print (coeffs[:, 0].shape)
    print (tage.shape)

    colors = ['red', 'orange', 'green', 'blue', 'purple']
    cmaps = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples']

    n_PC = 5
    fig, ax = plt.subplots(n_PC, 1, figsize = (12, 1*n_PC), gridspec_kw = {'hspace': 0.01})

    for i in range(n_PC):
        plot_SED(ax[i], wave, PCs[i], label = 'PC {}'.format(i + 1), color = colors[i])
        ax[i].set_yscale('linear')
        ax[i].set_ylabel('PC {}'.format(i))
    
        if i != n_PC - 1:
            ax[i].set_xticks([])

    ax[n_PC-1].set_xlabel('Wavelength [AA]')
    plt.savefig('PCs.png')
    # plt.show()


    fig2, ax2 = plt.subplots(n_PC, 1, figsize = (12, 24))
    # PC1
    for j in range(n_PC):
        ax2[j].scatter(tage, coeffs[:, j], s = 3, label = 'PC {}'.format(j+1), c = logzsol, cmap = cmaps[j])
        ax2[j].set_xscale('log')
        ax2[j].legend(loc = 'upper right', fontsize = 20)
        ax2[j].set_ylabel('PC coefficients', fontsize = 20)
        if j != n_PC - 1:
            ax2[j].set_xticks([])
    ax2[n_PC-1].set_xlabel('tage', fontsize = 20)
    plt.savefig('PCvsTage.png')
    plt.show()
else:
    tage = params


    print (wave.shape)
    print (params.shape)
    print (coeffs.shape)
    print (tage.shape)

    colors = ['red', 'orange', 'green', 'blue', 'purple']
    cmaps = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples']

    n_PC = 5
    fig, ax = plt.subplots(n_PC, 1, figsize = (12, 1*n_PC), gridspec_kw = {'hspace': 0.01})

    for i in range(n_PC):
        plot_SED(ax[i], wave, PCs[i], label = 'PC {}'.format(i + 1), color = colors[i])
        ax[i].set_yscale('linear')
        ax[i].set_ylabel('PC {}'.format(i))
    
        if i != n_PC - 1:
            ax[i].set_xticks([])

    ax[n_PC-1].set_xlabel('Wavelength [AA]')
    plt.savefig('PCs.png')
    # plt.show()


    fig2, ax2 = plt.subplots(n_PC, 1, figsize = (12, 24))
    # PC1
    for j in range(n_PC):
        ax2[j].scatter(tage, coeffs[:, j], s = 3, label = 'PC {}'.format(j+1), c = logzsol, cmap = cmaps[j])
        ax2[j].set_xscale('log')
        ax2[j].legend(loc = 'upper right', fontsize = 20)
        ax2[j].set_ylabel('PC coefficients', fontsize = 20)
        if j != n_PC - 1:
            ax2[j].set_xticks([])
    ax2[n_PC-1].set_xlabel('tage', fontsize = 20)
    plt.savefig('PCvsTage.png')
    plt.show()

