import numpy as np
import pandas as pd
import os, sys

from sklearn import preprocessing
from sklearn.decomposition import PCA as skPCA
from sklearn.decomposition import IncrementalPCA as skIncrementalPCA


def mask_SED(wave_filename, SED_filename, wave_min, wave_max, output_filename, chunk_size = 1000):
    '''
    Apply wavelength mask on SEDs without loading the entire file into memory.
    Read SEDs chunk by chunk and writes maksed results directly to disk.
    
    Parameters
    __________
    wave : array (N_lambda,)
        Wavelength grid
    SED_filename : str 
        The name of the file containing the input SED template file
    wave_min : float
        The lower limit of the wavelength 
    wave_max : float
        The upper limit of the wavelength
    output_filename : str
        The file for storing the masked results
    chunk_size : int
        The number of SEDs we want to mask each iteration


    Returns
    _______


    '''
    wave = np.genfromtxt(wave_filename)
    # Get wavelength mask

    wave_mask = np.where((wave >= wave_min) & (wave <= wave_max))[0]
    wave_output = wave[wave_mask]

    print ('\n*************************************************')
    print ('****************  Masking Process ***************')
    print ('*************************************************')

    print (f'\n Original wavelengths: shape {wave.shape}, '
           f'Min {wave.min():.2e} AA, Max {wave.max():.2e} AA')
    print ('f\n Ouput wavelengths: shape {wave_output.shape}, '
           f'Min {wave_output.min():.2e} AA, Max {wave_output.max():.2e} AA')
    print ('*************************************************')

    np.savetxt(wave_filename.replace('.csv', '_masked.csv'), wave_output)

    with open(SED_filename, 'r') as f_input, open(output_filename, 'w') as f_output:
        for line in f_input:
            row = np.fromstring(line, sep = ' ')
            row_masked = row[wave_mask]
            f_output.write(' '.join(f'{x:.10e}' for x in row_masked) + '\n')
        

    print (f'\n SEDs (masked) stored at:\n{output_filename}')
    print ('*************************************************')

def normalise_SED(wave, SED_filename, norm_type = 'std', eps = 0, chunk_size = 1000, output_filename = None):
    '''
    This perform a Memory-efficient SED normalisation.

    Parameters
    ----------
    wave : array (N_lambda,)
        wavelength grid
    SED_filename : str
        Instead of reading all SEDs at the sametime, we read the path to CSV file containing SEDs (N_s, N_lambda)
    norm_type : str
        'l1', 'l2', 'std', or '5500A'
    eps : float
        Small value to avoid division by zero
    chunk_size : int
        Number of rows to process at a time
    output_filename : str or None
        If given, normalised SEDs will be streamed to this file

    Returns
    _______

    SED_norm_factor : array
        Per-wavelength std (for 'std') or norms
    SED_norm_mean : array
        Per-wavelength mean
    '''
    print ('\n#################################################')
    print ('######### Memory-Efficient Norlisation ##########')
    print ('#################################################')

    print ('\n1st pass of the two-pass algorithm')
    print (' Computing mean and std across all spectra...')

    total = None
    total_sq = None
    n_rows = 0
    
    for chunk in np.array_split(np.loadtxt(SED_filename),
                               max(1, os.path.getsize(SED_filename)//(chunk_size*1000))):
        chunk = np.asarray(chunk, dtype = float)

        if total is None:
            total = np.zeros(chunk.shape[1])
            total_sq = np.zeros(chunk.shape[1])

        total += np.sum(chunk, axis = 0)
        total_sq += np.sum(chunk**2, axis = 0)
        n_rows += chunk.shape[0]

    mean_vals = total/n_rows
    std_vals = np.sqrt(total_sq/(n_rows-1) - mean_vals**2)
    
    print ('\n2nd pass of the two-pass algorithm')
    print (' Normalising spectra and writing output...')

    if output_filename:
        fout = open(output_filename, 'w')

    with open(SED_filename, 'r') as f:
        for line in f:
            row = np.fromstring(line, sep=' ')
            if norm_type == 'l1' or norm_type == 'l2':
                row_norm, _ = preprocessing.normalize(row,reshape(1, -1), norm = norm_type, return_norm = True)
                row_norm = row_norm.ravel()

            elif norm_type == '5500A':
                idx_5500A = np.argmin(abs(wave - 5500))
                factor = row[idx_5500A]
                row_norm = row/(factor + eps)
                    
            elif norm_type == 'std':
                row_norm = (row - mean_vals)/(std_vals + eps)

            else:
                raise Exception('Unknow norm_type (only l1, l2, 5500A, or std)')


            if output_filename:
                fout.write(' '.join(f'{val:.18e}' for val in row_norm) + '\n')

    if output_filename:
        fout.close()

    print ('\n#################################################')
    return std_vals if norm_type == 'std' else None, mean_vals if norm_type == 'std' else None, 






def incrementalPCA_SED(SED_filename, n_components = 20, output_dir = '../outputs/', chunk_size = 2000):
    '''
    Perform PCA (Incremental) on SED templates with reduced memory usage.

    Parameters
    ----------
    SED_filename : str
        Filename that contains input spectral energy distributions (2D array: [n_spectra, n_wavelengths])
    n_components : int
        Number of principal components
    output_dir : str
        Directory to save PCA results
    chunk_size: int
        Number of samples per batch for skIncrementalPCA

    Returns
    -------
    PCA_coeffs : np.ndarray
        PCA coefficients for each input spectrum
    PCA_PCs : np.ndarray
        Principal components (eigenvectors)
    ''' 

    print('\n#################################################')
    print('######## RUNNING Incremental PCA ON SPECTRA ######')
    print('#################################################')

    print(' Initialising IncrementalPCA')
    ipca = skIncrementalPCA(n_components = n_components, batch_size = chunk_size)

    print('\n First pass: fitting PCA in batches')
    reader = pd.read_csv(SED_filename, chunksize = chunk_size, header = None, sep = '\s+')
    for i, chunk in enumerate(reader):
        ipca.partial_fit(chunk.astype(float).values)
        print(f'  Fitted batch {i+1}')

    print('\n Second pass: transforming spectra')
    coeffs_list = []
    reader = pd.read_csv(SED_filename, chunksize = chunk_size, header = None, sep = '\s+')
    for i, chunk in enumerate(reader):
        coeffs = ipca.transform(chunk.astype(float).values)
        coeffs_list.append(coeffs)
        print(f'  Transformed batch {i+1}')

    PCA_coeffs = np.vstack(coeffs_list)
    PCA_PCs = ipca.components_

    os.makedirs(output_dir, exist_ok = True)
    np.savetxt(os.path.join(output_dir, 'PCA_coeffs.csv'), PCA_coeffs, fmt = '%.18e', delimiter = ' ')
    np.savetxt(os.path.join(output_dir, 'PCA_PCs.csv'), PCA_PCs, fmt = '%.18e', delimiter = ' ')

    return PCA_coeffs, PCA_PCs

def denormalise_SED(SED_norm, SED_std, SED_mean):
    # This de-normalise the normalised SEDs back to the original
    # This method can be applied at the very end as a postprocessing step
    # It is define (as simple as):
    #
    # SED_recon = SED_norm * SED_std + SED_mean,
    # where SED_norm is the normalised SED, SED_std is the normalisation factor for each SED
    #

    print ('\n#################################################')
    print ('############### Denormalise SED #################')
    print ('#################################################')
    return SED_norm*SED_std + SED_mean

def reconstruct_PCA_SED(coeffs, PCs, n_components):
    # This reconstruct the PCA SEDs back to the 'original SEDs'
    
    #######################################################################
    # input n_components cannot be larger than the total number of PCs or Featuers
    n_components = int(min(n_components, len(PCs), len(PCs[0])))
    print ('n_components used: ', n_components)
    SED_recon = np.zeros(len(PCs[0]))
    for i in range(n_components):
        SED_recon += coeffs[i]*PCs[i] 

    return SED_recon



###### TEST ######

def PCA_SED(SED, n_components = 20, output_dir = 'outputs/'):
    # The SED are dimentionally reduced to n_components
    # Instead of using the wavelengths as the vector for the spectra,
    # The spectra are projected onto the principal components (PCs)
    # The SED instead of represented as F_ij = f_ij * I_j,
    # They are described as S = c_ik* PC
    # notation: spectrum i-th, wavelentgh j-th.

    print ('\n#################################################')
    print ('########### RUNNING PCA ON SPECTRA ##############')
    print ('#################################################')


    # Initiailise the scikit-learn PCA package
    print (' Initialising sklearn PCA')
    PCA_fn = skPCA(n_components = n_components)

    # Fit the PCA on the SED templates
    print ('\n Fitting PCA on SEDs')
    PCA_fn.fit(SED)

    # Calculate the coefficient for each SED tempalte and get the eigen vectors (PCs)
    print ('\n Getting coeffs and PCs')
    PCA_coeffs = PCA_fn.transform(SED)
    PCA_PCs = PCA_fn.components_
    
    print ('\n Writing the output files')
    # Define the output filenames
    PCA_coeff_filename = 'PCA_coeffs.csv'
    PCA_PC_filename = 'PCA_PCs.csv'

    # Create the output directory
    os.makedirs(output_dir, exist_ok = True)


    np.savetxt(output_dir + PCA_coeff_filename, PCA_coeffs)
    print ('\n PCA coeffs stored at: ')
    print (output_dir + PCA_coeff_filename)
    
    np.savetxt(output_dir + PCA_PC_filename, PCA_PCs)
    print ('\n PCA PCs stored at: ')
    print (output_dir + PCA_PC_filename)

    print ('\n PCA RUN COMPLETE')
    print ('\n################################################')
    
    return PCA_coeffs, PCA_PCs
