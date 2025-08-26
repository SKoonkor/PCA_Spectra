# PCA_Spectra
This repository is for essential calculations for reconstruction of synthetic galaxy spectra using the Principal Component Analysis (PCA). Mainly, it makes use for the python-FSPS code for generating the simple stellar population (SSP) spectra
## The main structure of the code:
1. Use the Latin Hypercube Sampling method to generate the FSPS parameter grid, required for the SSP SED set.
2. Mask the SED to a specific redshift range essential for a given project (observation configurations)
3. Data preprocessing via L1, L2, or 5500A normalisation.
4. PCA application on the normalised SEDs.
5. SED reconstruction using the PC coefficients.



The output required for calculating galaxy SEDs:
1. Mean Spectra (the output from step 3 above).
2. First N PCs (Obtains from step 4).
3. PC coefficient of each PC.
