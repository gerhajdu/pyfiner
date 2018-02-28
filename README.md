# PyFiNeR: Fitting Near-infrared RR Lyrae light curves

 This routine implements the RR Lyrae near-infrared light curve fitting techniques described in Hajdu et al. (2018)

 The K-band light curve are fit with the Fourier transformed principal components coefficients; the H-band data is fit with the same light-curve shape; the shape of the J-band light curve is approximated with the K-band principal component amplitudes

 This routine was developed for:
 - Python 2.7
 - Numpy 1.12
 - Scipy 1.0
 - Matplotlib 2.1

## Installation

 Copy the files pyfiner, pyfiner.py, lcutil_pyfiner.py and coefs_pca.npz to the same directory in the system PATH

## Usage

 The routine requires six arguments to run; in order, these are:
 - COMMENT: a string that is passed printed out in the end; usually the name of the variable; can use quotation marks to carry additional information
 - PERIOD: the initial period estimate of the variable
 - FILE_K_LC: path of the K light curve file; must exist
 - FILE_J_LC: path of the J light curve file; must exist, but can be an empty file
 - FILE_H_LC: path of the H light curve file; must exist, but can be an empty file
 - OUTPUT.pdf: name of the output pdf file
