# PyFiNeR: Fitting Near-infrared RR Lyrae light curves

This routine implements the RR Lyrae near-infrared light curve fitting techniques described by
[Hajdu et al. (2018)](https://arxiv.org/abs/1804.01456).

The K-band light curve are fit with the Fourier transformed principal components coefficients;
the H-band data is fit with the same light-curve shape; the shape of the J-band light curve is
approximated with the K-band principal component amplitudes.

This routine was developed for:
 - `Python` 2.7+ or 3.6+
 - `Numpy` 1.12+
 - `Scipy` 1.0+
 - `Matplotlib` 2.1+

## Installation

Copy all files from the `bin` directory to the same directory in the system PATH.
If you get "ImportError: No module named builtins" error while using Python 2.7,
install the `future` package.

## Usage

The routine requires six arguments to run; in order, these are:
- COMMENT: a string that is passed printed out in the end; usually the name of the variable;
can use quotation marks to carry additional information
- PERIOD: the initial period estimate of the variable
- FILE_K_LC: path of the K light curve file; must exist
- FILE_J_LC: path of the J light curve file; must exist, but can be an empty file
- FILE_H_LC: path of the H light curve file; must exist, but can be an empty file
- OUTPUT.pdf: name of the output pdf file

Given these arguments, the routine does the fit and gives back the following output:
- COMMENT: the string from the argument list
- PERIOD: the final period 
- K_MAG: magnitude averaged mean K magnitude
- K_INT: intensity averaged mean K magnitude
- J_MAG_MEAN : magnitude averaged mean J magnitude, calculated with the mean offset
- J_INT_MEAN : intensity averaged mean J magnitude, calculated with the mean offset
- J_MAG_MEDI : magnitude averaged mean J magnitude, calculated with the median offset
- J_INT_MEDI : intensity averaged mean J magnitude, calculated with the median offset
- H_MAG_MEAN : magnitude averaged mean H magnitude, calculated with the mean offset
- H_INT_MEAN : intensity averaged mean H magnitude, calculated with the mean offset
- H_MAG_MEDI : magnitude averaged mean H magnitude, calculated with the median offset
- H_INT_MEDI : intensity averaged mean H magnitude, calculated with the median offset
- U1: amplitude of the first principal component of the K light curve
- U2: amplitude of the second principal component of the K light curve
- U3: amplitude of the third principal component of the K light curve
- U4: amplitude of the fourth principal component of the K light curve
- RMSE: the Root Mean Square Error of the K-band fit
- COST: the final value of the Huber cost function
- COSTN: the Uber cost function per K-band light curve cost
- N_FINAL: the number of K-band light curve points in the final fit
- N_INITIAL: the initial number of K-band light curve points

