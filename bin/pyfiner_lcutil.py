# Light curve utilities

from __future__ import print_function
import numpy as np
from os import path

# We read in the necessary constants contained in the coefs_pca.npz file
# This file contains the harmonic coefficients of the Fourier transform
# of the K light curve PCA eignevectors, the corresponding harmonic
# coefficients of the J-band fits, as well as the upper and lower limits
# of the principal component amplitudes, in the format of U1, U2/U1, U3/U1
# and U4/U1

try:
    npzfile = np.load(path.dirname(__file__)+'/pyfiner_coefs.npz')
    harmonic_coef_k = npzfile['k']
    harmonic_coef_j = npzfile['j']
    maxs            = npzfile['maxs']
    mins            = npzfile['mins']
except:
    print("Coefficient file not found")
    exit()

# The file we have read in has the information up to the 12 Fourier order
harmonic_order=12

# Simple function to calculate the phases of the timings of a light curve

def return_LC_phases(period, LC, epoch=0.0):
    return np.modf( (LC-epoch)/period )[0]


# Function to give back the sum of a Fourier series at given positions
# Can be used to give back the light-curve shape, if called with period=1.
# and positions (epochs) between 0 and 1

def return_harmonic_LC(order, period, coefs, intercept, positions):
    result=np.zeros(positions.size)
    for i in range(1,order+1):
        result=result +\
               coefs[i*2-2] *np.sin( i*2*np.pi*positions/period ) +\
               coefs[i*2-1] *np.cos( i*2*np.pi*positions/period)
    return result+intercept


# This function returns the residuals for a given light curve with
# amplitudes U1..U4
# Calls return_harmonic_LC, but passes the collapsed U*harmonic_coef
# to it
# This function is passed to the least_squares function of scipy.optimize
# as its first argument, which optimizes U1..U4 for the selected loss
# function (in our case, the Huber loss function)

def return_residuals(x, t, y):
    phases = return_LC_phases(x[2],t,x[1])
    U = np.array((x[3], x[3]*x[4], x[3]*x[5], x[3]*x[6]))
    return y - return_harmonic_LC(harmonic_order,x[2],(U*harmonic_coef_k.T).T.sum(axis=0),x[0],phases*x[2])


# Calculates the intensity averaged magnitude of a light curve
# The phases of the light curve have to be sampled evenly, but
# in a [0..1) way; should be also dense enough

def calc_int_average(lc):
    return -2.5*np.log10( np.average( (10.**(-lc/2.5)) ) )

