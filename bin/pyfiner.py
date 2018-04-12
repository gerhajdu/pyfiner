#!/usr/bin/env python
# PyFiNeR: Fitting Near-infrared RR Lyrae light curves

# This routine implements the RR Lyrae near-infrared light curve fitting techniques
# described in Hajdu et al. (2018)

# The K-band light curve are fit with the Fourier transformed principal components
# coefficients; the H-band data is fit with the same light-curve shape; the shape
# of the J-band light curve is approximated with the K-band principal component
# amplitudes

# This routine was developed for Python 2.7+ or 3.6+, Numpy 1.12+, Scipy 1.0+ and
# Matplotlib 2.1+

# We import the necessary packages and make sure that matplotlib is not in interactive
# mode

from __future__ import print_function
from builtins import range
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from scipy.optimize import least_squares
from sys import argv, exit
from pyfiner_lcutil import *
plt.ioff()


# We try to read in the arguments and open the light curves in the different bands
# We also calculate a basic shift to the times of the data

try:
    NAME      = argv[1]
    period    = np.float64(argv[2])
    K_LC_path = argv[3]
    J_LC_path = argv[4]
    H_LC_path = argv[5]
    saveto    = argv[6]
except:
    print("Usage: python",argv[0],"COMMENT PERIOD FILE_K_LC FILE_J_LC FILE_H_LC OUTPUT.pdf")
    print("COMMENT: a string passed through to the output, usually the name of the star")
    print("FILE_K_LC: path of the K light curve file")
    print("FILE_J_LC: path of the J light curve file, can be an empty file")
    print("FILE_H_LC: path of the H light curve file, can be an empty file")
    print("OUTPUT.pdf: name of the file to plot the results of the fit")

    exit()

try:
    LC_K=np.loadtxt(K_LC_path, unpack=True)
    jd_min=[np.min(LC_K[0,0])]
except:
    print("{:s}: K light curve file not found!".format(NAME))
    exit()

try:
    LC_J=np.loadtxt(J_LC_path, unpack=True)
    if LC_J.size==0:
        no_J = True
    else:
        no_J = False
        jd_min.append(np.min(LC_J[0]))
        
except:
    print("{:s}: J light curve file not found!".format(NAME))
    exit()

try:
    LC_H=np.loadtxt(H_LC_path, unpack=True)
    if LC_H.size==0:
        no_H = True
    else:
        no_H = False
        jd_min.append(np.min(LC_H[0]))
except:
    print("{:s}: H light curve file not found!".format(NAME))
    exit()

shift = np.int(np.min(jd_min))-2
LC_K[0]  = LC_K[0] - shift


# We make a first cut in the light curves at +/- 0.5 mag around the median value,
# while also discarding every magnitude brighter than 9 and dimmer than 20

first_mask = (LC_K[1] > 9) * (LC_K[1]<20)
mask =   ( LC_K[1] < np.median(LC_K[1][first_mask]+0.5) )\
       * ( LC_K[1] > np.median(LC_K[1][first_mask]-0.5) )

LC_K2 = LC_K[:,mask]


# With VVV photometry, it is possible that we won't have enough points ot continue
# with the fitting after the first cut; in this case, we exit here

if (LC_K2.shape[1]<10):
    print(NAME, "9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9", LC_K2.shape[1], LC_K.shape[1])
    exit()


# We start by defining the trusted regions
# We also define a few variables to store the costs and fits of the first round of fits

lower = np.array((-np.inf, 0.0, period-0.001, mins[0], mins[1], mins[2], mins[3]))
upper = np.array(( np.inf, 2.0, period+0.001, maxs[0], maxs[1], maxs[2], maxs[3]))

first_fits = np.zeros((5,7))
costs      = np.zeros(5)


# We fit the models starting from different phases, then choose the best one to
# be the starting point in the next step of the fitting; we also do a 3.5 sigma cut,
# where sigma is estimated from the Median Absolute Deviation

for i in range(5):
    parameters_to_fit =  np.asarray( (np.median(LC_K2[1]), 1.0+(-2+i)*period/5., period, 0.7, 0., 0., 0.))
    huber             = least_squares(return_residuals, x0=parameters_to_fit, args=(LC_K2[0], LC_K2[1]),
                                      loss='huber', f_scale=0.05, bounds=(lower, upper))
    first_fits[i] = huber.x
    costs[i]      = huber.cost

best_id       = np.argmin(costs)
residuals     = return_residuals(first_fits[best_id], LC_K2[0], LC_K2[1])
residuals_std = np.median(np.abs(residuals))*1.4826
mask_35sigma  = (np.abs(residuals) < 3.5*residuals_std )
LC_K3         = LC_K2[:,mask_35sigma]


# We refit the model on the remaining light-curve points, but starting from the phase
# determined in the previous step

parameters_to_fit = np.asarray( (np.median(LC_K3[1]), first_fits[best_id, 1], first_fits[best_id, 2],
                               first_fits[best_id, 3], first_fits[best_id, 4], first_fits[best_id, 5], first_fits[best_id, 6]))
huber = least_squares(return_residuals, x0=parameters_to_fit, args=(LC_K3[0], LC_K3[1]),
                      loss='huber', f_scale=0.05, bounds=(lower, upper))


# The vector U contains the fitted amplitudes of the individual principal components
U = np.array((huber.x[3], huber.x[3]*huber.x[4], huber.x[3]*huber.x[5], huber.x[3]*huber.x[6]))

residuals_K3_final = return_residuals(huber.x, LC_K3[0], LC_K3[1])
RMSE= np.sqrt( np.mean(residuals_K3_final**2))


# We start plotting the results
#

LC_K2_phases = return_LC_phases(huber.x[2], LC_K2[0], 0.0)
LC_K3_phases = LC_K2_phases[mask_35sigma]

plt.style.use('seaborn-dark-palette')

fig, ax1 = plt.subplots(3, sharex=True)
fig.set_size_inches(9,10)
fig.subplots_adjust(hspace=0.05)
ax1[0].invert_yaxis()
ax1[1].invert_yaxis()
ax1[2].invert_yaxis()

plt.xlim(-0.01,2.01)

centers_long = np.linspace(-0.01,2.01,202)
fitted_LC_final  = return_harmonic_LC(harmonic_order, 1.0, (U*harmonic_coef_k.T).T.sum(axis=0),
                                      huber.x[0], centers_long-huber.x[1]/huber.x[2])

K_int_mean   = calc_int_average(return_harmonic_LC(harmonic_order, 1.0, (U*harmonic_coef_k.T).T.sum(axis=0),
                                0.0, np.linspace(0.0,0.995,200)-huber.x[1]/huber.x[2]) + huber.x[0])

ax1[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
ax1[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
ax1[2].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
ax1[0].set_ylabel('$\mathrm{K}_\mathrm{S}$')
ax1[1].set_ylabel('H')
ax1[2].set_ylabel('J')

ax1[0].plot(LC_K2_phases,   LC_K2[1],'r.')
ax1[0].plot(LC_K2_phases+1, LC_K2[1],'r.')
ax1[0].plot(LC_K3_phases,   LC_K3[1],'b.')
ax1[0].plot(LC_K3_phases+1, LC_K3[1],'b.')
ax1[0].plot(centers_long, fitted_LC_final, 'r-', lw=2)

# We only calculate the average magnitudes and plot the J and H light curves if they have data

if no_J == True:
    J_mag_mean   = 99.999
    J_mag_median = 99.999
    J_int_mean   = 99.999
    J_int_median = 99.999
else:
    LC_J[0] = LC_J[0] - shift
    LC_J_phases        = return_LC_phases(huber.x[2], LC_J[0], 0.0)
    J_predicted_shape  = return_harmonic_LC(harmonic_order, 1.0,        (U*harmonic_coef_j.T).T.sum(axis=0),
                                            0.0, centers_long-huber.x[1]/huber.x[2])
    
    J_predicted_at_obs = return_harmonic_LC(harmonic_order, huber.x[2], (U*harmonic_coef_j.T).T.sum(axis=0),
                                            0.0, LC_J[0]-huber.x[1])
    J_mag_mean   = np.mean(LC_J[1]-J_predicted_at_obs)
    J_mag_median = np.median(LC_J[1]-J_predicted_at_obs)
    J_int_mean   = calc_int_average(return_harmonic_LC(harmonic_order, 1.0, (U*harmonic_coef_j.T).T.sum(axis=0),
                                    0.0, np.linspace(0.0,0.995,200)-huber.x[1]/huber.x[2]) + J_mag_mean)
    J_int_median = J_int_mean + J_mag_median - J_mag_mean

    ax1[2].plot(centers_long, J_predicted_shape+J_mag_mean, 'k-', lw=2)
    ax1[2].plot(centers_long, J_predicted_shape+J_mag_median, 'k--', lw=2)
    ax1[2].plot(LC_J_phases,  LC_J[1],'b.', ms=15.0)
    ax1[2].plot(LC_J_phases+1,LC_J[1],'b.', ms=15.0)


if no_H == True:
    H_mag_mean   = 99.999
    H_mag_median = 99.999
    H_int_mean   = 99.999
    H_int_median = 99.999
else:
    LC_H[0] = LC_H[0] - shift
    LC_H_phases        = return_LC_phases(huber.x[2], LC_H[0], 0.0)
    
    H_predicted_at_obs = return_harmonic_LC(harmonic_order, huber.x[2], (U*harmonic_coef_k.T).T.sum(axis=0),
                                            0.0, LC_H[0]-huber.x[1])
    H_mag_mean   = np.mean(  LC_H[1]-H_predicted_at_obs)
    H_mag_median = np.median(LC_H[1]-H_predicted_at_obs)
    H_int_mean   = H_mag_mean   - huber.x[0] + K_int_mean
    H_int_median = H_mag_median - huber.x[0] + K_int_mean

    ax1[1].plot(centers_long, fitted_LC_final-huber.x[0]+H_mag_mean,   'k-',  lw=2)
    ax1[1].plot(centers_long, fitted_LC_final-huber.x[0]+H_mag_median, 'k--', lw=2)
    ax1[1].plot(LC_H_phases,  LC_H[1],'b.', ms=15.0)
    ax1[1].plot(LC_H_phases+1,LC_H[1],'b.', ms=15.0)


if LC_K3.shape[1] > 0:
    cost=huber.cost
else:
    cost=9.9999

textstr = '{:10s}          P: {:9.7f}  N: {:d} ({:d})\n'.format(NAME, huber.x[2], LC_K3.shape[1], LC_K.shape[1])
textstr+= 'RMSE: {:6.4f}   cost: {:6.4f}   cost/N: {:8.6f}\n'.format(RMSE, cost, cost/LC_K3.shape[1])
textstr+= 'K = {:6.3f}  J = {:6.3f} ({:6.3f})  H = {:6.3f} ({:6.3f}) mag\n'.format(huber.x[0], J_mag_mean, J_mag_median, H_mag_mean, H_mag_median)
textstr+= 'K = {:6.3f}  J = {:6.3f} ({:6.3f})  H = {:6.3f} ({:6.3f}) int'.format(K_int_mean, J_int_mean, J_int_median, H_int_mean, H_int_median)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1[0].text(.02, 1.04,textstr,fontsize=11,
            horizontalalignment='left',verticalalignment='bottom',
            bbox=props,
            transform=ax1[0].transAxes)


# The figure is saved to a pdf with the name given as the 

plt.savefig(saveto, bbox_inches='tight', format='pdf')


# Finally, we print out the results of the fit

textstr = '{:7s} {:9.7f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} '.format(NAME, huber.x[2], huber.x[0], K_int_mean, J_mag_mean, J_int_mean)
textstr+= '{:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} '.format(J_mag_median, J_int_median, H_mag_mean, H_int_mean, H_mag_median, H_int_median)
textstr+= '{:9.6f} {:9.6f} {:9.6f} {:9.6f} {:6.4f} {:6.4f} '.format(U[0], U[1], U[2], U[3], RMSE, huber.cost)
textstr+= '{:8.6f} {:4d} {:4d}'.format(cost/LC_K3.shape[1], LC_K3.shape[1], LC_K.shape[1])


print(textstr)

