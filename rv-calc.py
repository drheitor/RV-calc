#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 09:51:08 2024

@author: heitor
"""

#turbospectrum
from __future__ import annotations
try:
    from scripts_for_plotting import *
except ModuleNotFoundError:
    import sys
    sys.path.append('/Users/heitor/Desktop/NLTE-code/TSFitPy/')
    from scripts_for_plotting import *
    
    
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import sys
import pandas as pd

from scipy.interpolate import interp1d
from scipy.optimize import minimize

from scipy.ndimage import gaussian_filter1d

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.optimize import minimize, differential_evolution


sns.set_style("white")
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 1.5})



turbospectrum_paths = {"turbospec_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/turbospectrum/exec-gf/",  # change to /exec-gf/ if gnu compiler
                       "interpol_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/scripts/model_interpolators/",
                       "model_atom_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/input_files/nlte_data/model_atoms/",
                       "departure_file_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/input_files/nlte_data/",
                       "model_atmosphere_grid_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/input_files/model_atmospheres/",
                       "line_list_path": "/Users/heitor/Desktop/NLTE-code/TSFitPy/input_files/linelists/linelist_for_fitting/"}


#---------------


def doppler_shift(wavelength, velocity):
    c = 299792.458  # speed of light in km/s
    return wavelength * np.sqrt((1 + velocity / c) / (1 - velocity / c))

def chi_squared(params, observed_wavelength, observed_flux, template_wavelength, template_flux):
    velocity = params[0]
    shifted_wavelength = doppler_shift(template_wavelength, velocity)
    interp_func = interp1d(shifted_wavelength, template_flux, kind='linear', fill_value='extrapolate')
    shifted_template_flux = interp_func(observed_wavelength)
    
    chi2 = np.sum((observed_flux - shifted_template_flux) ** 2)
    return chi2

def fit_spectrum(observed_wavelength, observed_flux, template_wavelength, template_flux):
    bounds = [(-1000, 1000)]  # search for velocity in the range -1000 km/s to 1000 km/s
    result = differential_evolution(chi_squared, bounds, args=(observed_wavelength, observed_flux, template_wavelength, template_flux), strategy='best1bin', disp=True)
    best_velocity = result.x[0]
    return best_velocity

def crop_spectrum(wavelength, flux, min_wavelength, max_wavelength):
    """
    Crop the spectrum to a given wavelength range.
    
    Parameters:
    - wavelength: array-like, the wavelengths of the spectrum
    - flux: array-like, the corresponding flux values
    - min_wavelength: float, the lower bound of the wavelength range
    - max_wavelength: float, the upper bound of the wavelength range
    
    Returns:
    - cropped_wavelength: array-like, the wavelengths within the specified range
    - cropped_flux: array-like, the flux values corresponding to the cropped wavelengths
    """
    # Convert to numpy arrays if not already
    wavelength = np.array(wavelength)
    flux = np.array(flux)
    
    # Create a mask for the specified wavelength range
    mask = (wavelength >= min_wavelength) & (wavelength <= max_wavelength)
    
    # Apply the mask to crop the spectrum
    cropped_wavelength = wavelength[mask]
    cropped_flux = flux[mask]
    
    return cropped_wavelength, cropped_flux



def turbo(teff,logg,met,vmic,lmin,lmax,FWHM):
    
    #teff = 5500
    #logg = 4.0
    #met = -1.0
    #vmic = 1.0
    
    #lmin = 4600
    #lmax = 5500

    ldelta = 0.01
    
    atmosphere_type = "1D"   # "1D" or "3D"
    nlte_flag = False
    
    elements_in_nlte = ["Fe", "Mg"]  # can choose several elements, used ONLY if nlte_flag = True
    element_abundances = {"Mg": 3.0, "O": 0.0}  # elemental abundances [X/Fe]; if not written solar scaled ones are used
    include_molecules = False  # way faster without them
    
    # plots the data, but can also save it for later use
    wavelength, flux = plot_synthetic_data(turbospectrum_paths, teff, logg, met, vmic, lmin, lmax, ldelta, atmosphere_type, nlte_flag, elements_in_nlte, element_abundances, include_molecules, resolution=0, macro=0, rotation=0, verbose=False)

    #convolution

    #FHWM (0.12A = 12 pixels) = 2.354 * sigma = 2.354 * 5.09 pixels
    #FWHM= 0.12
    
    pix = FWHM/ldelta

    sig= pix/2.354

    z = gaussian_filter1d(flux, sigma=sig)
    
    return wavelength, z




def add_star_data(filename, star_name, rv):
    # Open the file in append mode
    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        
        # Write the new star_name and rv
        writer.writerow([star_name, rv])

    # File is automatically closed when the 'with' block is exited



#/Users/heitor/Desktop/Astronomy/Lund/ASA-spec/codes

#---------------


#spec_path = '../out/helio/HD-219617_564l_OB2_helio.csv'  # Path to your FITS file

path='../out/helio/'

file_name = sys.argv[1]


spec_path = path + file_name


star_name = file_name


print(f"Measuring rv for: {star_name} spectrum")

spec = pd.read_csv(spec_path)

debug = 0


#-----------

print('Observed read')


# Example synthetic data
observed_wavelength = spec['wave']

#no noise 
observed_flux = spec['flux']


if debug == 1:
    
    plt.plot(observed_wavelength,observed_flux)
    
 

#limits 

#390
if max(observed_wavelength) < 4600 and min(observed_wavelength) > 3100:
    lmin = 4300
    lmax = 4450
   
    
#564l
if max(observed_wavelength) < 6800 and min(observed_wavelength) > 4500 :
    lmin = 4800
    lmax = 5200
    
    
#564u
if max(observed_wavelength) < 6800 and min(observed_wavelength) > 4800 :
    lmin = 5840
    lmax = 6000



observed_wavelength, observed_flux = crop_spectrum(observed_wavelength, observed_flux, lmin, lmax)


#if you wish to edith it yourself
# Generating synthetic data
#true_velocity = 51.724 # km/s

#noise 
#observed_flux = spec['flux'] + np.random.normal(0, max(spec['flux'])/1e1, len(observed_wavelength))  # observed flux
#flux_obs = np.sin(wavelength_obs/100) + np.random.normal(0, 0.05, len(wavelength_obs))  # observed flux


#-----------
#model


#model = pd.read_csv('data.csv')


#parameters rought numbers
teff = 5500
logg = 4.0
met = -1.0

vmic = 1.0


FWHM= 0.12

 
#-----------

print('Template')


#make template

template_wavelength, template_flux = turbo(teff,logg,met,vmic,lmin,lmax,FWHM)

template_flux = template_flux * np.median(observed_flux)


#-----------
#fit
print('fitting...')



#template_wavelength = doppler_shift(spec['wave'], true_velocity)  # template spectrum wavelengths
#template_flux  = spec['flux']  # template flux (without noise)

# Fit the spectrum
best_velocity = fit_spectrum(observed_wavelength, observed_flux, template_wavelength, template_flux)
print('\n')
print(f"Best fit radial velocity: {best_velocity} km/s")
print('\n')




best_velocity_r = -round(best_velocity,3)



#-----------
print('plotting...')


#plotting

fig, ax = plt.subplots(figsize=(15, 5))


# Plotting the results
#shifted_template_wavelength = doppler_shift(template_wavelength, best_velocity)
#shifted_template_flux = interp1d(shifted_template_wavelength, template_flux, kind='linear', fill_value='extrapolate')(observed_wavelength)


shifted_observed_wavelength = doppler_shift(observed_wavelength, -best_velocity)
shifted_observed_flux = interp1d(shifted_observed_wavelength, observed_flux, kind='linear', fill_value='extrapolate')(template_wavelength)


ax.plot(template_wavelength, template_flux, label='Template Spectrum', linestyle='dashed', linewidth=2.0, color='k')
ax.plot(observed_wavelength, observed_flux, label='Observed Spectrum', linestyle='dotted', color='k',alpha=0.3)


ax.plot(template_wavelength[2:-2], shifted_observed_flux[2:-2], label='Shifted Observed (Best Fit)', color='red')


ax.legend()
ax.set_xlabel('Wavelength (Angstroms)')
ax.set_ylabel('Intensity')
ax.set_title(f"radial velocity: {best_velocity_r} km/s")


plt.savefig('../fig/rv-figs/'+star_name+'_-rv.pdf')


#-----------
print('editing table...')



filename = 'stars_data.csv'

new_star_name = star_name
new_rv = best_velocity_r

add_star_data(filename, new_star_name, new_rv)





print('DONE')


#-----------






















#