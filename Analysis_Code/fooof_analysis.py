# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:05:02 2021

@author: ssshe
"""

# General imports
import numpy as np
# from scipy.io import loadmat, savemat
import scipy.io as sio
# from os.path import dirname, join as pjoin
import os
import matplotlib.pyplot as plt

# Import the FOOOF object, and custom helper & utility functions
from fooof import FOOOF, FOOOFGroup
from fooof.analysis import get_band_peak_fm, get_band_peak_fg
from fooof.utils import trim_spectrum
from fooof.data import FOOOFSettings
# Import the Bands object, which is used to define frequency bands
from fooof.bands import Bands
# Import plotting functions
from fooof.plts.spectra import plot_spectrum, plot_spectra
from fooof.plts.spectra import plot_spectrum_shading, plot_spectra_shading



# Change directory
os.chdir("C:\\Users\\ssshe\\Documents\\MathLab\\Analysis\\OrientWheel_Exo\\ML_Analysis")

# Location of saved data files
data_dir = os.path.join(os.getcwd(), 'fooof_data', 'byTargets_ML_v1') #single trial
# data_dir = os.path.join(os.getcwd(), 'fooof_AvgData', 'byTargets_ML_v1') #trial average

# Location to save fooof results
save_dir = os.path.join(os.getcwd(), 'fooof_data', 'byTargets_ML_v1', 'fooof_results') #single trial
# save_dir = os.path.join(os.getcwd(), 'fooof_AvgData', 'byTargets_ML_v1', 'fooof_results') #trial average

# Load the participant names
partname = sio.loadmat(os.path.join(data_dir, 'participants.mat'))
partlist = partname['nametmp']
# partlist = partname['nametmp'][2:3] #selecting a couple participants

# Define frequency bands of interest
bands = Bands({'theta' : [4, 7],
               'alpha' : [8, 14],
               'beta1' : [15, 22]})

# List electrode labels
elect = ['Oz','Pz','Cz','FCz','Fz','O1','O2','PO3','PO4','P7','P8','P5','P6',
         'P3','P4','CP5','CP6','CP1','CP2','C3','C4','FC5','FC6','FC1','FC2',
         'F7','F8','F3','F4','Fp1','Fp2']


for i, part in enumerate(partlist[:]):
    # Get list of data files in folder
    mat_fname = os.listdir(os.path.join(data_dir, part))
    
    # Make new folder is one does not exist
    if not(os.path.isdir(os.path.join(save_dir, part))):
        os.mkdir(os.path.join(save_dir, part))
        
        
    for j, ifile in enumerate(mat_fname):    
        # Load the mat file 
        data = sio.loadmat(os.path.join(data_dir, part, ifile))
        
        # Unpack data from dictionary, and squeeze numpy arrays
        freqs = np.squeeze(data['freqs']).astype('float')
        psds = np.squeeze(data['psd']).astype('float')
        # ^Note: this also explicitly enforces type as float (type casts to float64, instead of float32)
        #  This is not strictly necessary for fitting, but is for saving out as json from FOOOF, if you want to do that

        # Transpose power spectra, to have the expected orientation for FOOOF
        psds = psds.T
        
        # Initialize FOOOFGroup object
        fg = FOOOFGroup(max_n_peaks=7,peak_threshold=0.2,min_peak_height=0.22,peak_width_limits=[1,5])
        
        # Run FOOOF across all power spectra
        fg.fit(freqs, psds, [2, 50])
        
        # Fit the FOOOF model on all PSDs, and report
        # fg.report(freqs, psds, [2, 50])
        
        # Get all band peaks from a group of power spectrum models
        beta1 = get_band_peak_fg(fg, bands.beta1)
        alphas = get_band_peak_fg(fg, bands.alpha)
        thetas = get_band_peak_fg(fg, bands.theta)
        
        # Save out a specific FOOOF measure of interest - for example, slopes
        # exps = fg.get_params('aperiodic_params', 'exponent')
        # sio.savemat('exps.mat', {'exps' : exps})
        
        # Save out fooof results to json file
        #  There is a utility file to load this json file directly into Matlab
        savename = 'f_results_' + ifile[6:-4]
        fg.save(os.path.join(save_dir, part, savename), save_results=True)
        
        # Save out a specific FOOOF measure of interest
        savenameb = 'bands_' + ifile[6:-4] + '.mat'
        sio.savemat(os.path.join(save_dir, part, savenameb), {'beta1' : beta1,
                                                              'alphas' : alphas, 
                                                              'thetas' : thetas})
        
        # Clear variables
        del data, freqs, psds, savename, savenameb, ifile, fg, alphas, thetas, beta1 
    del mat_fname, part, j 
del i       
        


# Find the index of the worst model fit from the group
worst_fit_ind = np.argmax(fg.get_params('error'))
# Extract this model fit from the group
fm = fg.get_fooof(worst_fit_ind, regenerate=True)     

# Check results and visualize the extracted model
fm.plot() 
fm.print_results()
   


# Set whether to plot in log-log space
plt_log = True

# Plot the power spectrum
nfg.plot(plt_log)

plot_spectra(fg.freqs, fg.power_spectra, plt_log)        
     
# Plot full model, created by combining the peak and aperiodic fits        
plot_spectrum(fm.freqs, fm.fooofed_spectrum_, plt_log,
              label='Full Model', color='red')      
        
# Plot the peak fit: created by re-fitting all of the candidate peaks together
plot_spectrum(fm.freqs, fm._peak_fit, plt_log, color='green', label='Final Periodic Fit')


plot_spectrum(fm.freqs, fm._ap_fit, plt_log)

plot_spectrum(fm.freqs, (fm.power_spectrum - fm._ap_fit))


fg.fooofed_spectrum_
fg.power_spectra
fg._ap_fit


print(fm._ap_fit.__doc__)

