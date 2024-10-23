#! python3
# fooooooooof.py

import gc, os, sys
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from fooof import FOOOF
from fooof.utils.download import load_fooof_data
from fooof.plts.spectra import plot_spectra
from fooof.plts.spectra import plot_spectrum
from fooof.plts.annotate import plot_annotated_model

# ------------------------------------------------------------------------------------------------------
# Author: James Ives
# Email: james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
# Date: 14th October 2024
# 
# This script was written by James Ives and is released under the GNU General Public License v3.0. 
# This script is largely written with help from the fooof tutorial, which can be found:
# https://fooof-tools.github.io/fooof/auto_tutorials/index.html
#
# You are free to redistribute and/or modify this script under the terms of the GNU General Public 
# License as published by the Free Software Foundation, either version 3 of the License, or (at 
# your option) any later version.
# 
# This script is provided "as-is" without any warranty; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details: https://www.gnu.org/licenses/gpl-3.0.html
# 
# I am happy to collaborate on any projects related to this script. 
# Feel free to contact me at the email addresses provided.
# -----------------------------------------------------------------------------------------------------

# This script largely follows the fooof tutorial, which can be found here:
# https://fooof-tools.github.io/fooof/auto_tutorials/index.html
# This script expects that input data are a folder containing csv files. Each csv file should have
# two columns, the first is the frequency scale, and the second shows the corresponding power in
# each of the frequency bins.

# Currently the script excludes the first 4 values related to data under 1Hz, this is because the
# preprocessing scripts remove a lot of the data under 1Hz when removing low-frequency drift data.

# There is currently a bug. After too many files are process a bitmap error occurs suggesting that
# the memory has been filled while producing figures. To solve this either start with a smaller
# number of files with processed or comment out the sections that produce figures.
# Bug fix should come soon.

# Potential directory paths.
# directory = 'E:\\Birkbeck\\STREAM\\Epoched_full_spectral\\'
# directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4.2 FFT avg-segs 1s\\2s 50% overlap whole head'
# directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.3 Segmented Whole Head Averaged\\'
directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.4 Segmented Regional Averaged\\'
# directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.4 Segmented Regional Averaged\\Rest_Vid\\'
# directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.1.1 Full Whole Head Averaged'
# directory = 'E:\\Birkbeck\\STREAM\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.1.2 Full Regional Averaged'
# directory = 'E:\\Birkbeck\\STREAM INDIA\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.3 Segmented Whole Head Averaged\\'
# directory = 'E:\\Birkbeck\\STREAM INDIA\\Datasets\\2. Preprocessed\\2.4 FFT\\2.4.4 Segmented Regional Averaged\\'
results = []

# Conditions currently known
condition_list = ['Aud_gate', 'Between_Vid', 'Face_SSVEP_CL6Hz_FR7_5Hz', 'Face_SSVEP_CL7_5Hz_FR6Hz',\
                  'Face_SSVEP_FL6Hz_CR7_5Hz', 'Face_SSVEP_FL7_5Hz_CR6Hz','fast_erp', 'Gap',\
                  'Gap_Faces_English', 'Reading', 'Reading_English', 'Rest_Vid_Face_Onset',\
                  'Rest_Vid_Toy_Onset', 'Rest_Vid_EN_Face', 'Rest_Vid_NL_Face', 'Rest_Vid_SW_Face',\
                  'Rest_Vid_PL_Face','Rest_Vid_GM_Face', 'Rest_Vid_ES_Face', 'Rest_Vid_FR_Face']

region_list = ['central', 'frontal', 'occipital', 'parietal']
# Can be either "Whole head" or "Regional"
regional = 'Regional'

# Set up the figure object here to reuse and save space
fig = plt.figure()
ax = fig.add_subplot()

for file_index, filename in enumerate(os.listdir(directory)):

    # Used to skip files up to a specific number
    if file_index < 0:
        continue
    
    if not filename.endswith('.csv'):
        continue  # Skip to the next iteration if the filename doesn't end with '.csv'
    
    print(f'Processing: {filename}')

    # Check the conditions and regions
    condition = [c for c in condition_list if c in filename][0]
    region = [r for r in region_list if r in filename]
    if len(region) > 0:
        region = region[0]
    else:
        region = 'Whole_head'
        
    filepath = os.path.join(directory, filename)

    # Use context manager to read the CSV file
    with open(filepath, 'r') as file:
        data = pd.read_csv(file)

    # Extract the frequency and power spectrum data
    # Since implementing the clean_drifts script in matlab there is much less power in the <1Hz range
    # which throws this off, the first 5 values are all < 1Hz so are exlcuded here.
    freqs = data.iloc[4:, 0].values
    powers = data.iloc[4:, 1:].values

    # Loop over each electrode (column 2 to 21 in the CSV file)
    for electrode in range(powers.shape[1]):
        power_spectrum = powers[:, electrode]

        # Initialize the FOOOF model
        fm = FOOOF(peak_width_limits=[1, 6], max_n_peaks=6, min_peak_height=0.1, aperiodic_mode='fixed')

        # Fit the model to the power spectrum
        fm.fit(freqs, powers[:, electrode])

        # Extract the parameters
        aperiodic_params = fm.aperiodic_params_
        peak_params = fm.peak_params_
        r_squared = fm.r_squared_
        error = fm.error_

        peak_frequency = []
        peak_amplitude = []
        peak_width = []

        for i, peak in enumerate(peak_params):
            peak_frequency.append(peak[0])
            peak_amplitude.append(peak[1])
            peak_width.append(peak[2])
            
        # Save the results in a dictionary
        result = {
            'Filename': filename,
            'Condition': condition,
            'Region': region,
            'Offset': aperiodic_params[0],
            'Exponent': aperiodic_params[1],
            'R_squared': r_squared,
            'Error': error,
            'Peak_Freqs': peak_frequency,
            'Peak_Amp': peak_amplitude,
            'Peak_Width': peak_width
        }
        results.append(result)

        # Plot the original and the FOOOF model

        #plot_spectrum(freqs, powers[:, electrode], log_powers=True, label='Original Spectrum')
        fm.plot(ax=ax)
        f_name = filename[:-4]
        plt.title(f'FOOOF Model - {f_name} - Region {region}')
        plt.savefig(f'{directory}\\Figures\\FOOOF_Model_{f_name}_region_{region}.png')
#         plt.close('all')

        # Clear the figure from memory
#         plt.clf()
        plt.cla()
        gc.collect() 	 # Garbage collect to free up memory
    
    if file_index % 100 == 0:
        # Convert the results to a DataFrame
        results_df = pd.DataFrame(results)

        # Save the results to a new CSV file
        results_df.to_csv(f'{directory}\\FOOOF_Results\\FOOOF_Model_{regional}.csv', index=False)
            

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)

# Save the results to a new CSV file
results_df.to_csv(f'{directory}\\FOOOF_Results\\FOOOF_Model_{regional}.csv', index=False)

print("FOOOF analysis completed and results saved.")
        
