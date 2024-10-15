# STREAM-Analyses
## Background
Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com<br />
For: STREAM project 
-	https://www.bbk.ac.uk/school/psychological-sciences/research/scalable-transdiagnostic-early-assessment-of-mental-health
-	https://research.reading.ac.uk/stream/braintools/
<!-- end of the list -->
Date: 14th October 2024<br />
Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html<br />
Open to collaboration—feel free to contact me!<br />
## Goal
Process the EEG data from raw and complete some basic analyses (theta power/FOOOF).
## Expectations
This script expects that there will be a single folder with a subfolder for each new participant. No subfolders of folders containing participant data will be processed without changes to the code in step 1. Within each folder should be a time-stamped subfolder showing a data collection attempt. Within each time-stamped subfolder should be one set of EEG and one set of ET data.</p>
All EEG data is expected to be in a .easy format, ideally each .easy file will be paired with a corresponding .info file. All eyetracking (ET) data is expected to be in a further subfolder called “eye tracking”. This should be in a .mat format with all data in a Task Engine based uint8 format.</p>
Tracker and diary files may be present but will not be used within this pipeline.</p>
These scripts largely use EEGLAB and Task Engine, it is assumed that this has been installed. For a copy of Task Engine please email James. Most of the scripts are in Matlab but the fooof.py script is written in Python. The user is required to have Python installed (written in Python 3.10.9) with dependencies pandas, numpy, matplotlib and fooof.
## How to run these analyses
### Step 1
Converts EEG data from .easy into a .mat EEGLAB format. Converts ET data from a Task Engine uint8 format into a standard struct format. Where there are no events in the EEG data, these are reconstructed from the ET data.</p>
To run: read and change the appropriate variables (especially path names) within ET_EEG_Maintenance
### Step 2
Epochs EEG by task, preprocesses the EEG data in 7 stages, segments the data and saves a copy ready for FOOOF analysis, calculates an FFT and saves a copy ready for frequency power analysis.</p>
To run: 
-	Change paths and settings in getSettings as needed.
-	Change paths and settings in quick_bulk_preprocessing
-	Run quick_bulk_preprocessing
<!-- end of the list -->
Assumes that step 1 has run or data is in a similar format.
### Step 3
Run simple analyses that give the theta band frequency power and fooof results.</p>
To run: 
-	Change paths and settings in first_analyses.m/fooof.py
-	Run first_analyses.m and fooof.py
<!-- end of the list -->
Assumes that steps 1 and 2 have run or data is in a similar format.
