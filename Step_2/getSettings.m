function [ds] = getSettings

% ------------------------------------------------------------------------------------------------------
% Author: James Ives
% Email: james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% 
% This script was written by James Ives and is released under the GNU General Public License v3.0. 
% 
% You are free to redistribute and/or modify this script under the terms of the GNU General Public 
% License as published by the Free Software Foundation, either version 3 of the License, or (at 
% your option) any later version.
% 
% This script is provided "as-is" without any warranty; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details: https://www.gnu.org/licenses/gpl-3.0.html
% 
% I am happy to collaborate on any projects related to this script. 
% Feel free to contact me at the email addresses provided.
% -----------------------------------------------------------------------------------------------------

% The purpose of this function is the set all the settings and create a data structure that will be used later on.
% This version has been co-opted for STREAM. The original manages settings for a full preprocessing pipeline for EEG
% and eye tracking. Where results are 0, these are unused but required as part of the data structure.


% We'll store all the data, settings and other info in our datastructure (ds) in a struct so that it can be passed to functions easily and 
% they can use what they need from it.
ds = struct;


%% General
ds.settings.general.procEEG                         = 1; % Whether you want to process EEG or not 0 for no, 1 for yes
ds.settings.general.procET                          = 0; % Whether you want to process ET or not, 0 for no, 1 for yes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Paths
% First we're going to set the paths, it's expected that there will be a rootpath, within the rootpath there should be a data folder, within the
% data folder each participant should have it's own folder. Data will be saved in the raw folder and preprocessed folders
% General paths
ds.settings.paths.rootPath                          = fullfile('E:\Birkbeck\STREAM INDIA\Datasets');
[ds.settings.paths.rootPath]                        = checkPathEnd(ds.settings.paths.rootPath);

if ds.settings.general.procEEG
    % EEG Paths
    ds.settings.paths.dataPath                      = fullfile(ds.settings.paths.rootPath, '1. Raw\EEG\');
    ds.settings.paths.rawEEGPath                    = 0;
    ds.settings.paths.noTrialsFoundPath             = 0; % strcat(ds.settings.paths.rawEEGPath, 'No trials found/')
    ds.settings.paths.epochedEEGPath                = fullfile(ds.settings.paths.rootPath, '2. Preprocessed\2.1 Epoched_EEG\');
    ds.settings.paths.autoRejectEEGPath             = 0; % strcat(ds.settings.paths.epochedEEGPath, 'Auto rejected/')
    ds.settings.paths.preprocEEGPath                = 0;
    ds.settings.paths.concatEEGPath                 = 0;
    ds.settings.paths.icaEEGPath                    = 0;

    pathList = {ds.settings.paths.rawEEGPath; ds.settings.paths.noTrialsFoundPath; ds.settings.paths.epochedEEGPath; ...
        ds.settings.paths.autoRejectEEGPath; ds.settings.paths.preprocEEGPath; ds.settings.paths.concatEEGPath; ...
        ds.settings.paths.icaEEGPath};

    % Checks that the above folders have been created and if not creates them for you, so nothing crashes later on.
    checkAndCreateFolders(pathList);
    % Checks that there is a backslash at the end of each path, otherwise the folder structure could get confused.
    for i = 1:length(pathList); checkPathEnd(pathList{i}); end
end

if ds.settings.general.procET
    % ET Paths
    ds.settings.paths.rawETPath                     = strcat(ds.settings.paths.rootPath, '1.2 Raw_ET/');
    ds.settings.paths.epochedETPath                 = strcat(ds.settings.paths.rootPath, '2.2 Epoched_ET/');
    ds.settings.paths.fixationETPath                 = strcat(ds.settings.paths.rootPath, '3.2 Fixations_ET/');
    ds.settings.paths.concatETPath                  = strcat(ds.settings.paths.rootPath, '4.2 Concatenated_ET/');

    pathList = {ds.settings.paths.rawETPath; strcat(ds.settings.paths.rawETPath, 'No trials found/'); ds.settings.paths.epochedETPath; ds.settings.paths.fixationETPath; ...
        ds.settings.paths.concatETPath};

    % Checks that the above folders have been created and if not creates them for you, so nothing crashes later on.
    checkAndCreateFolders(pathList);
    % Checks that there is a backslash at the end of each path, otherwise the folder structure could get confused.
    for i = 1:length(pathList); checkPathEnd(pathList{i}); end
end

if ds.settings.general.procEEG & ds.settings.general.procET
    ds.settings.paths.lookingEEGPath                = strcat(ds.settings.paths.rootPath, '6.1 EEG_Filtered_with_ET/');
    ds.settings.paths.analysedDataPath              = strcat(ds.settings.paths.rootPath, '7.1 Analysed_Data/');

    pathList = {ds.settings.paths.lookingEEGPath; strcat(ds.settings.paths.lookingEEGPath, 'Filt_By_Null/'); ...
        strcat(ds.settings.paths.lookingEEGPath, 'Filt_By_Fix/')};

    % Checks that the above folders have been created and if not creates them for you, so nothing crashes later on.
    checkAndCreateFolders(pathList);
    % Checks that there is a backslash at the end of each path, otherwise the folder structure could get confused.
    for i = 1:length(pathList); checkPathEnd(pathList{i}); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Events
% If epochByEvents is set to auto then specify which event numbers signal the start and end of trials, these will be used to epoch the data later on.
ds.settings.epochByEvents                           = 'auto'; % Set this to manual if you would rather do it yourself
ds.settings.epochWithDINMarker                      = 0; % Sometimes the markers have "DIN"
ds.settings.onOffsetEventNumbers                    = [[58 59]; ...
                                                    [450 470]; [451 470]; [452 470]; [453 470]; ...
                                                    [60 61]; [84 85]; [82 83]; [80 81]; ...
                                                    [92 93]; [94 95]; [90 91]; ...
                                                    [21 19]; [16 19]; [11 19]; ...
                                                    [10 19]; [12 19]; [13 19]; [14 19]; ...
                                                    [15 19]; [17 19]; [18 19]; ...
                                                    [1100 1101]];
% Start and end of the trial only, these need to be a list of pairs e.g. [[U V]; [W X]; [Y Z]]
% You should name all trial types and this should match the number of event pairs you've put above. If you want repeats you can repeat and name them
% something different.

% Some of these are repeated because there is a Malawian, Indian and other counterpart, but either way you would know which is which.
ds.settings.eventNames                              = {'Between_Vid'; ...
                                                    'Face_SSVEP_FL6Hz_CR7_5Hz'; 'Face_SSVEP_FL7_5Hz_CR6Hz'; 'Face_SSVEP_CL6Hz_FR7_5Hz'; 'Face_SSVEP_CL7_5Hz_FR6Hz'; ...
                                                    'Gap'; 'Gap'; 'Gap'; 'Gap_Faces_English'; ...
                                                    'Reading'; 'Reading'; 'Reading_English'; ...
                                                    'Rest_Vid_Face_Onset'; 'Rest_Vid_Face_Onset'; 'Rest_Vid_Toy_Onset'; ...
                                                    'Rest_Vid_EN_Face'; 'Rest_Vid_NL_Face'; 'Rest_Vid_SW_Face'; 'Rest_Vid_PL_Face'; ...
                                                    'Rest_Vid_GM_Face'; 'Rest_Vid_ES_Face'; 'Rest_Vid_FR_Face'; ...
                                                    'Rocket'};

% Some trials are too short and are spread throughout the task list. Here we can merge these trials with pop_mergeset so that there is a longer
% continuous trial.
ds.settings.trialNumsToMerge                        = [[58 59]; ...
                                                       [60 61]; [84 85]; [82 83]; [80 81]];
ds.settings.trialNamesToMerge                       = {'Between_Vid'; ...
                                                    'Gap'; 'Gap_Faces_English'};

ds.settings.checkTrialLength                        = 0; % Whether we want to check trial lengths or not
ds.settings.maxEventDiscrepency                     = 0; % The maximum discrepencies allowed between trials within a trial type in seconds. This is used as a check when epoching the data.
ds.settings.minTrialLength                          = 10; % Measured in seconds, the minimum trial length, used when epoching the data
ds.settings.maxTrialLength                          = 16; % Measured in seconds, the maximum expected trial length, used when epoching the data.
ds.settings.expectedTrialLength                     = 15; % Measured in seconds, the expected trial length, used when epoching the data.
    
%% EEG preprocessing
% EEG preprocessing settings
if ds.settings.general.procEEG
    ds.settings.eegPreproc.expectedEEGSampleRate    = 500; %1000; % Checked during epoching
    ds.settings.eegPreproc.expectednbChannels       = 23; % Checked during epoching
    ds.settings.eegPreproc.filterOrder              = 6600; % The filter order used for high and low passes
    ds.settings.eegPreproc.highpass                 = 0.5; % Highpass cutoff
    ds.settings.eegPreproc.lowpass                  = 40; % Lowpass cutoff
    ds.settings.eegPreproc.notch                    = 50; % Notch filter Hz
    ds.settings.eegPreproc.noisyChanThreshold       = 0.25; % The threshold as a percentage for the number of noisy channels that can be rejected and interpolated
    ds.settings.eegPreproc.noisySegmentsThreshold   = 0.25; % The threshold as a percentage for the number of noisy segments (in 1 second chunks) that can be too noisy
end

%% ET preprocessing
% Eye tracking preprocessing settings
if ds.settings.general.procET
    ds.settings.et.etSampleRate                     = 120;
    ds.settings.et.opt.xres                         = 1920; % maximum value of horizontal resolution in pixels
    ds.settings.et.opt.yres                         = 1080; % maximum value of vertical resolution in pixels
    ds.settings.et.opt.missingx                     = -ds.settings.et.opt.xres; % missing value for horizontal position in eye-tracking data (example data uses -xres). used throughout functions as signal for data loss
    ds.settings.et.opt.missingy                     = -ds.settings.et.opt.yres; % missing value for vertical position in eye-tracking data (example data uses -yres). used throughout functions as signal for data loss
    ds.settings.et.opt.freq                         = 120; % sampling frequency of data (check that this value matches with values actually obtained from measurement!)
    ds.settings.et.opt.downsampFilter               = 0; % Always 0

    % Variables for the calculation of visual angle
    % These values are used to calculate noise measures (RMS and BCEA) of fixations. The may be left as is, but don't use the noise measures then.
    ds.settings.et.opt.scrSz                        = [50.9174 28.6411]; % screen size in cm
    ds.settings.et.opt.disttoscreen                 = 65; % distance to screen in cm.
end

%% EEG and ET preprocessing
% If preprocessing both
if ds.settings.general.procEEG & ds.settings.general.procET
    % This filter window factor is the factor that divides into both the EEG and ET sampling rate. This is used to search the data by a set window and
    % remove EEG where the infant isn't looking. For example, a filter window factor of 4 with EEG and ET sampling rates of 1000 and 120Hz means that you
    % are dividing each second by 4, 1000/4 = 250 (ms) (and this divides into both sampling frequencies). Increase the filter window factor to increase
    % the resolution. This number MUST divide into both sampling frequencies.
    ds.settings.eegPreproc.etFilterWindow           = 2;

    % Checks that the above divides into each of the sampling rates and throws an exception if not
    if mod(ds.settings.eegPreproc.expectedEEGSampleRate, ds.settings.eegPreproc.etFilterWindow) ~= 0 & mod(ds.settings.et.etSampleRate , ds.settings.eegPreproc.etFilterWindow) ~= 0
        ME = MException('MATLAB:filterFactorError', 'ET filter window factor does not divide into the EEG and/or ET sampling rate');
        throw(ME)
    end
end

end
