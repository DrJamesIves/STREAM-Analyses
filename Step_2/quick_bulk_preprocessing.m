function [success, eeg_data] = quick_bulk_preprocessing(eeg_data)

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

% This script just does a quick preprocessing with general things that we may want across all files. This script
% epochs the data (optional), corrects low frequecy drifts, performs a highpass, notch and lowpass filter, computes 
% and removes a robust average, searches for noisy and bridged (optional) channels, any found are interpolated, 
% then it searches for noisy segments and replaces them with zeros. Afterwards it segments the data and averages, 
% needed for FOOOF analyses and similar.

% Important: check the paths provided and that the settings used match the type of analysis that you would like to run
% Important: Paths related to the .loc and .elp file on line (approx) 151 need to be updated.
% Requires eeglab.

success = [];
dataset = 'Malawi';
epoch_eeg = 0;

%% Paths
% Add path to the theta scripts (including task engine).
addpath(genpath('E:\Birkbeck\Scripts\Stream\Theta'));
addpath('E:\Birkbeck\Scripts\Stream\');
addpath(genpath('E:\Birkbeck\Scripts\James Common\'));

%% Epoch EEG
% If you need to epoch the EEG you will need to manually change the getSettings file to match the location of the files and which tasks you would like
% to epoch
if epoch_eeg
    in = input('Have you changed getSettings file for the EEG trial epoching?', 's');
    ds = getSettings();
    epochEEG(ds);
end

switch dataset
    case 'Malawi'
        root = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\';
        root_path = fullfile(root, '2.1 Epoched_EEG\');
        outpath_full = fullfile(root, '2.2 Preprocessed_EEG\2.2.1 Full\');
        outpath_segmented = fullfile(root, '2.2 Preprocessed_EEG\2.2.2 Segmented\');
                rejected_path = fullfile(root, '2.2 Preprocessed_EEG\2.2.3 Rejected data\');
        outpath_fft_full = fullfile(root, '2.4 FFT\2.4.1 Full\');
        outpath_fft_segmented = fullfile(root, '2.4 FFT\2.4.2 Segmented\');
        outpath_fft_segmented_whole_head_avg = fullfile(root, '2.4 FFT\2.4.3 Segmented Whole Head Averaged\');
        outpath_fft_segmented_regional_avg = fullfile(root, '2.4 FFT\2.4.4 Segmented Regional Averaged\');

    case 'India'
        root = 'E:\Birkbeck\STREAM INDIA\Datasets\2. Preprocessed\';
        root_path = fullfile(root, '2.1 Epoched_EEG\');
        outpath_full = fullfile(root, '2.2 Preprocessed_EEG\2.2.1 Full\');
        outpath_segmented = fullfile(root, '2.2 Preprocessed_EEG\2.2.2 Segmented\');
                rejected_path = fullfile(root, '2.2 Preprocessed_EEG\2.2.3 Rejected data\');
        outpath_fft_full = fullfile(root, '2.4 FFT\2.4.1 Full\');
        outpath_fft_segmented = fullfile(root, '2.4 FFT\2.4.2 Segmented\');
        outpath_fft_segmented_whole_head_avg = fullfile(root, '2.4 FFT\2.4.3 Segmented Whole Head Averaged\');
        outpath_fft_segmented_regional_avg = fullfile(root, '2.4 FFT\2.4.4 Segmented Regional Averaged\');
end

% Grab the files
files = dir(root_path);
filename_extras = '';

% Segmentation settings
segment_length = 2;         % seconds
overlap_percent = 0.5;      % Percent
% Calculate step size for the epochs
step_size = segment_length * (1 - overlap_percent);

% Preprocessing settings
highpass = 0.1;             % Hz
lowpass = 48;               % Hz
notch = 50;                 % Hz

% Other settings
min_trial_length = 10;      % seconds
chan_rej_threshold = 0.25;  % percent; noisy channel rejection threshold
seg_rej_threshold = 0.25;   % percent; noisy segment rejection threshold

% Notes that will be included within the fft file as a read me for future users on the data structure.
fft_readme = '1st dimension, 1 = fft, 2 = fft scale, 3 = dBfft. 2nd dimension, electrodes. 3rd dimension, time';

% Channel labels
channel_labels = {'T7', 'P4', 'Cz', 'Pz', 'P3', 'P8', 'Oz', 'O2', 'T8', ...
                  'PO8', 'C4', 'F4', 'AF8', 'Fz', 'C3', 'F3', 'AF7', 'P7', ...
                  'PO7', 'Fpz', 'x', 'y', 'z'};

% Indices for brain regions
frontal_channels = {'F3', 'Fz', 'F4'};
central_channels = {'C3', 'Cz', 'C4'};
parietal_channels = {'P3', 'Pz', 'P4'};
occipital_channels = {'PO7', 'Oz', 'PO8'};

% Converting channels to indices
frontal_indices = find(ismember(channel_labels, frontal_channels));
central_indices = find(ismember(channel_labels, central_channels));
parietal_indices = find(ismember(channel_labels, parietal_channels));
occipital_indices = find(ismember(channel_labels, occipital_channels));

% Loop through all files found
for file = 1:length(files)
    % Grab the filename and make sure that it ends with .mat
    eeg_file = strcat(files(file).folder, '\', files(file).name);
    if ~endsWith(eeg_file, '.mat')
        continue
    end

    % Uncomment to only focus on one trial type
    % if ~contains(files(file).name, 'Rest_Vid')
    %     continue
    % end

    disp(fprintf('Processing ...%s', eeg_file))

    % Load the eeg file and grab the sampling rate
    load(eeg_file);
    Fs = EEG.srate;

    % Chech that the file is above the minimum length, if not reject and save
    if size(EEG.data, 2) < Fs * min_trial_length
        save(fullfile(rejected_path, files(file).name), 'EEG', 'dataInfo');
        continue
    end

    % Checks whether the correct set of channel labels was loaded, if Ch1 then it is default and needs
    % to be fixed, assumes that the channels are in the right order.
    if strcmp(EEG.chanlocs(1).labels, 'Ch1')
        for i = 1:length(EEG.chanlocs)
            EEG.chanlocs(i).labels = channel_labels{i};
        end
    end

    % Creates a copy of the EEG file otherwise the next step crashes
    EEG2 = EEG;

    % Add in channel coordinates
    EEG = pop_select(EEG2, 'channel',{'T7', 'P4', 'Cz', 'Pz', 'P3', 'P8', 'Oz', 'O2', 'T8', 'PO8', 'C4', 'F4', 'AF8', 'Fz','C3', 'F3','AF7', 'P7', 'PO7', 'Fpz'}); %'2-EXG1' '2-EXG2' '2-EXG3' '2-EXG4' '2-EXG5' '2-EXG6' '2-EXG7' '2-EXG8'
    EEG = pop_chanedit(EEG, 'lookup', 'E:\Birkbeck\Scripts\standard-10-5-cap385.elp','load',{'E:\Birkbeck\Scripts\Enobio20.loc' 'filetype' 'autodetect'});

    % If no channels are returned or no data then reject the file and save.
    if EEG.nbchan == 0 | isempty(EEG.data)
        save(fullfile(rejected_path, files(file).name), 'EEG', 'dataInfo');
        continue
    end

    %% Preprocessing
    %% 0. Clean out low frequency drifts from the data
    % Uses default settings, quite a few channels are overwhelmed with low-freq drifts that cause a lot of problems.
    EEG = clean_drifts(EEG);

    %% 1. Highpass the EEG
    % This highpasses the EEG signal using the highpass and filter order setting in ds.
    EEG = pop_eegfiltnew(EEG, [], highpass, 8250, 1, [], 0);

    %% 2. Notch filter the EEG
    % Uncomment if lowpass is above 50 and notch filter is needed.
    % This is used to remove a specific frequency, which usually represents the electrical frequency where the data
    % were recorded.
    % EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',(1:EEG.nbchan) ,'computepower',1,'linefreqs', notch,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',[],'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',3,'winstep',3);

    %% 3. Lowpass the EEG
    % This lowpasses the EEG signal using the highpass and filter order setting in ds.
    EEG = pop_eegfiltnew(EEG, [], lowpass, 8250, 0, [], 0);

    %% 4. Robust average reference
    % First identify and temporarily remove bad channels, compute and average and remove that average from all channels. dataInfo will be stored with the
    % saved data at the end with relevant info. clean_channels doesn't take data with multiple trials, so first we create a copy and feed in one trial
    % at a time. This may mean that some trials have more rejected channels than others, which is why we're keeping this info.

    % This will be used to store the data quality info for later analysis
    [dataInfo.chansRemovedFromRobustAvg, dataInfo.channelsRejectedForNoise, dataInfo.channelsRejectedForBridging, dataInfo.rejectedSections] = deal([]);

    % Identify noisy channels
    try
        [EEGOUT, removedChans] = clean_channels(EEG);
    catch
        % If there are too many noisy channels it will crash, this just moves that file to the auto reject folder and moves to the next iteration of the loop.
        % movefile(strcat(files(file).folder, '\', files(file).name), strcat(files(file).folder, '\Auto rejected\', files(file).name))
        disp('## Too many noisy channels for robust average')
        success = [success; 0];
        continue
    end

    % Average the clean channels and remove the average from all channels
    EEG.data = EEG.data-mean(EEGOUT.data,1);
    % Clear EEGOUT variable we no longer need and store the number of channels removed from the robust average for later
    clear EEGOUT
    dataInfo.chansRemovedFromRobustAvg = [dataInfo.chansRemovedFromRobustAvg, removedChans];
    
    %% 5. Check for bad channels, bridging and interpolate
    % Use a slightly narrower check to identify bad channels that are either noisy or bridged.
    try
        [~,chans2interp] = clean_channels(EEG,0.7,4,[],[],[],[]); noisyChans = find(chans2interp==1);
    catch
        % If there are too many noisy channels it will crash, this just moves that file to the auto reject folder and moves to the next iteration of the loop.
        % movefile(strcat(files(file).folder, '\', files(file).name), strcat(files(file).folder, '\Auto rejected\', files(file).name))
        disp('## Too many noisy channels after robust average')
        success = [success; 0];
        continue
    end

    dataInfo.channelsRejectedForNoise = [dataInfo.channelsRejectedForNoise, chans2interp]; % Store these electrodes

    % Uncomment to include bridge channel detection. eBridge will need to be downloaded.
    % eBridge to checked for bridged channels without bothering to check for bridging from noisy channels
    % [EB ED] = eBridge(EEG, {EEG.chanlocs(noisyChans).labels}, 'PlotMode', 0); % , 'EpLength', 250
    % 
    % % Find the bridged channels
    % y = zeros(1, EEG.nbchan);
    % for chan = 1:size(EB.Bridged.Labels, 2)
    %     x = strcmp({EB.Bridged.Labels{chan}},{EEG.chanlocs(:).labels}); y = y + x;
    % end
    % bridgedChans = find(y==1);
    % dataInfo.channelsRejectedForBridging.EB = EB; % Store electrode bridge data
    % dataInfo.channelsRejectedForBridging.ED = ED; % Store general bridging info

    % If detecting noisy and bridged channels they need to be included here.
    toReject = [noisyChans]; %; bridgedChans'];

    % Reject noisy and bridged channels. If there are less channels to interpolate than the rejection threshold, then interpolate them to replace them.
    % Below is commented code to not reject noisy channel sets but instead replace the bad channels with nan vectors.
    % If there are more than the threshold of rejected channels then just reject those channels and replace them with NaNs (this keeps the structure
    % of the channels, which is important for spatial analyses later on).
    if length(toReject) > EEG.nbchan * chan_rej_threshold
        disp('## Too many noisy and bridged channels to reject')
        success = [success; 0];
        save(fullfile(rejected_path, files(file).name), 'EEG', 'dataInfo');
        continue
        % nanVector = nan(1, EEG.pnts);
        % for r = 1:length(toReject)
        %     EEG.data(toReject(r), :) = nanVector;
        % end
    else
        % This interpolates channels using nearest neighbours.
        EEG = eeg_interp(EEG, toReject, 'spherical');
    end

    %% 6. Remove bad sections of continuous data
    % If 70% of channels are bad get rid of data
    [signal,mask]=clean_windows(EEG,0.7);
    % Get column index of bad data segments and mark as 0 within EEG data
    [~,indx]=find(mask==0);EEG.data(:,indx)=0;
    dataInfo.rejectedSections = [dataInfo.rejectedSections; {mask'}];

    % Reject if there are too many sections of noisy data as a percentage of the whole trial (currently set at 25%). Otherwise save as normal and continue
    % on with fft and other parts.
    if EEG.pnts - sum(mask) > EEG.pnts * seg_rej_threshold
        disp('## Too many noisy noisy segments to reject')
        success = [success; 0];
        save(fullfile(rejected_path, files(file).name), 'EEG', 'dataInfo');
        continue
    else
        fprintf(strcat('Saving ...\n'))
        save(fullfile(outpath_full, files(file).name), 'EEG', 'dataInfo');
        success = [success; 1];
    end
    
    [fft, ~] = run_fft_and_plot(EEG, eeg_file, outpath_fft_full, 0);
    save(fullfile(outpath_fft_full, files(file).name), 'fft', 'dataInfo', "fft_readme");

    %% Segment and epoch
    % If we are to segment and average, this segments the data into smaller (usually 1s) segments and averages before the data are put through the FFT
    % All data from this point continue in this way.
    % Segment samples
    segment_samples = Fs*segment_length;
    overlap_samples = segment_samples*overlap_percent;
    % Find the duration of the data in seconds
    data_duration = EEG.pnts / Fs;
    % Generate time points for epochs
    epoch_start_times = 0:step_size:(data_duration - segment_length);
    % Initialize a structure to hold the segmented data
    segmented_EEG = zeros(length(epoch_start_times), size(EEG.data,1), segment_samples);

    % Loop over each epoch start time and create epochs
    for i = 1:length(epoch_start_times)
        % Define the start and end time for each epoch
        start_time = epoch_start_times(i);
        start_sample = (epoch_start_times(i)*Fs)+1;
        end_sample = (start_time + segment_length)*Fs;

        segmented_EEG(i, :, :) = EEG.data(:, start_sample:end_sample);
    end
    % % Select the first second of data as the times
    times = EEG.times(1:segment_samples);

    % Resubmit the data back into the EEG data structure
    EEG.times = times;
    EEG.data = segmented_EEG;
    EEG.pnts = segment_samples;
    EEG.xmax = segment_length;
    comment = sprintf("This is data that has been segmented into %d second segments and averaged.", segment_length);
    EEG.info = [EEG.info; {comment}];
    EEG.comments = [EEG.comments; {comment}];

    % Save the segmented data
    save(fullfile(outpath_segmented, files(file).name), 'EEG', 'dataInfo');

    % Rerun fft and save
    [fft, ~] = run_fft_and_plot(EEG, eeg_file, outpath_fft_segmented, 1);
    save(fullfile(outpath_fft_segmented, files(file).name), 'fft', 'dataInfo', "fft_readme");

    % Reformat and save in a csv format with two columns for FOOOF analysis
    fft2 = fft;
    fft = squeeze(mean(fft, 2, 'omitnan'));
    freqs = squeeze(fft(2,:))'; % Assuming first column is frequency
    powers = squeeze(fft(1,:))'; % Assuming columns 2 to end are power data
    spectrum = [freqs, powers];

    % Write whole head fft
    writematrix(spectrum, fullfile(outpath_fft_segmented_whole_head_avg, [files(file).name(1:end-4), '.csv']))

    % Frontal electrodes only
    fft = fft2(:, frontal_indices, :);
    fft = squeeze(mean(fft, 2, 'omitnan'));
    freqs = squeeze(fft(2,:))'; % Assuming first column is frequency
    powers = squeeze(fft(1,:))'; % Assuming columns 2 to end are power data
    spectrum = [freqs, powers];
    writematrix(spectrum, fullfile(outpath_fft_segmented_regional_avg, [files(file).name(1:end-4), '_frontal.csv']))

    % Central electrodes only
    fft = fft2(:, central_indices, :);
    fft = squeeze(mean(fft, 2, 'omitnan'));
    freqs = squeeze(fft(2,:))'; % Assuming first column is frequency
    powers = squeeze(fft(1,:))'; % Assuming columns 2 to end are power data
    spectrum = [freqs, powers];
    writematrix(spectrum, fullfile(outpath_fft_segmented_regional_avg, [files(file).name(1:end-4), '_central.csv']))

    % Parietal electrodes only
    fft = fft2(:, parietal_indices, :);
    fft = squeeze(mean(fft, 2, 'omitnan'));
    freqs = squeeze(fft(2,:))'; % Assuming first column is frequency
    powers = squeeze(fft(1,:))'; % Assuming columns 2 to end are power data
    spectrum = [freqs, powers];
    writematrix(spectrum, fullfile(outpath_fft_segmented_regional_avg, [files(file).name(1:end-4), '_parietal.csv']))

    % Occipital electrodes only
    fft = fft2(:, frontal_indices, :);
    fft = squeeze(mean(fft, 2, 'omitnan'));
    freqs = squeeze(fft(2,:))'; % Assuming first column is frequency
    powers = squeeze(fft(1,:))'; % Assuming columns 2 to end are power data
    spectrum = [freqs, powers];
    writematrix(spectrum, fullfile(outpath_fft_segmented_regional_avg, [files(file).name(1:end-4), '_occipital.csv']))
    
    clear EEG EEG2 epoch_start_times fft fft2 freqs mask powers segmented_EEG signal spectrum times

end

end
