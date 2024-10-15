function first_analyses

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

% The function of this script is to run the first analyses and some other useful functions. 
% If the data have not yet been epoched or segmented and averaged (or the goal is to do this on raw data)
% this is the first two functions. Then calculates the fft, it can give the whole spectrum for all electrodes
% or averaged across the whole head, calculates theta 4-7Hz and can plot the spectral graphs. This script can
% also be used to copy eye tracking files from subfolders to a main folder.

%% Paths
addpath('E:\Birkbeck\Scripts\Stream\');
addpath('E:\Birkbeck\Scripts\Stream\Theta\');
addpath(genpath('E:\Birkbeck\Scripts\James Common\'));

root_path = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\';
output_averaged_spectral_folder = fullfile(root_path, 'Overview data\');
output_folder_ET = 'E:\Birkbeck\STREAM\Datasets\1. Raw\ET\';


folders = {output_averaged_spectral_folder; ...
    output_folder_ET};

for folder = 1:length(folders)
    if ~exist(folders{folder}, 'dir')
        mkdir(folders{folder});
    end
end

%% Setup
% Setting switches for what the script will do.
epoch_eeg = 0;              % Epochs and saves EEG, should have been done in preprocessing
segment_and_average = 0;    % Segments the epoched EEG into smaller (usually 1s) segments and averages, also should have been completed in preprocessing
average_head = 0;           % Averages across the whole head
save_full_spectrum = 0;     % Saves the full spectrum rather than just theta to a folder
process_theta = 1;          % Processes and saves theta
plot_spectral = 0;          % Plots full spectrum and saves
copy_ET_files = 0;          % Copies relevant ET files.

% Other settings
segment_length = 2;         % seconds
overlap_percent = 0.5;      % Percent
% Calculate step size for the epochs
step_size = segment_length * (1 - overlap_percent); % 1 second step with 50% overlap
max_hz  = 100;

%% Epoch EEG
if epoch_eeg; ds = getSettings(); epochEEG(ds); end

% To get these run the bulk preprocessing
files = dir(fullfile(root_path, '2.2 Preprocessed_EEG\'));
movedParticipants = {}; % List to keep track of moved participant files
all_avg_theta = [];
theta_table = table;

% Preallocate tables for social and non-social conditions
ppt_list = {};
condition_list = {};
region_list = {};
theta_power_list = [];

% Channel labels
channel_labels = {'T7', 'P4', 'Cz', 'Pz', 'P3', 'P8', 'Oz', 'O2', 'T8', ...
                  'PO8', 'C4', 'F4', 'AF8', 'Fz', 'C3', 'F3', 'AF7', 'P7', ...
                  'PO7', 'Fpz', 'x', 'y', 'z'};

% Indices for brain regions
frontal_channels = {'F3', 'Fz', 'F4'};
central_channels = {'C3', 'Cz', 'C4'};
parietal_channels = {'P3', 'Pz', 'P4'};
occipital_channels = {'PO7', 'Oz', 'PO8'};

% Converts channel locations to indices
frontal_indices = find(ismember(channel_labels, frontal_channels));
central_indices = find(ismember(channel_labels, central_channels));
parietal_indices = find(ismember(channel_labels, parietal_channels));
occipital_indices = find(ismember(channel_labels, occipital_channels));

%% Spectral loop
for file = 3:length(files)-1
    load(fullfile(files(file).folder, files(file).name));
    Fs = EEG.srate;

    % If we are to segment and average, this segments the data into smaller (usually 1s) segments and averages before the data are put through the FFT
    % All data from this point continue in this way.
    if segment_and_average
        % Segment samples
        segment_samples = Fs*segment_length;
        overlap_samples = segment_samples*overlap_percent;
        % Find the duration of the data in seconds
        data_duration = EEG.pnts / EEG.srate;
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

        % % Average the segmented data, ignoring the inserted nans
        data = squeeze(mean(segmented_EEG, 1, 'omitnan'));
        
        % Resubmit the data back into the EEG data structure
        EEG.times = times;
        EEG.data = data;
        EEG.pnts = segment_samples;
        EEG.xmax = segment_length;
        comment = sprintf("This is data that has been segmented into %d second segments and averaged.", segment_length);
        EEG.info = [EEG.info; {comment}];
        EEG.comments = [EEG.comments; {comment}];

        save(fullfile(output_avg_seg_folder, files(file).name), "EEG");
    end

    % Extract participant name from filename
    [~, name, ~] = fileparts(files(file).name);
    parts = strsplit(name, '_');
    participantName = parts{2};
    condition = parts{3}; % Extract condition (Face_Onset or Toy_Onset)

    if average_head 
        EEG.data = mean(EEG.data, 1, 'omitnan');
    end

    fft = [];
    for fftElec = 1:size(EEG.data, 1)
        [fftRes, fftHzScale, dBfft] = myFFT(squeeze(EEG.data(fftElec, :)),Fs,0,50);
        [minValue, closestIndex] = min(abs(fftHzScale-max_hz)); % We're filtering at 100Hz so no point keeping it.
        fftRes = fftRes(1:closestIndex);
        fftHzScale = fftHzScale(1:closestIndex);
        dBfft = dBfft(1:closestIndex);

        fft(1, fftElec, :) = fftRes;
        fft(2, fftElec, :) = fftHzScale;
        fft(3, fftElec, :) = dBfft;
    end
    
    if save_full_spectrum

        freqs = squeeze(fft(2,1,:)); % Assuming first column is frequency
        
        % Extract the frequency and power spectrum data
        if average_head
            % whole_head_fft = mean(fft(1,:,:), 2, 'omitnan');
            % powers = squeeze(whole_head_fft(1,:, :));
            powers = squeeze(fft(1,1, :)); % Assuming columns 2 to end are power data
        else
            powers = squeeze(fft(1,1:20, :))'; % Assuming columns 2 to end are power data
        end

        if segment_and_average
            % If segmenting and averaging to 1 second segments then the frequency resolution is already close to 0.5Hz steps
            spectrum = [freqs, powers];

            writematrix(spectrum, fullfile(output_avg_seg_spectral_folder, [files(file).name(1:end-4), '.csv']))
        else
            % New frequencies to downsample the current data to
            new_freqs = 0:0.5:25;
            new_freqs = new_freqs';

            % Interpolate the power data to match the new frequency resolution
            new_powers = interp1(freqs, powers, new_freqs, 'linear');

            spectrum = [new_freqs, new_powers];

            % save(fullfile(outputFullSpectralFolder, [files(file).name(1:end-4), '.mat']), 'spectrum');
            writematrix(spectrum, fullfile(output_full_spectral_folder, [files(file).name(1:end-4), '.csv']))
        end      
    end

    if process_theta
    % Extract and save theta data
    thetaIndex = find(fftHzScale >= 4 & fftHzScale <= 7);
    thetaData = fft(1:2, :, thetaIndex);

    if segment_and_average
        save(fullfile(output_avg_seg_theta_spectral_folder, [files(file).name(1:end-4), '_theta.mat']), 'thetaData');
    else
        save(fullfile(output_theta_spectral_folder, [files(file).name(1:end-4), '_theta.mat']), 'thetaData');
    end

    % Average theta
    thetaOnly = squeeze(thetaData(1, :, :));
    avg_theta = squeeze(mean(thetaOnly,  2, 'omitnan'));
    logTheta = squeeze(log10(avg_theta));

    % Calculate averages for each brain region
    frontal_avg = mean(logTheta(frontal_indices), 'omitnan');
    central_avg = mean(logTheta(central_indices), 'omitnan');
    parietal_avg = mean(logTheta(parietal_indices), 'omitnan');
    occipital_avg = mean(logTheta(occipital_indices), 'omitnan');

    % Store data in lists
    % ppt_list = [ppt_list; {participantName}];
    if strcmp(condition, 'Face')
        condition = {'Social'};
    elseif strcmp(condition, 'Toy')
        condition = {'Non-social'};
    else
        input('Condition not correct, check and reset', 's')
    end
    % theta_power_list = [theta_power_list; logTheta];

    % Append regional data to the table
    ppt_list = [ppt_list; {participantName}; {participantName}; {participantName}; {participantName}];
    condition_list = [condition_list; condition; condition; condition; condition];
    region_list = [region_list; {'Frontal'; 'Central'; 'Parietal'; 'Occipital'}];
    theta_power_list = [theta_power_list; frontal_avg; central_avg; parietal_avg; occipital_avg];
    end

    if plot_spectral
    % Plot the images and save them in the image folder
    fig = figure('Visible', 'off', 'Color', 'white'); % Set the figure to not be visible and background color to white
    if average_head
        plot(fftHzScale, squeeze(fft(1,1,:)));
    else
    for i = 1:size(fft, 2)
        plot(fftHzScale, squeeze(fft(1,i,:))); hold on;
    end
    end
    xlim([1 50]); title(files(file).name, 'Interpreter', 'none');

    % Save the figure
    [~, name, ~] = fileparts(files(file).name);

    if segment_and_average
        saveas(fig, fullfile(output_avg_seg_image_folder, [files(file).name(1:end-4), '.png'])); % Save as PNG file
    else
        saveas(fig, fullfile(output_image_folder, [files(file).name(1:end-4), '.png'])); % Save as PNG file
    end
    
    close(fig); % Close the figure to free up memory
     
    end
    
    if copy_ET_files

    % Check if participant's eyetracking files have already been moved
    if ismember(participantName, movedParticipants)
        continue; % Skip if the participant's files have already been moved
    end
    % Search and copy eyetracking files
    etFolder = fullfile('E:\Birkbeck\STREAM\data download MW', '**', 'decoded_eyetracking');
    etFiles = dir(fullfile(etFolder, ['ET_data_' participantName '.mat']));
    
    for etFile = 1:length(etFiles)
        copyfile(fullfile(etFiles(etFile).folder, etFiles(etFile).name), output_folder_ET);
    end

    movedParticipants{end+1} = participantName;
    end
end

if process_theta
    % save(fullfile(outputAveragedSpectralFolder, ['all_ppts_theta.mat']), 'all_avg_theta');

    % Create the table
    raw_dataset = table(ppt_list, condition_list, region_list, theta_power_list, ...
                        'VariableNames', {'ID', 'Condition', 'Region', 'AveragedLog10ThetaPower'});

    % Aggregate data by participant, condition, and region
    [uniqueKeys, ~, idx] = unique(raw_dataset(:, {'ID', 'Condition', 'Region'}), 'rows');
    averagedThetaPower = accumarray(idx, raw_dataset.AveragedLog10ThetaPower, [], @mean);

    % Create the final aggregated table
    theta_table = table(uniqueKeys.ID, uniqueKeys.Condition, uniqueKeys.Region, averagedThetaPower, ...
                           'VariableNames', {'ID', 'Condition', 'Region', 'AveragedLog10ThetaPower'});

    % Save the table
    save(fullfile(output_averaged_spectral_folder, 'theta_table.mat'), 'theta_table');
    % Save the table as a CSV file
    writetable(theta_table, fullfile(output_averaged_spectral_folder, 'theta_table.csv'));
end


end