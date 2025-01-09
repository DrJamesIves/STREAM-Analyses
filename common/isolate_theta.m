function isolate_theta

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!


% Simply takes in a csv with two columns, Hz scale in column 1 and power on the right, then segments out theta (4-7Hz)
% This is different from the script first_analyses as it takes in csvs where fft has already been calculated. So is quicker.

%% Paths
addpath('E:\Birkbeck\Scripts\Stream\');
addpath('E:\Birkbeck\Scripts\Stream\Theta\');
% addpath(genpath('E:\Birkbeck\Scripts\Stream\brainstorm-tools\'));
addpath(genpath('E:\Birkbeck\Scripts\James Common\'));

dataset = 'India';

switch dataset
    case 'Malawi'
        root = 'E:\Birkbeck\STREAM\Datasets\';
    case 'India'
        root = 'E:\Birkbeck\STREAM INDIA\Datasets\';
end

epoch_path = fullfile(root, '2. Preprocessed\2.1 Epoched_EEG');
preproc_path = fullfile(root, '2. Preprocessed\2.2 Preprocessed_EEG\2.2.1 Full');
data_path = fullfile(root, '3. Analysed\3.1 FFT\3.1.1 Full\');
output_averaged_spectral_folder = fullfile(root, '3. Analysed\3.2 Theta data\');
output_folder_ET = fullfile(root, '1. Raw\ET\');

% To save on clutter make the folders needed manually
folders = {output_averaged_spectral_folder; ...
    output_folder_ET};

for folder = 1:length(folders)
    if ~exist(folders{folder}, 'dir')
        mkdir(folders{folder});
    end
end

%% Setup
ds = getSettings(dataset);
lower_bound = 4;
upper_bound = 7;

% To get these tun the bulk preprocessing
files = dir(data_path);

all_avg_theta = [];
theta_table = table;

% Preallocate tables for social and non-social conditions
ppt_list = {};
condition_list = {};
region_list = {};
theta_power_list = [];
num_trials_list = [];

% Channel labels
channel_labels = {'T7', 'P4', 'Cz', 'Pz', 'P3', 'P8', 'Oz', 'O2', 'T8', ...
    'PO8', 'C4', 'F4', 'AF8', 'Fz', 'C3', 'F3', 'AF7', 'P7', ...
    'PO7', 'Fpz', 'x', 'y', 'z'};

% Indices for brain regions
frontal_channels = {'F3', 'Fz', 'F4'};
central_channels = {'C3', 'Cz', 'C4'};
parietal_channels = {'P3', 'Pz', 'P4'};
occipital_channels = {'PO7', 'Oz', 'PO8'};

frontal_indices = find(ismember(channel_labels, frontal_channels));
central_indices = find(ismember(channel_labels, central_channels));
parietal_indices = find(ismember(channel_labels, parietal_channels));
occipital_indices = find(ismember(channel_labels, occipital_channels));

fft = 0;

%% Spectral loop
for file = 3:length(files)-1
    load(fullfile(files(file).folder, files(file).name));

    % Find the indices in the fft scale corresponding to the boundaries of theta
    thetaIndex = find(fft(2,1,:) >= 4 & fft(2,1,:) <= 7);

    % Extract participant name from filename
    [~, name, ~] = fileparts(files(file).name);
    parts = strsplit(name, '_');
    participantName = parts{1};

    % Find and label the conditon, as there are multiple parts a direct search is needed
    switch parts{2}
        case {'Gap', 'Gap-merged', 'Reading', 'Reading-merged', 'Rocket', 'Rocket-merged'}
            condition = parts{2};
        case {'Aud', 'Between', 'fast'}
            condition = strjoin(parts(2:3), '_');
        case 'Face'
            condition = strjoin(parts(2:6), '_');
        case 'Rest'
            switch parts{4}
                case 'Face'
                    condition = 'Social_Rest_Vid';
                case 'Toy'
                    condition = 'Non-social_Rest_Vid';
            end
            % Differentiates merged and non-merged files
            if strcmp(parts{5}, 'Onset-merged')
                condition = [condition, '-merged'];
            % else
            %     condition = [condition, '_', parts{end}];
            end
    end

    if contains(condition, 'merged')
        switch condition
            case {'Between_Vid-merged', 'fast_erp-merged', 'Gap-merged'}
                temp = load(fullfile(epoch_path, files(file).name));
            case {'Social_Rest_Vid-merged', 'Non-social_Rest_Vid-merged', 'Rest_Vid_Face_Onset-merged', 'Rest_Vid_Toy_Onset-merged'}
                temp = load(fullfile(preproc_path, files(file).name));
            otherwise
                crashMe
        end
        num_trials = temp.EEG.num_trials_merged;
        clear temp
    else
        condition = [condition, '_', parts{end}];
        num_trials = 1;
    end

    whole_head_theta = mean(squeeze(mean(fft(1,:,thetaIndex), 2)),1);
    frontal_theta = mean(squeeze(mean(fft(1,frontal_indices,thetaIndex), 2)),1);
    central_theta = mean(squeeze(mean(fft(1,central_indices,thetaIndex), 2)),1);
    parietal_theta = mean(squeeze(mean(fft(1,parietal_indices,thetaIndex), 2)),1);
    occipital_theta = mean(squeeze(mean(fft(1,occipital_indices,thetaIndex), 2)),1);

    % Append regional data to the table
    ppt_list = [ppt_list; {participantName}; {participantName}; {participantName}; {participantName}; {participantName}];
    condition_list = [condition_list; condition; condition; condition; condition; condition];
    region_list = [region_list; {'Whole head'; 'Frontal'; 'Central'; 'Parietal'; 'Occipital'}];
    theta_power_list = [theta_power_list; whole_head_theta; frontal_theta; central_theta; parietal_theta; occipital_theta];
    num_trials_list = [num_trials_list; num_trials; num_trials; num_trials; num_trials; num_trials];
end

% Create the table
raw_dataset = table(ppt_list, condition_list, region_list, theta_power_list, num_trials_list, ...
    'VariableNames', {'ID', 'Condition', 'Region', 'AveragedThetaPower', 'Num merged trials'});

% Aggregate data by participant, condition, and region
[uniqueKeys, ~, idx] = unique(raw_dataset(:, {'ID', 'Condition', 'Region'}), 'rows');
averagedThetaPower = accumarray(idx, raw_dataset.AveragedThetaPower, [], @mean);

% Create the final aggregated table
theta_table = table(uniqueKeys.ID, uniqueKeys.Condition, uniqueKeys.Region, averagedThetaPower, num_trials_list, ...
    'VariableNames', {'ID', 'Condition', 'Region', 'AveragedThetaPower', 'NumTrials'});

% Save the table
save(fullfile(output_averaged_spectral_folder, 'theta_table.mat'), 'theta_table');
% Save the table as a CSV file
writetable(theta_table, fullfile(output_averaged_spectral_folder, 'theta_table.csv'));



end