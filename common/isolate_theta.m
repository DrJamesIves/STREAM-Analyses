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

root_path = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\2.4 FFT\';
data_path = fullfile(root_path, '2.4.1 Full\');


output_averaged_spectral_folder = fullfile(root_path, 'Overview data\');
output_folder_ET = 'E:\Birkbeck\STREAM\Datasets\1. Raw\ET\';

% To save on clutter make the folders needed manually
folders = {output_averaged_spectral_folder; ...
    output_folder_ET};

for folder = 1:length(folders)
    if ~exist(folders{folder}, 'dir')
        mkdir(folders{folder});
    end
end

%% Setup
ds = getSettings();
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
    condition = strcat(parts{2}, '_', parts{3},  '_', parts{4}, '_', parts{5}); % Extract condition (Face_Onset or Toy_Onset)

    whole_head_theta = mean(squeeze(mean(fft(1,:,thetaIndex), 2)),1);
    frontal_theta = mean(squeeze(mean(fft(1,frontal_indices,thetaIndex), 2)),1);
    central_theta = mean(squeeze(mean(fft(1,central_indices,thetaIndex), 2)),1);
    parietal_theta = mean(squeeze(mean(fft(1,parietal_indices,thetaIndex), 2)),1);
    occipital_theta = mean(squeeze(mean(fft(1,occipital_indices,thetaIndex), 2)),1);

    % Store data in lists
    % ppt_list = [ppt_list; {participantName}];
    if strcmp(condition, 'Rest_Vid_Face_Onset')
        condition = {'Social'};
    elseif strcmp(condition, 'Rest_Vid_Toy_Onset')
        condition = {'Non-social'};
    else
        input('Condition not correct, check and reset', 's')
    end
    % theta_power_list = [theta_power_list; logTheta];

    % Append regional data to the table
    ppt_list = [ppt_list; {participantName}; {participantName}; {participantName}; {participantName}; {participantName}];
    condition_list = [condition_list; condition; condition; condition; condition; condition];
    region_list = [region_list; {'Whole head'; 'Frontal'; 'Central'; 'Parietal'; 'Occipital'}];
    theta_power_list = [theta_power_list; whole_head_theta; frontal_theta; central_theta; parietal_theta; occipital_theta];
end

% Create the table
raw_dataset = table(ppt_list, condition_list, region_list, theta_power_list, ...
    'VariableNames', {'ID', 'Condition', 'Region', 'AveragedThetaPower'});

% Aggregate data by participant, condition, and region
[uniqueKeys, ~, idx] = unique(raw_dataset(:, {'ID', 'Condition', 'Region'}), 'rows');
averagedThetaPower = accumarray(idx, raw_dataset.AveragedThetaPower, [], @mean);

% Create the final aggregated table
theta_table = table(uniqueKeys.ID, uniqueKeys.Condition, uniqueKeys.Region, averagedThetaPower, ...
    'VariableNames', {'ID', 'Condition', 'Region', 'AveragedThetaPower'});

% Save the table
save(fullfile(output_averaged_spectral_folder, 'theta_table.mat'), 'theta_table');
% Save the table as a CSV file
writetable(theta_table, fullfile(output_averaged_spectral_folder, 'theta_table.csv'));



end