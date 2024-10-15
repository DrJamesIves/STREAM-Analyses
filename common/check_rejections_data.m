function [rejection_table] = check_rejections_data(folder)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function can be used after preprocessing to check the rejected data folder and generate stats
% on why different data were rejected using the data_info struct.

% Whatever your folder of interest is should be specified otherwise defaults.
if nargin == 0
    folder = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.3 Rejected data';
end

files = dir(fullfile(folder, '*.mat'));

% Settings
min_data_length = 10;       % Seconds

% Store sample wide stats
rejection_info = struct;
[rejection_info.name, rejection_info.chansRemovedDuringRobustAvg, rejection_info.chansInterpolatedForNoise, rejection_info.rejectedSegments] = deal([]);

for file = 1:length(files)

    % Load the file
    load(fullfile(files(file).folder, files(file).name));

    % Grab the name without extention for later
    f_name = files(file).name(1:end-4);

    Fs = EEG.srate;

    % Data are rejected with less than the minimum data are shown here.
    if size(EEG.data, 2) < Fs * min_data_length
        rejection_info.name = [rejection_info.name; {f_name}];
        rejection_info.chansRemovedDuringRobustAvg = [rejection_info.chansRemovedDuringRobustAvg; NaN];
        rejection_info.chansInterpolatedForNoise = [rejection_info.chansInterpolatedForNoise; NaN];
        rejection_info.rejectedSegments = [rejection_info.rejectedSegments; NaN];
        continue
    end

    % Sum the number of channels removed during robust average, number of interpolated noisy channels and number of rejected sections as a percentage
    rejection_info.name = [rejection_info.name; {f_name}];
    rejection_info.chansRemovedDuringRobustAvg = [rejection_info.chansRemovedDuringRobustAvg; sum(dataInfo.chansRemovedFromRobustAvg)];
    rejection_info.chansInterpolatedForNoise = [rejection_info.chansInterpolatedForNoise; sum(dataInfo.channelsRejectedForNoise)];
    if ~isempty(dataInfo.rejectedSections)
        rejection_info.rejectedSegments = [rejection_info.rejectedSegments; sum(dataInfo.rejectedSections{1,1})/length(dataInfo.rejectedSections{1,1})];
    else
        rejection_info.rejectedSegments = [rejection_info.rejectedSegments; NaN];
    end

    clear dataInfo EEG

end

table_headers = {'Name'; 'chansRemovedDuringRobustAvg'; 'chansInterpolatedForNoise'; 'rejectedSegmentsPercentage'};
rejection_table = table(rejection_info.name, rejection_info.chansRemovedDuringRobustAvg, rejection_info.chansInterpolatedForNoise, ...
    rejection_info.rejectedSegments);
rejection_table.Properties.VariableNames = table_headers;

%% Check rejections folders for small files specifically
% in_rej = 'E:\Birkbeck\STREAM INDIA\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.3 Rejected data';
% mw_rej = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.3 Rejected data';
% 
% rej_folders = {in_rej; mw_rej};
% rej_size = struct;
% rej_size.in = [];
% rej_size.mw = [];
% labels = {'in', 'mw'};
% 
% for folder = 1:length(rej_folders)
%     files = dir(rej_folders{folder});
% 
%     for file = 1:length(files)
%         eeg_file = strcat(files(file).folder, '\', files(file).name);
%         if ~endsWith(eeg_file, '.mat')
%             continue
%         end
% 
%         % Temporary measure to only focus on one trial type
%         % if ~contains(files(file).name, 'Rest_Vid')
%         %     continue
%         % end
% 
%         disp(fprintf('Processing ...%s', eeg_file))
% 
%         load(eeg_file);
% 
%         Fs = EEG.srate;
% 
%         if size(EEG.data, 2) < Fs * 10
%             rej_size.(labels{folder}) = [rej_size.(labels{folder}); EEG.pnts/Fs];
%             continue
%         end
%     end
% end

end