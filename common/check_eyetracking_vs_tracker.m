function check_eyetracking_vs_tracker

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% After running ET_EEG_Maintenance and preprocessing scripts, this scrip can be used to check
% the number of valid ET/EEG/tracker files.

raw_path = 'E:\Birkbeck\STREAM\data download MW';

% Paths reflect scripts for needed only, it is not advised to genpath into all scripts available.
addpath(genpath('E:\Birkbeck\Scripts\Stream\'));
addpath(genpath('E:\Birkbeck\Scripts\James Common\'));

subfolders = dir(raw_path);
subfolder_tracker = [];

for subfolderIndex = 3:length(subfolders)
    % Sets the subfolder name and if it starts with "._" it is a temporary file and ignored
    subfolder = subfolders(subfolderIndex).name;
    if ~startsWith(subfolder, 'MW')
        if ~startsWith(subfolder, 'IN')
            continue;
        end
    end

     % Sets the subfolder path and searches for subsubfolders
    fprintf('Processing: %s\n', subfolder)
    subfolderFullPath = strcat(subfolders(subfolderIndex).folder, '\', subfolders(subfolderIndex).name, '\');
    subsubfolders = dir(subfolderFullPath);

    subfolder_temp = [];

    % Loops through the subfolders within a participant folder and finds ET/EEG data
    % decodes this data and makes a note for later.
    for subsubfolderIndex = 3:length(subsubfolders)
        subsubfolder = subsubfolders(subsubfolderIndex).name;
        if startsWith(subsubfolder, '._')
            continue
        end

        decoded_et = 0;
        eeg = 0;
        tracker = 0;

        % Sets new subsubfolder path and records the number of subfolders
        subsubfolderFullPath = strcat(subfolderFullPath, subsubfolder, '\');

        if exist(fullfile(subsubfolderFullPath, 'decoded_eyetracking'))
            decoded_et = 1;
        end

        easy_files = dir(fullfile(subsubfolderFullPath, '**', '*.easy'));
        if length(easy_files) > 0
            eeg = 1;
        end

        tracker_files = dir(fullfile(subsubfolderFullPath, '**', 'tracker*.mat'));
        if length(tracker_files) > 0
            tracker = 1;
        end

        subfolder_temp = [subfolder_temp; {[subfolder, '\', subsubfolder]} decoded_et, eeg, tracker];
    end
    subfolder_tracker = [subfolder_tracker; {subfolder_temp}];
end

% Missing ET only
for subfolder = 1:length(subfolder_tracker)
    subfolder_to_test = vertcat(subfolder_tracker{subfolder});
    for subsubfolder = 1:size(subfolder_to_test, 1)
        if subfolder_to_test(subsubfolder, 2) == 0 & subfolder_to_test(subsubfolder, 3) == 1 & subfolder_to_test(subsubfolder, 4) == 1
            fprintf('EEG found without decoded ET but with a tracker, subfolder %s, subsubfolder %s.\n', num2str(subfolder), num2str(subsubfolder))
        end
           
    end
end

% Missing ET and tracker
for subfolder = 1:length(subfolder_tracker)
    subfolder_to_test = vertcat(subfolder_tracker{subfolder});
    for subsubfolder = 1:size(subfolder_to_test, 1)
        if subfolder_to_test(subsubfolder, 2) == 0 & subfolder_to_test(subsubfolder, 3) == 1 & subfolder_to_test(subsubfolder, 4) == 0
            fprintf('EEG found without decoded ET or a tracker, subfolder %s, subsubfolder %s.\n', num2str(subfolder), num2str(subsubfolder))
        end
           
    end
end

% Missing EEG only
count = 0;
file_recon = {};
for subfolder = 1:length(subfolder_tracker)
    subfolder_to_test = vertcat(subfolder_tracker{subfolder});
    for subsubfolder = 1:size(subfolder_to_test, 1)
        if subfolder_to_test{subsubfolder, 2} == 1 & subfolder_to_test{subsubfolder, 3} == 0 & subfolder_to_test{subsubfolder, 4} == 1
            % fprintf('EEG found without decoded ET or a tracker, subfolder %s, subsubfolder %s.\n', num2str(subfolder), num2str(subsubfolder))
            file_recon = [file_recon; {subfolder_to_test{subsubfolder, 1}}];
            count = count + 1;
        end
    end
end




end