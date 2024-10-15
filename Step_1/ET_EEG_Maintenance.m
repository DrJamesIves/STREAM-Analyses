function [] = ET_EEG_Maintenance()

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

%% Description

% This script has been co-opted so that it decodes the eye tracking (ET) data and removes it from the task engine ecosystem so it can be natively
% read in Matlab without extra scripts. It also takes in the EEG data, converts it to an EEGLAB format and adds in the events from the ET data.
% This has been added in because there are many EEG datasets without any events.

% To use this script the paths at the start of the script need to be replaced with local file paths.
% If you need a copy of task engine, email James Ives: james.white1@bbk.ac.uk or james.ernest.ives@gmail.com

%% Script
clear variables

% Paths reflect scripts for needed only, it is not advised to genpath into all scripts available.
addpath(genpath('E:\Birkbeck\Scripts\Stream\'));
addpath(genpath('E:\Birkbeck\Scripts\James Common\'));

% Settings
decode_et = 1;         % from uint8 to .mat, epochs ET data
convert_eeg = 1;       % from .easy to eeglab format
transfer_events = 1;   % Transfer events from ET to EEG
save_to_folders = 1;   % Save the files to a folder outside the raw data structure
save_data_log = 1;     % Save data log to show what has been converted/decoded and any notes
dataset = 'Malawi';

% Data structures
notes = {};
data_log = {};

% India/Malawi specifics
switch dataset
    case 'India'
        % Set India specific notes and raw_path
        notes = [notes; {"Intro - this data is from the STREAM India dataset PI Prof Emily Jones Birkbeck University, it has been decoded using Task Engine 2 custom scripts which were originally written by Luke Mason"}];
        root_path = 'E:\Birkbeck\STREAM INDIA\Datasets\';
        raw_path = fullfile(root_path, '0. Untouched\');
    case 'Malawi'
        % Set Malawi specific notes and raw_path
        notes = [notes; {"Intro - this data is from the STREAM Malawi dataset PI Prof Emily Jones Birkbeck University, it has been decoded using Task Engine 2 custom scripts which were originally written by Luke Mason"}];
        root_path = 'E:\Birkbeck\STREAM\Datasets\';
        raw_path = 'E:\Birkbeck\STREAM\data download MW\';
end

% The rest of the notes
notes = [notes; {'data - this is the eye tracking decoded raw data from the original recording. All task engine structures have been removed. No specialised code should be needed to view all of the data.'}];
notes = [notes; {'tab - is created using the teLogFilter script, which in this case has been set to look at just the fasterp trials.'}];
notes = [notes; {'eyeT_info - uses the events from tab to loop through the trials select the correct onset and offsets, calc the mid point of eyes in X and Y column, calculate the standard distance deviation and put all that in a table.'}];

%% Find data
% Set raw_outpath and get subfolders
raw_outpath = strcat(root_path, '1. Raw\');
eeg_output_path = fullfile(root_path, '1. Raw\EEG\');
et_output_path = fullfile(root_path, '1. Raw\ET\');
subfolders = dir(raw_path);

% Loop through the subfolders, Matlab on Windows finds two extra each time. Starting at 3 skips those.
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

    % Useful for recording info used later for stitching etc
    ppt_info = struct;
    ppt_info.ID = subfolder;
    ppt_info.num_subfolders = 0;
    ppt_info.et_found = [];
    ppt_info.decoded_et = {};
    ppt_info.eeg_found = [];
    ppt_info.decoded_eeg = {};
    ppt_info.notes = {};

    subfolder_et_data = [];
    subfolder_eeg_data = [];

    % Loops through the subfolders within a participant folder and finds ET/EEG data
    % decodes this data and makes a note for later.
    for subsubfolderIndex = 3:length(subsubfolders)
        subsubfolder = subsubfolders(subsubfolderIndex).name;
        if startsWith(subsubfolder, '._')
            continue
        end

        % Sets new subsubfolder path and records the number of subfolders
        ppt_info.num_subfolders = ppt_info.num_subfolders + 1;
        subsubfolderFullPath = strcat(subfolderFullPath, subsubfolder, '\');

        %% ET
        if decode_et
            try
                % Decodes the eyetracking while also removing task engine structures
                % Currently this only sorts two types of tasks as the end data has task specific names
                % decode_et_v2(subfolder, main_output_path)
                [et_data] = decode_et_v2(subfolder, subsubfolderFullPath, et_output_path);
                % If successful then records that it was successful, the path and the et_data
                et_found = 1;
                ppt_info.decoded_et = [ppt_info.decoded_et; {strcat(et_data.Path_Session, 'decoded_eyetracking\ET_Data_', et_data.ID, '.mat')}];
                subfolder_et_data = [subfolder_et_data; {et_data}];
            catch
                disp('Failed to decode the eyetracking, this may be because there is no eyetracking data or there is an issue with task engine.')
                et_found = 0;
            end
        else
            et_data = 0;
            et_found = 0;
        end

        % This will be important later to know how many valid files there are for this participant
        ppt_info.et_found = [ppt_info.et_found, et_found];

        %% EEG
        if convert_eeg
            % Get a list of any .easy files within the subfolder.
            easy_files = dir(fullfile(subsubfolderFullPath, '**', '*.easy'));
            if isempty(easy_files)
                eeg_found = 0;
            else
                % Convert .easy format into eeglab, while also pulling in the .info file as a comment/key
                eeg_data = convert_easy_to_eeglab(easy_files);

                % If successful store the eeg_data, path and mark that eeg was found in this subsubfolder
                subfolder_eeg_data = [subfolder_eeg_data; {eeg_data}];
                for i = 1:length(eeg_data)
                    ppt_info.decoded_eeg = [ppt_info.decoded_eeg; {eeg_data{i, 1}.new_filepath}];
                end
                eeg_found = 1;
            end
        else
            eeg_data = 0;
            eeg_found = 0;
        end
        % This will be important later to know how many valid files there are for this participant
        ppt_info.eeg_found = [ppt_info.eeg_found, eeg_found];
    end

    % Checks whether any et or eeg was found within all the subfolders checked for this participant.
    if sum(ppt_info.et_found) > 0; et_found = 1; end;
    if sum(ppt_info.eeg_found) > 0; eeg_found = 1; end;

    %% Check ET/EEG for splits
    if et_found && eeg_found && transfer_events
        % This section of the code finds similar eye tracking and eeg

        %% Check number of files for each, file sizes and timestamps
        % This became redundant later on but it could still be useful for future projects.

        % Get the number of files for each file type
        num_et = length(vertcat(subfolder_et_data{:}));
        num_eeg = length(vertcat(subfolder_eeg_data{:}));

        % Get the amount of data in seconds for the eyetracking datasets, both between the start-end events and overall
        event_et_durations = [];
        full_et_durations = [];
        for i = 1:num_et
            et = subfolder_et_data{i};
            num_events = size(et.Events); num_events = num_events(1);
            et_duration = et.Events.timestamp(num_events) - et.Events.timestamp(1);
            event_et_durations = [event_et_durations; et_duration];
            full_et_durations = [full_et_durations; et.ExternalData.Items.Duration];
        end

        % Get the amount of data in seconds for the EEG datasets
        full_eeg_durations = [];
        for i = 1:num_eeg
            eeg = subfolder_eeg_data{:}{i};
            full_eeg_durations = [full_eeg_durations; eeg.xmax];
        end
        % disp('Here')

        % Get the eyetracking timestamps from the subfolders
        et_timestamps = [];
        for i = 1:length(ppt_info.decoded_et)
            et_timestamp = split(ppt_info.decoded_et{i}, strcat(subfolder, '\')); et_timestamp = et_timestamp{2};
            et_timestamp = split(et_timestamp, '\'); et_timestamp = et_timestamp{1};
            et_timestamp = datetime(et_timestamp, 'InputFormat', 'yyyy-MM-dd''T''HHmmss');
            et_timestamps = [et_timestamps; et_timestamp];
        end

        % Get the EEG timestamps from the EEG file names
        eeg_timestamps = [];
        for i = 1:length(ppt_info.decoded_eeg)
            eeg_timestamp = split(ppt_info.decoded_eeg{i}, '\'); eeg_timestamp = eeg_timestamp(length(eeg_timestamp));
            eeg_timestamp = split(eeg_timestamp{:}, '_'); eeg_timestamp = eeg_timestamp(1);
            eeg_timestamp = datetime(eeg_timestamp, 'InputFormat', 'yyyyMMddHHmmss');
            eeg_timestamps = [eeg_timestamps; eeg_timestamp];
        end

        %% Check data to find a match
        % Checks each of the et_timestamps against the eeg timestamps, if they are less than an hour apart then they count as a match
        for et_t = 1:length(et_timestamps)
            match = find(abs(hours(eeg_timestamps - et_timestamps(et_t))) < 1);

            if isempty(match)
                ppt_info.notes = [ppt_info.notes; {sprintf("No matched EEG datasets for %s", ppt_info.decoded_et{et_t})}];
            else
                % Checks to see if the eeg_data has been matched previously. If so then continues.
                try
                    matched = subfolder_eeg_data{:}{match}.matched;
                catch
                    matched = [];
                    continue
                end

                % If the number of matches is 1 (it only copes with 1) and that eeg_data has not yet been matched then proceed
                if length(match) == 1 & subfolder_eeg_data{:}{match}.matched == 0
                    % Grab the relevant eeg and et data and recreate the eeg event structure
                    et_data = subfolder_et_data{et_t};
                    eeg_data = subfolder_eeg_data{:}{match};
                    eeg_data = create_event_structure(eeg_data, et_data);
                    % If successful (i.e there are events, then proceed to save)
                    if ~isempty(eeg_data.event)
                        eeg_data.matched = 1; EEG = eeg_data;
                        save(strcat(ppt_info.decoded_eeg{match}(1:end-4), '_added_events.mat'), "EEG")
                        f_name = split(ppt_info.decoded_eeg{match}, '\'); f_name = f_name(end);

                        ppt_info.notes = [ppt_info.notes; {sprintf("EEG events added for %s", ppt_info.decoded_eeg{match}(1:end-4))}];
                    else
                        ppt_info.notes = [ppt_info.notes; {sprintf("NO EEG events added for %s", ppt_info.decoded_eeg{match}(1:end-4))}];
                    end
                else
                    input('More than one match available', 's')
                end
            end
        end
    else
        % Save notes depending on whether events are being saved or not
        et_found = sum(ppt_info.et_found);
        eeg_found = sum(ppt_info.eeg_found);
        if transfer_events
            ppt_info.notes = [ppt_info.notes; {sprintf("%d ET valid datasets and %d EEG valid datasets found, not enough to match", et_found, eeg_found)}];
        else
            ppt_info.notes = [ppt_info.notes; {sprintf("%d ET valid datasets and %d EEG valid datasets found, no transfer requested", et_found, eeg_found)}];
        end
    end

    % Save ppt notes
    data_log = [data_log; {ppt_info}];

    % To save hassle can also save data to standalone folders rather than in the individual raw data subfolders
    if save_to_folders
        if eeg_found
            for i = 1:length(ppt_info.decoded_eeg)
                try
                    EEG = eeg_data{1,1};
                catch
                    try
                        EEG = eeg_data;
                    catch
                        input('Issue with EEG struct, may need a better solution')
                        crashMe
                    end
                end
                f_name = split(ppt_info.decoded_eeg{i}, '\'); f_name = f_name(end);
                if ~startsWith(f_name, 'IN') | ~startsWith(f_name, 'MW')
                    f_name = split(f_name, '_'); f_name = f_name(2);
                end
                save(strcat(eeg_output_path, f_name{:}), "EEG")
            end
        end
        % disp('ET saving currently disabled')
        if et_found
            save(strcat(et_output_path, et_data.ID, '.mat'), "et_data")
        end
    end

end

if save_data_log
    % Save the data log of changes made while converting to eeglab and decoding eyetracking
    f_name = strcat(root_path, 'Data_log-', string(datetime('now','Format','y-M-d_HH-mm-ss')), '.mat');
    save(f_name, "data_log")
end

end