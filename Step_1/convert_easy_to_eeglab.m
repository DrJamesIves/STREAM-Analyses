function [eeg_data] = convert_easy_to_eeglab(easy_files, eeg_output_path)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function takes in a subfolder within participant easy_files, it then searches the folders for .easy enobio files. It makes sure that the files
% are at least a megabyte and then loads them into an EEGLAB format and saves them as .mat files. This has the advantage of making the files a lot
% smaller without losing any data and also makes them quicker and easier to read into matlab. Also saves a quick output of what was converted in the
% root folder.

eeg_data = [];

% Settings for reading the .info file later
formatSpec = '%s';
delimiter = '\n';

%% Find info files
% Search through all subfolders where easy files were found to find the corresponding info files.
% Depreciated as they aren't strictly necessary and they get automatically loaded during pop_easy

% info_files = {};
% 
% % Get all unique folders/subfolders with easy files
% unique_folders = unique({easy_files.folder});
% 
% % Loop through each unique folder and search for .info files
% for i = 1:length(unique_folders)
%     % Find all .info files in this folder
%     info_files_in_folder = dir(fullfile(folder, '*.info'));
% 
%     % Append the results to the file_struct
%     info_files = vertcat(info_files{:}, info_files_in_folder);
% end

%% Load and convert easy to eeglab
% We are doing a standard conversion from easy to eeglab, but we are adding in the info file information as an extra key in the struct and in the
% comments. This is being put into the comments because eeglab operations may strip out non-standard keys.

% Loop through each subfolder
for i = 1:length(easy_files)
    filename = fullfile(easy_files(i).folder, easy_files(i).name);

    % Check if file name starts with "._", which is a temp file, and if its size is greater than 100KB
    if ~startsWith(easy_files(i).name, '._') && easy_files(i).bytes > 100000
        try
        % Extract the file name (without extension)
        [~, name, ~] = fileparts(filename);

        % Load the .easy file using pop_easy
        EEG = pop_easy(filename, 1);

        % Load the .easy file as a raw file to get the start time as a unix timestamp
        raw_EEG = load(filename);
        EEG.start_time = raw_EEG(1, end)/1000;

        % Used to match with an et file
        EEG.matched = 0;

        % Attempt to load in the corresponding info file
        try
            % Read in the info file
            info_filepath = strcat(easy_files(i).folder, '\', name, '.info');
            fileID = fopen(info_filepath);
         
            % Grab the info from the .info file, if fileID is -1 then no info file was found
            if fileID > -1
                info = textscan(fileID,formatSpec,'Delimiter',delimiter);
                fclose(fileID);
                info = info{1};

                % Redundancy in case the .info gets striped by eeglab functions
                EEG.comments = vertcat(EEG.comments, 'Information pulled from .info file', info);
                EEG.info = info;
            else
                EEG.comments = vertcat(EEG.comments, 'No .info file');
                EEG.info = 0;
            end
 
        catch
            fprintf('\nNo info file for %s\n', name)
        end

        % Check whether other eeg data has been saved for this participant and if so increase the counter
        file_counter = 1;
        temp_name = split(name, '_'); temp_name = temp_name{2};
        for j = 1:100
            if exist(fullfile(easy_files(i).name, [name, '_', num2str(j), '.mat']), 'file')
                file_counter = file_counter + 1;
            else 
                break
            end
        end

        raw_name = [temp_name, '_', num2str(file_counter)];
        name = [name, '_', num2str(file_counter)];

        % Store the new filename and location
        EEG.new_filepath = fullfile(easy_files(i).folder, [name '.mat']);

        % Store the EEG data in the structure that is returned
        eeg_data = [eeg_data; EEG];

        % Save data in .mat format in the same folder
        save(fullfile(easy_files(i).folder, [name '.mat']), 'EEG');

        % Print file name, number of points, and number of events
        fprintf('File: %s\tPoints: %d\tEvents: %d\n', name, EEG.pnts, length(EEG.event));

        % Place the output in the output struct for later.
        output{i, 1} = name; output{i, 2} = EEG.pnts; output{i,3} = length(EEG.event);

        catch ME

        % If an error occurs during loading, skip to the next file
        fprintf('Error loading file: %s. Skipping...\n', easy_files(i).name);
        fprintf('Error message: %s\n', ME.message);
        continue;

        end
    else
        % Ignore files that don't meet the criteria
        fprintf('Ignoring file: %s\n', easy_files(i).name);
    end
end

end
