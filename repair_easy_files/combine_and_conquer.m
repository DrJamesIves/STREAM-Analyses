function combine_and_conquer

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function concatenates multiple .easy files located in subfolders within a root directory.
% It handles time gaps between files by inserting blank (NaN) data arrays where necessary
% to ensure continuity of the EEG data, preserving the file size and structure.

% root and output paths need to be changed before use.

rootPath = 'D:\Ruihan_data\Raw data\';                      % Path to the directory containing subfolders with .easy files
outputFolder = 'D:\Ruihan_data\0. Concatenated output\';    % Directory for saving concatenated output files

% Uncomment below to automatically create the output folder if it doesn't exist
% if ~exist("outputFolder", 'dir'); mkdir(outputFolder); end

folders = dir(rootPath);

% Loop through each subfolder (skipping the first two entries which are '.' and '..')
for folder = 3:length(folders)

    % Get all files in the current folder, doesn't require all files to be easy files
    easy_files = dir(strcat(folders(folder).folder, '\', folders(folder).name));

    concat_data = [];

    % Loop through each .easy file in the folder
    for easy_file = 1:length(easy_files)
        easy_filename = strcat(easy_files(easy_file).folder, '\', easy_files(easy_file).name);

        % Skip hidden files or non-.easy files
        if startsWith(easy_files(easy_file).name, '.') || ~endsWith(easy_files(easy_file).name, '.easy')
            continue
        end

        new_easy = load_easy_file(easy_filename); % Load the current .easy file data with custom script

        if isempty(concat_data)
            % If this is the first file, just assign its data to the concatenated data
            concat_data = new_easy;
        else
            % Get the last timestamp from the concatenated data and the first timestamp of the new file
            last_timestamp = concat_data(end, end);
            new_start_timestamp = new_easy(1, end);

            % If there's a gap between the last file and the new file, create a NaN array to fill the gap
            % linspace generates equally spaced timestamps for the gap between the two files
            timestamps = linspace(last_timestamp, new_start_timestamp, ((new_start_timestamp - last_timestamp)/2)+1)';
            timestamps = timestamps(2:length(timestamps)-1);

            % Create a blank dataset (filled with NaNs) with the same number of columns as the data (25 columns)
            blank_dataset = nan(length(timestamps), 25);
            blank_dataset(:, 25) = timestamps;

            % Concatenate the existing data, the blank dataset, and the new data
            concat_data = [concat_data; blank_dataset; new_easy];
        end
    end

    % Save the concatenated data to the original folder and the output folder
    save_easy_file(strcat(folders(folder).folder, '\', folders(folder).name, '\', folders(folder).name, '_combined.easy'), concat_data);
    save_easy_file(strcat(outputFolder, folders(folder).name, '_combined.easy'), concat_data);

end