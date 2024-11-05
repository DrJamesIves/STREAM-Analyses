function filter_et_data_to_table

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date:1st November 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% The purpose of this function is to take in a root folder with et info files in, load
% them in and filter the resulting tab object for the relevant info for analysis.

% Expected that this is the epoched_info_mat folder from ET epoching
root = 'E:\Birkbeck\STREAM INDIA\Datasets\';
root_path = fullfile(root, '2. Preprocessed\2.1.1 Epoched_ET\epoched_info_mat');
out_path = fullfile(root, '3. Analysed\3.3 ET Overview Data\');
out_path_csv = fullfile(out_path, 'csv');
out_path_matlab = fullfile(out_path, 'mat');

% Check if these paths exist and if not make them
out_paths = {out_path; out_path_csv; out_path_matlab};
for o = 1:length(out_paths)
    if ~exist(out_paths{o}, 'dir')
        mkdir(out_paths{o})
    end
end

% Load all files
files = dir(root_path);

% Create an empty table for all participants and trials
eyeT_info_all = table;

% Loop through all the files
for file = 1:length(files)

    % Check that this is a .mat and not a hidden file (starts with .)
    if startsWith(files(file).name, '.') | ~endsWith(files(file).name, '.mat')
        continue
    end

    % Split the filename into its constituent parts.
    name_parts = split(files(file).name, '_');
    ppt_no = name_parts{3};
    attempt = name_parts{4};

    if strcmp(name_parts{5}, 'ns')
        condition = 'ns_50faces';
    elseif strcmp(name_parts{5}, 'face')
        condition = 'face_ssvep';
    else
        condition = name_parts{5};
    end

    % Add continuous suffix if used
    if strcmp(name_parts{6}, 'continuous')
        condition = [condition, '-continuous'];
    end

    % Create an empty table
    eyeT_info = table;

    % Determine which variable name the stim_onset and stim_offset has been given
    switch condition
        case {'audgate', 'audgate-continuous'}
            onset_name = 'soundonset';
            offset_name = 'soundoffset';
        case {'fasterp', 'fasterp-continuous'}
            onset_name = 'stim_onset';
            offset_name = 'stim_offset';
        case {'between', 'ns_50faces', 'reading', 'restingvideos2'}
            onset_name = 'movie_onset';
            offset_name = 'movie_offset';
        case {'face_ssvep', 'gap', 'rocket'}
            onset_name = 'presenter_trialonset';
            offset_name = 'presenter_trialoffset';
    end

    % Load the file, should give a tab file
    load(fullfile(files(file).folder, files(file).name))

    % Filter out the relevant columns
    eyeT_info.id = tab.id;
    eyeT_info.attempt = repmat(attempt, size(tab, 1), 1);
    eyeT_info.timestamp = tab.timestamp;
    eyeT_info.(onset_name) = tab.(onset_name);
    eyeT_info.(offset_name) = tab.(offset_name);
    eyeT_info.nan_percentage = tab.nan_percentage;
    eyeT_info.x_mean = tab.x_mean;
    eyeT_info.y_mean = tab.y_mean;
    eyeT_info.x_y_sdd = tab.x_y_sdd;
    eyeT_info.et_look_proportion = tab.et_looking_prop;

    switch condition
        case 'audgate'
            eyeT_info.trainno = tab.trainno;
        case {'fasterp', 'fasterp-continuous'}
            eyeT_info.trial_type = tab.trialtype;
            eyeT_info.face_orientation = tab.orientation;
            eyeT_info.face_ethnicity = tab.faceethnicity;
        case 'face_ssvep'
            eyeT_info.fs_left = tab.fs_left;
            eyeT_info.fs_right = tab.fs_right;
            eyeT_info.stim_type_left = tab.stim_type_left;
            eyeT_info.stim_type_right = tab.stim_type_right;
        case 'restingvideos2'
            eyeT_info.auditorygating = tab.auditorygatingvideochoice;
            eyeT_info.key = tab.key;

            % As we have to split by condition and save we will do this here rather than in the main loop
            % Matlab file
            eyeT_info2 = eyeT_info;
            eyeT_info = eyeT_info(startsWith(eyeT_info.key, 'FACE'), :);

            if ~isempty(eyeT_info)
                new_fname = [files(file).name(1:end-9), '_social_info'];

                save(fullfile(out_path_matlab, [new_fname, '.mat']), 'eyeT_info');
                writetable(eyeT_info, fullfile(out_path_csv, [new_fname, '.csv']));

                % Add to the "all" info table which can be saved and filtered later
                eyeT_info_all = [eyeT_info_all; ...
                    table(eyeT_info.id(1), ...
                    str2num(attempt), ...
                    {'Rest_Vid_Social'}, ...
                    mean(eyeT_info.nan_percentage, 'omitnan'), ...
                    mean(eyeT_info.x_mean, 'omitnan'), ...
                    mean(eyeT_info.y_mean, 'omitnan'), ...
                    mean(eyeT_info.x_y_sdd, 'omitnan'), ...
                    mean([eyeT_info.et_look_proportion{:}], 'omitnan'))];
            end

            % Now repeat for the non-social
            % Matlab file
            eyeT_info = eyeT_info2;
            eyeT_info = eyeT_info(startsWith(eyeT_info.key, 'TOY'), :);

            if ~isempty(eyeT_info)
                new_fname = [files(file).name(1:end-9), '_non-social_info'];

                save(fullfile(out_path_matlab, [new_fname, '.mat']), 'eyeT_info');
                writetable(eyeT_info, fullfile(out_path_csv, [new_fname, '.csv']));

                % Add to the "all" info table which can be saved and filtered later
                eyeT_info_all = [eyeT_info_all; ...
                    table(eyeT_info.id(1), ...
                    str2num(attempt), ...
                    {'Rest_Vid_Non-Social'}, ...
                    mean(eyeT_info.nan_percentage, 'omitnan'), ...
                    mean(eyeT_info.x_mean, 'omitnan'), ...
                    mean(eyeT_info.y_mean, 'omitnan'), ...
                    mean(eyeT_info.x_y_sdd, 'omitnan'), ...
                    mean([eyeT_info.et_look_proportion{:}], 'omitnan'))];
            end

    end

    if ~strcmp(condition, 'restingvideos2')
        % Matlab files
        save(fullfile(out_path_matlab, files(file).name), 'eyeT_info');
        writetable(eyeT_info, fullfile(out_path_csv, [files(file).name(1:end-4), '.csv']));

        % Add to the "all" info table which can be saved and filtered later
        eyeT_info_all = [eyeT_info_all; ...
            table(tab.id(1), ...
            str2num(attempt), ...
            {condition}, ...
            mean(tab.nan_percentage, 'omitnan'), ...
            mean(tab.x_mean, 'omitnan'), ...
            mean(tab.y_mean, 'omitnan'), ...
            mean(tab.x_y_sdd, 'omitnan'), ...
            mean([tab.et_looking_prop{:}], 'omitnan'))];
    end

end

% Set the headings on the all table
header={'id', 'attempt', 'task', 'nan_percentage','x_mean', 'y_mean', 'x_y_sdd', 'look_proportion'};
eyeT_info_all.Properties.VariableNames = header;

% Save info files on a trial by trial basis
% Matlab files
save(fullfile(out_path, 'overview_et_data.mat'), 'eyeT_info_all');
writetable(eyeT_info_all, fullfile(out_path, 'overview_et_data.csv'));

end
