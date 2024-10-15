function [data] = decode_et_v2(ppt_id, subfolder, main_output_path)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function takes in a subfolder expecting there to be eye tracking data present from a tobii eyetracker.
% Attempts to read in the data, and then will preprocess all the tasks which have been preprogrammed below. 
% If there are new tasks or issues then this will stop and wait for a response. This can be overridden.

% ###############################################################################
% This function expects that task engine is installed and the path is accessible.
% ###############################################################################

% If task Engine is unavailable, contact James Ives at:
% james.white1@bbk.ac.uk or
% james.ernest.ives@gmail.com

%% Settings
override_warnings = 0;
outPath = [subfolder, 'decoded_eyetracking\'];

% Tries to open the data with task engine, if it fails it is likely that either task engine can't be accessed
% by matlab or that there is no eye tracking data available. Check manually to see.
data = teData(subfolder);

%% Determine tasks
% These are the tasks that have been preprogrammed in, anything that is not in this list will need to be formatted
% before the below can run.
known_tasks = {'audgate'; 'between'; 'face_ssvep'; 'faceerp'; 'fasterp'; 'gap'; 'ns_50faces'; 'reading'; 'restingvideos2'; 'rocket'};
tasks_to_process = [];
unknown_tasks = [];

% Task list from the ET data
task_list = data.Tasks;

% Loop through the task list and categorize them
for i = 1:length(task_list)
    if ismember(task_list{i}, known_tasks)
        tasks_to_process = [tasks_to_process; task_list(i)];
    else
        unknown_tasks = [unknown_tasks; task_list(i)];
    end
end

% If there are unknown tasks we want the user to know about this
if ~override_warnings && ~isempty(unknown_tasks)
    fprintf('Warning: The following tasks are unknown:\n');
    fprintf('\t- %s\n', unknown_tasks{:});
    input('Press Enter to continue or stop the script and review the tasks...', 's');
    % Crashes the script, select "Pause on errors to see the variables at this stage."
    crashMe
end

% Filter the data using task engine
for task_index = 1:length(tasks_to_process)
    % Grab the task and filter the data log for that task
    task = tasks_to_process{task_index};
    tab = teLogFilter(data.Log, 'task', task, 'topic', 'trial_log_data');
    t_n = height(tab); %n of trials

    % If there are no trials then continue on
    if t_n == 0
        continue
    end

    % I haven't found a gap task with any data in, so if we find one we can use it to populate the stuff below.
    if strcmp(task, 'faceerp')
        disp('Found one!')
    end

    %% Format data
    % Variable names are not consistent between tasks, so this normalises what will be needed below.

    % onset
    % offset
    switch task
        case 'audgate'
            onset_name = 'soundonset';
            offset_name = 'soundoffset';
        case 'fasterp'
            onset_name = 'stim_onset';
            offset_name = 'stim_offset';
        case {'between', 'ns_50faces', 'reading', 'restingvideos2'}
            onset_name = 'movie_onset';
            offset_name = 'movie_offset';
        case {'face_ssvep', 'gap', 'rocket'}
            onset_name = 'presenter_trialonset';
            offset_name = 'presenter_trialoffset';
    end

    % New columns to be added to tab
    tab.nan_percentage = nan(t_n, 1);
    tab.x_mean = nan(t_n, 1);
    tab.y_mean = nan(t_n, 1);
    tab.x_y_sdd = nan(t_n, 1);

    % Raw data epoched by trial
    raw_et = table;
    raw_et.trial_index = linspace(1, t_n, t_n)';
    raw_et.raw_data = cell(t_n, 1);

    %% Loop Through Trials
    for i = 1:t_n
        % Adapted from a script from Teresa del Bianco
        % Select Onset and Offset of trial *1 & 2*
        onset = tab.(onset_name){i};
        offset = tab.(offset_name){i};

        % Get the raw eye tracking data from the eye tracking buffer
        try
            s1 = find(data.ExternalData('eyetracking').Buffer(:, 1) >= onset, 1);
            s2 = find(data.ExternalData('eyetracking').Buffer(:, 1) >= offset, 1);
        catch
            fprintf('No onset or offset for trial %d\n', i)
            continue
        end

        % Filter 1 trial from data and use the ET data from the raw to generate the mean x, mean y, nan proportion and standard deviation
        % of the eye tracking data.
        trial_onscreen = data.ExternalData('eyetracking').Buffer(s1:s2, :);

        % Calculate mid point of eyes in X and Y (column n based on https://docs.google.com/document/d/1-Ddb00b0nVeHw_JD-fexSB7D3x9aAqb4U4GXFyQoRq8/edit?usp=sharing)
        x_gaze = (trial_onscreen(:, 2)+trial_onscreen(:, 17))/2;
        y_gaze = (trial_onscreen(:, 3)+trial_onscreen(:, 18))/2;

        % Calculate Standard Distance Deviation (SDD - root mean squared distance from the centroid)
        x_m = mean(x_gaze, 'omitnan'); %centroid x
        y_m = mean(y_gaze, 'omitnan'); %centroid y

        x_sumSQ = sum((x_gaze-x_m).^2, 'omitnan'); %sum of squares
        y_sumSQ = sum((y_gaze-y_m).^2, 'omitnan');

        x_count = sum(x_gaze, 'omitnan'); %sum of cases (non NANs)
        y_count = sum(y_gaze, 'omitnan');

        x_d = x_sumSQ/x_count;
        y_d = y_sumSQ/y_count;

        sdd = sqrt(x_d+y_d);

        % Calculate the nan percentage of the time on screen
        nanP = round((sum(trial_onscreen(:,4)==0)/length(trial_onscreen))*100, 2); %gaze validity

        % Store in the table
        tab.x_mean(i) = x_m;
        tab.y_mean(i) = y_m;
        tab.x_y_sdd(i) = sdd;

        tab.nan_percentage(i) = nanP;

        % Store raw data
        raw_et.raw_data(i) = {trial_onscreen};
    end

    % Generate the task specific file name
    f_name = ['ET_data_', ppt_id, '_', task];

    %% Save Table to File

    path_list = [{outPath}; {main_output_path}; {fullfile(main_output_path, 'epoched_raw_mat')}; ...
        {fullfile(main_output_path, 'epoched_info_mat')}; {fullfile(main_output_path, 'epoched_info_csv')}; ...
        {fullfile(main_output_path, 'full_teData')}; {fullfile(main_output_path, 'full_data')}];

    for path = 1:length(path_list)
        if ~exist(path_list{path})
            mkdir(path_list{path})
        end
    end
   

    % If the data found is not empty then saves the task level information in .mat and .csv formats
    % both to the individual participant subfolder and to the main results folder
    if size(tab) ~= [0 0]
        % Save the raw data (split trial by trial), as these include variable length cells it doesn't make sense
        % to save these as csvs.
        % Matlab files
        save(fullfile(outPath, [f_name, '_raw_epoched.mat']), 'raw_et');
        save(fullfile(main_output_path, 'epoched_raw_mat', [f_name, '_raw_epoched.mat']), 'raw_et');

        % Save info files on a trial by trial basis
        % Matlab files
        save(fullfile(outPath, [f_name, '_info.mat']), 'tab');
        save(fullfile(main_output_path, 'epoched_info_mat', [f_name, '_info.mat']), 'tab');

        % csv files
        writetable(tab, fullfile(outPath, [f_name, '_info.csv']));
        writetable(tab, fullfile(main_output_path, 'epoched_info_csv', [f_name, '_info.csv']));

        % Store tab and epoched raw data to output file
        % output.(task).et_info = tab;
        % output.(task).raw_epochs = raw_et;
    end
end

%% Format the whole file without task engine
% So that the data can be shared without task engine, this then takes the entire dataset and removes the task engine elements.

save(fullfile(outPath, ['ET_data_', ppt_id, '_teData.mat']), "data");
save(fullfile(main_output_path, 'full_teData', ['ET_data_', ppt_id, '_teData.mat']), "data");

% registeredEvents, ExternalData are task engine classes, we want a version that
% doesn't use these classes in case people in future don't have task engine.
re = struct;
re_fieldnames = fieldnames(data.RegisteredEvents);
for f = 1:length(re_fieldnames)
    re.(re_fieldnames{f}) = data.RegisteredEvents.(re_fieldnames{f});
end
data.RegisteredEvents = re;

ex_items_paths = struct;
f_names = fieldnames(data.ExternalData.Items{1,1}.Paths);
for f = 1:length(f_names)
    ex_items_paths.(f_names{f}) = data.ExternalData.Items{1,1}.Paths.(f_names{f});
end

ex_items = struct;
f_names = fieldnames(data.ExternalData.Items{1,1});
remove = find(strcmp(f_names, 'Paths'));
f_names(remove) = [];
for f = 1:length(f_names)
    ex_items.(f_names{f}) = data.ExternalData.Items{1,1}.(f_names{f});
end
ex_items.Paths = ex_items_paths;

ex_path_path = struct;
f_names = fieldnames(data.ExternalData.Paths{1, 1});
for f = 1:length(f_names)
    ex_path.(f_names{f}) = data.ExternalData.(f_names{f});
end

ex_path = struct;
f_names = fieldnames(data.ExternalData.Paths{1, 1});
for f = 1:length(f_names)
    ex_path.(f_names{f}) = data.ExternalData.(f_names{f});
end

ex = struct;
f_names = fieldnames(data.ExternalData);
remove = find(strcmp(f_names, {'Items'}));
f_names(remove) = [];
remove = find(strcmp(f_names, {'Paths'}));
f_names(remove) = [];
for f = 1:length(f_names)
    ex.(f_names{f}) = data.ExternalData.(f_names{f});
end
ex.Items = ex_items;
ex.Paths = ex_path;

new_data = struct;
full_fieldnames = fieldnames(data);
remove = find(strcmp(full_fieldnames, 'ExternalData'));
full_fieldnames(remove) = [];
for f = 1:length(full_fieldnames)
    new_data.(full_fieldnames{f}) = data.(full_fieldnames{f});
end
new_data.ExternalData = ex;
data = new_data;

save(fullfile(outPath, ['ET_data_', ppt_id, '.mat']), "data");
save(fullfile(main_output_path, 'full_data', ['ET_data_', ppt_id, '.mat']), "data");

end