function remergeTrials(ds, merge_tasks, input_folder, delete_files)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 27th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!ds.settings.paths.epochedEEGPath

% The purpose of this function is to merge trials together that have been epoched. This is different to the continuous trials, here there are defined
% onset and offset markers but these trials are spread throughout the task list.

% This function expects ds from ds = getSettings(), which is mainly for paths, a cell array of names that should be found exactly within the file
% names, an input folder to focus on and a boolean value on whether to delete the constituent part files.

% Called only if running this separately
if nargin == 0
    addpath('C:\Users\james\Documents\MATLAB')
    eeglab
    dataset = 'India';
    ds = getSettings(dataset);
    % merge_tasks = ds.settings.trialNamesToMerge;
    % % Get the epoched filenames from the epoched files
    % epoched_files = ds.settings.paths.epochedEEGPath;
    % delete_files = 1;

    merge_tasks = {'Rest_Vid_Face_Onset'; 'Rest_Vid_Toy_Onset'};
    switch dataset
        case 'Malawi'
            input_folder = 'E:\Birkbeck\STREAM\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.1 Full';
        case 'India'
            input_folder = 'E:\Birkbeck\STREAM INDIA\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.1 Full';
    end
    delete_files = 0;
end

% Get the original participant numbers from the raw EEG folder
original_files = dir(strcat(ds.settings.paths.dataPath, '*.mat'));
original_names = cellfun(@(x) x(1:end-4), {original_files.name}, 'UniformOutput', false);
% original_names = vertcat(original_files.name);
% original_names = original_names(:, 1:end-4);

epoched_files = dir(fullfile(input_folder, '*.mat'));

% Loop through the epoched files
for name = 1:length(original_names)
    % Loop through the tasks
    for task = 1:length(merge_tasks)
        % Find the files to merge
        files_to_merge = find(contains({epoched_files.name}, [original_names{name}, '_', merge_tasks{task}, '_']));

        % If no files found to merge then skip to the next task/participant
        if isempty(files_to_merge)
            continue
        end

        % Get the file names, folders and naturally sort them
        merge_fnames = {epoched_files(files_to_merge).name};
        merge_folder = {epoched_files(files_to_merge(1)).folder};
        merge_fnames = natsortfiles(merge_fnames);

        % Load the first file and set it to EEG
        load(fullfile(merge_folder{:}, merge_fnames{1}));
        EEG_full = EEG;

        % Then loop through the others and concatenate them to the end of the EEG file
        for i = 2:length(merge_fnames)
            EEG_to_concat = load(fullfile(merge_folder{:}, merge_fnames{i}));
            EEG_to_concat = EEG_to_concat.EEG;
            EEG_full = pop_mergeset(EEG_full, EEG_to_concat);
        end

        % Delete the constituent file parts
        if delete_files
            for i = 1:length(merge_fnames)
                delete(fullfile(merge_folder{:}, merge_fnames{i}));
            end
        end

        % Save the new concatenated file
        EEG = EEG_full;
        EEG.num_trials_merged = length(files_to_merge);
        ds.dataInfo.num_trials_merged = length(files_to_merge);
        save(fullfile(merge_folder{:}, [merge_fnames{1}(1:end-6), '-merged.mat']), 'EEG', 'ds')

    end
end

end