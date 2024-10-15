function concat_EEG(filenames)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% The purpose of this function is to concatenate EEG files. The idea is that this is for files where the recording has been stopped and then
% restarted, guidance of things to think about is included below. 

%% Expected input
% eeg_filenames, expecting a vertically concatenated cell array of full filenames including extensions. E.g. [{eeg_a_1; eeg_a_2}, {eeg_b_1; eeg_b_2}]

%% Guidance
% The concatenation below concatenates EEGLAB formatted EEG using an EEGLAB process. EEGLAB will need to be installed. If non-EEGLAB formatted EEG are
% used it will fail.
%
% The concatenation below concatenates without taking into account timestamps etc. So, if the data stream was paused/stopped for a time and then
% restarted and you would like to respect the relative timestamps then you will need to create a blank dataset with data/times/pnts that fills the gap
% between the two. If you don't care about the time that was paused then this is fine. There will likely be edge effects but n-1 for the number of
% filenames provided.

%% Script
% Loop through the filenames found and load each internal cell within the array
for filename_set = 1:length(filenames)
    concatEEG = [];
    % Loop through and load each file within the filename set to be concatenated
    for filename = 1:length(filename_set)
        fprintf(strcat('Loading\t', filename_set{filename}, '\t\tfor concatenating EEG\n'))
        % For each of the files to be concatenated load the file and merge it into one dataset.
        EEG = load(filename, 'EEG'); EEG = EEG.EEG;
        if isempty(concatEEG)
            concatEEG = EEG;
        else
            concatEEG = pop_mergeset(concatEEG, EEG);
        end
    end
    % Save the concatenated dataset.
    EEG = concatEEG;
    save(strcat(filename(1:end-4), '_concatenated.mat'), 'EEG')
end

end