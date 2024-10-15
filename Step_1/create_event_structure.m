function [eeg_data] = create_event_structure(eeg_data, et_data)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function takes in EEG and eyetracking data, it then uses the unix timestamps in both to
% reconstruct the event structure in the eeg data. This is only for eeg data that is missing
% eeg events. It only works with corresponding eye tracking data, so assumes that appropriate
% checks have been made to ensure that the data corresponds to one another.

% Here we take the registered event names and codes so we can look up the codes later on.
et_labels = et_data.RegisteredEvents.Summary.Label;
et_event_codes = et_data.RegisteredEvents.Summary.eeg;

% Here we find the eeg latency in ms as the eeg start sample (from the .info file) and the timestamp of the 
% first event from et_data which is a syncing event, then minus one from the other, we also generate the 
% latency in eeg samples. It also ensures that these are set to a minimum of 1. 
et_start_time = et_data.Events.timestamp(1);
eeg_latency_ms = round((et_start_time - eeg_data.start_time)*1000);
eeg_latency_samples = round(eeg_latency_ms / (1000 / eeg_data.srate));

% Checks that eeg started after the et (which is how it should be based on matlab scripts)
if eeg_latency_ms > 0
    
    % Warn if latency exceeds 100,000 ms, which could indicate an issue
    if eeg_latency_ms > 100000
        disp(strcat("Latency check over 100,000 (actual value ", num2str(eeg_latency_ms), " for ", et_data.ID, " worth checking on to make sure it's valid"))
    end

    % Set up an empty event structure in the same EEGLAB format
    eeg_events = struct('type', [], 'latency', [], 'latency_ms', []);
    
    % Loop through each event in the ET data to compute corresponding EEG events
    for i = 1:size(et_data.Events, 1)
        % Find the EEG event code by matching ET event label with registered labels
        event_code = et_event_codes(find(strcmp(et_data.Events.data{i}, et_labels)));
        eeg_events(i).type = [num2str(event_code)];
        
        % Calculate latency in milliseconds relative to the first ET event
        eeg_events(i).latency_ms = round(((et_data.Events.timestamp(i) - et_start_time)*1000) + eeg_latency_ms);
        
        % Convert latency to EEG samples based on the sampling rate
        eeg_events(i).latency = round(eeg_events(i).latency_ms / (1000 / eeg_data.srate));
    end
    
    eeg_data.event = eeg_events;

end