function epochEEG(ds)

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

% This function loads data from the raw_EEG data folder, finds the onset and offset event codes as defined in the ds.settings struct, epochs the data and
% saves the data in the epoched_EEG folder. The script epochs the EEG based on each trial type an saves for each type, I.e. if you have social and
% non-social then the EEG file becomes two separate EEG files, labelled accordingly. The script can handle events where the recording was started
% after the start of the session or the recording was stopped before the end, providing these are within the min-max trial length.

% The script lets you know if it finds trials that are too long or too short. These likely represent trials that were cut off and there wasn't enough
% data, trials that were skipped (all the trial events are put in when skipping), trials where the wrong event has been sent at the start/end of a
% trial. In future I hope to be able to save these last trials but there are a lot of moving parts.

% The script will skip files where it finds that it has already epoched the trials. It does this by counting he number of epoched trials and comparing
% it against the expected number of trials total. Unfortunately, this means that if there are less than the expected number of trials in the recording
% then it will never skip that file. Unfortunately, this is likely to be a lot of trials.

%% Settings and init

if nargin == 0
    addpath('E:\Birkbeck\Scripts\Stream\Theta\');
    addpath('E:\Birkbeck\Scripts\James Common\');
    ds = getSettings;
end
% Finds all the raw EEG files
ds.dataInfo.files = dir(strcat(ds.settings.paths.dataPath, '*.mat'));

% Event prefix
if ds.settings.epochWithDINMarker
    prefix = 'DIN';
else
    prefix = '';
end

%% Endless events
% Are there events without a pre defined start and end point? I.e a bunch of events that just start, repeat and stop so it's unclear which event
% signifies the start and which the end? If so they have to be listed here along with events to ignore.
endless_events = 1;     % Whether to include these types of events
% Separate list for each trial type, within cells so they can be different lengths
endless_event_codes = {{'400', '401', '402'}'; ...
    {'310', '311', '312', '313', '314', '315', '316', '317', '318', '320', '321', '322', '330', '331', '340', '341', '351', '352'}'};
% Names for the above trials in the same order (these will be used for file names).
endless_event_names = {'Aud_gate'; 'fast_erp'};
% Event codes for syncing, skipping , attention getting, pausing and fixations among others.

event_codes_to_ignore = {'247', '248', '249', '250', '251', '252', '253', '254'}';

% These codes need neither be acknowledged or ignored
% 11XX codes are for the rocket task, these have a defined onset and offset
% 62-72 codes are for the gap task, these have a defined onset and offset
    % '1103', '1104', '1105', '1106', '1107', '1108', '1109', '1110', '1111', '1112', '1113', '1114', '1115', ...
    % '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72';

%% Find raw files
% Loops through all the raw EEG files
for file = 1:length(ds.dataInfo.files)

    % First check whether this epoching has already been done by checking whether the first of each of the trial types has been epoched. We will assume
    % that it has been done successfully if all trial types are present. If you would like to redo a file then delete the appropriate trial file.
    numCompleted = 0;
    for trialType = 1:length(ds.settings.eventNames)
        if exist(strcat(ds.settings.paths.epochedEEGPath, ds.dataInfo.files(file).name(1:end-4), '_', ds.settings.eventNames{trialType}, '_1.mat'))
            numCompleted = numCompleted + 1;
        end
    end

    % If there are one or more trials that are missing then attempt to epoch the data
    if numCompleted ~= length(ds.settings.eventNames)

        fprintf(strcat('Loading\t\t', ds.dataInfo.files(file).name, '\t\tfor epoching\n'))
        % Loads EEG file
        load(strcat(ds.dataInfo.files(file).folder, '\', ds.dataInfo.files(file).name));

        if isempty(EEG.event)
            continue
        end

        % We'll create a copy so we can use that later
        EEG2 = EEG;

        % This counter is used to count the number of trial types that have no trials. Later on, if it determines that there are no trials at all it moves
        % this file to a "No trials found" folder.
        trialTypeCounter = 0;

        %% Checks
        % Checks that the EEG sampling rate is the same as expected (specified in the master script)
        if EEG.srate ~= ds.settings.eegPreproc.expectedEEGSampleRate
            fprintf(strcat('Sampling rate was different from expected, resampling to correct rate'))
            if EEG.srate > ds.settings.eegPreproc.expectedEEGSampleRate
                EEG.data = resample(EEG.data, EEG.srate, ds.settings.eegPreproc.expectedEEGSampleRate);
            elseif EEG.srate < ds.settings.eegPreproc.expectedEEGSampleRate
                ME = MException('MATLAB:samplingError', 'Sampling rate is lower than expected, you could inpterpolate this data with interp1, but this script wont do that automatically.');
                throw(ME)
            end
        end

        if EEG.nbchan ~= ds.settings.eegPreproc.expectednbChannels
            ME = MException('MATLAB:channelError', 'Number of channels is not as expected.');
            throw(ME)
        end

        %% Search and epoching
        % Searches for each trial type individually
        for trialType = 1:length(ds.settings.eventNames)
            % Creates an onset/offset even space for each of the trial types. This will be filled later.
            ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType}) = [];
            ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType}) = [];

            % Finds all the onset and offset events for each trial type and stores in the dynamic struct.
            for i = 1:length(EEG2.event)
                event = replace(EEG2.event(i).type, ' ', '');
                if strcmp(event, strcat(prefix, num2str(ds.settings.onOffsetEventNumbers(trialType, 1)))) | strcmp(EEG2.event(i).type, strcat('D', num2str(ds.settings.onOffsetEventNumbers(trialType, 1))))
                    % Finds the onset/offset indices and latencies and storet them in dataInfo. Note that the mff file has incredibly precise latencies that aren't super
                    % useful to us when they are more specific than the sampling rate of the EEG setup. So when the latencies are taken they are rounded.
                    ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType}) = [ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType}); i, 0];
                    ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType}) = [ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType}); round(EEG2.event(i).latency), 0, 0];
                elseif strcmp(event, strcat(prefix, num2str(ds.settings.onOffsetEventNumbers(trialType, 2)))) | strcmp(EEG2.event(i).type, strcat('D', num2str(ds.settings.onOffsetEventNumbers(trialType, 2))))
                    % Checks to see if the latencies for this trial type are empty, which indicates that there wasn't a starting event found so it isn't a full trial.
                    % It also checks that there isn't already another value already in the end spot, this indicates that the end trigger has been sent twice. We want to
                    % ignore the second one.
                    if ~isempty(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})) & ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(end, 2) == 0
                        ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(end, 2) = i;
                        ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,2) = round(EEG2.event(i).latency);
                        % Gives us the duration, which will be used later
                        ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,3) = ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,2) - ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,1);
                    else
                        % If there is an end event before any start events it is likely that this means the recording was started after the screen task started, so this
                        % checks whether this event is close to the start of the recording and if it is and it meets the minimum trial length then the first data point is
                        % taken as the starting event. If it doesn't then there's not much point taking it anyway.
                        % if EEG2.event(i).latency < ds.settings.expectedTrialLength * EEG.srate & EEG2.event(i).latency > ds.settings.minTrialLength * EEG.srate
                        %     ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(end, 2) = 0;
                        %     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,2) = 1;
                        %     % Gives us the duration, which will be used later
                        %     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,3) = ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,2) - ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(end,1);
                        % end
                    end
                end
            end

            % Looks for orphaned trials, ones with a start and not an end. This can occur if the wrong end trigger is sent. So this takes those, and tries to find
            % any end trigger. It then checks to see if this corresponds to within 10% of the expected trial length. If so it assumes that this is correct and
            % creates a new trial
            % missedEndEvents = find([ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(:, 2)] == 0);
            % % If missing end events are found
            % if missedEndEvents > 0
            %     % Cycle through the end events
            %     for i = 1:length(missedEndEvents)
            %         % Grab the corresponding start event and latency
            %         missedOnsetEvent = ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(missedEndEvents(i), 1);
            %         missedOnsetLatency = ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i), 1);
            %
            %         % Get a list of codes from the onset event until the end, find the next end event and see if that matches
            %         % eventCodes = string(vertcat(EEG2.event(missedOnsetEvent:end).code));
            %         % endEventIndex = find(strcmp(string(strcat('DIN', num2str(ds.settings.onOffsetEventNumbers(:, 2)))), eventCodes), 1);
            %
            %         % Cycle through the events starting from the missed start event
            %         for j = missedOnsetEvent:length(EEG2.event)
            %             % This takes the event type readout and compares it to each of the end event and returns true if it finds any end event.
            %             if ~isempty(find(strcmp(EEG2.event(j).type, string(strcat('DIN', num2str(ds.settings.onOffsetEventNumbers(:, 2))))))) | ~isempty(find(strcmp(EEG2.event(j).type, string(strcat('D', num2str(ds.settings.onOffsetEventNumbers(:, 2)))))))
            %                 % This checks if that end event is within 10% either side of the expected trial length and if it is it takes  this event as the end event
            %                 if (missedOnsetLatency + (ds.settings.expectedTrialLength * EEG.srate * 0.9)) < round(EEG2.event(j).latency) < (missedOnsetLatency + (ds.settings.expectedTrialLength * EEG.srate * 1.1))
            %                     ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(missedEndEvents(i), 2) = j;
            %                     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),2) = round(EEG2.event(j).latency);
            %                     % Gives us the duration, which will be used later
            %                     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),3) = ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),2) - ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),1);
            %                 end
            %                 % If we get to the end of the data and there hasn't been an end event this suggests that the recording was stopped before the end event could be sent.
            %                 % So in this case take the last datapoint as the the end event and see if this is long enough
            %             elseif j == length(EEG2.event)
            %                 if (missedOnsetLatency + (ds.settings.expectedTrialLength * EEG.srate * 0.9)) < round(EEG2.pnts/EEG2.srate) < (missedOnsetLatency + (ds.settings.expectedTrialLength * EEG.srate * 1.1))
            %                     ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(missedEndEvents(i), 2) = j;
            %                     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),2) = round(EEG2.event(j).latency);
            %                     % Gives us the duration, which will be used later
            %                     ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),3) = ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),2) - ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(missedEndEvents(i),1);
            %                 end
            %             end
            %         end
            %     end
            % end

            % Any trials that are shorter than the minimum trial length are excluded in this loop
            if ds.settings.checkTrialLength
                for i = size(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType}), 1):-1:1
                    if ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(i,3) < ds.settings.minTrialLength * EEG.srate
                        fprintf(strcat('Warning: a trial in\t', ds.dataInfo.files(file).name, '\t', ds.settings.eventNames{trialType}, ' is under the minimum epoching length\n'))
                        ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(i,:) = [];
                        ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(i,:) = [];
                    elseif ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(i, 3) > ds.settings.maxTrialLength * EEG.srate
                        fprintf(strcat('Warning: a trial in\t', ds.dataInfo.files(file).name, '\t', ds.settings.eventNames{trialType}, ' is over the maximum epoching length\n'))
                        ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(i,:) = [];
                        ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(i,:) = [];
                    end
                end
            end

            % Does a check to make sure that there are some trials found, if not then it skips the next section.
            if ~isempty(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType}))

                for trial = 1:size(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType}), 1)
                    % We're going to reset the EEG variable so we can fill it with epochs to save
                    EEG.setname = ds.settings.eventNames{trialType};
                    [EEG.data, EEG.times, EEG.pnts, EEG.event] = deal([]);
                    EEG.trials = 0;
                    % minDur = round(min(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(:,3))); % Gives the minimum duration
                    % maxDur = round(max(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(:,3))); % Gives the maximum duration
                    %
                    % % This checks the range of durations found, and if there is more than a 1 second discrepency then it flags it to the user
                    % if maxDur-minDur > ds.settings.maxEventDiscrepency * EEG.srate;
                    %     resp = input(fprintf(strcat('Potential error - trial duration range is\t', num2str(maxDur-minDur), ' seconds which is above the expected limit. Hit enter to confirm this is okay or i to crash the script and inspect (if "Pause with errors" is on) -\t')), 's');
                    %     if resp == 'i'
                    %         hh
                    %     end
                    % % elseif minDur < ds.settings.minTrialLength * EEG.srate
                    % %     resp = input(fprintf(strcat('Potential error - trial duration is less than\t', num2str(ds.settings.minTrialLength), ' seconds. Hit enter to confirm this is okay or i to crash the script and inspect (if "Pause with errors" is on) -\t')), 's');
                    % %     if resp == 'i'
                    % %         hh
                    % %     end
                    % end

                    % Next now that we have found the indices for these events we'll go through and do the epoching.
                    % for i = 1:size(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType}), 1)
                    % We don't need to do this but it makes it much easier to read. We first set the epoch onset and offset, which we'll use to chunk the data later.
                    % The onset, offset, minDur and maxDur are rounded because the latencies are incredibly precise but we need these to be in samples so they need to be
                    % whole numbers.
                    onset = round(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(trial, 1));
                    offset = round(ds.dataInfo.onOffsetLatencies.(ds.settings.eventNames{trialType})(trial, 2));

                    if onset == 0 || offset == 0
                        continue
                    end

                    try
                        EEG.data = EEG2.data(:, onset:offset);
                    catch
                        continue
                    end
                    EEG.times = EEG2.times(1, onset:offset);

                    % With the events we want the information to be relative to the new epochs rather than to the whole original file, this section takes in the original
                    % events, stores the first row of information and removes that from every subsequent row.

                    % For the first event in the row, store its values
                    first_latency = round(EEG2.event(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(trial, 1)).latency);
                    first_init_index = first_latency;% EEG2.event(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(trial, 1)).description;
                    % first_init_time = EEG2.event(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(i, 1)).init_time;

                    % Create a temp variable we can edit
                    temp = EEG2.event(ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(trial, 1):ds.dataInfo.onOffsetEvents.(ds.settings.eventNames{trialType})(trial, 2));
                    first_init_time = temp(1).latency_ms;

                    % Modifies all subsequent rows
                    for j = 1:length(temp)
                        temp(j).latency = round(temp(j).latency) - first_latency + 1;
                        % temp(j).description = num2str(str2num(temp(j).description) - str2num(first_init_index) + 1);
                        temp(j).latency_ms = temp(j).latency_ms - first_init_time + 1;
                    end

                    EEG.event = temp;
                    clear temp
                    % end

                    EEG.trials = 1;
                    EEG.pnts = offset - onset + 1;
                    EEG.xmax = max(EEG.data, [], 'all');
                    EEG.xmin = min(EEG.data, [], 'all');

                    % Example of how you could epoch the data if you want to instead epoch specific times relative to events
                    % EEG3 = pop_epoch(EEG, { 'DIN241' 'DIN242' }, [-0.1 0.5], 'newname', 'Social4Hz', 'epochinfo', 'yes');

                    % Save the data
                    save(strcat(ds.settings.paths.epochedEEGPath, ds.dataInfo.files(file).name(1:end-4), '_', ds.settings.eventNames{trialType}, '_', num2str(trial), '.mat'), 'ds', 'EEG');
                end
            else
                % If there are no trials print this to the command window
                fprintf(strcat('No trials found for\t', ds.dataInfo.files(file).name, '\t-\t', ds.settings.eventNames{trialType}, '\n'))

                % This checks how many trial types have no trials, if it is every trial type being searched for then this file gets moved to the "No trials found"
                % folder, this is so that in future when this pipeline is run it skips this file.
                trialTypeCounter = trialTypeCounter + 1;
                % if trialTypeCounter == length(ds.settings.eventNames)
                %     movefile(strcat(ds.settings.paths.rawEEGPath, ds.dataInfo.files(file).name), strcat(ds.settings.paths.rawEEGPath, 'No trials found\', ds.dataInfo.files(file).name))
                % end
            end
        end

        % Clearing this info so we don't epoch based on previous files
        ds.dataInfo.onOffsetEvents = [];
        ds.dataInfo.onOffsetLatencies = [];

        if endless_events
            for index = 1:length(endless_event_codes)
                % Find where the events within the trial are.
                % Assuming binary_timeseries is your input timeseries
                % where 0 = exclude, 1 = include, 2 = flexible
                binary_timeseries = ismember(vertcat({EEG2.event(:).type}), endless_event_codes{index});
                ignore_timeseries = ismember(vertcat({EEG2.event(:).type}), event_codes_to_ignore);
                % Make the ignore events 2s so they are flexible
                ignore_timeseries = ignore_timeseries * 2;
                binary_timeseries = binary_timeseries + ignore_timeseries;

                % Identify all the segments of interest (1s and 2s)
                binary_timeseries_filtered = binary_timeseries;
                binary_timeseries_filtered(binary_timeseries_filtered == 2) = 1; % Treat 2 as 1 temporarily

                % Find the start and end indices considering the merged 1s and 2s
                start_indices = find(diff([0, binary_timeseries_filtered]) == 1);
                end_indices = find(diff([binary_timeseries_filtered, 0]) == -1);

                % Reevaluate start_indices and end_indices based on original binary_timeseries
                final_start_indices = [];
                final_end_indices = [];

                for i = 1:length(start_indices)
                    % Extract the segment from start to end
                    segment = binary_timeseries(start_indices(i):end_indices(i));

                    % Check if the segment contains at least one 1
                    if any(segment == 1)
                        % If the segment contains only 1s and 2s without interruption, treat as one block
                        if all(ismember(segment, [1 2]))
                            final_start_indices = [final_start_indices; start_indices(i)];
                            final_end_indices = [final_end_indices; end_indices(i)];
                        else
                            % If 2s are surrounded by 0s, consider them separately
                            % Find the actual start of 1s within this block
                            one_start = find(segment == 1, 1, 'first') + start_indices(i) - 1;
                            one_end = find(segment == 1, 1, 'last') + start_indices(i) - 1;

                            % Add this 1-sequence as a new block
                            final_start_indices = [final_start_indices; one_start];
                            final_end_indices = [final_end_indices; one_end];
                        end
                    end
                end

                onoffsets = [final_start_indices, final_end_indices];

                % % Then find the onsets and offsets of these events and class these as the trial on/offsets
                % % Initialize start_indices and end_indices with empty arrays
                % start_indices = []; end_indices = [];
                %
                % % Check if the first element is 1
                % if binary_timeseries(1) == 1; start_indices = 1; end
                % % Check if the last element is 1
                % % if binary_timeseries(end) == 1; end_indices = length(binary_timeseries); end
                %
                % % Find the start indices (transitions from 0 to 1)
                % start_indices = [start_indices; find(diff([0, binary_timeseries(2:end)]) == 1)];
                % % Need to add 1 to the start to give the correct indices, and recorrect the first one if it was 1
                % start_indices = start_indices + 1;
                % if start_indices(1) == 2; start_indices(1) = 1; end
                %
                % % Find the end indices (transitions from 1 to 0)
                % end_indices = [find(diff([binary_timeseries, 0]) == -1); end_indices];
                %
                % onoffsets = [start_indices', end_indices'];

                % Does a check to make sure that there are some trials found, if not then it skips the next section.
                if ~isempty(onoffsets)
                    for trial = 1:size(onoffsets, 1)
                        % We're going to reset the EEG variable so we can fill it with epochs to save
                        EEG.setname = endless_event_names{index};
                        [EEG.data, EEG.times, EEG.pnts, EEG.event] = deal([]);
                        EEG.trials = 0;

                        % Next now that we have found the indices for these events we'll go through and do the epoching.
                        % We don't need to do this but it makes it much easier to read. We first set the epoch onset and offset, which we'll use to chunk the data later.
                        % The onset, offset, minDur and maxDur are rounded because the latencies are incredibly precise but we need these to be in samples so they need to be
                        % whole numbers.
                        onset = EEG2.event(onoffsets(trial, 1)).latency;
                        offset = EEG2.event(onoffsets(trial, 2)).latency;
                        try
                            EEG.data = EEG2.data(:, onset:offset);
                        catch
                            continue
                        end
                        EEG.times = EEG2.times(1, onset:offset);

                        % With the events we want the information to be relative to the new epochs rather than to the whole original file, this section takes in the original
                        % events, stores the first row of information and removes that from every subsequent row.

                        % Create a temp variable we can edit
                        temp = EEG2.event(onoffsets(trial, 1):onoffsets(trial, 2));
                        first_init_time = temp(1).latency_ms;

                        % Modifies all subsequent rows
                        for j = 1:length(temp)
                            temp(j).latency = round(temp(j).latency) - onset + 1;
                            % temp(j).description = num2str(str2num(temp(j).description) - str2num(first_init_index) + 1);
                            temp(j).latency_ms = temp(j).latency_ms - first_init_time + 1;
                        end

                        EEG.event = temp;
                        clear temp
                        % end

                        EEG.trials = 1;
                        EEG.pnts = offset - onset + 1;
                        EEG.xmax = max(EEG.data, [], 'all');
                        EEG.xmin = min(EEG.data, [], 'all');

                        % Example of how you could epoch the data if you want to instead epoch specific times relative to events
                        % EEG3 = pop_epoch(EEG, { 'DIN241' 'DIN242' }, [-0.1 0.5], 'newname', 'Social4Hz', 'epochinfo', 'yes');

                        % Save the data
                        save(strcat(ds.settings.paths.epochedEEGPath, ds.dataInfo.files(file).name(1:end-4), '_', endless_event_names{index}, '_', num2str(trial), '.mat'), 'ds', 'EEG');
                    end
                end
            end
        end
    else
        fprintf(strcat('Skipping\t', ds.dataInfo.files(file).name, '\t\tepoching complete\n'))
    end
end
end
