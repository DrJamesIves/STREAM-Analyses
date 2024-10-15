function eeg_scroll_plot(EEG, startEEGLAB)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% This function takes in an eeglab formatted EEG dataset and plots it sequentially so you can quickly check through all your data.

% Set in seconds
pause_length = 0.1;

% Start EEGLAB
if startEEGLAB
    eeglab;
end

if nargin == 0
    % Load an EEG dataset (replace 'filename' with your EEG data file), this is an example
    load('E:\Birkbeck\STREAM\Datasets\2. Preprocessed\2.2 Preprocessed_EEG\2.2.1 Full\MW-0186_fast_erp_14.mat')
end

% Extract data
data = EEG.data; % Data is in the format [channels x points]

% Get the time vector in seconds
time = linspace(EEG.xmin, EEG.xmax, EEG.pnts);

% Normalize the data for each channel
data = bsxfun(@minus, data, mean(data, 2));
data = bsxfun(@rdivide, data, std(data, [], 2));

% Define the window size for scrolling (e.g., 2 seconds of data)
window_size = 2 * EEG.srate;

% Define the step size for scrolling (e.g., 0.5 second)
step_size = 0.5 * EEG.srate;

% Number of scrolling windows
num_windows = floor((EEG.pnts - window_size) / step_size) + 1;

% Create the figure
figure;
set(gcf, 'color', 'w')

for i = 1:num_windows
    % Calculate the indices for the current window
    idx_start = (i-1) * step_size + 1;
    idx_end = idx_start + window_size - 1;
    
    % Extract current window data and compute the scale
    current_window_data = data(:, idx_start:idx_end);
    data_max = [];
    data_min = [];

    for j = 1:size(current_window_data, 1)
        data_range = [max(current_window_data(j, :)) - min(current_window_data(j, :))];
    end
    scale_height = max(data_range);
    % data_min = [dmin(current_window_data(:));
    % scale_height = data_max - data_min;

    % Plot all electrodes in the current window
    clf; % Clear the figure
    hold on;
    for chan = 1:EEG.nbchan
        plot(time(idx_start:idx_end), data(chan, idx_start:idx_end) + chan * 10); % Offset each channel
    end
    hold off;
    
    % Set plot labels and limits
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(sprintf('EEG Data: Window %d/%d', i, num_windows));
    xlim([time(idx_start) time(idx_end)]);
    ylim([0 10 * (EEG.nbchan + 1)] ); % Adjust based on your data
    
    % Add scale as a text annotation or legend
    text(time(idx_start) + 0.1, 10 * EEG.nbchan, ...
        sprintf('Max range: %.2f μV', scale_height), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

    % Pause to create scrolling effect
    pause(pause_length);
end

end