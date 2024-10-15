function [fft, fig] = run_fft_and_plot(EEG, eeg_file, output_path, after_seg)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

%% Settings
fft = 0; fig = 0;
max_hz = 48;
Fs = EEG.srate;
filename = split(eeg_file, '\');
filename = filename(end);
filename = split(filename, '.');
filename = filename{1};

if after_seg
    % data = squeeze(mean(EEG.data, 1, 'omitnan'));
    data = EEG.data;
    t = 0:1/Fs:size(data, 2)-1/Fs; % Time vector
    % Define desired frequency resolution
    desired_resolution = 0.5; % Hz

    % Calculate the necessary segment length (N) to achieve the desired frequency resolution
    % Resolution = fs / N
    N = round(Fs / desired_resolution); % Segment length

    fft = [];
    for fftTrial = 1:size(data, 1)
        for fftElec = 1:size(data, 2)
            [fftRes, fftHzScale, dBfft] = myFFT(squeeze(data(fftTrial, fftElec, :))',Fs,0,50);
            [minValue, closestIndex] = min(abs(fftHzScale-max_hz)); % We're filtering at 100Hz so no point keeping it.
            fftRes = fftRes(1:closestIndex);
            fftHzScale = fftHzScale(1:closestIndex);
            dBfft = dBfft(1:closestIndex);

            fft(1, fftTrial, fftElec, :) = fftRes;
            fft(2, fftTrial, fftElec, :) = fftHzScale;
            fft(3, fftTrial, fftElec, :) = dBfft;
        end
    end

    fft = squeeze(mean(fft, 2, 'omitnan'));
else
    data = EEG.data;
    t = 0:1/Fs:EEG.pnts-1/Fs; % Time vector
    % Define desired frequency resolution
    desired_resolution = 0.1; % Hz

    % Calculate the necessary segment length (N) to achieve the desired frequency resolution
    % Resolution = fs / N
    N = round(Fs / desired_resolution); % Segment length

    %% Calc fft
    fft = [];
    for fftElec = 1:size(data, 1)
        [fftRes, fftHzScale, dBfft] = myFFT(squeeze(data(fftElec, :)),Fs,0,50);
        [minValue, closestIndex] = min(abs(fftHzScale-max_hz)); % We're filtering at 100Hz so no point keeping it.
        fftRes = fftRes(1:closestIndex);
        fftHzScale = fftHzScale(1:closestIndex);
        dBfft = dBfft(1:closestIndex);

        fft(1, fftElec, :) = fftRes;
        fft(2, fftElec, :) = fftHzScale;
        fft(3, fftElec, :) = dBfft;
    end
end

% Uncomment to plot the fft figure
% %% Plot the fft figure
% % Set the figure to not be visible and background color to white
% fig = figure('Visible', 'off', 'Color', 'white');
% % Cycle through all electrodes
% for i = 1:size(fft, 2)
%     plot(fftHzScale, squeeze(fft(1,i,:))); hold on;
% end
% % Set x limit and title
% xlim([1, max_hz]); title(filename, 'Interpreter', 'none');
% % Save image
% saveas(fig, fullfile(output_path, [filename, '.png'])); % Save as PNG file

% Uncomment below if you would instead like to use a pwelch PSD method
% %% Calc pwelch PSD
% 
% 
% % Grab the averaged data
% x = mean(data, 1, 'omitnan');
% 
% % Set other parameters for pwelch
% window = hamming(N); % Use a Hamming window of length N
% noverlap = N / 2; % 50% overlap
% nfft = N; % Number of FFT points
% 
% % Compute the Power Spectral Density using pwelch
% [pxx, f] = pwelch(x, window, noverlap, nfft, Fs);
% 
% % Plot the PSD
% fig = figure('Visible', 'off', 'Color', 'white');
% subplot(2,1,1)
% plot(f, 10*log10(pxx)); % Convert to decibels (dB)
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% title('Average Power Spectral Density');
% 
% % Now per electrode
% subplot(2,1,2)
% 
% for i = 1:size(EEG.data, 1)
%     % Compute the Power Spectral Density using pwelch
%     [pxx, f] = pwelch(EEG.data(i, :), window, noverlap, nfft, Fs);
%     plot(f, 10*log10(pxx)); % Convert to decibels (dB)
%     hold on;
% end
% 
% grid on;
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% % legend(EEG.chanlocs)
% title('PSD by electrode');
% 
% saveas(fig, fullfile(output_path, [filename, '_PSD.png'])); % Save as PNG file

end