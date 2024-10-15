function [fftFFR, HzScale, dBfft] = myFFT(data,Fs,toPlot,xlimVar,fileName)
%
% Compute and plot FFT
%
% Version 1 - September 2017
%
% Tim Schoof - t.schoof@ucl.ac.uk
% ----------------------------------
if nargin < 5
   fileName = 'file name not specified';
elseif nargin < 4
    xlimVar = 1200;
end
% zero-pad data
 data = [data, zeros(1,length(data))];
% compute FFT
%L= 32768*32;  %test with nfft = cst = 32768*16 for shorter fft
%NFFT = L;
L=length(data);   
NFFT = 2^nextpow2(L);       % Next power of 2 from length of y
Y = fft(data,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
fftFFR = 2*abs(Y(1:NFFT/2+1));
% single-sided amplitude spectrum re peak amplitude - not used for plotting
dBref = 1;
dBfft = 20*log10(fftFFR/dBref);
HzScale = [0:(Fs/2)/(length(fftFFR)-1):round(Fs/2)]'; % frequency 'axis'
% plot FFT (if required)
if toPlot == 1
    plot(HzScale,fftFFR,'LineWidth',2)
    xlim([0 xlimVar])
    set(0, 'DefaulttextInterpreter', 'none')
    title(['',fileName,''])
    ylabel('Spectral magnitude')
    xlabel('Frequency (Hz)')
end


