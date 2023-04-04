clear; clc; close all;
 %Load audio file
[x,fs] = audioread('Основной тон.wav');

% Convert to mono
x = mean(x, 2); % Take the mean of the left and right channels

% Design an FIR low-pass filter
N = 100; % Filter order
fc = 1000; % Cutoff frequency in Hz
Wn = fc/(fs/2); % Normalized cutoff frequency
b = fir1(N, Wn, 'low'); % FIR filter coefficients

% Apply the filter to the input signal
x = filter(b, 1, x);
%% 
% Compute the frequency response of the filter
[H, f] = freqz(b, 1, 1024, fs);

% Plot the magnitude response of the filter
figure;
plot(f, abs(H));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Response of FIR Low-Pass Filter');
grid on;

% Find the passband indices
passband_indices = find(f <= fc);

% Compute the passband ripple
passband_ripple = max(abs(abs(H(passband_indices)) - 1))

% Convert the passband ripple to decibels
passband_ripple_dB = 20*log10(passband_ripple)
%% 

% Apply window function
w = hamming(length(x));
xw = x.*w;

% Compute power spectrum
[S, f, T,Ps] = spectrogram(xw, [], [], [], fs, [],'onesided');

% Average the columns of the Ps matrix to obtain the PSD estimate
Pxx = mean(Ps, 2);

% Take logarithm of power spectrum
L = 10*log10(Pxx);

% Apply spectral smoothing
Ls = smoothdata(L,'gaussian',30);

% Select main pitch as the first maximum of the power logarithm in the range 70-450 Hz
[~, idx] = max(Ls(f >= 70 & f <= 450));
main_pitch = f(f >= 70 & f <= 450);
main_pitch = main_pitch(idx);
disp(['Main pitch: ', num2str(main_pitch), ' Hz'])

% Apply weighting function
f0 = main_pitch; % fundamental frequency

radius = 10; % radius of non-zero values in Hz

% Initialize weighting function to zeros
w = zeros(size(f));

% Find index of frequency closest to f0
[~, idx_f0] = min(abs(f - f0));

% Set values within radius to one and weights out of radius
w(idx_f0-radius:idx_f0+radius) = 1;
w(1:idx_f0-radius) = 0;
w(idx_f0+radius:end) = 0;
%w = 1./(1+(f./(0.2*f0)).^6);
Lsw = Ls.*w;

% Find spectral peaks
[pks,locs] = findpeaks(Lsw,'MinPeakDistance',round(fs/f0/4));

% Plot results
y0 = interp1(f,Lsw,f0);
y1=interp1(f,Ls,f0);

figure;
subplot(2,2,1);
plot(x);
xlabel('Time (samples)')
ylabel('Amplitude')
title('Original Signal');

subplot(2,2,2);
plot(f, Pxx);
xlim([0 4000])
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectrum');

subplot(2,2,3);
plot(f, Ls);
hold on
plot(f0,y1,'ro','MarkerSize',5,'LineWidth',1)
xlim([0 4000])
xlabel('Frequency (Hz)')
ylabel('Log Power (dB)')
title('Log Power Spectrum');

subplot(2,2,4);
plot(f,Lsw)
hold on
plot(f0,y0,'ro','MarkerSize',5,'LineWidth',1)
xlim([0 4000])
xlabel('Frequency (Hz)')
ylabel('Log Power (dB)')
title('Spectral Peaks')