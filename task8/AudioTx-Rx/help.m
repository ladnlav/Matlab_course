clc; clear;

[input_audio, sample_rate] = audioread('send_8.mp3');

% Create a low-pass filter
cutoff_frequency = 2000; % Choose a cutoff frequency (in Hz)
nyquist_frequency = sample_rate / 2;
normalized_cutoff = cutoff_frequency / nyquist_frequency;
filter_order = 6;
[b, a] = butter(filter_order, normalized_cutoff, 'low');

%compressed_audio = filter(b, a, input_audio);
compressed_audio=input_audio;
compressed_audio=(compressed_audio)*256;
compressed_audio=int8(floor(compressed_audio));

Y = abs(fft(compressed_audio)); % применение FFT
Y1=Y(1:(length(compressed_audio)/2));
f=(0:(length(compressed_audio)/2)-1);

% % построение графика амплитудного спектра
% figure;
% subplot(2,1,1);
% plot(f,Y1);
% title('Амплитудный спектр сигнала');
% xlabel('Частота (Гц)');
% ylabel('Амплитуда');