clear; clc; close all;

%% Упражнение 1. Fourier transform

S = zeros(128,1); %Инициализация исходного сигнала
S(5) = 1;
S(10) = 2;
S(60) = 3;

Signal = (ifft(S)); %Обратное преобразование Фурье

%figure;
%plot(abs(Signal));     %График исходного Signal
%hold on;
%title("abs(Signal)");

%F=fft(Signal);

%figure;
%plot(abs(F));     %График спектра
%hold on;
%title("Signal spectrum");

%% Упражнение 2. Noise generation

SNR=10;
[NoisedSignal,Noise]=NoiseGenerator(SNR,Signal);

%% Упражнение 3. Powers of signals

% Рассчитываем среднюю мощность исходного сигнала
P_Signal = PowerSignal(Signal);

% Pассчитываем среднюю мощность шума
P_Noise = PowerSignal(Noise);

% Pассчитываем среднюю мощность зашумленного сигнала
P_NoisedSignal = PowerSignal(NoisedSignal);

%% Упражнение 4. Parseval theorem

%Вычисляем спектры сигнала, шума и зашумленного сигнала
SignalSpec=fft(Signal);
NoiseSpec=fft(Noise);
NoisedSignalSpec=fft(NoisedSignal);

%Вычисляем мощности рассматриваемых сигналов, нормируя их на 128 -
%размерность преобразования Фурье
P_SignalSpec=PowerSignal(SignalSpec)/128;
P_NoiseSpec=PowerSignal(NoiseSpec)/128;
P_NoisedSignalSpec=PowerSignal(NoisedSignalSpec)/128;

% проверяем теорему Парсеваля
if abs(P_Signal - P_SignalSpec) / P_Signal < 0.001 && ...
   abs(P_Noise - P_NoiseSpec) / P_Noise < 0.001 && ...
   abs(P_NoisedSignal - P_NoisedSignalSpec) / P_NoisedSignal < 0.001
    disp('True');
else
    disp('False');
    disp(abs(P_SignalSpec) / P_Signal);
    disp(abs(P_Noise - P_NoiseSpec) / P_Noise);
    disp(abs(P_NoisedSignal - P_NoisedSignalSpec) / P_NoisedSignal);
end

%% Упражнение 5. Signal filtering
FilteredNoisedSignal=FilterSignal(NoisedSignal);

%% Упражнение 6. SNR comparison

% Разделение фильтрованного сигнала и шума
FilteredNoise = FilteredNoisedSignal - NoisedSignal;

% Расчет SNR для исходного зашумленного сигнала
SNR_NoisedSignal = 10 * log10(PowerSignal(NoisedSignal) / PowerSignal(Noise));

% Расчет SNR для отфильтрованного сигнала
SNR_FilteredNoisedSignal = 10 * log10(PowerSignal(FilteredNoisedSignal) / PowerSignal(FilteredNoise));

if SNR_FilteredNoisedSignal > SNR_NoisedSignal
    disp('Фильтрованный сигнал лучше, так как SNR выше');
else
    disp('Нефильтрованный сигнал лучше, так как SNR выше');
end

disp(['SNR нефильтрованного сигнала: ' num2str(SNR_NoisedSignal) ' дБ']);
disp(['SNR отфильтрованного сигнала: ' num2str(SNR_FilteredNoisedSignal) ' дБ']);

%% Выигрыш по SNR между фильтрованным и нефильтрованным сигналами

K=1000;
Gain = zeros(K, 1);
for itter = 1:K
    SNR=itter;
    [NoisedSignal,Noise]=NoiseGenerator(SNR,Signal);
    FilteredNoisedSignal=FilterSignal(NoisedSignal);

    % Разделение фильтрованного сигнала и шума
    FilteredNoise = FilteredNoisedSignal - NoisedSignal;

    % Расчет SNR для исходного зашумленного сигнала
    SNR_NoisedSignal = 10 * log10(PowerSignal(NoisedSignal) / PowerSignal(Noise));

    % Расчет SNR для отфильтрованного сигнала
    SNR_FilteredNoisedSignal = 10 * log10(PowerSignal(FilteredNoisedSignal) / PowerSignal(FilteredNoise));

    %Учет выигрыша
    Gain(itter)=abs(SNR_NoisedSignal-SNR_FilteredNoisedSignal);
end

figure;
plot(Gain);     %График Gain
hold on;
title("SNR Gain from input SNR");
xlabel('Input SNR');
ylabel('SNR Gain');
saveas(gcf, 'SNR Gain.png');
%% Функции

function P = PowerSignal(Signal)
    P =mean(abs(Signal).^2);
end

function [NoisedSignal, Noise] = NoiseGenerator(SNR, Signal)
    % вычисляем мощность сигнала
    SignalPower = PowerSignal(Signal);
    % вычисляем мощность шума
    NoisePower = SignalPower / 10^(SNR/10);
    % генерируем белый шум
    Noise = sqrt(NoisePower/2) * normrnd(0, 1, size(Signal)) + ...
            1i * sqrt(NoisePower/2) * normrnd(0, 1, size(Signal));
    % добавляем шум к сигналу
    NoisedSignal = Signal + Noise;
end

function FilteredNoisedSignal = FilterSignal(NoisedSignal)
    %Вычисляем спектр зашумлённого сигнала
    NoisedSignalSpec=fft(NoisedSignal);
    % Задаем фильтр
    filter = zeros(128, 1);
    filter(1:70) = 1;
    % Вычисляем отфильтрованный сигнал
    FilteredNoisedSignal = ifft(NoisedSignalSpec.* filter);
end
