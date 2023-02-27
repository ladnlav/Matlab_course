%% Task 2
clear; clc; close all;
load('Matlab_L3_6.mat');
%% Задание 1. Кадровая синхронизация
Stream=Bit_Stream;
Header = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
Special_bits=Header;

res_cor(1 : length(Stream)-length(Special_bits)) = 0;
for itter = 1 : length(Stream)-length(Special_bits)
    res_cor(itter) = sum(Special_bits.*Stream(itter + 1:itter+length(Special_bits) ))/length(Special_bits);
end

% находим пики корреляционной функции
[~, Indexes_of_frames] = findpeaks(res_cor, 'MinPeakHeight', 0.99999);
Start_Of_Frame_Position=Indexes_of_frames(1);

% Определяем количество кадров данных
Number_of_frame = length(Indexes_of_frames);

% Строим график корреляции и сохраняем его в файл
figure
plot(res_cor)
%xticks(0:10:10000); % установка значения для отображения на оси x
title('Корреляционный анализ битовой последовательности')
xlabel('Номер бита')
ylabel('Корреляция')
savefig('Frame_Corr.fig')

% Сохраняем переменные в файл
save('Frame_search.mat', 'Start_Of_Frame_Position', 'Number_of_frame')


%% Задание 2. Генератор псевдослучайной последовательности
%Функция Scrambler написана в файле "Scrambler.m"

Register = [1 0 0 1 0 1 0 1 0]; % начальное состояние регистра
sequence = Scrambler(Register); % генерация последовательности

%Функция циклической автокорреляции
N=length(sequence);
acf = zeros(1, N);
for itter = 1:N
    acf(itter) = sum(sequence.*circshift(sequence, itter-1))/N;
end
acf = acf / acf(1); % нормализация

figure;
plot(acf, 'LineWidth', 1.5);
title('Autocorrelation Function of Scrambler Output');
xlabel('Bit Offset');
ylabel('Autocorrelation');
saveas(gcf, 'ACF_Srambler.fig');

[~, max_index] = max(acf(2:end)); % поиск максимального значения автокорреляции
PN_Period = max_index; % вывод периода повторения