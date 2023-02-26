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

% находим максимум корреляционной функции
[~, Start_Of_Frame_Position] = max(res_cor); 

% Определяем количество кадров данных
Number_of_frame = sum(ismember(Stream, Header));

% Строим график корреляции и сохраняем его в файл
figure
plot(res_cor)
title('Корреляционный анализ битовой последовательности')
xlabel('Номер бита')
ylabel('Корреляция')
savefig('Frame_Corr.fig')

% Сохраняем переменные в файл
save('Frame_search.mat', 'Start_Of_Frame_Position', 'Number_of_frame')


%%