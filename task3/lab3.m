clear; clc; close all;

%% Задача 1: Анализ зашумленного звукового файла
%Пункт 1: Найти частоту гармоник
% загрузка аудиофайла
[y,Fs] = audioread('file6.wav');

Y = abs(fft(y)); % применение FFT
Y1=Y(1:(length(y)/2));
f=(0:(length(y)/2)-1);

% построение графика амплитудного спектра
figure;
subplot(2,1,1);
plot(f,Y1);
title('Амплитудный спектр сигнала');
xlabel('Частота (Гц)');
ylabel('Амплитуда');

% нахождение ярких гармоник
[~,idx] = sort(Y1,'descend');
top_harmonics = f(idx(1:3)); % выбор трех ярких гармоник

% Вывод результата
disp(['Частоты поврежденных гармоник: ', num2str(sort(top_harmonics)), ' Гц']);

%Пункт 2: Отфильтровать аудиофайл так, чтоб снизить помехи
% создание фильтра для ярких гармоник
freq_tol = 0; % ширина полосы частот вокруг гармоники, которую мы оставим в сигнале
harmonic_filter = zeros(size(Y));
for k = 1:length(top_harmonics)
    harmonic_idx = find(f >= (top_harmonics(k)-freq_tol) & f <= (top_harmonics(k)+freq_tol));
    harmonic_filter(harmonic_idx-1) = 1;
end

% применение фильтра к сигналу
y_filtered = real(ifft(y.*harmonic_filter));
y_filtered=y_filtered./max(y_filtered);
% вывод исходного и отфильтрованного сигнала
figure;
dur=length(y)/Fs;
t = linspace(0, dur, length(y));
subplot(2,1,1);
plot(t,y);
title('Исходный сигнал');
xlabel('Время, с');
ylabel('Амплитуда');

subplot(2,1,2);
plot(t,y_filtered);
title('Отфильтрованный сигнал');
xlabel('Время, с');
ylabel('Амплитуда');

%sound(y,Fs)
% Сохранение отфильтрованных данных в новом аудиофайле
audiowrite('filtered_audio_file.wav', y_filtered, Fs);

%% Задача 2: анализ клиппинг-эффекта.
clear; clc; close all;

%Пункт 1: Создайте звуковой файл из синусоиды 
% задаем параметры сигнала
fs = 2500;          % частота дискретизации, Гц
dur = 3;            % длительность, с
f0 = 1000;          % частота, Гц
A = 3;              % амплитуда, В

% создаем временной вектор
t = linspace(0, dur, dur*fs);

% создаем сигнал
x = A*sin(2*pi*f0*t);

%Пункт 2: Сымитируйте клиппинг-эффект 
% задаем максимальную амплитуду
max_amp = 2;

% обрезаем амплитуду сигнала
x_clip = clip(x, max_amp);

% создаем временные векторы
t_clip = t;

% создаем графики
subplot(2,1,1);
plot(t, x);
xlabel('Время, с');
ylabel('Амплитуда, В');
title('Исходный сигнал');

subplot(2,1,2);
plot(t_clip, x_clip);
xlabel('Время, с');
ylabel('Амплитуда, В');
title('Обрезанный сигнал');

% вычисляем спектры
X = abs(fft(x));
X_clip = abs(fft(x_clip));

% максимальная частота, Гц
fmax = fs/2;

% создаем вектор частот
fvec = linspace(0, fmax, length(X)/2);

% создаем графики спектров
figure();
subplot(2,1,1);
plot(fvec, X(1:length(X)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр исходного сигнала');

subplot(2,1,2);
plot(fvec, X_clip(1:length(X_clip)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр обрезанного сигнала');

audiowrite('xsin.wav', x./max(x), fs);
audiowrite('xsin_clipped.wav', x_clip./max(x_clip), fs);
%Проанализируйте, как это отразится на характеристиках сигнала
% Как видно на графиках спектров двух сигналов, после клиппинга появилась
% побочная гармоника на частоте примерно в 500 Гц. В исходном сигнале её
% не было. Из этого можно сделать вывод, что клиппинг привёл к искажению
% исходного сигнала, добавив в сигнал побочную гармонику. Кроме того,
% распределение энергии в спектре изменилось (часть перешла на побочную
% гармонику). Клиппинг может приводить к потере деталей и искажению звука. 
% Это может проявляться в виде неприятного звука или шума, который слышен 
% в тех местах, где обрезались вершины сигнала.

%% Задача 3: анализ влияния частоты дискретизации.
clear; clc; close all;

%Пункт 1: downsampling
% Загрузка аудиофайла и получение его характеристик
[x, fs] = audioread('task3.wav'); % загрузка аудиофайла

% Изменение частоты дискретизации вдвое
fs_downsampled = fs/2; % новая частота дискретизации

x_downsampled = zeros(round(length(x)/2), 2); 
x_downsampled(:,1) = x(1:2:end,1); % изменение частоты дискретизации на первом канале
x_downsampled(:,2) = x(1:2:end,2); % изменение частоты дискретизации на втором канале


% Проанализировать как это отразится на субъективном качестве и характеристиках сигнала
% Сохранение resampled data в новом аудиофайле
%sound(x_downsampled,fs_downsampled)
audiowrite('task3_downsampled.wav', x_downsampled, fs_downsampled);

%Проанализируйте, как это  отразится на субъективном качестве 
% и характеристиках сигнала 
% Субъективное восприятие: Узнать аудиозапись всё ещё можно. Однако сам 
% звук мелодии стал менее насыщенным, более приглушенным.

%Графики спектров сигналов
Y = abs(fft(x));
fvec = linspace(0, fs/2, length(Y)/2);

figure();
subplot(2,1,1);
plot(fvec, Y(1:length(Y)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр исходного сигнала');
xlim([0 fs_downsampled/2]);

Y_downsampled = abs(fft(x_downsampled));
fvec = linspace(0, fs_downsampled/2, length(Y_downsampled)/2);
subplot(2,1,2);
plot(fvec, Y_downsampled(1:length(Y_downsampled)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр сигнала с половинной частотой дискретизации');

% Анализ характеристик сигналов: как можно видеть на графиках все основные
% частоты остались прежними, их величина друг относительно друга осталась
% такой же. Однако изменилась абсолютная амплитуда сигналов на одних и тех
% же частотах. Вероятно, это и приводит к тому, что сигнал с половинной
% частотой дискретизации слышится более "приглушенным". 


% Пункт 2: upsampling

fs_upsampled = fs*2;
x_upsampled = zeros(2*length(x), 2);
x_upsampled (1:2:end,:) = x;
x_upsampled (2:2:end-2,:) = (x(1:end-1,:)+x(2:end,:))/2;
x_upsampled (end,:) = x(end,:);

% Сохранение resampled data в новом аудиофайле
audiowrite('task3_upsampled .wav', x_upsampled , fs_upsampled);

%Проанализируйте, как это  отразится на субъективном качестве 
% и характеристиках сигнала 
% Субъективное восприятие: Аудиозапись преобразилась. Звук стал более
% насыщенным и отчетливым. Если так можно сказать, то "яркость" каждой
% частоты увеличилась.

%Графики спектров сигналов
Y = abs(fft(x));
fvec = linspace(0, fs/2, length(Y)/2);

figure();
subplot(2,1,1);
plot(fvec, Y(1:length(Y)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр исходного сигнала');


Y_upsampled = abs(fft(x_upsampled));
fvec = linspace(0, fs_upsampled/2, length(Y_upsampled)/2);
subplot(2,1,2);
plot(fvec, Y_upsampled(1:length(Y_upsampled)/2));
xlabel('Частота, Гц');
ylabel('Амплитуда');
title('Спектр сигнала с удвоенной частотой дискретизации');
xlim([0 fs_upsampled/4]);

% Анализ характеристик сигналов: как можно видеть на графиках все основные
% частоты остались прежними, их величина друг относительно друга осталась
% такой же. Однако изменилась абсолютная амплитуда сигналов на одних и тех
% же частотах. Вероятно, это и приводит к тому, что сигнал с удвоенной
% частотой дискретизации слышится более отчётливо и ярко.

%% Дополнительная задача
clear; clc; close all;

% Загрузка аудиофайла
[x, fs] = audioread('budgie-chirping.wav');

% Построим спектрограмму
window_size = round(fs*0.0005); % размер окна 0,5 мс
noverlap = round(window_size/2); % перекрытие между окнами 50%
nfft = 2^nextpow2(window_size); % размер fft
[S, F, T,P] = spectrogram(x(:,1), window_size, noverlap, nfft, fs);

% Вычислим модуль мощности сигнала
Pmag = abs(P);

%Найдите максимальное значение модуля на спектрограмме вдоль оси времени 
maxVals = max(Pmag,[],1);

%Установим пороговое значение, при превышении которого будет звук "ч"
threshold = mean(maxVals) + 40*std(maxVals);

%Найдиём временные индексы, в которых максимальные значения модуля 
% превышают пороговое значение.
spikeTimes = T(maxVals > threshold);

%Построим спектрограмму и отметим время звука "ч" красными вертикальными
% линиями 
figure;
spectrogram(x(:,1), window_size, noverlap, nfft, fs,'yaxis');
hold on;
line(repmat(spikeTimes,2,1),repmat([0;max(F)],1,length(spikeTimes)),'Color','red');

%% Masker
% длительность сигнала-маскера в сек
masker_dur = 0.1; 
% амплитуда сигнала-маскера
%masker_amp = 0.05; 
mamp=max(x(:));
masker_t = linspace(0,masker_dur,round(masker_dur*fs));

%В качестве сигнала-маскера возьмём прямоугольный импульс
masker_signal = mamp * square(2*pi*600*masker_t);

%Определим временные индексы звука "ч"
outburst_indices = round(spikeTimes*fs);

for i = 1:length(outburst_indices)
    % Вычислим начальный и конечный индексы сигнала-маскера
    masker_start = max(1,outburst_indices(i)-round(length(masker_signal)/2));
    masker_end = min(length(x),outburst_indices(i)+round(length(masker_signal)/2)-1);
    
    % Добавим сигнал-маскер к исходному сигналу
    x(masker_start:masker_end) = x(masker_start:masker_end) + masker_signal(1:masker_end-masker_start+1);
end

x = x / max(abs(x)); % нормализация (опционально)
x(:,2) = x(:,1);

%sound(x,fs)
audiowrite('converted_audio.wav',x,fs); % Сохраняем новый сигнал в файл


%% Функции

%функция для создания эффекта-клиппинга
function y = clip(x, max_amp)
    y = x;
    y(y > max_amp) = max_amp;
    y(y < -max_amp) = -max_amp;
end