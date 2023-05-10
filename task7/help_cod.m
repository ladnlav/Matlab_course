clc; clear;
% Вводим сообщение на русском языке
message = 'I done this task! Im winner!';

% Преобразуем каждый символ в двоичный код
binary = dec2bin(message);

% Преобразуем двоичный код в массив из нулей и единиц
array = reshape(binary - '0', 1, []);

len2=length(array);
Register = [1 0 0 1 0 1 0 1 0];


seq=Scrambler(Register);
seq(seq==0)=-1;
len=length(seq);

totallen=len2;
maxdev=max_divisor_range(totallen);

Register = [1 0 0 1 0 1 0 1 0]; % начальное состояние регистра
scrseq=Scrambler(Register);

function d = max_divisor_range(n)
% Проверяем, является ли n положительным целым числом
if n < 1 || mod(n,1) ~= 0
    error('n должно быть положительным целым числом')
end
% Инициализируем d как 1
d = 1;
% Перебираем возможные делители от 2 до 32
for i = 2:32
    % Если i делит n без остатка и больше текущего значения d, то обновляем d как i
    if mod(n,i) == 0 && i > d
        d = i;
    end
end
end

function seq = Scrambler(Register)
    % Инициализация параметров
    %m = 8; % число бит регистра
    sequence_length = 2^7-1; % длина псевдослучайной последовательности
    sequence = zeros(1, sequence_length); % инициализация последовательности
    
    % Генерация псевдослучайной последовательности
    for i = 1:sequence_length
        sequence(i)=xor(Register(end), Register(end-1)); % выходной бит
        feedback = sequence(i); 

        Register = circshift(Register,1); % сдвиг регистра
        Register(1) = feedback; % обновление первого бита
    end
    seq=sequence.';
end