function sequence = Scrambler(Register)
    % Инициализация параметров
    %m = 8; % число бит регистра
    poly = [1 1 0 0 0 0 0 0 0 0]; % полином обратной связи
    sequence_length = 2^8 - 1; % длина псевдослучайной последовательности
    sequence = zeros(1, sequence_length); % инициализация последовательности
    
    % Генерация псевдослучайной последовательности
    for i = 1:sequence_length
        sequence(i)=xor(Register(end), Register(end-1)); % выходной бит
        feedback = mod(sum(Register(poly == 1)), 2); % вычисление обратной связи

        Register = circshift(Register,1); % сдвиг регистра
        Register(1) = feedback; % обновление первого бита
    end
end