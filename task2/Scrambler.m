function sequence = Scrambler(Register)
    % Инициализация параметров
    %m = 8; % число бит регистра
    sequence_length = 2^9 - 1; % длина псевдослучайной последовательности
    sequence = zeros(1, sequence_length); % инициализация последовательности
    
    % Генерация псевдослучайной последовательности
    for i = 1:sequence_length
        sequence(i)=xor(Register(end), Register(end-1)); % выходной бит
        feedback = sequence(i); 

        Register = circshift(Register,1); % сдвиг регистра
        Register(1) = feedback; % обновление первого бита
    end
end