function [de_bits] = demapping(IQ, Constellation)
% Make the different dictionary for BPSK, QPSK, 8PSK, 16QAM constellations
% calculate the Bit_depth for each contellation
[Dictionary, Bit_depth_Dict] = constellation_func(Constellation);

% Find the closest constellation point for each IQ value

distances =abs(IQ - Dictionary.');

[~, idx] = min(distances, [], 1);

% Convert the indices to binary representation
de_bits = int2bit(idx-1, Bit_depth_Dict);

% Reshape the bits into a row vector
de_bits = reshape(de_bits, 1, []);

end

