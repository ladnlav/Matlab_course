function [Eb_N0] = Eb_N0_convert(SNR, Constellation)
    [~, bit_depth] = constellation_func(Constellation);
    Rb=1;
    Eb_N0 = SNR / (bit_depth * Rb);

end

