close all; clear; clc;
%% Init parametrs of model
Length_Bit_vector = 1e5;
rng(321); % Fix the seed of the random number generator

Constellation = "16QAM"; % BPSK, QPSK, 8PSK, 16QAM

SNR = 30; % dB

%% Bit generator

Bit_Tx = generateBits(Constellation,Length_Bit_vector);

%% Mapping

IQ_TX = mapping(Bit_Tx, Constellation);

%plot_const(Constellation);
%% Channel
% Write your own function Eb_N0_convert(), which convert SNR to Eb/N0
Eb_N0 = Eb_N0_convert(SNR, Constellation);

% Use your own function of generating of AWGN from previous tasks
IQ_RX = Noise(SNR, IQ_TX);
%% Demapping
Bit_Rx = demapping(IQ_RX, Constellation);

%% Error check
% Write your own function Error_check() for calculation of BER
BER = Error_check(Bit_Tx, Bit_Rx);

%% Experimental BER(SNR) and BER(Eb/N0)
% Collect enough data to plot BER(SNR) and BER(Eb/N0) for each
% constellation.
% Compare the constalation. Describe the results
% You can use the cycle for collecting of data
% Save figure

constellations = ["BPSK", "QPSK", "8PSK", "16QAM"];

SNR = -50:0.05:50;
BERm_all=zeros(size(SNR,2),size(constellations,2));
EbN0_all=zeros(size(SNR,2),size(constellations,2));


for p = 1:length(constellations)
    Constellation = constellations{p};
    Bit_Tx = generateBits(Constellation,Length_Bit_vector);
    IQ_TX = mapping(Bit_Tx, Constellation);

    EbN0 = Eb_N0_convert(SNR, Constellation);
    BERm = zeros(size(SNR));
    MERm=  zeros(size(SNR));
    tic
    parfor i = 1:length(SNR)
        IQ_RX = Noise(SNR(i), IQ_TX);
        %IQ_RX= awgn(IQ_TX,SNR(i),'measured');
        Bit_Rx = demapping(IQ_RX, Constellation);
        BERm(i) = Error_check(Bit_Tx, Bit_Rx);
        MERm(i) = MER_my_func(IQ_RX, Constellation);
    end
    toc
    BERm_all(:,p)=BERm;
    EbN0_all(:,p)=EbN0;

    figure('Position', [100 0 1720 500]);
    subplot(1, 2, 1);
    
    lastNonZero = find(BERm ~= 0, 1, 'last');
    BERm(lastNonZero+1:end) = 4.000000000000000e-323;

    plot(SNR, BERm,'r','LineWidth',1.5);
    hold on;
    plot(EbN0, BERm,'b','LineWidth',1.5);
    hold on;
    EbN0_c = 10.^(EbN0/10);
    BERt = 1/2.*erfc(sqrt(EbN0_c));
    
    lastNonZero = find(BERt ~= 0, 1, 'last');  
    BERt(lastNonZero+1:end) = 4.000000000000000e-323;
    plot(EbN0, BERt,'g','LineWidth',1.5);
    hold on;
    legend('Зависимость BER от SNR','Экспериментальный BER от Eb/N0', ...
        'Теоретический BER от Eb/N0', 'Location','southwest');
    %set(gca, 'YScale', 'log');
    title(['Зависимости для ' Constellation]);
    grid on;
    xlabel('SNR/Eb_N0 (dB)');
    ylabel('BER');
    xlim([-50 50]);

%% Additional task. Modulation error ration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5MER_estimation = MER_my_func(IQ_RX, Constellation);
% Compare the SNR and MER_estimation from -50dB to +50dB for BPSK, QPSK,
% 8PSK and 16QAM constellation.
% Plot the function of error between SNR and MER for each constellation 
% Discribe the results. Make an conclusion about MER.
% You can use the cycle for collecting of data
% Save figure

    subplot(1, 2, 2);
    plot(SNR, abs(MERm-SNR),'magenta','LineWidth',1.5);
    xlabel('SNR (dB)');
    ylabel('Error(MERm-SNR)');
    title(['Ошибка между SNR и MER для ' Constellation]);
    grid on;

    name=Constellation+"_BER от SNR и от Eb_N0 и MER.png";
    saveas(gcf,name);
end

%% All in one
figure();
plot(EbN0_all(:,1), BERm_all(:,1),'r','LineWidth',2);
hold on;
plot(EbN0_all(:,2), BERm_all(:,2),'b','LineWidth',2);
hold on;
plot(EbN0_all(:,3), BERm_all(:,3),'g','LineWidth',2);
hold on;
plot(EbN0_all(:,4), BERm_all(:,4),'black','LineWidth',2);
legend("BPSK", "QPSK", "8PSK", "16QAM", 'Location','southwest');
set(gca, 'YScale', 'log');
xlabel('Eb/N0 (dB)');
ylabel('BER');
title('Зависимость BER от Eb/N0 для ' + strjoin(constellations,', '));
grid on;
hold on;
name="All in One.png";
saveas(gcf,name);
%% Functions

function plot_const(const_type)
% Plot constellation diagram for specified type (BPSK, QPSK, 8PSK, or 16-QAM)

% Get constellation points and depth
[points,~] = constellation_func(const_type);

% Map bit sequences to constellation points
num_points = length(points);
bits_m = de2bi(0:num_points-1, log2(num_points), 'left-msb')';
bits=bits_m(:)';
mapped_points = mapping(bits,const_type);

% Plot points on IQ plane
figure;
scatter(real(mapped_points), imag(mapped_points),100, 'g', 'filled');
hold on;
title(sprintf('%s Constellation Diagram', upper(const_type)));

% Label points with encoded bit sequence
for i = 1:num_points
    text(real(points(i)), imag(points(i))-0.15, ...
         sprintf('%s', num2str(bits_m(:,i))), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

xL = [min(real(points)) max(real(points))];
yL = [min(imag(points)) max(imag(points))];
xlim(xL + [-1 1]);
ylim(yL + [-1 1]);

% Add arrows to the x and y axes using xlim and ylim values
line([0 0], [0 100], 'Color', 'k');
line([100 0], [0 0], 'Color', 'k');
line([0 0], [0 -100], 'Color', 'k');
line([-100 0], [0 0], 'Color', 'k');
hold off;
grid on;
xlabel('In-Phase');
ylabel('Quadrature');
end

%FLOWBITS_GEN
function bits = generateBits(constellation,Length_Bit_vector)

    [~, Bit_depth_Dict] = constellation_func(constellation);

    bits = randi([0 1], 1, Length_Bit_vector);

    % Calculate the number of zeros needed to pad the bits vector
    num_zeros = Bit_depth_Dict - mod(length(bits), Bit_depth_Dict);
    if num_zeros == Bit_depth_Dict
        num_zeros = 0;
    end
    
    % Pad the bits vector with zeros
    bits = [bits zeros(1, num_zeros)];
end