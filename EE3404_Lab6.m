%% EE3404 Lab 6: QAM Demodulation | Aimee Nogoy | akn264

%% 1. Create a QAM constellation map
clear; clc; close all

%1 for QPSK, 2 for 16-QAM, 3 for 64-QAM

nbitsDim = input('Enter the modulation type (1 for QPSK, 2 for 16-QAM, or 3 for 64-QAM): ');
    if nbitsDim==1
        scale = 1/sqrt(2); %scale so avg power is 1/2 in 1D
        map1 = scale*[-1,1]';
    elseif nbitsDim==2
        scale = 1/sqrt(10);
        map1 = scale*[-3 -1 1 3]';
    else
        scale =1/sqrt(21);
        map1 = scale*[-7 -5 -3 -1 1 3 5 7]';
    end

%% 2. Gray code the constellation
%gray code mapping~~
switch nbitsDim
    case 1
        Igray = [0 1]';
    case 2
        Igray = [0 1 3 2]';
    case 3
        Igray = [0 1 3 2 6 7 5 4]';
    otherwise
        error('nbitsDim must 1, 2, or 3');
end
map1 = map1(Igray+1);

%% 3. Generate and modulate random bits
clc;
nsym = 1e5; 
modrate = 2*nbitsDim;
nbits = modrate*nsym;
bits = randi([0 1],nbits,1);

% sad = [bin2dec(bits(:,1:nbitsDim)),bin2dec(bits(:,(nbitsDim+1):end))];
% sym = map1(sad(:,1)+1) + 1i*map1(sad(:,2)+1);
%
switch nbitsDim
    case 1
        bi = bits(1:modrate:nbits); %inphase component
        bq = bits(2:modrate:nbits); %quadrature component 
        si = map1(bi+1); sq = map1(bq+1);
        sym = si+(1i*sq);
    case 2
        bi = bits(1:modrate:nbits) + bits(2:modrate:nbits);
        bq = bits(3:modrate:nbits) + bits(4:modrate:nbits);
        si = map1(bi+1); sq = map1(bq+1);
        sym = si+(1i*sq);
    case 3
        bi = bits(1:modrate:nbits) + bits(2:modrate:nbits) +...
            bits(3:modrate:nbits);
        bq = bits(4:modrate:nbits) + bits(5:modrate:nbits) +...
            bits(6:modrate:nbits);
        si = map1(bi+1); sq = map1(bq+1);
        sym = si+(1i*sq);
end

%% 4. Add Gaussian noise
% Add noise
EsN0 = 25; %in dB
wvar = 1/(4*(10^(EsN0/10)));
w = (randn(nsym,1) + 1i*randn(nsym,1))*sqrt(wvar);
y = sym+w;

%% 5. Perform the ML detection on the symbols

% get I and Q parts of y
z = [real(y), imag(y)];
% use repmat!!!
m = (1:(2^nbitsDim))';
mrep = repmat(map1',nsym,1);
Irep = repmat(z(:,1),1,length(m));
Qrep = repmat(z(:,2),1,length(m));
[~,a] = min(abs(Irep-mrep),[],2);
[~,b] = min(abs(Qrep-mrep),[],2);
mhat = [a,b];

%% 6. Extract the bits
%use dec2bin (coulda used it before)

RXbits = ([dec2bin(a-1),dec2bin(b-1)]);

%% 7. Measure the BER (bit error rate)
clc;
RXbits = str2num(RXbits(:));
bits = bits(:);
err = mean(RXbits ~= bits);
fprintf('The BER for the sample is %.2f %\n', err);

%% 7. Varying Es/N0 from -5 to 30 dB
clc;
clear; clc; close all;
EsN0Test = -5:30;
nsym = 1e5;
BER = zeros(length(EsN0Test), 3);
for nBitsDim = 1:3;
    switch nBitsDim,
        case 1
        % QPSK
        scale = 1 / sqrt(2);
        map1 = scale * [-1, 1]';
        Igray = [0, 1]';
        case 2
        % 16-QAM
        scale = 1 / sqrt(10);
        map1 = scale * [-3, -1, 1, 3]';
        Igray = [0, 1, 3, 2]';
        case 3
        % 64-QAM
        scale = 1 / sqrt(21);
        map1 = scale * [-7, -5, -3, -1, 1, 3, 5, 7]';
        Igray = [0, 2, 3, 1, 5, 7, 6, 4]';
    end % end switch
    
map1 = map1(Igray + 1);
bits = randi([0, 1], [nsym, 2 * nBitsDim]);
bits = num2str(bits);
Words = [bin2dec(bits(:, 1:nBitsDim)), bin2dec(bits(:, (nBitsDim + 1):end))];
sym = map1(Words(:, 1) + 1) + 1i * map1(Words(:, 2) + 1);
Wvar = repmat(1 ./ (4 * 10 .^ (EsN0Test / 10)), nsym, 1);
W = (randn(nsym, length(EsN0Test)) + 1i * randn(nsym, length(EsN0Test))) .*...
    sqrt(Wvar);
Y = repmat(sym, 1, length(EsN0Test)) + W;
for i = 1:length(EsN0Test)
    Icol = real(Y(:, i));
    Qcol = imag(Y(:, i));
    m = (1:(2^nBitsDim))';

    mrep = repmat(map1', nsym, 1);
    Irep = repmat(Icol, 1, length(m));
    Qrep = repmat(Qcol, 1, length(m));
    [~, a] = min(abs(Irep - mrep), [], 2);
    [~, b] = min(abs(Qrep - mrep), [], 2);
    RX_bits = [dec2bin(a - 1), dec2bin(b - 1)];
    RX_bits = RX_bits(:);
    bits = bits(:);
    BER(i, nBitsDim) = mean(RX_bits ~= bits);
    end % end for
end % end for


semilogy(EsN0Test, BER(:, 1), '-o', EsN0Test, BER(:, 2), '-go', EsN0Test, BER(:, 3));
grid on;
title('Bit Error Rate (BER) of 3 different modulation schemes', 'FontWeight', 'bold');
axis([min(EsN0Test), max(EsN0Test), 1e-4, 1]);
legend('QPSK', '16-QAM', '64-QAM');
xlabel('E_s/N_0 (dB)', 'FontWeight', 'bold');
ylabel('BER (%)', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');




