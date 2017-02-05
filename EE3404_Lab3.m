%% EE3404 Lab 3: QAM Modulator Simulation | Aimee Nogoy | akn264

%% Q1: generate random bits
clear; clc;
nbits = 2^17; %number of bits to process

%M = 16; %size of signal constellation
% Nsym = 2^15. Nb/Nsym = 4 bits/sym. hence M =16..?_?
%k = log2(M); %number of bits in each constellation pt

bits = randi([0 1],nbits,1)';

%% Q2: QAM mapping
%make sure bits are Gray coded
A = 1/sqrt(10); % for normalized energy
sym = zeros(1,nbits/4);
n = 1; %index var

for p = 1:length(sym)
    if bits(n)==0 %0000
        if bits(n+1)==0
            if bits(n+2)==0
                if bits(n+3)==0
                    sym(p) = -A + 3*A*1i;
                else %0001
                    sym(p) = -A + A*1i;
                end
            else 
                if bits(n+3)==0 %0010
                    sym(p) = -A - 3*A*1i;
                else %0011
                    sym(p) = -A - A*1i;
                end
            end
        else 
            if bits(n+2)==0 
                if bits(n+3)==0 %0100
                    sym(p) = A + 3*A*1i;
                else %0101
                    sym(p) = A + A*1i;
                end
            else
                if bits(n+3)==0 %0110
                    sym(p) = A - 3*A*1i;
                else %0111
                    sym(p) = A - A*1i;
                end
            end
        end        
    else 
        if bits(n+1)==0
            if bits(n+2)==0
                if bits(n+3)==0 %1000
                    sym(p) = -3*A + 3*A*1i;
                else %1001
                    sym(p) = -3*A + A*1i;
                end
            else 
                if bits(n+3)==0 %1010
                    sym(p) = -3*A - 3*A*1i;
                else %1011
                    sym(p) = -3*A - A*1i;
                end
            end
        else
            if bits(n+2)==0
                if bits(n+3)==0 %1100
                    sym(p) = 3*A + 3*A*1i;
                else %1101
                    sym(p) = 3*A + A*1i;
                end
            else 
                if bits(n+3)==0 %1110
                    sym(p) = 3*A - 3*A*1i;
                else %1111
                    sym(p) = 3*A - A*1i;
                end
            end
        end
    end
    n = n + 4; % move on to next 4 bits for next p
end
%tic; toc

%% Q3: Digital pulse shaping             

%Digital LPF and upsample
Nov = 4; %oversampling ratio
nfilt = 60; % filter length
wp = 1/Nov; %filter digital bandwidth
bfilt = Nov*fir1(nfilt,wp);% create filter

sym1 = upsample(sym, Nov); % insert zeros
sym1 = filter(bfilt,1,sym1); % filter
delay = (nfilt/2) + 1;
sym1 = sym1(delay:end);

%Plotting
np = 30; % number samples to plot
t0 = (0:length(sym)-1); 
t1 = (0:length(sym1)-1)/Nov;

figure
plot(t1(1:np*Nov),sym1(1:np*Nov),'-', t0(1:np), sym(1:np), 'o');
set(gca,'FontSize',16);
title('Original and Upsampled Signals (Re)')
legend('Interpolated (sym1)', 'Discrete (sym)');
axis([0, np, -3, 3])

figure
plot(t1(1:np*Nov),imag(sym1(1:np*Nov)),'-', t0(1:np), imag(sym(1:np)), 'o');
set(gca,'FontSize',16);
title('Original and Upsampled Signals (Imag)')
legend('Interpolated (sym1)', 'Discrete (sym)');
axis([0, np, -3, 3])

%% Q4: Plot the impulse response
close all
fsym = 10e6; 
Tsym = 1/fsym; 
Tsamp = Tsym/Nov;

n = 0:length(bfilt)-1;
tn = (n - nfilt/2)*Tsamp;
y = sinc(tn*fsym); %ideal sinc interpolation filter
plot(tn,bfilt,tn,y,'--')
grid on
xlabel('Time (s)'); ylabel('Amplitude')
legend('b[n]','ideal sinc')
%The graphs are similar shape, but the non ideal filter has a time delay.
% ^^ wrong

% It's like an estimate of the ideal filter
%% Q5: Plotting the frequency response
%compute the freq resp
nw = 128; %num points in the plot
[Hfilt,wfilt] = freqz(bfilt,1,nw); %get the freq resp
Hpow = 20*log10(abs(Hfilt)); %compute the mag squared in dB

plot(wfilt/(2*pi),Hpow)
grid on
axis([0,0.5,-90,20])
ylabel('Magnitude Squared (dB)'); xlabel('Frequency (MHz)')
title('Frequency Response of FIR Filter')                

%% Q7: Measure the PSD
%compute the PSD with the welch method
nfft = 1024;
utx = sym1;
[Su,fps] = pwelch(utx,[],[],nfft,(1/Tsamp),'centered');
plot(fps,10*log10(Su))
title('Power Spectral Density of U(t)')
xlabel('Frequency (Hz)');ylabel('PSD in dBm/Hz')
grid  on

%% Q8: Measure the excess bandwidth
%optional

%% Q9: Recover the transmitted samples
% RX the symbols
clc
doff = 1; %offset from 1 to Nov
urx = utx(doff:Nov:end); %downsample the TX samples
plot(urx,'o','Color',...
    [0,0.7,0.9],'MarkerFaceColor',[0,0.7,0.9]');
grid on
title('16-QAM TX/RX Symbols')




