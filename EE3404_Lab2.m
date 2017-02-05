%% EE3404 Lab 2: Power Spectrum and Frequency Modulation | Aimee Nogoy | akn264

%% Q1 Generate a random binary sequence
clear; close all; clc
fs0 = 10; %inital sample rate in MHz
nt0 = 1e3; %number time samples at the frequency fs0
%random binary sequence of nt0 samples. Each sample is +/-1

x0 = -1 + (1+1)*round(rand(1,nt0)); %row vec 

% or x0 = randi([0,1],[nt0,1]);
% x0(x0 == 0) = -1;   <-- finds zeros and sets them equal to -1
% sundeep way: 2*randi([0 1],nt0,1)-1

%% Q2 Upsample the signal
clc;
M = 8; %upconversion factor
nfilt = 80; %filter length
fs1 = M*fs0; %sample rate at the higher frequency in MHz

% design filter using FIR filter design
wp = 1/M; %sets the filter digital bandwidth. 1/M gives a cutoff of pi/M
bfilt = M*fir1(nfilt,wp);

% upsample
x1 = upsample(x0,M); %insert zeros(M-1 number of zeros between each sample)
x1 = filter(bfilt,1,x1); %filter

% shift signal to compensate for the delay
dly = 41; %approx half of the filter length 
x1 = x1(dly:end);


%% Q3 Plotting the interpolation
clc; close all
np = 20; %number samples to plot
dt0 = 1/fs0; 
dt1 = 1/fs1;

t0 = 0:dt0:((length(x0)*...
    dt0)-dt0); %time in microseconds
t1 = 0:dt1:((length(x1)*...
    dt1)-dt1);

plot(t1(1:np*M),x1(1:np*M),'-',t0(1:np),x0(1:np),'o');
set(gca,'FontSize',10);
title('Original and Upsampled Signals')
legend('Interpolated','Discrete');
xlabel('Time (\mus)');

%% Q4 Plot the PSD
clc
nfft = 512; %num freq points
stinx = fs1;
[Px,f] = pwelch(x1,[],[],nfft,stinx,'centered');

plot(f,10*log10(Px),'Color',[0,0.7,0.9])
hold on

% setting up the vertical lines
neg_f0 = [-fs0,-fs0];
pos_f0 = [fs0,fs0];
y1 = get(gca,'ylim'); %getting ylim for the vertical lines
plot(neg_f0,y1,'b','LineWidth',1.2)
plot(pos_f0,y1,'b','LineWidth',1.2)

title('Power Spectral Density Sx1(f)')
ylabel('PSD(dBm/MHz)'); xlabel('Frequency (MHz)');
grid
axis([-20,20,-75,0])

%% Q5 Upconversion
clc; close all
format long; clear Px f
fc = 20; %carrier frequency in MHz

xpt = x1 .* cos(2.*pi.*fc.*t1);

[Px,f] = pwelch(xpt,[],[],nfft,fc,'centered');
plot(f,10*log10(Px),'k')

title('Power Spectral Density Sxp(f)')
ylabel('PSD (dBm/MHz)'); xlabel('Frequency (MHz)');
grid on

%% Q6 Downconversion
clc; close all
% Downconvert
x2 = 2.*xpt.*exp(-2.*pi.*1i.*fc.*t1); %not sure what to put for t

% Downsample
bfiltd = bfilt/M;
x2 = filter(bfiltd,1,x2); %filter
x2 = x2(dly:M:end);

% Plot the first 100 samples
np = 100;
t0 = (0:np-1)/fs0;
plot(t0,[real(x2(1:np)); x0(1:np)]);
grid on;
set(gca,'FontSize',12);
legend('RX','TX');
title('Transmit and Receive Signals')
xlabel('Time (\mus)');

