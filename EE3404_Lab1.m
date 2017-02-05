%% EE3404 Lab 1: Plotting and Filtering Signals
%%  q1
clear; clc; close all
%Read the file and get dimensions
fn = 'music_orig.wav'; %u can use other files, even ur own mp3~~
[x,fs]=audioread(fn);

x1 = x(:,1);
x2 = x(:,2);
dt = 1/fs;
t = 0:dt:(length(x1)*dt)-dt; %b/c t and x have to be same length

%% q2
figure
plot(t,x1,'k',t,x2,'c'); 
xlabel('Time (seconds)'); ylabel('Amplitude')

%% ft graphs of channels (just for my reference)
[M1,f] = dtft(x1,dt);
figure
subplot(2,1,1)
plot(f,M1)
xlabel('Frequency (Hz)')
title('Spectrum of music orig channel 1')
[M2,f] = dtft(x2,dt);
subplot(2,1,2)
plot(f,M2,'g')
xlabel('Frequency (Hz)')
title('Spectrum of music orig channel 2')

%% q3
clear p
p = audioplayer(x,fs);
% nbits = 16;
% sound(x,fs,nbits) << not used. just experimented
start = 1;
stop = p.SampleRate * 8; %sample rate property listed in audioplayer
% and 8 is number of seconds i want song to be played
play(p,[start,stop]);

%% q4 lpf 1st order
clear p
fc = 500; %cutoff frequency in Hz
% tau = 1/fc; %do i need this var? idk
SampleRate = 48000; %or fs
T = 1/SampleRate;
Wn = (2*fc)/SampleRate; %0.0208; %normalized fc Wn = (2*fc)/SampleRate for Hz
[b,a] = butter(1,Wn); %spits out coefficients of tf (no idea if this is right tho)
y = filter(b,a,x);
yrmsTgt = 0.5; %target rms value for scaled signal
yscale = yrmsTgt * sqrt(2); %no clue what to do
y = yscale*y; %i dunno much about this

p = audioplayer(y,fs);
start = 1;
stop = p.SampleRate * 8;
play(p,[start,stop])

%Yes I can hear the difference. It sounds mellowed out. In the unfiltered
%version you can hear the string being plucked, but not in this version.
%This filtered signal sounds less crisp. 
%Oooh sounds like someone is playing the guitar and not letting the sound
%ring~~*~

%% again this is for my reference
[M1,f1] = dtft(y(:,1),dt);
[M2,f2] = dtft(y(:,2),dt);
figure
subplot(2,1,1)
plot(f1,M1)
xlabel('Frequency (Hz)')
title('Spectrum of filtered music orig channel 1')
subplot(2,1,2)
plot(f2,M2)
xlabel('Frequency (Hz)')
title('Spectrum of filtered music orig channel 2')

%% q5 higher order filter 
%try cheby1 filter
clear a b p y; clc;
n = 4; %filter order
Rp = 0.5; %passband ripple in dB @_@
Wp = (2*fc)/SampleRate; %"passband edge frequency"/Digital frequency;
% Wp => 1 is Nyquist
[b,a] = cheby1(n,Rp,Wp);
y = filter(b,a,x);
yrmsTgt2 = 0.75;
yscale = yrmsTgt2 * sqrt(2);
y = yscale*y; 

p = audioplayer(y,fs);
start = 1; stop = p.SampleRate * 8;
play(p,[start,stop])

% This filtered signal sounds really muffled yet echoey. The "real" guitar
% sound is there but it's not a full sound.

%% q6
% cheby1 bandpass. passband between 500 Hz and 3 kHz. fs = SampleRate =
% 48000
clear a b p y n Rp Wp; clc;
n = 4;
Rp = 0.5;
flo = 500; %low freq
fhi = 3e3; %hi freq
w1 = (2*flo)/SampleRate; w2 = (2*fhi)/SampleRate;
Wp = [w1 w2];
[b,a] = cheby1(n,Rp,Wp);
y =filter(b,a,x);
y = yscale*y; %using yscale from first filter

p = audioplayer(y,fs);
start = 1; stop = p.SampleRate * 8;
play(p,[start,stop])

%The music with the bandpass filter sounds kind of hollow. I
%feel like since it's missing lower frequencies, it's not as smooth and
%mellow. Less timbre/harmonic content. 


%The natural frequency of the A string on a guitar is 440 Hz and this was 
%filtered out since its below the low cutoff freqency.

[M4,f4] = dtft(y,dt);
figure
plot(f4,M4)
xlabel('Frequency (Hz)')
title('Spectrum of Bandpass Filtered signal')
