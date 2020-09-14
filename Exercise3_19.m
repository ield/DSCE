% speccos.m plot the spectrum of a cosine wave
f=100; phi=0;                % specify frequency and phase
time=2;                     % length of time
Ts=1/1000;                   % time interval between samples
t=Ts:Ts:time;               % create a time vector
x=cos(2*pi*f*t+phi);        % create cos wave
plotspec(x,Ts)              % draw waveform and spectrum

% 3.19.b
f1 = 100;
f2 = 150;
x=cos(2*pi*f1*t) + cos(2*pi*f2*t);
plotspec(x,Ts)

%% 3.19.c and d
% filternoise.m filter a noisy signal three ways

time=3;                          % length of time
Ts=1/10000;                      % time interval between samples
x=randn(1,time/Ts);              % generate noise signal
figure(1),plotspec(x,Ts)         % draw spectrum of input

% Cut frequency on 100 and 300 hz. The frequency band is determined by Ts.
% The maximum frequency is fs/2. Therefore, fs/2 corresponds to 1 when
% designing the filter. 
fs = 1/Ts;
freqs=[0 2*99/fs 2*100/fs 2*300/fs 2*301/fs 1];
amps=[0 0 1 1 0 0];
b=firpm(1000,freqs,amps);         % specify the LP filter
ylp=filter(b,1,x);               % do the filtering
figure(2),plotspec(ylp,Ts)       % plot the output spectrum


