% SNRcalculations.m: using linear filters to improve SNR
% Exercise done by Ignacio Esteban Lopez Delgado

time=3;                   % time
Ts=1/48000;               % sampling interval: Since fs = 48000 Hz

% Cut frequency on 3000 and 4000 hz. The frequency band is determined by Ts.
% The maximum frequency is fs/2. Therefore, fs/2 corresponds to 1 when
% designing the filter. 

%----------------------------------
freqs=[0 2*2999*Ts 2*3000*Ts 2*4000*Ts 2*4001*Ts 1];    % filter design, bandlimited
amps=[0 0 1 1 0 0];               % ...between 3K and 4K: WRONG
%----------------------------------
b=firpm(100,freqs,amps);          % BP filter
n=0.25*randn(1,time/Ts);          % generate white noise signal
%**********************************


% Calculate the PSD of n. Insert code.
% Plot it and check with the noise variance
% As a starting point use the code for plotting
% spectra you'll find below


%**********************************
x=filter(b,1,2*randn(1,time/Ts)); % do the filtering
y=filter(b,1,x+n);                % (a) filter the signal+noise
yx=filter(b,1,x);                 % or (b) filter signal 
yn=filter(b,1,n);                 % ...and noise separately
z=yx+yn;                          % add them
diffzy=max(abs(z-y))              % and make sure y = z

%Plot noise spectra
N=length(n);                         % length of the signal x
t=Ts*(1:N);                          % define time vector
ssf=(-N/2:N/2-1)/(Ts*N);             % frequency vector

fn=fftshift(fft(n(1:N)));
figure(4), subplot(2,1,1), plot(ssf,abs(fn))
xlabel('Noise spectra before filtering');

fnFiltered = fftshift(fft(yn(1:N)));
subplot(2,1,2), plot(ssf,abs(fnFiltered))
xlabel('Noise spectra after filtering');

%**********************************
% Obtain the theoretical values of SNR
% and check them with the estimates provided
% by the two following lines of code
%**********************************
snrinp=pow(x)/pow(n)              % SNR at input
snrout=pow(yx)/pow(yn)            % SNR at output

% Theoretical values: var(signal) / var(noise)
snrInpTh = var(x) / var(n);
snrOutTh = var(yx) / var(yn);
% The results observed show that the theoretical value is similar to the
%   ones obtaince

%**********************************

%% Repeat for a 2 kHz sinusoid with the same power than the noise
clear;
time=3;                   % time
Ts=1/48000;               % sampling interval: Since fs = 48000 Hz
f = 2000;
t=Ts:Ts:time;

% Filter design
freqs=[0 2*1950*Ts 2*1951*Ts 2*2050*Ts 2*2051*Ts 1];    % filter design, bandlimited
amps=[0 0 1 1 0 0];               % around 2K
%----------------------------------
b=firpm(100,freqs,amps);          % BP filter
n=0.25*randn(1,time/Ts);          % generate white noise signal

% From the noise power it is determined the power of the signal
x =(var(n)^2 / 2) * sin(2*pi*f*t); % With this amplitude the power is the same as the noise

%**********************************

% x=filter(b,1,2*randn(1,time/Ts)); % do the filtering
y=filter(b,1,x+n);                % (a) filter the signal+noise
yx=filter(b,1,x);                 % or (b) filter signal 
yn=filter(b,1,n);                 % ...and noise separately
z=yx+yn;                          % add them
diffzy=max(abs(z-y))              % and make sure y = z

%Plot noise spectra
N=length(n);                         % length of the signal x
t=Ts*(1:N);                          % define time vector
ssf=(-N/2:N/2-1)/(Ts*N);             % frequency vector

fn=fftshift(fft(n(1:N)));
figure(4), subplot(2,1,1), plot(ssf,abs(fn))
xlabel('Noise spectra before filtering');

fnFiltered = fftshift(fft(yn(1:N)));
subplot(2,1,2), plot(ssf,abs(fnFiltered))
xlabel('Noise spectra after filtering');

%**********************************
% Obtain the theoretical values of SNR
% and check them with the estimates provided
% by the two following lines of code
%**********************************
snrinp=pow(x)/pow(n)              % SNR at input
snrout=pow(yx)/pow(yn)            % SNR at output

% Theoretical values: var(signal) / var(noise)
snrInpTh = var(x) / var(n);
snrOutTh = var(yx) / var(yn);
% The results observed show that the theoretical value is similar to the
%   ones obtaince

% Plot the signal
fx=fftshift(fft(x(1:N)+n(1:N)));
figure(5), subplot(2,1,1), plot(ssf,abs(fx))
xlabel('magnitude spectrum of signal + noise')

fy=fftshift(fft(y(1:N)));
subplot(2,1,2), plot(ssf,abs(fy))
xlabel('magnitude spectrum after filtering')

%% How to plot spectra
N=length(x);                         % length of the signal x
t=Ts*(1:N);                          % define time vector
ssf=(-N/2:N/2-1)/(Ts*N);             % frequency vector

fx=fftshift(fft(x(1:N)+n(1:N)));
figure(5), subplot(2,1,1), plot(ssf,abs(fx))
xlabel('magnitude spectrum of signal + noise')

fy=fftshift(fft(y(1:N)));
subplot(2,1,2), plot(ssf,abs(fy))
xlabel('magnitude spectrum after filtering')