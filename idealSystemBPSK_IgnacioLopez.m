%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of Communication Systems and Equipment                           %
% Ignacio Esteban Lopez Delgado
% DCSE 2020/2021 MSTC-RF Track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clf;


%% Parameters
%Frame
header  = [1 0 1 0 0 1 1 1 0 1];

% Symbols and samples
M = 15;              % Samples per symbol=Oversamplig factor
ps = hamming(M);
% Modulation
fs=10000; Ts =1/fs;    % Sampling frequency
fc=2000;               % Carrier frequency
phi=0;                 % Carrier phase

% Channel
SNR= 50;               %dB


%%%%%%%%%%%%%%%%%%%%%%%%%% Tx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BaseBand Tx
load txStream.mat
L = length(txStream);

% Converting to bpsk
m = (txStream-0.5)*2;     % This way 0 are converted to -1 and 1 stay as 1.

me = zeros(L*M, 1);       % Transmitted sequence, with space for all data and zeros of the M
me(1:M:end) = m;          % M-1 zeros are inserted in between the data.

%shaping filter
ms = filter(ps, 1, me);    % Since it is a fir pulse, the denom = 1
%% Modulation
time = Ts*(0:length(ms)-1);
carrierTx = cos(2*pi*fc*time);
s = ms.*carrierTx';

%%%%%%%%%%%%%%%%%%%%%%%% Channel  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adding channel noise
%noise generation

% r = awgn(s, SNR);       % Function to add noise
pSignal = rms(s)^2;
pNoise = pSignal/(10^(SNR/10));                 % Noise power to achieve 50 dB ratio
noise=pNoise^0.5 * randn(1, length(s));      % Noise signal created scalated for snr = 50

r = s + noise';
% 10*log10(pSignal / rms(noise)^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Rx  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demodulation (modulation + filtering)
%cosine multiplication
carrierRx = cos(2*pi*fc*time + phi);
r2 = r .* carrierRx';

% LPF
order=30;                            % Order of the filter
fbe=[0 0.25 0.3 1]; damps=[1 1 0 0]; % LPF design. If the fs = 10000, 4000 is at .8. It is cur perfectly
b=firpm(order,fbe,damps);            % impulse response of LPF
xf = filter(b,1,r2);                 % LPF the demodulated signal


%% Matching filtering
x = filter(ps, 1, xf);    % Since it is a fir pulse, the denom = 1


%% Compressor
% 30 is taken looking at the plot of x. The filter takes some time to act.
% At 30 is the maximum (therefore the time at which the sample should be
% taken.
startSample = 30;
z = x(startSample:M:length(x));       

%% Decisor
% The signal z is a series of bits > and < 0. Therefore, the thresshold is
% placed in 0: + values are '1' and - values are '0'
rxStream = (sign(z)+1)/2;

%% frame synchronization. Unframe and message dispalay 
rxStream = rxStream(:)';
rh = xcorr(2*rxStream-1, 2*header-1); 
Rrh = rh( (length(rh)+1)/2:end );
[maxim, ind] = max(Rrh); headStart = ind(1); %headStart=9

%% Presentation layer
lengthMessagebin =  num2str( rxStream(headStart+10:headStart+10+(24-1)) );
lengthMessage    = bin2dec(lengthMessagebin(:)'); %lengthMessage=44;


rxMsg = rxStream(headStart+10+24:headStart+10+24+8*(lengthMessage)-1 );
disp(msg2text(rxMsg')); %ascii to text

