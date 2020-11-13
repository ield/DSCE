%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of Communication Systems and Equipment
% file: idealSystemQPSK.m
% System with QPSK modulation  and including: 
%   - PLL timing synchronization: timingRecoveryPLLQPSK.m
%   - PLL phase synchronization: phaseRecoveryPLLQPSK.m
%   - Disambiguation: ambiguityCorrection.m
%   - input data. txStream.mat
%   - msg2text. ascii codes to text 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clf; clc;

%% Parameters
%frame
header     = [1 0 1 0 0 1 1 1 0 1];    % frame Sychronization (included txStream file)
synchroHeader= zeros(1,500);           % Necessary for time synchronization
synchroHeader(1:2:end)=ones(1,250);
synchroHeader=repmat([1 1 0 0 1 0 0 1], 1, 65);
rng(0); synchroAmbiguity = round(rand(1,50));  % Remove the synchronization phase ambiguity

% Symbols
M = 15;              % Samples per symbol=Oversamplig factor

% Modulation
fs=10000; Ts =1/fs;            % Sampling frequency
fc=2000;                       % Carrier frequency
phi=pi/2; deltafc=1.2;           % Carrier phase & frequency error

% Channel
SNR=50;               %dB

%% Baseband Tx - mapping + shaping filter
% frame
load txStream.mat           % Loading the bits
txStream = [synchroHeader'; synchroAmbiguity'; txStream; zeros(100,1)]; %zeros are added                                                    %compensate the                                                      %filters

% QPSK mapping
%m = 2*txStream - 1;
mI = 2*txStream(1:2:end)-1;
mQ = 2*txStream(2:2:end)-1;
NSymbols = length(mI);   %Number of Symbols


%Shaping filtering - oversampling + pulse shaping
% Oversampling by M
meI=zeros(1, M*NSymbols);      
meI(1:M:end)= mI;        % Oversampling the message vector by M
meQ=zeros(1, M*NSymbols);      
meQ(1:M:end)= mQ;        % Oversampling the message vector by M

% Pulse shaping
ps=hamming(M);                     % hamming pulse of width M
msI=filter(ps,1,meI);              % convolve pulse shape with data
msQ=filter(ps,1,meQ);



%% Modulator
k = 0:M*NSymbols-1;                         % sampling index = length(ms)
c=exp( 1i*(2*pi*(fc+deltafc)*k*Ts+phi) );   % carrier at freq fc
s=real( c.*(msI + 1j*msQ) );     % modulation


%% Channel - adding noise
Ps    = sum(s.^2)/length(s);    % Power of the  transmitted signal

% generate noise signal
sigma = sqrt(Ps/(10^(SNR/10))); % Calculate the sigma of the noise
n = sigma*randn(size(s));
r = s + n;                  % received noisy signal


%% Demodulation - modulation + filering
c2  = exp( -1j*(2*pi*fc*k*Ts) );     % synchronized cosine for mixing
r2  = r.*c2;                         % demod received signal
r2I = real(r2);
r2Q = imag(r2);

% LPF
order=30;                              % Order of the filter
fbe = [0 0.25 0.3 1]; damps=[1 1 0 0]; % LPF design
b   = firpm(order,fbe,damps);          % impulse response of LPF
xfI = filter(b,1,r2I);                 % LPF the demodulated signal
xfQ = filter(b,1,r2Q);


%% Receiver - Matching filtering + Timing Recovery + Phase recovery and slicer +
%             + Phase Ambiguity Correction
%matching filter
xI=filter(ps,1,xfI); 
xQ=filter(ps,1,xfQ); 

% timing recovery
dseta = 1; Bn=25;
[ztI, ztQ] = timingRecoveryPLLQPSK(xI, xQ, M,dseta,Bn, fs);
%ztI = xI(30:M:end);   %by hand
%ztQ = xQ(30:M:end);

% Phase recovery and slicer
[zI, zQ] = phaseRecoveryPLLQPSK(fs, ztI, ztQ);

% Phase ambiguity correction
rxStream = ambiguityCorrection(zI, zQ, synchroAmbiguity);

% zI(zI == -1) = 0;
% zQ(zQ == -1) = 0;
% non-definitive code
% zI  = zeros( size(ztI) ); zI(ztI>0) = 1;
% zQ  = zeros( size(ztQ) ); zQ(ztQ>0) = 1;
% rxStream=zeros( 1, 2*length(zI) ); 
% rxStream(1:2:end) = zI;
% rxStream(2:2:end) = zQ;
% rxStream = [rxStream zeros(1, 1000)];






%% Unframe and Message Visualization - frame synchronization + Presentation Layer
rh = xcorr(2*rxStream-1, 2*header-1); Rrh = rh( (length(rh)+1)/2:end );
[maxim, ind] = max(Rrh); headStart = ind(1); 

% Presentation layer
lengthMessagebin = num2str( rxStream(headStart+10:headStart+10+(24-1)) );
lengthMessage    = bin2dec(lengthMessagebin(:)');
%lengthMessage=44;

rxMsg = rxStream(headStart+10+24:headStart+10+24+8*(lengthMessage)-1 );
disp(msg2text(rxMsg')); %ascii to text
