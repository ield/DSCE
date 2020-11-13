%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of Communication Systems and Equipment  
% file: idealSystemBPSKWithTimingRecovery
% System with BPSK modulation and including a PLL Timing Recovery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clf; clc


%% Parameters
%Frame
header  = [1 0 1 0 0 1 1 1 0 1];       %Frame synchronization (included in txStream file)
synchroHeader= repmat([1 1 1 0 0 0 0 1], 1, 65);           %Add this header for

% Symbols and samples
M = 15;                % Samples per symbol=Oversamplig factor

% Modulation
fs=10000; Ts =1/fs;    % Sampling frequency
fc=2000;               % Carrier frequency
phi=0;                 % Carrier phase

% Channel
SNR=50 ;               %dB


%% BaseBand Tx - mapping + shaping filter  
load txStream.mat
txStream = [synchroHeader'; txStream];

%mapping
m  = 2*txStream-1;    %mapping

% Generating the I and Q messages
mI = m(1:2:end);
mQ = m(2:2:end);

% figure;
% subplot(2, 1, 1)
% plot(mI, 'ro');
% subplot(2, 1, 2)
% plot(mQ, 'bo');

NSymbols = max(length(mI), length(mQ)); %number of symbols

%shaping filter
N=M*NSymbols;         %number of samples
meI = zeros(1, N);
meQ = zeros(1, N);

meI(1:M:end) = mI;     %expansion 
meQ(1:M:end) = mQ;
% figure;
% subplot(2, 1, 1)
% plot(meI, 'ro');
% subplot(2, 1, 2)
% plot(meQ, 'bo');

h  = hamming(M);     %filtering
msI = filter(h,1,meI);
msQ = filter(h,1,meQ);

% figure;
% subplot(2, 1, 1)
% plot(msI, 'r');
% subplot(2, 1, 2)
% plot(msQ, 'b');

%% Modulation
k = 0:N-1;                      % sampling index
deltaF = 0;                     % Frequency variation
phi = 0;                        % Phase variation

c=exp(1i*(2*pi*(fc+deltaF)*k*Ts + phi));    % carrier at freq fc with a small variation that must be recovered at the rx pll. The variation is introduced here. The rx will have frequency fc
s=real(c.*(msI + 1i*msQ));                  % modulation

% figure;
% plot(k*Ts, s, 'r');


%% Channel - adding noise  
n = randn(1,N);                   % noise generation
Ps = sum(s.^2)/N;                 % signal and noise powers
Pn = sum(n.^2)/N;
g = sqrt( Ps/(Pn*10^(SNR/10)) );  % noise gain
r = s + g*n;                      % corrupted signal

% figure;
% plot(k*Ts, r, 'r');


%% Demodulation (modulation + filtering)
k = 0:N-1;               % sampling index
c=exp(-1i*2*pi*fc*k*Ts);   % carrier at freq fc
r2I=real(c.*r);                        % modulation
r2Q=imag(c.*r); 

% figure;
% subplot(2, 1, 1)
% plot(r2I, 'r');
% subplot(2, 1, 2)
% plot(r2Q, 'b');

% LPF
order=30;                            % Order of the filter
fbe=[0 0.25 0.3 1]; damps=[1 1 0 0]; % LPF design
b=firpm(order,fbe,damps);            % impulse response of LPF
xfI = filter(b,1,r2I);                 % LPF the demodulated signal
xfQ = filter(b,1,r2Q);                 % LPF the demodulated signal

% figure;
% subplot(2, 1, 1)
% plot(xfI, 'r');
% subplot(2, 1, 2)
% plot(xfQ, 'b');

%% Rx  - Matching Filtering + Timing Recovery + Slicer/Decisor
%Matching Filtering
xI = filter(h,1,xfI);
xQ = filter(h,1,xfQ);

% figure;
% subplot(2, 1, 1)
% plot(xI, 'r');
% subplot(2, 1, 2)
% plot(xQ, 'b');

%Timing Recovery- Compressor
%z=x(30:M:end);   
dseta=1; Bn=25;
[zI, zQ] = timingRecoveryPLL(xI, xQ, M, dseta, Bn);
% figure;
% subplot(2, 1, 1)
% plot(zI, 'r');
% subplot(2, 1, 2)
% plot(zQ, 'b');

decI = zeros(size(zI));
decQ = zeros(size(zQ));
decI(zI<0) = 0;
decI(zI>=0) = 1;
decQ(zQ<0) = 0;
decQ(zQ>=0) = 1;

% figure;
% subplot(2, 1, 1)
% plot(decI, 'r');
% subplot(2, 1, 2)
% plot(decQ, 'b');

% Slicer/Decisor
rxStream = zeros(1, length(zI)+length(zQ));
rxStream(1:2:end) = decI;
rxStream(2:2:end) = decQ;
rxStream = [rxStream zeros(1, 100)];        % This is added so that all the message fits in the receiver: the length of the output of the pll is shorter than the message
figure;
plot(rxStream, 'r');


%% Unframe and Message Visualization - Frame synchronization + Presentation Layer   
%Unframe
rh = xcorr(2*rxStream-1, 2*header-1); Rrh = rh( (length(rh)+1)/2:end );
[maxim, ind] = max(Rrh); headStart = ind(1); %headStart=9

% Presentation layer
lengthMessagebin =  num2str( rxStream(headStart+10:headStart+10+(24-1)) );
lengthMessage    = bin2dec(lengthMessagebin(:)'); %lengthMessage=44;

rxMsg = rxStream(headStart+10+24:headStart+10+24+8*(lengthMessage)-1 );
disp(msg2text(rxMsg')); %ascii to text

