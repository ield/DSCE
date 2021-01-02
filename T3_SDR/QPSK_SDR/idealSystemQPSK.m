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
synchroHeaderType = 3;
header     = [1 0 1 0 0 1 1 1 0 1];     % frame Sychronization (included txStream file)
% To have zero average, we have to use the synchro header for qpsk
switch synchroHeaderType                 % Necessary for time synchronization
    case 0
        synchroHeader=repmat([1 0], 1, 260);
    case 1
        synchroHeader=repmat([1 0 0 1], 1, 130);
    otherwise
        synchroHeader=repmat([1 1 0 0 1 0 0 1], 1, 65);
end
rng(0);
synchroAmbiguity = round(rand(1,50));  % Remove the synchronization phase ambiguity

% Symbols
M = 9;              % Samples per symbol=Oversamplig factor

% Modulation
fs=1e6; Ts =1/fs;            % Sampling frequency
% fc=2000;                       % Carrier frequency
phi=pi/8; deltafc=0.25;           % Carrier phase & frequency error

% Channel
SNR=50;               %dB

%% Baseband Tx - mapping + shaping filter
% frame
load txStream.mat           % Loading the bits

% Randomize the transmitted sequence: only the information bits.
txStream_rand = txStream(1+8+10+24:end-8);  % Select the message sequence
rand_seq = round(rand(1,length(txStream_rand)));    % Create the random sequence to perform the xor operation
txStream_rand = xor(txStream_rand, rand_seq');   % Perform the xor operation
txStream(1+8+10+24:end-8) = txStream_rand;      % Include the new values 

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
%R = fs/M   BW(hamming) = 2*2.6*R = 0.5 MHz
%           BW(cos) = 2*(1+alpha)*R = 0.333 MHz
% ps=hamming(M); % hamming pulse of width M
alpha = 0.5;
ps = srrc(0.5, alpha, 8);
msI=filter(ps,1,meI);              % convolve pulse shape with data
msQ=filter(ps,1,meQ);


%% Trasnmitting with Pluto
% It is used as a model test_Tx_carrier
Tt=20;      % Transmission duration (s)

fc=1e9;     % Carrier frequency (GHz)

% Select and configure the SDR
tx=sdrtx('Pluto');
tx.CenterFrequency = fc;
tx.Gain = -10;
tx.BasebandSampleRate = fs;

% Set up the signal to transmit
N=5000; % Transmit buffer length
% Real part of the signal
xtransI = zeros(1, N);
xtransI(1:length(msI)) = msI;
% Complex part of the signal
xtransQ = zeros(1, N);
xtransQ(1:length(msQ)) = msQ;
xtrans=xtransI' + 1i*xtransQ';

%Transmit
transmitRepeat(tx,xtrans);
% If we see xtrans, we see that the signal is the head, thesigna, the tail,
% and the end with zeros
%% Receiving with Pluto
pause(5);

rx=sdrrx('Pluto');
rx.CenterFrequency = fc;
rx.BasebandSampleRate = fs;
rx.GainSource = 'Manual';
rx.Gain = 0;

Ny=20000;
Nf = 1;     % Number of received buffers
y=zeros(Ny*Nf,1);

for ic=1:Nf
    [y(Ny*(ic-1)+(1:Ny)),valid,overflow]=rx();
    if overflow; disp(ic); end
end

release(rx);
release(tx);

%% Selecting the receiving section
% plot(real(y))
% Watching the plot of y we see that of a given random period, the starting
% point is in 7757, and the end in 7757+5000.
% However, it is also necessary to delete the last 0, so this is done just
% selecting the length of the sent message
% starting_point = 7757;


% This can also be done  by searching the starting pòint with some
% algorithm. It is looked for 25 0s in the derivative, and that point is
% where the chain of 0s starts. Therefore, it is looked for that chain of
% 0s.
starting_point = searchStartingPoint(y) + 175;

wanted_y = y(starting_point:starting_point + length(msI))/rms(y);
xfI = real(wanted_y);
xfQ = imag(wanted_y);

subplot(2, 1, 1); plot(xfI);
subplot(2, 1, 2); plot(xfQ);
% Note distortions:
% - The amplitude varies linearly with time.It increases and then
% decreases. The opposite happens in the imaginary plane.
% - There is no (almost) DC component but the part of the signal is not
% balanced in any of the channels, the amplitude in the -1 almost doubles
% the amplitude in the +1
% - The rest is almost the same
%% Receiver - Matching filtering + Timing Recovery + Phase recovery and slicer +
%             + Phase Ambiguity Correction
%matching filter
xI=filter(ps,1,xfI); 
xQ=filter(ps,1,xfQ); 

% timing recovery
figure(1); clf;
dseta = 1; Bn=25;
[ztI, ztQ] = timingRecoveryPLLQPSK(xI, xQ, M,dseta,Bn);

% Phase recovery and slicer
figure(2); clf;
[zI, zQ] = phaseRecoveryPLLQPSK(ztI, ztQ);

% Phase ambiguity correction
rxStream = ambiguityCorrection(zI, zQ, synchroAmbiguity);


%% Unframe and Message Visualization - frame synchronization + Presentation Layer
rh = xcorr(2*rxStream-1, 2*header-1); Rrh = rh( (length(rh)+1)/2:end );
[maxim, ind] = max(Rrh); headStart = ind(1); 

% Presentation layer
lengthMessagebin = num2str( rxStream(headStart+10:headStart+10+(24-1)) );
lengthMessage    = bin2dec(lengthMessagebin(:)');
%lengthMessage=44;

rxMsg = rxStream(headStart+10+24:headStart+10+24+8*(lengthMessage)-1 );
rxMsg = xor(rxMsg, rand_seq);   %Recover the message multiplying again with the xor
disp(msg2text(rxMsg')); %ascii to text
