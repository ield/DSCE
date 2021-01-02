fs=6e5;
fc=2.2e9;

% Open transmitter and configure it

tx=sdrtx('Pluto');

% Set transmitter parameters:
% Carrier frequency
% Sampling rate
% Transmitter gain. -10 is adequate for connecting the spectrum analyser

tx.CenterFrequency = fc;
tx.BasebandSampleRate = fs;
tx.Gain = 0;        % Gain goes from -89.75 dB to 0 dB
% tx.FrequencyCorrection = 5.45;     %Nacho
tx.FrequencyCorrection = 2;     %Diego
N = 5000;           % Length of transmission buffer
x = ones(N,1);      % Signal to be transmitted
Tt = 20;            % Transmission duration (s)

x=complex(x);       % PLUTO requires complex (I/Q) signals

transmitRepeat(tx,x);  % Transmission buffer is repeated until release

pause(Tt)           % Release after Tt seconds
release(tx);
