fs=6e5;         % Sampling frequency: it determines the bw
fc=1e9;         % Carrier frequency. some frequencies are forbidden because they have other uses.

% Open transmitter and configure it

tx=sdrtx('Pluto');

% Set transmitter parameters:
% Carrier frequency
% Sampling rate
% Transmitter gain. -10 is adequate for connecting the spectrum analyser

tx.CenterFrequency = fc;
tx.BasebandSampleRate = fs;
tx.Gain = -10;  %dB

N = 5000;         % Length of transmission buffer
Nb= 1000;         % Number of buffers to be transmitted
L = Nb*N;       % Length of transmitted signal
x=ones(L,1);      % Signal

for bufferCount = 1:Nb
    buffer = x((bufferCount-1)*N+(1:N));  % Select signal block to be transmitted
    buffer=complex(buffer);               % PLUTO requires complex (I/Q) signals
    tx(buffer);
end
release(tx);        % Close transmitter

