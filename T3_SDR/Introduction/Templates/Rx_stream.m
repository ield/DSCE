fs=10e5;
fc=2.4e9;
% Open transmitter and configure it

rx=sdrrx('Pluto');

% Set receiver parameters:
% Carrier frequency
% Sampling rate
% AGC
% Receiver gain in dB

rx.CenterFrequency = fc;
rx.BasebandSampleRate = fs;
% The SDR has a very wide range so it has a agc to control the amplitude of
% the received signal.
rx.GainSource = 'Manual';       % Disable AGC for manual setting of gain.
rx.Gain = 0;

Ny=20000;       % Length of receive buffer
Nb=5;           % Number of signal buffers to be received
y=zeros(Ny*Nb,1); % Received signal

for ic=1:Nb
    [y(Ny*(ic-1)+(1:Ny)),valid,overflow]=rx();  % Store buffers at proper: it captures as many samples as there are in the vector
                                                % position in received signal
    if overflow; disp(ic); end
end

release(rx);

% If you connect your PLUTO to another one configured as transmitter
% you can plot the received signal as follows

figure(1);
subplot(211); plot(real(y));
subplot(212); plot(imag(y));
