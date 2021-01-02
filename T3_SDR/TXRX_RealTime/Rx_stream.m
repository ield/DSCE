fs=6e5;
fc=2.2e9;
% Open transmitter and configure it

rx=sdrrx('Pluto');

% Set receiver parameters:
% Carrier frequency
% Sampling rate
% AGC
% Receiver gain in dB

rx.CenterFrequency = fc;
rx.BasebandSampleRate = fs;
rx.GainSource = 'Manual';       % Disable AGC for manual setting of gain
rx.Gain = 55;                    % The gain goes from -4 to 71 dB

Ny=20000;       % Length of receive buffer
Nb=5;           % Number of signal buffers to be received
y=zeros(Ny*Nb,1); % Received signal

for ic=1:Nb
    [y(Ny*(ic-1)+(1:Ny)),valid,overflow]=rx();  % Store buffers at proper
                                                % position in received signal
    if overflow; disp(ic); end
end

release(rx);

% If you connect your PLUTO to another one configured as transmitter
% you can plot the received signal as follows
y_agc = agc(y);

figure(1);
subplot(211); plot(real(y_agc));
subplot(212); plot(imag(y_agc));

delta_f0 = 13e3;
y_carr = carrierRecovery(y_agc, fc, 1/fs, delta_f0);
%% Test agc
t = [1:10000]/1e7;
x = 2000*cos(2*pi*500e3*t) + 1i*200*sin(2*pi*500e3*t);
s = agc(x);

figure;
subplot(2, 1, 1)
plot(x)
subplot(2, 1, 2)
plot(real(s))
%%
% The transmitter is 4
[~,locs] = findpeaks(real(y));
diff_sample = diff(locs);
% mean(diff_sample);
freq_diff = fs/mean(diff_sample);

% To know if the frequency shift is positive or negative it is plotted the
% Lisajoux. If the shift is positive (counter clock wise) the it means that
% the frequency of the receiver (mine) is lower than the transmitter.
% Therefore, the transmitter must reduce (negative) the ppm
begin = 500;
samples_freq = y(begin:begin+round(0.8*mean(diff_sample)));
figure;
plot(samples_freq(1:end-1));
hold on;
plot(samples_freq(end), 'o');

% My first result is
ppm_shift = freq_diff/fc * 1e6;
fprintf('The shift to be done is %f ppm', ppm_shift);

