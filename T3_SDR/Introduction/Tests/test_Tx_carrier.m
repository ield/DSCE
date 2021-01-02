%
% Signal menu:
% 0: Carrier only
% 1: Cosine with period Nx samples
% 2: Two complex sinusoids with period Nx samples and different amplitudes
% 3: Pseudo-random QPSK with length N samples, 10 samples per symbol
%    and Hanning pulse shaping.
%

Sselect=0;  % Signal type selection
Tt=20;      % Transmission duration (s)

fs=1e6;
fc=1e9;

tx=sdrtx('Pluto');
tx.CenterFrequency = fc;
tx.Gain = -10;
tx.BasebandSampleRate = fs;

% Try to reduce the frequency offset as
% measured by the analyser. Units are ppm ref to carrier


% When it is selected 1 GHz, the peak is at 999.98898 MHz Therefore, the
% frequency is corrected to We want a correction of 1000000-9999889.98
% After a first approach the result is measured in 999.99438
correction_frequency = fc - 0.999982930e9;
% tx.FrequencyCorrection = (correction_frequency)/fc*1e6    % Carrier correction in ppm
% It is corrected to 0.99999991 GHz
new_error = (0.99999991/fc)*1e6;
fprintf('Old error = %f ppm. New error = %f ppm\n', tx.FrequencyCorrection, new_error);
% The result is not very accurate, it moves from test to test, as the
% temperature changes.

N=5000; % Transmit buffer length
switch Sselect
    case 0
        x=ones(N,1);
    case 1
        Nx=10;
        x=0.5*(exp(1i*2*pi/Nx*(1:N)')+exp(-1i*2*pi/Nx*(1:N)'));
        % f/fs = 1/10. Therefore, f = fs/10 = 1e5 Hz (expected)
    case 2
        Nx=10;
        x=0.5*(exp(1i*2*pi/Nx*(1:N)')+exp(-1i*2*pi/Nx*(1:N)')/sqrt(10));
    case 3
        x=base_QPSK(N,10,0).';
end
x=complex(x);
transmitRepeat(tx,x);
pause(Tt)
release(tx);

fread = 999.9811e6;       % Carrier frequency read in the analyzer (with spam of 10 kHz
fprintf('The carrier frequency is %f MHz\n', fread/1e6);
% See case 1 to determine the expected carrier
fc_saw = 999.88133e6;
fprintf('Expected fcos = %f kHz. Seen fcos = %f kHz\n', fs/10/1e3, (fread-fc_saw)/1e3);
% Expectre bw
% bw = 2/Tb for a quared signal. However, this is not a pure squared signal. 
% There the bw is scales by a bw factor of 1.8


