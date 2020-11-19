fs=10e5;
fc=2.4e9;

tx=sdrtx('Pluto');
tx.CenterFrequency = fc;
tx.Gain = -10;
tx.BasebandSampleRate = fs;

Nf=40;          % Number of receive buffers
Nx=100;         % Period of transmit sinusoid (samples)

N=5000;
x=exp(1i*2*pi/Nx*(1:N)');   % They fit exactly 50p eriods in the buffer: so the last sample is consistant with the first.
x=complex(x);

transmitRepeat(tx,x);

pause(5);               % Let the transmitter stabilize

rx=sdrrx('Pluto');
rx.CenterFrequency = fc;
rx.BasebandSampleRate = fs;
rx.GainSource = 'Manual';
rx.Gain = 0;

Ny=20000;
y=zeros(Ny*Nf,1);

for ic=1:Nf
    [y(Ny*(ic-1)+(1:Ny)),valid,overflow]=rx();
    if overflow; disp(ic); end
end

release(rx);
release(tx);

%% Plot the received signal
load rx.mat;        % When working from home.

figure(1);
subplot (121); plot(y); grid on; axis('square'); axis('equal');
xlabel('in-phase'); ylabel('quadrature');
subplot(222); plot(real(y)); title('in-phase');
subplot(224); plot(imag(y)); title('quadrature');

% Sometimes Pluto fails.
%   There can be discontinuities in amplitude, so there are different
%       circles.
%   The amplitude is not constant so there is some phase error (the circle
%       gets thicker and thinner in some places).
%   The amplitde is not the same in the inphase and quadrature planes
%   The phase difference is not perfect 90º, 
%   Try to have constant amplitude: as if using an AGC

%% Amplitude correction
% We are going to create s'1(t) = A(t) sin (wc t) as the Hilbert transform
% of s1(t)
lenY = length(y);

s1 = real(y);
s1_prime = imag(hilbert(s1));       % According to matlab help it is the imaginary part the interesting one
A_1 = sqrt(s1.^2 + s1_prime.^2);      % A(t), variation of amplitude with time. The error that must be corrected
s1_corrected = s1 ./ A_1;

% figure;
% subplot(2, 1, 1);
% plot(s1);
% hold on;
% plot(s1_prime);
% subplot(2, 1, 2)
% plot(s1_corrected);

% Since we saw that there is also need for phase correction in channel 2,

s2 = imag(y);
s2_prime = imag(hilbert(s2));       % According to matlab help it is the imaginary part the interesting one
A_2 = sqrt(s2.^2 + s2_prime.^2);      % A(t), variation of amplitude with time. The error that must be corrected
s2_corrected = s2 ./ A_2;

% Now we have that s1_corrected is m1(t), the original message. Now we have
% to correct the phase error

%% Phase correction
% We have now that m2(t) = sin(wct-a) = while mi(t) = cos(wct) = s1_corrected. Therefore,
% if we do <m1, m2> and we do the average of the result, we have that the
% average is -0.5sin(alpha)
m1 = s1_corrected;
m2 = s2_corrected;
average_scalar_product = mean(m1.*m2);
alpha = asin(2*average_scalar_product);

% Also knwoing that m2 = (s2(t) + s1(t)sina) / cosa
s2_corrected = m2*cos(alpha)-m1*sin(alpha);

% Now we have all the signals corrected, so we have to plot them
%% Plot and anayze results
y_corrected = s1_corrected + 1i*s2_corrected;


figure(2);
subplot (121); plot(y_corrected); grid on; axis('square'); axis('equal');
xlabel('in-phase'); ylabel('quadrature');
subplot(222); plot(real(y_corrected)); title('in-phase');
subplot(224); plot(imag(y_corrected)); title('quadrature');



