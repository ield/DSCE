function [zI, zQ] = phaseRecoveryPLLQPSK(fs, ztI, ztQ)
%function [zsI, zsQ] = phaseRecoveryPLLQPSK(zI, ZQ),
%QPSK phase synchronization
%
%Input parameters
% zI, ZQ: unsynchronized phase symbols
% dseta: PLL parameter. Damping factor 
% Bn: PLL parameter. Noise BandWidth
%
%Output parameters
% zsI, zsQ: synchronized phase  symbols
%

%PLL constants
fs = 500; T=1/fs;  %sampling frequency
Bn = 25;            %noise Bandwidth(Hz)
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2
accError = 0;

zt= ztI + 1j*ztQ;
zs=zeros(size(zt)); zs(1)=zt(1);
phi = 0;
erroracc = zeros(size(zt));
M1=0;M2=0;
for ii=1:length(zt)-1
  %slicer
  % Step 1. Given the input, you decide in which quadrant of the
  % constellation is the signal
  zI(ii) = sign(real(zs(ii)));
  zQ(ii) = sign(imag(zs(ii)));
  z = zI(ii) + 1j*zQ(ii);
  
  %error
  % Step 2. You calculate the error between the estimate and the input (the
  % angle) see the notes on this part.
  error = atan2(imag(zs(ii)*conj(z)), real(zs(ii)*conj(z)));
    
  %PLL loop filter
%   v = error*k1 + error*k2 + M1; %loop filter
%   phi = phi + v;   %w0 can be cancelled
%   
%   M1 = error*k2 + M1;  %memories update
%   M2 = phi;
  
  %
  accError = accError + error;
  phi = phi + k1*error + k2*accError;
  
  %correction update
  % Step 3. Correct the next estimation with this step.
  zs(ii+1) = zt(ii+1)*exp(-1j*phi);                    %correction
  erroracc(ii)=phi;
end
% zI = real(zs);
% zQ = imag(zs);


% plot results
figure;
subplot(221); plot(ztI, ztQ,'b.')   % plot input constellation diagram
title('Input constellation diagram'); axis square; grid
subplot(222); plot(zI, zQ,'b.')   % plot output constellation diagram
title('Output constellation diagram'); axis square; grid
subplot(212); plot(erroracc)
title('Phase recovery PLL'); ylabel('errors')



