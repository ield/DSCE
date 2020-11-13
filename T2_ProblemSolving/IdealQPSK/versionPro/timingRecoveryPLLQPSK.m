function [zI, zQ]=timingRecoveryPLL(xI, xQ, M, dseta, Bn, fs)
%function xs=timingRecoveryPLL(x),
%Timing synchronization for the ideal system
%
%Input parameters
% x: Matching filter output
% M: Oversampling factor or samples per Symbol
% dseta: PLL parameter. Damping factor 
% Bn: PLL parameter. Noise BandWidth
%
%Output parameters
% zI, zQ: Decimated output samples (one per symbol)
%

x = sqrt( xI.^2 + xQ.^2 ); 

%PLL constants
fs = 500; T=1/fs;  %sampling frequency
%Bn = 25;            %noise Bandwidth(Hz)
kp = 1;            %phase gain
k0 = 1;            %VCO gain
%dseta = 1;         %damping factor

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2

% clock recovery algorithm
lside=50;
tnow=lside*M+1;      % initialize variables
tau=0; 
xs=zeros(size(x));  
tausave=zeros(size(x)); tausave(1)=tau; %for plotting
i=0;   %loop index
%mu=0.01;                           % gradient algorithm stepsize
delta=0.1;                          % time for derivative
M1=0; M2=0;                         % PLL memories
accError = 0;
while tnow<length(x)-lside         % run iteration
  i=i+1;
  %TED
  xs(i)    = interpsinc(x,tnow+tau,lside);       % interp value at tnow+tau
  zI(i)    = interpsinc(xI,tnow+tau,lside);
  zQ(i)    = interpsinc(xQ,tnow+tau,lside);
  x_deltap = interpsinc(x,tnow+tau+delta, lside);% value to right
  x_deltam = interpsinc(x,tnow+tau-delta,lside); % value to left
  dx       = x_deltap-x_deltam;                  % numerical derivative
  
  %tau=tau+mu*dx*xs(i);   % gradient alg. update: OP
  
  error = dx*xs(i);                 %PLL alg. update:OP 
  
  %PLL version 1
  %v   = k1*error + k2*error + M1;
  %tau = tau +v;
  %M1  = k2*error + M1; %memories update
  %M2  = tau;
  
  %PLL version 2
  accError = accError + error;
  tau = tau + k1*error + k2*accError;
  
  tnow=tnow+M; 
  tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,2,1), plot(xI, xQ,'.')        % plot constellation diagram
title('Timing sychro. Input constellation diagram'); axis square; grid
subplot(2,2,2), plot(zI, zQ,'.')        % plot constellation diagram
title('Timing sychro. Output constellation diagram'); axis square; grid
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

