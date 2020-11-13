function [zI, zQ] = timingRecoveryPLL(xI, xQ, M, dseta, Bn)
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
% xs: Decimated output samples (one per symbol)
%

%PLL constants
fs = 10000; T=1/fs;  %sampling frequency OJO Que en el main es otra
%Bn = 25;            %noise Bandwidth(Hz)
kp = 1;            %phase gain
k0 = 1;            %VCO gain
%dseta = 1;         %damping factor

%Computing the P-I loop-filter constants; k1 and k2
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2

% clock recovery algorithm
lside=20;        %one sided length of data to interpolate
tnow=2*lside+1;  % initialize variables
tau=0;
x = sqrt(xI.^2 + xQ.^2);
xs = zeros(size(x));
xsI=zeros(size(xI));
xsQ=zeros(size(xQ));
tausave=zeros(size(xI)); tausave(1)=tau; %for plotting
i=0;         %loop index
%mu=0.01;    % gradient algorithm stepsize
delta=0.1;   % time for derivative
M1=0; M2=0; accError=0; % PLL memories
while tnow<length(xQ)-lside  % run iteration
  i=i+1;
  xsI(i)    = interpsinc(xI,tnow+tau,lside);        % interp value at tnow+tau
  xsQ(i)    = interpsinc(xQ,tnow+tau,lside);        % interp value at tnow+tau
  xs(i)    = interpsinc(x,tnow+tau,lside);        % interp value at tnow+tau
  
  x_deltap = interpsinc(x,tnow+tau+delta, lside); % value to right
  x_deltam = interpsinc(x,tnow+tau-delta,lside);  % value to left
  
  dx       = x_deltap-x_deltam;                   % numerical derivative

  %tau=tau+mu*dx*xs(i);     % gradient alg. update: OP-TED
%   error = atan((dxQ*xsQ(i))/(dxI*xsI(i)));
%   error = sqrt((dxQ*xsQ(i))^2+(dxI*xsI(i))^2);
  error = dx*xs(i);                   % PLL alg. update: OP-TED

  %PLL version 1
  v     = error*k1 + error*k2 + M1;   %loop filter
  tau   = tau + v;
  M1    = error*k2 + M1;  %memories update
  M2    = tau;
  
  %PLL version 2
  %accError = accError + error;
  %tau = tau + k1*error + k2*accError;

  
  tnow=tnow+M; 
  tausave(i)=tau;      % save for plotting
end

% plot results
figure;
subplot(2,1,1), plot(xsI(1:i-2),'b.');
hold on;
plot(xsQ(1:i-2),'r.');        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

zI = xsI(1:i-2);
zQ = xsQ(1:i-2);
