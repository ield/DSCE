% clockrecOP.m:  clock recovery maximizing output power
% find tau to optimize J(tau)=|x(kT/m+tau)|^2

% prepare transmitted signal
n=5000;                            % number of data points
m=2;                               % oversampling factor
constel=2;                         % 2-pam constellation
beta=0.5;                          % rolloff parameter for srrc
l=50;                              % 1/2 length of pulse shape (in symbols)
chan=[1];                          % T/m "channel"
toffset=-0.3;                      % initial timing offset
pulshap=srrc(l,beta,m,toffset);    % srrc pulse shape
s=pam(n,constel,1);                % random data sequence with var=1
% s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,6,1);                    % 2-PAM

sup=zeros(1,n*m);                  % upsample the data by placing...
sup(1:m:end)=s;                    % ... p zeros between each data point
hh=conv(pulshap,chan);             % ... and pulse shape
r=conv(hh,sup);                    % ... to get received signal

r = r + 0.5*randn(size(r));     % Add noise (12.12)

matchfilt=srrc(l,beta,m,0);        % filter = pulse shape
x=conv(r,matchfilt);               % convolve signal with matched filter

% run clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.1;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*m            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  tau=tau+mu*dx*xs(i);              % alg update (energy)
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%
Sol12_10a = 'For 2-PAM, the maximum mu is 1.7. Higher mu imply not convergence'
Sol12_10b = 'Higher modulations need smaller mu'
Sol12_10b1 = 'For 4-PAM, the maximum mu is 0.2. Higher mu imply not convergence'
Sol12_10b2 = 'For 6-PAM, the maximum mu is 0.075. Higher mu imply not convergence'

Sol12_12a = 'For 2-PAM, var = 0.1, the maximum mu is 1.6.'
Sol12_12b = 'For 2-PAM, var = 0.3, the maximum mu is 0.7.'
Sol12_12c = 'For 2-PAM, var = 0.5, it is too noisy.'
Sol12_12d = 'For 2-PAM, var = 0.7, it is too noisy.'

%% Adapt code to PLL

% prepare transmitted signal
n=5000;                            % number of data points
m=2;                               % oversampling factor
constel=2;                         % 2-pam constellation
beta=0.5;                          % rolloff parameter for srrc
l=50;                              % 1/2 length of pulse shape (in symbols)
chan=[1];                          % T/m "channel"
toffset=-0.3;                      % initial timing offset
pulshap=srrc(l,beta,m,toffset);    % srrc pulse shape
s=pam(n,constel,1);                % random data sequence with var=1
% s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,6,1);                    % 2-PAM

sup=zeros(1,n*m);                  % upsample the data by placing...
sup(1:m:end)=s;                    % ... p zeros between each data point
hh=conv(pulshap,chan);             % ... and pulse shape
r=conv(hh,sup);                    % ... to get received signal

r = r + 0.1*randn(size(r));     % Add noise (12.12)

matchfilt=srrc(l,beta,m,0);        % filter = pulse shape
x=conv(r,matchfilt);               % convolve signal with matched filter

% run clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.1;                            % algorithm stepsize
delta=0.1;                          % time for derivative

% PLL parameters
fs = 500; T=1/fs;    %sampling frequency, see line 73
Bn = 5;        %noise Bandwidth(Hz). Addapted to the sampling time, so that it is proportional to the values of the costas loop cover
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

M1=0;M2=0; accError=0;

%Computing the P-I loop-filter constants; k1 and k2. How to do it is
%explained in slide 22
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

while tnow<length(x)-l*m            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  
  % tau=tau+mu*dx*xs(i);              % alg update (energy)
  error = dx*xs(i);
  
  v = error*k1 + error*k2 + M1; %loop filter
  tau=tau+v;                    % alg update: DD  
  
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
  
  M1 = error*k2 + M1;  %memories update
  M2 = tau;
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%% Normal code and unagreement between m in tx and rx

% prepare transmitted signal
n=5000;                            % number of data points
m=2;                               % oversampling factor
constel=2;                         % 2-pam constellation
beta=0.5;                          % rolloff parameter for srrc
l=50;                              % 1/2 length of pulse shape (in symbols)
chan=[1];                          % T/m "channel"
toffset=-0.3;                      % initial timing offset
pulshap=srrc(l,beta,m,toffset);    % srrc pulse shape
s=pam(n,constel,1);                % random data sequence with var=1
% s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,6,1);                    % 2-PAM

sup=zeros(1,n*m);                  % upsample the data by placing...
sup(1:m:end)=s;                    % ... p zeros between each data point
hh=conv(pulshap,chan);             % ... and pulse shape
r=conv(hh,sup);                    % ... to get received signal

r = r + 0.3*randn(size(r));     % Add noise (12.12)

matchfilt=srrc(l,beta,m,0);        % filter = pulse shape
x=conv(r,matchfilt);               % convolve signal with matched filter


% clockrecperiod.m: resample to change period It is easy to understant if
% fac=1
% This effect is analogous to have a frequency difference, so it can be
% corrected if it is not too low
fac=1.001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal

% run clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.8;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*m            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  tau=tau+mu*dx*xs(i);              % alg update (energy)
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%
Sol12_25 = 'For 2-PAM, and var = 0.1, the maximum fac is 1.001.'

%% Adapt code to PLL and unagreement between mm in tx and rx
clear;

% prepare transmitted signal
n=5000;                            % number of data points
m=2;                               % oversampling factor
constel=2;                         % 2-pam constellation
beta=0.5;                          % rolloff parameter for srrc
l=50;                              % 1/2 length of pulse shape (in symbols)
chan=[1];                          % T/m "channel"
toffset=-0.3;                      % initial timing offset
pulshap=srrc(l,beta,m,toffset);    % srrc pulse shape
s=pam(n,constel,1);                % random data sequence with var=1
% s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,6,1);                    % 2-PAM

sup=zeros(1,n*m);                  % upsample the data by placing...
sup(1:m:end)=s;                    % ... p zeros between each data point
hh=conv(pulshap,chan);             % ... and pulse shape
r=conv(hh,sup);                    % ... to get received signal

% r = r + 0.1*randn(size(r));     % Add noise (12.12)

matchfilt=srrc(l,beta,m,0);        % filter = pulse shape
x=conv(r,matchfilt);               % convolve signal with matched filter

% clockrecperiod.m: resample to change period
fac=1.05; z=zeros(size(x));   % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal

% run clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
delta=0.1;                          % time for derivative

% PLL parameters
fs = 500; T=1/fs;    %sampling frequency, see line 73
Bn = 25;        %noise Bandwidth(Hz). Addapted to the sampling time, so that it is proportional to the values of the costas loop cover
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

M1=0;M2=0; accError=0;

%Computing the P-I loop-filter constants; k1 and k2. How to do it is
%explained in slide 22
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

while tnow<length(x)-l*m            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  
  dx=x_deltap-x_deltam;             % numerical derivative
  
  % tau=tau+mu*dx*xs(i);              % alg update (energy)
  error = dx*xs(i);
  
  v = error*k1 + error*k2 + M1; %loop filter
  tau=M2;                    % alg update: DD  
  
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
  
  M1 = error*k2 + M1;  %memories update
  M2 = M2 + v;
  
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')