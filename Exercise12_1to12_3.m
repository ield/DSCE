% clockrecDD.m: clock recovery minimizing 4-PAM cluster variance
% to minimize J(tau) = (Q(x(kT/M+tau))-x(kT/M+tau))^2

% prepare transmitted signal
clear;

n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
% s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,2,1);                    % 2-PAM
s=pam(n,6,12);                   % 6-PAM


sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal

r = r + 0.3*randn(size(r));     % Add noise (12.2)

matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.2;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  
%   qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
%   qx=quantalph(xs(i),[-1,1]);  % quantize to alphabet
  qx=quantalph(xs(i),[-5, -3,-1,1,3, 5]);  % quantize to alphabet
  
  
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%% Sol 12_1
Sol12_1_a = 'For mu >= 1.4 some errors start to appear'
Sol12_1_b = '2-PAM: For mu >= 7.8 some errors start to appear'
Sol12_1_b = '6-PAM: For mu >= 0.4 some errors start to appear'

%% Sol 12_3
Sol12_3 = 'Yes, with noise, it is necessary to reduce the value of mu';


%% Adapt code to PLL
clear;

n=10000;                         % number of data points
m=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,m,toffset);  % srrc pulse shape with timing offset
s=pam(n,4,5);                    % random data sequence with var=5
% s=pam(n,2,1);                    % 2-PAM
% s=pam(n,6,12);                   % 6-PAM


sup=zeros(1,n*m);                % upsample the data by placing...
sup(1:m:n*m)=s;                  % ... m-1 zeros between each data point
hh=conv(pulshap,chan);           % ... and pulse shape
r=conv(hh,sup);                  % ... to get received signal

r = r + 0.3*randn(size(r));     % Add noise (12.2)

matchfilt=srrc(l,beta,m,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

% clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=0.1;                          % time for derivative

% PLL parameters
T=m;    %sampling frequency, see line 73
Bn = 0.005;            %noise Bandwidth(Hz). Addapted to the sampling time, so that it is proportional to the values of the costas loop cover
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

M1=0;M2=0; accError=0;

%Computing the P-I loop-filter constants; k1 and k2. How to do it is
%explained in slide 22
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

while tnow<length(x)-2*l*m          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
%   qx=quantalph(xs(i),[-1,1]);  % quantize to alphabet
%   qx=quantalph(xs(i),[-5, -3,-1,1,3, 5]);  % quantize to alphabet
  
  error = dx*(qx-xs(i));
  v = error*k1 + error*k2 + M1; %loop filter
  
  tau=tau+v;         % alg update: DD  
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
  
  M1 = error*k2 + M1;  %memories update
  M2 = tau;
end

% plot results
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

