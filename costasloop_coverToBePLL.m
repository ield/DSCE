%% Signal generation
% pulrecsig.m: create pulse shaped received signal

N=10000; M=20; Ts=.0001;   % # symbols, oversampling factor
time=Ts*N*M; t=Ts:Ts:time; % sampling interval & time vector
m=pam(N,4,5);              % 4-level signal of length N
mup=zeros(1,N*M); 
mup(1:M:N*M)=m;            % oversample by integer length M
ps=hamming(M);             % blip pulse of width M
s=filter(ps,1,mup);        % convolve pulse shape with data
fc=1000; phoff=-1.0;       % carrier freq. and phase
c=cos(2*pi*fc*t+phoff);    % construct carrier
rsc=s.*c;                  % modulated signal (small carrier)
rlc=(s+1).*c;              % modulated signal (large carrier)

%% PLL parameters
f0 = 50;           %exponential signal frequency (Hz):  e^(j*2*pi*f0*t)
fs = 500; T=Ts;  %sampling frequency
Bn = 25;            %noise Bandwidth(Hz)
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

M1=0;M2=0; accError=0;

%Computing the P-I loop-filter constants; k1 and k2. How to do it is
%explained in slide 22
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

%% Costas Loop modified
% costasloop.m simulate costas loop with input from pulrecsig.m
r=rsc;                                % rsc from pulrecsig.m
fl=500; ff=[0 .01 .02 1]; fa=[1 1 0 0];
h=firpm(fl,ff,fa);                    % LPF design
mu=.003;                              % algorithm stepsize
f0=1000;                              % freq. at receiver
theta=zeros(1,length(t)); theta(1)=0; % estimate vector
zs=zeros(1,fl+1); zc=zeros(1,fl+1);   % buffers for LPFs
for k=1:length(t)-1                   % z contains past inputs
  zs=[zs(2:fl+1), 2*r(k)*sin(2*pi*f0*t(k)+theta(k))];
  zc=[zc(2:fl+1), 2*r(k)*cos(2*pi*f0*t(k)+theta(k))];
  lpfs=fliplr(h)*zs'; lpfc=fliplr(h)*zc'; % output of filters
  
%   theta(k+1)=theta(k)-mu*lpfs*lpfc;   % algorithm update
  error = -1*lpfs*lpfc;                % -1 because it is a function to maximize 
  
  v = error*k1 + error*k2 + M1; %loop filter
  theta(k+1) = theta(k) + v;   %w0 can be cancelled
  
  M1 = error*k2 + M1;  %memories update
  M2 = theta(k+1);

end

plot(t,theta),
title('Phase Tracking via the Costas Loop')
xlabel('time'); ylabel('phase offset')
