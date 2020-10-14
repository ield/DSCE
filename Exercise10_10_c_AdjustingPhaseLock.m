% Code to solve the optimization problem in 10_10_c inspired in pllsd.m

Ts=1/10000; time=2; t=0:Ts:time-Ts;    % time interval & vector
fc=100; alphaoff=-0.8;                    % carrier freq. and phase

mu=-0.39;                               % algorithm stepsize
beta=zeros(1,length(t)); beta(1)=0;  % initialize estimates
fl=25; h=ones(1,fl)/fl;                % averaging coefficients
z=zeros(1,fl); f0=fc;                  % buffer for avg
for k=1:length(t)-1                     % run algorithm
  jd = 0.25*sin(2*alphaoff-2*beta(k));
  z=[z(2:fl), jd];                 % z contains past inputs
  beta(k+1)=beta(k)+mu*fliplr(h)*z'; % update = z convolve h. Mu is positive because it is wanted to maximize J.
end

plot(t,beta)                             % plot estimated phase
title('Phase Tracking via SD cost')
xlabel('time'); ylabel('phase offset')


%% Option 2: not considering alpha
% pllsd.m: phase tracking minimizing SD

Ts=1/10000; time=1; t=0:Ts:time-Ts;    % time interval & vector
fc=100; alpha=-0.8;                    % carrier freq. and phase
rp=cos(2*pi*fc*t+alpha);             % simplified rec signal
mu=.2437;                               % algorithm stepsize
beta=zeros(1,length(t)); beta(1)=0;  % initialize estimates
fl=25; h=ones(1,fl)/fl;                % averaging coefficients
z=zeros(1,fl); f0=fc;                  % buffer for avg
for k=1:length(t)-1                    % run algorithm
  filtin=(rp(k)-cos(2*pi*f0*t(k)+beta(k)))*sin(2*pi*f0*t(k)+beta(k));
  z=[z(2:fl), filtin];                 % z contains past inputs
  beta(k+1)=beta(k)-mu*fliplr(h)*z'; % update = z convolve h
end

plot(t,beta)                             % plot estimated phase
title('Phase Tracking via SD cost')
xlabel('time'); ylabel('phase offset')

