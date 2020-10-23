% agcvsfading.m: compensating for fading with an AGC
n=50000;                           % # steps in simulation

% Generating the sin function
fs = 8e3;
t = 1/fs:1/fs:n/fs;
r=5*cos(2*pi*1000*t);              % generate the signal

% Fading profile
% The specifications state that the signal must be 1, then reduced 16 db
% and then back to 1.
n1 = 15000;
n2 = 35000;
env1 = ones(1, n1);                 % the fading profile
att = 16;                           % Attenuation of the faded part in dB
env2 = 10^(-att/20)*ones(1, n2-n1);
env3 = ones(1, n-n2);
env = [env1 env2 env3];

r=r.*env;                          % apply to raw input r[k]
ds=0.5;                            % desired power of output
a=zeros(1,n); a(1)=1;              % initialize AGC parameter
s=zeros(1,n);                      % initialize outputs
mu=0.001;                           % algorithm stepsize
for k=1:n-1
  s(k)=a(k)*r(k);                  % normalize by a to get s
  a(k+1)=a(k)-mu*(s(k)^2-ds);      % adaptive update of a(k)
end


% draw agcgrad.eps
subplot(3,1,2), plot(a,'g')        % plot AGC values
axis([0,length(r),0,1.5])
title('Adaptive gain parameter')
subplot(3,1,1), plot(r,'r')        % plot inputs and outputs
axis([0,length(r),-7,7])
title('Input r(k)')
subplot(3,1,3),plot(s,'b')
axis([0,length(r),-7,7])
title('Output s(k)')
xlabel('iterations')


