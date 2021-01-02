function [s] = agc(r)
% Amplitude gain control for task 3 based on schematic on page 2

n=length(r);                           % number of steps in simulation

ds=5;                           % desired power of output
mu=0.001;                          % algorithm stepsize
lenavg=10;                         % length over which to average
a=zeros(n,1); a(1)=1/max(abs(r));              % initialize AGC parameter
s=zeros(n,1);                      % initialize outputs
avec=zeros(1,lenavg);              % vector to store terms for averaging
for k=1:n-1
  s(k)=a(k)*r(k);                  % normalize by a(k)
%   avec=[sign(a(k))*(abs(s(k))^2-ds),avec(1:lenavg-1)];  % incorporate new update into avec
  a(k+1)=a(k)-mu*(1-abs(s(k))^2);       % average adaptive update of a(k)
end

% draw agcgrad.eps
subplot(3,1,1)
plot(a)              % plot AGC values
title('Adaptive gain parameter')
subplot(3,1,2)
plot(r,'r')          % plot inputs and outputs
axis([0,10^4,-5,5])
title('Input r(k)')
subplot(3,1,3)
plot(s,'b')
axis([0,10^4,-5,5])
title('Output s(k)')
xlabel('iterations')

end

