clear;
% LMSequalizer.m find a LMS equalizer f for the channel b
b=[1 1 -0.8 -0.3 1 1];             % define channel
m=1000; s=pam(m,2,1);        % binary source of length m
r=filter(b,1,s);             % output of channel
n=25; f=zeros(n,1);           % initialize equalizer at 0. The length of the equalizer cannot be too large because it conflicts with mu so there can be divergence.
mu=.01;     % This constant is important because if the length of the equalizer is too large and it converges too fast there can be inestabilities.
delta=floor((length(b)+n)/2);% n-delta >= 0;  (length(b)+delta)/2 = delta
for i=n+1:m                  % iterate
  rr=r(i:-1:i-n+1)';         % vector of received signal
  e=s(i-delta)-rr'*f;        % calculate error
  f=f+mu*e*rr;               % update equalizer coefficients
end
y=filter(f,1,r);                   % equalizer is a filter
dec=sign(y);                       % quantize and find errors
err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))

plotEqualizers(b, f, y, dec);
