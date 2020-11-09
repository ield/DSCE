clear;
% LSequalizer.m find a LS equalizer f for the channel b
b=[1 1 -0.8 -0.3 1 1];             % define channel
m=1000; s=sign(randn(1,m));        % binary source of length m
r=filter(b,1,s);                   % output of channel
n=20;                               % length of equalizer - 1
delta=floor((length(b)+n)/2)       % n-delta >= 0  (length(b)+delta)/2 = delta
p=length(r)-delta;
R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
S=s(n+1-delta:p-delta)';           % and vector S
f=inv(R'*R)*R'*S                   % calculate equalizer f, as in diapo 17
Jmin=S'*S-S'*R*inv(R'*R)*R'*S      % Jmin for this f and delta
y=filter(f,1,r);                   % equalizer is a filter
dec=sign(y);                       % quantize and find errors
err=0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))

%Why is the delay related to the length of the equalizer? (see line 7)

fprintf('a. For n >= 10 the error is very small\n');
fprintf('b. floor((length(b)+n)/2) the error is drastically reduced, the delta is at the center of the convolution and the greates values of the equalizer are centered.\n');
fprintf('c. Around 100\n');

plotEqualizers(b, f, y, dec);