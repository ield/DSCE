clear;
% LSequalizer.m find a LS equalizer f for the channel b
b=[1 1 -0.8 -0.3 1 1];             % define channel
m=1000; s=sign(randn(1,m));        % binary source of length m
r=filter(b,1,s);                   % output of channel
n=15;                               % length of equalizer - 1
delta=5;                           % n-delta >= 0
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
fprintf('b. delta = 5 the error is drastically reduced, in 5 and 6 the delta is at the center of the F vector\n');
fprintf('c. Around 100\n');

figure;
subplot(1, 3, 1);
stem(0:length(b)-1, b);

subplot(1, 3, 2);
stem(0:length(f)-1, f)

subplot(1, 3, 3);
prod = conv(b, f);
stem(0:length(prod)-1, prod);
fprintf('e. The multipliation is the convolution of the time responses of the channel. The result must be the highest possible delta centered in the delay\n');

figure;
subplot(1, 3, 1);
plot(abs(fft(b)));

subplot(1, 3, 2);
plot(abs(fft(f)));

subplot(1, 3, 3);
plot(abs(fft(prod)));
fprintf('e. The multipliation in frequency domain. The result close to a constant (the transformation of the delta)\n');

