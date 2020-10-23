% specnoise.m plot the spectrum of a noise signal
time=1;                     % length of time
Ts=1/10000;                 % time interval between samples
% x=-1 + 2*rand(1,time/Ts);    % 3.4.a Ts points of noise for time seconds. The signal is adapted to be in [-1 1]
% x=sign(-1 + 2*rand(1,time/Ts));    % 3.4.b
% x=3^0.5*sign(-1 + 2*rand(1,time/Ts));    % 3.4.c It is multiplied by sqrt(3) so that var(x) = 3
x=3^0.5*randn(1, time/Ts);
plotspec(x,Ts)              % call plotspec to draw spectrum