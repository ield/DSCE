% specsquare.m plot the spectrum of a square wave
f=1000;                       % "frequency" of square wave
phi = pi/4;
time=10;                    % length of time
Ts=1/10000;                  % time interval between samples
t=Ts:Ts:time;               % create a time vector
% x=sign(cos(2*pi*f*t));      % square wave = sign of cos wave
% x=5*exp(-t);                % 3.3.a
x=sin(2*pi*f*t + phi);        % 3.3.e
plotspec(x,Ts)              % call plotspec to draw spectrum
