% Representation of a sin wave. It is seen when there is alliasing. It is
% when the frequancy is smaller than the nyquist frequency: the sin
% corresponds to the one at another frequency because of sampling errors.
f=50;                       % "frequency" of square wave
time=100;                    % length of time
Ts=1/60;                  % time interval between samples
t=Ts:Ts:time;               % create a time vector
x=cos(2*pi*f*t);        
plotspec(x,Ts)              % call plotspec to draw spectrum