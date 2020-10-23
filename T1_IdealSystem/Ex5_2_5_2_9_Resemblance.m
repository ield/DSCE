%% Author: ield
% Description: Plot made to see the resembplance point of s1 and s2 for
% exercise 5.2.9 in task 5.2
f1 = 200;
f2 = 400;
fc = 2000;

fs = 8000;
T = 0.02;
t = 1/fs:1/fs:T;

s1 = cos(2*pi*f1*t);
s1Phase = 0.25*(cos(2*pi*(f1+20)*t)+cos(2*pi*(f1-20)*t));

s2 = -sin(2*pi*f2*t);
s2Phase = -0.25*(sin(2*pi*(f1+20)*t)+sin(2*pi*(f1-20)*t));


plot(t, s1-s1Phase);
hold on;
plot(t, s2-s2Phase);
