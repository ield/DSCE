%% 6.05. Creation of x(f)
% Creation of a traingle signal
T = 100;  
fs = 1; 
t = -(T/fs)/2:1/fs:(T/fs)/2;
x = (sawtooth(2*pi*0.01*(t-50),1/2)+1)*0.5;  
plot(x)
% plot(t,x) 

% Replacing signal in x axis
x = [x zeros(1, 349*2+1) x];
t = -450:1:450;
% plot(t, x)

y = conv(x, x);
t = -900:1:900;
plot(t, y)

