function [s] = carrierRecovery(r, f0, Ts, delta_f0)
% Amplitude gain control for task 3 based on schematic on page 2

n = length(r);
phi = zeros(size(r));

s = zeros(size(r)); s(1)=exp(1i*(2*pi*f0*Ts*1+phi(1)));   
v = zeros(size(r));

fs = 1/Ts;
factor = 3/4*pi*3/2*fs;
k1 = delta_f0/factor;
k1 = 0.07;
k2 = k1^2;

mem_ek2 = 0;
mem_w = 0;

for k=1:n-1
  e = angle(r(k)*conj(s(k)));
  
  v(k) = e*k1 + e*k2+mem_ek2;
  mem_ek2 = e*k2;
  
  w = v(k) + mem_w;
  mem_w = w;
  
  s(k + 1) = exp(1i * w);
end

% draw
figure;
subplot(3,1,1)
plot(v/(2*pi))              % plot AGC values
title('Frequency')
subplot(3,1,2)
plot(real(s),'r')          % plot inputs and outputs
xlim([0 1000])
subplot(3,1,3)
plot(real(r),'b')
xlim([0 1000])

end

