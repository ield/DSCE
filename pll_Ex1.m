%PLL Example 1
%complex exponential PLL
clear;
%PLL parameters
f0 = 50;           %exponential signal frequency (Hz):  e^(j*2*pi*f0*t)
fs = 500; T=1/fs;  %sampling frequency
Bn = 25;            %noise Bandwidth(Hz)
kp = 1;            %phase gain
k0 = 1;            %VCO gain
dseta = 1;         %damping factor

%signal perturbations
theta_in  = pi; %phase offset (radians)
deltaf0   = 0;  %frequency offset (Hz)
Pn        = 0;  %Noise power (in wats)

%Computing the P-I loop-filter constants; k1 and k2. How to do it is
%explained in slide 22
wn=Bn*T/( dseta+1/(4*dseta) );       %natural frequency 
k1 = 4*dseta*wn/(1+2*dseta*wn+wn^2); %k1=0.1479;
k2 = 4*wn^2/(1+2*dseta*wn+wn^2);     %k2=0.0059;

%Acquistion theoretical results
f_pullin = 2*pi*sqrt(2*dseta)*Bn;   % Diapo 15. Pull-in frequency
tacq = 2/Bn;                        % Diapo 15. Acquisition time: time it takes to switch to tracking
disp(['f_pullin (teórico): ' num2str(f_pullin) ' Hz']);
disp(['tacq (teórico): ' num2str(tacq) ' seg']);


%discrete-time simulation 
N         = 1000; %iterations
theta(1) = 0;
w0        = 2*pi*f0*T;      %frequency (rad/sample)
deltaw0   = 2*pi*deltaf0*T; %frequency offset (rad/sample)

s_inVector   = zeros(1, N);
error_vector = zeros(1, N);
Real_exp_in  = zeros(1, N);
Real_exp_est = zeros(1, N);
M1=0;M2=0; accError=0;
for n=1:N-1,    
  %error
  exp_in = exp( j*((w0+deltaw0)*n + theta_in) ); %input signal with phase and freq offset 
  s_in  =  exp_in + (0.5*sqrt(Pn)*randn + j*0.5*sqrt(Pn)*randn);    %signal plus noise
  s_est =  exp(j*w0*n+theta(n));
  s_mod = s_in*conj(s_est); 
  error= angle(s_mod);

  %data for plots
  s_inVector(n+1)   = real(s_in); 
  Real_exp_in(n+1)  = real(exp_in);
  Real_s_est(n+1)   = real(s_est);
  error_vector(n+1) = error;
  
  %PLL version 1
  v = error*k1 + error*k2 + M1; %loop filter
  theta(n+1) = theta(n) + v + w0;   %w0 can be cancelled
  
  M1 = error*k2 + M1;  %memories update
  M2 = theta(n+1);
  
  %PLL version 2
  %accError  = accError + error;
  %theta(n+1)= theta(n) + k1*error + k2*(error+ accError)+w0;
  
end;

subplot(311)
plot(0:N-1, s_inVector); axis tight; grid
title('Re\{s_{in}\}')

subplot(312)
plot(0:N-1, Real_exp_in, 'b'); hold on
plot(0:N-1, Real_s_est, 'r'); axis tight; grid; hold off
title('Re\{exp_{in}\} (blue) and Re\{s_{est}\} (red)')

subplot(313)
plot(0:N-1, error_vector); grid; axis tight
xlabel('samples'); title('error')



