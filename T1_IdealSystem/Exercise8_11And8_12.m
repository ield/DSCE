% pulseshape.m: applying a pulse shape to a text string
str='Transmit this text string';         % message to be transmitted
m=letters2pam(str); N=length(m);         % 4-level signal of length N
M=10; mup=zeros(1,N*M); mup(1:M:N*M)=m;  % oversample by M
ps=hamming(M);                           % blip pulse of width M
x=filter(ps,1,mup);                      % convolve pulse shape with data

% Adding noise to the signal (8.12)
% It is seen that the higher the noise the higher the errors. However, the
% errors are reduced with prx = ps
x = x + 1*randn(size(x));

% Pulse at the receiver is modified (8.11)
prx = ps;
% prx = sin(0.1*pi*(0:M-1));               % No aparent errors, but some errors when adding noise
% prx = cos(0.1*pi*(0:M-1));               % All errors

y=xcorr(x,prx);                          % correlate pulse with received signal
z=y(N*M:M:2*N*M-1)/(pow(prx)*M);         % downsample to symbol rate and normalize
mprime=quantalph(z,[-3,-1,1,3])';        % quantize to +/-1 and +/-3 alphabet
pam2letters(mprime)                      % reconstruct message
sum(abs(sign(mprime-m)))                 % calculate number of errors

