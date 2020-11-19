function s=base_QPSK(Lt,T,Nz)

if nargin<3; Nz=0; end   %Duracion en simbolos de la cola de ceros
Ar=0;
if nargin<2; T=11; end  % Duración del símbolo
if nargin<1; Lt=5000; end % Duracion de la secuencia
N=floor(Lt/T)-Nz; % Numero de símbolos
Ln = T*N; % Duracion de la secuencia de informacion

rng(1,'twister');
a = randi(4,[1,N]);
vA=exp(1i*[1,3,5,7]*pi/4);

p=hanning(T)';

b=zeros(1,T*N);
b(1+T*(0:N-1))=vA(a);
c=filter(p,1,b);

r = c + Ar*(randn(size(c))+1i*randn(size(c)))/sqrt(2);
D = (T+rem(T,2))/2;
x = filter(p(end:-1:1),1,r)/sqrt(2*T);
Dr = 2*D-1;
s=[c,zeros(1,Lt-Ln)];


%{
figure(1)
subplot(211); plot(1:4*T,[real(c(1:4*T))',imag(c(1:4*T))']);
subplot(223); plot(r(D:T:end),'.'); axis ('square'); grid on

figure(2)
subplot(211); plot(1:4*T,[real(x(1:4*T))',imag(x(1:4*T))']);
subplot(223); plot(x(Dr:T:end),'.'); axis ('square'); grid on
%}

