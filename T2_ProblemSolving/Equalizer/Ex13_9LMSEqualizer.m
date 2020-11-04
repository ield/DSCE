clear;
% LMSequalizer.m find a LMS equalizer f for the channel b
tests = 30;
error = zeros(20);      % Final vector of errors
b=[1 1 -0.8 -0.3 1 1];             % define channel
m=1000; s=pam(m,2,1);        % binary source of length m
r=filter(b,1,s);             % output of channel
mu=.01;                         % Stepsize
minError = 1000;
% Tests after f is optimized
s1=sign(randn(1,m));
r1=filter(b,1,s1);

for n = 1:tests
    f=zeros(n,1);           % initialize equalizer at 0
    e_fin = zeros(1, (m-(n+1)));
    for delta = 1:n               %delya <= n+2 for lin 20
        for ii=n+1:m                  % iterate
          rr=r(ii:-1:ii-n+1)';         % vector of received signal
          e=s(ii-delta)-rr'*f;        % calculate error
          f=f+mu*e*rr;               % update equalizer coefficients
        end
        
        % Try equalizer with new signal

        for ii=n+1:m                  % iterate
          rr=r1(ii:-1:ii-n+1)';         % vector of received signal
          e_fin(ii)=s1(ii-delta)-rr'*f;        % calculate error
        end

        error(n, delta) = sum(abs(e_fin));
        minError = min(error(n, delta), minError);
        
    end
end
[n, delta] = find(error==minError);
fprintf('There is no delta giving 0 error, the minimum error is %f obtained for n = % i and delta = %i.\n', minError, n, delta);
fprintf('The answer here is more random');

%Plot the errors
n = 1:tests;
delta = 1:tests;
[n, delta] = meshgrid(n, delta);
surf(n, delta, error);
xlabel('n')
ylabel('delta')



