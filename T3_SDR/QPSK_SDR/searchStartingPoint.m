function [startingPoint] = searchStartingPoint(y)   %y is already realor imaginary
% The starting point is defined as the point in which there are 50
% consecutive 0s. It is looked for the point where this chain starts. 
Iy = diff(real(y));
% Iy = abs(Iy)/max(abs(Iy));
Iy = abs(Iy)/rms(Iy);
Iy = Iy>0.01;
plot(Iy)


% With only 25 zeros is enough
numZeros = 25;

for ii = 1:length(Iy)
    if(Iy(ii:ii+numZeros) == zeros(1,25))
        startingPoint = ii;
        return;
    end
end
end

