%% 6.28. Creation of J(x)
x = -5:1:5;

j_x = (1 - abs(x-2)).^2;

plot(x, j_x);