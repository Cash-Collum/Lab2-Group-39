clc;
clear;
close all;

n = 0;
L = 1;
x = 0.123825;
H = 91.0871;
T0 = 15.9830;
t = linspace(0,1000,10000);
alpha = 4.8191 * 10^-5;

for n = 1:10
lambda = ((2*n - 1) * pi) / (2 * L);
bn = ((-2 * H) / (lambda*L)) * ((sin(L * lambda) / lambda) - L * cos(L * lambda));
sum = bn .* sin(lambda .* x) .* exp(-1 .* lambda.^2 .* alpha .* t);
end;

u = T0 + (H * x) + sum;

plot(t,u, 'r', lineWidth=1.5);


