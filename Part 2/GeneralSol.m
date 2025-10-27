clc;
clear;
close all;

%% Task 1
%t=1s
N = 10;
L = .2413;
x = .2159;
H = 91.0871;
T0 = 15.9830;
t = 1;
k = 130;
ro = 2810;
cp = 960;
alpha = k / (ro * cp);
total1 = 0;
total2= 0;
start = T0 + (H * x);
TOTAL1 = zeros(1,N);
TOTAL2 = zeros(1,N);
lambda1 = zeros(1,N);
bn1 = zeros(1,N);
sum1 = zeros(1,N);
den1 = 0;
den = 0;

for n = 1:N
lambda1(n) = ((2*n - 1) * pi) / (2 * L);
bn1(n) = ((-1)^(n) * (4*H*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
den1 = ((-1)^(1) * (4*H*L)) / (2*1-1);
den = 2/((2*1-1)*pi);
sum1(n) = bn1(n) * sin(lambda1(n)*x) * exp(-lambda1(n)^2 * alpha * t);
total1 = sum1(n) + total1;
TOTAL1(n) = total1;
end

u1 = start + TOTAL1;




% t = 1000s
t = 1000;
for n = 1:N
lambda2 = ((2*n - 1) * pi) / (2 * L);
bn2 = ((-1)^(n) * (4*H*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
sum2 = bn2 * sin(lambda2*x) * exp(-lambda2^2 * alpha * t);
total2 = sum2 + total2;
TOTAL2(n) = total2;
end

u2 = start + TOTAL2;

n = linspace(1,10,10);
t1 = [start,u1(1)];
t2 = [start,u2(1)];
m = linspace (0,1,2);

figure;

plot(n, u1, 'b'); hold on;
plot(n,u2, 'r');
plot(m,t1, 'b');
plot(m,t2, 'r');
xlabel('Number of iterations');
ylabel('Temperature approximation (C)');
legend('t = 1s', 't = 1000s');
grid on;

%% Task 2





