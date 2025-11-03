clc;
clear;
close all;

%% Task 2
%Variables
N = 10;
L = .2413;
x = .2159;
T0 = 15.9830;
k_array = [130,115,16.2];
ro_array = [2810,8500,8000];
cp_array = [960,380,500];
Name_array = {'ALuminum25', 'Aluminum30', 'Brass25', 'Brass30', 'Steel22'};


%Steady State finding


lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 14000;
total = zeros(1,tf);

u = zeros(1,8);
Fo = zeros(1,5);

HEXP = [55.399,78.553,104.987,150.169,287.308];
%
for h = 1:length(HEXP)

  if h<=2 %aluminum alpha
            alpha = (k_array(1)/(ro_array(1)*cp_array(1)));
        elseif h<=4 %brass Alpha
            alpha = (k_array(2)/(ro_array(2)*cp_array(2)));
        else %steel alpha
            alpha = (k_array(3)/(ro_array(3)*cp_array(3)));
        end

alpha_adjusted = [1.2048, 1.2048,1.7802,1.7802,1.1745];
alpha_adjusted = alpha_adjusted * (10^-5);

start = T0 + (HEXP(h) * x);

%figure; hold on;

%for x = linspace(.0381,L,8)

x = zeros(1,8);
x = linspace(.0381,L,8);

for q = 1:8
    START(q) = T0 + (HEXP(h) * x(q));

    for t = 1:tf
        TOTAL = 0;
        TOtAL = 0;
        for n = 1:N

            Lambda = ((2*n - 1) * pi) / (2 * L);
            Bn = ((-1)^(n) * (4*HEXP(h)*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
            Sum = Bn * sin(Lambda*x(q)) * exp(-Lambda^2 * alpha_adjusted(h) * t);
            TOtAL = Sum + TOtAL;
        end
        total(t) = TOtAL;
    end
u = START(q) + total;

time = linspace(1,tf,tf);

Difference = diff(total);
tx = 0;
for j = 1:length(Difference)
    if Difference(j) >= .001
        tss = time(j);
    end
end
figure;
plot(time(1:tss),u(1:tss), 'b'); hold on;
plot(time(tss:length(time)), u(tss:length(u)), 'r');
xline(tss);
legend('Transient State', 'Steady State', 'Steady State Time', location='south');
titleNum = Name_array{h};
title(titleNum);
ylim([0,100]);
xlim([0,tf]);
xline(tss);
text(4000, 90, ['Tss = ', num2str(tss), ' seconds'], 'HorizontalAlignment', 'right', 'FontSize', 10);

Fo(h) = (alpha_adjusted(h) * tss) / (L * L);
plot(time,u, 'r'); hold on;
legend('H_Experimental');

end






