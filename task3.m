clc;
clear;
close all;

Aluminum25V = readmatrix("Aluminum_25V_240mA");
Aluminum30V = readmatrix("Aluminum_30V_290mA");
Brass25V = readmatrix("Brass_25V_237mA");
Brass30V = readmatrix("Brass_30V_285mA");
Steel22V = readmatrix("Steel_22V_203mA");

Cases = {Aluminum25V,Aluminum30V,Brass25V,Brass30V,Steel22V};
CaseNames = { ...
    'Aluminum 25V 240mA', ...
    'Aluminum 30V 290mA', ...
    'Brass 25V 237mA', ...
    'Brass 30V 285mA', ...
    'Steel 22V 203mA'};

%% experimental data overlayed with H_an

for j=1:length(Cases)

    %time for each case
    t = Cases{j}(:,1);

    %plotting exp data
    figure(j)
    hold on
    ylabel('Temperature (°C)')
    xlabel('Time (s)')
    xlim([0 1000])
    title(['H_{exp} over Experimental data for ', CaseNames{j}])
    %for loop for experimental data
    for i=2:8
        Exper = Cases{j};
        plot(t,Exper(:,i),'b')
    end
    hold off
end

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

for n = 1:N
lambda1(n) = ((2*n - 1) * pi) / (2 * L);
bn1(n) = ((-1)^(n) * (4*H*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
sum1(n) = bn1(n) * sin(lambda1(n)*x) * exp(-lambda1(n)^2 * alpha * t);
total1 = sum1(n) + total1;
TOTAL1(n) = total1;
end

u1 = start + TOTAL1;

%% Task 2 overlay

%Model 1A

lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 1000;
total = zeros(1,tf);
Total = zeros(1,tf);

u = zeros(1,8);
HAN = [91.0871,132.0763,101.61,146.7295,554.0675];
HEXP = [55.399,78.553,104.987,150.169,287.308];

for h = 1:length(HAN)

        %different alpha values for each case/material
        if h<=2 %aluminum alpha
            alpha = (k_array(1)/(ro_array(1)*cp_array(1)));
        elseif h<=4 %brass Alpha
            alpha = (k_array(2)/(ro_array(2)*cp_array(2)));
        else %steel alpha
            alpha = (k_array(3)/(ro_array(3)*cp_array(3)));
        end

start = T0 + (HAN(h) * x);
START = T0 + (HEXP(h) * x);

for x = linspace(.0381,L,8)
    for t = 1:tf
        TOTAL = 0;
        TOtAL = 0;
        for n = 1:N
            lambda(n) = ((2*n - 1) * pi) / (2 * L);
            bn(n) = ((-1)^(n) * (4*HAN(h)*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
            sum(n) = bn(n) * sin(lambda(n)*x) * exp(-lambda(n)^2 * alpha * t);
            TOTAL = sum(n) + TOTAL;

            Lambda = ((2*n - 1) * pi) / (2 * L);
            Bn = ((-1)^(n) * (4*HEXP(h)*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
            Sum = Bn * sin(Lambda*x) * exp(-Lambda^2 * alpha * t);
            TOtAL = Sum + TOtAL;
        end
        total(t) = TOtAL;
        Total(t) = TOTAL;
    end
%u = start + total;
U = START + Total;

time = linspace(1,1000,1000);

%figure for H_exp
figure(h)
hold on
plot(time,U, 'r');
legend('Experimental Data','','','','','','','','H_{exp}','Location','best')
hold off

%saving each figure
saveas(gcf, ['Hexp_Case_' CaseNames{h} '.png']);


end

end



