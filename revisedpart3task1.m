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

%% experimental data

for j=1:length(Cases)

    %time for each case
    t = Cases{j}(:,1);

    %plotting exp data
    figure(j)
    hold on
    ylabel('Temperature (°C)')
    xlabel('Time (s)')
    xlim([0 1000])
    title(['\alpha against \alpha_{adj} for ', CaseNames{j}])
    %for loop for experimental data
    for i=2:9
        Exper = Cases{j};
        plot(t,Exper(:,i),'b')
    end
    hold off
end

%% Constants
N=10;
L = .2413;
x = .2159;
H = 91.0871;
T0 = [15.983,15.741,14.539,13.916,9.6274];
t = 1;
k_array = [130,115,16.2];
ro_array = [2810,8500,8000];
cp_array = [960,380,500];
total1 = 0;
total2= 0;
start = T0 + (H * x);
TOTAL1 = zeros(1,N);
TOTAL2 = zeros(1,N);
lambda1 = zeros(1,N);
bn1 = zeros(1,N);
sum1 = zeros(1,N);

%% Task 2 overlay

%Model 1B

lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 1000;
total = zeros(1,tf);
Total = zeros(1,tf);

u = zeros(1,8);
HAN = [91.0871,132.0763,101.61,146.7295,554.0675];
HEXP = [55.399,78.553,104.987,150.169,287.308];

for h = 1:length(HEXP)
    if h<=2 
        alpha = (k_array(1)/(ro_array(1)*cp_array(1)));
    elseif h<=4
        alpha = (k_array(2)/(ro_array(2)*cp_array(2)));
    else
        alpha = (k_array(3)/(ro_array(3)*cp_array(3)));
    end

    Hexp = HEXP(h);
    T_0 = T0(h);

    START = T_0 + (Hexp * x);

    x = linspace(.0381,L,8);

for q=1:8
    START(q) = T_0 + (Hexp * x(q));

    for t = 1:tf
        TOTAL = 0;
        for n = 1:N
            lambda(n) = ((2*n - 1) * pi) / (2 * L);
            bn(n) = ((-1)^(n) * (4*Hexp*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
            sum(n) = bn(n) * sin(lambda(n)*x(q)) * exp(-lambda(n)^2 * alpha * t);
            TOTAL = sum(n) + TOTAL;
        end
      
        Total(t) = TOTAL;
    end
U = START(q) + Total;

time = linspace(1,1000,1000);

%figure for H_exp
figure(h)
hold on
plot(time,U, 'r');
hold off

%saving each figure
%saveas(gcf, ['Hexp_Case_' CaseNames{h} '.png']);


end

end


%% Part 3 task 1
% 
% %Model III

lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 1000;
total = zeros(1,tf);
Total = zeros(1,tf);

u = zeros(1,8);
HAN = [91.0871,132.0763,101.61,146.7295,554.0675];
HEXP = [55.399,78.553,104.987,150.169,287.308];

for h = 1:length(HEXP)

    %different alpha values for each case/material
    if h<=2 %aluminum alpha
        alpha = (k_array(1)/(ro_array(1)*cp_array(1)));
    elseif h<=4 %brass Alpha
        alpha = (k_array(2)/(ro_array(2)*cp_array(2)));
    else %steel alpha
        alpha = (k_array(3)/(ro_array(3)*cp_array(3)));
    end

    %varying alpha
    if h<=2 %aluminum best fit percentage
        alpha_adj = (alpha*0.35);
    elseif h<=4 %brass best fit percentage
        alpha_adj = (alpha*0.4);
    else %steel best fit percentage
        alpha_adj = (alpha*0.7);
    end
    
    Hexp = HEXP(h);
    T_0 = T0(h);

    START = T_0 + (Hexp * x);

    x = linspace(.0381,L,8);

for q=1:8
    START(q) = T_0 + (Hexp * x(q));

    for t = 1:tf
        TOTAL = 0;
        for n = 1:N
            lambda(n) = ((2*n - 1) * pi) / (2 * L);
            bn(n) = ((-1)^(n) * (4*Hexp*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
            sum(n) = bn(n) * sin(lambda(n)*x(q)) * exp(-lambda(n)^2 * alpha_adj * t);
            TOTAL = sum(n) + TOTAL;
        end
      
        Total(t) = TOTAL;
    end
U = START(q) + Total;

time = linspace(1,1000,1000);

%figure for H_exp
figure(h)
hold on
plot(time,U, 'g');
hold off

%saving each figure
%saveas(gcf, ['Hexp_Case_' CaseNames{h} '.png']);


end

end