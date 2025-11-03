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
    xlim([0 20000])
    title(['\alpha_{adj} for ', CaseNames{j}])
    %for loop for experimental data
    for i=2:9
        Exper = Cases{j};
        plot(t,Exper(:,i),'b')
    end
    hold off
end

%% Constants
N=50;
L = .2413;
H = 91.0871;
T0 = [15.983,15.741,14.539,13.916,9.227];
t = 1;
x = 0.03495 + (0:7)*0.0254; 
k_array = [130,115,16.2];
ro_array = [2810,8500,8000];
cp_array = [960,380,500];
total1 = 0;
total2= 0;
TOTAL1 = zeros(1,N);
TOTAL2 = zeros(1,N);
lambda1 = zeros(1,N);
bn1 = zeros(1,N);
sum1 = zeros(1,N);




%% Part 3 task 1
% 
% %Model III

lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 30000;
total = zeros(1,tf);
Total = zeros(1,tf);
START = zeros(1,8);
time = linspace(1,30000,30000);

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
    pct = [1.5,1.1,3.15]; % These values are chosen based off of the time required to reach steady state. Varying alpha does not effect Temperature magnitude, only rate of H.T. (Based off of 8th thermocouple)


    if h==1 %aluminum best fit percentage
        alpha_adj = (alpha*pct(1));
    elseif h==2 
        alpha_adj = (alpha*pct(1)) %Adjusted diffusivity is within normal ranges
    elseif h==3 %brass best fit percentage
        alpha_adj = (alpha*pct(2))
    elseif h==4 
        alpha_adj = (alpha*pct(2)); %Adjusted diffusivity is within normal ranges
    else %steel best fit percentage
        alpha_adj = (alpha*pct(3))%Adjusted diffusivity is within normal ranges but a little bit on the high side.
    end
  
    Hexp = HEXP(h);
    T_0 = T0(h);



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



%figure for H_exp
figure(h)
hold on
plot(time,U, 'g');
legend('Experimental Data','','','','','','','','\alpha_{adj}','Location','best')
if h<=4
    xlim([0 4000])
else
    xlim([0 12000])
end
hold off

end

%saving each figure
%saveas(gcf, ['Alpha_comparison_' CaseNames{h} '.png']);

end
