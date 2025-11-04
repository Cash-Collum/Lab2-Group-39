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

% for j=1:length(Cases)
% 
%     %time for each case
%     t = Cases{j}(:,1);
% 
%     %plotting exp data
%     figure(j)
%     hold on
%     ylabel('Temperature (°C)')
%     xlabel('Time (s)')
%     xlim([0 20000]) 
%     title(['H_{exp} over Experimental Data for ', CaseNames{j}])
%     %for loop for experimental data
%     for i=2:9
%         Exper = Cases{j};
%         plot(t,Exper(:,i),'b')
%     end
%     hold off
% end

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

%% Task 3 overlay

%Model 1b

lambda = zeros(1,N);
bn = zeros(1,N);
sum = zeros(1,N);
tf = 10000;
Total = zeros(1,tf);

HAN = [91.0871,132.0763,101.61,146.7295,554.0675]; %H_an for task 2
HEXP = [55.399,78.553,104.987,150.169,287.308];

for h = 1:length(HEXP)
    if h<=2 
        alpha = (k_array(1)/(ro_array(1)*cp_array(1)));
    elseif h<=4
        alpha = (k_array(2)/(ro_array(2)*cp_array(2)));
    else
        alpha = (k_array(3)/(ro_array(3)*cp_array(3)));
    end

    H = HEXP(h);
    T_0 = T0(h);


START = zeros(1,8);
time = linspace(1,30000,30000);

% for q=1:8
%     START(q) = T_0 + (H * x(q));
% 
%     for t = 1:length(time)
%         TOTAL = 0;
%         for n = 1:N
%             lambda(n) = ((2*n - 1) * pi) / (2 * L);
%             bn(n) = ((-1)^(n) * (4*H*L)) / (2*n-1) * (2 / ((2*n-1) * pi * pi));
%             sum(n) = bn(n) * sin(lambda(n)*x(q)) * exp(-lambda(n)^2 * alpha * time(t));
%             TOTAL = sum(n) + TOTAL;
%         end
% 
%         Total(t) = TOTAL;
%     end
% U = START(q) + Total;



%figure for H_exp
% err = 0.02 *ones(size(U));
% 
% figure(h)
% hold on
% plot(time,U, 'r'); %Turn on if you want Original Overlay
% errorbar(time(1:20:end), U(1:20:end), err(1:20:end), 'r')
% legend('Experimental Data','','','','','','','','H_{exp}','Location','best')
% hold off
% if h<=4
%     xlim([0 4000]);
% else
%     xlim([0 12000]);
% end
%     %saving each figure
%     saveas(gcf, ['Hexp_Case_' CaseNames{h} '.png']);
% end

% Only calculate for 8th thermocouple
q = 8;
START = T_0 + (H * x(q));
Total = zeros(1,length(time));

for t_idx = 1:length(time)
    TOTAL = 0;
    for n = 1:N
        lambda_n = ((2*n - 1) * pi) / (2 * L);
        bn_n = ((-1)^n * (4*H*L)) / (2*n-1) * (2 / ((2*n-1) * pi^2));
        TOTAL = TOTAL + bn_n * sin(lambda_n * x(q)) * exp(-lambda_n^2 * alpha * time(t_idx));
    end
    Total(t_idx) = TOTAL;
end

U = START + Total;

% Figure for overlay
figure(h)
hold on

% Plot experimental data for 8th thermocouple
Exper = Cases{h};
plot(Exper(:,1), Exper(:,9), 'b')  % blue line

% Plot experimental error bars (every 20th point)
err = 2 * ones(size(Exper(:,9)));
errorbar(Exper(1:20:end,1), Exper(1:20:end,9), err(1:20:end), 'g.', 'MarkerSize',8, 'CapSize',3)

% Plot model prediction for 8th thermocouple
plot(time, U, 'r', 'LineWidth',1.5)

legend('Experimental Data','Experimental Error','Model Prediction','Location','best')
xlabel('Time (s)')
ylabel('Temperature (°C)')
title(['H_{exp} vs Model for ', CaseNames{h}, ' (Thermocouple 8)'])

if h<=4
    xlim([0 4000]);
else
    xlim([0 12000]);
end

hold off

% Save figure
saveas(gcf, ['Hexp_Case_' CaseNames{h}, '_TC8.png']);


end