clc;
clear;
close all;


%% Aluminum 25V experimental data
Aluminum25V = readmatrix("Aluminum_25V_240mA");
t = Aluminum25V(:,1);

figure(1)
hold on
for i=2:8
    plot(t,Aluminum25V(:,i),'b')
end;
ylabel('Temperature (°C)')
xlabel('Time (s)')
title('Temperature vs. Time for Aluminum 25V 240mA')
hold off

%% Aluminum 30V experimental data
Aluminum30V = readmatrix("Aluminum_30V_290mA");
t = Aluminum30V(:,1);

figure(2)
hold on
for i=2:8
    plot(t,Aluminum30V(:,i),'b')
end;
ylabel('Temperature (°C)')
xlabel('Time (s)')
title('Temperature vs. Time for Aluminum 30V 290mA')
hold off

%% Brass 25V experimental data
Brass25V = readmatrix("Brass_25V_237mA");
t = Brass25V(:,1);

figure(3)
hold on
for i=2:8
    plot(t,Brass25V(:,i),'b')
end;
ylabel('Temperature (°C)')
xlabel('Time (s)')
title('Temperature vs. Time for Brass 25V 237mA')
hold off

%% Brass 30V experimental data
Brass30V = readmatrix("Brass_30V_285mA");
t = Brass30V(:,1);

figure(4)
hold on
for i=2:8
    plot(t,Brass30V(:,i),'b')
end;
ylabel('Temperature (°C)')
xlabel('Time (s)')
title('Temperature vs. Time for Brass 30V 285mA')
hold off

%% Steel 22V experimental data
Steel22V = readmatrix("Steel_22V_203mA");
t = Steel22V(:,1);

figure(5)
hold on
for i=2:8
    plot(t,Steel22V(:,i),'b')
end;
ylabel('Temperature (°C)')
xlabel('Time (s)')
title('Temperature vs. Time for Steel 22V 203mA')
hold off
