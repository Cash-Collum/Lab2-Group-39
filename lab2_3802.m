clear; close all; clc;

Aluminum25V = readmatrix("Aluminum_25V_240mA");
Aluminum30V = readmatrix("Aluminum_30V_290mA");
Brass25V = readmatrix("Brass_25V_237mA");
Brass30V = readmatrix("Brass_30V_285mA");
Steel22V = readmatrix("Steel_22V_203mA");

Aluminum25V = dataProcessor(Aluminum25V);
Aluminum30V = dataProcessor(Aluminum30V);
Brass25V = dataProcessor(Brass25V);
Brass30V = dataProcessor(Brass30V);
Steel22V = dataProcessor(Steel22V);


numchannels = width(Aluminum25V);
%measured in meters
spaceBtwn = 0.0127; 
startSpace = 0.0762; 
endSpace = 0.0254; 
bandHeater = 0.0254; 
diameter = 0.0254; 
x0 = 0.0349; 

channelLocations = zeros(1,8);

for i = 1:numchannels
    channelLocations(i) = startSpace + spaceBtwn*i;
end

plotting(Aluminum25V,x0,channelLocations,"Aluminum-25V-240mA");
plotting(Aluminum30V,x0,channelLocations,"Aluminum-30V-290mA");
plotting(Brass25V,x0,channelLocations,"Brass-25V-237mA");
plotting(Brass30V,x0,channelLocations,"Brass-30V-285mA");
plotting(Steel22V,x0,channelLocations,"Steel-22V-203mA");

%%Find H_analytical
Aluminumk = 130;
Brassk = 115;
Steelk = 16.2;

AluminumV = [25, 30];
BrassV = [25,30];
SteelV = 22;

AluminumI = [.240, .290];
BrassI = [.237,.285];
SteelI = .203;

A = .0005067;

HSteel = (SteelV .* SteelI) ./ (Steelk * A);
HAluminum = (AluminumV .* AluminumI) ./ (Aluminumk * A);
HBrass = (BrassV .* BrassI) ./ (Brassk * A);

L = .2159;
To = 14.539;
B2T0 = 13.916;
A1To = 15.983;
A2To = 15.741;
STo = 9.627;
x = linspace(0, L, 10);
xTherm = [.034925, 0.047625, 0.060325, 0.073025, 0.085725, 0.098425, 0.111125, 0.123825];



figure();
plot(x,HBrass(1) * x + To); hold on; grid on;
scatter(xTherm, HBrass(1).*xTherm + To);
title('Brass-25V-237mA');
legend('Brass at 25 V and 237mA', 'Location of thermocouples');
xlabel('Location(m)');
ylabel('Temperature(deg C)');

figure();
plot(x,HBrass(2) * x + B2T0); hold on; grid on;
scatter(xTherm, HBrass(2).*xTherm + B2T0);
title('Brass-30V-285mA');
legend('Brass at 30 V and 285mA', 'Location of thermocouples');
xlabel('Location(m)');
ylabel('Temperature(deg C)');

figure();
plot(x,HAluminum(1) * x + A1To); hold on; grid on;
scatter(xTherm, HAluminum(1).*xTherm + A1To);
title('Aluminum-25V-240mA');
legend('Aluminum at 25 V and 240mA', 'Location of thermocouples');
xlabel('Location(m)');
ylabel('Temperature(deg C)');

figure();
plot(x,HAluminum(2) * x + A2To); hold on; grid on;
scatter(xTherm, HAluminum(2).*xTherm + A2To);
title('Aluminum-30V-290mA');
legend('Aluminum at 30 V and 290mA', 'Location of thermocouples');
xlabel('Location(m)');
ylabel('Temperature(deg C)');

figure();
plot(x,HSteel * x + STo); hold on; grid on;
scatter(xTherm, HSteel.*xTherm + STo);
title('Steel-22V-203mA');
legend('Steel at 25 V and 237mA', 'Location of thermocouples');
xlabel('Location(m)');
ylabel('Temperature(deg C)');




