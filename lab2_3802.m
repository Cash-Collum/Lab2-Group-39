clear; close all; clc;

Aluminum25V = readmatrix("Aluminum_25V_240mA");
Aluminum30V = readmatrix("Aluminum_30V_290mA");
Brass25V = readmatrix("Brass_25V_237mA");
Brass30V = readmatrix("Brass_30V_285mA");
Steel22V = readmatrix("Steel_22V_203mA");

Aluminum25V = Aluminum25V(~any(isnan(Aluminum25V), 2), :);
Aluminum30V = Aluminum30V(~any(isnan(Aluminum25V), 2), :);
Brass25V = Brass25V(~any(isnan(Brass25V), 2), :);
Brass30V = Brass30V(~any(isnan(Brass30V), 2), :);
Steel22V = Steel22V(~any(isnan(Steel22V), 2), :);


numchannels = width(Aluminum25V);
spaceBtwn = 0.0127; %m
startSpace = 0.0762; %m
endSpace = 0.0254; %m
bandHeater = 0.0254; %m
diameter = 0.0254; %m
x0 = 0.0349; %m

channelLocations = zeros(1,8);

for i = 1:numchannels-1
    channelLocations(i) = startSpace + spaceBtwn*i;
end

plotting(Aluminum25V,x0,channelLocations);
plotting(Aluminum30V,x0,channelLocations);
plotting(Brass25V,x0,channelLocations);
plotting(Brass30V,x0,channelLocations);
plotting(Steel22V,x0,channelLocations);






