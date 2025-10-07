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
spaceBtwn = 0.0127; %m
startSpace = 0.0762; %m
endSpace = 0.0254; %m
bandHeater = 0.0254; %m
diameter = 0.0254; %m
x0 = 0.0349; %m

channelLocations = zeros(1,7);

for i = 1:numchannels-1
    channelLocations(i) = startSpace + spaceBtwn*i;
end

plotting(Aluminum25V,x0,channelLocations,"Aluminum-25V-240mA");
plotting(Aluminum30V,x0,channelLocations,"Aluminum-30V-290mA");
plotting(Brass25V,x0,channelLocations,"Brass-25V-237mA");
plotting(Brass30V,x0,channelLocations,"Brass-30V-285mA");
plotting(Steel22V,x0,channelLocations,"Steel-22V-203mA");






