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


[AlSteady1,AlInit1,Alx11,Aly11,Alx21,Aly21] = experimental(Aluminum25V,x0,channelLocations,"Aluminum-25V-240mA");
[AlSteady2,AlInit2,Alx12,Aly12,Alx22,Aly22] = experimental(Aluminum30V,x0,channelLocations,"Aluminum-30V-290mA");
[BSteady1,BInit1,Bx11,By11,Bx21,By21] = experimental(Brass25V,x0,channelLocations,"Brass-25V-237mA");
[BSteady2,BInit2,Bx12,By12,Bx22,By22] = experimental(Brass30V,x0,channelLocations,"Brass-30V-285mA");
[SSteady,SInit,Sx1,Sy1,Sx2,Sy2] = experimental(Steel22V,x0,channelLocations,"Steel-22V-203mA");

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
analytical.x = linspace(0, L, 10);
analytical.xTherm = [.034925, 0.047625, 0.060325, 0.073025, 0.085725, 0.098425, 0.111125, 0.123825];

Linearbrass1 = HBrass(1) * analytical.x + To;
Rawbrass1 = HBrass(1).*analytical.xTherm + To;
Linearbrass2 = HBrass(2) * analytical.x + B2T0;
Rawbrass2 = HBrass(2).*analytical.xTherm + B2T0;
LinearAluminum1 = HAluminum(1) * analytical.x + A1To;
RawAluminum1 = HAluminum(1).*analytical.xTherm + A1To;
LinearAluminum2 = HAluminum(2) * analytical.x + A2To;
RawAluminum2 = HAluminum(2).*analytical.xTherm + A2To;
LinearSteel1 = HSteel * analytical.x + STo;
RawSteel1 = HSteel.*analytical.xTherm + STo;

plotting(AlSteady1,AlInit1,channelLocations,Alx11,Aly11,Alx21,Aly21,analytical,LinearAluminum1,RawAluminum1,"Aluminum-25V-240mA");
plotting(AlSteady2,AlInit2,channelLocations,Alx12,Aly12,Alx22,Aly22,analytical,LinearAluminum2,RawAluminum2,"Aluminum-30V-290mA");
plotting(BSteady1,BInit1,channelLocations,Bx11,By11,Bx21,By21,analytical,Linearbrass1,Rawbrass1,"Brass-25V-237mA");
plotting(BSteady2,BInit2,channelLocations,Bx12,By12,Bx22,By22,analytical,Linearbrass2,Rawbrass2,"Brass-30V-285mA");
plotting(SSteady,SInit,channelLocations,Sx1,Sy1,Sx2,Sy2,analytical,LinearSteel1,RawSteel1,"Steel-22V-203mA");

