clc; close all;
addpath(genpath('./'))

%% Subclasses
cd0 = CircleDomain();
cd1 = CircleDomain(123.213432);
%cd2 = CircleDomain('wadwafa'); 

cod0 = CosineDomain();
cod1 = CosineDomain(123.22);
cod2 = CosineDomain(2.22,4.22);
cod3 = CosineDomain(1.22,2.22,3);
%cod4 = CircleDomain('wadwafa'); 

pd0 = PolygonDomain();

spd0 = SmoothPolygonDomain();
spd1 = SmoothPolygonDomain(2,3,0.2);


%% User-implemented classes


%% Main class, fallbacks
d0 = Domain();

%d4 = Domain(@(),456); % fallack currently expects user to imput max radius, find a nice way to compute this