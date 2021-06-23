clc; close all;
addpath(genpath('./'))

%% Subclasses
cd0 = CircleDomain();
cd1 = CircleDomain(123.213432);
%cd2 = CircleDomain('wadwafa'); % TODO: Bad parameter doesn't error

cod0 = CosineDomain();
cod1 = CosineDomain(123.22);
cod2 = CosineDomain(2.22,4.22);
cod3 = CosineDomain(1.22,2.22,3.22); % -- Correctly rounds down 3rd parameter
cod4 = CircleDomain('wadwafa'); % TODO: Bad parameter doesn't error

%% User-implemented classes


%% Main class, fallbacks
d0 = Domain();

d1 = Domain(1.23); % Transform params
d2 = Domain(1.23,4.56);
d3 = Domain(1.23,4.56,7.89);

d4 = Domain(@(),456); % fallack currently expects user to imput max radius, find a nice way to compute this