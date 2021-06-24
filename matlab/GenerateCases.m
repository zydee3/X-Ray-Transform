clc; close all;
addpath(genpath('./'))

%% Domain.build
d0 = Domain.build();

dc0 = Domain.build('circle', 'radius', 4, 'origin', [100000,100000000000000]);
dcs0 = Domain.build('cosine', 'radius', 4, 'amplitude', 1, 'cycles', 4);
dp0 = Domain.build('polygon', 'radius', 4, 'sides', 8);
dsp0 = Domain.build('smoothpoly', 'radius', 3, 'sides', 3, 'bevelRadius', 0.3);

%% Metric.build
m0 = Metric.build();

me0 = Metric.build('euclid');
mh0 = Metric.build('hyperbolic', 'radius', 4);
mp0 = Metric.build('poly', 'coeffs', [12,0.43,5,2,-1,1]);
ms0 = Metric.build('sphere', 'radius', 4);

%% Metric (do not use)

%m1 = Metric();
%m2 = Metric(@(x,y) exp(x.^2 + y.^2));
    % Same as: m2 = Metric.build('poly', 'coeffs', [1,0,1,0,0,0]);

%% RiemannSuface
rs0 = RiemannSurface();
rs1 = RiemannSurface(dc0,me0);