clc; clear; close all;
addpath(genpath('./'))

% test cmnds from the cmnd bar ig

m0 = Metric(@(x,y) x.^2 + y.^2);
m1 = EuclidMetric();
R = 2;
m2 = Metric(@(x,y) (4*R^4)./(R^2 - (x.^2 + y.^2)).^2);
m3 = Metric();

d0 = CircleDomain(1);
d1 = d0.transform(2,0.5,0.5);
d1 = CosineDomain();
d1.plot

r0 = RiemannSurface(m0,d1);
r1 = RiemannSurface(d0,m2);

%{
[X,Y] = meshgrid(-5:0.1:5);
tic
m1.metricVals(X,Y);
toc
tic
m3.metricVals(X,Y);
toc
%}
%{
figure; hold on
fplot(@(x) -cos(x))
fplot(@(x) dderiv(@(z) cos(z),x))
%}