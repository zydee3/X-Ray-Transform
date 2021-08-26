clc; close all; clear;

dom = polygonDomain(radius = 15*sqrt(2), sides = 4); dom.theta = pi/4; 
%dom = cosineDomain(radius = 4);
%dom = circleDomain(radius = 1.4);

   %dom.originY = 0.5;

%met = sphereMetric(radius=1);
%met = cylinderMetric(radius=2);
%met = euclidMetric;
%met = eaton180Metric;
%met = eaton360Metric;
met = eatonMazeMetric(seed=12345647554341324356);

testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';
testsurf.stepSize = 0.02;
testsurf.geoDur = 300;

func = @(x,y) zeros(size(x));

Figureizer.figure(testsurf,func);
%[minB,maxB] = dom.getBoundingBox;
%met.plot(100,minB,maxB);
%Figureizer.plotGeoParallels(testsurf,0.7 ,linspace(0,pi,160)+pi/2+pi/4);
Figureizer.plotGeoNormals(testsurf,0.5 ,300);
%Figureizer.plotGeoFan(testsurf,0,20)
%Figureizer.plotGeo(testsurf,0.1,1.1,-pi/2);