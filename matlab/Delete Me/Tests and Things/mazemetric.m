clc; close all; clear;

dom = polygonDomain(radius = 15*sqrt(2), sides = 4); dom.theta = pi/4; 
met = eatonMazeMetric(seed=12345647554341324356);

testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';
testsurf.stepSize = 0.02;
testsurf.geoDur = 300;

Figureizer.figure(testsurf);
Figureizer.plotGeoNormals(testsurf,0.5 ,100);
%Figureizer.plotGeoNormals(testsurf,0.5-pi/4 ,100);
Figureizer.plotGeoNormals(testsurf,0.5-pi/2 ,100);

        
        