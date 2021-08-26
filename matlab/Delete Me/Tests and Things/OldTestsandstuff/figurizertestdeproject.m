close all; clc; clear; clear classes

%dom = polygonDomain(radius = 2, sides = 4);
%dom = cosineDomain(radius = 4);
dom = circleDomain(radius = 15);
    %dom.originY = 3;
%met = sphereMetric(radius=1);
met = cylinderMetric(radius=2);
%met = capsuleMetric(radius=2);
%met = euclidMetric;

testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';
%testsurf.stepSize = pi/4;
testsurf.geoDur = 7*pi;

%func = @(x,y) zeros(size(x));
func = @(x,y) mod(floor(x*4)+floor(y*4),2);
    
sett = Figureizer.settings
%sett.penColor = [0,0,0]
%sett.funcColormap = summer


Figureizer.figureDeprojected(testsurf,func);
%Figureizer.plotGeoFanDeprojected(testsurf, -2,20);
Figureizer.plotGeoNormalsDeprojected(testsurf, 0.8,20);
%Figureizer.plotGeoDeprojected(testsurf, 0,-0.5,pi-0.5);


%abeta = testsurf.StoBeta(linspace(0,40,200),0.0001);
%Figureizer.plotBdrPointDeprojected(testsurf,abeta);


Figureizer.figure(testsurf,func);
%Figureizer.plotGeoFan(testsurf,-2,20);
Figureizer.plotGeoNormals(testsurf, 0.8,20);
%testsurf.figureJacobiRadiate(0, 1, linspace(-pi,pi, 100), enablePlotConjugates = true, enableClamped = false, enableAbsed = true );

%testsurf.plotGeo(0,1,pi+0.1);
