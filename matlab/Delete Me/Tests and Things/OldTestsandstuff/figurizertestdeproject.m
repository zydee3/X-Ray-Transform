close all; clc; clear; %clear classes

%dom = aabDomain(minB = [-1,-1],maxB = [5,1]);   %dom.theta = 1
%dom = unboundedCircleDomain(radius = 1);
%dom = cosineDomain(radius = 1);
dom = circleDomain(radius = 10);
    %dom.originX = -0.3;
    dom.originY = 10;
    
met = sphereMetric(radius=1);
%met = cylinderMetric;
%met = capsuleMetric;
%met = euclidMetric;

testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';
%testsurf.stepSize = pi/4;
testsurf.geoDur = 4*pi + testsurf.stepSize;

%func = @(x,y) zeros(size(x));
func = @(x,y) mod(floor(x*4)+floor(y*4),2) + 4*double((~dom.isInside(x,y)));
%func = @(x,y) sqrt(0.3.^2 - min(x.^2 + y.^2,0.3.^2));
    
sett = Figureizer.settings;
sett.domainThickness = 1;
sett.domainResolution = 200;
sett.gridResolution = 200;
%sett.penColor = [0,0,0]
%sett.domainColor = uint8([0,50,020])
%sett.funcColormap = summer


%sett.betaSpacing = 'default';



sett.betaSpacing = 'euclid';
Figureizer.figure(testsurf,func);
axis xy off
%Figureizer.plotGeoFan(testsurf,-2,20);
sett.penColor = uint8([256,0,0]);
Figureizer.plotGeoNormals(testsurf, 0.3,50);
%sett.penColor = uint8([256,100,200]);
Figureizer.plotGeoNormals(testsurf, -1,50);

%testsurf.figureJacobiRadiate(0, 1, linspace(-pi,pi, 100), enablePlotConjugates = true, enableClamped = false, enableAbsed = true );

%testsurf.plotGeo(0,1,pi+0.1);



sett.betaSpacing = 'arclength';
Figureizer.figureDeprojected(testsurf,func);
axis xy off
%Figureizer.plotGeoFanDeprojected(testsurf, -2,20);
sett.penColor = uint8([256,0,0]);
Figureizer.plotGeoNormalsDeprojected(testsurf, 0.3,50);
%sett.penColor = uint8([256,100,200]);
Figureizer.plotGeoNormalsDeprojected(testsurf, -1,50);
%Figureizer.plotGeoDeprojected(testsurf, 0,-0.5,pi-0.5);


%abeta = testsurf.StoBeta(linspace(0,40,200),0.0001);
%Figureizer.plotBdrPointDeprojected(testsurf,abeta);