%% initialization

% variables 

r = 0.5;
R = 1;

stepSizes = linspace(0.0001,0.1,1000);

% begin by initializing a surface
dc = circleDomain(radius = R);
ms = euclidMetric;
rsurf = RiemannSurface(dc,ms);

% initialize some betas an alphas to use
    beta = 0;%[0,3];
    alpha = [0];%linspace(-pi,pi,20)*0.5;
    [betaG,alphaG] = ndgrid(beta,alpha); % alternatively use meshgrid (doesnt work with griddedinterpolant or scatteredinterpolant)

% plot function under domain
    numsteps = 5
    rsurf.stepSize = stepSizes(end);
    rsurf.geoDur = stepSizes(end)*numsteps;

    figure, hold on, axis equal
        [minB,maxB] = dc.getBoundingBox();
        %{
        [VX,VY] = ndgrid(dc.aabbgrid(250));
        pl = pcolor(VX,VY,func0(VX,VY));
        pl.EdgeColor = 'none';
        %}
        [plotX,plotY,plotTh] = rsurf.BAtoXYTh(betaG,alphaG);
        rsurf.plotGeo(plotX,plotY,plotTh);
        dc.plotAlNormal;
     
 eeVals = zeros(1); 
 ieVals = zeros(1); 
 rkVals = zeros(1); 
 oVals = zeros(1); 
                                     
% run test
for h = stepSizes
    rsurf.stepSize = h;
    rsurf.geoDur = h*numsteps;
    %disp(h)
    
    %smoot
    func0 = @(x,y) sinh(x+3);
    xray0 = @(beta,alpha) (cosh(1+3)-cosh((1-h*numsteps)+3)); 
    %dicontinuou
    %func0 = @(x,y) double(x<0.9);
    %xray0 = @(beta,alpha) max(rsurf.geoDur-0.1,0); 
    %bumpy
    %func0 = @(x,y) double(x<0.9);
    %xray0 = @(beta,alpha) max(rsurf.geoDur-0.1,0); 
    
  
%% EE        
        rsurf.stepType = 'EE';

        i0Data = rsurf.I0(betaG,alphaG, func0);
        %i0Data = (func0(1,0))*h;

        M = abs(xray0(betaG,alphaG) - i0Data);
        eeVals = [eeVals,M];
        
%% IE           
        rsurf.stepType = 'IE';

        i0Data = rsurf.I0(betaG,alphaG, func0);
        %i0Data = (func0(1,0)+func0(1-h,0))*h/2;

        M = abs(xray0(betaG,alphaG) - i0Data);
        ieVals = [ieVals,M];
        
%% RK4       
    rsurf.stepType = 'RK4';

        i0Data = rsurf.I0(betaG,alphaG, func0);
        %i0Data = (func0(1,0)+4*func0(1-h/2,0)+func0(1-h,0))*h/6;


        M = abs(xray0(betaG,alphaG) - i0Data);
        rkVals = [rkVals,M];

    
end

figure, hold on
axis equal
plot(log(stepSizes),log(eeVals(2:end)),'r*')   
plot(log(stepSizes),log(ieVals(2:end)),'g*')
plot(log(stepSizes),log(rkVals(2:end)),'b*')
%title('dont worry, its a loglog, but i couldnt figure out how to get hold on to work with that')
     


        