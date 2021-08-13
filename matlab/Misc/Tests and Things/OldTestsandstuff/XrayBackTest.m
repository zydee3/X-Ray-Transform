dom = circleDomain(radius=1.5);
met = euclidMetric();
 
func = @(x,y) (  double((x-0.7).^2+(y*2-0.5).^2 <= 0.25) - double((x+0.7).^2+(y - 0.1).^2 <= 0.25) - double((x+0.7).^2+(y - 0.5).^2 <= 0.1)  );


beta = linspace(0,2*pi,200);
alpha = linspace(-pi,pi,200)*0.5;

ss = 0.01;
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=ss, geoDur=20);

figure, hold on
[minB,maxB] = dom.getBoundingBox();

[VX,VY] = ndgrid(minB(1):0.04:maxB(1),minB(2):0.04:maxB(2));
figgg = pcolor(VX,VY,func(VX,VY));
figgg.EdgeColor = 'none';

%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;


[betaI,alphaI] = ndgrid(beta,alpha);
tic
   xrayData = surf.I0(betaI,alphaI, func);
toc

figure
pl = pcolor(betaI,alphaI,xrayData);
pl.EdgeColor = 'none';

[betaI,alphaI] = ndgrid(beta,alpha);
xray = griddedInterpolant(betaI,alphaI, xrayData);


tic 
    funcBack = surf.I0star(xray, VX,VY, 100);
toc


figure; hold on;
pl = pcolor(VX,VY,funcBack);
pl.EdgeColor = 'none';

dom.plotAlNormal;
