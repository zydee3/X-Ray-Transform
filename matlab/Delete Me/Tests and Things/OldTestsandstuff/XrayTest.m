

dom = polygonDomain(radius = 7, sides = 4); dom.theta = pi/4; 
met = eatonMazeMetric();%(radius=1.7);

func = @(x,y) double(x.^2+y.^2 < 0.5^2);


ss = pi * 0.01;
surf0 = RiemannSurface(dom,met, stepType='RK4', stepSize=ss, geoDur=20);
[beta,alpha] = meshgrid(linspace(0,2*pi,100),linspace(-pi,pi,100)*0.5 *0.3);




figure, hold on
[minB,maxB] = dom.getBoundingBox();

[VX,VY] = meshgrid(minB(1):0.01:maxB(1),minB(2):0.01:maxB(2));
figgg = pcolor(VX,VY,func(VX,VY));
figgg.EdgeColor = 'none';

%surf.plot;
surf0.plotGeoFan(0);
dom.plotAlNormal;


surf0.stepType = 'RK4';
tic
   xray1 = surf0.I0(beta,alpha, func);
toc

surf0.stepType = 'IE';
tic
    xray2 = surf0.I0(beta,alpha, func);
toc

surf0.stepType = 'EE';
tic
    xray3 = surf0.I0(beta,alpha, func);
toc


figure; hold on;
pl = surf(beta,alpha,xray1);
pl.EdgeColor = 'none';
colorbar


figure; hold on;
pl = surf(beta,alpha,abs(xray2-xray3));
pl.EdgeColor = 'none';
colorbar
