dom = circleDomain(radius=1.5);
met = euclidMetric();
%{
im = InMap([0,0,0,0,0;...
            0,0,-2,0,0;...
            0,0,0,1,0;...
            0,2,0,0,0;
            0,0,0,0,0;]);
%}
func = @(x,y) (  exp(-1./(1-min((x/2.4).^2+(y/1.5-0.5).^2,1))).*(sin(x)+sin(4*y)) + exp(-1./(1-min((x).^2+(y+1.3).^2,1))).*(sin(10*x)+sin(10*y) )  );
im = InMap(func);
im = im.transform(0,0,0.5,0.5,0.3); 
func = im.eval;



ss = 0.01;
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=ss, geoDur=20);
[beta,alpha] = meshgrid(linspace(0,2*pi,200),linspace(-pi,pi,200)*0.5);



figure, hold on
[minB,maxB] = dom.getBoundingBox();

[VX,VY] = meshgrid(minB(1):0.01:maxB(1),minB(2):0.01:maxB(2));
figgg = pcolor(VX,VY,func(VX,VY));
figgg.EdgeColor = 'none';

%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;



tic
   xray = surf.I0(beta,alpha, func);
toc


figure
pl = pcolor(beta,alpha,xray);
pl.EdgeColor = 'none';



tic 
    funcBack = surf.Backproject(xray, beta,alpha, 100);
toc


figure
pl = pcolor(funcBack);
pl.EdgeColor = 'none';