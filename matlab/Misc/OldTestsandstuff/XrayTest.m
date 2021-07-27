

dom = circleDomain(radius=1.5);
met = euclidMetric();%(radius=1.7);
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



ss = 0.2;
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=ss, geoDur=20);
[beta,alpha] = meshgrid(linspace(0,2*pi,250),linspace(-pi,pi,250)*0.5);




figure, hold on
[minB,maxB] = dom.getBoundingBox();

[VX,VY] = meshgrid(minB(1):0.01:maxB(1),minB(2):0.01:maxB(2));
figgg = pcolor(VX,VY,func(VX,VY));
figgg.EdgeColor = 'none';

%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;



tic
   xray1 = surf.I0(beta,alpha, func);
toc

surf.stepType = 'IE';
tic
    xray2 = surf.I0(beta,alpha, func);
toc

surf.stepType = 'EE';
tic
    xray3 = surf.I0(beta,alpha, func);
toc


%{
surf.stepType = 'EE';
surf.stepSize = 0.00001;
tic
    xrayComp = surf.I0(beta,alpha, im.eval);
toc
%}

figure
pl = pcolor(beta,alpha,xray2);
pl.EdgeColor = 'none';


figure; hold on;
subplot(2,2,1);
pl = pcolor(beta,alpha,abs(xray1-xrayComp));
pl.EdgeColor = 'none';
title('RK4 error',['step size: ', string(ss)]) 
colorbar

subplot(2,2,2);
pl = pcolor(beta,alpha,abs(xray2-xrayComp));
pl.EdgeColor = 'none';
title('IE error') 
colorbar

subplot(2,2,3);
pl = pcolor(beta,alpha,abs(xray3-xrayComp));
pl.EdgeColor = 'none';
title('EE error') 
colorbar


subplot(2,2,4);
pl = pcolor(beta,alpha,abs(xray3-xrayComp)-abs(xray2-xrayComp));
pl.EdgeColor = 'none';
title('EE error - IE error') 
colorbar

%}