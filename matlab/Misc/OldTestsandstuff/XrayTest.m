

dom = circleDomain(radius=1.5);
met = euclidMetric();%(radius=1.7);
%{
im = InMap([0,0,0,0,0;...
            0,0,-2,0,0;...
            0,0,0,1,0;...
            0,2,0,0,0;
            0,0,0,0,0;]);
%}
im = InMap(@(x,y) (exp(-1./(1-min((x/2.4).^2+(y/1.5-0.5).^2,1) )).*(sin(x)+0.5*sin(8*y))));
im = im.transform(0,0,0.5,0.5,0.3);  

ss = 0.01;
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=ss, geoDur=20);
[beta,alpha] = meshgrid(linspace(0,2*pi,250),linspace(-pi,pi,250)*0.5);
tic
   xray1 = surf.I0(beta,alpha, im);
toc

figure, hold on
[minB,maxB] = dom.getBoundingBox();
im.plot(minB,maxB);
%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;

%{.
surf.stepType = 'IE';
%surf.stepSize = 0.1;
tic
    xray2 = surf.I0(beta,alpha, im);
toc

%{
surf.stepType = 'IE';
surf.stepSize = 0.00001;
tic
    xrayComp = surf.I0(beta,alpha, im);
toc
%}
figure
pl = pcolor(beta,alpha,xrayComp);
pl.EdgeColor = 'none';

figure; hold on;
subplot(1,3,1);
pl = pcolor(beta,alpha,abs(xray1-xrayComp));
pl.EdgeColor = 'none';
title('RK4 error',['step size: ', string(ss)]) 
colorbar

subplot(1,3,2);
pl = pcolor(beta,alpha,abs(xray2-xrayComp));
pl.EdgeColor = 'none';
title('IE error') 
colorbar

subplot(1,3,3);
pl = pcolor(beta,alpha,abs(xray1-xrayComp)-abs(xray2-xrayComp));
pl.EdgeColor = 'none';
title('RK4 error - IE error') 
colorbar

%}