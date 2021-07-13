

dom = circleDomain(radius=2);
met = hyperbolicMetric(radius=5);

im = InMap([0,0,0,0,0;...
            0,0,-2,0,0;...
            0,0,0,1,0;...
            0,2,0,0,0;
            0,0,0,0,0;]);
    im = im.transform(0,0,0.7,0.7,0);       
    
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=20);

[beta,alpha] = meshgrid(linspace(0,2*pi,256),linspace(-pi,pi,256)*0.5);

tic
xray = XrayI0(im,surf,beta,alpha);
toc



figure, hold on
subplot(1,2,1);
im.plot;
%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;


subplot(1,2,2);
pl = pcolor(beta,alpha,xray);
pl.EdgeColor = 'none';