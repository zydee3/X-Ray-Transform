%initialize surface
testdom = circleDomain(radius = 3);
testmetric = euclidMetric();

testsurf = RiemannSurface(testdom,testmetric, ...
                stepType = 'RK4', stepSize = 0.05);


%just using a cheap vanishing function (defined at the bottom of this script)
beta = linspace(0,2*pi,200);
alpha = linspace(-pi,pi,200)*0.5;
[betaI,alphaI] = ndgrid(beta,alpha);

i0Data = testsurf.I0(betaI,alphaI, @(x,y) f(x-0.5,y));
i1Data = testsurf.I1(betaI,alphaI, @(x,y) dxf(x-0.5,y), @(x,y) dyf(x-0.5,y));


% figures
figure, hold on
[minB,maxB] = testdom.getBoundingBox();

[VX,VY] = ndgrid(minB(1):0.01:maxB(1),minB(2):0.01:maxB(2));
pl = pcolor(VX,VY,f(VX-0.5,VY));
pl.EdgeColor = 'none';

%surf.plot;
testsurf.plotGeoFan(0);
testdom.plotAlNormal;

figure;
pl = pcolor(betaI,alphaI,i0Data);
pl.EdgeColor = 'none';

figure;
pl = pcolor(betaI,alphaI,i1Data);
pl.EdgeColor = 'none';



function out = dyf(X,Y)
    out = (4*Y.*(X.^2+Y.^2-1)).*ind(X,Y);
end

function out = dxf(X,Y)
    out = (4*X.*(X.^2+Y.^2-1)).*ind(X,Y);
end

function out = f(X,Y)
    out = ((X.^2+Y.^2-1).^2).*ind(X,Y);
end

function out = ind(X,Y)
    out = X.^2+Y.^2 <= 1;
end