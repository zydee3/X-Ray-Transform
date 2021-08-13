randElemX = rand(1,100000);
randElemY = rand(1,100000);

%{
tic
me0.metricValsCurv(randElemX, randElemY);
toc
tic
mh0.metricValsCurv(randElemX, randElemY);
toc
tic
mp0.metricValsCurv(randElemX, randElemY);
toc
tic
ms0.metricValsCurv(randElemX, randElemY);
toc
%}


xI = ones(1,50000);
yI = ones(1,50000);
thI = rand(1,50000);


met = hyperbolicMetric('radius',5);
domain = cosineDomain('radius',3,'amplitude',1);
domain.theta = 2;
riemsuf = RiemannSurface(domain, met);
riemsuf.stepSize = 0.005;
riemsuf.stepType = 'RK4';
riemsuf.geoDur = 20;


met = build_metric_negative(5);
params.type = 'cos';
params.a = 3;
params.b = 1;
params.xOrig = 0;
params.yOrig = 0;
params.th0 = 2;
domain = build_boundary(params);


tic
geodesics(domain, xI, yI, thI, 0.005, met, 'RK4');
toc
tic
riemsuf.geodesic(xI, yI, thI);
toc
tic
geodesics(domain, xI, yI, thI, 0.005, met, 'RK4');
toc
tic
riemsuf.geodesic(xI, yI, thI);
toc
tic
geodesics(domain, xI, yI, thI, 0.005, met, 'RK4');
toc
tic
riemsuf.geodesic(xI, yI, thI);
toc
