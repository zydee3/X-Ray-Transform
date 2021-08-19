vals = @(x,y) 1./(x.^2+y.^2+1);
B = linspace(0,2*pi,256);
A = linspace(-pi,pi,256)*0.5;



%--------------------------------------------------------------------------

dom = cosineDomain();
met = hyperbolicMetric(radius=5);

im = InMap(vals);

surf = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=10);

[beta,alpha] = meshgrid(B,A);

tic
    I0_1 = XrayI0(im,surf,beta,alpha);
toc

%--------------------------------------------------------------------------

params.type = 'cos';
params.a = 2;
params.b = 0.111;
params.n = 4;
params.xOrig = 0;
params.yOrig = 0;
params.rmax = 2;
params.th0 = 0;

domain = build_boundary(params);
metric = build_metric_negative(5);


    I0_2 = geoI0(domain, vals, B, A, 0.01, metric, 'RK4');


%--------------------------------------------------------------------------
