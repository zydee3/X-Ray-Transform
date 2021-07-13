dom = circleDomain(radius=2);
met = sphereMetric(radius=4);
surf = RiemannSurface(dom,met, stepType='EE', stepSize=0.01, geoDur=20);

xI = rand(1,65536);
yI = rand(1,65536);
thI = rand(1,65536);

tic
for (i = 1:100)
    surf.geoStep(xI,yI,thI);
end    
time = toc;
time = time/100