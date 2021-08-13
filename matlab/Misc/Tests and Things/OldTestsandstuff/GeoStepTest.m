dom = circleDomain(radius=2);
met = sphereMetric(radius=4);
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=20);

tests = 500;
size = [256,256];

time1 = zeros(1,tests);
time2 = zeros(1,tests);

for (i = 1:tests)
    xI = rand(size);
    yI = rand(size);
    thI = rand(size);
tic    
    surf.geoStep(xI,yI,thI);
time1(i) = toc;
end   
%{
surf = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=20);

for (i = 1:tests)
    xI = rand(size);
    yI = rand(size);
    thI = rand(size);
tic    
    surf.geoStep(xI,yI,thI);
time2(i) = toc;
end   
%}
figure; hold on
plot(time1,0,'r.')
plot(time2,0.1,'b.')
ylim([-0.1,0.2])