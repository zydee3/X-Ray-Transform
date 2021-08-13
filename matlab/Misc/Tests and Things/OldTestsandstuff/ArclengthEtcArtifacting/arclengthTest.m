dom = polygonDomain();
%dom = cosineDomain();
testsurf = RiemannSurface(dom);
testsurf.stepType = 'RK4'

circ = testsurf.arcLength(0,pi, 0.01);

sI = linspace(0,circ,20);
betas = testsurf.StoBeta(sI(1:end), 0.01);

figure; hold on;
testsurf.plot;
dom.plotBdrPoint(betas);

S = testsurf.BetatoS(betas, 0.01);


% plot errors |S - F(beta)|
figure
plot(abs(betas-S))

% plot difference between points
figure; hold on;
plot(abs(S(1:end-1) - S(2:end)))
[xs,ys,~] = testsurf.BAtoXYTh(betas,0);
%plot(sqrt( (xs(1:end-1)-xs(2:end)).^2 + (ys(1:end-1)-ys(2:end)).^2 ) , 'r')
yline(sI(2)-sI(1))
