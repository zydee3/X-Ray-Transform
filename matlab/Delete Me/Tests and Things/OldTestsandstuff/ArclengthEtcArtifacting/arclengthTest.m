dom = polygonDomain();
met = gaussMetric(weights=3,widths=0.3,xOffsets=-1);
%dom = cosineDomain();
testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';

circ = testsurf.BetatoS(3*pi/2, 0.01);

sI = linspace(0,circ,50);
betas = testsurf.StoBeta(sI(1:end), 0.01);

figure; hold on;
testsurf.plot;
dom.plotBdrPoint(betas);
dom.plotOrigin;

S = testsurf.BetatoS(betas, 0.001);

%{
% plot errors |S - F(beta)|
figure
plot(abs(betas-S))

% plot difference between points
figure; hold on;
plot(abs(S(1:end-1) - S(2:end)))
[xs,ys,~] = testsurf.BAtoXYTh(betas,0);
%plot(sqrt( (xs(1:end-1)-xs(2:end)).^2 + (ys(1:end-1)-ys(2:end)).^2 ) , 'r')
yline(sI(2)-sI(1))
%}
