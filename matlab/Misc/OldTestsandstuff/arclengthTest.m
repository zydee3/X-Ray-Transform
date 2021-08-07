dom = polygonDomain();
testsurf = RiemannSurface(dom);

circ = testsurf.arcLength(0,2*pi, 0.00001);

betas = testsurf.StoBeta(linspace(0,circ,50), 0.0005);

figure; hold on;
testsurf.plot;
dom.plotBdrPoint(betas);

S = testsurf.BetatoS(betas, 0.00001);


% plot errors |S - F(beta)|
figure
plot(abs(betas-S))

% plot difference between points
figure; hold on;
plot(abs(S(1:end-1) - S(2:end)))
[xs,ys,~] = testsurf.BAtoXYTh(betas,0);
plot(sqrt( (xs(1:end-1)-xs(2:end)).^2 + (ys(1:end-1)-ys(2:end)).^2 ) , 'r')

% plot difference between points using euclidean distance
figure
