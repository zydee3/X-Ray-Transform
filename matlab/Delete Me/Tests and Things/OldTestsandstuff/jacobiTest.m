dom = circleDomain(radius = 1);
met = gaussMetric(weights=  [2],...
                  widths=   [0.2],...
                  xOffsets= [0],...
                  yOffsets= [0]);
        
%dom = cosineDomain();
testsurf = RiemannSurface(dom,met);
testsurf.stepType = 'RK4';

figure; hold on;
testsurf.figureJacobiRadiate(-1, 0,...
                              linspace(-pi/8,pi/8, 400), enablePlotConjugates = false, enableClamped = false, enableAbsed = true );
met.plot(0.1);
%{
testsurf.plotConjugates(ones(1,50) * -1,...
                        ones(1,50) * 0,...
                        linspace(-pi/4,pi/4, 50));
%}

