function [lgt,dxlgt,dylgt] = metricValues4(x, y, met)

R2 = (met.R)^2;

lgt = log(4*R2*R2) - 2*log(R2 - x.*x - y.*y);

dxlgt = 4.*x./(R2 - x.*x - y.*y);
dylgt = 4.*y./(R2 - x.*x - y.*y);

end