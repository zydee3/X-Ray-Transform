function [metric,lg,dxlg,dylg,curv] = build_metric_negative(R)
%BUILD_METRIC_NEGATIVE builds symbolic expression for an isotropic metric 
% of constant negative curvature
% 
% Author: Francois Monard
% Date: 2-28-2015

type = 4;

if nargin == 0,
    R = 1;
end

%%%%% define the metric
metric.type = type;
metric.metVal = str2func(['metricValues',num2str(type)]);
metric.metValCurv = str2func(['metricValuesCurvature',num2str(type)]);
metric.R = R;


lg   = @(x,y) log(4*R^4) - 2*log(R*R - x.*x - y.*y);
dxlg = @(x,y) 4.*x./(R*R - x.*x - y.*y);
dylg = @(x,y) 4.*y./(R*R - x.*x - y.*y);
curv = @(x,y) -1/R^2;

metric.lg = lg;
metric.dxlg = dxlg;
metric.dylg = dylg;
metric.curv = curv;

end