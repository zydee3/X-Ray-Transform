function [xO,yO,thO,lg] = geoStep(xI,yI,thI,metric,h,type)
%GEOSTEP geodesic step
% 
% inputs: - xI,yI,thI: input position 
%         - metric: structure describing the metric
%         - h: stepsize
%         - type: timestepper 'EE', 'IE' or 'RK4'
% outputs - xO,yO,thO: output position
%         - lg: values of the (conformal) metric at the current points
%
% Author: Francois Monard
% Date: 06-23-2016

switch type
    case 'EE' % Explicit Euler
        [lg,dxlg,dylg] = metric.metVal(xI, yI, metric);
        cth = cos(thI); sth = sin(thI);
        
        hh = exp(-.5*lg)*h;
        
        xO = xI + hh.*cth;
        yO = yI + hh.*sth;
        thO = thI + .5*hh.*(cth.*dylg - sth.*dxlg);
        
    case 'IE' % Improved Euler
        
        % predictor
        [lg,dxlg,dylg] = metric.metVal(xI, yI, metric);
        cth = cos(thI); sth = sin(thI);
        
        k1x = exp(-.5*lg).*cth;
        k1y = exp(-.5*lg).*sth;
        k1th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h * k1x;
        yO = yI + h * k1y;
        thO = thI + h * k1th;
        
        % corrector
        [lg,dxlg,dylg] = metric.metVal(xO, yO, metric);
        cth = cos(thO); sth = sin(thO);
        
        k2x = exp(-.5*lg).*cth;
        k2y = exp(-.5*lg).*sth;
        k2th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h * (k1x + k2x)/2;
        yO = yI + h * (k1y + k2y)/2;
        thO = thI + h * (k1th + k2th)/2;
        
    case 'RK4' % Runge-Kutta 4
        
        % first slope
        [lg,dxlg,dylg] = metric.metVal(xI, yI, metric);
        cth = cos(thI); sth = sin(thI);
        
        k1x = exp(-.5*lg).*cth;
        k1y = exp(-.5*lg).*sth;
        k1th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h/2 * k1x;
        yO = yI + h/2 * k1y;
        thO = thI + h/2 * k1th;
        
        % second slope
        [lg,dxlg,dylg] = metric.metVal(xO, yO, metric);
        cth = cos(thO); sth = sin(thO);
        
        k2x = exp(-.5*lg).*cth;
        k2y = exp(-.5*lg).*sth;
        k2th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h/2 * k2x;
        yO = yI + h/2 * k2y;
        thO = thI + h/2 * k2th;
        
        % third slope
        [lg,dxlg,dylg] = metric.metVal(xO, yO, metric);
        cth = cos(thO); sth = sin(thO);
        
        k3x = exp(-.5*lg).*cth;
        k3y = exp(-.5*lg).*sth;
        k3th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h * k3x;
        yO = yI + h * k3y;
        thO = thI + h * k3th;
        
        % fourth slope
        [lg,dxlg,dylg] = metric.metVal(xO, yO, metric);
        cth = cos(thO); sth = sin(thO);
        
        k4x = exp(-.5*lg).*cth;
        k4y = exp(-.5*lg).*sth;
        k4th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);
        
        xO = xI + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
        yO = yI + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
        thO = thI + h/6 * (k1th + 2*k2th + 2*k3th + k4th);        
        
    otherwise 
        error('wrong timestepper')
end