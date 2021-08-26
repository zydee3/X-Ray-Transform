classdef eatonGeneralMetric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        theta
    end
    
    methods      
        function obj = eatonGeneralMetric(args)
            arguments
                args.theta (1,1) {mustBeNumeric} = pi/4
            end
            
            obj.theta = args.theta;
        end
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            
            th = obj.theta;
            r = sqrt(X.*X+Y.*Y);
            p = (2*th)/(pi+th);
            
            bool = r < 1;
            r = r(bool); X = X(bool); Y = Y(bool);
            lgt(bool)   = real(log((2./r-1).^p));
            
            dxlgt(bool) = real(p*2*X./((r-2).*r.^2));
            dylgt(bool) = real(p*2*Y./((r-2).*r.^2));            
        end        
        
        
    end
end

