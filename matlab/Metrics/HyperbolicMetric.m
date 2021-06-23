classdef HyperbolicMetric < Metric
    
    properties
        radius = 2
    end
    
    methods
        function obj = HyperbolicMetric(radius)
            if (nargin == 0), radius = obj.radius; 
            elseif(~isnumeric(radius))
                error("Bad input arguments.")
            end
            
            obj.radius = radius;
            
            R2 = radius*radius;
            obj.lg   = @(x,y) log((4*R2*R2) ./ diffsq(x.*x,y.*y,R2)); % @(x,y) log(4*R^4) - 2*log(x.*x + y.*y + R*R);
            obj.dxlg = @(x,y) 4.*x./(R2 - x.*x - y.*y);
            obj.dylg = @(x,y) 4.*y./(R2 - x.*x - y.*y);
            obj.curv = @(x,y) -ones(size(x))/R2;
        end
        
        %{
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            its fine for now ig
        end
        %}
        
    end
end

function out = diffsq(x,y,z) % minor optimization 
    s = (z - (x+y));
    out = s.*s;
end

