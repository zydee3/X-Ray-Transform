classdef SphereMetric < Metric
    
    properties
        radius = 2
    end
    
    methods
        function obj = SphereMetric(radius)
            if (nargin == 0), radius = obj.radius; 
            elseif(~isnumeric(radius))
                error("Bad input arguments.")
            end
            
            obj.radius = radius;
            
            R2 = radius*radius;
            obj.lg   = @(x,y) log((4*R2*R2) ./ sumsq(x.*x,y.*y,R2)); % @(x,y) log(4*R^4) - 2*log(x.*x + y.*y + R*R);
            obj.dxlg = @(x,y) -4.*x./(x.*x + y.*y + R2);
            obj.dylg = @(x,y) -4.*y./(x.*x + y.*y + R2);
            obj.curv = @(x,y) ones(size(x))/R2;
        end
        
        %{
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            its fine for now ig
        end
        %}
        
    end
end

function out = sumsq(x,y,z) % minor optimization 
    s = (x + y + z);
    out = s.*s;
end

