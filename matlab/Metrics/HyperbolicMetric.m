classdef hyperbolicMetric < Metric
    
    properties
        radius
    end
    
    methods
        function obj = hyperbolicMetric(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
        end
        
        
        function out = lg(obj,x,y)
            R2 = obj.radius*obj.radius;
            out = log((4*R2*R2) ./ diffsq(x.*x,y.*y,R2));
        end
        
        function out = dxlg(obj,x,y)
            out = 4.*x./(obj.radius*obj.radius - x.*x - y.*y);
        end
        
        function out = dylg(obj,x,y)
            out = 4.*y./(obj.radius*obj.radius - x.*x - y.*y);
        end
        
        function out = curv(obj,x,y)
            out = -ones(size(x))/(obj.radius*obj.radius);
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

