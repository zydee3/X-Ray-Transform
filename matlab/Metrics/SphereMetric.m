classdef sphereMetric < Metric
    
    properties
        radius
    end
    
    methods      
        function obj = sphereMetric(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
        end
        
        function out = lg(obj,x,y)
            R2 = obj.radius*obj.radius;
            out = log((4*R2*R2) ./ sumsq(x.*x,y.*y,R2));
        end
        
        function out = dxlg(obj,x,y)
            out = -4.*x./(x.*x + y.*y + obj.radius*obj.radius);
        end
        
        function out = dylg(obj,x,y)
            out = -4.*y./(x.*x + y.*y + obj.radius*obj.radius);
        end
        
        function out = curv(obj,x,y)
            out = ones(size(x))/(obj.radius*obj.radius);
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

