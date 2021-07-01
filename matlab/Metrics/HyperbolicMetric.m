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
            s = (R2 - (x.*x + y.*y));
            R2os = R2./s; %lol this saves exactly 1 flop.
            out = log(4*R2os.*R2os);
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
       
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            R2 = obj.radius*obj.radius;
            %s = (R2 - (X.*X + Y.*Y));
            
            lgt = log((4*R2*R2) ./ ((R2 - (X.*X + Y.*Y)).*(R2 - (X.*X + Y.*Y))));
            dxlgt = 4.*X./(R2 - (X.*X + Y.*Y));
            dylgt = 4.*Y./(R2 - (X.*X + Y.*Y));
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            R2 = obj.radius*obj.radius;
            s = (R2 - (X.*X + Y.*Y));
            
            lgt = log((4*R2*R2) ./ (s.*s));
            dxlgt = 4*X./s;
            dylgt = 4*Y./s;
            curvt = -ones(size(X))/R2;
        end
        
        
    end
end
