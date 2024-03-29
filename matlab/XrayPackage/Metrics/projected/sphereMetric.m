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
            s = (R2 + x.*x + y.*y);
            out = log(4 * R2 * R2 / (s.*s));
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
        
       
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            R2 = obj.radius*obj.radius;

            lgt = log((4*R2*R2) ./ (exp(2*log(R2 + X.*X + Y.*Y)) ));
            dxlgt = -4*X./(R2 + X.*X + Y.*Y);
            dylgt = -4*Y./(R2 + X.*X + Y.*Y);
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            R2 = obj.radius*obj.radius;
            %s = (R2 + X.*X + Y.*Y);
            
            lgt = log((4*R2*R2) ./ (exp(2*log(R2 + X.*X + Y.*Y)) ));
            dxlgt = -4*X./(R2 + X.*X + Y.*Y);
            dylgt = -4*Y./(R2 + X.*X + Y.*Y);
            curvt = ones(size(X))/R2;
        end
        
        
        function [xO,yO,zO] = deproject(obj,X,Y)
            r = obj.radius;
            s = 1./(r*r+X.*X+Y.*Y);
            xO = 2*r*r*X.*s;
            yO = 2*r*r*Y.*s;
            zO = r*(r*r-X.*X-Y.*Y).*s;
        end   
        
    end
end

