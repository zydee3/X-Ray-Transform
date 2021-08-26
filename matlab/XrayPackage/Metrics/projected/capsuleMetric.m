classdef capsuleMetric < Metric
    
    properties
        radius
        height
    end
    
    methods      
        function obj = capsuleMetric(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
        end
           
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            R2 = 1;
            m2 = (X.*X + Y.*Y);
            bool = m2 > R2;
            
            lgt(bool)   = -log(m2(bool));
            dxlgt(bool) = -2*X(bool)./m2(bool);
            dylgt(bool) = -2*Y(bool)./m2(bool);
            
            lgt(~bool) = log((4*R2*R2) ./ (exp(2*log(R2 + m2(~bool))) ));
            dxlgt(~bool) = -4*X(~bool)./(R2 + m2(~bool));
            dylgt(~bool) = -4*Y(~bool)./(R2 + m2(~bool));
        end
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            curvt = lgt;
            R = 1;
            m2 = (X.*X + Y.*Y);
            bool = m2 > R2;
            
            lgt   = -log(m2(bool));
            dxlgt = -2*X(bool)./m2(bool);
            dylgt = -2*Y(bool)./m2(bool);
            
            lgt(~bool) = log((4*R2*R2) ./ (exp(2*log(R2 + m2(~bool))) ));
            dxlgt(~bool) = -4*X(~bool)./(R2 + m2(~bool));
            dylgt(~bool) = -4*Y(~bool)./(R2 + m2(~bool));
            curvt(~bool) = ones(size(X(~bool)))/R2;
        end
        
        function [xO,yO,zO] = deproject(obj,X,Y)
            m = (X.*X + Y.*Y);
            
            R = 1;
            bool = m > R*R;
            
            xO = zeros(size(X));
            yO = xO;
            zO = xO;
            
            xO(bool) = R*X(bool)./sqrt(m(bool));
            yO(bool) = R*Y(bool)./sqrt(m(bool));
            zO(bool) = -log(sqrt(m(bool)));
            
            s = 1./(R*R+m(~bool));
            xO(~bool) = 2*R*R*X(~bool).*s;
            yO(~bool) = 2*R*R*Y(~bool).*s;
            zO(~bool) = R*(R*R-m(~bool)).*s;
        end   

        
    end
end

