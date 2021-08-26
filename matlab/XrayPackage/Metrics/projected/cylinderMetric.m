classdef cylinderMetric < Metric
    
    properties
        radius
    end
    
    methods      
        function obj = cylinderMetric(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
        end
           
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            m2 = (X.*X + Y.*Y);
            lgt   = -log(m2);
            dxlgt = -2*X./m2;
            dylgt = -2*Y./m2;
        end
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            m2 = (X.*X + Y.*Y);
            lgt   = -log(m2);
            dxlgt = -2*X./m2;
            dylgt = -2*Y./m2;
            curvt = zeros(size(X));
        end
        
        function [xO,yO,zO] = deproject(obj,X,Y)
            m = sqrt(X.*X + Y.*Y);
            xO = X./m;
            zO = Y./m;
            yO = log(m);
        end   
        
        function [xO,yO] = aabbspace(obj,domain,resoX,resoY)
            arguments
                obj
                domain (1,1) {Domain.mustBeDomain}
                resoX (1,1) {mustBeNumeric} = 250
                resoY (1,1) {mustBeNumeric} = resoX
            end
                        
            [minB,maxB] = domain.getBoundingBox; 
            origin = [domain.originX,domain.originY];
            minB = logTemp(minB+origin);   maxB = logTemp(maxB+origin);
            xO = invlogTemp( (linspace(minB(1),maxB(1),resoX)));
            yO = invlogTemp( (linspace(minB(2),maxB(2),resoY)));
            %[xO,yO] = meshgrid(xO,yO);
            
            % "signed log"
            function out = logTemp(V)
                R = 0.00001; % in a way, determines the length of the cylinder rendered
                out = zeros(size(V));
                bool = V<0;
                out(bool) = log(R./(-V(bool)+R));
                out(~bool) = log((V(~bool)+R)/R);
            end
            
            function out = invlogTemp(V)
                R = 0.00001; % in a way, determines the length of the cylinder rendered
                out = zeros(size(V));
                bool = V<0;
                out(bool) = -exp(-V(bool))+1;
                out(~bool) = exp(V(~bool))-1;
                out = out * R;
            end
        end 
        
        
    end
end

