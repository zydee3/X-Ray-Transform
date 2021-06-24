classdef ellipseDomain < Domain
    
    properties
        radiusA
        radiusB
    end
    
    methods
        function obj = ellipseDomain(args)
            arguments
                args.radiusA (1,1) {mustBeNumeric} = 3;
                args.radiusB (1,1) {mustBeNumeric} = 1;
            end
            
            obj.radiusA = args.radiusA;
            obj.radiusB = args.radiusB;
            obj.rMax = max(obj.radiusA, obj.radiusB); 
        end    
        
        function out = bdr(obj,th)
            a = obj.radiusA;
            b = obj.radiusB;
            asin = a*sin(th);
            bcos = b*cos(th);
            out = a*b ./ sqrt( bcos.*bcos + asin.*asin );
        end

        function out = dbdr(obj,th)
            a = obj.radiusA;
            b = obj.radiusB;
            out = obj.bdr(th).^3 .* sin(2*th) .* (b*b-a*a)/(2*(a*b)*(a*b));
        end

        function out = ddbdr(obj,th)
            a = obj.radiusA;
            b = obj.radiusB;
            out = 3*dbdr(th).^2./bdr(th) + bdr(th).^3 .* cos(2*th) .* (b*b-a*a)./(a*b)^2;
        end
    end    
    
end