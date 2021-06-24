classdef polygonDomain < Domain
    
    properties
        radius = 2 
        sides = 4
    end
    
    methods
        function obj = polygonDomain(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2;
                args.sides (1,1) {mustBeInteger} = 2;
            end
            
            obj.radius = args.radius;
            obj.sides = args.sides;
            obj.rMax = obj.radius;        
        end 
        
        function out = bdr(obj,th)
            pion = pi/obj.sides;
            th = mod(th,2*pion)-pion;
            out = obj.radius*cos(pion) ./cos(th);
        end

        function out = dbdr(obj,th)
            pion = pi/obj.sides;
            th = mod(th,2*pion)-pion;
            out = obj.radius * cos(pion) * sin(th)./cos(th).^2;
        end

        function out = ddbdr(obj,th)
            pion = pi/obj.sides;
            th = mod(th,2*pion)-pion;
            out = obj.radius * cos(pion) * (sin(th).^2+1)./cos(th).^3;
        end
    end    
    
end


