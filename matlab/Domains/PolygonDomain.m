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
            %obj.rMax = obj.radius;        
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
        
        
        function [minB,maxB] = getBoundingBox(obj) 
            th0 = -obj.theta;
            r = obj.radius;
            pion = pi/obj.sides; 
            maxB = r*[cos(mod(th0+pion, 2*pion) - pion),...
                    cos(mod(th0+pi/2-pion, 2*pion) - pion)];
            minB = r*[cos(mod(th0+pion+pi, 2*pion) - pion - pi),...
                    cos(mod(th0-pion-pi/2, 2*pion) - pion-pi)];
        end
        
        
        function minR = getMinRadius(obj) 
            minR = cos(pi/obj.sides) * obj.radius;
        end
    end    
    
end


