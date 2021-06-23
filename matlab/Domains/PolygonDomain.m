classdef PolygonDomain < Domain
    
    properties
        radius = 2
        sides = 4
    end
    
    methods
        function obj = PolygonDomain(radius, sides)
            arguments
                radius (1,1) {mustBeNumeric} = 2;
                sides (1,1) {mustBeNumeric} = 2;
            end
            
            obj.radius = radius;
            obj.sides = floor(sides);
                        
            obj.bdr = @(th) bdr(th,obj.radius,obj.sides);
            obj.dbdr = @(th) dbdr(th,obj.radius,obj.sides);
            obj.ddbdr = @(th) ddbdr(th,obj.radius,obj.sides);
            obj.rMax = radius;         
        end    
    end    
    
end


function r = bdr(th,r,s)
    pion = pi/s;
    th = mod(th,2*pion)-pion;
    r = r*cos(pion) ./cos(th);
end

function r = dbdr(th,r,s)
    pion = pi/s;
    th = mod(th,2*pion)-pion;
    r = r*cos(pion) *sin(th)./cos(th).^2;
end

function r = ddbdr(th,r,s)
    pion = pi/s;
    th = mod(th,2*pion)-pion;
    r = r*cos(pion) *(sin(th).^2+1)./cos(th).^3;
end