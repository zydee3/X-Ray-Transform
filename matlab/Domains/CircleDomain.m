classdef CircleDomain < Domain
        
    properties
        radius
    end
    
    methods
        
        function obj = CircleDomain(radius)  
            arguments
                radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = radius;
            
            obj.bdr = @(th) radius; % !! TODO: Obj dont change dynamically, write a set method and reconstruct from there? !!
            obj.dbdr = @(th) 0;
            obj.ddbdr = @(th) 0;
            obj.rMax = radius;            
        end   
         
        
    end    
    
end