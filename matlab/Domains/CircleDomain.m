classdef CircleDomain < Domain
        
    properties
        radius = 2
    end
    
    methods
        
        function obj = CircleDomain(radius)  
            if (nargin == 0), radius = obj.radius; end % !! TODO: Better way of doing this? !!
            
            obj.radius = radius;
            
            obj.bdr = @(th) radius; % !! TODO: Obj dont change dynamically, write a set method and reconstruct from there? !!
            obj.dbdr = @(th) 0;
            obj.ddbdr = @(th) 0;
            obj.rMax = radius;            
        end   
         
        
    end    
    
end