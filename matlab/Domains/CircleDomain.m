classdef circleDomain < Domain
        
    properties
        radius
    end
    
    methods
        
        function obj = circleDomain(args)  
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
            obj.rMax = obj.radius;                     
        end   
        
        function out = bdr(obj,th)
            out = obj.radius * ones(size(th));
        end

        function out = dbdr(obj,th)
            out = zeros(size(th));
        end

        function out = ddbdr(obj,th)
            out = zeros(size(th));
        end
         
        
    end    
    
end