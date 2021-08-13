classdef circleDomain < Domain
        
    properties
        radius (1,1) {mustBeNumeric} = 2
    end  
    
    methods
        
        function obj = circleDomain(args)  
            arguments
                args.radius (1,1) {mustBeNumeric} = 2
            end
            
            obj.radius = args.radius;
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
         
        
        function [minB,maxB] = getBoundingBox(obj) 
            r = obj.radius;
            maxB = [r,r];
            minB = -maxB;
        end  
        
        function minR = getMinRadius(obj) 
            minR = obj.radius;
        end 
        
    end    
    
end