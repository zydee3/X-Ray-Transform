classdef splineDomain < Domain
    %SPLINEDOMAIN The domain to end all domains
    
    properties
        verticies = [];
        bevelRadii = [];
    end
    
    methods
        function obj = splineDomain(args)
            arguments
                args.verticies (:,2) {mustBeNumeric} = 2;
                args.bevelRadii (:,2) {mustBeNumeric} = 4;
            end
            
            obj.updateMinMax();
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

%{
piecewise scheme
input:
verticies [th_n,r_n], bevel radii [br]
piece.
it must be that th_{n+1} > th_n

%}