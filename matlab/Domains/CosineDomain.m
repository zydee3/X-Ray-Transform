classdef cosineDomain < Domain
    
    properties
        radius
        amplitude
        cycles
    end
    
    methods
        function obj = cosineDomain(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2;
                args.amplitude (1,1) {mustBeNumeric} = 0.111;
                args.cycles (1,1) {mustBeInteger} = 4;
            end
            
            obj.radius = args.radius;
            obj.amplitude = args.amplitude;
            obj.cycles = args.cycles;
            obj.rMax = obj.radius + obj.amplitude; 
        end    
        
        function out = bdr(obj,th)
            out = obj.radius + obj.amplitude * cos(obj.cycles*th);
        end

        function out = dbdr(obj,th)
            c = obj.cycles
            out = -c * obj.amplitude * sin(c*th);
        end

        function out = ddbdr(obj,th)
            c = obj.cycles
            out = -c*c * obj.amplitude * cos(c*th);
        end
    end    
    
end