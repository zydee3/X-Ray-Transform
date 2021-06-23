classdef CosineDomain < Domain
    
    properties
        radius = 2
        amplitude = 0.111
        cycles = 4
    end
    
    methods
        function obj = CosineDomain(radius, amplitude, cycles)
            if (nargin < 3), cycles = obj.cycles; end % !! TODO: Better way of doing this? !!
            if (nargin < 2), amplitude = obj.amplitude; end
            if (nargin == 0), radius = obj.radius; end
            
            obj.radius = radius;
            obj.amplitude = amplitude;
            obj.cycles = floor(cycles);
            n = obj.cycles;
                        
            obj.bdr = @(th) radius + amplitude*cos(n*th);
            obj.dbdr = @(th) -n*amplitude*sin(n*th);
            obj.ddbdr = @(th) -n*n*amplitude*cos(n*th);
            obj.rMax = radius + amplitude;         
        end    
    end    
    
end