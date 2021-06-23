classdef CosineDomain < Domain
    
    properties
        radius
        amplitude
        cycles
    end
    
    methods
        function obj = CosineDomain(radius, amplitude, cycles)
            arguments
                radius (1,1) {mustBeNumeric} = 2;
                amplitude (1,1) {mustBeNumeric} = 0.111;
                cycles (1,1) {mustBeInteger} = 4;
            end
            
            obj.radius = radius;
            obj.amplitude = amplitude;
            obj.cycles = cycles;
                        
            obj.bdr = @(th) radius + amplitude*cos(cycles*th);
            obj.dbdr = @(th) -cycles*amplitude*sin(cycles*th);
            obj.ddbdr = @(th) -cycles*cycles*amplitude*cos(cycles*th);
            obj.rMax = radius + amplitude;         
        end    
    end    
    
end