classdef RiemannSurface
    %RIEMANNSURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        metric
        
        innerNorm
    end
    
    methods
        
        %%%-- !! TODO: yo make this cleaner !!
        
        function obj = RiemannSurface(argA, argB)
            if (nargin == 0)        % Initialize generic surface
                obj = RiemannSurface(EuclidMetric(), Domain());
            elseif (nargin == 1)    % Determine args and initialize partially generic surface
                if (isa(argA,'Metric'))
                    obj = RiemannSurface(argA, Domain());
                elseif (isa(argA,'Domain'))  
                    obj = RiemannSurface(argA, EuclidMetric());
                end
            elseif (nargin == 2)    % Initialize specified surface
                
                if (isa(argA,'Metric') && isa(argB,'Domain'))
                    obj.metric = argA;
                    obj.domain = argB;
                    return;
                elseif (isa(argA,'Domain') && isa(argB,'Metric'))  
                    obj.domain = argA;
                    obj.metric = argB;
                    return;
                end
                error('Constructor arguments must be a Metric and Domain.');
                
            end    
        end
        
        
        
        function plot(obj) % make nicer?
            holdBool = ishold;
            hold on;
            
            dom = obj.domain;
            r = dom.rMax;
            xo = dom.originX;
            yo = dom.originY;
            
            met = obj.metric;
            pr = r * 1.1;
            
            mp = met.plot((-pr:0.2:pr) + xo, (-pr:0.2:pr) + yo);
            dp = dom.plot();
            
            %dp.Color = [0.9,0.9,0.9];
            
            if (~holdBool), hold off; end;
        end
        
        
    end
end


