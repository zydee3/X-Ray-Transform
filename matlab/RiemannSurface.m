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
        
        function obj = RiemannSurface(domain, metric)
            arguments
                domain (1,1) {Domain.mustBeDomain} = CircleDomain()
                metric (1,1) {Metric.mustBeMetric} = EuclidMetric()
            end
            
            obj.domain = domain;
            obj.metric = metric;
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


