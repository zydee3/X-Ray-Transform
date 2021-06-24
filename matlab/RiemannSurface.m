classdef RiemannSurface
    %RIEMANNSURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        metric
        
        innerNorm
    end
    
    methods
        
        
        function obj = RiemannSurface(domain, metric)
            arguments
                domain (1,1) {Domain.mustBeDomain} = circleDomain()
                metric (1,1) {Metric.mustBeMetric} = euclidMetric()
            end
            
            obj.domain = domain;
            obj.metric = metric;
        end
        
        
        
        function plot(obj) % make nicer?
            holdBool = ishold;
            hold on;
            
            dom = obj.domain;
            r = dom.rMax;
            x0 = dom.originX;
            y0 = dom.originY;
            
            met = obj.metric;
            pr = r * 1.1;
            
            mp = met.plot((-pr:0.2:pr) + x0, (-pr:0.2:pr) + y0);
            dp = dom.plot();
            
            %dp.Color = [0.9,0.9,0.9];
            
            if (~holdBool), hold off; end;
        end
        
        
    end
end


