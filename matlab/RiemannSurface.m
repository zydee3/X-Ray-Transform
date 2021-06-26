classdef RiemannSurface
    %RIEMANNSURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        metric
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
            [minB,maxB] = dom.getBoundingBox();
            psizeX = 0.01 * (maxB(1) - minB(1));
            psizeY = 0.01 * (maxB(2) - minB(2));
            x0 = dom.originX;
            y0 = dom.originY;
            
            met = obj.metric;
            
            met.plot((minB(1):psizeX:maxB(1)) + x0, (minB(2):psizeY:maxB(2)) + y0);
            dom.plotAlNorm();
            
            %consider removing this
                axis equal
                border = max(maxB(1)-minB(1),maxB(2)-minB(2)) * 0.05;
                xlim([minB(1),maxB(1)] + border * [-1,1] + x0); 
                ylim([minB(2),maxB(2)] + border * [-1,1] + y0)
            
            %dp.Color = [0.9,0.9,0.9];
            
            if (~holdBool), hold off; end;
        end
        
        
    end
end


