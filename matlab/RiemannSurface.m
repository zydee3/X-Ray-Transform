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
        
        
        function [xO,yO,thO] = geodesic(obj,xI,yI,thI)
            met = obj.metric;
            dom = obj.domain;
            
            h = 0.01;
            NMAX = floor(100/h);

            % initialize (x,y,th)

            ngeo = length(xI);

            xO = zeros(ngeo, NMAX); xO(:,1) = xI;
            yO = zeros(ngeo, NMAX); yO(:,1) = yI;
            thO = zeros(ngeo,NMAX); thO(:,1) = thI;
            
            t = 1;
            insidepoints = dom.isInside(xI,yI);

            while any(insidepoints) && (t < NMAX)
                IPidx = find(insidepoints);
                noIPidx = find(~insidepoints);

                % move the inside points forward
                [xO(IPidx,t+1), yO(IPidx,t+1), thO(IPidx,t+1)] =  ...
                    obj.geoStep(xO(IPidx,t), yO(IPidx,t), thO(IPidx,t));

                % keep everything fixed for points that reached the boundary
                xO(noIPidx, t+1) = xO(noIPidx, t);
                yO(noIPidx, t+1) = yO(noIPidx, t);
                thO(noIPidx, t+1) = thO(noIPidx, t);

                % update inside points and march time forward
                insidepoints(IPidx) = dom.isInside(xO(IPidx, t+1),yO(IPidx, t+1));
                t = t+1;    

                % debug visualisation
            %     clf    
            %     plot(exp(1i*[0:2*pi/100:2*pi, 0]), 'b'); hold on; axis equal
            %     plot(x(:,1:t)',y(:,1:t)','r'); 
            %     pause

            end
            
            xO = xO(:, 1:t-1); yO = yO(:, 1:t-1); 
            thO = thO(:, 1:t-1); 
            
            xO = xO';
            yO = yO';
            thO = thO';
        end    
        
        function plotGeoRadiate(obj, x,y)
            holdBool = ishold;
            hold on;
            obj.plot()
            
            
            thI = linspace(0,2*pi, 80);
            xI = ones(size(thI)) * x;
            yI = ones(size(thI)) * y;
            
            [xO,yO,~] = obj.geodesic(xI,yI,thI);
            plot(xO,yO,'r')
            plot(x,y,'r*')
            
            if (~holdBool), hold off; end;
        end    
        
        
        function plotGeoParallels(obj, x,y,th)
            holdBool = ishold;
            hold on;
            obj.plot()
            
            
            off = linspace(-4,4, 80);
            xI = -sin(th) * off + x;
            yI = cos(th) * off + y;
            thI = ones(size(off)) * th;
            
            [xO,yO,~] = obj.geodesic(xI,yI,thI);
            plot(xO,yO,'r')
            plot(x,y,'r*')
            plot(xI,yI,'r*','MarkerSize',1)
            
            if (~holdBool), hold off; end;
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
    
    methods (Access = 'protected')
        function [xO,yO,thO] = geoStep(obj,xI,yI,thI)
            h = 0.01;
            
            met = obj.metric;
        	[lg,dxlg,dylg] = met.metricVals(xI, yI);
            cth = cos(thI); sth = sin(thI);

            hh = exp(-0.5*lg)*h;

            xO = xI + hh.*cth;
            yO = yI + hh.*sth;
            thO = thI + .5*hh.*(cth.*dylg - sth.*dxlg);
        end    
    end    
end


