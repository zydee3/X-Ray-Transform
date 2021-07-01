classdef RiemannSurface
    %RIEMANNSURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        metric
        
        stepType
        stepSize
        geoDur
    end
    
    methods
        
        
        function obj = RiemannSurface(domain, metric, args)
            arguments
                domain (1,1) {Domain.mustBeDomain} = circleDomain()
                metric (1,1) {Metric.mustBeMetric} = euclidMetric()
                args.stepType (1,:) {mustBeText} = 'EE'
                args.stepSize (1,1) {mustBeNumeric} = 0.01
                args.geoDur (1,1) {mustBeNumeric} = 100
            end
            
            obj.stepType = args.stepType;
            obj.stepSize = args.stepSize;
            obj.geoDur = args.geoDur;
            
            obj.domain = domain;
            obj.metric = metric;
        end
        
        
        
        
        
        function [xO,yO,thO] = geodesic(obj,xI,yI,thI)
            dom = obj.domain;
            
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            NMAX = floor(obj.geoDur/obj.stepSize);

            % initialize (x,y,th)

            ngeo = length(xI);

            xO = zeros(ngeo, NMAX); xO(:,1) = xI;
            yO = zeros(ngeo, NMAX); yO(:,1) = yI;
            thO = zeros(ngeo,NMAX); thO(:,1) = thI;
            
            t = 1;
           
            insidepoints = ones(1,ngeo);%;dom.isInsideR2(xI,yI,minR2);

            while any(insidepoints) && (t ~= NMAX)
                
                IPidx = find(insidepoints); 
                noIPidx = find(~insidepoints);

                % move the inside points forward
                [xO(IPidx,t+1), yO(IPidx,t+1), thO(IPidx,t+1)] =  ...
                    obj.geoStep(xO(IPidx,t), yO(IPidx,t), thO(IPidx,t));

                insidepoints(IPidx) = dom.isInsideR2(xO(IPidx, t),yO(IPidx, t),minR2);
                
                % keep everything fixed for points that reached the boundary
                xO(noIPidx, t+1) = xO(noIPidx, t);
                yO(noIPidx, t+1) = yO(noIPidx, t);
                thO(noIPidx, t+1) = thO(noIPidx, t);

                % march time forward
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
        
        
              
        
        
        function plotGeo(obj, X,Y,Th)       
            [xO,yO,~] = obj.geodesic(X,Y,Th);
            
            plot(xO,yO,'r')
        end  
        
        
        function plotGeoFan(obj, beta)
            
            dom = obj.domain;
            ra = dom.bdr(beta - dom.theta);
            x = cos(beta) * ra + dom.originX;
            y = sin(beta) * ra + dom.originY;
            
            plot(x,y,'r*')
            
            holdBool = ishold;
            hold on;
            
            thI = linspace(0.5*pi, 1.5*pi, 40) + dom.alNormal(beta) + dom.theta;
            xI = ones(size(thI)) * x;
            yI = ones(size(thI)) * y;
            
            obj.plotGeo(xI,yI,thI);            
            
            if (~holdBool), hold off; end
        end  
        
        
        function plotGeoRadiate(obj, x,y)
            plot(x,y,'r*')
            holdBool = ishold;
            hold on;
            
            thI = linspace(0,2*pi, 81);
            xI = ones(size(thI)) * x;
            yI = ones(size(thI)) * y;
            
            inside = find(obj.domain.isInside(xI,yI));
            xI = xI(inside); yI = yI(inside); thI = thI(inside);
            
            obj.plotGeo(xI,yI,thI);            
            
            if (~holdBool), hold off; end;
        end    
        
        
        function plotGeoParallels(obj, x,y,th)
            
            plot(x,y,'r*')
            holdBool = ishold;
            hold on;
            
            off = linspace(-4,4, 81);
            xI = -sin(th) * off + x;
            yI = cos(th) * off + y;
            thI = ones(size(off)) * th;
            
            inside = find(obj.domain.isInside(xI,yI));
            xI = xI(inside); yI = yI(inside); thI = thI(inside);

            obj.plotGeo(xI,yI,thI);           
            plot(xI,yI,'r.','MarkerSize',5)
            
            if (~holdBool), hold off; end;
        end  
        
        
        function plotGeoCircle(obj, x,y,r,th)
            plot(x,y,'r*')
            holdBool = ishold;
            hold on;
            
            
            off = linspace(0,2*pi, 81);
            xI = sin(off) * r + x;
            yI = cos(off) * r + y;
            thI = pi*0.5 -off + th;
            
            inside = find(obj.domain.isInside(xI,yI));
            xI = xI(inside); yI = yI(inside); thI = thI(inside);
          
            obj.plotGeo(xI,yI,thI);
            plot(xI,yI,'r.','MarkerSize',5)
            
            if (~holdBool), hold off; end;
        end  
        
        function plot(obj) % make nicer?
            
            dom = obj.domain;
            [minB,maxB] = dom.getBoundingBox();
            psizeX = 0.01 * (maxB(1) - minB(1));
            psizeY = 0.01 * (maxB(2) - minB(2));
            x0 = dom.originX;
            y0 = dom.originY;
            
            met = obj.metric;
            
            met.plot((minB(1):psizeX:maxB(1)) + x0, (minB(2):psizeY:maxB(2)) + y0);      
            
            holdBool = ishold;
            hold on;
            
            dom.plotAlNormal();
            
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
            h = obj.stepSize;
            type = obj.stepType;
            met = obj.metric;

            switch type
                case 'EE' % Explicit Euler
                    [lg,dxlg,dylg] = met.metricVals(xI, yI);
                    cth = cos(thI); sth = sin(thI);

                    hh = exp(-0.5*lg)*h;

                    xO = xI + hh.*cth;
                    yO = yI + hh.*sth;
                    thO = thI + .5*hh.*(cth.*dylg - sth.*dxlg);
                case 'IE' % Improved Euler
                    % predictor
                    [lg,dxlg,dylg] = met.metricVals(xI, yI);
                    cth = cos(thI); sth = sin(thI);

                    k1x = exp(-.5*lg).*cth;
                    k1y = exp(-.5*lg).*sth;
                    k1th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h * k1x;
                    yO = yI + h * k1y;
                    thO = thI + h * k1th;

                    % corrector
                    [lg,dxlg,dylg] = met.metricVals(xO, yO);
                    cth = cos(thO); sth = sin(thO);

                    k2x = exp(-.5*lg).*cth;
                    k2y = exp(-.5*lg).*sth;
                    k2th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h * (k1x + k2x)/2;
                    yO = yI + h * (k1y + k2y)/2;
                    thO = thI + h * (k1th + k2th)/2;

                case 'RK4' % Runge-Kutta 4

                    % first slope
                    [lg,dxlg,dylg] = met.metricVals(xI, yI);
                    cth = cos(thI); sth = sin(thI);

                    k1x = exp(-.5*lg).*cth;
                    k1y = exp(-.5*lg).*sth;
                    k1th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h/2 * k1x;
                    yO = yI + h/2 * k1y;
                    thO = thI + h/2 * k1th;

                    % second slope
                    [lg,dxlg,dylg] = met.metricVals(xO, yO);
                    cth = cos(thO); sth = sin(thO);

                    k2x = exp(-.5*lg).*cth;
                    k2y = exp(-.5*lg).*sth;
                    k2th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h/2 * k2x;
                    yO = yI + h/2 * k2y;
                    thO = thI + h/2 * k2th;

                    % third slope
                    [lg,dxlg,dylg] = met.metricVals(xO, yO);
                    cth = cos(thO); sth = sin(thO);

                    k3x = exp(-.5*lg).*cth;
                    k3y = exp(-.5*lg).*sth;
                    k3th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h * k3x;
                    yO = yI + h * k3y;
                    thO = thI + h * k3th;

                    % fourth slope
                    [lg,dxlg,dylg] = met.metricVals(xO, yO);
                    cth = cos(thO); sth = sin(thO);

                    k4x = exp(-.5*lg).*cth;
                    k4y = exp(-.5*lg).*sth;
                    k4th = .5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg);

                    xO = xI + h/6 * (k1x + 2*k2x + 2*k3x + k4x);
                    yO = yI + h/6 * (k1y + 2*k2y + 2*k3y + k4y);
                    thO = thI + h/6 * (k1th + 2*k2th + 2*k3th + k4th);        
                otherwise 
                    error('wrong timestepper')
            end
        end    
    end    
end


