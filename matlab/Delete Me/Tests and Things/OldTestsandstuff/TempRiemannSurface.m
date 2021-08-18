classdef TempRiemannSurface
    % class intended to compare effects of proposed changes
    
    properties
        domain (1,1) {Domain.mustBeDomain} = circleDomain()
        metric (1,1) {Metric.mustBeMetric} = euclidMetric()
        
        stepType (1,:) {mustBeText} = 'RK4'
        stepSize (1,1) {mustBeNumeric} = 0.01
        geoDur (1,1) {mustBeNumeric} = 100
    end
    
    methods      
        function obj = TempRiemannSurface(domain, metric, args)
            %RIEMANNSURFACE Constructs an instance of RiemannSurface
            arguments
                domain (1,1) {Domain.mustBeDomain} = circleDomain()
                metric (1,1) {Metric.mustBeMetric} = euclidMetric()
                args.stepType (1,:) {mustBeText} = 'RK4'
                args.stepSize (1,1) {mustBeNumeric} = 0.01
                args.geoDur (1,1) {mustBeNumeric} = 100
            end
            
            obj.stepType = args.stepType;
            obj.stepSize = args.stepSize;
            obj.geoDur = args.geoDur;
            
            obj.domain = domain;
            obj.metric = metric;
        end
        
%--------------------------------------------------------------------------
%%                           Misc Helpers                                  
%--------------------------------------------------------------------------        

        function [xO,yO,thO] = BAtoXYTh(obj, Beta,Alpha)
            % BATOXYTH A method to convert from a (Beta,Alpha)
            % characterization of geodesics to a (X,Y,Th) one.
            dom = obj.domain;
            ra = dom.bdr(Beta - dom.theta);
            xO = cos(Beta) .* ra + dom.originX;
            yO = sin(Beta) .* ra + dom.originY;
            
            thO = pi + Alpha + dom.alNormal(Beta) + dom.theta;
        end    
                    
        
        function [betaO] = StoBeta(obj, S, stepSize)
            
            h = stepSize; % perhaps replace this with obj.stepSize
            
            sze = size(S);
            S = S(:);
            if (~issorted(S))
                S = sort(S);
                warning('input lengths are not already sorted. expect unpredictable ordering')                
            end
            if any(S<0), error('arclengths cannot be negative'), end
            
            
            betaO = zeros(sze);
            spart = 0;
            betapart = 0;
            
            for i = 1:length(S) % loop over lengths in order
                
                sD = 1;
                
                while spart < S(i)
                    sD = obj.arcLengthStep(betapart, h);
                    spart = spart + sD;
                    betapart = betapart + h;
                end
                % lin interpolate result (removes some artifacting)
                betaO(i) = betapart + h * (S(i)-spart+sD)./sD - h;
                                
                %spart = S(i);
            end    
            
            betaO = mod(betaO,2*pi);
        end    
                  
        
        function [sO] = BetatoS(obj, Beta, stepSize)
            
            h = stepSize; % perhaps replace this with obj.stepSize
            
            sze = size(Beta);
            Beta = Beta(:);
            if (~issorted(Beta))
                Beta = sort(Beta);
                warning('input lengths are not already sorted. expect unpredictable ordering')
            end
            if any(Beta<0), error('beta cannot be negative'), end
            
            
            sO = zeros(sze);
            spart = 0;
            betapart = 0;
            
            for i = 1:length(Beta) % loop over betas in order
                
                sD = 0;
                
                while betapart < Beta(i)
                    sD = obj.arcLengthStep(betapart, h);
                    betapart = betapart + h;
                    spart = spart + sD;
                end
                sO(i) = spart + sD .*(Beta(i)-betapart + h)./h - sD;
                
            end    
            
        end    
 
        
                
        
%--------------------------------------------------------------------------
%%                          Misc Computers                                 
%--------------------------------------------------------------------------        
                
        function [xO,yO,thO] = geodesic(obj,X,Y,Th)
            %GEODESIC Solves for the geodesic associated with the metric
            %given X,Y,Th
            %   Much of the geodesic function is associated with
            %   obj.geoStep.
            %   Function continues to compute until the geodesic has exited
            %   the domain or until it acheives the length described by
            %   geoDur.
            %   All geodesic functions use stepType, stepSize, and geoDur
            %   to determine the method that the geodesic travels.
            
            arguments
                obj
                X {mustBeNumeric}
                Y {mustBeNumeric}
                Th {mustBeNumeric}
            end
            
            % Reshape inputs
            tsze = size(X);
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            dom = obj.domain;
            
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

            % initialize (x,y,th)

            ngeo = length(X);

            xO = zeros(ngeo, NMAX); xO(:,1) = X;
            yO = xO;                yO(:,1) = Y;
            thO = xO;               thO(:,1) = Th;
            
            t = 1;
           
            insidepoints = ones(1,ngeo);%;dom.isInsideR2(X,Y,minR2);
            IPidx = find(insidepoints); 
            noIPidx = find(~insidepoints);
            
            while any(insidepoints) && (t ~= NMAX)

                % move the inside points forward
                [xO(IPidx,t+1), yO(IPidx,t+1), thO(IPidx,t+1)] =  ...
                    obj.geoStep(xO(IPidx,t), yO(IPidx,t), thO(IPidx,t));
          
                % keep everything fixed for points that reached the boundary
                xO(noIPidx, t+1) = xO(noIPidx, t);
                yO(noIPidx, t+1) = yO(noIPidx, t);
                thO(noIPidx, t+1) = thO(noIPidx, t);

                % march time forward
                t = t+1;    
                insidepoints(IPidx) = dom.isInsideR2(xO(IPidx, t),yO(IPidx, t),minR2);
                IPidx = find(insidepoints); 
                noIPidx = find(~insidepoints);

            end
            
            xO = xO(:, 1:t); yO = yO(:, 1:t); 
            thO = thO(:, 1:t); 
            
            xO = geoReshape(xO,tsze); % reshape to agree with inputs
            yO = geoReshape(yO,tsze);
            thO = geoReshape(thO,tsze);
        end    
   
                
        
%--------------------------------------------------------------------------
%%                               XRay                                      
%--------------------------------------------------------------------------     
        
        function [uO] = I_(obj,Beta,Alpha, integrand)
            %Integrates along an integrand that takes position and
            %direction arguments
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

            % initialize (x,y,th,int)

            ra = dom.bdr(Beta - dom.theta);
            x = cos(Beta) .* ra + dom.originX;
            y = sin(Beta) .* ra + dom.originY;

            th = pi + Alpha + dom.alNormal(Beta) + dom.theta;
            uO = zeros(size(Beta));

            t = 1;

            insidePoints = ones(size(Beta));
            IPidx = find(ones(size(Beta)));%;dom.isInsideR2(X,Y,minR2);

            while (t ~= NMAX) && (~isempty(IPidx))

                % move the inside points forward
                [x(IPidx), y(IPidx), th(IPidx), uO(IPidx)] =  ...
                    obj.geoStepI_(x(IPidx),y(IPidx),th(IPidx), uO(IPidx), integrand);

                % march time forward
                t = t+1;    
                insidePoints(IPidx) = dom.isInsideR2(x(IPidx),y(IPidx),minR2);
                IPidx = find(insidePoints); 
            end

            
            uO = uO * obj.stepSize; %reintroduce factored values
            switch obj.stepType
                case 'IE'
                    uO = uO/2;
                case 'RK4'
                    uO = uO/6;    
            end    
            
            
        end    
                
        
       
        
%--------------------------------------------------------------------------
%%                              Ploters                                    
%--------------------------------------------------------------------------        
                             
        function plotGeo(obj, X,Y,Th) 
            %PLOTGEO Plots the geodesic paths described by obj.metric using
            %obj.geodesic.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                Th {mustBeNumeric} = []
            end


            % reshape for plotting
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            [xO,yO,~] = obj.geodesic(X,Y,Th);
            
            % transpose for plotting
            xO = geoDeshape(xO)';   yO = geoDeshape(yO)';
            plot(xO,yO,'r');
        end  
        
        
        function plotGeoFan(obj, beta)
            %PLOTGEOFAN Plots a fan of geodesics at a point on the domain using
            %obj.plotGeo.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                beta (1,1) {mustBeNumeric} = 0
            end

            dom = obj.domain;
            ra = dom.bdr(beta - dom.theta);
            x = cos(beta) * ra + dom.originX;
            y = sin(beta) * ra + dom.originY;
            
            plot(x,y,'r*')
            
            holdBool = ishold;
            hold on;
            
            Th = linspace(0.5*pi, 1.5*pi, 40) + dom.alNormal(beta) + dom.theta;
            X = ones(size(Th)) * x;
            Y = ones(size(Th)) * y;
            
            obj.plotGeo(X,Y,Th);            
            
            if (~holdBool), hold off; end
        end  
        
        
        function plotGeoRadiate(obj, x,y)
            %PLOTGEORADIATE Plots a star of geodesics from a point in the domain 
            %using obj.plotGeo.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                x (1,1) {mustBeNumeric} = 0
                y (1,1) {mustBeNumeric} = 0
            end

            plot(x,y,'r*')
            holdBool = ishold;
            hold on;
            
            Th = linspace(0,2*pi, 81);
            X = ones(size(Th)) * x;
            Y = ones(size(Th)) * y;
            
            inside = find(obj.domain.isInside(X,Y));
            X = X(inside); Y = Y(inside); Th = Th(inside);
            
            obj.plotGeo(X,Y,Th);            
            
            if (~holdBool), hold off; end;
        end    
        
        
        function plotGeoParallels(obj, th)
            %PLOTGEOPARALLELS Plots an array of geodesics from points in the domain 
            %using obj.plotGeo.
            %   Geodesics are all initially parallel, th controls their direction.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                th (1,1) {mustBeNumeric} = 0
            end

            off = linspace(0,2*pi, 100);
            dom = obj.domain;
            ra = dom.bdr(off - dom.theta);
            X = cos(off) .* ra + dom.originX;
            Y = sin(off) .* ra + dom.originY;

            obj.plotGeo(X,Y,th);        
            
            holdBool = ishold;
            hold on;
            
            plot(X,Y,'r.','MarkerSize',5)
            
            if (~holdBool), hold off; end;
        end  
        
        
        function plotGeoNormals(obj, th)
            %PLOTGEONORMALS Plots an array of geodesics from points in the domain 
            %using obj.plotGeo.
            %   Geodesics initially point relative to the direction of the 
            %   inner normal of the domain, th controls their direction.
            %   This method is to be used in tandom with obj.plot.

            arguments
                obj
                th (1,1) {mustBeNumeric} = 0
            end

            plot(x,y,'r*')
            holdBool = ishold;
            hold on;
            
            off = linspace(0,2*pi, 100);
            dom = obj.domain;
            ra = dom.bdr(off - dom.theta);
            X = cos(off) * ra + dom.originX;
            Y = sin(off) * ra + dom.originY;
            Th = off - dom.theta + th;

            obj.plotGeo(X,Y,Th);    
            
            plot(X,Y,'r.','MarkerSize',5)
            
            if (~holdBool), hold off; end;
        end  

        
        function figureJacobiRadiate(obj, x,y, Th, args)
            
            arguments
                obj
                x {mustBeNumeric} = 0
                y {mustBeNumeric} = 0
                Th {mustBeNumeric} = linspace(0,2*pi,40)
                args.enableClamped (1,1) {mustBeNumericOrLogical} = 0
                args.enableAbsed (1,1) {mustBeNumericOrLogical} = 0
            end
            
            geos = prod(size(Th));

            dom = obj.domain;
            [minB,maxB] = dom.getBoundingBox;
                        
            X = ones(size(Th)) * x(1);
            Y = ones(size(Th)) * y(1);
            
            inside = find(dom.isInside(X(1:end-1),Y(1:end-1)));
            X = X(inside)';   Y = Y(inside)';   Th = Th(inside)';
            
            
            [xO,yO,~, aO,bO] = obj.geodesicJacobiAB(X,Y,Th);
            xO = geoDeshape(xO)';   yO = geoDeshape(yO)';
            aO = geoDeshape(aO)';   bO = geoDeshape(bO)';
            
            if args.enableAbsed()
                aO = abs(aO);
                bO = abs(bO);
            end

            figure;
            % left figure, plot a
            subplot(1,2,1); axis equal; hold on;
                dom.plotAlNormal;
                for i= 1:geos-1
                    s = pcolor(repmat(xO(:,i),1,2),...
                               repmat(yO(:,i),1,2),...
                               repmat(aO(:,i),1,2));
                   s.FaceAlpha=0; s.EdgeColor='interp'; s.LineWidth = 2;
                end
                title('geodesic flow a') 
                xlim([minB(1),maxB(1)]);
                ylim([minB(1),maxB(1)]);
                if args.enableClamped, caxis([-2,2]); end
            
            subplot(1,2,2); axis equal; hold on;
            % right figure, plot b
                dom.plotAlNormal;
                for i= 1:geos-1
                    s = pcolor(repmat(xO(:,i),1,2),...
                               repmat(yO(:,i),1,2),...
                               repmat(bO(:,i),1,2));
                   s.FaceAlpha=0; s.EdgeColor='interp'; s.LineWidth = 2;
                end
                title('geodesic flow b') 
                xlim([minB(1),maxB(1)]);
                ylim([minB(1),maxB(1)]);
                if args.enableClamped, caxis([-2,2]); end
        end    
              
        
        function plotConjugates(obj, X,Y,Th)
            %PLOTGEONORMALS Plots an array of geodesics from points in the domain 
            %using obj.plotGeo.
            %   Geodesics initially point relative to the direction of the 
            %   inner normal of the domain, th controls their direction.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                Th {mustBeNumeric} = []
            end

            [cX,cY] = obj.findConjugates(X,Y,Th);
            
            plot(cX,cY,'g*','MarkerSize',5);
            
        end  

        
        function plot(obj) 
            %PLOT Plots the domain on top of the metric using
            %domain.plotAlNormal and metric.plot.
            %   The region of the metric plotted is determined by domain.getBoundingBox
            
            dom = obj.domain;
            [minB,maxB] = dom.getBoundingBox();
            psizeX = 0.01 * (maxB(1) - minB(1));
            psizeY = 0.01 * (maxB(2) - minB(2));
            x0 = dom.originX;
            y0 = dom.originY;
            
            met = obj.metric;
            
            met.plot(0.1,minB,maxB);      
            
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
        
        
        
%--------------------------------------------------------------------------
%%                            Geosteppers                                  
%--------------------------------------------------------------------------        
        
        function [xO,yO,thO] = geoStep(obj,X,Y,Th)
            %GEOSTEP Steps a geodesic forward in accordance to the metric
            %and stepping methodds.
            
            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;

            switch type
                case 'EE' % Explicit Euler  ------------------------------------------------- Step
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + .5*hh.*(cth.*dylg - sth.*dxlg);
                case 'IE' % Improved Euler  ------------------------------------------------- Step
                    % predictor
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);

                    % corrector
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    elg = exp(-.5*lg);
                    
                    xO = X + hovr *    (xO  + elg.*cth);
                    yO = Y + hovr *    (yO  + elg.*sth);
                    thO = Th + hovr *  (thO + .5*elg.*(cth.*dylg - sth.*dxlg));   

                case 'RK4' % Runge-Kutta 4  ------------------------------------------------- Step
                    
                    % first slope
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    k1x = elg.*cth;
                    k1y = elg.*sth;
                    k1th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    % second slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k1x, Y+hovr*k1y);
                    cth = cos(Th+h/2*k1th); sth = sin(Th+h/2*k1th);

                    elg = exp(-.5*lg);
                    
                    k2x = elg.*cth;
                    k2y = elg.*sth;
                    k2th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    % third slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k2x, Y+hovr*k2y);
                    cth = cos(Th+h/2*k2th); sth = sin(Th+h/2*k2th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);

                    % fourth slope
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    hovr = h/6;
                    elg = exp(-.5*lg);
                    
                    xO = X + hovr * (k1x + 2*(k2x +xO) + elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO) + elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO) + 0.5*elg.*(cth.*dylg - sth.*dxlg));          
                otherwise 
                    error('wrong timestepper')
            end
        end    
                
        
        function [xO,yO,thO, uO] = geoStepI_(obj,X,Y,Th, U,integrand)
            
            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;

            switch type
                case 'EE' % Explicit Euler  ------------------------------------------------- I0
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    uO = U + integrand(X,Y);

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + 0.5 * hh.*(cth.*dylg - sth.*dxlg);
                    
                case 'IE' % Improved Euler  ------------------------------------------------- I0
                    % predictor
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);
                                        
                    % corrector
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    elg = exp(-.5*lg);
                                
                    %uO = U + (integrand(X, Y) +...
                    %          integrand(X+h*xO, Y+h*yO)); % multipy by h/2, 
                    
                    xO = X + hovr *    (xO  + elg.*cth);
                    yO = Y + hovr *    (yO  + elg.*sth);
                    thO = Th + hovr *  (thO + 0.5*elg.*(cth.*dylg - sth.*dxlg));
                    
                    %uO = U + (integrand(X, Y) + integrand(X + h*elg.*cth,Y + h*elg.*sth)); % multipy by h/2, 
                    uO = U + (integrand(X, Y) + integrand(xO,yO)); % multipy by h/2, 
                                        
                case 'RK4' % Runge-Kutta 4 ------------------------------------------------- I0

                    % first slope
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    k1x = elg.*cth;
                    k1y = elg.*sth;
                    k1th = .5*elg.*(cth.*dylg - sth.*dxlg);
                    
                    % second slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k1x, Y+hovr*k1y);
                    cth = cos(Th+h/2*k1th); sth = sin(Th+h/2*k1th);

                    elg = exp(-.5*lg);
                    
                    k2x = elg.*cth;
                    k2y = elg.*sth;
                    k2th = .5*elg.*(cth.*dylg - sth.*dxlg);
                    
                    % third slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k2x, Y+hovr*k2y);
                    cth = cos(Th+h/2*k2th); sth = sin(Th+h/2*k2th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);

                    % fourth slope
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    
                    %{
                    uO = U + (integrand(X, Y) +...
                            2*integrand(X + hovr*k1x, Y + hovr*k1y) +...
                            2*integrand(X + hovr*k2x, Y + hovr*k2y) +...
                              integrand(X + h*xO, Y + h*yO)); % multipy by h/6  
                    %}
                    hovr = h/6;                    
                    elg = exp(-.5*lg);
                    
                    xO = X + hovr * (k1x + 2*(k2x +xO) + elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO) + elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO) + 0.5*elg.*(cth.*dylg - sth.*dxlg));                  
                    %{
                    uO = U + (integrand(X, Y) +...
                            2*integrand(X + h/2*k1x, Y + h/2*k1y) +...
                            2*integrand(X + h/2*k2x, Y + h/2*k2y) +...
                              integrand(xO, yO)) * hovr; % multipy by h/6  
                    %}
                    %{.
                    uO = U + (integrand(X, Y) +...
                            4*integrand((X + xO)/2, (Y + yO)/2) +...
                              integrand(xO, yO)); % multipy by h/6  
                    %}
                otherwise 
                    error('wrong timestepper')
            end
        end    
        
        
        
  
        
    end 
end


