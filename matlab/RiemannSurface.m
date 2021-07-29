classdef RiemannSurface
    %RIEMANNSURFACE Contains domain and metric information, as well as
    %various ODE solvers.
    
    properties
        domain
        metric
        
        stepType
        stepSize
        geoDur
    end
    
    
    methods      
        function obj = RiemannSurface(domain, metric, args)
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
%%                               Misc
%--------------------------------------------------------------------------        

        function [xO,yO,thO] = BAtoXYTh(obj, Beta,Alpha)
            % BATOXYTH A method to convert from a (Beta,Alpha)
            % characterization of geodesics to a (X,Y,Th) one.
            dom = obj.domain;
            ra = dom.bdr(Beta - dom.theta);
            xO = cos(Beta) * ra + dom.originX;
            yO = sin(Beta) * ra + dom.originY;
            
            thO = pi + Alpha + dom.alNormal(Beta) + dom.theta;
        end    
        
        %{.
        function [betaO,alphaO] = XYThtoBA(obj, X,Y,Th)
            % XYTHTOBA A very inefficent method to convert from a (X,Y,Th)
            % characterization of geodesics to a (Beta,Alpha) one.
            %   It is expected that the given (X,Y) is inside the domain.
            %   Trapped geodesics are returned as NaN.
            %
            %   Hey, never use this.
            dom = obj.domain;
            origX = dom.originX;
            origY = dom.originX;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            % begin by computing beta (project geodesics)

            NMAX = floor(obj.geoDur/obj.stepSize);
            
            % initialize betaO,xn,yn,thn
            betaO = NaN(size(X));
            xn = zeros(size(X));   yn = zeros(size(X));   thn = zeros(size(X));
            
            t = 1;
            IPidx = find(ones(size(X))); %TODO: consider changing this?  %;dom.isInsideR2(X,Y,minR2);

            while (t ~= NMAX) && any(IPidx)
                
                % move the inside points forward (assume the point is inside at the begining of every loop)
                [xn(IPidx),yn(IPidx),thn(IPidx)] = obj.geoStep(X(IPidx),Y(IPidx),Th(IPidx));
      
                % march time forward, update X,Y,Th if they are still inside
                t = t+1;    
                IPidx = find(dom.isInsideR2(xn,yn,minR2)); % TODO dont do this
                                
                
                X(IPidx) = xn(IPidx);   Y(IPidx) = yn(IPidx);   Th(IPidx) = thn(IPidx);
            end

            % search for intersections + interpolate results (TODO: optimize)
            noIPidx = find(~dom.isInsideR2(xn,yn,minR2));   
            X = X(noIPidx);   Y = Y(noIPidx); 
            xn = xn(noIPidx);   yn = yn(noIPidx); 
            
            arcTin = atan2(Y-origY, X-origX);
            arcTout = atan2(yn-origY, xn-origX);
            
            valout = sqrt((X-origX).^2 + (Y-origY).^2) - dom.bdr(arcTout);
            valin  = sqrt((xn-origX).^2 + (yn-origY).^2) - dom.bdr(arcTin);
            
            betaO(noIPidx) = (arcTin - arcTout) .* valout./(valout-valin) + arcTout;
            
            % compute alpha
            alphaO(noIPidx) = pi + atan2(Y-yn,X-xn) - dom.alNormal(betaO(noIPidx));
            alphaO(noIPidx) = mod(alphaO(noIPidx) + pi, 2*pi) - pi;
        end    
        %}
        
        %{
        function [betaO] = arcLengthBeta(obj, stepSize, maxSteps, tol)
            % ARCLENGTHBETA Generates a list of betas distributed uniformly 
            % along the domain using newtons method and surface properties
            % defined by the metric.
            %   Very inefficient, it's best to pre-generate these values of
            %   beta only once and to store their computation in a variable
            %   for re-use.
            %
            % Scheme: 
            %   - Begin with Beta_0 = 0, 
            %   - Apply secant method with arbitrary IC Beta_0, Beta_0 + stepSize
            %       to find the zeros of the function parameterized by Beta_1 "(arclength from Beta_0 to Beta_1) - stepSize"
            %   - repeat with x until the max number of steps has been reached
            %       or until Beta_n has completed a full cycle
            
            h = obj.stepSize;   hovr = h/2;   pi2 = 2*pi;
            met = obj.metric;
            dom = obj.domain;
            
            % Initialize betas, step counter
            betaO = zeros(1,maxSteps);
            step = 1;
            
            while (step ~= maxSteps) && ~(betaO(step) >= pi2)
                %evaluate next beta
                len = @(b) () - stepSize;
                phi = @(b1,b2) 
                
                %update vals to prepare for the next step
                step = step + 1;
            end    
            
            betaO = betaO(1:step-1);
        end    
        %}
        
%--------------------------------------------------------------------------
%%                             Computers
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
            
            
            % Reshape inputs
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            dom = obj.domain;
            
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            NMAX = floor(obj.geoDur/obj.stepSize);

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
            
            sze = size(xO);
            sze = [size(X),sze(end)];
            
            xO = reshape(xO,sze); % reshape to agree with inputs
            yO = reshape(yO,sze);
            thO = reshape(thO,sze);
        end    
        
        
        function [xO,yO,thO, aO,bO] = geodesicJacobiAB(obj,X,Y,Th)
            % Reshape inputs
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            dom = obj.domain;
            
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            NMAX = floor(obj.geoDur/obj.stepSize);

            % initialize (x,y,th)

            ngeo = length(X);

            xO = zeros(ngeo, NMAX); xO(:,1) = X;
            yO = xO;                yO(:,1) = Y;
            thO = xO;               thO(:,1) = Th;
            
            % initialize ab, abdot
            ab = ones(ngeo*2,NMAX);   ab(1:end/2,1) = 0;
            abdot = flip(ab);

            t = 1;
           
            insidepoints = ones(1,ngeo);%;dom.isInsideR2(X,Y,minR2);
            IPidx = find(insidepoints);
            noIPidx = find(~insidepoints);
            
            while ~isempty(IPidx) && (t ~= NMAX)

                % move the inside points forward
                [xO(IPidx,t+1),yO(IPidx,t+1),thO(IPidx,t+1), ab([IPidx,IPidx+ngeo],t+1),abdot([IPidx,IPidx+ngeo],t+1)] =  ...
                    obj.geoStepJacobi2(xO(IPidx,t),yO(IPidx,t),thO(IPidx,t), ab([IPidx,IPidx+ngeo],t),abdot([IPidx,IPidx+ngeo],t));
          
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
                        
            % split ab into aO, bO
            bO = ab(1:end/2,:);
            aO = ab(end/2+1:end,:);
            
            aO = aO(:, 1:t); bO = bO(:, 1:t); 
            
            % reshape to agree with inputs
            sze = size(xO);
            sze = [size(X),sze(end)];
            
            xO = reshape(xO,sze);
            yO = reshape(yO,sze);
            thO = reshape(thO,sze);
            aO = reshape(aO,sze);
            bO = reshape(bO,sze);
            
        end
        
        
        function [uO] = I0(obj,Beta,Alpha, integrand)
            %I0 Solves for the Xray transfrom of the given integrand along
            %the geodesics described by Beta and Alpha.
            %   The I0 function is associated with obj.geoStepI0.
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize);

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
                    obj.geoStepI0(x(IPidx),y(IPidx),th(IPidx), uO(IPidx), integrand);

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
        
        
        function [uO] = I1(obj,Beta,Alpha, integrandU,integrandV)
            %I1 Solves for the Xray transfrom of the given integrand along
            %the geodesics described by Beta and Alpha.
            %   The I0 function is associated with obj.geoStepI0.
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize);

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
                    obj.geoStepI1(x(IPidx),y(IPidx),th(IPidx), uO(IPidx), integrandU,integrandV);

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
        
        %{.
        function [intO] = Backproject(obj, xray, Beta,Alpha, resolution)
            
            %Note: Beta,Alpha should be the values used to generate the
            %xray.
            
            dom = obj.domain;
            
            % initialize function grid
            [minB,maxB] = dom.getBoundingBox;
            intO = zeros(resolution,resolution);
            gCount = zeros(resolution,resolution);

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize);

            % initialize (x,y,th,int)
            ra = dom.bdr(Beta - dom.theta);
            x = cos(Beta) .* ra + dom.originX;
            y = sin(Beta) .* ra + dom.originY;

            th = pi + Alpha + dom.alNormal(Beta) + dom.theta;

            t = 1;

            insidePoints = ones(size(Beta));
            IPidx = find(ones(size(Beta)));%;dom.isInsideR2(X,Y,minR2);

            while (t ~= NMAX) && (~isempty(IPidx))

                % move the inside points forward
                [x(IPidx), y(IPidx), th(IPidx)] =  ...
                    obj.geoStep(x(IPidx),y(IPidx),th(IPidx));

                % place xray values onto function grid (TODO: interpolate)
                inX = clamp(1, round((x(IPidx) - minB(1)) / (maxB(1)-minB(1)) * resolution), resolution );
                inY = clamp(1, round((y(IPidx) - minB(2)) / (maxB(2)-minB(2)) * resolution), resolution );
                in = sub2ind([resolution,resolution], inY,inX);
                                
                intO(in) = intO(in) + xray(IPidx);
                gCount(in) = gCount(in) + 1;

                % march time forward
                t = t+1;    
                insidePoints(IPidx) = dom.isInsideR2(x(IPidx),y(IPidx),minR2);
                IPidx = find(insidePoints); 
            end 
            
            intO = intO./gCount;
        end   
        %}
        
        function [xO,yO] = findConjugates(obj, X,Y,Th)
            %FINDCONUGATES Identifies conjugate points characterized by given geodesics
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize);

            % initialize (b,bdot, xO,yO)
            b = zeros(size(X));   bdot = ones(size(X));
            xO = [];   yO = [];
            
            t = 1;
            insidePoints = ones(size(X));
            IPidx = find(insidePoints); %TODO: consider changing this?  %;dom.isInsideR2(X,Y,minR2);
            while (t ~= NMAX) && (~isempty(IPidx))
                
                % move the inside points forward (assume the point is inside at the begining of every loop)
                [xn, yn, thn, bn, bdotn] = obj.geoStepJacobi(X,Y,Th, b,bdot);

                % search for zeros + interpolate results (TODO: optimize)
                % (Note: broken into cases to avoid redundantly identifying
                % or failing to identify conjugate points.)
                
                %case: b_{n+1} = 0
                    zidx = find(bn == 0);
                    xO = [xO, xn(zidx)]; 
                    yO = [yO, yn(zidx)];
                %case: b_n > 0, b_{n+1} < 0
                    zidx = find(b.*bn < 0);
                    lt = bn(zidx)./(bn(zidx)-b(zidx));
                    xO = [xO, (X(zidx)-xn(zidx)).*lt+xn(zidx)];
                    yO = [yO, (Y(zidx)-yn(zidx)).*lt+yn(zidx)];

                     
                % march time forward, update X,Y,Th,b,bdot 
                t = t+1;    
                insidePoints = dom.isInsideR2(xn,yn,minR2);
                IPidx = find(insidePoints);
                
                X = xn(IPidx);   Y = yn(IPidx);   Th = thn(IPidx);
                b = bn(IPidx);   bdot = bdotn(IPidx);
            end
        end  
       
       
       
%--------------------------------------------------------------------------
%%                                Ploters
%--------------------------------------------------------------------------        
                             
        function plotGeo(obj, X,Y,Th) 
            %PLOTGEO Plots the geodesic paths described by obj.metric using
            %obj.geodesic.
            %   This method is to be used in tandom with obj.plot.
            
            % reshape for plotting
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            [xO,yO,~] = obj.geodesic(X,Y,Th);
            
            % transpose for plotting
            xO = xO(:,:)';   yO = yO(:,:)';
            plot(xO,yO,'r')
        end  
        
        
        function plotGeoFan(obj, beta)
            %PLOTGEOFAN Plots a fan of geodesics at a point on the domain using
            %obj.plotGeo.
            %   This method is to be used in tandom with obj.plot.
            
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

        
        function plotJacobiRadiate(obj, x,y)
            
            count = 40;
            
            Th = linspace(0,2*pi, count);
            X = ones(size(Th)) * x;
            Y = ones(size(Th)) * y;
            
            inside = find(obj.domain.isInside(X(1:end-1),Y(1:end-1)));
            X = X(inside)';   Y = Y(inside)';   Th = Th(inside)';
            
            
            [xO,yO,~, aO,bO] = obj.geodesicJacobiAB(X,Y,Th);
            
            for i= 1:count-1
                hold on
                s = pcolor(repmat(xO(i,:)',1,2),...
                         repmat(yO(i,:)',1,2),...
                         repmat(abs(bO(i,:))',1,2));
               s.FaceAlpha=0; s.EdgeColor='interp'; s.LineWidth = 2;
            end
        end    
              
        
        function plotConjugates(obj, X,Y,Th)
            %PLOTGEONORMALS Plots an array of geodesics from points in the domain 
            %using obj.plotGeo.
            %   Geodesics initially point relative to the direction of the 
            %   inner normal of the domain, th controls their direction.
            %   This method is to be used in tandom with obj.plot.
            
            [cX,cY] = obj.findConjugates(X,Y,Th);
            
            plot(cX,cY,'g*','MarkerSize',5);
            
        end  

        
                
        function plot(obj) % make nicer?
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
            
            met.plot(minB,maxB);      
            
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
%%                             Geosteppers
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
                
        
        function [xO,yO,thO, uO] = geoStepI0(obj,X,Y,Th, U,integrand)
            
            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;

            switch type
                case 'EE' % Explicit Euler  ------------------------------------------------- I0
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + 0.5 * hh.*(cth.*dylg - sth.*dxlg);
                    
                    uO = U + integrand(xO,yO);
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
                                
                    uO = U + (integrand(X, Y) +...
                              integrand(X+h*xO, Y+h*yO)); % multipy by h/2, 
                    
                    xO = X + hovr *    (xO  + elg.*cth);
                    yO = Y + hovr *    (yO  + elg.*sth);
                    thO = Th + hovr *  (thO + 0.5*elg.*(cth.*dylg - sth.*dxlg));
                    
                    %uO = U + (integrand(X, Y) + integrand(X + h*elg.*cth,Y + h*elg.*sth)) * hovr; % multipy by h/2, 
                    %uO = U + (integrand(X, Y) + integrand(xO,yO)) * hovr; % multipy by h/2, 
                                        
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

                    hovr = h/6;
                    elg = exp(-.5*lg);
                    %{.
                    uO = U + (integrand(X, Y) +...
                            2*integrand(X + h/2*k1x, Y + h/2*k1y) +...
                            2*integrand(X + h/2*k2x, Y + h/2*k2y) +...
                              integrand(X + h*xO, Y + h*yO)); % multipy by h/6  
                    %}
                    
                    xO = X + hovr * (k1x + 2*(k2x +xO) + elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO) + elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO) + 0.5*elg.*(cth.*dylg - sth.*dxlg));                  
                    %{
                    uO = U + (integrand(X, Y) +...
                            2*integrand(X + h/2*k1x, Y + h/2*k1y) +...
                            2*integrand(X + h/2*k2x, Y + h/2*k2y) +...
                              integrand(xO, yO)) * hovr; % multipy by h/6  
                    %{.
                    uO = U + (integrand(X, Y) +...
                            4*integrand((X + xO)/2, (Y + yO)/2) +...
                              integrand(xO, yO)) * h/6; % multipy by h/6  
                    %}
                otherwise 
                    error('wrong timestepper')
            end
        end    
        
        
        function [xO,yO,thO, uO] = geoStepI1(obj,X,Y,Th, U,integrandU,integrandV)

            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;

            switch type
                case 'EE' % Explicit Euler  ------------------------------------------------- I1
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + 0.5 * hh.*(cth.*dylg - sth.*dxlg);
                    
                    uO = U + (integrandU(xO,yO).*cth + integrandV(xO,yO).*sth) .* hh;
                case 'IE' % Improved Euler  ------------------------------------------------- I1
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
                    
                    uO = U + ((integrandU(X,Y).*cth + integrandV(X,Y).*sth) +...
                              (integrandU(X+h*xO, Y+h*yO).*cth + integrandV(X+h*xO, Y+h*yO).*sth) ); % multipy by h/2, 
                    
                    xO = X + hovr *    (xO  + elg.*cth);
                    yO = Y + hovr *    (yO  + elg.*sth);
                    thO = Th + hovr *  (thO + .5*elg.*(cth.*dylg - sth.*dxlg));   
                                        
                case 'RK4' % Runge-Kutta 4 ------------------------------------------------- I1

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
                   
                    uO = U + ((integrandU(X, Y).*cth + integrandV(X, Y).*sth) +... 
                            2*(integrandU(X + h/2*k1x, Y + h/2*k1y).*cth + integrandV(X + h/2*k1x, Y + h/2*k1y).*sth) +...
                            2*(integrandU(X + h/2*k2x, Y + h/2*k2y).*cth + integrandV(X + h/2*k2x, Y + h/2*k2y).*sth) +...
                              (integrandU(X + h*xO, Y + h*yO).*cth + integrandV(X + h*xO, Y + h*yO).*sth)  ); % multipy by h/6  
                    
                    xO = X + hovr * (k1x + 2*(k2x +xO) + elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO) + elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO) + 0.5*elg.*(cth.*dylg - sth.*dxlg));                  
                    
                otherwise 
                    error('wrong timestepper')
            end
        end    
 
        
        function [xO,yO,thO, bO,bdotO] = geoStepJacobi(obj,X,Y,Th, B,Bdot)
            
            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;
            
            switch type
                case 'EE' % Explicit Euler
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + 0.5 * hh.*(cth.*dylg - sth.*dxlg);
                    
                    bO = B + h * Bdot;
                    bdotO = Bdot - h .* curv .* B;

                case 'IE' % Improved Euler

                    % predictor
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);
                    k1bdot = - curv .* B; 

                    % corrector
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    xO = X + hovr * (xO + exp(-.5*lg).*cth);
                    yO = Y + hovr * (yO + exp(-.5*lg).*sth);
                    thO = Th + hovr * (thO + 0.5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg));
                    bO = B + hovr * (Bdot + Bdot + h*k1bdot);
                    bdotO = Bdot + hovr * (k1bdot + curv .* (h*Bdot - B));

                case 'RK4' % Runge-Kutta 4
                    
                    % first slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);

                    k1x = elg.*cth;
                    k1y = elg.*sth;
                    k1th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    k1bdot = - curv .* B;         

                    % second slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+hovr*k1x, Y+hovr*k1y);
                    cth = cos(Th+hovr*k1th); sth = sin(Th+hovr*k1th);

                    elg = exp(-.5*lg);

                    k2x = elg.*cth;
                    k2y = elg.*sth;
                    k2th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    k2b = Bdot + hovr * k1bdot;
                    k2bdot = - curv .* (B + hovr * Bdot);

                    % third slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+hovr*k2x, Y+hovr*k2y);
                    cth = cos(Th+hovr*k2th); sth = sin(Th+hovr*k2th);

                    elg = exp(-.5*lg);

                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = 0.5*elg.*(cth.*dylg - sth.*dxlg);
                    
                    k3b = Bdot + hovr * k2bdot;
                    bdotO = - curv .* (B + hovr * k2b);
                    
                    % fourth slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    elg = exp(-.5*lg);
                    
                    hovr = h/6;
                    
                    xO = X + hovr * (k1x + 2*(k2x + xO)  +elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO)  +elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO)  +0.5*elg.*(cth.*dylg - sth.*dxlg));
                    bO = B + h/6 * (Bdot + 2*(k2b + k3b) + Bdot+h*bdotO);
                    bdotO = Bdot + h/6 * (k1bdot + 2*(k2bdot + bdotO) + -curv.*(B+h*k3b));
        
                otherwise 
                    error('wrong timestepper')
            end
        end
        
        
        function [xO,yO,thO, bO,bdotO] = geoStepJacobi2(obj,X,Y,Th, B,Bdot)
            
            h = obj.stepSize;
            hovr = h/2;
            type = obj.stepType;
            met = obj.metric;
            
            switch type
                case 'EE' % Explicit Euler
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th); sth = sin(Th);

                    hh = exp(-0.5*lg)*h;

                    xO = X + hh.*cth;
                    yO = Y + hh.*sth;
                    thO = Th + 0.5 * hh.*(cth.*dylg - sth.*dxlg);
                    
                    bO = B + h * Bdot;
                    bdotO = Bdot - h .* curv .* B;

                case 'IE' % Improved Euler

                    % predictor
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);
                    
                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = .5*elg.*(cth.*dylg - sth.*dxlg);
                    k1bdot = - curv .* B; 

                    % corrector
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+h*xO, Y+h*yO);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    xO = X + hovr * (xO + exp(-.5*lg).*cth);
                    yO = Y + hovr * (yO + exp(-.5*lg).*sth);
                    thO = Th + hovr * (thO + 0.5*exp(-.5*lg).*(cth.*dylg - sth.*dxlg));
                    bO = B + hovr * (Bdot + Bdot + h*k1bdot);
                    bdotO = Bdot + hovr * (k1bdot + curv .* (h*Bdot - B));

                case 'RK4' % Runge-Kutta 4
                    
                    % first slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X, Y);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th); sth = sin(Th);

                    elg = exp(-.5*lg);

                    k1x = elg.*cth;
                    k1y = elg.*sth;
                    k1th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    k1bdot = - curv .* B;         

                    % second slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+hovr*k1x, Y+hovr*k1y);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th+hovr*k1th); sth = sin(Th+hovr*k1th);

                    elg = exp(-.5*lg);

                    k2x = elg.*cth;
                    k2y = elg.*sth;
                    k2th = .5*elg.*(cth.*dylg - sth.*dxlg);

                    k2b = Bdot + hovr * k1bdot;
                    k2bdot = - curv .* (B + hovr * Bdot);

                    % third slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+hovr*k2x, Y+hovr*k2y);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th+hovr*k2th); sth = sin(Th+hovr*k2th);

                    elg = exp(-.5*lg);

                    xO = elg.*cth;
                    yO = elg.*sth;
                    thO = 0.5*elg.*(cth.*dylg - sth.*dxlg);
                    
                    k3b = Bdot + hovr * k2bdot;
                    bdotO = - curv .* (B + hovr * k2b);
                    
                    % fourth slope
                    [lg,dxlg,dylg,curv] = met.metricValsCurv(X+h*xO, Y+h*yO);
                    curv = repmat(curv, 2,1);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    elg = exp(-.5*lg);
                    
                    hovr = h/6;
                    
                    xO = X + hovr * (k1x + 2*(k2x + xO)  +elg.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO)  +elg.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO)  +0.5*elg.*(cth.*dylg - sth.*dxlg));
                    bO = B + h/6 * (Bdot + 2*(k2b + k3b) + Bdot+h*bdotO);
                    bdotO = Bdot + h/6 * (k1bdot + 2*(k2bdot + bdotO) + -curv.*(B+h*k3b));
        
                otherwise 
                    error('wrong timestepper')
            end
        end   
        
                
                
        
    end
    
    
end


