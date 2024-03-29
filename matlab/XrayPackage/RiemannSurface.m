classdef RiemannSurface
    %RIEMANNSURFACE Contains domain and metric information, as well as
    %various ODE solvers.
    
    properties
        domain (1,1) {Domain.mustBeDomain} = circleDomain()
        metric (1,1) {Metric.mustBeMetric} = euclidMetric()
        
        stepType (1,:) {mustBeText} = 'RK4'
        stepSize (1,1) {mustBeNumeric} = 0.01
        geoDur (1,1) {mustBeNumeric} = 100
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
                    
        
        function [betaO] = StoBeta(obj, S, betaStepSize)
            
            h = betaStepSize; % perhaps replace this with obj.stepSize
            
            sze = size(S);
            S = reshape(S,[],1);
            if (~issorted(S))
                S = sort(S);
                warning('input lengths are not already sorted. expect unpredictable ordering')                
            end
            if any(S<0), error('arclengths cannot be negative'), end
                        
            betaO = NaN(sze);
            spart = 0;
            betapart = 0;
                        
            i = 0;
            sD = 1;
            while i ~= length(S) % loop over lengths in order
                
                i = i+1;
                
                while spart < S(i)
                    sD = obj.arcLengthStep(betapart, h);
                    spart = spart + sD;
                    betapart = betapart + h;
                end
                % lin interpolate result (removes some artifacting)
                betaO(i) = betapart + h * (S(i)-spart+sD)./sD - h;
                if betapart >= 2*pi, break; end
            end    
            
            betaO = mod(betaO,2*pi);
        end                      
               
        function [sO] = BetatoS(obj, Beta, stepSize)
            
            h = stepSize; % perhaps replace this with obj.stepSize
            
            sze = size(Beta);
            Beta = reshape(Beta,[],1);
            if (~issorted(Beta))
                Beta = sort(Beta);
                warning('input lengths are not already sorted. expect unpredictable ordering')
            end
            if any(Beta<0), error('beta cannot be negative'), end
            
            
            sO = NaN(sze);
            spart = 0;
            betapart = 0;
            sD = 0;            
            
            for i = 1:length(Beta) % loop over betas in order
                
                while betapart < Beta(i)
                    sD = obj.arcLengthStep(betapart, h);
                    betapart = betapart + h;
                    spart = spart + sD;
                end
                sO(i) = spart + sD .*(Beta(i)-betapart + h)./h - sD;
                
            end    
            
        end    
         
        
        function [betaO,alphaO] = geodesicEnd(obj, X,Y,Th)
            % Formerly known as XYThtoBA, code transfered into geodesicFoot
            [betaO,alphaO] = obj.geodesicFoot(X,Y,Th + pi);
            alphaO = mod(alphaO-pi/2,2*pi)-pi/2;
        end    
                
        function [betaO,alphaO] = geodesicFoot(obj, X,Y,Th)
            sze = size(X);
            Th = Th + pi;
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            % begin by computing beta (project geodesics)

            NMAX = floor(obj.geoDur/obj.stepSize) + 1;
            
            % initialize xn,yn,thn (ideally, these become the first steps that are outside the domain)
            xn = zeros(sze);   yn = zeros(sze);   thn = zeros(sze);
            
            T = 1;
            insidePoints = zeros(sze);
            IPidx = find(~insidePoints); 
            
            while (T ~= NMAX) && any(IPidx)
                
                % move the inside points forward, assume all values are
                % inside for the first step
                [xn(IPidx), yn(IPidx), thn(IPidx)] = obj.geoStep(X(IPidx),Y(IPidx),Th(IPidx));
      
                % march time forward, update X,Y,Th only if they are still inside
                T = T+1;    
                insidePoints(IPidx) = dom.isInsideR2(xn(IPidx),yn(IPidx),minR2);
                IPidx = find(insidePoints); 
                
                X(IPidx) = xn(IPidx);   Y(IPidx) = yn(IPidx);   Th(IPidx) = thn(IPidx);
            end

            % search for the points that exited (just in case some points are trapped)
            noIPidx = find(~dom.isInsideR2(xn,yn,minR2));   
            
            %compute/interpolate scattering relation from the exited points
            origX = dom.originX;
            origY = dom.originY;
            betaO = NaN(sze); alphaO = NaN(sze);
            %[betaO(noIPidx),alphaO(noIPidx), ~] = dom.exitInterp(X(noIPidx)-origX, Y(noIPidx)-origY,...
            %                                                    xn(noIPidx)-origX,yn(noIPidx)-origY);
           
            betaO(noIPidx) = atan2(yn(noIPidx)-origY,xn(noIPidx)-origX);
            alphaO(noIPidx) = thn(noIPidx) - dom.alNormal(betaO(noIPidx)) + pi - dom.theta;
                                                            
        end            
               
        function [betaScattO,alphaScattO] = scatteringRelation(obj, Beta,Alpha)
            [X,Y,Th] = obj.BAtoXYTh(Beta,Alpha);
            [betaScattO,alphaScattO] = geodesicEnd(obj, X,Y,Th);
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
            X = reshape(X,[],1);   Y = reshape(Y,[],1);   Th = reshape(Th,[],1);
            
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
            
            xO = xO(:, 1:t)'; yO = yO(:, 1:t)'; 
            thO = thO(:, 1:t)'; 
            
        end    
        
        
        function [xO,yO,thO, aO,bO] = geodesicJacobiAB(obj,X,Y,Th)

            arguments
                obj
                X {mustBeNumeric}
                Y {mustBeNumeric}
                Th {mustBeNumeric}
            end

            % Reshape inputs
            tsze = size(X);

            X = reshape(X,[],1);   Y = reshape(Y,[],1);   Th = reshape(Th,[],1);
            
            dom = obj.domain;
            
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

            % initialize (x,y,th)

            ngeo = length(X);

            xO = zeros(ngeo, NMAX); xO(:,1) = X;
            yO = xO;                yO(:,1) = Y;
            thO = xO;               thO(:,1) = Th;
            
            % initialize ab, abdot
            ab = NaN(ngeo*2,NMAX);   ab(1:end/2,1) = 0;   ab(end/2+1:end,1) = 1;
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
            
            xO = xO(:, 1:t)'; yO = yO(:, 1:t)'; 
            thO = thO(:, 1:t)'; 
                        
            % split ab into aO, bO
            bO = ab(1:end/2,:);
            aO = ab(end/2+1:end,:);
            
            aO = aO(:, 1:t)'; bO = bO(:, 1:t)'; 
            
            % reshape to agree with inputs
            
        end
        
        function [xO,yO] = findConjugates(obj, X,Y,Th)
            %FINDCONUGATES Identifies conjugate points characterized by given geodesics

            arguments
                obj
                X {mustBeNumeric}
                Y {mustBeNumeric}
                Th {mustBeNumeric}
            end
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

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
                    xtemp = xn(zidx);   ytemp = yn(zidx);
                    xO = [xO; xtemp(:)];   yO = [yO; ytemp(:)];
                %case: b_n > 0, b_{n+1} < 0
                    zidx = find(b.*bn < 0);
                    lt = bn(zidx)./(bn(zidx)-b(zidx));
                          
                    xtemp = (X(zidx)-xn(zidx)).*lt+xn(zidx);   ytemp = (Y(zidx)-yn(zidx)).*lt+yn(zidx);
                    xO = [xO; xtemp(:)];   yO = [yO; ytemp(:)];
                     
                % march time forward, update X,Y,Th,b,bdot 
                t = t+1;    
                insidePoints = dom.isInsideR2(xn,yn,minR2);
                IPidx = find(insidePoints);
                
                X = xn(IPidx);   Y = yn(IPidx);   Th = thn(IPidx);
                b = bn(IPidx);   bdot = bdotn(IPidx);
            end
            
        end  
 
        
                
%--------------------------------------------------------------------------
%%                               XRay                                      
%--------------------------------------------------------------------------     
        
        function [uDataO] = I0(obj, Beta,Alpha, integrand)
            %I0 Solves for the Xray transfrom of the given integrand along
            %the geodesics described by Beta and Alpha.
            %   The I0 function is associated with obj.geoStepI0.
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

            % initialize (x,y,th,int)

            ra = dom.bdr(Beta - dom.theta);
            x = cos(Beta) .* ra + dom.originX;
            y = sin(Beta) .* ra + dom.originY;

            th = pi + Alpha + dom.alNormal(Beta) + dom.theta;
            uDataO = zeros(size(Beta));

            t = 1;

            insidePoints = ones(size(Beta));
            IPidx = find(ones(size(Beta)));%;dom.isInsideR2(X,Y,minR2);

            while (t ~= NMAX) && (~isempty(IPidx))

                % move the inside points forward
                [x(IPidx), y(IPidx), th(IPidx), uDataO(IPidx)] =  ...
                    obj.geoStepI0(x(IPidx),y(IPidx),th(IPidx), uDataO(IPidx), integrand);

                % march time forward
                t = t+1;    
                insidePoints(IPidx) = dom.isInsideR2(x(IPidx),y(IPidx),minR2);
                IPidx = find(insidePoints); 
            end

            
            uDataO = uDataO * obj.stepSize; %reintroduce factored values
            switch obj.stepType
                case 'IE'
                    uDataO = uDataO/2;
                case 'RK4'
                    uDataO = uDataO/6;    
            end    
           
            
        end    

        function [uDataO] = I1(obj, Beta,Alpha, integrandU,integrandV)
            %I1 Solves for the Xray transfrom of the given integrand along
            %the geodesics described by Beta and Alpha.
            %   The I0 function is associated with obj.geoStepI0.
            
            dom = obj.domain;

            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;

            NMAX = floor(obj.geoDur/obj.stepSize) + 1;

            % initialize (x,y,th,int)

            ra = dom.bdr(Beta - dom.theta);
            x = cos(Beta) .* ra + dom.originX;
            y = sin(Beta) .* ra + dom.originY;

            th = pi + Alpha + dom.alNormal(Beta) + dom.theta;
            uDataO = zeros(size(Beta));

            t = 1;

            insidePoints = ones(size(Beta));
            IPidx = find(ones(size(Beta)));%;dom.isInsideR2(X,Y,minR2);

            while (t ~= NMAX) && (~isempty(IPidx))

                % move the inside points forward
                [x(IPidx), y(IPidx), th(IPidx), uDataO(IPidx)] =  ...
                    obj.geoStepI1(x(IPidx),y(IPidx),th(IPidx), uDataO(IPidx), integrandU,integrandV);

                % march time forward
                t = t+1;    
                insidePoints(IPidx) = dom.isInsideR2(x(IPidx),y(IPidx),minR2);
                IPidx = find(insidePoints); 
            end

            
            %reintroduce factored values
            switch obj.stepType
                case 'IE'
                    uDataO = uDataO/2* obj.stepSize;
                case 'RK4'
                    uDataO = uDataO/6* obj.stepSize;    
            end    
            
            
        end   
                           
        function [fDataO] = I0star(obj, xray, X,Y, geosPer)
            % Xray must be a 2 parameter function defined on [0,2pi]x[-pi,pi].
            % 
            sze = size(X);
            fDataO = zeros(sze);
            
            dom = obj.domain;
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            IPidx = find(dom.isInsideR2(X,Y,minR2));
            os = ones(size(IPidx));
            
            ths = linspace(0,2*pi,geosPer+1);
            ths = ths(1:end-1);
            
            for (i = 1:geosPer)
                [beta,alpha] = obj.geodesicFoot(X(IPidx),Y(IPidx),ths(i)*os);
                fDataO(IPidx) = fDataO(IPidx) + xray(mod(beta,2*pi),mod(alpha+pi,2*pi)-pi);
            end    
            
            fDataO = fDataO * 2*pi/geosPer;
            
        end
        
        function [fDataO] = I0perpstar(obj, xray, X,Y, geosPer)

            sze = size(X);
            
            dom = obj.domain;
            minR2 = dom.getMinRadius;
            minR2 = minR2*minR2;
            
            met = obj.metric;
            
            IPidx = find(dom.isInsideR2(X,Y,minR2));
            os = ones(size(IPidx));
            
            ths = linspace(0,2*pi,geosPer+1);
            ths = ths(1:end-1);
            sth = sin(ths);  cth = cos(ths);
            
            u = zeros(sze);   v = zeros(sze);
                        
            for (k = 1:geosPer)
                [beta,alpha] = obj.geodesicFoot(X(IPidx),Y(IPidx),ths(k)*os);
                
                val = xray(mod(beta,2*pi),mod(alpha+pi/2,pi)-pi/2);
                
                u(IPidx) = u(IPidx) - val * sth(k);
                v(IPidx) = v(IPidx) + val * cth(k);
            end    
            
            % multiplication by g^{1/2}
            lgt = met.metricVals(X,Y);
            u = u .* exp(0.5*lgt);
            v = v .* exp(0.5*lgt);

            % divergence
            val = divergence(Y,X,v,u); % !!!!!TODO!!!!! figure out why divergence(X,Y,u,v) doesnt work

            % multiplication by g^{-1} and integral over angles
            fDataO = val .* exp(-lgt) * 2*pi/geosPer; 
            
            % crop image TODO remove this or do this differently somewhere else
            IPidx = ~dom.isInsideR2(X*1.0526,Y*1.0526,minR2);
            fDataO(IPidx) = 0;
        end
        
        
        
        function [fDataO, betaO,alphaO] = geoA_precomp(~, Beta,Alpha,func, BetaScatt,AlphaScatt, sign)
            po2 = pi/2;
            
            bscatt = mod(BetaScatt, 2*pi);
            ascatt = mod(AlphaScatt+po2,pi)-po2;
            
            betaO = [Beta, Beta];
            alphaO = [Alpha, Alpha+pi];
            fDataO = [func(Beta,Alpha), sign*func(bscatt, ascatt)];
        end      
      
        function [fDataO] = geoHilbert(~, Beta,Alpha,func, nalo2)
            % F is handled as an anonymous function
            % for each input beta, descretizes the function along alpha in order to perform FFT
            % output is interpolated samples of fft output
            
            sze = size(Beta);
            ubt = unique(reshape(Beta,[],1))';
            
            al = linspace(-pi/2,3*pi/2,nalo2*2);
            Alpha = mod(Alpha+pi/2,2*pi)-pi/2;
            
            fDataO = zeros(sze);
            
            H = 2*nalo2*[-1j*ones(1,nalo2), 1j*ones(1,nalo2)];
            
            for bt = ubt
                alind = (bt == Beta); 
                 
                fDataO(alind) = interp1(al, ifft(H .* fft( func(mod(ones(1,nalo2*2)*bt,2*pi),al), nalo2*2 ) ) /(nalo2*2),...
                                         Alpha(alind), 'linear');
            end
            fDataO = real(fDataO);

        end   
        
        function [fDataO, betaO,alphaO] = geoAstar_precomp(~, Beta,Alpha,func, BetaScatt,AlphaScatt, sign)            
            po2 = pi/2;
            
            bscatt = mod(BetaScatt, 2*pi);
            ascatt = mod(AlphaScatt+po2+pi,2*pi)-po2; % !!!!!----!!!!! figure out why mod(AlphaScatt+po2,2*pi)-po2; doesnt work
                                                   % edit: TOOD!!!! figure out what changed so that its no longer mod(-AlphaScatt+po2,2*pi)-po2; doesnt work
            fDataO = func(Beta,Alpha) + sign*func(bscatt, ascatt); 
            betaO = Beta;
            alphaO = Alpha;
        end  
        
        function [fDataO, Beta,Alpha] = geoR_precomp(~, Beta,Alpha,func, BetaScatt,AlphaScatt, nalo2)
            
            % A operator
            po2 = pi/2;
            
            bscatt = mod(BetaScatt, 2*pi);
            ascatt = mod(AlphaScatt+po2,pi)-po2;
            
            intF = scatteredInterpolant(reshape([Beta, Beta],                              [],1),...
                                        reshape([Alpha, Alpha+pi],                         [],1),...
                                        reshape([func(Beta,Alpha), -func(bscatt, ascatt)], [],1) );
             
            
            % Hilbert + A* operator
            
            ubt = unique(reshape(Beta,[],1))'; % TODO: add a threshold for very near fibres to be considered as the same fibre
            fDataO = zeros(size(Beta));
            al = linspace(-pi/2,3*pi/2,nalo2*2);
            H = 2*nalo2*[-1j*ones(1,nalo2), 1j*ones(1,nalo2)];
                        
            for bt = ubt
                alind = (bt == Beta); 
                
                HAbt = ifft(H .* fft( intF(mod(ones(1,nalo2*2)*bt,2*pi),al), nalo2*2 ) ) /(nalo2*2);
                
                fDataO(alind) = interp1(al, HAbt, Alpha(alind), 'linear') +...
                                interp1(al, HAbt, Alpha(alind), 'linear'); % TODO: this is incorrect, should be ascatt + pi on associated bscatt fibres, consider interpolating for these values
            end
           
            
            fDataO = real(fDataO)/4;
            
        end    
   
        function [fDataO, Beta,Alpha] = geoR_simple(~, Beta,Alpha,func, nalo2)
            % Covers only the case that func(Beta,Alpha+pi) = -func(Beta,Alpha)
            % A operator
            po2 = pi/2;
            
            fba = func(Beta,Alpha);
            
            intF = scatteredInterpolant(reshape([Beta, Beta],       [],1),...
                                        reshape([Alpha, Alpha+pi],  [],1),...
                                        reshape([fba, -fba],        [],1) );                    
                                    
            % Hilbert + A* operator
            
            ubt = unique(reshape(Beta,[],1))';
            fDataO = zeros(size(Beta));
            al = linspace(-pi/2,3*pi/2,nalo2*2);
            H = 2*nalo2*[-1j*ones(1,nalo2), 1j*ones(1,nalo2)];
            
            for bt = ubt
                alind = (bt == Beta); 
                
                HAbt = ifft(H .* fft( intF(mod(ones(1,nalo2*2)*bt,2*pi),al), nalo2*2 ) ) /(nalo2*2);
                
                fDataO(alind) = interp1(al, HAbt, Alpha(alind), 'linear') - ... % TODO figure out why + doesnt work
                                interp1(al, HAbt, Alpha(alind) + pi, 'linear');
            end
            
            fDataO = real(fDataO)/4;
            
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


            % reshape for plotting (not totally neccessary)
            X = reshape(X,[],1);   Y = reshape(Y,[],1);   Th = reshape(Th,[],1);

            [xO,yO,~] = obj.geodesic(X,Y,Th);

            plot(xO,yO,'r');
        end          
        
        function plotGeoFan(obj, beta, numgeos)
            %PLOTGEOFAN Plots a fan of geodesics at a point on the domain using
            %obj.plotGeo.
            %   This method is to be used in tandom with obj.plot.
            
            arguments
                obj
                beta (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 0
            end

            dom = obj.domain;
            ra = dom.bdr(beta - dom.theta);
            x = cos(beta) * ra + dom.originX;
            y = sin(beta) * ra + dom.originY;
            
            plot(x,y,'r*')
            
            holdBool = ishold;
            hold on;
            
            Th = linspace(0.5*pi, 1.5*pi, numgeos) + dom.alNormal(beta) + dom.theta;
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
        
        function plotGeoNormals(obj, th,numgeos)
            %PLOTGEONORMALS Plots an array of geodesics from points in the domain 
            %using obj.plotGeo.
            %   Geodesics initially point relative to the direction of the 
            %   inner normal of the domain, th controls their direction.
            %   This method is to be used in tandom with obj.plot.

            arguments
                obj
                th (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 0                
            end
            
            beta = linspace(0,2*pi, numgeos+1); beta = beta(1:end-1);
            dom = obj.domain;
            ra = dom.bdr(beta - dom.theta);
            X = cos(beta) .* ra + dom.originX;
            Y = sin(beta) .* ra + dom.originY;
            Th = th + dom.alNormal(beta) + pi;
            
            plot(X,Y,'r.','MarkerSize',5)
            holdBool = ishold;
            hold on;  
            
            obj.plotGeo(X,Y,Th)
                        
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
                args.enablePlotConjugates (1,1) {mustBeNumericOrLogical} = 0
            end
            
            dom = obj.domain;
            [minB,maxB] = dom.getBoundingBox;
                        
            X = ones(size(Th)) * x(1);
            Y = ones(size(Th)) * y(1);
            
            inside = find(dom.isInside(X(1:end),Y(1:end)));
            X = X(inside);   Y = Y(inside);   Th = Th(inside);
            geos = length(Th);
            
            [xO,yO,~, aP,bP] = obj.geodesicJacobiAB(X,Y,Th);
            aO = aP;   bO = bP;
            
            if args.enableAbsed()
                aP = abs(aP);
                bP = abs(bP);
            end

            figure;
            % left figure, plot a
            subplot(1,2,1); axis equal; hold on;
                dom.plotAlNormal;
                for i= 1:geos
                    s = pcolor(repmat(xO(:,i),1,2),...
                               repmat(yO(:,i),1,2),...
                               repmat(aP(:,i),1,2));
                   s.FaceAlpha=0; s.EdgeColor='interp'; s.LineWidth = 2;
                end
                title('geodesic flow a') 
                xlim([minB(1),maxB(1)]);
                ylim([minB(1),maxB(1)]);
                if args.enableClamped, caxis([-2,2]); end
            
            subplot(1,2,2); axis equal; hold on;
            % right figure, plot b
                dom.plotAlNormal;
                
                for i= 1:geos
                    s = pcolor(repmat(xO(:,i),1,2),...
                               repmat(yO(:,i),1,2),...
                               repmat(bP(:,i),1,2));
                   s.FaceAlpha=0; s.EdgeColor='interp'; s.LineWidth = 2;
                end
                
                title('geodesic flow b') 
                xlim([minB(1),maxB(1)]);
                ylim([minB(1),maxB(1)]);
                if args.enableClamped, caxis([-2,2]); end
                
                
                if args.enablePlotConjugates
                    
                    %case: b_n > 0, b_{n+1} < 0
                    zidx = find([(bO(1:end-1,:).*bO(2:end,:)) < 0; zeros(1,geos)]);
                    t = bO(zidx)./(bO(zidx)-bO(zidx+1));
                    X = (xO(zidx+1)-xO(zidx)).*t+xO(zidx);
                    Y = (yO(zidx+1)-yO(zidx)).*t+yO(zidx);
                    %X = xO(zidx); Y = yO(zidx);
                    %case: b_{n+1} = 0
                    zidx = find([bO(2:end,:)==0; zeros(1,geos)]);
                    X = [X,xO(zidx)]; Y = [Y,yO(zidx)];
                    
                    
                    plot(X,Y,'g*','MarkerSize',5);
                end
                    
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
            
            met.plot(200,minB,maxB);      
            
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

                    elg1 = exp(-.5*lg);
                    
                    xO = elg1.*cth;
                    yO = elg1.*sth;
                    thO = .5*elg1.*(cth.*dylg - sth.*dxlg);
                                                            
                    % corrector
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    elg2 = exp(-.5*lg);
                    
                    uO = U + ((integrandU(X,Y).*cth + integrandV(X,Y).*sth) .* elg1 +...
                              (integrandU(X+h*xO, Y+h*yO).*cth + integrandV(X+h*xO, Y+h*yO).*sth) .* elg2 ); % multipy by h/2, 
                    
                    xO = X + hovr *    (xO  + elg2.*cth);
                    yO = Y + hovr *    (yO  + elg2.*sth);
                    thO = Th + hovr *  (thO + .5*elg2.*(cth.*dylg - sth.*dxlg));   
                                        
                case 'RK4' % Runge-Kutta 4 ------------------------------------------------- I1

                    % first slope
                    [lg,dxlg,dylg] = met.metricVals(X, Y);
                    cth = cos(Th); sth = sin(Th);

                    elg1 = exp(-.5*lg);
                    
                    k1x = elg1.*cth;
                    k1y = elg1.*sth;
                    k1th = .5*elg1.*(cth.*dylg - sth.*dxlg);
                    
                    % second slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k1x, Y+hovr*k1y);
                    cth = cos(Th+h/2*k1th); sth = sin(Th+h/2*k1th);

                    elg2 = exp(-.5*lg);
                    
                    k2x = elg2.*cth;
                    k2y = elg2.*sth;
                    k2th = .5*elg2.*(cth.*dylg - sth.*dxlg);
                    
                    % third slope
                    [lg,dxlg,dylg] = met.metricVals(X+hovr*k2x, Y+hovr*k2y);
                    cth = cos(Th+h/2*k2th); sth = sin(Th+h/2*k2th);

                    elg3 = exp(-.5*lg);
                    
                    xO = elg3.*cth;
                    yO = elg3.*sth;
                    thO = .5*elg3.*(cth.*dylg - sth.*dxlg);

                    % fourth slope
                    [lg,dxlg,dylg] = met.metricVals(X+h*xO, Y+h*yO);
                    cth = cos(Th+h*thO); sth = sin(Th+h*thO);

                    hovr = h/6;
                    elg4 = exp(-.5*lg);
                   
                    uO = U + ((integrandU(X, Y).*cth + integrandV(X, Y).*sth) .* elg1 +... 
                            2*(integrandU(X + h/2*k1x, Y + h/2*k1y).*cth + integrandV(X + h/2*k1x, Y + h/2*k1y).*sth) .* elg2 +...
                            2*(integrandU(X + h/2*k2x, Y + h/2*k2y).*cth + integrandV(X + h/2*k2x, Y + h/2*k2y).*sth) .* elg3 +...
                              (integrandU(X + h*xO, Y + h*yO).*cth + integrandV(X + h*xO, Y + h*yO).*sth) .* elg4  ); % multipy by h/6  
                    
                    xO = X + hovr * (k1x + 2*(k2x +xO) + elg4.*cth);
                    yO = Y + hovr * (k1y + 2*(k2y + yO) + elg4.*sth);
                    thO = Th + hovr * (k1th + 2*(k2th + thO) + 0.5*elg4.*(cth.*dylg - sth.*dxlg));                  
                    
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
        
                
        function [sDO] = arcLengthStep(obj, Beta, stepSize)
            % ARCLENGTHSTEP
            % Unlike the geostep functions, this method does not update
            % values, rather evaluates the difference between steps.
            h = stepSize; % perhaps replace this with obj.stepSize
            hovr = h/2;
            met = obj.metric;
            dom = obj.domain;
            xoff = dom.originX; yoff = dom.originY;
            type = obj.stepType;
            th0 = dom.theta;
            
            switch type
                case 'EE' % Explicit Euler  ------------------------------------------------- ArcStep
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                                        
                    sDO = h * exp(0.5*met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db);
                    
                case 'IE' % Trapazoid rule  ------------------------------------------------- ArcStep
                    % left sample
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                                        
                    sDO = exp(0.5*met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db);
                    
                    % right sample
                    Beta = Beta + h;
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                    
                    sDO = hovr * (0.5*exp(met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db) + sDO);
                case 'RK4' % Simpsons Rule  ------------------------------------------------- ArcStep
                    % left sample
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                                        
                    k1sD = exp(0.5*met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db);
                    
                    % middle sample
                    Beta = Beta + hovr;
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                    
                    sDO = exp(0.5*met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db);
                    
                    % right sample                    
                    Beta = Beta + hovr;
                    cth = cos(Beta); sth = sin(Beta);
                    b = dom.bdr(Beta-th0);
                    db = dom.dbdr(Beta-th0);
                    
                    sDO = h/6 * (0.5*exp(met.lg(cth.*b+xoff, sth.*b+yoff)) .* sqrt(b.*b + db.*db) + 4*sDO + k1sD);
                    
                otherwise 
                    error('wrong timestepper')
            end
        end       
        
%--------------------------------------------------------------------------        
    end 
    
    methods (Static)
       function mustBeSurface(obj)     
           %MUSTBEMETRIC Errors if the passed object is not an instance of 
           %RiemannSurface or of a subclass of RiemannSurface.
           if (~isa(obj,'RiemannSurface'))
               error("Value must be a RiemannSurface.");
           end
       end  
       
    end
end


