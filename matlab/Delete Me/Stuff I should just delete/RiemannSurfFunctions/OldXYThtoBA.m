        function [betaO,alphaO] = XYThtoBA(obj, X,Y,Th)
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