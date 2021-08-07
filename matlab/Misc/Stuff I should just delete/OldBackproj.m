  function [intO] = backproject(obj, xray, Beta,Alpha, resolution)
            
            %Note: Beta,Alpha should be the values used to generate the
            %xray.
            %belonged in riemann
            
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