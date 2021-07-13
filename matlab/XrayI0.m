function int = XrayI0(inmap,surface, beta,alpha) 
%XRAY Summary of this function goes here
% -- beta and alph to be implemented in the style of a meshgrid


    dom = surface.domain;

    minR2 = dom.getMinRadius;
    minR2 = minR2*minR2;

    NMAX = floor(surface.geoDur/surface.stepSize);

    % initialize (x,y,th,int)

    % this block can be optimized but is not a bottleneck
        ra = dom.bdr(beta - dom.theta);
        x = cos(beta) .* ra + dom.originX;
        y = sin(beta) .* ra + dom.originY;
        
    th = pi + alpha + dom.alNormal(beta) + dom.theta;
    int = zeros(size(beta));

    t = 1;

    insidepoints = ones(size(beta));%;dom.isInsideR2(xI,yI,minR2);

    while (t ~= NMAX) && any(any(insidepoints))

        IPidx = find(insidepoints); 

        % move the inside points forward
        [x(IPidx), y(IPidx), th(IPidx)] =  ...
            surface.geoStep(x(IPidx), y(IPidx), th(IPidx));
        %calculate int
        int(IPidx) = int(IPidx) + inmap.eval(x(IPidx), y(IPidx));

        % march time forward
        t = t+1;    
        insidepoints(IPidx) = dom.isInsideR2(x(IPidx),y(IPidx),minR2);
    end

    int = int * surface.stepSize;
    
end

