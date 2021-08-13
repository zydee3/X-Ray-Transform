function int = XrayI0(inmap,surface, Beta,Alpha) 
%XRAY Summary of this function goes here
% -- beta and alph to be implemented in the style of a meshgrid


    dom = surface.domain;

    minR2 = dom.getMinRadius;
    minR2 = minR2*minR2;

    NMAX = floor(surface.geoDur/surface.stepSize);

    % initialize (x,y,th,int)

    % this block can be optimized but is not a bottleneck
        ra = dom.bdr(Beta - dom.theta);
        x = cos(Beta) .* ra + dom.originX;
        y = sin(Beta) .* ra + dom.originY;
        
    th = pi + Alpha + dom.alNormal(Beta) + dom.theta;
    int = zeros(size(Beta));

    t = 1;

    insidepoints = ones(size(Beta));%;dom.isInsideR2(xI,yI,minR2);

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

    int = int' * surface.stepSize; % transpose for plotting/agree with inputs
    
end

