x = 0.3*ones(1,3);
y = linspace(-1,1,3);
th = 1*ones(1,3);

[b,a] = surf.geodesicFoot(x,y,th - pi);

figure; hold on;
surf.plot;
surf.plotGeo(x,y,th);

    dom = surf.domain;
    ra = dom.bdr(b - dom.theta);
        size(b)
        size(ra)
    x = cos(b) .* ra + dom.originX;
    y = sin(b) .* ra + dom.originY;
    a
surf.plotGeo(x,y,a+dom.alNormal(b));