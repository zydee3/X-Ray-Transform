classdef SmoothPolygonDomain < Domain
    % Twice-differentiable version of PolygonDomain
    
    properties
        bevelRadius = 0
        radius = 2
        sides = 4
    end
    
    methods
        function obj = SmoothPolygonDomain(radius, sides, bevelRadius)
            arguments
                radius (1,1) {mustBeNumeric} = 2;
                sides (1,1) {mustBeInteger} = 4;
                bevelRadius (1,1) {mustBeNumeric} = 0;
            end
            
            obj.radius = radius;
            obj.sides = sides;
            obj.bevelRadius = max(0,min(bevelRadius,pi/sides));
                        
            obj.bdr = @(th) bdr(th,obj.radius,obj.sides,obj.bevelRadius);
            obj.dbdr = @(th) dbdr(th,obj.radius,obj.sides,obj.bevelRadius);
            obj.ddbdr = @(th) ddbdr(th,obj.radius,obj.sides,obj.bevelRadius);
            obj.rMax = radius;         
        end    
    end    
    
end



function outr = bdr(th,r,s,br)
    pion = pi/s;
    ct = mod(th+pion,2*pion)-pion; % c
    bt = mod(th,2*pion)-pion; % b
    
    outr = (abs(ct) > br) .* r*cos(pion)./cos(bt) %!!!!!!!!!! vectorize this less bad, (use find?)
    if (br == 0), return, end;
    
    bR = mod(br,2*pion)-pion;
    q = (ct-br).*(ct+br);
    secbr = sec(bR);
    tanbr = tan(bR);
    
    constA = (br*(tanbr*tanbr+secbr*secbr)-tanbr)/(8*br*br*br);
    constB = tanbr/(2*br);
    
    outr = outr + (abs(ct) <= br) .* r*cos(pion)*secbr .* (constA.*q.*q + constB.*q + 1);
    %!!!!!!!!!!
end

function r = dbdr(th,r,s,br) % VVVV TODO!! VVVV
    pion = pi/s;
    th = mod(th,2*pion)-pion;
    r = r*cos(pion) *sin(th)./cos(th).^2;
end

function r = ddbdr(th,r,s,br)
    pion = pi/s;
    th = mod(th,2*pion)-pion;
    r = r*cos(pion) *(sin(th).^2+1)./cos(th).^3;
end