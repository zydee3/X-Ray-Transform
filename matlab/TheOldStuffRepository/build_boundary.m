function domain = build_boundary(params)
%BUILD_BOUNDARY: builds a boundary function parameterized by an angle theta
% - bdr is the function r(\theta)
% - dbdr is its derivative
% - ddbdr is its second derivative
% - (nx,ny) are the components of the outward normal vector
%
% Author: Francois Monard
% Date: 2-28-2015

switch params.type
    case 'circle'
        bdr = @(th) params.a;
        dbdr = @(th) 0;
        ddbdr = @(th) 0;
        rmax = params.a;
        
%    case 'shifted_circle'
%        bdr = @(th) 0.4*sin(th) + sqrt((.6)^2 - (.4*cos(th)).^2);
%        dbdr = @(th) 0.4*cos(th).* (1. + 0.4*sin(th)./sqrt((.6)^2 - (.4*cos(th)).^2));
%        rmax = 1.;
        
    case 'ellipse'
        a = params.a;
        b = params.b;
        th0 = params.th0; % angle of great axis
        bdr = @(th) a*b ./ sqrt( (b*cos(th-th0)).^2 + (a*sin(th-th0)).^2 );
        dbdr = @(th) bdr(th).^3 .* sin(2*(th-th0)) .* (b^2-a^2)./2/(a*b)^2;
        ddbdr = @(th) 3*dbdr(th).^2./bdr(th) + bdr(th).^3 .* cos(2*(th-th0)) .* (b^2-a^2)./(a*b)^2;
        rmax = max(params.a, params.b);
        
    case 'cos'
        n=4;
        
        a = params.a;
        b = params.b;
        th0 = params.th0;
        bdr = @(th) a + b*cos(n*(th-th0));
        dbdr = @(th) -n*b*sin(n*(th-th0));
        ddbdr = @(th) -n*n*b*cos(n*(th-th0));
        rmax = params.a + params.b;
        
    otherwise 
        error('wrong boundary type');
        
end

% components of the non-normalized normal
nx = @(th) bdr(th).*cos(th) + dbdr(th).*sin(th);
ny = @(th) bdr(th).*sin(th) - dbdr(th).*cos(th);


domain.bdr = bdr;
domain.dbdr = dbdr;
domain.ddbdr = ddbdr;
domain.alnorm = @(th) angle(nx(th) + 1i*ny(th));
domain.rmax = rmax;
domain.xOrig = params.xOrig;
domain.yOrig = params.yOrig;

end
