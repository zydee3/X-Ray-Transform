function gRT = geoI0(domain, f, betas, alphas, h, metric, step)
%GEOI0 geodesic Xray Transform of a function using geoStep
% 
% inputs: - domain: contains information about the boundary
%         - f : function to be imaged, either defined on a grid or 
% a function handle of two variables
%         - betas, alphas: discretization of the ingoing boundary
%         - h : stepsize 
%         - metric : structure describing the metric
%         - step : 'EE', 'IE' or 'RK4'
%
% Author: Francois Monard
% Date: 10-19-2017

disp('geoI0.')

if isa(f,'function_handle')
    n = length(alphas);  
    accessFunc = @accessVal;  
elseif isa(f,'double')
    n = size(f,1);
    accessFunc = @interp4pt;
else
    error('Wrong data type for f.');
end

bdr = domain.bdr;
alnorm = domain.alnorm;
rmax = domain.rmax;
xOrig = domain.xOrig;
yOrig = domain.yOrig;

NMAX = 10*ceil(1/h);

% allocate memory for Radon transform
nbt = length(betas);
nal = length(alphas);
gRT = zeros(nal,nbt);

% loop over the boundary point
msg = ['k = 1 / ',num2str(nbt)]; disp(msg);

for k=1:nbt
    
    if ~mod(k,10)
        disp(char(repmat(8,1,length(msg)+2)));
        msg = ['k = ',num2str(k),' / ',num2str(nbt)]; disp(msg);
    end
        
    %-- point where to shoot geodesics from
    xb = xOrig + bdr(betas(k)) * cos(betas(k)); 
    yb = yOrig + bdr(betas(k)) * sin(betas(k));
    
    x = xb * ones(nal,1);
    y = yb * ones(nal,1);
    theta = alnorm(betas(k)) + pi + alphas';
    
    values = zeros(nal,1);
    insidepoints = ones(nal,1); %-- indices of points that haven't exitted
    idx = find(insidepoints);
    
    itcount = 0;
    
    while any(insidepoints) && itcount <= NMAX
        
        itcount = itcount + 1;
        
        %-- move points along geodesics
        [x(idx),y(idx),theta(idx)] = geoStep(x(idx),y(idx),theta(idx), metric, h, step);
        
        %-- update values
        values(idx) = values(idx) + ...
            accessFunc(f, x(idx), y(idx), n, rmax, xOrig, yOrig);
        
        %-- update what points are still inside the domain
        insidepoints = ((x-xOrig).*(x-xOrig) + (y-yOrig).*(y-yOrig) ...
            <= bdr(angle( x-xOrig + 1i*(y-yOrig) )).^2);% + h);
        
        idx = find(insidepoints);
    end
    
    gRT(:,k) = values;
end

gRT = h*gRT;
    
 
end

function values = accessVal(f,x,y,n,rmax,xOrig,yOrig)

values = f(x,y);

%-- debug visualization 
% if any(imag(values))
%     coco = find(imag(values)~=0)
%     disp(imag(values(coco)));
%     disp(x(coco))
%     disp(y(coco))
% end

end

function values = interp4pt(f,x,y,n,rmax,xOrig,yOrig)
%-- bilinear interpolation

dx = 2*rmax/n;

xrs = (x - xOrig + rmax-dx/2)/dx;
yrs = (y - yOrig + rmax-dx/2)/dx;

i0 = floor(yrs);
j0 = floor(xrs);
xbar = xrs - j0;
ybar = yrs - i0;            
i0 = max(min(i0 + 1,n-1),1);
j0 = max(min(j0 + 1,n-1),1);

values = ( f( n*(j0-1) + i0) .* (1.-xbar) .* (1.-ybar) ...
    + f( n*(j0-1) + i0+1) .* ybar .* (1.-xbar) ...
    + f( n*(j0+1-1) + i0) .* xbar .* (1.-ybar) ...
    + f( n*(j0+1-1) + i0+1) .* xbar .* ybar );

end