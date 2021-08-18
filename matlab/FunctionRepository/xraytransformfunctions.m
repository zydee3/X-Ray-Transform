% Script stores/ generates some function handles to use

% variables:
r = 0.5;
R = 1;



% Case: Euclidean Metric, Circle Domain /w radius R, integrand func
% Characteristic function of a disk
func0 = @(x,y) x.^2+y.^2 < r^2;
xray0 = @(beta,alpha) 2.*sqrt( r.*r - R.*R.*sin(...
                               max( min(abs(alpha),asin(r/R)) ,0)...
                                             ).^2 );
                                         
% Half sphere
func0 = @(x,y) sqrt(r.^2 - min(x.^2 + y.^2,r.^2));
xray0 = @(beta,alpha) pi*0.5*( r.*r - R.*R.*sin(...
                       max( min(abs(alpha),asin(r/R)) ,0)...
                                     ).^2 );                                       
                                                                                             