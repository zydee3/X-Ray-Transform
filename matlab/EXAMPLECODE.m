% To generate many generic cases, run the script "GenerateCases.m"


%% Metric
%% Methods universally available to Metrics
% obj.plot                    -- plots exp(lg/2) onto an available figure
% obj.plotALL                 -- plots lg,dxlg,dylg,curv onto a new figure (intended for debugging)

% obj.lg,dxlg,dylg,curv (x,y) -- metric function, derivatives
% obj.getHandles              -- returns function handles of lg,dxlg,dylg,curv
% obj.metricVals        (X,Y) -- computes lg,dxlg,dylg over some inputs
% obj.metricValsCurv    (X,Y) -- computes lg,dxlg,dylg,curv over some inputs

%% Static Metric Methods
% Metric.mustBeMetric(obj)    -- errors if the obj is not a metric
% Metric.build('type',params) -- constucts a sublcass with parameters



%% Domain
%% Properties universally available to Domains
% originX,originY,theta           -- transforms domain
%% Methods universally available to Domains
% obj.plot                        -- plots bdr onto an available figure
% obj.plotAABB                    -- plots axis aligned bounding box as given by getBoundingBox (intended for debugging)
% obj.plotOrigin                  -- plots a dot at originX, originY (intended for debugging)
% obj.plotAlNormal                -- plots bdr with outer normals
% obj.plotALL                     -- runs plotAABB, plotOrigin, plotAlNormal (intended for debugging)

% obj.isInside                    -- tests if a point is inside the domain
% obj.isInsideR2                  -- slightly optimized version of the above method (first tests if a point is within the minimum of bdr)

% obj.bdr,dbdr,ddbdr(th)          -- boundry function, derivatives
% obj.alNormal                    -- returns angle representing the angle of outer normal against the x axis
% obj.getHandles                  -- returns function handles of lg,dxlg,dylg,curv
% obj.getBoundingBox              -- returns vector extents the axis oriented bounding box of the domain [minB,maxB] 
% obj.getMinRadius                -- returns the minimum value of bdr
% obj.transform                   -- helper method to set values of originX,originY,theta

%% Static Domains Methods
% Domain.mustBeDomain(obj)    -- errors if the obj is not a metric
% Domain.build('type',params) -- constucts a sublcass with parameters



%% RiemannSurface
%% Properties universally available to RiemannSurfaces
% domain                        -- domain, defaults to a generic circleDomain
% metric                        -- metric, defaults to a euclidMetric

% stepType                      -- 'EE', 'IE', 'RK4'
% stepSize                      -- h
% geoDur                        -- duration that the geodesic is allowed to exist before the geodesic method assumes it is trapped
%% Methods universally available to RiemannSurfaces
% obj.plot                      -- plots domain.plotAlNormal over metric.plot
% obj.plotGeo                   -- plots geodesics directly according to the geodesic method
% obj.plotGeoFan                -- plots a geodesic fan on the boundry given a beta (currently sometimes bugged?)
% obj.plotGeoRadiate            -- plots a family of geodesics that intersect a point
% obj.plotGeoParallels          -- plots a family of parallel geodesics
% obj.plotGeoCircle             -- plots a family of geodesics that lie on a circle

% obj.geodesic                  -- computes geodesics given initial points and directions

%% Static RiemannSurface Methods
% geoStep                       -- given a x_n,y_n,theta_n, steps to a next x_{n+1},y_{n+1},theta_{n+1}






%% To construct a sphere metric:
%(Domain and metric are coded in basically the same way)

R = 2;

ms0 = sphereMetric;
ms0.radius = R; % (we can modify parameters by name after creation)

%equivically 
ms1 = sphereMetric('radius', R);

%equivically
ms2 = Metric.build('sphere','radius', R); % For this implementation, the class must be named [lowercased name]Metric
                                        % This method just calls the above constructor
%{
ms0.plotALL
ms1.plotALL
ms2.plotALL
%}
                                        
%you can also construct the same metric directly from a function handle expression of log(g), 
%but this is less efficient and can be prone to error
ms3 = Metric( @(x,y) log(4*R^4) - 2*log(x.*x + y.*y + R*R) );
%ms3.plotALL

%% The order and number of parameters referenced doesn't matter 
dcos0 = cosineDomain; % a domain constructed with the default parameters
dcos1 = cosineDomain('amplitude', 3, 'radius', 7); % doesnt reference 'cycles'
dcos2 = Domain.build('cosine','radius', 7, 'amplitude', 3, 'origin', [100,200]); % same domain as above, but shifted somewhat

%{
figure, dcos0.plot
figure, dcos1.plot
figure, dcos2.plot
%}





%% To plot geodesics
% First plot the surface and set hold to on
surf = RiemannSurface(dcos1, ms1);
surf.plot; hold on
% Then, use plot function or something
surf.plotGeoParallels(2,2,4);
%surf.plotGeoFan(1);
