%% Naming Convention
% Any capitalized arguments (ie: (X,Y,Th) in many of the geosteppers) are to indicate arguments that
% share a dimension, and are always used in places where arguments are handled in a manner to produce 
% results in parallel.
%   - This convention was changed from (xI,yI,thI) to reflect to the naming
%     convention that can be found in MatLab's documentation.
%   - Lowercased arguments are sometimes still vectorized or require 
%     particular dimensions but are not typically intended to be used in a 
%     vectorized manner.
%   - This convention is broken in riemannsurface.geoStepJacobi2, where the dimensions of B
%     and Bdot must be twice the height of the other arguments (see riemannsurface.geodesicJacobiAB).

% In Plotters:
% new "plot" methods roughly reflects the behaviour of the native plot function
%   - respects hold on/hold off
%   - only generates a new figure when there isn't already one
%   - typically doesnt set axis properties or limits
%       - This convention is broken in riemannsurface.plot
% new "figure" methods always adhears to the behaviour of the native figure function
%   - always generates a new figure
%   - sometimes with titles and subplots
%   - wholly intended for debugging/ visualization 

% In RiemannSurface:
% Computers are named after what they return.
%   - obj.I0,1 returns the XrayI0,1 trasformations only.
%   - obj.geodesic returns a list of points along the geodesic
%   - obj.geodesicJacobiAB returns both the geodesic path (as does obj.geodesic) and the Jacobi a,b functions
%   - etc...
%   - This is similarly true with the associated geosteppers.





%% Metric -----------------------------------------------------------------
% - Defaults to euclidean metric

%-##-NonStatic Metric Methods-##-
%--Plotters--
% obj.plot                    -- plots exp(lg/2) onto an available figure
% obj.figureALL               -- plots lg,dxlg,dylg,curv onto a new figure (intended for debugging)

%--Characteristic Functions, Misc--
% obj.lg,dxlg,dylg,curv (x,y) -- metric function, derivatives, curvature of a surface described by the metric
% obj.getHandles              -- returns function handles of lg,dxlg,dylg,curv
% obj.metricVals        (X,Y) -- computes lg,dxlg, and dylg over some inputs
% obj.metricValsCurv    (X,Y) -- computes lg,dxlg,dylg, and curv over some inputs

%-##-Static Metric Methods-##-
% Metric.mustBeMetric(obj)    -- errors if the obj is not a metric



%% Domain -----------------------------------------------------------------
% - Defaults to circle domain with radius 2

%-##-Properties-##-
% originX,originY,theta           -- transforms domain

%-##-NonStatic Domain Methods-##-
%--Plotters--
% obj.plot                        -- plots bdr onto an available figure
% obj.plotAlNormal                -- plots bdr with outer normals
% obj.plotAABB                    -- plots axis aligned bounding box as given by getBoundingBox (intended for debugging)
% obj.plotOrigin                  -- plots a dot at originX, originY (intended for debugging)
% obj.plotBdrPoint                -- plots a dots along the boundary at given betas (intended for debugging)
% obj.plotIsInsideTest            -- plots a bunch of colour coded points to visually confirm that isInsideR2 is working as intended (intended for debugging)
% obj.plotALL                     -- runs plotAABB, plotOrigin, plotAlNormal, plotIsInsideTest (intended for debugging)

%--Characteristic Functions, Misc--
% obj.isInside                    -- tests if a point is inside the domain
% obj.isInsideR2                  -- slightly optimized version of the above method (first tests if a point is within the minimum of bdr)

% obj.bdr,dbdr,ddbdr(th)          -- boundry function, derivatives
% obj.alNormal                    -- returns angle representing the angle of outer normal against the x axis
% obj.getHandles                  -- returns function handles of lg,dxlg,dylg,curv
% obj.getMinRadius                -- returns the minimum value of bdr
% obj.transform                   -- helper method to set values of originX,originY,theta
% obj.getBoundingBox              -- returns vector extents the axis oriented bounding box of the domain [minB,maxB] 

% obj.aabbgrid                    -- generates inputs to be used in meshgrid or ndgrid using obj.getBoundingBox  

%-##-Static Domains Methods-##-
% Domain.mustBeDomain(obj)    -- errors if the obj is not a metric



%% RiemannSurface ---------------------------------------------------------
% - Defaults with euclidean metric and circle domain, RK4 steppers.

%-##-Properties-##-
% domain                        -- domain, defaults to a generic circleDomain
% metric                        -- metric, defaults to a euclidMetric

% stepType                      -- 'EE', 'IE', 'RK4'
% stepSize                      -- h, distance the geodesic can travel each step as according to the metric
% geoDur                        -- duration that the geodesic is allowed to exist before the geodesic method assumes it is trapped
%-##-NonStatic RiemannSurface Methods-##-
%--Plotters--
% obj.plot                      -- plots domain.plotAlNormal over metric.plot

% obj.plotGeo                   -- plots geodesics given as X,Y,Th
% obj.plotGeoFan                -- plots a geodesic fan on the boundry given a beta
% obj.plotGeoRadiate            -- plots a family of geodesics that intersect a point
% obj.plotGeoParallels          -- plots a family of parallel geodesics
% obj.plotGeoNormals            -- plots a family of geodesics that lie on the domain in a direction relative to the inward normals

% obj.figureJacobiRadiate       -- creates a new figure while plotting geodesic flow a,b functions onto geodesics
% obj.plotConjugates            -- plots conjugate points identified along given geodesics

%--Computers, Misc--
% obj.I0                        -- computes the xray transfomration of a scalar field 
% obj.I1                        -- computes the xray transfomration of a vector field 
% obj.I0star                    -- returns the original function of a xray transformation at a set of given points

% obj.geodesic                  -- computes geodesics given initial points and directions
% obj.geodesicJacobiAB          -- same as obj.geodesic, but also computes the jacobi a,b functions

% obj.findConjugates            -- searches for conjugate points along geodesics

% obj.StoBeta                   -- computes betas from arclength parameterization of domain [0,2pi] -> [0,totalLength]
% obj.BetatoS                   -- undoes obj.StoBeta 
% obj.BAtoXYTh                  -- converts betas and alphas to X,Y, and Theta according to the domain
% obj.geodesicEnd               -- projects geodescs forward onto the domain
% obj.geodesicFoot              -- same as geodesic end, but sends geodesics in the opposite direction of theta
% obj.scatteringRelation        -- same as geodesic end, but inputs are given as betas,alphas

%--Geosteppers--
% obj.geoStep                       -- given a x_n,y_n,theta_n, steps to a next x_{n+1},y_{n+1},theta_{n+1}
% obj.geoStepI0                     -- same as obj.geoStep, but also integrates a function
% obj.geoStepI1                     -- same as obj.geoStep, but also integrates a vector field
% obj.geoStepJacobi                 -- same as obj.geoStep, but solves the jacobi DE
% obj.geoStepJacobi2                -- same as obj.geoStepJacobi, but parallelized to handle computing the a and b functions at the same time (see obj.geodesicJacobiAB for implementation)

% obj.arcLengthStep                 -- computes approximate distance traveled along the domain from stepping beta some amount



%% MiscFunctions ---------------------------------------------------------
% geoDeshape             -- Converts from the typical output of a riemannsurface.geodesic___ method to a 2D array indexed as (nth geodesic, nth step)
% geoReshape             -- Undoes geoDeshape