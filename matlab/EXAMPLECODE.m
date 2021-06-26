%(Domain and metric are coded in basically the same way)
%% To construct a sphere metric:
R = 7;

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
ms3 = Metric( @(x,y) log(4*R^4) - 2*log(x.*x + y.*y + R*R) )
%ms3.plotALL

%% The order and number of parameters referenced doesn't matter 
dcos0 = cosineDomain % a domain constructed with the default parameters
dcos1 = cosineDomain('amplitude', 3, 'radius', 7); % doesnt reference 'cycles'
dcos2 = Domain.build('cosine','radius', 7, 'amplitude', 3, 'origin', [100,200]); % same domain as above, but shifted somewhat


figure, dcos0.plot
figure, dcos1.plot
figure, dcos2.plot



%% Methods universally available to Metrics
% obj.plot                    -- plots exp(lg/2) onto an available figure
% obj.plotALL                 -- plots lg,dxlg,dylg,curv onto a new figure (intended for debugging)
% obj.lg,dxlg,dylg,curv (x,y) -- is what it is
% obj.getHandles              -- returns function handles of lg,dxlg,dylg,curv
% obj.metricVals        (X,Y) -- computes lg,dxlg,dylg over some inputs
% obj.metricValsCurv    (X,Y) -- computes lg,dxlg,dylg,curv over some inputs

%% Static Metric Methods
% Metric.mustBeMetric(obj)    -- errors if the obj is not a metric
% Metric.build('type',params) -- constucts a sublcass with parameters


%% Properties universally available to Domains
% originX,originY,theta          -- transforms domain
%% Methods universally available to Domains
% obj.plot                       -- plots bdr onto an available figure
% obj.bdr,dbdr,ddbdr,alNorm (th) -- is what it is
% obj.getHandles                 -- returns function handles of lg,dxlg,dylg,curv
% obj.transform                  -- method to set values of originX,originY,theta

%% Static Metric Domains
% Domain.mustBeDomain(obj)    -- errors if the obj is not a metric
% Domain.build('type',params) -- constucts a sublcass with parameters

