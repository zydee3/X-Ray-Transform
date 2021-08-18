clc; close all;
addpath(genpath('./'))

%% Subclasses
em0 = EuclidMetric();
%em1 = EuclidMetric('error') %--Correctly errors

sm0 = SphereMetric();
sm1 = SphereMetric(1);
sm2 = SphereMetric(1.5376532);
%sm3 = SphereMetric('error') %--Correctly errors

h0 = HyperbolicMetric();
h1 = HyperbolicMetric(1);
h2 = HyperbolicMetric(1.5376532);
%h3 = HyperbolicMetric('error') %--Correctly errors

pm0 = PolyMetric();
pm1 = PolyMetric([12,12,12,12,12,12]);
pm2 = PolyMetric([0.1,0.2,0.3,0.4,0.5,0.6]);
%pm3 = PolyMetric('lool'); %--Correctly errors
%pm4 = PolyMetric([1,2]); %--Correctly errors

% TODO: Changing parameters doesn't update functions
%       - Set methods are awkward to implement for users, so maybe don't do
%       that
%       - Maybe just make parameters private set or allow this to keep
%       happening


%% User-implemented subclasses


%% Main class, fallbacks
m0 = Metric();
m1 = Metric(sm1);

mS.lg = 'lol'; % TODO: error checks for stuct -> metric
mS.dxlg = 'waf';
mS.dylg = 2131231;
mS.curv = 'awd';
m2 = Metric(mS);

%mG.thisisastruct = 1000; % TODO: more error checks for stuct -> metric
%m3 = Metric(mG);

R=4;
m4A = Metric( @(x,y) (4*R^4)./(x.*x+y.*y+R^2).^2 ); % TODO: confirm that numeric fallback is accurate 
    % -- Correctly warns about inefficeincy
m4B = SphereMetric(R);
%m4A.plotALL
%m4B.plotALL

m5A = Metric( @(x,y) (4*R^4)./(-x.*x-y.*y+R^2).^2 );
m5B = HyperbolicMetric(R);
%m5A.plotALL
%m5B.plotALL

coeffs = [1,0,1,0,0,0];
m6A = Metric( @(x,y) exp(coeffs(1)*x.^2 + coeffs(2)*x.*y + coeffs(3)*y.^2 + coeffs(4)*x + coeffs(5)*y + coeffs(6)) ); 
m6B = PolyMetric(coeffs);
m6A.plotALL
m6B.plotALL

clear mS mG R
