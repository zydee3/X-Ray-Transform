vals = @(x,y) 1./(x.^2+y.^2+1);
B = linspace(0,2*pi,128);
A = linspace(-pi,pi,256)*0.5;

%--------------------------------------------------------------------------

dom = cosineDomain();
met = hyperbolicMetric(radius=5);

im = InMap(vals);

surf = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=10);


tic
    I0_1 = zeros(length(B),length(A));
    
    msg = ['bpart = 1 / ',num2str(length(B))]; disp(msg);

    for bpart = 1:length(B)
        if ~mod(bpart,10)
            disp(char(repmat(8,1,length(msg)+2)));
            msg = ['bpart = ',num2str(bpart),' / ',num2str(length(B))]; disp(msg);            
        end
    
        [beta,alpha] = meshgrid(bpart,A);
        I0_1(bpart,:) = XrayI0(im,surf,beta,alpha);
    end
toc

%--------------------------------------------------------------------------

params.type = 'cos';
params.a = 2;
params.b = 0.111;
params.n = 4;
params.xOrig = 0;
params.yOrig = 0;
params.rmax = 2;
params.th0 = 0;

domain = build_boundary(params);
metric = build_metric_negative(5);

tic
    I0_2 = geoI0(domain, vals, B, A, 0.01, metric, 'RK4');
toc

%--------------------------------------------------------------------------




figure, hold on
subplot(1,2,1);
pl = pcolor(beta,alpha,I0_1);
pl.EdgeColor = 'none';

subplot(1,2,2);
pl = pcolor(beta,alpha,I0_2);
pl.EdgeColor = 'none';

figure, hold on
im.plot;
%surf.plot;
surf.plotGeoFan(0);
dom.plotAlNormal;