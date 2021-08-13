%% initialization

    % begin by initializing a surface (see EXAMPLECODE_initializeSurface)
        dc = circleDomain(radius = 1.3);
        ms = sphereMetric(radius = 2);
        rsurf = RiemannSurface(dc,ms, stepSize = 0.05, stepType = 'EE');


    % initialize some functions
        func0 = @(x,y) -(((x-0.2).^2 + abs(y-0.1)) < 0.5);
        %{
            [x,y] = ndgrid(-5:0.8:5);
            z = sin(x.^2 + y.^2) ./ (x.^2 + y.^2) .* ((x.^2 + y.^2) < 0.7);
        func0 = griddedInterpolant(x,y,z);
        %}

    % plot function under domain
        figure, hold on, axis equal
        [minB,maxB] = dc.getBoundingBox();

        [VX,VY] = ndgrid(dc.aabbgrid(250));
        pl = pcolor(VX,VY,func0(VX,VY));
        pl.EdgeColor = 'none';

        %rsurf.plotGeoFan(0);
        dc.plotAlNormal;




%% I0

    % initialize some betas an alphas to use
        beta = linspace(0,2*pi,150);
        alpha = linspace(-pi,pi,150)*0.5;
        [betaG,alphaG] = ndgrid(beta,alpha); % alternatively use meshgrid (doesnt work with griddedinterpolant or scatteredinterpolant)

    % run I0
        tic
        i0Data = rsurf.I0(betaG,alphaG, func0);
        toc

    % plot XrayI0
        figure;
        pl = pcolor(betaG,alphaG,i0Data);
        pl.EdgeColor = 'none';


%% I0 backprojeciton (I0star)

    % convert the transformation to a useable function
        xrayI0 = griddedInterpolant(betaG,alphaG,i0Data);

    % initialize points to use (we could just re-use the points used to plot the original function)
        [VX,VY] = ndgrid(dc.aabbgrid(100));

    % run I0star (this could take a moment)
        tic
        funcData = rsurf.I0star(xrayI0, VX,VY, 50);
        toc

    % plot reconstructed function
        figure, hold on, axis equal
        pl = pcolor(VX,VY,funcData);
        pl.EdgeColor = 'none';