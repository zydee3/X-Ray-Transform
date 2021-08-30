clc, close all, clear


%% initialization

    % begin by initializing a surface (see EXAMPLECODE_initializeSurface)
        dc = circleDomain(radius = 1.3);
            %dc.exitInterpType = 'slinear';
        ms = sphereMetric(radius = 2);
        rsurf = RiemannSurface(dc,ms, stepSize = 0.05, stepType = 'RK4');


    % initialize some functions
        
        func0 = @(x,y)  exp(-((x-0).^2 + (y-0.5).^2) / 0.1) - exp(-((x-0.3).^2 + (y-0.5).^2) / 0.1);
        %func0 = @(x,y) double(((x-0.2).^2 + abs(y-0.1)) < 0.5); 
        
        
        %{
            [x,y] = ndgrid(-5:0.8:5);
            z = sin(x.^2 + y.^2) ./ (x.^2 + y.^2) .* ((x.^2 + y.^2) < 0.7);
        func0 = griddedInterpolant(x,y,z);
        %}

    % plot function under domain
        %{.
        figure, hold on, axis equal
        [minB,maxB] = dc.getBoundingBox();

        [VX,VY] = ndgrid(dc.aabbspace(250));
        pl = pcolor(VX,VY,func0(VX,VY));
        pl.EdgeColor = 'none';

        rsurf.plotGeoFan(0);
        dc.plotAlNormal;
        %}




%% I0

    % initialize some betas an alphas to use
        beta = linspace(0,2*pi,200);
        alpha = linspace(-pi,pi,200)*0.5;
        [betaG,alphaG] = ndgrid(beta,alpha); % alternatively use meshgrid (doesnt work with griddedinterpolant or scatteredinterpolant)

    % run I0
        tic
        i0Data = rsurf.I0(betaG,alphaG, func0);
        toc
     
    % convert the transformation to a useable function
        I0f = griddedInterpolant(betaG,alphaG,i0Data);

    % plot XrayI0
        figure;
        %pl = pcolor(betaG,alphaG,i0Data);
        pl = pcolor(betaG,alphaG,i0Data);
        pl.EdgeColor = 'none';
        title('I0')

%{.
%% I0 Inversion (I0perp*, geoA*, Hilbert, geoA of I0)
        
    % initialize points to transform over and precompute scattering relation
        beta = linspace(0,2*pi,200);
        alpha = linspace(-pi,pi,200+2)*0.5;   alpha = alpha(2:end-1);
        [VBeta,VAlpha] = ndgrid(beta, alpha);    
        [VBetaS,VAlphaS] = rsurf.scatteringRelation(VBeta,VAlpha);

    % geoA_precomp
        disp('A_-');
        [aData,betaO,alphaO] = rsurf.geoA_precomp(VBeta,VAlpha,I0f, VBetaS,VAlphaS, -1);
        
        % sorting step so that data is still in ndgrid format
        [alphaO, aind] = sort(alphaO(1,:));
        alphaO = repmat(alphaO,[height(betaO),1]);
        aData = aData(:,aind);
        
        AI0f = griddedInterpolant(betaO,alphaO,aData);
            %{.
            figure, hold on
            %pl = pcolor(VX,VY,funcData);
            pl = pcolor(betaO,alphaO,aData);
            pl.EdgeColor = 'none';
            title('AI0f');
            %}.
        
    % geoHilbert
        disp('H');    
        hData = rsurf.geoHilbert(betaO,alphaO,AI0f, 500);
        HAI0f = griddedInterpolant(betaO,alphaO,hData);
            
            %{
            figure, hold on
            %pl = pcolor(VX,VY,funcData);
            pl = pcolor(betaO,alphaO,hData);
            pl.EdgeColor = 'none';
            title('HAI0f');
            %}

    % geoAstar_precomp
        disp('A_+^*');   
        [astarData,betaO,alphaO] = rsurf.geoAstar_precomp(VBeta,VAlpha,HAI0f, VBetaS,VAlphaS, 1);
        
        % sorting step so that data is in ndgrid format
        [alphaO, aind] = sort(alphaO(1,:));
        alphaO = repmat(alphaO,[height(betaO),1]);
        astarData = astarData(:,aind);
        
        AstarHAI0f = griddedInterpolant(betaO,alphaO,astarData);
            
            %{.
            figure, hold on
            %pl = pcolor(VX,VY,funcData);
            pl = pcolor(betaO,alphaO,astarData);
            pl.EdgeColor = 'none';
            title('A*HAI0f');
            %}.
      
    % I0perpstar
        % initialize points to reconstruct at (we could just re-use the points used to plot the original function)
        [VX,VY] = ndgrid(dc.aabbspace(100));
        
        disp('Backproject I0_perp^*');   
        rsurf.stepType = 'RK4';
        tic
            %funcData = rsurf.I0perpstar(AstarHAI0f, VX,VY, 40)/(8*pi); % !!!!TODO!!!! why is this division by 8 as opposed to division by 2?
        toc


     %plot reconstructed function
        %{
        figure, hold on, axis equal
        pl = pcolor(VX,VY,funcData);
        pl.EdgeColor = 'none';
        title('reconstruction');
        %}.
     %plot error
        %{.
        figure, hold on, axis equal
        pl = pcolor(VX,VY,abs(funcData-func0(VX,VY)));
        pl.EdgeColor = 'none';
        title('error');
        %}.
%}     
        
%% I0 Inversion (I0perp*, geoR of I0)      
        beta = linspace(0,2*pi,200);
        alpha = linspace(-pi,pi,200+2)*0.5;   alpha = alpha(2:end-1);
        [VBeta,VAlpha] = ndgrid(beta, alpha);    
        [VBetaS,VAlphaS] = rsurf.scatteringRelation(VBeta,VAlpha);

        disp('R');
        [aData,betaO,alphaO] = rsurf.geoR_precomp(VBeta,VAlpha,I0f, VBetaS,VAlphaS, 100);
        
        AI0f = griddedInterpolant(betaO,alphaO,aData);


                figure, hold on, axis equal
        pl = pcolor(betaO,alphaO,aData);
        pl.EdgeColor = 'none';
        title('R');