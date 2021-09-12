%% clear everything and set up plotting settings
clc, close all, clear

sett = Figureizer.settings;
sett.swapSettings(FigureizerSettings.classic);
sett.gridResolution = 80;

%% initialization

    % begin by initializing a surface (see EXAMPLECODE_initializeSurface)
        dc = circleDomain(radius = 1.3);
        %dc = smoothpolyDomain(radius = 1.3, sides = 5, bevelRadius = 0.2);
            %dc.exitInterpType = 'squad';
            dc.theta = 0.2 + pi/2;
            dc.originY = 0;
            dc.originX = -1;
        ms = sphereMetric(radius = 1.7);
        rsurf = RiemannSurface(dc,ms, stepSize = 0.01, stepType = 'RK4');


    % initialize some functions
        
        func0 = @(x,y)  1*exp(-1./(1-min((x-0).^2*9 + (y-0.3).^2*9,1)) + 1) + ...
                        0.7*exp(-1./(1-min((x+0.2).^2*3 + (y-0.1).^2*3,1)) + 1) - ...
                        1*exp(-1./(1-min((x-0.3).^2*9 + (y-0.4).^2*9,1)) + 1) - ...
                        0.7*exp(-1./(1-min((x+0.2).^2*13 + (y+0.1).^2*13,1)) + 1);
        %func0 = @(x,y)  1*exp(-((x-0).^2 + (y-0.3).^2) / 0.2) - exp(-((x-0.3).^2 + (y-0.5).^2) / 0.1);
        %func0 = @(x,y) double(((x-0.2).^2 + abs(y-0.1)) < 0.5); 
        
        
        %{
            [x,y] = ndgrid(-5:0.8:5);
            z = sin(x.^2 + y.^2) ./ (x.^2 + y.^2) .* ((x.^2 + y.^2) < 0.7);
        func0 = griddedInterpolant(x,y,z);
        %}

    % plot function under domain
        %{.
        Figureizer.figure(rsurf,func0);
        Figureizer.plotGeoFan(rsurf,0,30);
        %}




%% I0

    % initialize some betas an alphas to use
        beta = linspace(0,2*pi,200);
        alpha = linspace(-pi,pi,200)*0.5;
        [betaG,alphaG] = ndgrid(beta,alpha); % alternatively use meshgrid (doesnt work with griddedinterpolant or scatteredinterpolant)

    % run I0
        disp('I0');
        tic
            i0Data = rsurf.I0(betaG,alphaG, func0);
        toc
     
    % convert the data to a useable function
        I0f = griddedInterpolant(betaG,alphaG,i0Data);

    % plot XrayI0
        Figureizer.figureGSpace(rsurf,I0f);
        title('I0')

%{.
        
%% initialize points to transform over and precompute scattering relation
        beta = linspace(0,2*pi,120);
        alpha = linspace(-pi,pi,120+2)*0.5;   alpha = alpha(2:end-1);
        [VBeta,VAlpha] = ndgrid(beta, alpha);
        [VBetaS,VAlphaS] = rsurf.scatteringRelation(VBeta,VAlpha);
        

%% I0 Inversion (I0perp*, geoR of I0)      

        % R
            disp('R');
            tic
                %[aData,betaO,alphaO] = rsurf.geoR_precomp(VBeta,VAlpha,I0f, VBetaS,VAlphaS, 200);
                [aData,betaO,alphaO] = rsurf.geoR_simple(VBeta,VAlpha,I0f, 200);
            toc
            
            RI0f = griddedInterpolant(betaO,alphaO,aData);
            
            Figureizer.figureGSpace(rsurf,RI0f);
            title('RI0f');
        
        % I0perpstar
            % initialize points to reconstruct at (we could just re-use the points used to plot the original function)
            [VX,VY] = dc.aabbspace(100);
            [VX,VY] = ndgrid(VX,VY);

            disp('Backproject I0_perp^*');   
            rsurf.stepType = 'RK4';
            tic
                funcData = rsurf.I0perpstar(RI0f, VX,VY, 30)/(2*pi);
            toc

            reconstruct = griddedInterpolant(VX,VY,funcData);
            error = griddedInterpolant(VX,VY,abs(funcData-func0(VX,VY)));

         %plot reconstructed function
            %{.
            Figureizer.figure(rsurf,reconstruct);
            title('reconstruction');
            %}.
         %plot error (disable forceMidZero for just this)
            %{.
            %sett.forceMidZero = false;
            Figureizer.figure(rsurf,error);
            %sett.forceMidZero = true;
            title('error');
            %}.
%}

%{.
%% I0 Inversion (I0perp*, geoA*, Hilbert, geoA of I0)

    % geoA_precomp
        disp('A_-');
        tic
            [aData,betaO,alphaO] = rsurf.geoA_precomp(VBeta,VAlpha,I0f, VBetaS,VAlphaS, -1);
        toc
        
        % sorting step so that data is still in ndgrid format
        [alphaO, aind] = sort(alphaO(1,:));
        alphaO = repmat(alphaO,[height(betaO),1]);
        aData = aData(:,aind);
        
        AI0f = griddedInterpolant(betaO,alphaO,aData);
        
            %{
            Figureizer.figureGSpace(rsurf,AI0f,alphaRange = [-pi,3*pi]/2);
            title('AI0f');
            %}
        
    % geoHilbert
        disp('H');    
        tic
            hData = rsurf.geoHilbert(betaO,alphaO,AI0f, 200);
        toc
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
        tic
            [astarData,betaO,alphaO] = rsurf.geoAstar_precomp(VBeta,VAlpha,HAI0f, VBetaS,VAlphaS, 1);
        toc 
        
        % sorting step so that data is in ndgrid format
        [alphaO, aind] = sort(alphaO(1,:));
        alphaO = repmat(alphaO,[height(betaO),1]);
        astarData = astarData(:,aind);
        
        AstarHAI0f = griddedInterpolant(betaO,alphaO,astarData/4);
            
            %{.
            Figureizer.figureGSpace(rsurf,AstarHAI0f);
            title('A*HAI0f');
            %}.
      
    % I0perpstar
        % initialize points to reconstruct at (we could just re-use the points used to plot the original function)
        [VX,VY] = dc.aabbspace(100);
        [VX,VY] = ndgrid(VX,VY);
        
        disp('Backproject I0_perp^*');   
        rsurf.stepType = 'RK4';
        tic
            funcData = rsurf.I0perpstar(AstarHAI0f, VX,VY, 30)/(2*pi);
        toc
        
        reconstruct = griddedInterpolant(VX,VY,funcData);
        error = griddedInterpolant(VX,VY,abs(funcData-func0(VX,VY)));

     %plot reconstructed function
        %{.
        Figureizer.figure(rsurf,reconstruct);
        title('reconstruction');
        %}.
     %plot error (disable forceMidZero for just this)
        %{.
        %sett.forceMidZero = false;
        Figureizer.figure(rsurf,error);
        %sett.forceMidZero = true;
        title('error');
        %}.
%}     

