clc, close all, clear


%% initialization

    % begin by initializing a surface (see EXAMPLECODE_initializeSurface)
        R = pi/3;
        dc = circleDomain(radius = R);
        ms = euclidMetric;
        rsurf = RiemannSurface(dc,ms, stepSize = 0.1, stepType = 'RK4');


    % initialize some functions
        
        func_uncropped = @(x,y) 1;
        func_cropped = @(x,y) double((x.^2 + y.^2) < R.^2); 
        truexray = @(beta,alpha) 2.*sqrt( R*R - R*R*sin(...
                               max( min(abs(alpha),asin(1)) ,0)...
                                             ).^2 );



%% I0

    % initialize some betas an alphas to use
        beta = linspace(0,2*pi,200);
        alpha = linspace(-pi,pi,200)*0.5;
        [betaG,alphaG] = ndgrid(beta,alpha); % alternatively use meshgrid (doesnt work with griddedinterpolant or scatteredinterpolant)

    % run I0
        tic
        i0Data0 = rsurf.I0(betaG,alphaG, func_uncropped);
        toc
                
    % run I0 on the cropped integrand
        tic
        i0Data1 = rsurf.I0(betaG,alphaG, func_cropped);
        toc  
       
    % run I0_interp, slinear
        dc.exitInterpType = 'slinear';
        rsurf.domain=dc;
        tic
        i0Data2 = rsurf.I0_scatt(betaG,alphaG, func_uncropped);
        toc
        
    % run I0_interp, slinearC
        dc.exitInterpType = 'slinearC';
        rsurf.domain=dc;
        tic
        i0Data3 = rsurf.I0_scatt(betaG,alphaG, func_uncropped);
        toc    
               
    % construct figures
        figure; sgtitle('I0 errors, alpha on [-pi/2,pi/2]')
            subplot(1,4,1);
            pl = pcolor(betaG,alphaG, abs(i0Data0 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('original I0');

            subplot(1,4,2);
            pl = pcolor(betaG,alphaG, abs(i0Data1 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 with cropped integrand')

            subplot(1,4,3);
            pl = pcolor(betaG,alphaG, abs(i0Data2 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 slinear')
            
            subplot(1,4,4);
            pl = pcolor(betaG,alphaG, abs(i0Data3 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 slinearB')

            
        
        figure; sgtitle('I0 errors, alpha restricted')
            alind = find(i0Data2(1,:)>rsurf.stepSize);
            betaG = betaG(:,alind); alphaG = alphaG(:,alind);
            i0Data0 = i0Data0(:,alind);
            i0Data1 = i0Data1(:,alind);
            i0Data2 = i0Data2(:,alind);
            i0Data3 = i0Data3(:,alind);
            
            subplot(1,4,1);
            pl = pcolor(betaG,alphaG, abs(i0Data0 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('original I0');

            subplot(1,4,2);
            pl = pcolor(betaG,alphaG, abs(i0Data1 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 with cropped integrand')

            subplot(1,4,3);
            pl = pcolor(betaG,alphaG, abs(i0Data2 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 slinear')
            
            subplot(1,4,4);
            pl = pcolor(betaG,alphaG, abs(i0Data3 - truexray(betaG,alphaG)));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 slinearB')
            
            %{
        figure; sgtitle('I0 errors, alpha on [-pi/8,pi/8]')
            hal = height(alphaG);
            alind = floor(3*hal/8):floor(5*hal/8);
        
            subplot(1,3,1);
            pl = pcolor(betaG(:,alind),alphaG(:,alind), abs(i0Data0(:,alind) - truexray(betaG(:,alind),alphaG(:,alind))));
            pl.EdgeColor = 'none';
            colorbar;
            title('original I0');

            subplot(1,3,2);
            pl = pcolor(betaG(:,alind),alphaG(:,alind), abs(i0Data1(:,alind) - truexray(betaG(:,alind),alphaG(:,alind))));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 with cropped integrand')

            subplot(1,3,3);
            pl = pcolor(betaG(:,alind),alphaG(:,alind), abs(i0Data2(:,alind) - truexray(betaG(:,alind),alphaG(:,alind))));
            pl.EdgeColor = 'none';
            colorbar;
            title('I0 with interpolated exits')            
%}
           
