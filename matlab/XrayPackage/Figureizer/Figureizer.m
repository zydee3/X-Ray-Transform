classdef Figureizer
    %FIGUREIZER Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        settings = FigureizerSettings;
    end    
    
    methods (Static)
%--------------------------------------------------------------------------
%%                       Helper Stuff                                      
%--------------------------------------------------------------------------          

        function betaO = spaceBetas(surface)
            sett = Figureizer.settings;
             
            dom = surface.domain;
             
            n = sett.domainResolution + 1;
            
            switch sett.betaSpacing
                case 'default'
                    betaO = linspace(0,2*pi,n);
                case 'arclength'    
                    len = surface.BetatoS(2*pi,0.0005);
                    betaO = surface.StoBeta(linspace(0,len,n),0.0005);
                case 'euclid'    
                    tsurf = RiemannSurface(dom);
                    len = tsurf.BetatoS(2*pi,0.0005);
                    betaO = tsurf.StoBeta(linspace(0,len,n),0.0005); 
                otherwise
                    error('wrong spacing method')
            end
            
        end

%--------------------------------------------------------------------------
%%                   Deprojected Plotters                                  
%--------------------------------------------------------------------------     
        
        function fig = figureDeprojected(surface,func)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
                func (1,1) = @(x,y) NaN(size(x));
            end
            
            sett = Figureizer.settings;
            
            dom = surface.domain;
            met = surface.metric;
            
            [px, py] = met.aabbspace(dom, sett.gridResolution);
            [px, py] = ndgrid(px,py);
            [dpx,dpy,dpz] = met.deproject(px,py);
            
            fig = figure; hold on;
            ax = gca;   
                      
            %plot surface
            s = surf(dpx,dpy,dpz,func(px,py));
            light(Position = [1 1 1], Style = 'infinite', Color = 'white');
            material([0.9 1 0])
            s.EdgeColor = 'none';
            s.FaceColor = 'interp';
            s.FaceLighting = 'gouraud';
            
            %plot domain
            Figureizer.plotDomainDeprojected(surface)
            
            
            colormap(sett.funcColormap);
            
            ax.ClippingStyle = 'rectangle';
            ax.Clipping = 'off';
            axis equal;
            
            view([45,20])
            if (sett.forceMidZero)
                cmax = max(ax.CLim);
                ax.CLim = [-cmax, cmax];
            end
        end
        
        function plotDomainDeprojected(surface)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
            end
            
            sett = Figureizer.settings;

            dom = surface.domain;
            met = surface.metric;
            
            ox = dom.originX;
            oy = dom.originY;
                        
            th = Figureizer.spaceBetas(surface);
            domr = dom.bdr(th - dom.theta);
            domx = cos(th).*domr + ox;   domy = sin(th).*domr + oy;
            [ddomx,ddomy,ddomz] = met.deproject(domx,domy);
            plot3(ddomx,ddomy,ddomz,sett.domainStyle,Color = sett.domainColor, LineWidth = sett.domainThickness)
            %compute alnormal spacing
            ath = th;
            %compute alnormals
            ant = dom.alNormal(ath) + dom.theta;
            [dax,day,daz] = met.deproject(domx + 0.01*cos(ant),domy + 0.01*sin(ant));
            
            %compute sizes for positioning camera and rescaling alnormal
            adx = [min(ddomx(1:ceil(end/40):end),[],'all'),max(ddomx(1:ceil(end/40):end),[],'all')]; 
            ady = [min(ddomy(1:ceil(end/40):end),[],'all'),max(ddomy(1:ceil(end/40):end),[],'all')]; 
            adz = [min(ddomz(1:ceil(end/40):end),[],'all'),max(ddomz(1:ceil(end/40):end),[],'all')]; 
            ab = max(abs([adx(2)-adx(1), ady(2)-ady(1), adz(2)-adz(1)]) ) * 0.03;
            
            %rescale alnormals
            mag = sqrt((dax-ddomx).^2+(day-ddomy).^2+(daz-ddomz).^2);
            dax = ab*(dax-ddomx)./mag + ddomx;
            day = ab*(day-ddomy)./mag + ddomy;
            daz = ab*(daz-ddomz)./mag + ddomz;
            %plot alnormals
            plot3([ddomx; dax], [ddomy; day], [ddomz; daz],...
                                     sett.domainStyle, Color = sett.domainColor, LineWidth = sett.domainThickness);
        end
        
        
        function plotGeoDeprojected(surface, X,Y,Th, style,args)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                Th {mustBeNumeric} = []
                style {mustBeText} = '-'
                args.penColor = Figureizer.settings.penColor
            end
            
            % reshape for plotting (not totally neccessary)
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            [xp,yp,~] = surface.geodesic(X,Y,Th);
            
            Figureizer.plotDeprojected(surface,xp,yp,style,penColor = args.penColor);
        end
        
        function plotGeoFanDeprojected(surface, beta, numgeos)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                beta (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 0
            end

            dom = surface.domain;
            met = surface.metric;
            ra = dom.bdr(beta - dom.theta);
            x = cos(beta) * ra + dom.originX;
            y = sin(beta) * ra + dom.originY;
            
            Th = linspace(0.5*pi, 1.5*pi, numgeos) + dom.alNormal(beta) + dom.theta;
            X = ones(size(Th)) * x;
            Y = ones(size(Th)) * y;
            
            Figureizer.plotGeoDeprojected(surface,X,Y,Th);            
            
        end
        
        function plotGeoNormalsDeprojected(surface, th, numgeos)

            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                th (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 0                
            end
            
            beta = linspace(0,2*pi, numgeos+1); beta = beta(1:end-1);
            dom = surface.domain;
            ra = dom.bdr(beta - dom.theta);
            X = cos(beta) .* ra + dom.originX;
            Y = sin(beta) .* ra + dom.originY;
            Th = th + dom.alNormal(beta) + pi;        
            
            Figureizer.plotGeoDeprojected(surface,X,Y,Th);            
        end

        function plotGeoParallelsDeprojected(surface, th, Beta)

            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                th (1,1) {mustBeNumeric} = 0
                Beta {mustBeNumeric} = linspace(pi,1.25*pi,40)     
            end
            
            dom = surface.domain;
            ra = dom.bdr(Beta - dom.theta);
            X = cos(Beta) .* ra + dom.originX;
            Y = sin(Beta) .* ra + dom.originY;
            Th = th;
            
            Figureizer.plotPoint(surface,X,Y);
            holdBool = ishold;
            hold on;         
            
            Figureizer.plotGeoDeprojected(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
        end

        
        
        function plotDeprojected(surface, X,Y, style,args)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                style {mustBeText} = Figureizer.settings.penStyle
                args.penColor = Figureizer.settings.penColor
                args.lineWidth = Figureizer.settings.penThickness
            end
            
            [dx,dy,dz] = surface.metric.deproject(X,Y);
            plot3(dx,dy,dz,style,Color = args.penColor,LineWidth = args.lineWidth)
            
        end
        
        function plotBdrPointDeprojected(surface, Beta)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                Beta {mustBeNumeric} = []
            end
            
            [xO,yO,~] = surface.BAtoXYTh(Beta,zeros(size(Beta)));
            
            plotDeprojected(surface, xO,yO);
        end
         
%--------------------------------------------------------------------------
%%                    Standard Plotters                                    
%--------------------------------------------------------------------------          
        
        function fig = figure(surface,func)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
                func (1,1) = @(x,y) NaN(size(x));
            end
            
            sett = Figureizer.settings;

            fig = figure; hold on;
            ax = gca;
            
            dom = surface.domain;
            met = surface.metric;
            
            [px,py] = dom.aabbspace(sett.gridResolution);
            [px, py] = ndgrid(px,py);
            ox = dom.originX;
            oy = dom.originY;
                       
            %plot function
            %[px,py] = meshgrid(gx,gy);
            s = pcolor(px,py,func(px,py));
            s.EdgeColor = 'none';
            s.FaceColor = 'interp';
            
            %plot domain
            Figureizer.plotDomain(surface);
            
            axis equal;            
            colormap(sett.funcColormap);
            
            %set lims
            [minB,maxB] = dom.getBoundingBox;
            xlim([minB(1),maxB(1)] + dom.originX);
            ylim([minB(2),maxB(2)] + dom.originY);
            
            if (sett.forceMidZero)
                cmax = max(ax.CLim);
                ax.CLim = [-cmax, cmax];
            end
        end
        
        function plotDomain(surface)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
            end
            
            sett = Figureizer.settings;

            dom = surface.domain;
            met = surface.metric;
            
            ox = dom.originX;
            oy = dom.originY;
            
            th = Figureizer.spaceBetas(surface);
            domr = dom.bdr(th - dom.theta);
            domx = cos(th).*domr + ox;   domy = sin(th).*domr + oy;
            Figureizer.plot(surface,domx,domy,sett.domainStyle,penColor = sett.domainColor, lineWidth = sett.domainThickness)
            %compute alnormal spacing
            ath = th;
            %compute alnormals
            ant = dom.alNormal(ath) + dom.theta;
            ax = domx + cos(ant);   ay = domy + sin(ant);
            
            %compute sizes for positioning camera and rescaling alnormal
            adx = [min(domx(1:ceil(end/40):end),[],'all'),max(domx(1:ceil(end/40):end),[],'all')]; 
            ady = [min(domy(1:ceil(end/40):end),[],'all'),max(domy(1:ceil(end/40):end),[],'all')]; 
            ab = max(abs([adx(2)-adx(1), ady(2)-ady(1), adx(2)-ady(1)]) ) * 0.03;
            
            %rescale alnormals
            ax = ab*(ax-domx) + domx;
            ay = ab*(ay-domy) + domy;
            %plot alnormals
            Figureizer.plot(surface,[domx; ax], [domy; ay],...
                     sett.domainStyle, penColor = sett.domainColor, lineWidth = sett.domainThickness);
        end
       
        
        function plotGeo(surface, X,Y,Th, style,args)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface} 
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                Th {mustBeNumeric} = []
                style {mustBeText} = '-'
                args.penColor = Figureizer.settings.penColor
            end

            X = X(:);   Y = Y(:);   Th = Th(:);
            
            [xO,yO,~] = surface.geodesic(X,Y,Th);

            Figureizer.plot(surface,xO,yO,style,penColor = args.penColor);
        end   
 
        function plotGeoFan(surface, beta, numgeos)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface} 
                beta (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 1
            end

            dom = surface.domain;
            ra = dom.bdr(beta - dom.theta);
            x = cos(beta) * ra + dom.originX;
            y = sin(beta) * ra + dom.originY;
            
            Th = linspace(0.5*pi, 1.5*pi, numgeos) + dom.alNormal(beta) + dom.theta;
            X = ones(1,numgeos) * x;
            Y = ones(1,numgeos) * y;
            
            Figureizer.plotGeo(surface,X,Y,Th);            
            
        end   
        
        function plotGeoNormals(surface, th, numgeos)

            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                th (1,1) {mustBeNumeric} = 0
                numgeos (1,1) {mustBeNumeric} = 0                
            end
            
            beta = linspace(0,2*pi, numgeos+1); beta = beta(1:end-1);
            dom = surface.domain;
            ra = dom.bdr(beta - dom.theta);
            X = cos(beta) .* ra + dom.originX;
            Y = sin(beta) .* ra + dom.originY;
            Th = th + dom.alNormal(beta) + pi + dom.theta;

            Figureizer.plotGeo(surface,X,Y,Th);            
        end
        
        function plotGeoParallels(surface, th, Beta)

            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                th (1,1) {mustBeNumeric} = 0
                Beta {mustBeNumeric} = linspace(pi,1.25*pi,40)     
            end
            
            dom = surface.domain;
            ra = dom.bdr(Beta - dom.theta);
            X = cos(Beta) .* ra + dom.originX;
            Y = sin(Beta) .* ra + dom.originY;
            Th = th;
            
            Figureizer.plotPoint(surface,X,Y);
            holdBool = ishold;
            hold on;         
            
            Figureizer.plotGeo(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
        end
        
        
        function plot(surface, X,Y, style,args)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                style {mustBeText} = Figureizer.settings.penStyle
                args.penColor = Figureizer.settings.penColor
                args.lineWidth = Figureizer.settings.penThickness
            end
            
            plot(X,Y,style,Color = args.penColor, LineWidth = args.lineWidth);
        end    

        function plotBdrPoint(surface, Beta)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                Beta {mustBeNumeric} = []
            end
            
            [xO,yO,~] = surface.BAtoXYTh(Beta,zeros(size(Beta)));
                        
            plotPoint(surface, xO,yO);
        end

%--------------------------------------------------------------------------              
%%                    Geo-Space Plotters                                   
%--------------------------------------------------------------------------

        function fig = figureGSpace(surface, func, args)
            
            arguments
                surface
                func (1,1) = @(x,y) NaN(size(x));
                args.betaRange (1,2) {mustBeNumeric} = [0,2*pi];
                args.alphaRange (1,2) {mustBeNumeric} = [-pi,pi]/2;
            end
            
            sett = Figureizer.settings;

            fig = figure; hold on;
            ax = gca;   

            betaS  = linspace(args.betaRange(1), args.betaRange(2), sett.gridResolution);
            alphaS = linspace(args.alphaRange(1),args.alphaRange(2),sett.gridResolution);
            
            [betaG,alphaG] = ndgrid(betaS,alphaS);
            
            s = pcolor(betaG,alphaG,func(betaG,alphaG));
            s.EdgeColor = 'none';
            s.FaceColor = 'interp';
            colormap(sett.funcColormap);
            
            xlim(args.betaRange);
            ylim(args.alphaRange);
            
            if (sett.forceMidZero)
                cmax = max(ax.CLim);
                ax.CLim = [-cmax, cmax];
            end
            
        end


%--------------------------------------------------------------------------
    end
end