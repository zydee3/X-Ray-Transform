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

%--------------------------------------------------------------------------
%%                   Deprojected Plotters
%--------------------------------------------------------------------------     
        
        function fig = figureDeprojected(surface,func)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
                func (1,1) {mustBe2DHandle} = @(x,y) zeros(size(x));
            end
            
            sett = Figureizer.settings;
            
            dom = surface.domain;
            met = surface.metric;
            
            [px, py] = met.aabbspace(dom, sett.gridResolution);
            [px, py] = meshgrid(px,py);
            [dpx,dpy,dpz] = met.deproject(px,py);
            
            fig = figure; hold on;
            
            ox = dom.originX;
            oy = dom.originY;
                       
            %plot surface
            surf(dpx,dpy,dpz,func(px,py),...
                FaceColor = 'interp', EdgeColor = 'none', FaceLighting = 'gouraud');
            light(Position = [1 1 1], Style = 'infinite', Color = 'white');
            material([0.7 0.3 0])
            
            %plot domain
            n = 250;
            th = linspace(0,2*pi,n);
            domr = dom.bdr(th - dom.theta);
            domx = cos(th).*domr + ox;   domy = sin(th).*domr + oy;
            [ddomx,ddomy,ddomz] = met.deproject(domx,domy);
            plot3(ddomx,ddomy,ddomz,Color = sett.domainColor)
            %compute alnormals
            ant = dom.alNormal(th) + dom.theta;
            [dax,day,daz] = met.deproject(domx + 0.01*cos(ant),domy + 0.01*sin(ant));
            
            %compute sizes for positioning camera and rescaling alnormal
            adx = [min(ddomx(1:ceil(end/40):end),[],'all'),max(ddomx(1:ceil(end/40):end),[],'all')]; 
            ady = [min(ddomy(1:ceil(end/40):end),[],'all'),max(ddomy(1:ceil(end/40):end),[],'all')]; 
            adz = [min(ddomz(1:ceil(end/40):end),[],'all'),max(ddomz(1:ceil(end/40):end),[],'all')]; 
            ab = max(abs([adx(2)-adx(1), ady(2)-ady(1), adx(2)-ady(1)]) ) * 0.2;
            
            %rescale alnormals
            mag = sqrt((dax-ddomx).^2+(day-ddomy).^2+(daz-ddomz).^2) * 4;
            dax = ab*(dax-ddomx)./mag + ddomx;
            day = ab*(day-ddomy)./mag + ddomy;
            daz = ab*(daz-ddomz)./mag + ddomz;
            %plot alnormals
            plot3([ddomx; dax],...
                 [ddomy; day],...
                 [ddomz; daz],...
                 Color = sett.domainColor);
            
            
            colormap(sett.funcColormap);
            ax = gca;   
            ax.ClippingStyle = 'rectangle';
            ax.Clipping = 'off';
            axis equal;
            
            % position camera
            switch sett.plotCenterType
                case 'default'
                    %\do nothing lol
                case 'domain'
                    xlim(adx + [-ab,ab]);
                    ylim(ady + [-ab,ab]);
                    zlim(adz + [-ab,ab]);
            end    
            view([45,20])
            
        end

        
        function plotGeoDeprojected(surface, X,Y,Th, style,args)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                Th {mustBeNumeric} = []
                style {mustBeText} = Figureizer.settings.linestyle
                args.penColor = Figureizer.settings.penColor
            end
            
            % reshape for plotting (not totally neccessary)
            X = X(:);   Y = Y(:);   Th = Th(:);
            
            [xp,yp,~] = surface.geodesic(X,Y,Th);
            [dx,dy,dz] = surface.metric.deproject(xp,yp);
            
            plot3(dx,dy,dz,style,Color = args.penColor);
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
            
            [dx,dy,dz] = met.deproject(x,y);
            plot3(dx,dy,dz,'r*')
            
            holdBool = ishold;
            hold on;
            
            Th = linspace(0.5*pi, 1.5*pi, numgeos) + dom.alNormal(beta) + dom.theta;
            X = ones(size(Th)) * x;
            Y = ones(size(Th)) * y;
            
            Figureizer.plotGeoDeprojected(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
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
            
            Figureizer.plotPointDeprojected(surface,X,Y);
            holdBool = ishold;
            hold on;         
            
            Figureizer.plotGeoDeprojected(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
        end
        
        
        function plotPointDeprojected(surface, X,Y, style,args)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                style {mustBeText} = '*'
                args.penColor = Figureizer.settings.penColor
            end
            
            [dx,dy,dz] = surface.metric.deproject(X,Y);
            plot3(dx,dy,dz,style,Color = args.penColor)
            
        end
        
        function plotBdrPointDeprojected(surface, Beta)
            
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                Beta {mustBeNumeric} = []
            end
            
            [xO,yO,~] = surface.BAtoXYTh(Beta,zeros(size(Beta)));
            
            plotPointDeprojected(surface, xO,yO);
        end
 
        
        
%--------------------------------------------------------------------------
%%                   Standard Plotters
%--------------------------------------------------------------------------          
        
        function fig = figure(surface,func)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface};
                func (1,1) {mustBe2DHandle} = @(x,y) NaN(size(x));
            end
            
            sett = Figureizer.settings;

            fig = figure; hold on;
            
            dom = surface.domain;
            met = surface.metric;
            
            [px,py] = dom.aabbspace(sett.gridResolution);
            [px, py] = meshgrid(px,py);
            x0 = dom.originX;
            y0 = dom.originY;
                       
            %plot function
            %[px,py] = meshgrid(gx,gy);
            s = pcolor(px,py,func(px,py));
            s.EdgeColor = 'none';
            
            %plot domain
            dom.plotAlNormal
            
            axis equal;            
            colormap(sett.funcColormap);

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

            plot(xO,yO,style,Color = args.penColor);
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
            
            plot(x,y,'*',Color = Figureizer.geoColor)
            
            holdBool = ishold;
            hold on;
            
            Th = linspace(0.5*pi, 1.5*pi, numgeos) + dom.alNormal(beta) + dom.theta;
            X = ones(1,numgeos) * x;
            Y = ones(1,numgeos) * y;
            
            Figureizer.plotGeo(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
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
            
            Figureizer.plotPoint(surface,X,Y);
            holdBool = ishold;
            hold on;         
            
            Figureizer.plotGeo(surface,X,Y,Th);            
            
            if (~holdBool), hold off; end
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
        
        
        function plotPoint(surface, X,Y, style,args)
            arguments
                surface (1,1) {RiemannSurface.mustBeSurface}
                X {mustBeNumeric} = []
                Y {mustBeNumeric} = []
                style {mustBeText} = '*'
                args.penColor = Figureizer.settings.penColor
            end
            
            plot(X,Y,style,Color = args.penColor);
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

    end
end

