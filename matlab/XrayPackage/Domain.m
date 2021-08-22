classdef Domain
    %DOMAIN Parent class of all other domains
    
    properties
        originX (1,1) {mustBeNumeric} = 0
        originY (1,1) {mustBeNumeric} = 0
        theta (1,1) {mustBeNumeric} = 0
        exitInterpType (1,:) {mustBeText} = 'last'
    end
    
    properties (Access = 'private')
        bdrVAR = @(th) ones(size(th));
    end    
        
    
    
    methods
        
%--------------------------------------------------------------------------
%%                 Constructor, Characteristic Functions
%--------------------------------------------------------------------------  

        function out = bdr(obj,th)
            %BDR Computes a position of a point of the boundry in polar
            %coordinates.
            %   Method uses the function handle obj.bdrVAR if not 
            %   overwritten by a subclass. 
            out = obj.bdrVAR(th);
        end
        
        function out = dbdr(obj,th)
            %DBDR Derivative of obj.bdr.
            %   Method computes the derivative numerically from obj.bdr if not
            %   overwritten by a subclass.
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = deriv(@(t) obj.bdr(t), th);
        end
        
        function out = ddbdr(obj,th)
            %DDBDR Second derivative of obj.bdr.
            %   Method computes the derivatives numerically from obj.bdr if not
            %   overwritten by a subclass.
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = dderiv(@(t) obj.bdr(t), th);
        end
        
                
        function obj = Domain(bdr)
        	%DOMAIN Constructs an instance of Domain.           
        if (strcmp(class(obj), 'Domain')) %%% -- !! TODO: hack to avoid subclass calling superclass constructor, should prolly change !!    
            if (nargin == 1 && ...
                isa(bdr,'function_handle') && nargin(lg) == 1)
                
                obj.bdrVAR = bdr;
            end
            
        end,end   
   
    
%--------------------------------------------------------------------------
%%                                  Misc
%--------------------------------------------------------------------------      

       
        function minR = getMinRadius(obj) 
            %GETMINRADIUS Returns the minimum value of obj.bdr on the
            %interval [0,2pi].
            %   Commputation is numerical and potentially inaccurate if not
            %   overwritten by a subclass.
            %   See obj.isInsideR2, getBoundingBox.
            
            %warning('Method "getMinRadius" is slow, consider overriding.'); % see comment for getBoundingBox
            minR = min(obj.bdr(linspace(0,2*pi, 1000)));
        end 
        
        
        function out = alNormal(obj, Th)
            %ALNORMAL Computes the angle of the outer normal of the domain.
            Th = Th - obj.theta;
            out = 0.5*pi - atan2(obj.bdr(Th).*cos(Th) + obj.dbdr(Th).*sin(Th),...
                        (obj.bdr(Th).*sin(Th) - obj.dbdr(Th).*cos(Th))); 
        end  
        
                              
        function bool = isInside(obj, X, Y, ~)
            %ISINSIDE Determines if a sequence of points is inside the domain.
            
            %%{
            X = X - obj.originX;
            Y = Y - obj.originY;
            r = obj.bdr(atan2(Y,X) - obj.theta);
            bool = (X).*(X) + (Y).*(Y) <= r.*r;            
            %}
            %{
            x = obj.originX;
            y = obj.originY;
            bool = (X-x).*(X-x) + (Y-y).*(Y-y) <= obj.bdr(atan2(Y-y,X-x) - obj.theta).^2;            
            %}
        end  
                
        function bool = isInsideR2(obj, X, Y, r)
            %ISINSIDE An optimized approach of determining if a sequence of 
            %points is inside the domain that takes a minimum radius as an 
            %additional argument.
            %   See obj.isInside.
            
            X = X - obj.originX;
            Y = Y - obj.originY;
            XY2 = X.*X + Y.*Y;
            bool = XY2 <= r;
            
            if (any(~bool))
                r = obj.bdr(atan2(Y(bool),X(bool)) - obj.theta);
                bool(bool) = XY2(bool) <= r.*r;
            end
        end
        
        function [betaO,alphaO, tO] = exitInterp(obj, xin,yin, xout,yout)
            % constructs betaO,alphaO from an interpolant described by a
            % pair of points.

            switch obj.exitInterpType
                case 'last'
                    tO = 1;
                    betaO = atan2(yout,xout);
                    
                case 'slinear'
                    mIn = sqrt(xin.*xin + yin.*yin);
                    mOut = sqrt(xout.*xout + yout.*yout);
                    [btIn,btOut] = atan2near(yin, xin, yout, xout);
                    
                    funcin = -mIn+obj.bdr(btIn);
                    tO = (funcin)./( mOut-obj.bdr(btOut) +funcin);

                    betaO = mod((btOut-btIn).*tO + btIn,2*pi);
                case 'slinearB'
                    weight = 0.4;
                    xint = (xout-xin)*weight+xin; yint = (yout-yin)*weight+yin;
                    
                    mIn = sqrt(xint.*xint + yint.*yint);
                    mOut = sqrt(xout.*xout + yout.*yout);
                    [btIn,btOut] = atan2near(yint, xint, yout, xout);
                    funcin = -mIn+obj.bdr(btIn);
                    
                    tO = (funcin)./( mOut-obj.bdr(btOut) +funcin);

                    betaO = mod((btOut-btIn).*tO + btIn,2*pi);    
                    alphaO = -obj.alNormal(betaO) + atan2(yout-yin,xout-xin) + pi;
                    
                    tO = (tO.*(1-weight) + weight);
                    return;
                        
                case 'slinearC'
                    weight = 0.273;
                    xint = (xout-xin)*weight+xin; yint = (yout-yin)*weight+yin;
                    xout = (xin-xout)*weight+xout; yout = (yin-yout)*weight+yout;
                    
                    mIn = sqrt(xint.*xint + yint.*yint);
                    mOut = sqrt(xout.*xout + yout.*yout);
                    [btIn,btOut] = atan2near(yint, xint, yout, xout);
                    funcin = -mIn+obj.bdr(btIn);
                    
                    tO = (funcin)./( mOut-obj.bdr(btOut) +funcin);

                    betaO = mod((btOut-btIn).*tO + btIn,2*pi);    
                    alphaO = -obj.alNormal(betaO) + atan2(yout-yin,xout-xin) + pi;
                    
                    tO = (tO.*(1-2*weight) + weight);
                    return;
                    
                otherwise 
                    error('wrong interpolation method')
            end
            alphaO = -obj.alNormal(betaO) + atan2(yout-yin,xout-xin) + pi;
        end
        
        
        function [minB,maxB] = getBoundingBox(obj) 
            %GETBOUNDINGBOX Returns two vectors representing the lower
            %right and upper left corners of the domain's axis-aligned
            %bounding box relative to the domain's origin.
            %   Commputation is numerical and potentially inaccurate if not
            %   overwritten by a subclass.
            %   See obj.getMinRadius.
            
            
            %warning('Method "getBoundingBox" is slow, consider overriding.'); % consider implementing a faster numerical update to this method, something with newtons itter could work
            rSamples = obj.bdr(linspace(0,2*pi, 1000) - obj.theta);
            eSamples = cos(linspace(0,2*pi, 1000)).* rSamples;
            minB = [min(eSamples),0];
            maxB = [max(eSamples),0];
            eSamples = sin(linspace(0,2*pi, 1000)).* rSamples;
            minB(2) = min(eSamples);
            maxB(2) = max(eSamples);
        end 

        function [xO,yO] = aabbgrid(obj,resoX,resoY)
            %GRIDAABB generates an input for meshgrid or ngrid using obj.getBoundingBox
            arguments
                obj
                resoX (1,1) {mustBeNumeric} = 250
                resoY (1,1) {mustBeNumeric} = resoX
            end

            [minB,maxB] = obj.getBoundingBox; 
            xO = linspace(minB(1),maxB(1),resoX);
            yO = linspace(minB(2),maxB(2),resoY);
        end

        
        function obj = transform(obj, xoff, yoff, rot)
            %TRANSFORM A method to set the properties originX, originY,
            %theta all at once.
            arguments
                obj
                xoff {mustBeNumeric}
                yoff {mustBeNumeric}
                rot {mustBeNumeric}
            end

            obj.originX = xoff;
            obj.originY = yoff;
            obj.theta = rot;
        end    
        
        
        function [bdr,dbdr,ddbdr] = getHandles(obj)
            %METRICVALSCURV Returns bdr, dbdr, and ddbdr as function
            %handles.
            bdr = @(t) obj.bdr(t);
            dbdr = @(t) obj.dbdr(t);
            ddbdr = @(t) obj.ddbdr(t);
        end 
           
        
%--------------------------------------------------------------------------
%%                                Plotters
%--------------------------------------------------------------------------              
    
        function out = plot(obj)
            %PLOT Displays the boundry using obj.bdr.
            
            n = 500;
            th = linspace(0,2*pi,n);
            th0 = obj.theta;
            x0 = obj.originX;
            y0 = obj.originY;
            r = obj.bdr(th - th0);
            pointX = cos(th) .* r + x0;
            pointY = sin(th) .* r + y0;
            
            out = plot(pointX,pointY,'b');
        end        
  
        
        function out = plotBdrPoint(obj, Beta)
            %PLOT Displays the boundry using obj.bdr.
            
            arguments
                obj
                Beta {mustBeNumeric} = []
            end

            Beta = Beta(:);

            n = 500;
            th0 = obj.theta;
            x0 = obj.originX;
            y0 = obj.originY;
            r = obj.bdr(Beta - th0);
            pointX = cos(Beta) .* r + x0;
            pointY = sin(Beta) .* r + y0;
            
            out = plot(pointX,pointY,'r*');
        end        
        
        function out = plotAABB(obj)
            %PLOTAABB Displays the axis aligned bounding box of the domain
            %using obj.getBoundingBox.
           
            [minB,maxB] = obj.getBoundingBox();
            pointX = [maxB(1),maxB(1),minB(1),minB(1),maxB(1)] + obj.originX;
            pointY = [maxB(2),minB(2),minB(2),maxB(2),maxB(2)] + obj.originY;
            out = plot(pointX,pointY,'b');
        end
               
        function out = plotOrigin(obj)
            %PLOTORIGIN Plots a point at the origin of the domain.
            out = plot(obj.originX,obj.originY,'b*');
        end 
        
        function out = plotAlNormal(obj)
            %PLOTALNORMAL Plots the domain and the outer normals of the domain
            %described by obj.alNormal.
            
            holdBool = ishold;
            hold on;
            
            n = 250;
            th = linspace(0,2*pi,n);
            th0 = obj.theta;
            x0 = obj.originX;
            y0 = obj.originY;
            r = obj.bdr(th - th0);
            an = obj.alNormal(th) + th0;
            pointX = cos(th) .* r + x0;
            pointY = sin(th) .* r + y0;
            [minB,maxB] = obj.getBoundingBox();
            als = 0.012*sum(abs(maxB-minB));
                        
            plot([pointX; pointX + als*cos(an)],...
                 [pointY; pointY + als*sin(an)],...
                 'b');
            
            plot(pointX,pointY,'b');
            if (~holdBool), hold off; end;
        end
                
        function out = plotIsInsideTest(obj)
            %PLOTISINSIDETEST Plots a random collection of points within 
            %obj.getBoundingBox to test IsInside and IsInsideR2 methods.
            %   Red points indicate points outside the domain, while green
            %   points indicate inside points.
                        
            holdBool = ishold;
            hold on;
            

            [minB,maxB] = obj.getBoundingBox();
            minR2 = obj.getMinRadius();
            minR2 = minR2*minR2;
            pointX = rand(1,5000)*(maxB(1)-minB(1)) + minB(1) + obj.originX;
            pointY = rand(1,5000)*(maxB(2)-minB(2)) + minB(2) + obj.originY;       
            inside = obj.isInsideR2(pointX, pointY, minR2);
            
            
            plot(pointX(find(inside)),pointY(find(inside)),'g.');
            plot(pointX(find(~inside)),pointY(find(~inside)),'r.');
            if (~holdBool), hold off; end;
        end        
        
        function out = plotALL(obj)
            %PLOTALL A convinent function that runs a menagerie of other
            %plotting functions.
            %   Intended to test that everything is working as it should.
            
            holdBool = ishold;
            hold on;
            
            %obj.plotIsInsideTest();
            obj.plotAABB();
            obj.plotAlNormal();
            obj.plotOrigin();
            
            if (~holdBool), hold off; end;
        end
        
       
    end
    
    methods (Static)
        
        function mustBeDomain(obj)
           %MUSTBEDOMAIN Errors if the passed object is not an instance of 
           %Domain or of a subclass of Domain.
           if (~isa(obj,'Domain'))
               error("Value must be a Domain.");
           end
        end     
        
    end
    
end