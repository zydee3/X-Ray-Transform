classdef Domain
    %DOMAIN Parent class of all other domains
    
    properties
        originX = 0; % Handle these values after construction
        originY = 0;
        theta = 0;
    end
    
    properties (Access = 'private')
        bdrVAR = @(th) ones(size(th));
    end    
    
    properties (SetAccess = 'protected')
        %rMax; % now unused
    end    
    
    
    
    methods
        function out = bdr(obj,th)
            out = obj.bdrVAR(th);
        end
        
        function out = dbdr(obj,th)
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = deriv(@(t) obj.bdr(t), th);
        end
        
        function out = ddbdr(obj,th)
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = dderiv(@(t) obj.bdr(t), th);
        end
        
        
        
        function obj = Domain(bdr)
        if (strcmp(class(obj), 'Domain')) %%% -- !! TODO: hack to avoid subclass calling superclass constructor, should prolly change !!    
            if (nargin == 1 && ...
                isa(bdr,'function_handle') && nargin(lg) == 1)
                
                obj.bdrVAR = bdr;
                %obj.rMax = rMax;
            end
            
        end,end   
    
        
    
        function [bdr,dbdr,ddbdr] = getHandles(obj)
            bdr = @(t) obj.bdr(t);
            dbdr = @(t) obj.dbdr(t);
            ddbdr = @(t) obj.ddbdr(t);
        end 
        
        
        function out = alNormal(obj, th)
            th = th - obj.theta;
            out = 0.5*pi - atan2(obj.bdr(th).*cos(th) + obj.dbdr(th).*sin(th),...
                        (obj.bdr(th).*sin(th) - obj.dbdr(th).*cos(th))); 
        end  
        
        
        function [minB,maxB] = getBoundingBox(obj) 
            %warning('Method "getBoundingBox" is slow, consider overriding.'); % consider implementing a faster numerical update to this method, something with newtons itter could work
            rSamples = obj.bdr(linspace(0,2*pi, 1000) - obj.theta);
            eSamples = cos(linspace(0,2*pi, 1000)).* rSamples;
            minB = [min(eSamples),0];
            maxB = [max(eSamples),0];
            eSamples = sin(linspace(0,2*pi, 1000)).* rSamples;
            minB(2) = min(eSamples);
            maxB(2) = max(eSamples);
        end 
        
        
        function minR = getMinRadius(obj) 
            %warning('Method "getMinRadius" is slow, consider overriding.'); % see comment for getBoundingBox
            minR = min(obj.bdr(linspace(0,2*pi, 1000)));
        end 
              
        
        function bool = isInside(obj, X, Y, ~)
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
            X = X- obj.originX;
            Y = Y - obj.originY;
            XY2 = X.*X + Y.*Y;
            bool = XY2 <= r;
            
            if (any(~bool))
                i = find(~bool);
                r = obj.bdr(atan2(Y(i),X(i)) - obj.theta);
                bool(i) = XY2(i) <= r.*r;
            end
        end
        
    
        
        function obj = transform(obj, argA, argB, argC) % !! TODO: This bit of code is redudnant with the constructor, fix !!
            if (nargin == 1 && isnumeric(argA))   
                obj.theta = argA;
                return;   
            elseif (nargin == 2 && isnumeric(argA) && isnumeric(argB))   
                obj.originX = argA;
                obj.originY = argB;
                return;
            elseif (isnumeric(argA) && isnumeric(argB) && isnumeric(argC))       
                obj.theta = argA;
                obj.originX = argB;
                obj.originY = argC;
                return;
            end
            error("Bad input arguments.")
        end    
            
        
    
        function out = plot(obj) % !! TODO: make nicer lool!!
            %plot Displays boundry
            
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
        
        
        function out = plotAABB(obj)
            %plot Displays the axis aligned bounding box of the domain
           
            [minB,maxB] = obj.getBoundingBox();
            pointX = [maxB(1),maxB(1),minB(1),minB(1),maxB(1)] + obj.originX;
            pointY = [maxB(2),minB(2),minB(2),maxB(2),maxB(2)] + obj.originY;
            out = plot(pointX,pointY,'b');
        end
        
        
        function out = plotOrigin(obj) % !! TODO: make nicer lool!!
            %plot Displays the axis aligned bounding box of the domain
            out = plot(obj.originX,obj.originY,'b*');
        end
        
        function out = plotAlNormal(obj) % !! TODO: make nicer lool!!
            %plot Plots with alNorm
            
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
        
        
        function out = plotIsInsideTest(obj) % !! TODO: make nicer lool!!
            %plot Plots random points to test isInside
            
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
           if (~isa(obj,'Domain'))
               error("Value must be a Domain");
           end
       end     

       
       function obj = build(varargin)
            obj = Domain();
            if (nargin == 0), return;

            elseif (nargin == 1)
                varargin = varargin{1};
                if (isa(varargin,'struct'))
                    if (isfield(varargin, 'type') && isfield(varargin,'args')) % check if this is a parsed struct
                        if (~strcmp(varargin.type, 'default'))
                            name = lower(varargin.type) + "Domain"; 
                            if (exist(name,'class') == 8)
                                argNames = fieldnames(varargin.args);
                                celin = cell(1,length(argNames) * 2+1);
                                celin{1} = name;
                                for i = 1:(numel(argNames))
                                    celin{2*i} = argNames{i};
                                    celin{2*i+1} = varargin.args.(argNames{i});
                                end
                                obj = feval(celin{:});

                            end
                        
                        end
                    else    
                        obj.bdr   = varargin.bdr;
                        obj.dbdr  = varargin.dbdr;
                        obj.ddbdr = varargin.ddbdr;
                    end
                    obj.theta = varargin.theta;
                    obj.originX = varargin.originX;
                    obj.originY = varargin.originY;
                    return;
                end
            elseif (nargin > 1) 
                
            end    


            p = inputParser; % it parses inputs
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            p.FunctionName = 'Domain.build';

            isAString = @(in) isstring(in) || ischar(in);
            isANumber = @(in) isnumeric(in);
            isA2Vector = @(in) isnumeric(in) && all(size(in) == [1,2]);
            isA1Handle = @(in) isa(in,'function_handle') && nargin(in) == 1;
           addOptional(p,'type','default',isAString);

            addParameter(p,'theta',obj.theta,isANumber);
            addParameter(p,'origin',[obj.originX,obj.originY],isA2Vector);

            addParameter(p,'bdr',@(t) obj.bdr(t),isA1Handle);
            addParameter(p,'dbdr',@(t) obj.dbdr(t),isA1Handle);
            addParameter(p,'ddbdr',@(t) obj.ddbdr(t),isA1Handle);

            parse(p,varargin{:})
            r = p.Results;
            stin.type = r.type;
            stin.theta = r.theta;
            stin.originX = r.origin(1);
            stin.originY = r.origin(2);
            stin.args = p.Unmatched;

            obj = Domain.build(stin);
        end    
        
       
        
    end
    
end