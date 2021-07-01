classdef Metric
    %METRIC Parent class of other metrics
    %   - Default implementation of a metric.

    properties (Access = 'private')
        lgVAR = @(x,y) zeros(size(x));
    end    
    
    methods  
        
        function out = lg(obj,x,y)
            out = obj.lgVAR(x,y);
        end
        
        function out = dxlg(obj,x,y)
            out = deriv(@(x0) obj.lg(x0,y), x);
        end
        
        function out = dylg(obj,x,y)
            out = deriv(@(y0) obj.lg(x,y0), y);
        end
        
        function out = curv(obj,x,y)
            out = -(dderiv(@(x0) obj.lg(x0,y), x) + ...
                    dderiv(@(y0) obj.lg(x,y0), y)) .* ...
                   (0.5*exp(-obj.lg(x,y)));
        end
        
        
        
        function obj = Metric(lg)
        if (strcmp(class(obj), 'Metric')) 
            if (nargin == 1 && isa(lg,'function_handle') && nargin(lg) == 2)
                obj.lgVAR = lg;
            end
              
        end,end 
               
    
    
        function [lg,dxlg,dylg,curv] = getHandles(obj)
            lg = @(x,y) obj.lg(x,y);
            dxlg = @(x,y) obj.dxlg(x,y);
            dylg = @(x,y) obj.dylg(x,y);
            curv = @(x,y) obj.curv(x,y);
        end    
    
        
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            %metricVals Precompute values over pairs of (x,y)
            %   
            
            if (nargin ~= 3)
                error("Incorrect input arguments.")
            elseif (~(ismatrix(X) && ismatrix(Y)))
                error("X,Y must be vectors.")
            elseif (size(X) ~= size(Y))
                error("The sizes of X,Y do not agree.")
            end

            %warning("Method is slow, consider implementing a faster override of metricVals or use a subclass of Metric.")
            
            [nr,nc] = size(X);

            xtmp = reshape(X,nr*nc,1);
            ytmp = reshape(Y,nr*nc,1);
        
            lgt =   reshape(obj.lg(xtmp,ytmp)  , nr,nc);
            dxlgt = reshape(obj.dxlg(xtmp,ytmp), nr,nc);
            dylgt = reshape(obj.dylg(xtmp,ytmp), nr,nc);
        end
        
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            %metricValsCurv Precompute values
            %   
            
            if (nargin ~= 3)
                error("Incorrect input arguments.")
            elseif (~(ismatrix(X) && ismatrix(Y)))
                error("X,Y must be vectors.")
            elseif (size(X) ~= size(Y))
                error("The sizes of X,Y do not agree.")
            end

            %warning("Method is slow, consider implementing a faster override of metricVals or use a subclass of Metric.")
            
            [nr,nc] = size(X);

            xtmp = reshape(X,nr*nc,1);
            ytmp = reshape(Y,nr*nc,1);
            
            [lgt,dxlgt,dylgt] = metricVals(obj, X, Y);        
            curvt = reshape(obj.curv(xtmp,ytmp), nr,nc);
        end
        
        
        
        function out = plot(obj,X,Y) % -- !! TODO: make nicer to work more like plot, replace surf or smthn !!
            %plot Plots gl
            %   !! Constructor kept calling "display", so I had to rename this !!
            %  
            
            [X0,Y0] = meshgrid(-5:0.2:5);
            if (nargin == 3 && (isvector(X) && isvector(Y))),
                [X0,Y0] = meshgrid(X,Y);
            elseif (nargin ~= 1), warning("Incorrect input arguments, plotting over generic domain.");
            end
            
            Z = exp(metricVals(obj,X0,Y0)*0.5); 
            m = mean(Z(1:5:end));
            s = sqrt(std(Z(1:5:end))); % TODO should implement percentiles or something
            
            Z = clamp(m-s,Z,m+s); % clamp output % !! TODO make a parameter for this !!
            out = contour(X0,Y0,Z,20);
        end
        
        
        
        function plotALL(obj,X,Y)
            %plot Plots gl, dxgl,dygl, curv
            %   - intended for debuging
            %  
            
            [X0,Y0] = meshgrid(-5:0.2:5);
            if (nargin == 3 && (isvector(X) && isvector(Y))),
                [X0,Y0] = meshgrid(X,Y);
            elseif (nargin ~= 1), warning("Incorrect input arguments, plotting over generic domain.");
            end
            
            [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj,X0,Y0);
            
            %%{
            m = -10; M = 10;
            lgt = clamp(m,lgt,M); % clamp output
            dxlgt = clamp(m,dxlgt,M); % !! TODO make a parameter for this !!
            dylgt = clamp(m,dylgt,M);
            curvt = clamp(m,curvt,M);
            %%}
            
            figure;
            subplot(2,2,1);
            p = surf(X0,Y0,lgt,'EdgeColor','none');
            view(2); 
            subplot(2,2,2);
            p = surf(X0,Y0,dxlgt,'EdgeColor','none');
            view(2); 
            subplot(2,2,3);
            p = surf(X0,Y0,dylgt,'EdgeColor','none');
            view(2); 
            subplot(2,2,4);
            p = surf(X0,Y0,curvt,'EdgeColor','none');
            view(2); 
        end
        
        
        
    end    
    
    methods (Static)
       function mustBeMetric(obj)           
           if (~isa(obj,'Metric'))
               error("Value must be a Metric");
           end
       end    
       
       
       
       function obj = build(varargin)
           obj = Metric();
            if (nargin == 0), return; 
            
            elseif (nargin == 1)
                varargin = varargin{1};
                if (isa(varargin,'struct'))
                    if (isfield(varargin, 'type') && isfield(varargin,'args')) % check if this is a parsed struct
                        if (~strcmp(varargin.type, 'default'))
                            name = lower(varargin.type) + "Metric"; 
                            if (exist(name,'class') == 8)

                                argNames = fieldnames(varargin.args);
                                celin = cell(1,length(argNames) * 2+1);
                                celin{1} = name;
                                for i = 1:numel(argNames)
                                    celin{i+1} = argNames{i};
                                    celin{i+2} = varargin.args.(argNames{i});
                                end

                                obj = feval(celin{:});

                            end
                        else    
                            obj.lg   = varargin.lg;
                            obj.dxlg = varargin.dxlg;
                            obj.dylg = varargin.dylg;
                            obj.curv = varargin.curv;
                        end
                        return;
                    end

                end
            elseif (nargin > 1) 
                  
            end
            p = inputParser; % it parses inputs
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            p.FunctionName = 'Metric.build';

            isAString = @(in) isstring(in) || ischar(in);
            isA2Handle = @(in) isa(in,'function_handle') && nargin(in) == 2;

            addOptional(p,'type','default',isAString);

            if (isa(varargin,'cell'))
                parse(p,varargin{:})
            else
                parse(p,varargin)
            end    
            r = p.Results;
            stin.type = r.type;
            stin.args = p.Unmatched;

            obj = Metric.build(stin);
        end    
    end  
    
    
end



