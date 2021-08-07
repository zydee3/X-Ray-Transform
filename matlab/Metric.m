classdef Metric
    %METRIC Parent class of other metrics
    %   - Default implementation of a metric.

    properties (Access = 'private')
        lgVAR = @(x,y) zeros(size(x));
    end    
    
    methods  

%--------------------------------------------------------------------------
%%                 Constructor, Characteristic Functions
%--------------------------------------------------------------------------          
        
        function out = lg(obj,x,y)
            %LG Computes the log of g.
            %   Method uses the function handle obj.lgVAR if not 
            %   overwritten by a subclass. 
            out = obj.lgVAR(x,y);
        end
        
        function out = dxlg(obj,x,y)
            %DXLG Derivative of obj.bdr with respect to x.
            %   Method computes the derivative numerically from obj.lg if not
            %   overwritten by a subclass.
            out = deriv(@(x0) obj.lg(x0,y), x);
        end
        
        function out = dylg(obj,x,y)
            %DYLG Derivative of obj.bdr with respect to y.
            %   Method computes the derivative numerically from obj.lg if not
            %   overwritten by a subclass.
            out = deriv(@(y0) obj.lg(x,y0), y);
        end
        
        function out = curv(obj,x,y)
            %CURV Curvature described by the metric at the point (x,y).
            %   Method computes the curvature numerically from obj.lg if not
            %   overwritten by a subclass.
            out = -(dderiv(@(x0) obj.lg(x0,y), x) + ...
                    dderiv(@(y0) obj.lg(x,y0), y)) .* ...
                   (0.5*exp(-obj.lg(x,y)));
        end
        
               
        function obj = Metric(lg)
            %METRIC Constructs an instance of Metric.
        if (strcmp(class(obj), 'Metric')) 
            if (nargin == 1 && isa(lg,'function_handle') && nargin(lg) == 2)
                obj.lgVAR = lg;
            end
              
        end,end 
                  
        
%--------------------------------------------------------------------------
%%                                  Misc
%--------------------------------------------------------------------------  

        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            %METRICVALS Precomputes values of lg,dxlg,dylg over an array of
            %points (X,Y).
            %   Method is intended to be a vectorized, faster computation
            %   of obj.lg, obj.dxlg, and obj.dylg
            %   Method just uses obj.lg, obj.dxlg, and obj.dylg if not
            %   overwritten by a subclass.
            
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
            %METRICVALSCURV Precomputes values of lg,dxlg,dylg and curv over an array of
            %points (X,Y).
            %   Method is intended to be a vectorized, faster computation
            %   of obj.lg, obj.dxlg, obj.dylg and obj.curv.
            %   Method just uses obj.lg, obj.dxlg, obj.dylg and obj.curv if not
            %   overwritten by a subclass.
            %   See obj.metricVals
            
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
        

        function [lg,dxlg,dylg,curv] = getHandles(obj)
            %METRICVALSCURV Returns lg, dxlg, dylg and curv as function
            %handles.
            lg = @(x,y) obj.lg(x,y);
            dxlg = @(x,y) obj.dxlg(x,y);
            dylg = @(x,y) obj.dylg(x,y);
            curv = @(x,y) obj.curv(x,y);
        end    
    
                               
%--------------------------------------------------------------------------
%%                                Plotters
%--------------------------------------------------------------------------              

        function out = plot(obj, lbound,ubound, reso)
            %PLOT Plots a contour plot of exp(lg/2) given a lower and upper
            %bound.
            
            arguments
                obj
                lbound (1,2) {mustBeNumeric} = [-5,-5]
                ubound (1,2) {mustBeNumeric} = [5,5]
                reso (1,1){mustBeNumeric} = 0.2
            end
            
            [X0,Y0] = meshgrid(lbound(1):reso:ubound(1),...
                               lbound(2):reso:ubound(2));
            
            Z = exp(metricVals(obj,X0,Y0)*0.5); 
            m = mean(Z(1:5:end));
            s = sqrt(std(Z(1:5:end))); % TODO should implement percentiles or something
            
            Z = clamp(m-s,Z,m+s); % clamp output % !! TODO make a parameter for this !!
            out = contour(X0,Y0,Z,20);
        end
                
        
        function figureALL(obj,X,Y)
            %PLOTALL A convinient function that creates a new figure and 
            %plots lg, dxlg, dylg and curv. 
            %   Intended to test that everything is working as it should.
            
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
            title('lg') 
            subplot(2,2,2);
            p = surf(X0,Y0,dxlgt,'EdgeColor','none');
            view(2); 
            title('dxlg') 
            subplot(2,2,3);
            p = surf(X0,Y0,dylgt,'EdgeColor','none');
            view(2); 
            title('dylg') 
            subplot(2,2,4);
            p = surf(X0,Y0,curvt,'EdgeColor','none');
            view(2); 
            title('curv') 
        end
        
        
        
    end    
    
    methods (Static)
       function mustBeMetric(obj)     
           %MUSTBEMETRIC Errors if the passed object is not an instance of 
           %Metric or of a subclass of Metric.
           if (~isa(obj,'Metric'))
               error("Value must be a Metric");
           end
       end    
              
       
       function obj = build(varargin) %TODO: Delete this grabage.
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



