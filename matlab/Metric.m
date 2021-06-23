classdef Metric
    %METRIC Parent class of other metrics
    %   - Default implementation of a metric.
    %   - Constructor defaults to euclidean, but accepts a function handle
    %   "g" so that the metric is conformal to the euclidean metric by 
    %   e^{log(g)}.
    %   - This class provides a display function and numerical precompute 
    %   methods. 
    
    properties        
        lg   = @(x,y) zeros(size(x));
        dxlg = @(x,y) zeros(size(x));
        dylg = @(x,y) zeros(size(x));
        curv = @(x,y) zeros(size(x));   
    end
    
    
    methods
        
        function obj = Metric(argA, argB)
        if (strcmp(class(obj), 'Metric')) %%% -- !! TODO: hack to avoid subclass calling superclass constructor, should prolly change !! 
            %METRIC Construct an instance of this class
            %   
                if (nargin == 0), return; end
            
                if (nargin == 1)   
                    if (isa(argA,'Metric') || isa(argA,'struct')) % Casts, we shouldn't need this.
                        obj.lg   = argA.lg;
                        obj.dxlg = argA.dxlg;
                        obj.dylg = argA.dylg;
                        obj.curv = argA.curv;    
                        return;
                        
                        
                        
                    elseif (isa(argA,'function_handle') && nargin(argA) == 2) % argA = g, computes the rest
                        warning("Metrics constructed by anonymous functions are very slow. Consider using or creating a subclass of metric.")
                        g = argA;
                        obj.lg   = @(x,y) log(g(x,y));
                        obj.dxlg = @(x,y) deriv(@(x0) log(g(x0,y)), x);
                        obj.dylg = @(x,y) deriv(@(y0) log(g(x,y0)), y);
                        obj.curv = @(x,y) -(dderiv(@(x0) log(g(x0,y)), x) + ...
                                            dderiv(@(y0) log(g(x,y0)), y)) ./ ...
                                           (2*g(x,y));  
                        return;
                     
                        
                    end    
                elseif (nargin == 2)                    % argA = type, argB = params
                    
                    return
                end
                
                
                error("Bad input arguments.")                    
        end,end 
               
    
    
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

            warning("Method is slow, consider implementing a faster override of metricVals or use a subclass of Metric.")
            
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

            warning("Method is slow, consider implementing a faster override of metricValsCurv or use a subclass of Metric.")
            
            if (nargin ~= 3)
                error("Incorrect input arguments.")
            elseif (~(ismatrix(X) && ismatrix(Y)))
                error("X,Y must be vectors.")
            elseif (size(X) ~= size(Y))
                error("The sizes of X,Y do not agree.")
            end

            warning("Method is slow, consider implementing a faster override of metricVals or use a subclass of Metric.")
            
            [nr,nc] = size(X);

            xtmp = reshape(X,nr*nc,1);
            ytmp = reshape(Y,nr*nc,1);
        
            lgt =   reshape(obj.lg(xtmp,ytmp)  , nr,nc);
            dxlgt = reshape(obj.dxlg(xtmp,ytmp), nr,nc);
            dylgt = reshape(obj.dylg(xtmp,ytmp), nr,nc);
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
            
            Z = metricVals(obj,X0,Y0); 
            Z = min(7,max(Z,-7)); % clamp output % !! TODO make a parameter for this !!
            out = pcolor(X0,Y0,Z);
            
            out.EdgeColor = 'none';
            out.FaceColor = 'interp';
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
            
            m = -7; M = 7;
            lgt = min(M,max(lgt,m)); % clamp output
            dxlgt = min(M,max(dxlgt,m)); % !! TODO make a parameter for this !!
            dylgt = min(M,max(dylgt,m));
            curvt = min(M,max(curvt,m));
            
            figure;
            subplot(2,2,1)
            p = pcolor(X0,Y0,lgt);
                p.EdgeColor = 'none';
                p.FaceColor = 'interp';
            subplot(2,2,2)
            p = pcolor(X0,Y0,dxlgt);
                p.EdgeColor = 'none';
                p.FaceColor = 'interp';
            subplot(2,2,3)
            p = pcolor(X0,Y0,dylgt);
                p.EdgeColor = 'none';
                p.FaceColor = 'interp';
            subplot(2,2,4)
            p = pcolor(X0,Y0,curvt);
                p.EdgeColor = 'none';
                p.FaceColor = 'interp';
        end
        
        
        
    end
end



