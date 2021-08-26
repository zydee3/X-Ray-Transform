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
                                       
            if (~strcmp(class(obj), 'Metric')) % determine if method has been overwritten
                out = obj.metricVals(x, y);
            else
                out = obj.lgVAR(x,y);
            end    
            
        end
               
        function out = dxlg(obj,x,y)
            [~,out] = obj.metricVals(x,y);
        end
        
        function out = dylg(obj,x,y)
            [~,~,out] = obj.metricVals(x,y);
        end
        
        function out = curv(obj,x,y)
            [~,~,~,out] = obj.metricValsCurv(x,y);
        end
        
        function obj = Metric(lg)
            %METRIC Constructs an instance of Metric.
        if (strcmp(class(obj), 'Metric')) 
            if (nargin == 1 && isa(lg,'function_handle') && nargin(lg) == 2)
                obj.lgVAR = lg;
            end
              
        end,end 

    
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
            
            [nr,nc] = size(X);

            xtmp = reshape(X,nr*nc,1);
            ytmp = reshape(Y,nr*nc,1);
        
            lgt =   reshape(obj.lgVAR(xtmp,ytmp)            , nr,nc);
            dxlgt = reshape(deriv(@(x0) obj.lgVAR(x0,Y), X) , nr,nc);
            dylgt = reshape(deriv(@(y0) obj.lgVAR(X,y0), Y) , nr,nc);
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
            
            [nr,nc] = size(X);
            
            [lgt,dxlgt,dylgt] = metricVals(obj, X, Y);     
            curvt = reshape(-(dderiv(@(x0) obj.metricVals(x0,Y), X) + ...
                              dderiv(@(y0) obj.metricVals(X,y0), Y)) .* ...
                             (0.5*exp(-obj.metricVals(X,Y))),        nr,nc);
        end    
  
%--------------------------------------------------------------------------
%%                                  Misc
%--------------------------------------------------------------------------             
                          
        function [xO,yO,zO] = deproject(obj,X,Y)
            error('Metric must define a deproject method in order for deprojected plotters to function. Note that not all metrics can have such a method accurately defined.')            
        end
        
        function [xO,yO] = aabbspace(obj,domain,resoX,resoY)
            arguments
                obj
                domain (1,1) {Domain.mustBeDomain}
                resoX (1,1) {mustBeNumeric} = 250
                resoY (1,1) {mustBeNumeric} = resoX
            end
            
            [xO,yO] = domain.aabbspace(resoX,resoY);
        end

%--------------------------------------------------------------------------
%%                                Plotters
%--------------------------------------------------------------------------              

        function out = plot(obj, reso, minB,maxB)
            %PLOT Plots a contour plot of exp(lg/2) given a lower and upper
            %bound.
            
            arguments
                obj
                reso (1,1){mustBeNumeric} = 200
                minB (1,2) {mustBeNumeric} = [-5,-5]
                maxB (1,2) {mustBeNumeric} = -minB
            end
            
            [X0,Y0] = meshgrid(linspace(minB(1),maxB(1),reso),...
                               linspace(minB(2),maxB(2),reso));
            
            Z = exp(metricVals(obj,X0,Y0)*0.5); 
            out = contour(X0,Y0,Z,20);
        end
        
        function figureALL(obj, reso, lbound,ubound)
            %PLOTALL A convinient function that creates a new figure and 
            %plots lg, dxlg, dylg and curv. 
            %   Intended to test that everything is working as it should.
            
            arguments
                obj
                reso (1,1){mustBeNumeric} = 0.2
                lbound (1,2) {mustBeNumeric} = [-5,-5]
                ubound (1,2) {mustBeNumeric} = [5,5]
            end
            [X0,Y0] = meshgrid(lbound(1):reso:ubound(1),...
                               lbound(2):reso:ubound(2));
            
            [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj,X0,Y0);

            
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
        
%--------------------------------------------------------------------------              
                        
    end    
    
    methods (Static)
       function mustBeMetric(obj)     
           %MUSTBEMETRIC Errors if the passed object is not an instance of 
           %Metric or of a subclass of Metric.
           if (~isa(obj,'Metric'))
               error("Value must be a Metric.");
           end
       end    
             
    end  
    
    
end



