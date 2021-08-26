classdef constcurveMetric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kappa
    end
    
    methods
        function obj = constcurveMetric(args)
            arguments
                args.kappa (1,1) {mustBeNumeric} = 2
            end
            
            obj.kappa = args.kappa;
        end
        
        function out = lg(obj,x,y)
            out = -2*log(1 + obj.kappa*(x.*x + y.*y));
        end
        
        function out = dxlg(obj,x,y)
            k = obj.kappa;
            out = -4*k*x./(1+k*(x.*x + y.*y));
        end
        
        function out = dylg(obj,x,y)
            k = obj.kappa;
            out = -4*k*y./(1+k*(x.*x + y.*y));        
        end
        
        function out = curv(obj,~,~)
            out = 4*obj.kappa;
        end
       
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            k = obj.kappa;
            
            dxlgt = -4*k*X./(1+k*(X.*X + Y.*Y));
            dylgt = -4*k*Y./(1+k*(X.*X + Y.*Y));
            lgt = -2*log(1 + k*(X.*X + Y.*Y));
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            k = obj.kappa;
            
            dxlgt = -4*k*X./(1+k*(X.*X + Y.*Y));
            dylgt = -4*k*Y./(1+k*(X.*X + Y.*Y));
            lgt = -2*log(1 + k*(X.*X + Y.*Y));
            curvt = 4*k*ones(size(X));
        end
        
        %{
        function [xO,yO,zO] = deproject(obj,X,Y)
            k = 1/(4*obj.kappa); % = R^2
            if (k < 0); error('Constant-negative surfaces cannot be represented in R3.'); end
            
            s = 1./(k+X.*X+Y.*Y);
            xO = 2*k*X.*s;
            yO = 2*k*Y.*s;
            zO = sqrt(k)*(k-X.*X-Y.*Y).*s;
        end 
        %}
    end
end

