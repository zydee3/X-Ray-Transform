classdef polynomialMetric < Metric

    properties
        coeffs (1,6) {mustBeNumeric}
    end    
    
    methods
        function obj = polynomialMetric(args)
            arguments
                args.coeffs (1,6) {mustBeNumeric} = [0,0,0,0,0,1]
            end
            
            obj.coeffs = args.coeffs;
        end
        
        function out = lg(obj,x,y)
            c = obj.coeffs;
            out = x.*(c(1)*x + c(2)*y + c(4)) ...
                  + y.*(c(3)*y + c(5)) + c(6);
        end
        
        function out = dxlg(obj,x,y)
            c = obj.coeffs;
            out = 2*c(1)*x + c(2)*y + c(4);
        end
        
        function out = dylg(obj,x,y)
            c = obj.coeffs;
            out = c(2)*x + 2*c(3)*y + c(5);
        end
        
        function out = curv(obj,x,y)
            c = obj.coeffs;
            out = -exp(-obj.lg(x,y)).*(c(1) + c(3));
        end
        
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            c = obj.coeffs;
            lgt = X.*(c(1)*X + c(2)*Y + c(4)) ...
                  + Y.*(c(3)*Y + c(5)) + c(6);
            dxlgt = 2*c(1)*X + c(2)*Y + c(4);
            dylgt = c(2)*X + 2*c(3)*Y + c(5);
        end
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            c = obj.coeffs;
            lgt = X.*(c(1)*X + c(2)*Y + c(4)) ...
                  + Y.*(c(3)*Y + c(5)) + c(6);
            dxlgt = 2*c(1)*X + c(2)*Y + c(4);
            dylgt = c(2)*X + 2*c(3)*Y + c(5);
            curvt = -exp(-lgt).*(c(1) + c(3));
        end
        

    end
end


