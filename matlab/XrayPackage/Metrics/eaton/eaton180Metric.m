classdef eaton180Metric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function out = lg(obj,x,y)
            out = log(2./sqrt(x.*x+y.*y)-1);
        end
        
        function out = dxlg(obj,x,y)
            out = -2./(2./sqrt(x.*x+y.*y)-1).*x.*(x.*x+y.*y).^(-1.5);
        end
        
        function out = dylg(obj,x,y)
            out = -2./(2./sqrt(x.*x+y.*y)-1).*y.*(x.*x+y.*y).^(-1.5);
        end
        
        function out = curv(obj,x,y)
            out = 1./(2-sqrt(x.*x+y.*y)).^3;
        end
       
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            
            r = sqrt(X.*X+Y.*Y);
            
            bool = r < 1;
            r = r(bool);

            lgt(bool) = log(2./r-1);
            dxlgt(bool) = -2*exp(-lgt(bool)).*X(bool)./r.^3;
            dylgt(bool) = -2*exp(-lgt(bool)).*Y(bool)./r.^3;
        end
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            
            r = sqrt(X.*X+Y.*Y);
            
            bool = r < 1;
            r = r(bool);

            lgt(bool) = real(log(2./r-1));
            dxlgt(bool) = real(-2*exp(-lgt(bool)).*X(bool)./r.^3);
            dylgt(bool) = real(-2*exp(-lgt(bool)).*Y(bool)./r.^3);
            curvt(bool) = real(1./(2-r).^3);
                        
        end        
        
    end
end

