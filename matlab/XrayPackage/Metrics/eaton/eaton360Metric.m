classdef eaton360Metric < Metric
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
            nt = ((sqrt(r.^(-2)+1/27)+1./r).^(1/3) - (sqrt(r.^(-2)+1/27)-1./r).^(1/3)).^2;
            fnt = 3*nt.^(1/2)+nt.^(-1/2);
            dnt = -4./r.^2./fnt;

            lgt(bool)   = real(2*log(nt));
            dxlgt(bool) = real(2.*X(bool)./(nt.*r).*dnt);
            dylgt(bool) = real(2.*Y(bool)./(nt.*r).*dnt);            
        end
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            lgt = zeros(size(X));
            dxlgt = lgt;
            dylgt = lgt;
            
            r = sqrt(X.*X+Y.*Y);
            
            bool = r < 1;
            r = r(bool);
            nt = ((sqrt(r.^(-2)+1/27)+1./r).^(1/3) - (sqrt(r.^(-2)+1/27)-1./r).^(1/3)).^2;
            fnt = 3*nt.^(1/2)+nt.^(-1/2);
            dnt = -4./r.^2./fnt;
            ddnt = (8./r.^3 - .5*(3*nt.^(-1/2)-nt.^(-3/2)).*dnt.^2)./fnt;

            lgt(bool)   = real(2*log(nt));
            dxlgt(bool) = real(2.*X(bool)./(nt.*r).*dnt);
            dylgt(bool) = real(2.*Y(bool)./(nt.*r).*dnt);
            curvt(bool) = real(-1./nt.^2 .* (ddnt./nt + dnt./(r.*nt) - (dnt/nt).^2));
                        
        end
        
        
    end
end

