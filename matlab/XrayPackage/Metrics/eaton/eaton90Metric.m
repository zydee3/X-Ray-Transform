classdef eaton90Metric < Metric
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
            
            r2 = X.*X+Y.*Y;
            r = sqrt(r2);
            
            bool = r < 1;
            r2 = r2(bool); r = r(bool); X = X(bool); Y = Y(bool);
            
            sq3 = sqrt(3);
            cb2 = 2^(1/3);
            cb3 = 3^(1/3);
            
            pt = (9*r + sqrt(r2.*(81-48*r2.*r2))).^(1/3);
            qt = sqrt(cb2^5*r./(cb3*pt)+cb2*pt./(cb3*cb3*r));
            nt = qt+sqrt(-qt.*qt+4./(r.*qt));
            lgt(bool)   = 2*real(log(nt));
                        
            % dx
            dr = X./r;
            dpt = X.*(6*sq3+(54-96*r2.*r2)./sqrt(27-16*r2.*r2))./(2*sq3*r.*pt.*pt);
            dqt = (cb2^5./(pt.*pt)-cb2./(cb3*r2)).*(pt.*dr-r.*dpt)./(2*cb3*qt);
            dxlgt(bool) = 2*real( (dqt+(qt.*dqt+2*(dr.*qt+dqt.*r)./(qt.*qt.*r2))./(qt-nt))./nt );
            % dy
            dr = Y./r;
            dpt = Y.*(6*sq3+(54-96*r2.*r2)./sqrt(27-16*r2.*r2))./(2*sq3*r.*pt.*pt);
            dqt = (cb2^5./(pt.*pt)-cb2./(cb3*r2)).*(pt.*dr-r.*dpt)./(2*cb3*qt);
            dylgt(bool) = 2*real( (dqt+(qt.*dqt+2*(dr.*qt+dqt.*r)./(qt.*qt.*r2))./(qt-nt))./nt );
 
        end
              
    end
end

