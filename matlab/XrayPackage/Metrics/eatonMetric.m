classdef eatonMetric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
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
            srxy2u2 = 2./ sqrt( max(min(X.*X+Y.*Y, zeros(size(X))+4 ),zeros(size(X))) );
            
            lgt = log(srxy2u2 - 1);
            dxlgt = -X.*srxy2u2.^3./(4*srxy2u2 - 4);
            dylgt = -Y.*srxy2u2.^3./(4*srxy2u2 - 4);
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            
            srxy2 = sqrt( max(min(X.*X+Y.*Y, zeros(size(X))+4 ),zeros(size(X))) );
            
            lgt = log(2./srxy2-1);
            dxlgt = -2.*X./(srxy2.^3.*(2./srxy2-1));
            dylgt = -2.*Y./(srxy2.^3.*(2./srxy2-1));

            curvt = 1./(2-srxy2).^3;
                        
        end
        
        
    end
end

