classdef polyMetric < Metric

    properties
        coeffs (1,6) {mustBeNumeric}
    end    
    
    methods
        function obj = polyMetric(args)
            arguments
                args.coeffs (1,6) {mustBeNumeric} = [0,0,0,0,0,1]
            end
            
            obj.coeffs = args.coeffs;
        end
        
        function out = lg(obj,x,y)
            c = obj.coeffs;
            out = c(1)*x.^2 + c(2)*x.*y + c(3)*y.^2 + c(4)*x + c(5)*y + c(6);
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
        
        %{
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)

        end
        %}

    end
end


