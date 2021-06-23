classdef PolyMetric < Metric

    properties
        coeffs = [0,0,0,0,0,1]
    end    
    
    methods
        function obj = PolyMetric(coeffs)
            if (nargin == 0), coeffs = obj.coeffs; 
            elseif(~isnumeric(coeffs) || ~all(size(coeffs)==[1,6]))
                error('Bad input arguments.');
            end % !! TODO: Better way of doing this? !!
            
            obj.coeffs = coeffs;
            
            obj.lg = @(x,y) coeffs(1)*x.^2 + coeffs(2)*x.*y + coeffs(3)*y.^2 + coeffs(4)*x + coeffs(5)*y + coeffs(6);
            obj.dxlg = @(x,y) 2*coeffs(1)*x + coeffs(2)*y + coeffs(4);
            obj.dylg = @(x,y) coeffs(2)*x + 2*coeffs(3)*y + coeffs(5);
            obj.curv = @(x,y) -exp(-obj.lg(x,y)).*(coeffs(1) + coeffs(3));
        end
        
        %{
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)

        end
        %}

    end
end


