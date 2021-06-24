classdef euclidMetric < Metric

    methods
        function out = lg(obj,x,y)
            out = zeros(size(x));
        end
        
        function out = dxlg(obj,x,y)
            out = zeros(size(x));
        end
        
        function out = dylg(obj,x,y)
            out = zeros(size(x));
        end
        
        function out = curv(obj,x,y)
            out = zeros(size(x));
        end
        
         
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            s = zeros(size(X));
            lgt   = s;
            dxlgt = s;
            dylgt = s;
        end

    end
end

