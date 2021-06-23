classdef EuclidMetric < Metric

    methods
        function obj = EuclidMetric()
            obj.lg   = @(x,y) 0.;
            obj.dxlg = @(x,y) 0.;
            obj.dylg = @(x,y) 0.;
            obj.curv = @(x,y) 0.;   
        end
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            s = size(X);
            lgt   = zeros(s);
            dxlgt = zeros(s);
            dylgt = zeros(s);
        end

    end
end

