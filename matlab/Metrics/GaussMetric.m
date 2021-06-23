classdef GaussMetric < Metric

    methods
        function obj = GaussMetric()
            obj.lg   = @(x,y) 0.;
            obj.dxlg = @(x,y) 0.;
            obj.dylg = @(x,y) 0.;
            obj.curv = @(x,y) 0.;   
        end
        
        %{
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)

        end
        %}

    end
end

