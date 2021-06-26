classdef myexample1Metric < Metric
    methods
        function out = lg(obj,x,y)
            out = x.^2 - y.^2;
        end
    end
end

