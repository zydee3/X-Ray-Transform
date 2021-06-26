classdef myexample1Metric < Metric
    methods
        function out = lg(obj,x,y)
            out = square(x) + square(y);
        end
    end
end

function out = square(X)
    out = 2 * (sin(X) > 0) - 1;
end