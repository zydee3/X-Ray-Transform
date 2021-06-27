classdef myexample1Metric < Metric
    methods
        function out = lg(obj,x,y)
            %out = 5 * exp(-1./(1-min(abs(x),1).^2));
            out = -2 + 2 * exp(-1./(1-min(max(x*5,0),1).^2));
        end
    end
end
