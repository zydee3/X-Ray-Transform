classdef EXAMPLECODE_myCustomMetric0 < Metric
    % To define a new metric, define a new class that inherits from the Metric class.
    % This is enough for our metric to be functional, although it will
    % default to the euclidean metric if we don't override any methods.
    
    methods
        
        % The simplest meaningful implementation of a metric additionally requiers a definition of the metricVals function. 
        % The remaining functions are then *autocompleted* by the Metric class unless also overwritten.
        % Implementation of these methods (especially metricVals and metricValsCurv) should be vectorized for many inputs.
        
        % Of course, implementing definitions that contradict one another or that contradict the super class implementations leads to unpredictable and unwanted behaviour.
        
        % Ideally, a metric is defined with all of metricVals and metricValsCurv and perhaps lg,dxlg,dylg, and curv.
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)            
            lgt = X-Y;
            dxlgt = 1;
            dylgt = -1;
        end
        
        
    end
    
end

