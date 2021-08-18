classdef gutterMetric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sig
        scale
        yOffset
    end
    
    methods
        function obj = gutterMetric(args)
            arguments
                args.sig (1,1) {mustBeNumeric} = 1
                args.scale (1,1) {mustBeNumeric} = 1
                args.yOffset (1,1) {mustBeNumeric} = 0
            end
            
            obj.sig = args.sig;
            obj.scale = args.scale;
            obj.yOffset = args.yOffset;
        end
        
        function out = lg(obj,x,y)
            out = obj.scale*exp(-(y-obj.yOffset).^2/(2*obj.sig^2));
        end
        
        function out = dxlg(obj,~,y)
            out = zeros(size(y)); 
        end
        
        function out = dylg(obj,x,y)
            out = -(y-obj.yOffset)/obj.sig^2 .* obj.lg(x,y);     
        end
        
        function out = curv(obj,~,~)
            s = obj.sig;
            out = -.5*exp(-obj.lg(x,y)) .* ((y-y0).^2 - s^2) .* obj.lg(x,y) /s^4;
        end
       
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            s = obj.sig;
            ybar = (Y - obj.yOffset) / s;
            k = obj.scale;
            
            lgt = k .* exp(- ybar.^2 /2 );
            dxlgt = zeros(size(X));
            dylgt = - (ybar./s) .* lgt;
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            s = obj.sig;
            ybar = (Y - obj.yOffset) / s;
            k = obj.scale;
            
            lgt = k .* exp(- ybar.^2 /2 );
            dxlgt = zeros(size(X));
            dylgt = - (ybar./s) .* lgt;

            curvt = - .5 * exp(-lgt) .* (ybar.^2 - 1)/s^2 .* lgt;
        end
        
        
    end
end

