classdef gaussMetric < Metric

    properties
        weights
        widths
        xOffsets
        yOffsets
    end
    
    methods
        function obj = constcurveMetric(args)
            arguments
                args.weights (1,:) {mustBeNumeric} = 1
                args.widths (1,:) {mustBeNumeric} = 1
                args.xOffsets (1,:) {mustBeNumeric} = 0
                args.yOffsets (1,:) {mustBeNumeric} = 0
            end
            
            obj.weights = args.weights;
            obj.widths = args.widths;
            obj.xOffsets = args.xOffsets;
            obj.yOffsets = args.yOffsets;
        end
        
        function out = lg(obj,x,y)
            xC = obj.xOffsets(:)';
            yC = obj.yOffsets(:)';
            k = obj.weights(:)';
            
            nb = length(sigma); % number of bumps
            [nr,nc] = size(x); % number of points is nr*nc
            np = nr*nc; % number of points

            sigs = ones(nr*nc,1) * obj.widths(:)';

            x = (x(:)*ones(1,nb) - ones(np,1)*xC)./sigs; % reshape transform X and Y
            y = (y(:)*ones(1,nb) - ones(np,1)*yC)./sigs;
            
            out = (ones(nr*nc,1)* k) .* exp(-0.5 * (x.*x + y.*y)); % compute stuffs, reshape back
            out = reshape( sum(out, 2), nr, nc);
        end
        
        %{
        
        function out = dxlg(obj,x,y)
            % TODO: implement 
            out = 0;
        end
        
        function out = dylg(obj,x,y)
            % TODO: implement 
            out = 0;      
        end
        
        function out = curv(obj,~,~)
            out = 4*obj.kappa;
        end
        %}
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            xC = obj.xOffsets;
            yC = obj.yOffsets;
            k = obj.weights;
            
            nb = length(sigma); % number of bumps
            [nr,nc] = size(X); % number of points is nr*nc
            np = nr*nc; % number of points

            sigs = ones(nr*nc,1) * obj.widths;

            X = X(:); % reshape X and Y
            Y = Y(:);
            X = (X*ones(1,nb) - ones(np,1)*xC)./sigs; % transform X and Y
            Y = (Y*ones(1,nb) - ones(np,1)*yC)./sigs;
            
            tmp = (ones(nr*nc,1)* k) .* exp(-0.5 * (X.*X + Y.*Y)); % compute stuffs, reshape back
            dxlgt = reshape( sum( -X./sigs .* tmp, 2 ), nr, nc);
            dylgt = reshape( sum( -Y./sigs .* tmp, 2 ), nr, nc);
            lgt = reshape( sum(tmp, 2), nr, nc);
        end
        
        
        function [lgt,dxlgt,dylgt,curvt] = metricValsCurv(obj, X, Y)
            xC = obj.xOffsets;
            yC = obj.yOffsets;
            k = obj.weights;
            
            nb = length(sigma); % number of bumps
            [nr,nc] = size(X); % number of points is nr*nc
            np = nr*nc; % number of points

            sigs = ones(nr*nc,1) * obj.widths;

            X = X(:); % reshape X and Y
            Y = Y(:);
            X = (X*ones(1,nb) - ones(np,1)*xC)./sigs; % transform X and Y
            Y = (Y*ones(1,nb) - ones(np,1)*yC)./sigs;
            
            tmp = (ones(nr*nc,1)* k) .* exp(-0.5 * (X.*X + Y.*Y)); % compute stuffs, reshape back
            dxlgt = reshape( sum( -X./sigs .* tmp, 2 ), nr, nc);
            dylgt = reshape( sum( -Y./sigs .* tmp, 2 ), nr, nc);
            lgt = reshape( sum(tmp, 2), nr, nc);

            curvt = -.5 * exp(-lgt) .* reshape( sum( (X.^2 + Y.^2 - 2)./sigs.^2 .* tmp, 2 ), nr, nc);
        end
        
    end    
end

