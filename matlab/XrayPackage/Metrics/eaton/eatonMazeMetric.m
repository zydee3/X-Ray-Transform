classdef eatonMazeMetric < Metric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        randomItter = 7
    end
    
    properties
        seed = 12345678
        metrics = {euclidMetric,...euclidMetric,...
                   eaton90Metric,...
                   eaton180Metric,...
                   ...eatonGeneralMetric(theta = 3*pi/2),...                   
                   eaton360Metric...
                   }
    end
    
    methods
        
        function obj = eatonMazeMetric(args)
            arguments
                args.seed (1,1) {mustBeNumeric} = 12345678
            end
            
            obj.seed = args.seed;
        end
        
        function [lgt,dxlgt,dylgt] = metricVals(obj, X, Y)
            cs = 1;%obj.cellRadius;
            rx = round(X/(2*cs)); ry = round(Y/(2*cs));
            
            % get type index
            mets = obj.metrics;
            numType = length(mets);
            typeInd = floor(prandom(obj.seed,rx,ry)*numType)+1; %mod(rx+ry,numType)+1;%
            % modularize position
            xin = X - rx*2;
            yin = Y - ry*2;
            % remove weird values, initialize lgt,dxlgt,dylgt
            lgt = zeros(size(X));
            bool = (xin.*xin + yin.*yin) < 0.000001;
            typeInd(bool) = xin(bool);
            lgt(bool) = NaN(1,nnz(bool));
            dxlgt = lgt;
            dylgt = lgt;
            % access metricvals
            for (ti = 1:numType)
                bool = (ti == typeInd);
                [lgt(bool),dxlgt(bool),dylgt(bool)] = mets{ti}.metricVals(xin(bool),yin(bool));
            end
            
            %[testx,testy] = meshgrid(-10:10);
            %pcolor(prandom(obj.seed,testx,testy))
        end
                
        
        
    end
    
end




function out = prandom(seed, X,Y)
    w = 131071;
    X = X+w; Y = Y+w;
    for i = 1:eatonMazeMetric.randomItter
        seed = mod(w + seed - 524281,1048576);
        Y=mod(X.*Y + seed,65536);
        X=mod(X.*Y + seed,65536); % lotta arbitrarily chosen numbers you got here
    end
    out = (X+Y) /(65536+65536);
end

