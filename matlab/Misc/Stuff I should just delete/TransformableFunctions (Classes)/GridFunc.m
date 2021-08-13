classdef GridFunc < TFunc
     
    properties
        interpType 
    end
    
    methods
        function obj = GridFunc(values, args)  
            %INMAP Constructs an instance of GridFunc
            arguments
                values = @(x,y) zeros(size(x));   
                args.offset (1,2) {mustBeNumeric} = [0,0]
                args.scale (1,2) {mustBeNumeric} = [1,1]
                args.rotation (1,1) {mustBeNumeric} = 0
                args.interpType (1,:) {mustBeText} = 'linear'
            end
            obj@TFunc(values, offset=args.offset, scale=args.scale, rotation=args.rotation); 
            
            obj.interpType = args.interpType;
        end
    end
    
    
    methods (Access = 'protected')
        
        function out = evaluate(obj, X,Y)
            %EVAL Evaluates the function at an array of points (X,Y).
            
            tmat = obj.tmatrix;

            x = X*tmat(1) + Y*tmat(2) + tmat(3); % matrix multiplication (when handled correctly, only these 6 indicies are necessarily unknown, and the resulting vector will always have a 1 in the third place)
            y = X*tmat(4) + Y*tmat(5) + tmat(6);
            % z = ones(size(X)) is what would be written here if we cared abt this index 

            v = obj.valuesA;
            s = size(v);
            o = (s+1)*0.5;
            switch obj.interpType
                case 'nearest'
                    out = nearest2(v, min(max(1, x + o(1)), s(2)),...
                                      min(max(1, y + o(2)), s(1)));   
                case 'linear'
                    out = lininterp2(v, min(max(1, x + o(1)), s(2)),...
                                        min(max(1, y + o(2)), s(1)));
                otherwise
                    error('Wrong interpolation method')
            end
            
        end
        
    end    
       
end    


% grid re-sample functions
%   x taken from [1,length(vals)], y from [1,height(vals)]

function out = nearest2(vals,x,y)
    s = size(vals);
    out = vals(sub2ind(s,round(y),round(x)));
end

function out = lininterp2(vals, x,y) 
    s = size(vals);
    vfy = vals( sub2ind(s,ceil(y),floor(x))) .* (y-floor(y))  +  vals(sub2ind(s,floor(y),floor(x))) .* (1-y+floor(y));
    out = vals( sub2ind(s,ceil(y), ceil(x))) .* (y-floor(y))  +  vals(sub2ind(s,floor(y), ceil(x))) .* (1-y+floor(y));
    out = (out-vfy) .* (x-floor(x)) + vfy;
end