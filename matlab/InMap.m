classdef InMap
    %INFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = 'protected')
        values = [1];
        type = 'i';
        
        tmatrix = eye(3);
    end
    
    properties (Dependent)
        eval 
    end
    
    
    methods
        function obj = InMap(values)  
            %INMAP Constructs an instance of InMap
            
            if (nargin == 1)
                obj.values = values; 

                if isa(values,'function_handle')
                    obj.type = 'h';
                elseif isnumeric(values)
                    obj.type = 'i';
                else
                    error('Wrong data type for values.');
                end 
            end
        end
        
        function obj = transform(obj, xoff,yoff, xscale,yscale, rot)
            %TRANSFORM Takes arguments to rotate, scale and offset the
            %function.
            arguments
                obj
                xoff (1,1) {mustBeNumeric} = 0
                yoff (1,1) {mustBeNumeric} = 0
                xscale (1,1) {mustBeNumeric} = 1
                yscale (1,1) {mustBeNumeric} = 1
                rot (1,1) {mustBeNumeric} = 0
            end
            
            obj.tmatrix = inv([cos(rot)*xscale, sin(rot)*xscale,0;...
                              -sin(rot)*yscale, cos(rot)*yscale,0;...
                               xoff,            yoff,           1]);
        end
        
        function obj = bake(obj, xres,yres, xmin,ymin,xmax,ymax)
            %BAKE Converts the type of InMap from a function handle to an
            %image.
            
            arguments
                obj
                xres (1,1) {mustBeNumeric} = 0.2
                yres (1,1) {mustBeNumeric} = 0.2
                xmin (1,1) {mustBeNumeric} = -2
                ymin (1,1) {mustBeNumeric} = -2
                xmax (1,1) {mustBeNumeric} = 2
                ymax (1,1) {mustBeNumeric} = 2                
            end
            
            [X,Y] = meshgrid(xmin:1/xres:xmax,ymin:1/yres:ymax);
            obj.values = obj.values(obj.eval(X,Y));
            obj.type = 'i';
            obj.tmatrix = inv([xscale, 0,0;...
                               0, yscale,0;...
                               0.5*(xmin+xmax),0.5*(ymin+ymax),1]);
        end         
        
        function out = plot(obj, lbound,ubound)
            %PLOT Plots InMap on the axis aligned box defined by a lower
            %and upper bound.
            %   Uses obj.eval.
            
            arguments
                obj
                lbound (1,2) {mustBeNumeric} = [-2,-2]
                ubound (1,2) {mustBeNumeric} = [2,2]
            end
            
            [X,Y] = meshgrid(lbound(1):0.01:ubound(1),...
                             lbound(2):0.01:ubound(2));
                         
            Z = obj.eval(X,Y);
            hold on;
            out = pcolor(X,Y,Z);
            regT = inv(obj.tmatrix);
            plot(regT(3),regT(6),'r*');
            out.EdgeColor = 'none';
            axis square;
            xlim([-2,2]); ylim([-2,2]);
        end
        
        function value = get.eval(obj)
           value = @(X,Y) obj.evaluate(X,Y);
        end
        
    end
    
    methods (Access = 'protected')
        
        function out = evaluate(obj, X,Y)
            %EVAL Evaluates the function at an array of points (X,Y).
            
            tmat = obj.tmatrix;

            x = X*tmat(1) + Y*tmat(2) + tmat(3); % matrix multiplication (when handled correctly, only these 6 indicies are necessarily unknown, and the resulting vector will always have a 1 in the third place)
            y = X*tmat(4) + Y*tmat(5) + tmat(6);
            % z = ones(size(X)) is what would be written here if we cared abt this index 
            switch obj.type
                case 'h'
                    out = obj.values(x,y);
                case 'i'
                    v = obj.values;
                    s = size(v);
                    o = (s+1)*0.5;
                    out = lininterp2(v, clamp(1, x + o(1), s(2)),... % maybe include other edge cases using mod/force zero
                                        clamp(1, y + o(2), s(1)));
            end
            
        end
        
    end    
end


% image re-sample functions
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
