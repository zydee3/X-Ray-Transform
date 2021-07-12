classdef InMap
    %INFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = 'protected')
        values = [1];
        type = 'i';
        
        tmatrix = eye(3);
    end
    
    
    methods
        function obj = InMap(values)  
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
            
            obj.tmatrix = inv([cos(rot)*xscale, sin(rot)*xscale,0;...
                              -sin(rot)*yscale, cos(rot)*yscale,0;...
                               xoff,            yoff,           1]);
        end
        
        function obj = bake(obj, xres,yres, xmin,ymin,xmax,ymax)
            
            [X,Y] = meshgrid(xmin:1/xres:xmax,ymin:1/yres:ymax);
            obj.values = obj.values(obj.eval(X,Y));
            obj.type = 'i';
            obj.tmatrix = inv([xscale, 0,0;...
                               0, yscale,0;...
                               0.5*(xmin+xmax),0.5*(ymin+ymax),1]);
        end 
        
        function out = eval(obj, X,Y)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            tmat = obj.tmatrix;

            x = X*tmat(1) + Y*tmat(2) + tmat(3); % matrix multiplication
            y = X*tmat(4) + Y*tmat(5) + tmat(6);
            
            switch obj.type
                case 'h'
                    out = obj.values(x,y);
                case 'i'
                    v = obj.values;
                    s = size(v);
                    o = (s+1)*0.5;
                    out = lininterp2(v, clamp(1, x + o(1), s(1)),... % maybe include other edge cases using mod/force zero
                                        clamp(1, y + o(2), s(2)));
            end
            
        end
           
        
        
        function out = plot(obj)
            [X,Y] = meshgrid(-6:0.1:6);
            Z = obj.eval(X,Y);
            hold on;
            out = pcolor(X,Y,Z);
            plot(obj.tmatrix(3),obj.tmatrix(6),'r*');
            out.EdgeColor = 'none';
            axis square;
            xlim([-6,6]); ylim([-6,6]);
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
    vcy = vals( sub2ind(s,ceil(y), ceil(x))) .* (y-floor(y))  +  vals(sub2ind(s,floor(y), ceil(x))) .* (1-y+floor(y));
    vfy = vals( sub2ind(s,ceil(y),floor(x))) .* (y-floor(y))  +  vals(sub2ind(s,floor(y),floor(x))) .* (1-y+floor(y));
    out = (vcy-vfy) .* (x-floor(x)) + vfy;
end
