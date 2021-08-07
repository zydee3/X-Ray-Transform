classdef Map
    % class: Map
    % Author: Francois Monard - 10-27-2016
    
    properties
        a
        n % gridsize of definition. n=0 if function handles
        dom % domain
    end
    methods
        function obj=Map(A)
            % MAP: creates a connection, whose components are either
            % function handles or defined on a grid
            
            obj.a = A; 
            
            if isa(A,'function_handle')
                obj.n = 0;
            elseif isa(A,'double')
                obj.n = size(A,1);
            else
                error('Wrong data type for A.');
            end            
            
            %obj.dom = domain;
        end
        
        function A=val(obj,x,y)
            % MAP.VAL: builds the values of the attenuation 
            % at given gridpoints on SM
            
            if obj.n == 0, % regular function evaluation
                A = obj.a(x,y);
            else % bilinear interpolation of values at gridpoints
                dx = 10/obj.n;
                
                xrs = (x - 6 + 5-dx/2)/dx;
                yrs = (y - 6 + 5-dx/2)/dx;
                
                i0 = floor(yrs);
                j0 = floor(xrs);
                xbar = xrs - j0;
                ybar = yrs - i0;            
                i0 = max(min(i0 + 1,obj.n-1),1);
                j0 = max(min(j0 + 1,obj.n-1),1);
                
                A = ( obj.a( obj.n*(j0-1) + i0) .* (1.-xbar) .* (1.-ybar) ...
                    + obj.a( obj.n*(j0-1) + i0+1) .* ybar .* (1.-xbar) ...
                    + obj.a( obj.n*(j0+1-1) + i0) .* xbar .* (1.-ybar) ...
                    + obj.a( obj.n*(j0+1-1) + i0+1) .* xbar .* ybar );
            end            
        end
    end
end
         

            
