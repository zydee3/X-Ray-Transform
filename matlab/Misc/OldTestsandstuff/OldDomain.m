classdef Domain
    %DOMAIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bdr     
        dbdr
        ddbdr
        
        rMax
        originX
        originY
    end
    
    methods
        function obj = Domain(argA, argB) % -- !! i want to rework this constructor to accept more flavours of input
            
            % -- !! Todo, accept origin as an arg !!
            obj.originX = 0;
            obj.originY = 0;
            
            
            if (nargin == 0)        % Defaults to unit circle 
                params.a = 1;
                obj = Domain('circle', params);
                return;
                
                
                
            elseif (nargin == 1)   
                if (isa(argA,'function_handle') && nargin(argA) == 1) % argA = r(\theta)
                    obj.bdr = argA;
                    % compute the other variables 
                    return;
                    
                    
                    
                elseif (isa(argA,'struct')) % converts from struct
                    obj.bdr   = argA.bdr;
                    obj.dbdr  = argA.dbdr;
                    obj.ddbdr = argA.ddbdr;
                    obj.rMax  = argA.rMax;
                    if (isfield(argB, 'originX')), obj.originX = argB.origin; else, obj.originX = [0,0]; end  
                    if (isfield(argB, 'originY')), obj.originY = argB.origin; else, obj.originY = [0,0]; end 
                    return;
                end
                
                
                
            elseif (nargin == 2)   % argA = type, argB = params (argB.a,b,c) (pre-built boundries)
                if (isa(argB,'struct') && (isstring(argA) || ischar(argA)))
                    
                    if (isfield(argB, 'originX')), obj.originX = argB.origin; else, obj.originX = [0,0]; end  
                    if (isfield(argB, 'originY')), obj.originY = argB.origin; else, obj.originY = [0,0]; end 
                    if (~isfield(argB, 'th0')), argB.th0 = 0; end 
                        
                    switch argA % th0 = offset
                        case 'circle' % param.a = radius
                            obj.bdr = @(th) argB.a;
                            obj.dbdr = @(th) 0;
                            obj.ddbdr = @(th) 0;
                            obj.rMax = argB.a;
 
                            
                            
                        case 'elipse' % param.a = radiusA, param.b = radiusB
                            a = argB.a;
                            b = argB.b;
                            th0 = argB.th0;
                            obj.bdr = @(th) a*b ./ sqrt( (b*cos(th-th0)).^2 + (a*sin(th-th0)).^2 );
                            obj.dbdr = @(th) bdr(th).^3 .* sin(2*(th-th0)) .* (b^2-a^2)./2/(a*b)^2;
                            obj.ddbdr = @(th) 3*dbdr(th).^2./bdr(th) + bdr(th).^3 .* cos(2*(th-th0)) .* (b^2-a^2)./(a*b)^2;
                            obj.rMax = max(argB.a, argB.b);
  
                            
                            
                        case 'cos' % param.a = offest, param.b = amplitude, param.c = cycle count 
                            a = argB.a;
                            b = argB.b;
                            n = argB.c;
                            th0 = argB.th0;
                            obj.bdr = @(th) a + b*cos(n*(th-th0));
                            obj.dbdr = @(th) -n*b*sin(n*(th-th0));
                            obj.ddbdr = @(th) -n*n*b*cos(n*(th-th0));
                            obj.rMax = argB.a + argB.b;         
                            
                            
                            
                        otherwise 
                            error('wrong boundary type');
                    end,return;
                end
            
                
            elseif (nargin > 2)    
                
            end
            
            
            error("Bad input arguments.")
        end      
            
        
        
        function out = plot(obj) % !! TODO: make nicer lool!!
            %plot Displays boundry
            
            n = 100;
            th = 2*pi*(1:n+1)/n;
            r = obj.bdr(th);
            pointX = cos(th) .* r;
            pointY = sin(th) .* r;
            
            out = plot(pointX,pointY);
            
        end
        
        
    end
end

