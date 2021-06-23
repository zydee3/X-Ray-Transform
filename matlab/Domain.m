classdef Domain
    %DOMAIN Parent class of all other domains
    
    properties
        bdr = @(th) 1;
        dbdr = @(th) 0;
        ddbdr = @(th) 0;
        
        rMax = 1;
        
        originX = 0; % Handle these values after construction
        originY = 0;
        theta = 0;
    end
    
    properties (Dependent)
        alNorm
    end
    
    methods (Static)
        
       function mustBeDomain(obj)
           if (~isa(obj,'Domain'))
               error("Value must be a Domain");
           end
       end    
       
    end    
    
    methods
        function obj = Domain(argA, argB, argC) % -- !! TODO: i want to rework this constructor to accept more flavours of input
        if (strcmp(class(obj), 'Domain')) %%% -- !! TODO: hack to avoid subclass calling superclass constructor, should prolly change !!    
            % -- !! Todo, accept origin as an arg !!

            if (nargin == 0), return; end
            
            if (nargin == 1)   
                if (isnumeric(argA)) % argA = \theta
                    obj.theta = argA;
                    return;            
  
                    
                    
                elseif (isa(argA,'struct')) % converts from struct
                    % !! TODO: Exception for not having vvvv these vvvv fields !!
                    obj.bdr   = argA.bdr;
                    obj.dbdr  = argA.dbdr;
                    obj.ddbdr = argA.ddbdr;
                    obj.rMax  = argA.rMax;
                                            
                    if (isfield(argB, 'originX')), obj.originX = argB.origin; else, obj.originX = [0,0]; end  
                    if (isfield(argB, 'originY')), obj.originY = argB.origin; else, obj.originY = [0,0]; end 
                    if (isfield(argB, 'theta')), obj.theta = argB.theta; else, obj.originY = [0,0]; end 
                    return;
                end
                
                
                
            elseif (nargin == 2) 
                if (isnumeric(argA) && isnumeric(argB)) % argA,B = originX,Y
                    obj.originX = argA;
                    obj.originY = argB;
                    return; 
                    
                    
                    
                elseif (isa(argA,'function_handle') && nargin(argA) == 1 && isnumeric(argB)) % argA = r(\theta), argB is an upper bound of r
                    obj.bdr   = argA;
                    obj.dbdr  = @(x) deriv(argA, x);
                    obj.ddbdr = @(x) dderiv(argA,x);
                    obj.rMax  = argB;
                    return;    
                end    
                
            
                
            elseif (nargin == 3)    
                if (isnumeric(argA) && isnumeric(argB) && isnumeric(argC)) % argA = theta, argB,C = originX,Y
                    obj.theta = argA;
                    obj.originX = argB;
                    obj.originY = argC;
                    return; 
                end    
            end
            
            
            error("Bad input arguments.")
            
                       
            
        end,end   
        
        
        
        function obj = transform(obj, argA, argB, argC) % !! TODO: This bit of code is redudnant with the constructor, fix !!
            if (nargin == 1 && isnumeric(argA))   
                obj.theta = argA;
                return;   
            elseif (nargin == 2 && isnumeric(argA) && isnumeric(argB))   
                obj.originX = argA;
                obj.originY = argB;
                return;
            elseif (isnumeric(argA) && isnumeric(argB) && isnumeric(argC))       
                obj.theta = argA;
                obj.originX = argB;
                obj.originY = argC;
                return;
            end
            error("Bad input arguments.")
        end    
            
        
        
        function aln = get.alNorm(obj) % !! TODO: Bad implementation, initialize this in the constructor !!
            aln = @(th) angle(obj.bdr(th).*cos(th) + obj.dbdr(th).*sin(th) + 1i *(obj.bdr(th).*sin(th) - obj.dbdr(th).*cos(th))); 
        end    
        
        
        
        function out = plot(obj) % !! TODO: make nicer lool!!
            %plot Displays boundry
            
            n = 200;
            th = 2*pi*(1:n+1)/n;
            r = obj.bdr(th);
            pointX = cos(th) .* r;
            pointY = sin(th) .* r;
            
            out = plot(pointX,pointY);
        end
        
               
        
    end
end