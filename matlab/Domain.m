classdef Domain < handle
    %DOMAIN Parent class of all other domains
    
    properties
        originX = 0; % Handle these values after construction
        originY = 0;
        theta = 0;
    end
    
    properties (Access = 'private')
        bdrVAR = @(th) ones(size(th));
    end    
    
    properties (SetAccess = 'protected')
        rMax;
    end    
    
    
    methods
        function out = bdr(obj,th)
            out = obj.bdrVAR(th);
        end
        
        function out = dbdr(obj,th)
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = dderiv(obj.bdrVAR, th);
        end
        
        function out = ddbdr(obj,th)
            warning('Default derivatives are very inefficient, consider overriding with explicit implementations.')
            out = dderiv(obj.bdrVAR, th);
        end
        
        
        
        function obj = Domain(bdr, rMax)
        if (strcmp(class(obj), 'Domain')) %%% -- !! TODO: hack to avoid subclass calling superclass constructor, should prolly change !!    
            if (nargin == 2 && ...
                isa(bdr,'function_handle') && nargin(lg) == 1 && ...
                isnumeric(rMax) && all(size(rMax)==1))
            
                obj.bdrVAR = bdr;
                obj.rMax = rMax;
            end
            
        end,end   
        
    
        function [bdr,dbdr,ddbdr] = getHandles(obj)
            bdr = @(t) obj.bdr(t);
            dbdr = @(t) obj.dbdr(t);
            ddbdr = @(t) obj.ddbdr(t);
        end 
        
        
        function out = alNorm(obj)
            aln = @(th) angle(obj.bdr(th).*cos(th) + obj.dbdr(th).*sin(th) + 1i *(obj.bdr(th).*sin(th) - obj.dbdr(th).*cos(th))); 
        end  
    
        
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
            
        
    
        function out = plot(obj) % !! TODO: make nicer lool!!
            %plot Displays boundry
            
            n = 200;
            th = 2*pi*(1:n+1)/n;
            th0 = obj.theta;
            x0 = obj.originX;
            y0 = obj.originY;
            r = obj.bdr(th);
            pointX = cos(th - th0) .* r + x0;
            pointY = sin(th - th0) .* r + y0;
            
            out = plot(pointX,pointY);
        end
        
               
        
    end
    
    methods (Static)
       function mustBeDomain(obj)
           if (~isa(obj,'Domain'))
               error("Value must be a Domain");
           end
       end     

       
       function obj = build(varargin)
            obj = Domain();
            if (nargin == 0), return;

            elseif (nargin == 1)
                varargin = varargin{1};
                if (isa(varargin,'struct'))
                    if (isfield(varargin, 'type') && isfield(varargin,'args')) % check if this is a parsed struct
                        if (~strcmp(varargin.type, 'default'))
                            name = lower(varargin.type) + "Domain"; 
                            if (exist(name,'class') == 8)
                                argNames = fieldnames(varargin.args);
                                celin = cell(1,length(argNames) * 2+1);
                                celin{1} = name;
                                for i = 1:(numel(argNames))
                                    celin{2*i} = argNames{i};
                                    celin{2*i+1} = varargin.args.(argNames{i});
                                end
                                obj = feval(celin{:});

                            end
                        
                        end
                    else    
                        obj.bdr   = varargin.bdr;
                        obj.dbdr  = varargin.dbdr;
                        obj.ddbdr = varargin.ddbdr;
                        obj.rMax  = varargin.rMax;
                    end
                    obj.theta = varargin.theta;
                    obj.originX = varargin.originX;
                    obj.originY = varargin.originY;
                    return;
                end
            elseif (nargin > 1) 
                
            end    


            p = inputParser; % it parses inputs
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            p.FunctionName = 'Domain.build';

            isAString = @(in) isstring(in) || ischar(in);
            isANumber = @(in) isnumeric(in);
            isA2Vector = @(in) isnumeric(in) && all(size(in) == [1,2]);
            isA1Handle = @(in) isa(in,'function_handle') && nargin(in) == 1;

            addOptional(p,'type','default',isAString);

            addParameter(p,'theta',obj.theta,isANumber);
            addParameter(p,'origin',[obj.originX,obj.originY],isA2Vector);

            addParameter(p,'bdr',@(t) obj.bdr(t),isA1Handle);
            addParameter(p,'dbdr',@(t) obj.dbdr(t),isA1Handle);
            addParameter(p,'ddbdr',@(t) obj.ddbdr(t),isA1Handle);
            addParameter(p,'rMax',obj.rMax,isANumber);

            parse(p,varargin{:})
            r = p.Results;
            stin.type = r.type;
            stin.theta = r.theta;
            stin.originX = r.origin(1);
            stin.originY = r.origin(2);
            stin.args = p.Unmatched;

            obj = Domain.build(stin);
        end    
    end
    
end