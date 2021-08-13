        function obj = build(varargin) %TODO: Delete this grabage.
            
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

            parse(p,varargin{:})
            r = p.Results;
            stin.type = r.type;
            stin.theta = r.theta;
            stin.originX = r.origin(1);
            stin.originY = r.origin(2);
            stin.args = p.Unmatched;

            obj = Domain.build(stin);
        end    
       