       function obj = build(varargin) %TODO: Delete this grabage.
           obj = Metric();
            if (nargin == 0), return; 
            
            elseif (nargin == 1)
                varargin = varargin{1};
                if (isa(varargin,'struct'))
                    if (isfield(varargin, 'type') && isfield(varargin,'args')) % check if this is a parsed struct
                        if (~strcmp(varargin.type, 'default'))
                            name = lower(varargin.type) + "Metric"; 
                            if (exist(name,'class') == 8)

                                argNames = fieldnames(varargin.args);
                                celin = cell(1,length(argNames) * 2+1);
                                celin{1} = name;
                                for i = 1:numel(argNames)
                                    celin{i+1} = argNames{i};
                                    celin{i+2} = varargin.args.(argNames{i});
                                end

                                obj = feval(celin{:});

                            end
                        else    
                            obj.lg   = varargin.lg;
                            obj.dxlg = varargin.dxlg;
                            obj.dylg = varargin.dylg;
                            obj.curv = varargin.curv;
                        end
                        return;
                    end

                end
            elseif (nargin > 1) 
                  
            end
            p = inputParser; % it parses inputs
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            p.FunctionName = 'Metric.build';

            isAString = @(in) isstring(in) || ischar(in);
            isA2Handle = @(in) isa(in,'function_handle') && nargin(in) == 2;

            addOptional(p,'type','default',isAString);

            if (isa(varargin,'cell'))
                parse(p,varargin{:})
            else
                parse(p,varargin)
            end    
            r = p.Results;
            stin.type = r.type;
            stin.args = p.Unmatched;

            obj = Metric.build(stin);
        end   