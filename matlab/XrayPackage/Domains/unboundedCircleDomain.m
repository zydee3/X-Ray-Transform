classdef unboundedCircleDomain < circleDomain
    %UNBOUNDEDDOMAIN Summary of this class goes here
    %   Detailed explanation goes here
    methods
        
        function bool = isInside(~, ~, ~, ~)
            bool = true;
        end
        
        function bool = isInsideR2(~, ~, ~, ~)
            bool = true;
        end
        
    end
end

