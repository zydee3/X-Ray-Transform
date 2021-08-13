classdef splineDomain < Domain
    %SPLINEDOMAIN The domain to end all domains
    
    properties %(Access = 'private')
        intervals;
    end
    
    properties 
        cycles;
    end
    
    methods
        function obj = splineDomain(args)
            arguments
                args.verticies (:,3) {mustBeNumeric} = 4;
                args.cycles (1,1) {mustBeInteger} = 1;
            end
            v = args.verticies;
            
            int = zeros(height(v),4);
            int(:,1) = v(:,1);
            int(:,2) = 0.5 * pi + atan2(v(:,2).*sin(v(:,1)) - v([2:end,1],2).*sin(v([2:end,1],1)),...
                                        v(:,2).*cos(v(:,1)) - v([2:end,1],2).*cos(v([2:end,1],1)) );
            int(:,3) = v(:,2) .* cos(int(:,1) - int(:,2));
            
            obj.intervals  = int;
            obj.cycles = args.cycles;
        end
        
        
        function out = bdr(obj,th) % lol this is awful
            int = obj.intervals;
            
            S = size(th);
            th = th(:)';
            
            th = mod(th, 2*pi/obj.cycles);
            %[~,index] = max( th' > repelem(int(:,1)',length(th),1) , [],2);
                           
            %out = (int(index,3)' .* sec(th - int(index,2)'));
            out = sum((int(:,3) .* sec(repelem(th,height(int),1) - int(:,2))) .* (th >= int(:,1) & th <= int([2:end,1],1)));
            out = reshape(out,S);
        end
        
        
        
        
    end
end

%{
piecewise scheme
input:
[th_n, r_n, br_n; 
th_{n+1}, r_{n+1}, br_{n+1}]
* it must be that th_{n+1} > th_n, 
* clamp br_n under min(|th_n-th_{n-1}|, |th_n-th_{n+1}|)
* wrap last vertex+1 with the first vertex

storage:
(for every interval [th_n,th_{n+1}])
[th_n, sec shift, sec coeff, br_n (for the vertex at n if br_n>0)], [poly coeffs]

%}