classdef aabDomain < Domain
    
    properties
        relativeB (1,2) {mustBeNumeric,mustBePositive} = [1,2]
    end
    
    methods
        function obj = aabDomain(args)
            arguments
                args.maxB (1,2) {mustBeNumeric} = [1,2];
                args.minB (1,2) {mustBeNumeric} = [0,0];
            end
            obj.relativeB = abs(args.maxB-args.minB)/2;
            
            origin = min(args.maxB,args.minB) + obj.relativeB;
            obj.originX = origin(1);
            obj.originY = origin(2);
               
        end 
        
        function out = bdr(obj,th)
            b = obj.relativeB;
            out = min( b(1)*abs(sec(th)), b(2)*abs(csc(th)) );
        end
        
        function out = dbdr(obj,th)
            b = obj.relativeB;
            cth = cos(th); sth = sin(th);
            w = b(1)*abs(1./cth);
            h = b(2)*abs(1./sth);
            
            bool = (w<h);
            
            out = zeros(size(th));
            ind = find(bool);
            out(ind) = w(ind).*sth(ind)./cth(ind);
            ind = find(~bool);
            out(ind) = -h(ind).*cth(ind)./sth(ind);
        end
        
        %{
        function out = ddbdr(obj,th)
            pion = pi/obj.sides;
            th = mod(th,2*pion)-pion;
            out = obj.radius * cos(pion) * (sin(th).^2+1)./cos(th).^3;
        end
        %}
        
        function [minB,maxB] = getBoundingBox(obj) 
            th0 = obj.theta;
            b = obj.relativeB;

            if (mod(th0,2*pi) == 0)
                maxB = b;
                minB = -b;
                return;
            end
                
            mr = sqrt(b(1).*b(1) + b(2).*b(2));
            at = atan2(b(2),b(1));
            
            uy = mr * max(abs(sin(th0 + at)), abs(sin(th0 - at)));
            ux = mr * max(abs(cos(th0 + at)), abs(cos(th0 - at)));
            
            maxB = [ux,uy];
            minB = -maxB;
        end
        
        
        
        function minR = getMinRadius(obj) 
            minR = min(obj.relativeB(1),obj.relativeB(2));
        end
        
        function maxR = maxRadius(obj) 
            bx = obj.relativeB(1);   by = obj.relativeB(2);
            maxR = sqrt(bx.*bx + by.*by);
        end
    end    
    
end


