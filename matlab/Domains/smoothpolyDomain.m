classdef smoothpolyDomain < Domain
    % Twice-differentiable version of PolygonDomain
    
    properties
        bevelRadius = 0
        radius = 2
        sides = 4
    end
    
    methods
        function obj = smoothpolyDomain(args)
            arguments
                args.radius (1,1) {mustBeNumeric} = 2;
                args.sides (1,1) {mustBeInteger} = 4;
                args.bevelRadius (1,1) {mustBeNumeric} = 0;
            end
            
            obj.radius = args.radius;
            obj.sides = args.sides;
            obj.bevelRadius = max(0,min(args.bevelRadius,pi/obj.sides));
            obj.rMax = obj.radius;
              
        end
        
        
        function outr = bdr(obj,th)
            r = obj.radius;
            br = obj.bevelRadius;
            
            pion = pi/obj.sides;
            ct = mod(th+pion,2*pion)-pion; % c
            bt = mod(th,2*pion)-pion; % b

            outr = (abs(ct) > br) .* r*cos(pion)./cos(bt); %!!!!!!!!!! vectorize this less bad, (use find?)
            if (br == 0), return, end;

            bR = mod(br,2*pion)-pion;
            q = (ct-br).*(ct+br);
            secbr = sec(bR);
            tanbr = tan(bR);

            constA = (br*(tanbr*tanbr+secbr*secbr)-tanbr)/(8*br*br*br);
            constB = tanbr/(2*br);

            outr = outr + (abs(ct) <= br) .* r*cos(pion)*secbr .* (constA.*q.*q + constB.*q + 1);
            %!!!!!!!!!!
        end


        %{
        function out = dbdr(obj,th) % VVVV TODO!! VVVV

        end

        function out = ddbdr(obj,th)

        end
        %}
    end    
    
end





