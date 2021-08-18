classdef TFunc
    %INFUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = 'protected')
        valuesA
        tmatrix = eye(3);
    end
    
    properties (Dependent)
        eval 
    end
    
    
    methods
        function obj = TFunc(values, args)  
            %INMAP Constructs an instance of TFunc
            arguments
                values = @(x,y) zeros(size(x));   
                args.offset (1,2) {mustBeNumeric} = [0,0]
                args.scale (1,2) {mustBeNumeric} = [1,1]
                args.rotation (1,1) {mustBeNumeric} = 0
            end
            
            obj = obj.transform(args.offset(1),args.offset(2), args.scale(1),args.scale(2), args.rotation);
            obj.valuesA = values; 
        end
        
        
        
        function obj = transform(obj, xoff,yoff, xscale,yscale, rot)
            %TRANSFORM Takes arguments to rotate, scale and offset the
            %function.
            arguments
                obj
                xoff (1,1) {mustBeNumeric} = 0
                yoff (1,1) {mustBeNumeric} = 0
                xscale (1,1) {mustBeNumeric} = 1
                yscale (1,1) {mustBeNumeric} = 1
                rot (1,1) {mustBeNumeric} = 0
            end
            
            obj.tmatrix = inv([cos(rot)*xscale, sin(rot)*xscale,0;...
                              -sin(rot)*yscale, cos(rot)*yscale,0;...
                               xoff,            yoff,           1]);
        end
               
        %{ 
       function obj = bake(obj, xres,yres, xmin,ymin,xmax,ymax)
            %BAKE Converts the type of InMap from a function handle to an
            %image.
            
            arguments
                obj
                xres (1,1) {mustBeNumeric} = 0.2
                yres (1,1) {mustBeNumeric} = 0.2
                xmin (1,1) {mustBeNumeric} = -2
                ymin (1,1) {mustBeNumeric} = -2
                xmax (1,1) {mustBeNumeric} = 2
                ymax (1,1) {mustBeNumeric} = 2                
            end
            
            [X,Y] = meshgrid(xmin:1/xres:xmax,ymin:1/yres:ymax);
            obj.values = obj.values(obj.eval(X,Y));
            obj.tmatrix = inv([xscale, 0,0;...
                               0, yscale,0;...
                               0.5*(xmin+xmax),0.5*(ymin+ymax),1]);
        end  
        %}
        
        function out = plot(obj, lbound,ubound, resolution)
            %PLOT Plots InMap on the axis aligned box defined by a lower
            %and upper bound.
            %   Uses obj.eval.
            
            arguments
                obj
                lbound (1,2) {mustBeNumeric} = [-2,-2]
                ubound (1,2) {mustBeNumeric} = [2,2]
                resolution (1,1) {mustBeNumeric} = 100
            end
            
            [X,Y] = meshgrid(linspace(lbound(1),ubound(1),resolution),...
                             linspace(lbound(2),ubound(2),resolution));
                         
            Z = obj.eval(X,Y);
            hold on;
            out = pcolor(X,Y,Z);
            regT = inv(obj.tmatrix);
            plot(regT(3),regT(6),'r*');
            out.EdgeColor = 'none';
            
        end
        
        
        function value = get.eval(obj)
           value = @(X,Y) obj.evaluate(X,Y);
        end
        
    end
    
    
    methods (Access = 'protected')
        
        function out = evaluate(obj, X,Y)
            %EVAL Evaluates the function at an array of points (X,Y).
            
            tmat = obj.tmatrix;

            x = X*tmat(1) + Y*tmat(2) + tmat(3); % matrix multiplication (when handled correctly, only these 6 indicies are necessarily unknown, and the resulting vector will always have a 1 in the third place)
            y = X*tmat(4) + Y*tmat(5) + tmat(6);
            % z = ones(size(X)) is what would be written here if we cared abt this index 
            
            out = obj.valuesA(x,y);

        end        
        
        
    end    
    
end


