classdef FigureizerSettings < handle
    
    properties
        funcColormap = uint8([linspace(256,100,256);...
                      linspace(256,100,256);...
                      linspace(256,100,256)]')
        deprojSaturation = 0.5    
                  
        domainColor = uint8([0, 0, 255])
        domainStyle = '-'
        domainThickness = 0.8
        
        penColor = uint8([255, 0, 0])
        penStyle = '*'
        penThickness = 0.6
        
        gridResolution = 150;
        plotCenterType = 'default' % or use 'domain'
        
        domainResolution = 150;
        betaSpacing = 'default' % or use 'arclength' or 'euclid'
        
    end
    
    methods (Static)
        function sett = classic()
            sett = FigureizerSettings;
        end
        
        function sett = autumn()
            sett = FigureizerSettings;
        end
    end
    
    methods
        function swapSettings(obj,newSettings)
            
        end
        
        function reset(obj)
            obj.swapSettings(FigureizerSettings)
        end
        
    end
    
end

