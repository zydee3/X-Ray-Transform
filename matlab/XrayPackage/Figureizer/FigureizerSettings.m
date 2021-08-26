classdef FigureizerSettings < handle
    
    properties
        funcColormap = uint8([linspace(150,256,256)',...
                      linspace(150,256,256)',...
                      linspace(150,256,256)'])
                  
        domainColor = uint8([0, 0, 255])
        penColor = uint8([255, 0, 0])
        linestyle = '-'
        pointstyle = '*'
        
        gridResolution = 500;
        plotCenterType = 'default' % or use 'domain'
        
        domainResolution = 250;
        betaspacing = 'default' % or use 'arclength' or 'euclid'
        alnormalspacing = 'default' % or use 'arclength', 'betas' or 'euclid'
        
    end
    
    methods (Static)
        function sett = default()
            sett = FigureizerSettings;
        end
        
        function sett = summer()
            sett = FigureizerSettings;
        end
    end
    
    methods
        function swapSettings(newSettings)
            
        end
        
    end
    
end

