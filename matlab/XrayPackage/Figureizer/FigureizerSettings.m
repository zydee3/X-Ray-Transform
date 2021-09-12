classdef FigureizerSettings < handle
    
    properties
        funcColormap = uint8([linspace(256,100,256);...
                      linspace(256,100,256);...
                      linspace(256,100,256)]')
        deprojSaturation = 0.5    
        forceMidZero = false;
                  
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
            sett.forceMidZero = false;
            sett.funcColormap = turbo;
            sett.domainColor = uint8([74, 95, 169]);
            sett.penColor = uint8([238, 42, 40]);
        end
        
        function sett = heat()
            sett = FigureizerSettings;
            sett.forceMidZero = true;
            sett.funcColormap = buildColorGradient([255,79, 37;...
                                                    255,178,49;...
                                                    255,226,143;...
                                                    255,255,255;...
                                                    192,228,255;...
                                                    122,193,255;...
                                                    0,  169,255],128);
            sett.domainColor = uint8([0, 0, 0]);
            sett.penColor = uint8([100, 100, 100]);
        end
        
        function sett = gray()
            sett = FigureizerSettings;
            sett.forceMidZero = false;
            sett.funcColormap = uint8([linspace(256,100,256);...
                                      linspace(256,100,256);...
                                      linspace(256,100,256)]');
            sett.domainColor = uint8([0, 0, 0]);
            sett.penColor = uint8([100, 100, 100]);
        end
        
        function sett = lines()
            sett = FigureizerSettings;
            sett.forceMidZero = false;
            sett.funcColormap = lines(64);%repmat(repelem([0,0,0;1,1,1],2,1),64,1) .* lines;
            sett.domainColor = uint8([255, 255, 255]);
            sett.penColor = uint8([200, 200, 200]);
        end
    end
    
    methods
        function swapSettings(obj,newSettings)
            arguments
                obj
                newSettings% (1,1) {FigureizerSettings.mustBeSett};
            end
            
            for fn = fieldnames(obj)'
                obj.(fn{1}) = newSettings.(fn{1});
            end
        end
        
        
        function clear(obj)
            obj.swapSettings(FigureizerSettings)
        end
        
    end
    
end

