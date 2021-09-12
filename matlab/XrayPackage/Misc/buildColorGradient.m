function out = buildColorGradient(keyColors, steps)
    arguments
        keyColors (:,3)
        steps (1,1) = 10;
    end
    
    keyColors = double(keyColors);
    
    stepPoints = linspace(0,1, steps);
    keyPoints = linspace(0,1, height(keyColors));
    
    interpX = interp1(keyPoints,keyColors(:,1),stepPoints);
    interpY = interp1(keyPoints,keyColors(:,2),stepPoints);
    interpZ = interp1(keyPoints,keyColors(:,3),stepPoints);
    
    out = uint8(reshape([interpX;interpY;interpZ]',[],3));
end