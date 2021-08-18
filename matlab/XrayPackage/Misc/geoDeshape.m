function geosO =  geoDeshape(geos)
% Reshapes a matrix of geodesics from (time index, array of geos) to (geos index, time index) 

    inDim = size(geos); %size(1) = numSteps, size(2:end) = sizeGeos;
    
    geosO = zeros([prod(inDim(2:end)),inDim(1)]);
    
    for (i = 1:inDim(1))
        geosO(:,i) = geos(i,:);
    end    
end
