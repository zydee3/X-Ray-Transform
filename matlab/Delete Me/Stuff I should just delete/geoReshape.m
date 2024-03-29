function geosO =  geoReshape(geos, dimensions)
    % Reshapes a matrix of geodesics from (geos index, time index) to (time index, array of geos)

    inDim = size(geos); %inDim(1) = numGeos, inDim(2) = numSteps;

    if (~isvector(dimensions)), error('Dimensions should be a vector representing the size of an array.'); end
    if (inDim(1) ~= prod(dimensions)), error('Reshape should not change the number of geodesics.'); end

    geosO = zeros([inDim(2),dimensions]);
    
    for (i = 1:inDim(2))
        geosO(i,:) = geos(:,i);
    end   
end