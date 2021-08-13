function mustBeIntegerGreaterThan2(obj)
           %MUSTBEMETRIC Errors if the passed object is not an integer
           %greater than 2
    if (~isinteger(obj)), error('Value must be an integer.'), end
    if (obj > 2), error('Value must be geater than 2.'), end

end
