function mustBe2DHandle(obj)
           %MUSTBE2DHandle Errors if the passed object is not a function
           %handle that takes two arguments
    if (~isa(obj,'function_handle')), error('Value must be a function.'), end
    if (nargin(obj) ~= 2), error('Function must take 2 arguments.'), end

end
