function dxn = deriv(f, xn) % numerical 1st derivative function, 4th order convergent with eps
    eps = 2^-10; % Heuristically chosen epsilon 
    dxn = (8*(f(xn+eps) - f(xn-eps)) + f(xn-2*eps) - f(xn+2*eps))/(12*eps);
end