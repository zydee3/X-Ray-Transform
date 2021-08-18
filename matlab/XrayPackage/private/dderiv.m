
function ddxn = dderiv(f, xn) % numerical 2nd derivative function, 5th order convergent with eps
    eps = 2^-8; % Heuristically chosen epsilon 
    ddxn = ((f(xn+3*eps)+f(xn-3*eps))/90 + 3*(f(xn+eps)+f(xn-eps) - (f(xn+2*eps)+f(xn-2*eps))/10)/2 - 49*f(xn)/18)/(eps*eps);
end

%{
function ddxn = dderiv(f, xn) % numerical 2nd derivative function, 3th order convergent with eps
    eps = 2^-9; % Heuristically chosen epsilon 
    ddxn = (8*(f(xn+eps)+f(xn-eps))/3 - (f(xn+2*eps)+f(xn-2*eps))/6 - 5*f(xn))/(2*eps*eps);
end
%}

