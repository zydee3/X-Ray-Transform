function [xO] = iterate2(phi, xn1,xn2, Nmax)
    if (Nmax == 0); xO = xn2; return; end
    xO = iterate2(phi, xn2,phi(xn1,xn2), Nmax - 1);
end