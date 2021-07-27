function x = iterate(phi, xn, Nmax)
    if (Nmax == 0); x = xn; return; end
    x = iterate(phi, phi(xn), Nmax - 1);
end