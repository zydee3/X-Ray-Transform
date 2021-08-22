function [out0,out1] = atan2near(Y0,X0, Y1,X1)
%ATAN2PAIR computes atan2 so that the output arguments are within pi of
%each other
%   Values lie on [-pi,2pi]
    out0 = atan2(Y0,X0);
    out1 = atan2(Y1,X1);
    pi2 = 2*pi;
    
    ind = (abs(out0-out1-pi2) < pi);
    out1(ind) = out1(ind)+pi2;
    ind = (abs(out0+pi2-out1) < pi);
    out0(ind) = out0(ind)+pi2;

end

