isize = 50;
rsize = 256;
values = peaks(isize);

gf = GridFunc(values, interpType = 'nearest', offset = [isize,isize]/2);
m = Map(values);

[gx,gy] = ndgrid(linspace(1,isize,isize));
gi = griddedInterpolant(gx,gy,values);
gi.Method='nearest';

%[sx,sy] = ndgrid(linspace(1,isize,rsize));

sx = rand(rsize)*50;
sy = rand(rsize)*50;


%figure
tic
e=m.val(sx,sy);
toc
%s=surf(sx,sy,e);
%s.EdgeColor = 'none';


tic
e=gf.eval(sx,sy);
toc



%figure
tic
e=gi(sx,sy);
toc
%s=surf(sx,sy,e);
%s.EdgeColor = 'none';


