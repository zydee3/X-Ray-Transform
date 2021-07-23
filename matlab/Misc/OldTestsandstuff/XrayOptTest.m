dom = circleDomain(radius=2);
met = hyperbolicMetric(radius=3);

im = InMap([0,0,0,0,0;...
            0,0,-2,0,0;...
            0,0,0,1,0;...
            0,2,0,0,0;
            0,0,0,0,0;]);
    im = im.transform(0,0,0.7,0.7,0);       
    
surface = RiemannSurface(dom,met, stepType='RK4', stepSize=0.01, geoDur=20);





pow = 8;
pow2 = 2^pow;
times = zeros(pow+1);

B = linspace(0,2*pi,1);
A = linspace(-pi,pi,pow2*100)*0.5;
vals = zeros(pow2);

for bpart = 0:0
    for apart = 0:pow
        
        vals = zeros(pow2);
        
        bsize = 1;
        asize = (pow2/(2^apart));
        
        tic
        for bpos = (1:(2^bpart))-1
            for apos = (1:(2^apart))-1

                [beta,alpha] = meshgrid(B((1:bsize)+bsize*bpos),...
                                        A((1:asize)+asize*apos));   
                                                     
                 vals((1:bsize)+bsize*bpos,...
                      (1:asize)+asize*apos) = XrayI0(im,surface,beta,alpha);
                
            end
        end
        times(bpart+1,apart+1) = toc
        
        

        pause(1);
    end
end

s=surf(times);
s.EdgeColor = 'none';
s.FaceColor = 'interp';

