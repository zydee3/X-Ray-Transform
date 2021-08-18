clc; close; clear;

repeat = 2000;

Times = zeros(1,10);
numbsA = rand(250);
numbsB = rand(250);
result = zeros(250);

for i = 1:repeat
    times = zeros(1,10);
   
    tic
    result = numbsA .* numbsB;
    times(2) = toc;
   
    tic
    result = numbsA + numbsB;
    times(1) = toc;

    tic
    result = numbsA - numbsB;
    times(3) = toc;

    tic
    result = numbsA ./ numbsB;
    times(4) = toc;

    tic
    result = exp(numbsB .* log(numbsA));
    times(6) = toc;
    
    tic
    result = numbsA .^ numbsB;
    times(5) = toc;
    
    Times = Times + times;
end
Times