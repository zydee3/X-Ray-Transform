function TheExpThing()
runs(3) = runs(3) + 1;

A = double(rand(1000))*100;

tic
    A = exp(n * log(A));
totalTime(3) = totalTime(3) + toc;
avgTime(3) = totalTime(3)/runs(3);

clear A
end