function WithPower()
runs(1) = runs(1) + 1;

A = double(rand(1000))*100;

tic
    A = A.^n;
totalTime(1) = totalTime(1) + toc;
avgTime(1) = totalTime(1)/runs(1);

clear A
end