function WithTimes()
runs(2) = runs(2) + 1;

A = double(rand(1000))*100;

tic
for i = 1:n
    A = A.*A;
end
totalTime(2) = totalTime(2) + toc;
avgTime(2) = totalTime(2)/runs(2);

clear A i
end