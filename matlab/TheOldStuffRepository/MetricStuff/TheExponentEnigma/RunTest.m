
for n = 1:10
   [runs,totalTime,avgTime] = InitializeTest();
   run(n,runs,totalTime,avgTime) 
   
end


% 1 -- ^n,   2 -- ****,  3 -- exp(n*log(n))

function A=generateA()
    A = rand(1000);
end

function avg = run(n,runs,totalTime,avgTime)
    for i = 1:100
        clc
        i
        pause(0.01)
        [n,runs,totalTime,avgTime] = WithTimes(n,runs,totalTime,avgTime);
        pause(0.01)
        [n,runs,totalTime,avgTime] = WithPower(n,runs,totalTime,avgTime);
        pause(0.01)
        [n,runs,totalTime,avgTime] = TheExpThing(n,runs,totalTime,avgTime);
    end
    avg = avgTime;
end




function [n,runs,totalTime,avgTime] = WithTimes(n,runs,totalTime,avgTime)
runs(2) = runs(2) + 1;

A = generateA();

tic
for i = 1:n
    A = A.*A;
end
totalTime(2) = totalTime(2) + toc;
avgTime(2) = totalTime(2)/runs(2);

clear A i
end


function [n,runs,totalTime,avgTime] = WithPower(n,runs,totalTime,avgTime)
runs(1) = runs(1) + 1;

A = generateA();

tic
    A = A.^n;
totalTime(1) = totalTime(1) + toc;
avgTime(1) = totalTime(1)/runs(1);

clear A
end


function [n,runs,totalTime,avgTime] = TheExpThing(n,runs,totalTime,avgTime)
runs(3) = runs(3) + 1;

A = generateA();

tic
    A = exp(n * log(A));
totalTime(3) = totalTime(3) + toc;
avgTime(3) = totalTime(3)/runs(3);

clear A
end