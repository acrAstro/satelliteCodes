function [A,B] = hcwmatrices(meanMotion, numInputs)

n = meanMotion;

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0];
if numInputs == 3
    B = [zeros(3); eye(3)];
elseif numInputs == 2
    B = [zeros(4,2); eye(2)];
else
    error('Improper input matrix defined');
end
end