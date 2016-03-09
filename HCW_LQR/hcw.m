function xdot = hcw(t,x,K,n,numInputs)

[A,B] = hcwmatrices(n,numInputs);
if isempty(K)
    u = zeros(numInputs,1);
else
    u = -K*x;
end
xdot = A*x + B*u;
end