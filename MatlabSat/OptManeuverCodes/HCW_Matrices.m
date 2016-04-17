function [A,B] = HCW_Matrices(n,nu,mass)

a = 3*n^2;
b = 2*n;
c = -2*n;
d = n^2;

A11 = zeros(3,3); A12 = eye(3);
A21 = [a 0 0; 0 0 0; 0 0 d];
A22 = [0 b 0; c 0 0; 0 0 0];
A = [A11,A12; A21,A22];

switch nu
    case 2
        B = 1/mass.*[zeros(4,2); eye(2)];
    case 3
        B = 1/mass.*[zeros(3,3); eye(3)];
    case 4
        B = 1/mass.*[zeros(4,4); eye(2), -eye(2)];
    case 6
        B = 1/mass.*[zeros(3,6); eye(3), -eye(3)];
end

end