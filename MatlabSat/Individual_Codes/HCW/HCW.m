function DX = HCW(t,X,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All units are in meters and meters/sec, instead of mu = 3.986e5 km^3/s^3,
% it is 3.986e14 m^3/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is provided 'as-is' for use with AOE 5234: Orbital Mechanics.
% Numerical accuracy is guaranteed only to within the bounds specified by
% The MathWorks Inc.
%
% Author: Andrew Rogers, 10 February 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A = [zeros(3,3) eye(3)
    3*n^2 0 0 0 2*n 0
    0 0 0 -2*n 0 0
    0 0 -n^2 0 0 0];

DX = A*X;
end