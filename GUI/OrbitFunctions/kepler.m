function [E,iter] = kepler(M,e)
% Function solves Kepler's equation:
% M = n*(t-t_0) = E-e*sin(E)
% Using Newton-Raphson iteration
% AC Rogers, 21 February 2013
% Inputs:
%           M    = mean anomaly
%           e    = eccentricity
% Outputs:
%           E    = eccentric anomaly
%           iter = number of iterations
format long
tol = 1e-10;
iter = 1;
E_n = M; f_n = E_n-e*sin(E_n) - M;
while (abs(f_n) > tol)
    E_n = E_n - f_n/(1-e*cos(E_n));
    f_n = E_n - e*sin(E_n) - M;
    iter = iter + 1;
end

E = E_n;
format short
end