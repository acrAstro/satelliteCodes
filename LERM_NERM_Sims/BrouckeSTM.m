function X = BrouckeSTM()



end

function [r,rdot,thetadot,thetaddot] = polarParams(T,mu,a,ecc,omega)

t0 = T(1);
t = T(2);
dt = t-t0;
n = sqrt(mu/a^3);
MeanAnom = n*dt;
[EccAnom,~] = kepler(MeanAnom,ecc);
f = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(EccAnom/2));
theta = omega + f;

r = a*(1-ecc^2)/(1+ecc*cos(theta));
rdot = sqrt(mu*a*(1-ecc^2))*(1+ecc*cos(theta)^2)/(a^2*(1-ecc^2)^2);
thetadot = ecc*mu*sin(theta)/(sqrt(mu*a*(1-ecc^2)));
thetaddot = 2*mu*ecc*(1+ecc*cos(theta))^3*sin(theta)/(a^3*(1-ecc^2)^3);
end

function [E,iter] = kepler(M,ecc)
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
for ii = 1:length(M)
    E_n = M(ii);
    f_n = E_n-ecc.*sin(E_n) - M(ii);
    while (abs(f_n) > tol)
        E_n = E_n - f_n./(1-ecc.*cos(E_n));
        f_n = E_n - ecc.*sin(E_n) - M(ii);
        iter = iter + 1;
    end
    E(ii) = E_n;
end
format short
end