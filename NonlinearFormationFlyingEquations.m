function DX = NonlinearFormationFlyingEquations(t,X,mu,a,ecc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is provided 'as-is' for use with AOE 5234: Orbital Mechanics.
% Numerical accuracy is guaranteed only to within the bounds specified by
% The MathWorks Inc.
%
% Author: Andrew Rogers, 10 February 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time varying coefficients and solution to Kepler's equation
n = sqrt(mu/a^3);
M = n.*t;
M = M(:);
[EccAnom,~] = kepler(M,ecc);
f = 2.*atan(sqrt((1+ecc)./(1-ecc)).*tan(EccAnom./2));
r = (a.*(1-ecc.^2))./(1+ecc.*cos(f));
fdot = sqrt(mu.*a.*(1-ecc.^2)).*(1+ecc.*cos(f)).^2./(a.^2.*(1-ecc.^2).^2);
rdot = ecc.*sin(f).*sqrt(mu.*a.*(1-ecc.^2))./(a.*(1-ecc.^2));
fddot = -2.*rdot.*fdot./r;

% Deputy radius
x = X(1); y = X(2); z = X(3);
xd = X(4); yd = X(5); zd = X(6);
rd = ((r+x).^2+y.^2+z.^2).^(1/2);

% State space form of equations of motion (NERM)
dx1 = xd;
dx2 = yd;
dx3 = zd;
dx4 = 2.*fdot.*yd + fddot.*y + fdot.^2.*x + mu./r.^2 - mu.*(x+r)./rd.^3;
dx5 = -2.*fdot.*xd - fddot.*x + fdot.^2.*y - mu.*y./rd.^3;
dx6 = -mu.*z./rd.^3;

DX = [dx1 dx2 dx3 dx4 dx5 dx6]';
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
tol = 1e-12;
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