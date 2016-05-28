classdef RelativeMotion < handle
%% RelativeMotion class
%
% This is a simple class to propagate the linear and nonlinear equations of
% relative motion for spacecraft formation flying.
%
% Author: Andrew Rogers
% Date:   11 April 2016
%
    properties
    mu          % Gravitational parameter of Earth
    kepElems    % Initial Kepler elements of chief
    t0          % Initial time
    numPeriod   % Number of orbital periods
    period      % Period of orbit, (in seconds)
    time        % Time vector
    relState    % Relative initial conditions
    XL          % Linear relative state history
    XN          % Nonlinear relative state history
    intOptions  % Integrator options
    end
    
    methods
        % Constructor method
        function sff = RelativeMotion(initStruct)
            sff.mu = initStruct.params{1};
            sff.kepElems = initStruct.params{2};
            sff.relState = initStruct.params{3};
            sff.t0 = initStruct.timeParams{1};
            sff.numPeriod = initStruct.timeParams{2};
            sff.intOptions = initStruct.options;
            sff.makeTimeVector();
        end
        
        % Makes the time vectory
        function sff = makeTimeVector(sff)
            sff.period = 2*pi/sqrt(sff.mu/sff.kepElems(1)^3);
            tf = sff.numPeriod*sff.period;
            sff.time = linspace(sff.t0,tf,200*sff.numPeriod);
        end
        
        % Propagates the linear equations of relative motion
        function sff = propagateLinear(sff)
            [~,sff.XL] = ode113(@LinearFormationFlyingEquations,sff.time,sff.relState,sff.intOptions,sff.mu,sff.kepElems(1),sff.kepElems(2));
        end
        
        % Propagates the nonlinear equations of relative motion
        function sff = propagateNonlinear(sff)
            [~,sff.XN] = ode113(@NonlinearFormationFlyingEquations,sff.time,sff.relState,sff.intOptions,sff.mu,sff.kepElems(1),sff.kepElems(2));
        end
    end
end

function DX = NonlinearFormationFlyingEquations(t,X,mu,a,ecc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function DX = LinearFormationFlyingEquations(t,X,mu,a,ecc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% State matrix for LERM
A = [0,0,0,1,0,0;
      0,0,0,0,1,0;
      0,0,0,0,0,1;
      (fdot^2+2*mu/r^3) fddot 0 0 2*fdot 0;
      -fddot (fdot^2-mu/r^3) 0 -2*fdot 0 0
      0 0 -mu/r^3 0 0 0];
DX = A*X;
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