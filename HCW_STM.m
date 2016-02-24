function x = HCW_STM(x0,n,t0,tf,steps)

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


t = linspace(t0,tf,steps);
x = zeros(length(t),6);
x(1,:) = x0;
for ii = 1:length(t)-1
    x(ii+1,:) = [4-3*cos(n*t(ii)) 0 0 1/n*sin(n*t(ii)) 2/n - 2/n*cos(n*t(ii)) 0;
        -6*n*t(ii)+6*sin(n*t(ii)) 1 0 -2/n+2/n*cos(n*t(ii)) 4/n*sin(n*t(ii))-3*t(ii) 0;
        0 0 cos(n*t(ii)) 0 0 1/n*sin(n*t(ii));
        3*n*sin(n*t(ii)) 0 0 cos(n*t(ii)) 2*sin(n*t(ii)) 0;
        -6*n+6*n*cos(n*t(ii)) 0 0 -2*sin(n*t(ii)) -3+4*cos(n*t(ii)) 0;
        0 0 -n*sin(n*t(ii)) 0 0 cos(n*t(ii))]*x0;
end
end