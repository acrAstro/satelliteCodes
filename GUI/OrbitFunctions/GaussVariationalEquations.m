function DX = GaussVariationalEquations(t,X,req,mu,J2flag,flag)

% This function contains the necessary procedures for computing Gauss'
% Variation of Parameters. Equations derived from Battin (1987), Schaub and
% Junkins (2009) and Vallado (2007). This function also contains the
% differential equation used to describe the relative motion in cartesian
% coordinates of a chief and deputy. GVE give the chief's state, and the
% formation flying equations give the deputy's. J2 is added to provide the
% disturbance vector field but can be turned off using the J2flag option.
% The flag option calls angle arguments in either radians ('r') or degrees
% ('d'). The vector 'u' will take the place of any control vector field.

% Global options to be used throughout analysis
% global mu;

% Extract the Orbital Elements from the input vector to make typing/
% reading the equations easier.
a = X(1); ecc = X(2); inc = X(3); 
W = X(4);  w = X(5); f = X(6);
% p = a*(1-ecc^2);

% rho = [x y z]';

% More parameters calculated for convenience
p = a*(1-ecc^2);
r = p/(1+ecc*cos(f));
h = sqrt(mu*p);
th = w+f;

% The acceleration vector is written in terms of the position of the
% satellite, so convert the orbital elements into the R and V vectors using
% oe2rv subroutine.
els = [a ecc inc W w f];
[R,V] = oe2rv(els,mu,flag);
rx = R(1); ry = R(2); rz = R(3);
Rot = ijk_to_LVLH([R V]);
% Rdep = R + Rot^(-1)*rho;
% rdx = Rdep(1); rdy = Rdep(2); rdz = Rdep(3);

if strcmp(J2flag,'yes')
    AC = Rot*j2vec(rx,ry,rz,req,mu);
%     AJ2d = Rot*j2vec(rdx,rdy,rdz);
%     AD = AJ2d-AC;
elseif strcmp(J2flag,'no')
    AC = zeros(3,1);
%     AJ2d = zeros(3,1);
%     AD = AJ2d-AC;
else
end
ar = AC(1); at = AC(2); ah = AC(3);
% Gauss' Variational Equations
de(1) = (2*a^2/h)*(ecc*sin(f)*ar+(p/r)*at);
de(2) = (1/h)*(p*sin(f)*ar+((p+r)*cos(f)+r*ecc)*at);
de(3) = (r*cos(th)/h)*ah;
de(4) = r*sin(th)/(h*sin(inc))*ah;
de(5) = (1/(h*ecc))*(-p*cos(f)*ar+(p+r)*sin(f)*at)-(r*sin(th)...
    *cos(inc))/(h*sin(inc))*ah;
de(6) = h/r^2+1/(ecc*h)*(p*cos(f)*ar-(p+r)*sin(f)*at);

DX = de(:);
end

function A = j2vec(rx,ry,rz,req,mu)

r = norm([rx ry rz]);
J2 = 1082.63e-6;
J2Factor = -3/2*J2.*(mu./r.^2).*(req./r).^2;
A =  J2Factor*[(1-5.*(rz./r).^2).*(rx./r);
    (1-5.*(rz./r).^2).*(ry./r);
    (3-5.*(rz./r).^2).*(rz./r)];
end

function [r,v] = oe2rv(elements,mu,flag)
% This function converts the Keplerian elements to the Cartesian state
% vector. From Fundamentals of Astrodynamics with Applications,...
% David Vallado, 2007., pp.126-127, Algorithm 10.
%
% AC Rogers 14 Feb 2013
%
% Inputs:
% p   = parameter
% ecc = eccentricity
% inc = inclination (in radians)
% O   = right-ascension of the ascending node
% w   = argument of periapsis
% nu  = true anomaly
%
% Outputs:
% r   = 3x1 vector of current Cartesian positions (ECI frame)
% v   = 3x1 vector of current Cartesian velocities (ECI frame)

% Set default value of mu to that of Earth if no value is given
if nargin < 2
    mu = 3.986e5;
end
% Determine if vector of orbit elements is the correct size and return an
% error if it is not
if length(elements) ~= 6
    error('Orbit elements vector is not the correct size.')
end
a = elements(1); ecc = elements(2); inc = elements(3);
O = elements(4); w = elements(5);  f = elements(6);
p = a*(1-ecc^2);

rpqw = [p*cos(f)/(1+ecc*cos(f));
    p*sin(f)/(1+ecc*cos(f));
    0];
vpqw = [-sqrt(mu/p)*sin(f);
    sqrt(mu/p)*(ecc+cos(f));
    0];
map = rotation3(-O,flag)*rotation1(-inc,flag)*rotation3(-w,flag);

r = map*rpqw; v = map*vpqw;
if length(r)~= 3, error('Wrong size Position vector')
end
if length(v)~= 3, error('Wrong size Velocity vector')
end

end
