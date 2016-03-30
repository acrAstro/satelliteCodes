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

if strcmp(flag,'r') == 1
    rpqw = [p*cos(f)/(1+ecc*cos(f));
        p*sin(f)/(1+ecc*cos(f));
        0];
    vpqw = [-sqrt(mu/p)*sin(f);
        sqrt(mu/p)*(ecc+cos(f));
        0];
    map = rotation3(-O,flag)*rotation1(-inc,flag)*rotation3(-w,flag);
elseif strcmp(flag,'d') == 1
    rpqw = [p*cosd(f)/(1+ecc*cosd(f));
        p*sind(f)/(1+ecc*cosd(f));
        0];
    vpqw = [-sqrt(mu/p)*sind(f);
        sqrt(mu/p)*(ecc+cosd(f));
        0];
    map = rotation3(-O,flag)*rotation1(-inc,flag)*rotation3(-w,flag);
else
end
r = map*rpqw; v = map*vpqw;
if length(r)~= 3, error('Wrong size Position vector')
end
if length(v)~= 3, error('Wrong size Velocity vector')
end

end

function [R,V] = HillFrameToInertial(R1,V1,X)
X = X(:);
Rot = ijk_to_LVLH([R1(:);V1(:)]);

R = Rot'*X(1:3);
V = Rot'*X(4:6);
end

function RSW = ijk_to_LVLH(ic)

% Converts the inertial frame to the orbit frame (Local-Vertical,
% Local-Horizontal). Page 163, 172 of Fundamentals of Astrodynamics with
% Applications, David Vallado, 2007.

rvec = [ic(1) ic(2) ic(3)]'; vvec = [ic(4) ic(5) ic(6)]';
r_hat = (1/norm(rvec)).*rvec;
w_hat = (1/(norm(cross(rvec,vvec)))).*cross(rvec,vvec);
s_hat = cross(w_hat,r_hat);

RSW = [r_hat(:), s_hat(:), w_hat(:)]';

% Make sure it's a good rotation matrix. Needs full rank and unity
% determinant
if round(det(RSW)) ~= 1 || rank(RSW) ~= 3
    error('Inappropriately-sized matrix')
else
end
end

function R2 = rotation2(theta,flag)
% This function performs a 2-rotation, takes arguments in radians
if strcmp(flag,'d') == 1
    theta = theta*180/pi;
    R2 = [cosd(theta) 0 -sind(theta)
        0 1 0
        sind(theta) 0 cosd(theta)];
elseif strcmp(flag,'r') == 1
    R2 = [cos(theta) 0 -sin(theta)
        0 1 0
        sin(theta) 0 cos(theta)];
else
end
end

function R3 = rotation3(theta,flag)
% This function performs a 3-rotation, takes arguments in radians
if strcmp(flag,'d') == 1
    theta = theta*180/pi;
    R3 = [cosd(theta) sind(theta) 0;
        -sind(theta) cosd(theta)  0;
        0                 0    1];
elseif strcmp(flag,'r') == 1
    R3 = [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta)  0;
        0                 0    1];
else
end
end