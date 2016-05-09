function [mag,phi,beta] = thrustVector(U)

x = -U(1); y = -U(2); z = -U(3);

% Euclidean norm
mag = sqrt(x^2+y^2+z^2);
% Elevation angle
phi = acosd(x/mag);
% Azimuth angle
beta = atan2d(y,z);

end