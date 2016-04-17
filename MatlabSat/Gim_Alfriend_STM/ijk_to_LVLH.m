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