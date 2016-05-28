function [KepElems] = rv2oe(ic,mu)
% Function imports a 1x6 state vector of ephemeris data and computes the
% corresponding initial Keplerian state vector for use with Gauss'
% Equations, follows algorithm 9 (pp 120-121) in Fundamentals of ...
% Astrodynamics with Applications, David Vallado, 2007.

% Inputs:   ic         = Cartesian state vector [x,y,z,xd,yd,zd]
% Outputs:  KepElems   = 1x6 vector of Keplerian Elements
%
% Inertial unit vectors, IJK
KK = [0 0 1];
% Radial and velocity vectors
rvec = [ic(1) ic(2) ic(3)]; vvec = [ic(4) ic(5) ic(6)];
r = norm(rvec); v = norm(vvec);
% Angular momentum
hvec = cross(rvec,vvec); h = norm(hvec);
nvec = cross(KK,hvec); n = norm(nvec);
% Laplace-Runge-Lenz vector (eccentricity vector)
eccvec = (1/mu)*((v^2-mu/r).*rvec - (dot(rvec,vvec)).*vvec);
ecc = norm(eccvec);
% Specific mechanical energy
Energy = v^2/2 - mu/r;
% Compute semi-major axis, semi-minor axis and semi-latus rectum
if ecc ~= 1.0
    a = -mu/(2*Energy);
else
    a = inf;
end
b = a*sqrt(1-ecc^2);
% Inclination
inc = acosd((hvec(3))/h);
% Right-ascension of the ascending node
if nvec(2) > 0
    W = acosd((nvec(1))/n);
else
    W = 360 - acosd((nvec(1))/n);
end
% Argument of periapsis
if eccvec(3) > 0
    w = acosd((dot(nvec,eccvec))/(n*ecc));
else
    w = 360 - acosd((dot(nvec,eccvec))/(n*ecc));
end
% True anomaly
if dot(rvec,vvec) > 0
    f = acosd(dot(eccvec,rvec)/(ecc*r));
else
    f = 360 - acosd(dot(eccvec,rvec)/(ecc*r));
end

% Special cases
% Elliptic Equatorial
if inc == 0 && ecc ~= 0
    if eccvec(2) > 0
        w_true = acosd(eccvec(1)/ecc);
    else
        w_true = 360 - acosd(eccvec(1)/ecc);
    end
    KepElems = [a,ecc,inc,W,w_true,f]';
    % Circular Inclined
elseif inc ~= 0 && ecc == 0
    if rvec(3) > 0
        u_l = acosd(dot(nvec,rvec)/(n*r));
    else
        u_l = 360 - acosd(dot(nvec,rvec)/(n*r));
    end
    KepElems = [a,ecc,inc,W,w,u_l]';
    % Circular Equatorial
elseif inc == 0 && ecc == 0
    if rvec(2) > 0
        lam = acosd(rvec(1)/r);
    else
        lam = 360 - acosd(rvec(1)/r);
    end
    KepElems = [a,ecc,inc,W,w,lam]';
    % All other orbits
else
    KepElems = [a,ecc,inc,W,w,f]';
end
end