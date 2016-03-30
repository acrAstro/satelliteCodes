function DX = orbitEquation(t,x,mu,J2,Req)
% Main function to integrate the orbital equation of motion with the
% influence of J2
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%
% State vector
r = x(1:3,1);
v = x(4:6,1);
R = norm(r); % Magnitude of R

% Return J2 acceleration
aj2 = J2Vector(mu,r,J2,Req);

% Differential equation in state space form
rdot = v;
vdot = -mu/R^3.*r + aj2;
DX = [rdot(:); vdot(:)];
end

function aj2 = J2Vector(mu,r,J2,req)
% J2 vector
x = r(1); y = r(2); z = r(3);
R = norm(r);
gamma = -3/2*J2*(mu/R^2)*(req/R)^2;
j2x = (1 - 5*(z/R)^2)*x/R;
j2y = (1 - 5*(z/R)^2)*y/R;
j2z = (3 - 5*(z/R)^2)*z/R;

aj2 = gamma.*[j2x; j2y; j2z];
end