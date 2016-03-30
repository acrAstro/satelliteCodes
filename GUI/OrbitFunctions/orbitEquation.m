function DX = orbitEquation(t,x,mu,J2,Req)

r = x(1:3,1);
v = x(4:6,1);
R = norm(r);

aj2 = J2Vector(mu,r,J2,Req);

rdot = v;
vdot = -mu/R^3.*r + aj2;
DX = [rdot(:); vdot(:)];
end

function aj2 = J2Vector(mu,r,J2,req)

x = r(1); y = r(2); z = r(3);
R = norm(r);
gamma = -3/2*J2*(mu/R^2)*(req/R)^2;
j2x = (1 - 5*(z/R)^2)*x/R;
j2y = (1 - 5*(z/R)^2)*y/R;
j2z = (3 - 5*(z/R)^2)*z/R;

aj2 = gamma.*[j2x; j2y; j2z];
end

% function ad = DragVector(Cd,mass