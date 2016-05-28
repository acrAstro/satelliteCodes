function Phi_Broucke = BrouckeSTM(tspan,mu,a,ecc,omega)

n = sqrt(mu/a^3);
eta = sqrt(1-ecc^2);
p = a*(1-ecc^2);
t0 = tspan(1);
Phi_Broucke = zeros(6,6,length(tspan));
for ii = 1:length(tspan)
    dt = tspan(ii)-t0;
    [r,~,theta,~,~,EccAnom] = polarParams(dt,mu,a,ecc,omega);
    f = theta - omega;
    
    R = Rmatrix(dt,mu,a,ecc,r,f,EccAnom,n,eta,p);
    RI = RInverseMatrix(dt,mu,a,ecc,r,f,EccAnom,n,eta,p);
    Phi_Broucke(:,:,ii) = R*RI;
    
end
end

function R = Rmatrix(dt,mu,a,ecc,r,f,EccAnom,n,eta,p)
R = zeros(6,6);
Rxy(:,1) = [r/a - 3*n*dt*ecc*sin(f)/(2*eta);
    -3/2*n*dt*eta*(a/r);
    -n*ecc*sin(f)/(2*eta) - 3/2*dt*ecc*cos(f)*(n*a/r)^2;
    3/2*dt*ecc*sin(f)*(n*a/r)^2 - 3/2*eta*(n*a/r)];

Rxy(:,2) = [-cos(f);
    (1+r/p)*sin(f);
    n*sin(f)*eta*(a/r)^2
    n*eta*(1+r/p)*(a/r)^2*cos(f) + ecc*n*sin(f)^2/sqrt((1-ecc)^3)];

%     Rxy(:,3) = [ecc*sin(f)/eta;
%         eta*(a/r);
%         ecc*n*cos(f)*(a/r)^2;
%         -ecc*n*sin(f)*(a/r)^2];
Rxy(:,3) = [sin(f);
    r/p*cos(f)*(2+ecc*cos(f));
    sqrt(mu*p)/r^2*cos(f);
    -sqrt(mu*p/r^2)*ecc*sin(f)*(1+r^2/p^2)];

Rxy(:,4) = [0;
    r/a;
    0;
    ecc*n*sin(f)/eta];

z1 = sqrt(a*n/mu)*r*cos(f); z2 = sqrt(a*n/(mu*(1-ecc^2))*r*sin(f));
z1Dot = -a*sqrt(n)*sin(EccAnom)/r; z2Dot = a*sqrt(n)*cos(EccAnom)/r;

Rz = [z1, z2
    z1Dot, z2Dot];
R(1:4,1:4) = Rxy; R(5:6,5:6) = Rz;
end

function RI = RInverseMatrix(dt,mu,a,ecc,r,f,EccAnom,n,eta,p)
RI = zeros(6,6);
Imatrix = [2*(a/r)^2*(1+p/r), (a/r)*(3*cos(f) + ecc*(2+cos(f)^2));
    -2*(a/r)^2*ecc*sin(f), -ecc*sin(f)*a/p*(ecc+cos(f)*(2+ecc*cos(f)));
    2*ecc*sin(f)/(n*eta), eta/n*sin(f);
    2*(a/r)*eta/n, eta/n*(cos(f)+r/p*(ecc+cos(f)))];

Jmatrix = [-a/p*sin(f)*(ecc^2+(1+ecc*cos(f))*(3+ecc*cos(f))), -a^2/(p*r)*ecc*sin(f)*(2+p/r+r/p);
    (a/p)*ecc*sin(f)^2*(2+ecc*cos(f)), (a/p)^2*(1+ecc*cos(f) - ecc^2*(ecc*cos(f)^3+2*cos(f)^2-1));
    r*eta/(p*n)*(cos(f)+ecc*(cos(f)^2-2)), (ecc*cos(f)-1)/(n*eta)*(1+r/p);
    -eta/n*sin(f)*(1+r/p), -ecc*sin(f)/(n*eta)*(1+r/p)];

Kmatrix = [n*ecc/eta*(a/r)*(1+p/r), n/sqrt((1-ecc^2)^3)*(a/r)^2*(1+p/r);
    -n*ecc^2*sin(f)/eta*(a/r)^2, -n*ecc*sin(f)/sqrt((1-ecc^2)^3)*(a/r)^2;
    ecc^2*sin(f)/(1-ecc^2), ecc*sin(f)/((1-ecc^2)^2);
    a*ecc/r, a/((1-ecc^2)*r)];


z1 = sqrt(a*n/mu)*r*cos(f); z2 = sqrt(a*n/(mu*(1-ecc^2))*r*sin(f));
z1Dot = -a*sqrt(n)*sin(EccAnom)/r; z2Dot = a*sqrt(n)*cos(EccAnom)/r;

RIxy = transpose([Imatrix, Jmatrix + 3*dt*Kmatrix]);
RIz = [z2Dot, -z2;
    -z1Dot, z1];

RI(1:4,1:4) = RIxy; RI(5:6,5:6) = RIz;

end

function [r,rdot,theta,thetadot,thetaddot,EccAnom] = polarParams(dt,mu,a,ecc,omega)

% t0 = tspan(1);
% t = tspan(2);
% dt = t-t0;
n = sqrt(mu/a^3);
MeanAnom = n*dt;
[EccAnom,~] = kepler(MeanAnom,ecc);
f = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(EccAnom/2));
theta = omega + f;

r = a*(1-ecc^2)/(1+ecc*cos(theta));
rdot = sqrt(mu*a*(1-ecc^2))*(1+ecc*cos(theta)^2)/(a^2*(1-ecc^2)^2);
thetadot = ecc*mu*sin(theta)/(sqrt(mu*a*(1-ecc^2)));
thetaddot = 2*mu*ecc*(1+ecc*cos(theta))^3*sin(theta)/(a^3*(1-ecc^2)^3);
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
tol = 1e-10;
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