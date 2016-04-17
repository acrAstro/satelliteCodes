function Sigma = SigmaMatrix(J2,elems,Re,mu)

% System_matrix, Sigma in osculating element with perturbation by J2
%
% input :
%    elems(1)= a
%    elems(2)= theta
%    elems(3)= i
%    elems(4)= q1
%    elems(5)= q2
%    elems(6)= Omega
%
% output :
%    6x6 system matrix
%    Sigma = A + (3*J2*Re^2)*B
%
gamma   = 3*J2*Re^2;
a       = elems(1);
argLat  = elems(2);
inc     = elems(3);
q1      = elems(4);
q2      = elems(5);
% RAAN   = elems(6);

% Evaluations from the inputs
% Hamiltonian = -mu/(2*a);
% n = sqrt(mu/a^3);
p = a*(1 - q1^2 - q2^2);
% h = sqrt(mu*p);
R = p/(1 + q1*cos(argLat) + q2*sin(argLat));
Vr = sqrt(mu/p)*(q1*sin(argLat) - q2*cos(argLat));
Vt = sqrt(mu/p)*(1 + q1*cos(argLat) + q2*sin(argLat));

% Transformation Matrix A
A11 = R/a;
A12 = R*Vr/Vt;
A13 = 0;
A14 = -(2*a*R*q1/p) - (R^2/p)*cos(argLat);
A15 = -(2*a*R*q2/p) - (R^2/p)*sin(argLat);
A16 = 0;

A21 = -(1/2)*Vr/a;
A22 = sqrt(mu/p)*((p/R)-1);
A23 = 0;
A24 = (Vr*a*q1/p) + sqrt(mu/p)*sin(argLat);
A25 = (Vr*a*q2/p) - sqrt(mu/p)*cos(argLat);
A26 = 0;

A31 = 0;
A32 = R;
A33 = 0;
A34 = 0;
A35 = 0;
A36 = R*cos(inc);

A41 = -(3/2)*Vt/a;
A42 = -Vr;
A43 = 0;
A44 = (3*Vt*a*q1/p) + 2*sqrt(mu/p)*cos(argLat);
A45 = (3*Vt*a*q2/p) + 2*sqrt(mu/p)*sin(argLat);
A46 = Vr*cos(inc);

A51 = 0;
A52 = 0;
A53 = R*sin(argLat);
A54 = 0;
A55 = 0;
A56 = -R*cos(argLat)*sin(inc);

A61 = 0;
A62 = 0;
A63 = Vt*cos(argLat) + Vr*sin(argLat);
A64 = 0;
A65 = 0;
A66 = (Vt*sin(argLat) - Vr*cos(argLat))*sin(inc);

A = [ A11  A12  A13  A14  A15  A16;
    A21  A22  A23  A24  A25  A26;
    A31  A32  A33  A34  A35  A36;
    A41  A42  A43  A44  A45  A46;
    A51  A52  A53  A54  A55  A56;
    A61  A62  A63  A64  A65  A66 ];


% Transformation Matrix B
B11 = 0;
B12 = 0;
B13 = 0;
B14 = 0;
B15 = 0;
B16 = 0;

B21 = 0;
B22 = 0;
B23 = 0;
B24 = 0;
B25 = 0;
B26 = 0;

B31 = 0;
B32 = 0;
B33 = 0;
B34 = 0;
B35 = 0;
B36 = 0;

B41 = 0;
B42 = 0;
B43 = -Vt*sin(inc)*cos(inc)*sin(argLat)^2/(p*R);
B44 = 0;
B45 = 0;
B46 = Vt*sin(inc)^2*cos(inc)*sin(argLat)*cos(argLat)/(p*R);

B51 = 0;
B52 = 0;
B53 = 0;
B54 = 0;
B55 = 0;
B56 = 0;

B61 = 0;
B62 = Vt*sin(inc)*cos(inc)*sin(argLat)/(p*R);
B63 = 0;
B64 = 0;
B65 = 0;
B66 = Vt*sin(inc)*cos(inc)^2*sin(argLat)/(p*R);

B = [ B11  B12  B13  B14  B15  B16;
    B21  B22  B23  B24  B25  B26;
    B31  B32  B33  B34  B35  B36;
    B41  B42  B43  B44  B45  B46;
    B51  B52  B53  B54  B55  B56;
    B61  B62  B63  B64  B65  B66 ];


% System Matrix AA_B

Sigma = A + gamma*B;

end