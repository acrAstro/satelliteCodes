Req     = 6378.1363e3; % Radius of Earth (meters)
mu      = 3.986004415e14; % Gravitational parameter (m^3/s^2)
J2      = 1082.629e-6; % J2 coefficient
tol     = 1e-13; % tolerance for transcendental root finding
safetyAltitude = 50e3;
samples = 3;
B = [zeros(3,3); eye(3)];
umax = 0.1;
umin = 0;

% Valid descriptions are 'Classical'; 'Nonsingular'
chiefOrbitDescription = 'Classical';
% Valid descriptions are 'Cartesian'; 'Relative Classical'; 'Relative Nonsingular'
deputyOrbitDescriptionInit = 'Cartesian';
% Valid descriptions are 'Cartesian'; 'Relative Classical'; 'Relative
% Nonsingular';
deputyOrbitDescriptionFinal = 'Cartesian';

[chiefOrbitDescription,deputyOrbitDescriptionInit,deputyOrbitDescriptionFinal]...
    = checkDescriptors(chiefOrbitDescription,deputyOrbitDescriptionInit,deputyOrbitDescriptionFinal);

% Note: Classical with Cartesian is the most stable currently

method = chiefOrbitDescription;
switch method
    case 'Classical'
        a = 6678e3;
        ecc = 0.0;
        inc = 28*pi/180;
        raan = pi/4;
        w = 0;
        M0 = 0;
        n = sqrt(mu/a^3);
        Elements = [a ecc inc raan w M0]';
    case 'Nonsingular'
         a = 8494.549e3;
         th = 170.003*pi/180;
         inc = 69.988*pi/180;
         q1 = 9.420e-2;
         q2 = 3.407e-2;
         raan = 45.006*pi/180;
         n = sqrt(mu/a^3);
         Elements = [a th inc q1 q2 raan]';
         ecc = sqrt(q1^2 + q2^2);         
end
period = 2*pi/n;
method = deputyOrbitDescriptionInit;
switch method
    case 'Cartesian'
        eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));
        x0 = 100;
        y0 = 100;
        z0 = -100;
        xd0 = 0;
        yd0 = eccFactor*x0;
        zd0 = 0;
        RelInitState = [x0 xd0 y0 yd0 z0 zd0]';
    case 'Relative Classical'
        da = 0;
        de = 0;
        di = 0.2*pi/180;
        dO = 0.2*pi/180;
        dw = 0;
        dM = 0*pi/180;
        RelInitState = [da de di dO dw dM]';
    case 'Relative Nonsingular'
        da = -103.624;
        dth = -1.104e-3;
        di = 7.7076e-4;
        dq1 = 4.262e-5;
        dq2 = -9.708e-6;
        dO = 3.227e-3;
        RelInitState = [da dth di dq1 dq2 dO]';
end

method = deputyOrbitDescriptionFinal;
switch method
    case 'Cartesian'
        eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));
        x0 = 0;
        y0 = 0;
        z0 = 0;
        xd0 = 0;
        yd0 = eccFactor*x0;
        zd0 = 0;
        RelFinalState = [x0 xd0 y0 yd0 z0 zd0]';
    case 'Relative Classical'
        da = 0;
        de = 0;
        di = 0.2*pi/180;
        dO = 0.2*pi/180;
        dw = 0;
        dM = 0*pi/180;
        RelFinalState = [da de di dO dw dM]';
    case 'Relative Nonsingular'
        da = -103.624;
        dth = -1.104e-3;
        di = 7.7076e-4;
        dq1 = 4.262e-5;
        dq2 = -9.708e-6;
        dO = 3.227e-3;
        RelFinalState = [da dth di dq1 dq2 dO]';
end

t0 = 0; numPeriod = 3; numSteps = 100; dt = 1;
tf = 300;

initStruct.descriptor = 'GimAlfriendSTM';
initStruct.params = {Req,mu,J2,tol,safetyAltitude};
initStruct.maneuverParams = {samples,B,umax,umin};
initStruct.timeParams.t0 = t0;
initStruct.timeParams.dt = dt;
initStruct.timeParams.tf = tf;
initStruct.initChiefDescription = chiefOrbitDescription;
initStruct.initDeputyDescription = deputyOrbitDescriptionInit;
initStruct.finalDeputyDescription = deputyOrbitDescriptionFinal;
initStruct.RelInitState = RelInitState(:);
initStruct.RelFinalState = RelFinalState(:);
initStruct.Elements = Elements(:);

