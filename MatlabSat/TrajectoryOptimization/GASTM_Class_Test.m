clear; close all; clc; asv;

Req     = 6378.1363e3; % Radius of Earth (meters)
mu      = 3.986004415e14; % Gravitational parameter (m^3/s^2)
J2      = 1082.629e-6; % J2 coefficient
tol     = 1e-12; % tolerance for transcendental root finding
safetyAltitude = 50e3;
samples = 3;
B = [zeros(3,3); eye(3)];

% Valid descriptions are 'Classical'; 'Nonsingular'
chiefOrbitDescription = 'Nonsingular';
% Valid descriptions are 'Cartesian'; 'Relative Classical'; 'Relative Nonsingular'
deputyOrbitDescription = 'Relative Classical';

[chiefOrbitDescription,deputyOrbitDescription] = checkDescriptors(chiefOrbitDescription,deputyOrbitDescription);

% Note: Classical with Cartesian is the most stable currently

method = chiefOrbitDescription;
switch method
    case 'Classical'
        a = 6678e3;
        ecc = 0.01;
        inc = 45*pi/180;
        raan = pi/4;
        w = pi/6;
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

method = deputyOrbitDescription;
switch method
    case 'Cartesian'
        eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));
        x0 = 4000;
        y0 = 4000;
        z0 = 1000;
        xd0 = 0;
        yd0 = eccFactor*x0;
        zd0 = 0;
        RelInitState = [x0 xd0 y0 yd0 z0 zd0]';
    case 'Relative Classical'
        da = 0;
        de = 0;
        di = 0.2*pi/180;
        dO = 0*pi/180;
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

t0 = 0; numPeriod = 3; numSteps = 100;

initStruct.params = {Req,mu,J2,tol,t0,numPeriod,safetyAltitude,numSteps,samples,B};
initStruct.initChiefDescription = chiefOrbitDescription;
initStruct.initDeputyDescription = deputyOrbitDescription;
initStruct.RelInitState = RelInitState(:);
initStruct.Elements = Elements(:);

GA = GimAlfriendSTM(initStruct);
GA.propagateState();
GA.plotOrbit();
GA.makeDiscreteMatrices();