clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All units are in meters and meters/sec, instead of mu = 3.986e5 km^3/s^3,
% it is 3.986e14 m^3/s^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is provided 'as-is' for use with AOE 5234: Orbital Mechanics.
% Numerical accuracy is guaranteed only to within the bounds specified by
% The MathWorks Inc.
%
% Author: Andrew Rogers, 10 February 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravitational parameter
mu = 3.986e14;

%% Set up Chief orbit with orbital elements
caseNum   = 'Case 2';
if strcmp(caseNum,'Case 1') == 1
    a      = 7000e3;
    ecc    = 0.0;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
elseif strcmp(caseNum,'Case 2') == 1
    a      = 15000e3;
    ecc    = 0.4;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
elseif strcmp(caseNum,'Case 3') == 1
    a      = 45000e3;
    ecc    = 0.8;
    inc    = 45*pi/180;
    raan   = pi/6;
    argper = pi/6;
    f0     = 0;
else
    a      = 6678e3;
    ecc    = 0.025;
    inc    = 28.5*pi/180;
    raan   = 0;
    argper = 0;
    f0     = 0;
end

n = sqrt(mu/a^3);

%% Set up deputy initial conditions
% Bounded relative motion constraint factor (dimensionless)
eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));

% Relative offsets from Chief (meters)
x0 = -300;
y0 = -300;
z0 = 0;

% Relative velocities from Chief (meters/second)
xd0 = 0;
yd0 = eccFactor*x0;
zd0 = 0;

% Initial condition vector for deputy satellite
XT_INIT = [x0 y0 z0 xd0 yd0 zd0]';

% Initial time
t0 = 0;
% Number of orbital periods
numPeriod = 1;
% Integrator options
options = odeset('RelTol',3e-12,'AbsTol',1e-15,'Stats','on');

initStruct.params{1} = mu;
initStruct.params{2} = [a ecc inc raan argper f0]';
initStruct.params{3} = XT_INIT;
initStruct.timeParams{1} = t0;
initStruct.timeParams{2} = numPeriod;
initStruct.options = options;

States = {RelativeMotion(initStruct),RelativeMotion(initStruct)};
States{1}.propagateLinear();
States{2}.propagateNonlinear();

figure
hold on
grid on
plot3(States{1}.XL(:,1),States{1}.XL(:,2),States{1}.XL(:,3),'k','linewidth',2)
plot3(States{2}.XN(:,1),States{2}.XN(:,2),States{2}.XN(:,3),'ro','markersize',6)
axis tight
