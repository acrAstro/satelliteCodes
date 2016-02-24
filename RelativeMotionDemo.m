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

%% Set up analysis time
% Mean motion
n = sqrt(mu/a^3);

% Initial time
t0 = 0;

% Compute orbital period
period = 2*pi/n;

% Number of orbital periods
numPeriod = 1;

% Total analysis time, given in seconds
tf = numPeriod*period;

% Number of analysis time-steps (defaults to 100 if not specified)
N = 200;

% Time vector for ODE solver
t = linspace(t0,tf,N);

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

% Integrator options
options = odeset('RelTol',3e-12,'AbsTol',1e-15,'Stats','on');

%% Numerically integrate equations of motion

% Linear Equations (LERM)
[T1,XL] = ode113(@LinearFormationFlyingEquations,t,XT_INIT,options,mu,a,ecc);

% Nonlinear Equations (NERM)
[T2,XN] = ode113(@NonlinearFormationFlyingEquations,t,XT_INIT,options,mu,a,ecc);

% Plot nice pretty picture
figure(1)
hold on
grid on
plot3(XL(:,1),XL(:,2),XL(:,3),'k','LineWidth',2)
plot3(XN(:,1),XN(:,2),XN(:,3),'bo','MarkerSize',6)
plot3(0,0,0,'r.','MarkerSize',26)
plot3(x0,y0,z0,'m.','MarkerSize',26)
leg1 = legend('Linear','Nonlinear','Chief','$X_0$','Location','Best');
title1 = title('Linear and Nonlinear Relative Motion Simulation');
xl = xlabel('Radial, $x$, m');
yl = ylabel('In-Track, $y$, m');
zl = zlabel('Cross-Track, $z$, m');
set([title1 xl yl zl leg1],'interpreter','latex','fontsize',12)
axis tight