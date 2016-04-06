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

% Semi-major axis of chief orbit
a = 7000e3;

% Mean motion
n = sqrt(mu/a^3);

% Compute orbital period
period = 2*pi/n;

% Number of orbital periods
numPeriod = 1;
disp(numPeriod)

% Total analysis time, given in seconds
tf = numPeriod*period;

% Number of analysis time-steps
N = 100*numPeriod;

% Time Vector for ODE solver
t0 = 0;
t = linspace(t0,tf,N);

% Bounded relative motion constraint factor (dimensionless)
eccFactor = -2*n;

% Relative offsets from Chief (meters)
x0 = -300;
y0 = -300;
z0 = 100;

% Relative velocities from Chief (meters/second)
xd0 = 0;
yd0 = eccFactor*x0;
zd0 = -0.2;

% Initial condition vector for deputy satellite
XT_INIT = [x0 y0 z0 xd0 yd0 zd0]';

tic
% Hill-Clohessy-Wiltshire State Transition Matrix
X = HCW_STM(XT_INIT,n,t0,tf,N);
time1 = toc;

% Integrator options
options = odeset('RelTol',3e-12,'AbsTol',1e-15);
% Hill-Clohessy-Wiltshire equations direct numerical simulation
[T,X1] = ode113(@HCW,t,XT_INIT,options,n);
time2 = toc;

STMTIME = time1;
DNSTIME = time2-time1;
fprintf('State Transition Simulation Time (seconds) = \n');
disp(STMTIME)
fprintf('\n');
fprintf('Numerical Integration Solution Time (seconds) = \n');
disp(DNSTIME)

% Plot nice pretty picture
figure(1)
hold on
grid on
plot3(X(:,1),X(:,2),X(:,3),'k','LineWidth',2)
plot3(X1(:,1),X1(:,2),X1(:,3),'bo','LineWidth',2)
plot3(0,0,0,'r.','MarkerSize',26)
plot3(x0,y0,z0,'m.','MarkerSize',26)
leg1 = legend('STM','DNS','Chief','$X_0$','Location','Best');
title1 = title('HCW State Transition and Direct Numerical Simulation');
xl = xlabel('Radial, $x$, m');
yl = ylabel('In-Track, $y$, m');
zl = zlabel('Cross-Track, $z$, m');
set([title1 xl yl zl leg1],'interpreter','latex','fontsize',12)
axis tight