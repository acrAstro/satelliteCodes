clear all; close all; clc;
%% Constants to be used for this program

Req     = 6378.1363e3; % Radius of Earth (meters)
mu      = 3.986004415e14; % Gravitational parameter (m^3/s^2)
J2      = 1082.629e-6; % J2 coefficient
%J2      = 0; % Set J2 to zero if you don't want the perturbed STM
tol     = 1e-12; % tolerance for transcendental root finding

%% Chief Initial Conditions
% Here you can specify classical orbital elements, nonsingular elements or
% Cartesian variables (not yet though)

% Classical orbital elements
a = 45000e3;
ecc = 0.8;
inc = 45*pi/180;
raan = pi/4;
w = pi/6;
M0 = 0;

% Check if orbit would intersect Earth
isInsideEarth(a,ecc,Req,0);

kepElemsInit = [a ecc inc raan w M0]';
eta = sqrt(1 - ecc^2);
ChiefElemsNSMean = COE_to_Nonsingular(kepElemsInit,tol);
ChiefElemsNSMean = ChiefElemsNSMean';
[DJ2mat,ChiefOsc] = MeanToOsculatingElements(J2,ChiefElemsNSMean,Req,mu);
[R0,V0] = oe2rv(kepElemsInit,mu,'r');

% Nonsingular orbital elements- mean
% a = 8494.549e3;
% th = 170.003*pi/180;
% inc = 69.988*pi/180;
% q1 = 9.420e-2;
% q2 = 3.407e-2;
% raan = 45.006*pi/180;

% ChiefElemsNSMean = [a th inc q1 q2 raan]';
% [DJ2mat,ChiefOsc] = MeanToOsculatingElements(J2,ChiefElemsNSMean,Req,mu);

%% Time
n = sqrt(mu/a^3);
t0 = 0;
period = 2*pi/n;
numPeriod = 3;
tf = numPeriod*period;
% tf = 800;
N = 301;
t = linspace(t0,tf,N*numPeriod);
%% Deputy Initial Conditions
% Here you can specify whether you want initial orbital element differences
% or XYZ differences

% Cartesian
eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));
x0 = 4000;
y0 = 4000;
z0 = 10000;
xd0 = 0;
yd0 = eccFactor*x0;
zd0 = 0;
X0GA = [x0; xd0; y0; yd0; z0; zd0];
X0 = [x0; y0; z0; xd0; yd0; zd0];

% Elements- COEs
% da = 103.624;
% de = 0.0001;
% di = 0.2*pi/180;
% dO = 0.2*pi/180;
% dw = 0;
% dM = 0.01*pi/180;
% 
% DepElemsInit = kepElemsInit + [da; de; di; dO; dw; dM];
% DepElemsInitNS =  COE_to_Nonsingular(DepElemsInit,tol);
% [DJ2Dep,DepOsc] = MeanToOsculatingElements(J2,DepElemsInitNS,Req,mu);
% deltaElems = DepOsc - ChiefOsc;
% X0 = SigmaMatrix(J2,ChiefOsc,Req,mu)*deltaElems;
% x0 = X0(1); y0 = X0(3); z0 = X0(5);

% Elements- Nonsingular
% da = -103.624;
% dth = -1.104e-3;
% di = 7.7076e-4;
% dq1 = 4.262e-5;
% dq2 = -9.708e-6;
% dO = 3.227e-3;

% DepElemsInitNS = ChiefElemsNSMean + [da; dth; di; dq1; dq2; dO];
% [DJ2Dep,DepOsc] = MeanToOsculatingElements(J2,DepElemsInitNS,Req,mu);
% deltaElems = DepOsc - ChiefOsc;
% X0 = SigmaMatrix(J2,ChiefOsc,Req,mu)*deltaElems;
% x0 = X0(1); y0 = X0(3); z0 = X0(5);

%% Compute the STM

Xt = zeros(6,length(t));
PhiJ2 = PHI_GA_STM(t,J2,ChiefOsc,ChiefElemsNSMean,Req,mu,tol);

for ii = 1:length(t)
    Xt(:,ii) = PhiJ2(:,:,ii)*X0GA;
end
options = odeset('RelTol',1e-12,'AbsTol',1e-15);
[t,Xl] = ode45(@LinearFormationFlyingEquations,t,X0,options,mu,a,ecc);

Xt = Xt';

x = Xt(:,1);
y = Xt(:,3);
z = Xt(:,5);

figure(1)
hold on
grid on
plot3(x,y,z,'k','LineWidth',2)
plot3(Xl(:,1),Xl(:,2),Xl(:,3),'ro','MarkerSize',6)
plot3(x0,y0,z0,'r.','MarkerSize',30)
plot3(x(end),y(end),z(end),'b.','MarkerSize',30)
leg1 = legend('GA-STM','LERM','$\textbf{X}_0$','$\textbf{X}_f$','Location','Best');
xl = xlabel('Radial, $x$, m');
yl = ylabel('In-track, $y$, m');
zl = zlabel('Cross-track, $z$, m');
title1 = title('Gim-Alfriend State Transition Matrix vs LERM');
set([leg1 xl yl zl],'interpreter','latex','fontsize',12)
set(title1,'interpreter','latex','fontsize',14)
axis tight

figure(4)
subplot(311)
hold on
grid on
plot(n/(2*pi).*t,x - Xl(:,1),'k','LineWidth',2)
axis tight
title1 = title('Relative Error Between GA-STM and LERM');
yl1 = ylabel('Radial, $e_x$, m');
% leg1 = legend('Deputy','Location','Best');
subplot(312)
hold on
grid on
plot(n/(2*pi).*t,y - Xl(:,2),'k','LineWidth',2)
axis tight
yl2 = ylabel('In-track, $e_y$, m');
subplot(313)
hold on
grid on
plot(n/(2*pi).*t,z - Xl(:,3),'k','LineWidth',2)
axis tight
xl = xlabel('Time, $n$-Orbits');
yl3 = ylabel('Cross-track, $e_z$, m');
set([xl yl1 yl2 yl3],'interpreter','latex','fontsize',10)
set(title1,'interpreter','latex','fontsize',14)
