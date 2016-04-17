clear all; close all; clc;

mu = 3.986e14;
a = 7000e3;
ecc = 0.01;
inc = 28*pi/180;
RAAN = 45*pi/180;
ArgPer = 0*pi/180;
f0 = 0*pi/180;
n = sqrt(mu/a^3);
eccFactor = -n*(2+ecc)/(sqrt((1+ecc)*(1-ecc)^3));

x0 = 100;
y0 = 100;
z0 = 100;
vx0 = 0;
vy0 = eccFactor*x0;
vz0 = 0;

X0 = [x0 y0 vx0 vy0 z0 vz0]';
t0 = 0;
period = 2*pi/n;
numPeriod = 1;
tf = numPeriod*period;
tspan = linspace(t0,tf,numPeriod*150);
Phi = BrouckeSTM(tspan,mu,a,ecc,ArgPer);

for ii = 1:length(tspan)
    X(:,ii+1) = Phi(:,:,ii)*X0;
end

figure
hold on
grid on
plot3(X(1,:),X(2,:),X(5,:),'k','linewidth',2)
axis tight
