clear; close all; clc; asv;

mu = 3.986e6;
J2 = 1082.63e-6;
Req = 6378.137;

a = 6678;
ecc = 0.01;
inc = 28.5*pi/180;
raan = 45*pi/180;
argPer = 45*pi/180;
f0 = 0*pi/180;

kepElems = [a ecc inc raan argPer f0]';
n = sqrt(mu/a^3);
period = 2*pi/n;
numPeriod = 5;
t0 = 0;
tf = numPeriod*period;

t = linspace(t0,tf,100*numPeriod);

orbit = TwoBodyOrbit(kepElems,J2,mu,Req,t0,numPeriod);
orbit.setInitialConditions();
orbit.propagateOrbit();
orbit.plotOrbit();

