clear; close all; clc; asv;
% This is a simple test script to debug the TwoBodyOrbit class
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%

parameterization = 'OE';

% Parameters
mu = 3.986e6;
J2 = 1082.63e-6;
Req = 6378.137;

% Initial Kepler elements
a = 6678;
ecc = 0.01;
inc = 48.5*pi/180;
raan = 45*pi/180;
argPer = 45*pi/180;
f0 = 0*pi/180;
kepElems = [a ecc inc raan argPer f0]';
safetyAltitude = 75;
% Total time of flight
numPeriod = 3;
t0 = 0;

% Instantiate TwoBodyOrbit class
orbit = TwoBodyOrbit(kepElems,J2,mu,Req,t0,numPeriod,safetyAltitude,parameterization);
% Compute initial conditions
orbit.setInitialConditions();
% Propagate orbit
orbit.propagateOrbit();
% Plot orbit
orbit.plotOrbit();

