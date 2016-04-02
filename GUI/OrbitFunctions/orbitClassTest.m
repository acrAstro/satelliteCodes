clear; close all; clc; asv;
% This is a simple test script to debug the TwoBodyOrbit class
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%

parameterization = 'RV';

% Parameters
mu = 3.986e5;
J2 = 1082.63e-6;
Req = 6378.137;

% Initial Kepler elements
a = 45000;
ecc = 0.8;
inc = 68.4*pi/180;
raan = 45*pi/180;
argPer = 30*pi/180;
f0 = 0*pi/180;
kepElems = [a ecc inc raan argPer f0]';
safetyAltitude = 75;
% Total time of flight
numPeriod = 1;
t0 = 0;

% Initialization structure for orbit class
initStruct.kepElems         = kepElems;
initStruct.params           = {J2,mu,Req,t0,numPeriod,safetyAltitude};
initStruct.Parameterization = parameterization;

% Instantiate TwoBodyOrbit class
orbit = TwoBodyOrbit(initStruct);
% Compute initial conditions
orbit.setInitialConditions();
% Propagate orbit
orbit.propagateOrbit();
% Plot orbit
orbit.plotOrbit();
% % Plot orbital elements
orbit.plotOrbitalElements();

