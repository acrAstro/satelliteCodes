clear; close all; clc;
% This is a simple test script to debug the TwoBodyOrbit class
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%

parameterization = 'OE';

% Parameters
mu = 3.986e5;
J2 = 1082.63e-6;
Req = 6378.137;
mass = 15;

% Initial Kepler elements
a = 6878;
ecc = 0.01;
inc = 50.4*pi/180;
raan = 45*pi/180;
argPer = -30*pi/180;
f0 = 0*pi/180;
kepElems = [a ecc inc raan argPer f0]';
safetyAltitude = 75;
% Total time of flight
t0 = 0;
dt = 1;
tf = 800;

% Initialization structure for orbit class
chiefStruct.kepElems         = kepElems;
chiefStruct.params           = {J2,mu,Req,safetyAltitude,mass};
chiefStruct.timeParams       = {t0,dt,tf};
chiefStruct.Parameterization = parameterization;

% Instantiate TwoBodyOrbit class
orbit = TwoBodyOrbit(chiefStruct);
% Propagate orbit
orbit.propagateOrbit();
% Plot orbit
orbit.plotOrbit();
% % Plot orbital elements
orbit.plotOrbitalElements();

