clear; close all; clc; asv; addPaths();

HCW_Config_Script();

mfocp = OptimalManeuver(initStruct);
mfocp.minimizeFuel();

meocp = OptimalManeuver(initStruct);
meocp.minimizeEnergy_Constrained();

inputStruct.states.states = {mfocp.X,meocp.X};
inputStruct.states.statesq = {mfocp.Xq,meocp.Xq};
inputStruct.states.points = {mfocp.X0,mfocp.Xf};
inputStruct.controls.controls = {mfocp.U,meocp.U};
inputStruct.controls.controlsq = {mfocp.Uq,meocp.Uq};
inputStruct.times = {mfocp.T,meocp.T};
inputStruct.id = {'hcw','hcw'};
inputStruct.lines.linestates = {'k-','b-'};
inputStruct.lines.linemods = {'linewidth','linewidth'};
inputStruct.lines.linesizes = [2,2];
inputStruct.lines.linestatesq = {'r','r'};
inputStruct.lines.linemodsq = {'linewidth','linewidth'};
inputStruct.lines.linesizesq = [2,2];
inputStruct.lines.points = {'k.','b.'};
inputStruct.lines.pointmods = {'markersize','markersize'};
inputStruct.lines.pointsizes = [25,25];
inputStruct.legends = {'HCW Min Fuel','HCW Min Energy','$X_0$','$X_f$','Min Fuel Thrust','Min Energy Thrust'};
inputStruct.title = 'Relative Trajectory';
inputStruct.labels = {'X, m','Y, m','Z, m'};
inputStruct.bounds = 'tight';

plotMotion = OrbitPlotter(inputStruct);
plotMotion.plot3DOrbit();