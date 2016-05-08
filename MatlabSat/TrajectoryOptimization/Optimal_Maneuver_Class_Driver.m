clear; close all; clc; asv; addPaths();

HCW_Config_Script();

mfocp = OptimalManeuver(initStruct);
mfocp.minimizeFuel();

meocp = OptimalManeuver(initStruct);
meocp.minimizeEnergy_Constrained();

inputStruct.states = {mfocp.X,meocp.X,mfocp.X0,mfocp.Xf,mfocp.Xq,meocp.Xq};
inputStruct.controls = {mfocp.U,meocp.U,[],[],mfocp.Uq,meocp.Uq};
inputStruct.times = {mfocp.T,meocp.T,[],[],mfocp.T,meocp.T};
inputStruct.id = {'hcw','hcw','0','0','q','q'};
inputStruct.lines = {'k-','b-','k.','b.','r-','r-'};
inputStruct.lineMods = {'linewidth','linewidth','markersize','markersize','linewidth','linewidth'};
inputStruct.lineSizes = [2,2,25,25,2,2];
inputStruct.legends = {'HCW Min Fuel','HCW Min Energy','$X_0$','$X_f$','Min Fuel Thrust','Min Energy Thrust'};
inputStruct.title = 'Relative Trajectory';
inputStruct.labels = {'X, m','Y, m','Z, m'};
inputStruct.bounds = 'tight';

plotMotion = OrbitPlotter(inputStruct);
plotMotion.plot3DOrbit();