clear; close all; clc; asv; addPaths();

HCW_Config_Script();

mfocp = OptimalManeuver(initStruct);
mfocp.minimizeFuel();

meocp = OptimalManeuver(initStruct);
meocp.minimizeEnergy_Constrained();

inputStruct.states.states = {mfocp.X,meocp.X,mfocp.X0,mfocp.Xf};
inputStruct.states.statesq = {mfocp.Xq,meocp.Xq};
inputStruct.controls.controls = {mfocp.U,meocp.U};
inputStruct.controls.controlsq = {mfocp.Uq,meocp.Uq};
inputStruct.times = {mfocp.T,meocp.T};
inputStruct.id = {'hcw','hcw','.','.'};
inputStruct.lines.linestates = {'k-','b-','k.','b.'};
inputStruct.lines.linemods = {'linewidth','linewidth','markersize','markersize'};
inputStruct.lines.linesizes = [2,2,25,25];
inputStruct.lines.linestatesq = {'r','r'};
inputStruct.lines.linemodsq = {'linewidth','linewidth'};
inputStruct.lines.linesizesq = [2,2];
inputStruct.legends = {'HCW Min Fuel','HCW Min Energy','$X_0$','$X_f$','Min Fuel Thrust','Min Energy Thrust'};
inputStruct.title = 'Relative Trajectory';
inputStruct.labels = {'X, m','Y, m','Z, m'};
inputStruct.bounds = 'tight';
inputStruct.shuttleFlag = 'no';

plotMotion = OrbitPlotter(inputStruct);
plotMotion.plot3DOrbit();

figure
subplot(211)
hold on
grid on
plot(mfocp.T,mfocp.U(1,:),'k','LineWidth',2);
plot(meocp.T,meocp.U(1,:),'r','Linewidth',2);
plot([mfocp.T(1),mfocp.T(end)],[mfocp.umax,mfocp.umax],'k--','lineWidth',2)
plot([mfocp.T(1),mfocp.T(end)],[-mfocp.umax,-mfocp.umax],'k--','lineWidth',2)
axis([mfocp.T(1),mfocp.T(end),-2*umax,2*umax])
subplot(212)
hold on
grid on
plot(mfocp.T,mfocp.U(2,:),'k','LineWidth',2);
plot(meocp.T,meocp.U(2,:),'r','Linewidth',2);
plot([mfocp.T(1),mfocp.T(end)],[mfocp.umax,mfocp.umax],'k--','lineWidth',2)
plot([mfocp.T(1),mfocp.T(end)],[-mfocp.umax,-mfocp.umax],'k--','lineWidth',2)
axis([mfocp.T(1),mfocp.T(end),-2*umax,2*umax])