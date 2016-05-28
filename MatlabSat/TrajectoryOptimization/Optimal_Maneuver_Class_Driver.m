clear; close all; clc; asv; addPaths();

HCW_Config_Script();

mf = ConvexSFFManeuver(initStruct);
mf.minimizeFuel();

meu = ConvexSFFManeuver(initStruct);
meu.minimizeEnergy();

mec = ConvexSFFManeuver(initStruct);
mec.minimizeEnergy_Constrained();

inputStruct.states.states = {mf.X,meu.X,mec.X,mf.X0,mf.Xf};
inputStruct.states.statesq = {mf.Xq,meu.Xq,mec.Xq};
inputStruct.controls.controls = {mf.U,meu.U,mec.U};
inputStruct.controls.controlsq = {mf.Uq,meu.Uq,mec.Uq};
inputStruct.times = {mf.T,meu.T,mec.T};
inputStruct.id = {'hcw','hcw','hcw','.','.'};
inputStruct.lines.linestates = {'k-','b-','g-','k.','b.'};
inputStruct.lines.linemods = {'linewidth','linewidth','linewidth','markersize','markersize'};
inputStruct.lines.linesizes = [2,2,2,25,25];
inputStruct.lines.linestatesq = {'r','r','r'};
inputStruct.lines.linemodsq = {'linewidth','linewidth','linewidth'};
inputStruct.lines.linesizesq = [2,2,2];
inputStruct.legends = {'HCW Min Fuel','HCW Min Energy','HCW Min Energ Con','$X_0$','$X_f$','Min Fuel Thrust','Min Energy Thrust'};
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
plot(mf.T,mf.U(1,:),'k','LineWidth',2);
plot(meu.T,meu.U(1,:),'b','Linewidth',2);
plot(mec.T,mec.U(1,:),'g','Linewidth',2);
plot([mf.T(1),mf.T(end)],[mf.umax,mf.umax],'k--','lineWidth',2)
plot([mf.T(1),mf.T(end)],[-mf.umax,-mf.umax],'k--','lineWidth',2)
axis([mf.T(1),mf.T(end),-2*umax,2*umax])
subplot(212)
hold on
grid on
plot(mf.T,mf.U(2,:),'k','LineWidth',2);
plot(meu.T,meu.U(2,:),'b','Linewidth',2);
plot(mec.T,mec.U(2,:),'g','Linewidth',2);
plot([mf.T(1),mf.T(end)],[mf.umax,mf.umax],'k--','lineWidth',2)
plot([mf.T(1),mf.T(end)],[-mf.umax,-mf.umax],'k--','lineWidth',2)
axis([mf.T(1),mf.T(end),-2*umax,2*umax])