clear; close all; clc; asv; addPaths();

HCW_Config_Script();

mt = ConvexSFFManeuver(initStruct);
mt.findFeasibleTime(536,1);

mf = ConvexSFFManeuver(initStruct);
mf.minimizeFuel();

mec = ConvexSFFManeuver(initStruct);
mec.minimizeEnergy_Constrained();

meu = ConvexSFFManeuver(initStruct);
meu.minimizeEnergy();

inputStruct.states.states = {mt.X,mf.X,mec.X,meu.X,mf.X0,mf.Xf};
inputStruct.states.statesq = {mt.Xq,mf.Xq,mec.Xq,meu.Xq};
inputStruct.controls.controls = {mt.U,mf.U,mec.U,meu.U};
inputStruct.controls.controlsq = {mt.Uq,mf.Uq,mec.Uq,meu.Uq};
inputStruct.times = {mt.T,mf.T,mec.T,meu.T};
inputStruct.id = {'hcw','hcw','hcw','hcw','.','.'};
inputStruct.lines.linestates = {'r--','k-','b-','g-','k.','b.'};
inputStruct.lines.linemods = {'linewidth','linewidth','linewidth','linewidth','markersize','markersize'};
inputStruct.lines.linesizes = [2,2,2,2,25,25];
inputStruct.lines.linestatesq = {'r','r','r','r'};
inputStruct.lines.linemodsq = {'linewidth','linewidth','linewidth','linewidth'};
inputStruct.lines.linesizesq = [2,2,2,2];
inputStruct.legends = {'Min Time','Min Fuel','Constrained Min Energy','Unconstrained Min Energy','$X_0$','$X_f$','Thrust'};
inputStruct.title = 'Relative Trajectory';
inputStruct.labels = {'X, m','Y, m','Z, m'};
inputStruct.bounds = 'tight';
inputStruct.shuttleFlag = 'no';

plotMotion = OrbitPlotter(inputStruct);

figure
subplot(211)
hold on
grid on
plot(mt.T(1:size(mt.U,2)),mt.U(1,:),'r','linewidth',2);
plot(mf.T,mf.U(1,:),'k','LineWidth',2);
plot(mec.T,mec.U(1,:),'b','Linewidth',2);
plot(meu.T,meu.U(1,:),'g','Linewidth',2);
plot([mf.T(1),mf.T(end)],[mf.umax,mf.umax],'k--','lineWidth',2)
plot([mf.T(1),mf.T(end)],[-mf.umax,-mf.umax],'k--','lineWidth',2)
axis([mf.T(1),mf.T(end),-2*umax,2*umax])
subplot(212)
hold on
grid on
plot(mt.T(1:size(mt.U,2)),mt.U(2,:),'r','linewidth',2);
plot(mf.T,mf.U(2,:),'k','LineWidth',2);
plot(mec.T,mec.U(2,:),'b','Linewidth',2);
plot(meu.T,meu.U(2,:),'g','Linewidth',2);
plot([mf.T(1),mf.T(end)],[mf.umax,mf.umax],'k--','lineWidth',2)
plot([mf.T(1),mf.T(end)],[-mf.umax,-mf.umax],'k--','lineWidth',2)
axis([mf.T(1),mf.T(end),-2*umax,2*umax])