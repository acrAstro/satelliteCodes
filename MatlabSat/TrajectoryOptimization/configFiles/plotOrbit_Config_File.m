inputStruct.states = {X1,X2,orig};
inputStruct.times = {T1,T2,0};
inputStruct.id = {'hcw','hcw2','0'};
inputStruct.lines = {'kd-','r*-','k.'};
inputStruct.lineMods = {'linewidth','linewidth','markersize'};
inputStruct.lineSizes = [2,2,20];
inputStruct.legends = {'HCW','HCW2','Origin'};
inputStruct.title = 'Relative Motion';
inputStruct.labels = {'X, m','Y, m','Z, m'};
inputStruct.bounds = 'tight';

plotMotion = OrbitPlotter(inputStruct);
plotMotion.plot3DOrbit();