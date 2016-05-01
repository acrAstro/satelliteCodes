clear; close all; clc; asv;

mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
t0 = 0;
numPeriod = 1;
period = 2*pi/n;
tf = period*numPeriod;
numSteps = 100;
numInputs = 0;

dt = (tf-t0)/(numSteps*numPeriod);

x0 = 100;
y0 = 100;
z0 = -100;
vx0 = 0;
vy0 = -2*n*x0;
vz0 = 0;
X01 = [x0 y0 z0 vx0 vy0 vz0]';

x0 = 100;
y0 = -100;
z0 = 100;
vx0 = 0;
vy0 = -2*n*x0;
vz0 = 0;
X02 = [x0 y0 z0 vx0 vy0 vz0]';

initStruct1.params = {mu,a,X01,numInputs};
initStruct1.timeParams = {t0,dt,tf};

hcw = HCW(initStruct1);
hcw.propagateState();
X1 = hcw.X;
T1 = hcw.time;

initStruct2.params = {mu,a,X02,numInputs};
initStruct2.timeParams = {t0,dt,tf};

hcw = HCW(initStruct2);
hcw.propagateState();
X2 = hcw.X;
T2 = hcw.time;

orig = [0,0,0];

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