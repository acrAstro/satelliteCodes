clear; close all; clc; asv;

mu = 3.986e14;
a = 7000e3; n = sqrt(mu/a^3);
nu = 2;
mass = 1;
Lb = -0.03;
Ub = 0.03;
initStruct.params = {mu,a,nu,mass,Ub,Lb};

t0 = 0;
dt = 1;
tf = 800;
initStruct.timeParams = {t0,dt,tf};

x0 = 0;
y0 = 0;
z0 = -20;
vx0 = 0;
vy0 = 0;
vz0 = 0;

xf = 100;
yf = 0;
zf = 20;
vxf = 0;
vyf = -2*n*xf;
vzf = 0;

X0 = [x0 y0 z0 vx0 vy0 vz0]';  initStruct.X0 = X0;
Xf = [xf yf zf vxf vyf vzf]';  initStruct.Xf = Xf;

hcw = hcwOpt(initStruct);
hcw.energyOptimalTransfer();
hcw.plotTransfer();