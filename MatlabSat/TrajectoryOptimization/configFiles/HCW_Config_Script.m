mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
t0 = 0;
numPeriod = 1;
period = 2*pi/n;
dt = 1;
tf = 300;
numSteps = 100;
B = [zeros(3,3); eye(3)];
% B = [zeros(4,2); eye(2)];
samples = [];

umax = 0.075;
umin = -0.075;

eccFactor = -2*n;

x0 = 100;
y0 = 100;
z0 = -100;
vx0 = 0;
vy0 = eccFactor*x0;
vz0 = 0;
X0 = [x0 y0 z0 vx0 vy0 vz0]';

xf = 0;
yf = 0;
zf = 0;
vxf = 0;
vyf = 0;
vzf = 0;
Xf = [xf yf zf vxf vyf vzf]';

initStruct.descriptor = 'HCW';
initStruct.params = {mu,a,X0};
initStruct.terminalCondition = Xf;
initStruct.timeParams.t0 = t0;
initStruct.timeParams.dt = dt;
initStruct.timeParams.tf = tf;
initStruct.maneuverParams = {samples,B,umax,umin};