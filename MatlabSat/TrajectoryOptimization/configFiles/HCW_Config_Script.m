mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
t0 = 0;
numPeriod = 1;
period = 2*pi/n;
tf = period*numPeriod;
numSteps = 100;
B = [zeros(3,3); eye(3)];
samples = [];

umax = 0.1;

dt = (tf-t0)/(numSteps*numPeriod);

x0 = 100;
y0 = 100;
z0 = -100;
vx0 = 0;
vy0 = -2*n*x0;
vz0 = 0;
X0 = [x0 y0 z0 vx0 vy0 vz0]';

initStruct.descriptor = 'HCW';
initStruct.params = {mu,a,X0};
initStruct.timeParams.t0 = t0;
initStruct.timeParams.dt = dt;
initStruct.timeParams.tf = tf;
initStruct.maneuverParams = {samples,B,umax};