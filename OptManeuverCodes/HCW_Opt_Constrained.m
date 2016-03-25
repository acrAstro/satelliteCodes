clear all; close all; clc;

mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
nu = 6;
nx = 6;
mass = 20;
[Ac,Bc] = HCW_Matrices(n,nu,mass);

period = 2*pi/n;
dt = 1;
numPeriod = 1/4;
tf = numPeriod*period;
Nsim = round(tf/dt);
sysc = ss(Ac,Bc,[],[]);
sysd = c2d(sysc,dt,'zoh');

A = sysd.a;
B = sysd.b;

Lb = 0;
Ub = 0.5;

x0 = 0;
y0 = 0;
z0 = 0;
vx0 = 0;
vy0 = 0;
vz0 = 0;

xf = 400;
yf = 0;
zf = 0;
vxf = 0;
vyf = -2*n*xf;
vzf = 0;

X0 = [x0 y0 z0 vx0 vy0 vz0]';
Xf = [xf yf zf vxf vyf vzf]';

u = sdpvar(nu,Nsim);
uSlack = sdpvar(nu,Nsim);

X = sdpvar(nx,Nsim+1);
x1 = sdpvar(nx,1);
x2 = sdpvar(nx,1);
constraints = [];
objective = 0;
constraints = [constraints, X(:,1) == X0];

for ii = 1:Nsim
    objective = objective + sum(uSlack(1:nu,ii));
    constraints = [constraints, X(:,ii+1) == A*X(:,ii) + B*u(:,ii)];
    constraints = [constraints, u(1:nu,ii) <= uSlack(1:nu,ii)];
    constraints = [constraints, Lb.*ones(nu,1) <= u(1:nu,ii) <= Ub.*ones(nu,1)];
end

constraints = [constraints, X(:,Nsim+1) == Xf];
options = sdpsettings('solver','gurobi','saveyalmipmodel',1,'verbose',3);

inputParams = {x1,x2};
outputSolutions = {u,objective};

controller = optimizer(constraints,objective,options,inputParams,outputSolutions);

[solutions,~] = controller{{X0,Xf}};

u = solutions{1};
x = zeros(6,Nsim);
x(:,1) = x0;
for ii = 1:length(Nsim)
    x(:,ii+1) = A*x(:,ii) + B*u(:,ii);
end

figure
hold on
grid on
plot3(x(1,:),x(2,:),x(3,:),'k','linewidth',2)
axis tight
