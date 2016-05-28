clear all; close all; clc;

mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
nu = 3;
nx = 6;
mass = 1;
[Ac,Bc] = HCW_Matrices(n,nu,mass);

period = 2*pi/n;
dt = 2;
% numPeriod = 1/10;
% tf = numPeriod*period;
tf = 300;
% Nsim = round(tf/dt);
% Time = 0:dt:tf;
sysc = ss(Ac,Bc,[],[]);
sysd = c2d(sysc,dt,'zoh');

A = sysd.a;
B = sysd.b;

umin = -0.03;
umax = 0.03;

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

X0 = [x0 y0 z0 vx0 vy0 vz0]';
Xf = [xf yf zf vxf vyf vzf]';

varhi = tf; varlo = 0;
while diagnostics ~= 0
    N = round((varlo + varhi)/2);
    u = sdpvar(nu,N);
    uslack = sdpvar(nu,N);
    
    x = sdpvar(6,N+1);
    
    x1 = sdpvar(6,1);
    x2 = sdpvar(6,1);
    
    %% Build constraints and cost function
    constraints = [];
    objective = 0;
    
    constraints = [constraints, x(:,1) == X0];
    for kk = 1:N
        for jj = 1:nu
            objective = objective + uslack(jj,kk);
        end
        constraints = [constraints, x(:,kk+1) == A*x(:,kk) + B*u(:,kk)];
        for ii = 1:nu
            constraints = [constraints, u(ii,kk) <= uslack(ii,kk)];
            constraints = [constraints, -u(ii,kk) <= uslack(ii,kk)];
            constraints = [constraints, umin <= u(ii,kk) <= umax];
        end
    end
    constraints = [constraints, x(:,N+1) == Xf];
    
    options = sdpsettings('solver','gurobi','saveyalmipmodel',1,'verbose',0);
    
    parameters_in = {x1,x2};
    solutions_out = {u,x,objective};
    
    controller = optimizer(constraints,objective,options,parameters_in,...
        solutions_out);
    
    [solutions,flag] = controller{{X0,Xf}};
    flag
    U = solutions{1};
    X = solutions{2};
    objectiveOpt = solutions{3};
    disp([N,objectiveOpt])
    if round(objectiveOpt) == 0
        varhi = N;
        Nsim = N;
        Xs = X;
        Us = U;
    else
        varlo = N;
    end
    iter = iter + 1;
end
Time = 0:dt:Nsim-1;
% [Xq,Uq] = quivThrust(Time(1:end-1),transpose(X),transpose(U),nu,1);

switch nu
    case 2
        figure
        hold on
        grid on
        plot3(X(1,:),X(2,:),X(3,:),'k','linewidth',2)
        axis tight
        plot3(x0,y0,z0,'b.','MarkerSize',25)
        plot3(xf,yf,zf,'r.','MarkerSize',25)
%         quiver3(Xq(:,1),Xq(:,2),Xq(:,3),zeros(size(Uq(:,1))),Uq(:,1),Uq(:,2),'r','LineWidth',2);
        title1 = title('Trajectory in Space');
        xl = xlabel('Radial, $x$, m');
        yl = ylabel('In-track, $y$, m');
        zl = zlabel('Cross-track, $z$, m');
        leg1 = legend('$x(t)$','$x_0$','$x_f$','Thrust Plume','location','best');
        set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
        view([-90,90])
        
        figure
        subplot(211)
        hold on
        grid on
        plot(Time(1:end-1),U(1,:),'k','linewidth',2)
        plot([Time(1), Time(end)],[umax, umax],'k--','LineWidth',2)
        plot([Time(1), Time(end)],[umin, umin],'k--','LineWidth',2)
        axis([Time(1),Time(end),2*umin,2*umax])
        subplot(212)
        hold on
        grid on
        plot(Time(1:end-1),U(2,:),'k','linewidth',2)
        plot([Time(1), Time(end)],[umax, umax],'k--','LineWidth',2)
        plot([Time(1), Time(end)],[umin, umin],'k--','LineWidth',2)
        axis([Time(1),Time(end),2*umin,2*umax])
    case 3
        
        figure
        hold on
        grid on
        plot3(X(1,:),X(2,:),X(3,:),'k','linewidth',2)
        axis tight
        plot3(x0,y0,z0,'b.','MarkerSize',25)
        plot3(xf,yf,zf,'r.','MarkerSize',25)
%         quiver3(Xq(:,1),Xq(:,2),Xq(:,3),Uq(:,1),Uq(:,2),Uq(:,3),'r','LineWidth',2);
        title1 = title('Trajectory in Space');
        xl = xlabel('Radial, $x$, m');
        yl = ylabel('In-track, $y$, m');
        zl = zlabel('Cross-track, $z$, m');
        leg1 = legend('$x(t)$','$x_0$','$x_f$','Thrust Plume','location','best');
        set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
        view([-90,90])
        
        
        figure
        subplot(311)
        hold on
        grid on
        plot(Time(1:end-1),U(1,:),'k','linewidth',2)
        plot([Time(1), Time(end)],[umax, umax],'k--','LineWidth',2)
        plot([Time(1), Time(end)],[umin, umin],'k--','LineWidth',2)
        axis([Time(1),Time(end),2*umin,2*umax])
        subplot(312)
        hold on
        grid on
        plot(Time(1:end-1),U(2,:),'k','linewidth',2)
        plot([Time(1), Time(end)],[umax, umax],'k--','LineWidth',2)
        plot([Time(1), Time(end)],[umin, umin],'k--','LineWidth',2)
        axis([Time(1),Time(end),2*umin,2*umax])
        subplot(313)
        hold on
        grid on
        plot(Time(1:end-1),U(3,:),'k','linewidth',2)
        plot([Time(1), Time(end)],[umax, umax],'k--','LineWidth',2)
        plot([Time(1), Time(end)],[umin, umin],'k--','LineWidth',2)
        axis([Time(1),Time(end),2*umin,2*umax])
end
