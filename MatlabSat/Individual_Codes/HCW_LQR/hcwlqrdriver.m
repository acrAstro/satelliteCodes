clear all; close all; clc;

mu = 3.986e14;
a = 7000e3;
n = sqrt(mu/a^3);
period = 2*pi/n;
numPeriod = 1;

numInputs = 3;

[A,B] = hcwmatrices(n,numInputs);

Q = eye(6);
R = 1e10.*eye(numInputs);

[K,P,E] = lqr(A,B,Q,R);

t0 = 0;
tf = numPeriod*period;

time = linspace(t0,tf,300);

X0 = [100 -100 100 0 -2*n*100 0]';
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[Time,X1] = ode45(@hcw,time,X0,options,K,n,numInputs);
[~,X2] = ode45(@hcw,time,X0,options,[],n,numInputs);

for ii = 1:length(time)
    U(:,ii) = -K*(X1(ii,:)');
end

Time = 1/period.*Time;

figure
hold on
grid on
plot3(X1(:,1),X1(:,2),X1(:,3),'k-','LineWidth',2)
plot3(X2(:,1),X2(:,2),X2(:,3),'r-','LineWidth',2)
plot3(X0(1),X0(2),X0(3),'rd','MarkerSize',8);
plot3(0,0,0,'b.','MarkerSize',25)
leg1 = legend('LQR','Free','$X_0$','$X_f$','Location','Best');
xl = xlabel('Radial');
yl = ylabel('In-track');
zl = zlabel('Cross-track');
title1 = title(['LQR Rendezvous with ' num2str(numInputs) ' Inputs']);
set([leg1 xl yl zl title1],'interpreter','latex','fontsize',12)
axis tight

figure
subplot(311)
hold on
grid on
plot(Time,X1(:,1),'k','LineWidth',2)
yl1 = ylabel('Radial');
title1 = title(['LQR State History with ' num2str(numInputs) ' Inputs']);
axis tight
subplot(312)
hold on
grid on
plot(Time,X1(:,2),'k','LineWidth',2)
yl2 = ylabel('In-track');
axis tight
subplot(313)
hold on
grid on
plot(Time,X1(:,3),'k','LineWidth',2)
yl3 = ylabel('Cross-track');
xl = xlabel('Time, Orbits');
set([title1 yl1 yl2 yl3 xl],'interpreter','latex','fontsize',12)
axis tight


if numInputs == 3
    figure
    subplot(311)
    hold on
    grid on
    plot(Time,U(1,:),'k','LineWidth',2)
    yl1 = ylabel('Radial');
    title1 = title(['LQR Control with ' num2str(numInputs) ' Inputs']);
    axis tight
    subplot(312)
    hold on
    grid on
    plot(Time,U(2,:),'k','LineWidth',2)
    yl2 = ylabel('In-track');
    axis tight
    subplot(313)
    hold on
    grid on
    plot(Time,U(3,:),'k','LineWidth',2)
    yl3 = ylabel('Cross-track');
    xl = xlabel('Time, Orbits');
    set([title1 yl1 yl2 yl3 xl],'interpreter','latex','fontsize',12)
    axis tight
elseif numInputs == 2
    figure
    subplot(211)
    hold on
    grid on
    plot(Time,U(1,:),'k','LineWidth',2)
    yl1 = ylabel('In-track');
    title1 = title(['LQR Control with ' num2str(numInputs) ' Inputs']);
    axis tight
    subplot(212)
    hold on
    grid on
    plot(Time,U(2,:),'k','LineWidth',2)
    yl2 = ylabel('Cross-track');
    xl = xlabel('Time, Orbits');
    set([title1 yl1 yl2 xl],'interpreter','latex','fontsize',12)
    axis tight
else
end