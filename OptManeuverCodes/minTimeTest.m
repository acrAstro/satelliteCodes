clear; close all; clc;
cvx_solver sdpt3
a = 6678e3;
mu = 3.986e14;
mm = sqrt(mu/a^3);
mass = 20;

Ac = [0      0     0     1    0 0;
      0      0     0     0    1 0;
      0      0     0     0    0 1;
      3*mm^2 0     0     0 2*mm 0;
      0      0     0 -2*mm    0 0;
      0      0 -mm^2     0    0 0];

Bc = 1/mass.*[0 0 0;
              0 0 0;
              0 0 0;
              1 0 0;
              0 1 0;
              0 0 1];

Cc = [eye(3) zeros(3,3)];
Dc = [];
Ts = 1;
sysc = ss(Ac,Bc,Cc,Dc);
sysd = c2d(sysc,Ts);
[A,H,~,~] = ssdata(sysd);
B = [H -H];

umax = 0.1;

x0 = -100;
y0 = 0;
z0 = 100;
vx0 = 0;
vy0 = 0;
vz0 = 0;
X0 = [x0 y0 z0 vx0 vy0 vz0]';

xf = 100;
yf = 0;
zf = 0;
vxf = 0;
vyf = -2*mm*xf;
vzf = 0;
Xf = [xf yf zf vxf vyf vzf]';


n = 6;
nu = size(B,2);

varlo = 50;
varhi = 1000;
N = varhi;
err = 1;
% cvx_quiet(true);

while (varhi ~= varlo)
%     N = round((varlo + varhi)/2);
    cvx_begin
        variable X(n,N+1);
        variable U(nu,N);
    subject to
        X(:,2:N+1) == A*X(:,1:N) + B*U;
        X(:,1) == X0;
        X(:,N+1) == Xf;
%         X(1:3,round((N+1)/2)) == xdesi;
        U >= 0;
        U <= umax;
    cvx_end
    if ~isempty(strfind(cvx_status,'Solved'))
        varhi = N;
    else
        varlo = N + 1;
    end
    if (varhi - varlo > 1)
        N = varlo + round((varhi - varlo)/2);
    else
        N = varlo;
    end
    
end

N = varhi;
cvx_begin
        variable Xs(n,N+1);
        variable Us(nu,N);
    subject to
        Xs(:,2:N+1) == A*Xs(:,1:N) + B*Us;
        Xs(:,1) == X0;
        Xs(:,N+1) == Xf;
%         X(1:3,round((N+1)/2)) == xdesi;
        Us >= 0;
        Us <= umax;
    cvx_end

figure(1)
hold on
grid on
plot3(Xs(1,:),Xs(2,:),Xs(3,:),'k','LineWidth',2)
plot3(X0(1),X0(2),X0(3),'r.','MarkerSize',25)
plot3(Xf(1),Xf(2),Xf(3),'b.','MarkerSize',25)
% plotShuttle(0,0,0,0,-pi/2,pi/2,0.0487,1e-3,[1,1,0.5])
leg1 = legend('$X(t)$','$X_0$','$X_f$','Location','Best');
axis tight
set(leg1,'interpreter','latex','fontsize',12)
view([-90 90])

tu = 0:Ts:(N-1);
figure(2)
subplot(211)
hold on
grid on
stairs(tu,Us(1,:)-Us(3,:))
axis tight
subplot(212)
hold on
grid on
stairs(tu,Us(2,:)-Us(4,:))
axis tight