classdef hcwOpt < handle
    % This class constructs a fully-actuated or under-actuated, minimum-fuel
    % transfer for the HCW equations.
    properties
        mu                  % gravitational parameter (m^3/s^2)
        a                   % semi-major axis of chief orbit (m)
        n                   % mean motion of chief orbit (rad/s)
        period              % period of chief orbit (s)
        nu                  % number of inputs (must be 2 or 3)
        mass                % mass of deputy (kg)
        Lb                  % lower bound on inputs
        Ub                  % upper bound on inputs
        Ac                  % continuous HCW state matrix
        Bc                  % continuous input matrix
        A                   % discrete HCW state matrix
        B                   % discrete input matrix
        dt                  % time-step (s)
        t0                  % initial time (just set it to zero, HCW are LTI) given in s
        tf                  % time-of-flight (s)
        Time                % Time vector
        Nsim                % number of simulation steps
        X0                  % initial state, [x0 y0 z0 vx0 vy0 vz0]^T in m, m/s
        Xf                  % final state, [xf yf zf vxf vyf vf]^T in m, m/s
        X                   % state history, 6 x Nsim+1 matrix
        U                   % input history, nu x Nsim matrix
        Xq                  % states that use the quiver function to plot plumes
        Uq                  % inputs that use the quives function to plot plumes
        optimalObjective    % optimal objective
        varhi               % upper bound on feasibility
        varlo               % lower bound on feasibility
        Xs                  % minimum time state solution
        Us                  % minimum time control solution
        tu                  % minimum time vector
        err                 % bisection error
        feas                % status of feasibility (1=yes,0=n)
    end
    
    methods
        % This is the constructor function, as you can see, the basic
        % minimum fuel transfer is defined simply by the 6 input
        % parameters, boundary conditions, and the time-of-flight
        function hcw = hcwOpt(initStruct)
            hcw.mu      = initStruct.params{1};
            hcw.a       = initStruct.params{2};
            hcw.nu      = initStruct.params{3};
            hcw.mass    = initStruct.params{4};
            hcw.Ub      = initStruct.params{5};
            hcw.Lb      = initStruct.params{6};
            hcw.n       = sqrt(hcw.mu/hcw.a^3);
            hcw.t0      = initStruct.timeParams{1};
            hcw.dt      = initStruct.timeParams{2};
            hcw.tf      = initStruct.timeParams{3};
            [hcw.Ac,hcw.Bc] = HCW_Matrices(hcw.n,hcw.nu,hcw.mass);
            sysc = ss(hcw.Ac,hcw.Bc,[],[]);
            sysd = c2d(sysc,hcw.dt,'zoh');
            hcw.A = sysd.a;
            hcw.B = sysd.b;
            hcw.makeTimeVector();
            hcw.X0 = initStruct.X0;
            hcw.Xf = initStruct.Xf;
            try
                exist(initStruct.feasParams,'var')
            catch ME
                if strcmp(ME.identifier,'MATLAB:nonExistentField') == 1
                    hcw.varhi = [];
                    hcw.varlo = [];
                    hcw.err = [];
                else
                    hcw.varhi = initStruct.feasParams{1};
                    hcw.varlo = initStruct.feasParams{2};
                    hcw.err = initStruct.feasParams{3};
                end
            end
        end
        
        % This function makes the time vector (probably not necessary, but
        % all of my other codes have something similar, and they need it
        function hcw = makeTimeVector(hcw)
            hcw.Time = hcw.t0:hcw.dt:hcw.tf;
            hcw.Nsim = round(hcw.tf/hcw.dt);
            hcw.period = 2*pi/hcw.n;
        end
        
        % This function defines the optimal control problem and then solves
        % it and stores the output in the object fields we didn't use
        % before. NOTE: the user MUST have the free Matlab toolbox CVX
        % found at:
        %
        % http://cvxr.com/cvx/
        %
        % I also use the Gurobi optimizer for the linear program, found at:
        %
        % http://www.gurobi.com/
        %
        % and the SDPT3 toolbox found at:
        %
        % http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
        %
        function hcw = feasibility(hcw)
            [hcw.tu,hcw.Xs,hcw.Us] = minTimeFeasibility(hcw.A,hcw.B,...
                hcw.varhi,hcw.varlo,hcw.err,hcw.X0,hcw.Xf);
        end
        
        % This function defines the optimal control problem and then solves
        % it and stores the output in the object fields we didn't use
        % before. NOTE: the user MUST have the free Matlab toolbox YALMIP
        % found at:
        %
        % http://users.isy.liu.se/johanl/yalmip/
        %
        % I also use the Gurobi optimizer for the linear program, found at:
        %
        % http://www.gurobi.com/
        %
        function hcw = fuelOptimalTransfer(hcw)
            % Define the input sdpvars and slack sdpvars
            u = sdpvar(hcw.nu,hcw.Nsim);
            uslack = sdpvar(hcw.nu,hcw.Nsim);
            % The state variables
            x = sdpvar(6,hcw.Nsim+1);
            % The initial conditions
            x1 = sdpvar(6,1);
            x2 = sdpvar(6,1);
            % Build constraints and cost function
            constraints = [];
            objective = 0;
            % The constraints are concatenated pointwise in time
            constraints = [constraints, x(:,1) == hcw.X0];
            for kk = 1:hcw.Nsim
                % The cost funtion is additive pointwise in time
                for jj = 1:hcw.nu
                    objective = objective + uslack(jj,kk);
                end
                constraints = [constraints, x(:,kk+1) == hcw.A*x(:,kk) + hcw.B*u(:,kk)];
                for ii = 1:hcw.nu
                    constraints = [constraints, u(ii,kk) <= uslack(ii,kk)];
                    constraints = [constraints, -u(ii,kk) <= uslack(ii,kk)];
                    constraints = [constraints, hcw.Lb <= u(ii,kk) <= hcw.Ub];
                end
            end
            constraints = [constraints, x(:,hcw.Nsim+1) == hcw.Xf];
            % Settings for the optimization
            options = sdpsettings('solver','gurobi','saveyalmipmodel',1,'verbose',3);
            % The optimization problem is parameterized by the boundary
            % conditions, so we include those here
            parameters_in = {x1,x2};
            % We want to return the optimal inputs, the state history, and
            % the optimal value of the objective function, so we state that
            % here
            solutions_out = {u,x,objective};
            % The controller object constructs the problem and transforms
            % it for use with whatever optimizer you choose (Gurobi!!!)
            controller = optimizer(constraints,objective,options,parameters_in,...
                solutions_out);
            % Solve the transfer problem using the boundry conditions
            [solutions,~] = controller{{hcw.X0,hcw.Xf}};
            % Extract the relevant information from the output
            hcw.U = solutions{1};
            hcw.X = solutions{2};
            hcw.optimalObjective = solutions{3};
            % Decimate the states and inputs so we don't have tons of
            % arrows pointing everywhere
            [hcw.Xq,hcw.Uq] = quivThrust(hcw.Time(1:end-1),transpose(hcw.X),transpose(hcw.U),hcw.nu,10);
        end
        
        % This function defines the optimal control problem and then solves
        % it and stores the output in the object fields we didn't use
        % before. NOTE: the user MUST have the free Matlab toolbox YALMIP
        % found at:
        %
        % http://users.isy.liu.se/johanl/yalmip/
        %
        % I also use the Gurobi optimizer for the linear program, found at:
        %
        % http://www.gurobi.com/
        %
        function hcw = energyOptimalTransfer(hcw)
            % Define the input sdpvars
            u = sdpvar(hcw.nu,hcw.Nsim);
            % The state variables
            x = sdpvar(6,hcw.Nsim+1);
            % The initial conditions
            x1 = sdpvar(6,1);
            x2 = sdpvar(6,1);
            % Build constraints and cost function
            constraints = [];
            objective = 0;
            % The constraints are concatenated pointwise in time
            constraints = [constraints, x(:,1) == hcw.X0];
            for kk = 1:hcw.Nsim
                % The cost funtion is additive pointwise in time
                for jj = 1:hcw.nu
                    objective = objective + u(jj,kk)^2;
                end
                constraints = [constraints, x(:,kk+1) == hcw.A*x(:,kk) + hcw.B*u(:,kk)];
                for ii = 1:hcw.nu
                    constraints = [constraints, hcw.Lb <= u(ii,kk) <= hcw.Ub];
                end
            end
            constraints = [constraints, x(:,hcw.Nsim+1) == hcw.Xf];
            % Settings for the optimization
            options = sdpsettings('solver','gurobi','saveyalmipmodel',1,'verbose',3);
            % The optimization problem is parameterized by the boundary
            % conditions, so we include those here
            parameters_in = {x1,x2};
            % We want to return the optimal inputs, the state history, and
            % the optimal value of the objective function, so we state that
            % here
            solutions_out = {u,x,objective};
            % The controller object constructs the problem and transforms
            % it for use with whatever optimizer you choose (Gurobi!!!)
            controller = optimizer(constraints,objective,options,parameters_in,...
                solutions_out);
            % Solve the transfer problem using the boundry conditions
            [solutions,~] = controller{{hcw.X0,hcw.Xf}};
            % Extract the relevant information from the output
            hcw.U = solutions{1};
            hcw.X = solutions{2};
            hcw.optimalObjective = solutions{3};
            % Decimate the states and inputs so we don't have tons of
            % arrows pointing everywhere
            [hcw.Xq,hcw.Uq] = quivThrust(hcw.Time(1:end-1),transpose(hcw.X),transpose(hcw.U),hcw.nu,10);
        end
        
        function hcw = plotTransfer(hcw,color)
            if nargin < 2 || isempty(color)
                plotColor = 'k';
            else
                plotColor = color;
            end 
            % This function plots the relative trajectory and the control
            % forces
            switch hcw.nu
                case 2
                    figure(1)
                    hold on
                    grid on
                    plot3(hcw.X(1,:),hcw.X(2,:),hcw.X(3,:),plotColor,'linewidth',2)
                    axis tight
                    plot3(hcw.X0(1),hcw.X0(2),hcw.X0(3),'b.','MarkerSize',25)
                    plot3(hcw.Xf(1),hcw.Xf(2),hcw.Xf(3),'r.','MarkerSize',25)
                    quiver3(hcw.Xq(:,1),hcw.Xq(:,2),hcw.Xq(:,3),zeros(size(hcw.Uq(:,1))),hcw.Uq(:,1),hcw.Uq(:,2),'r','LineWidth',2);
                    title1 = title('Trajectory in Space');
                    xl = xlabel('Radial, $x$, m');
                    yl = ylabel('In-track, $y$, m');
                    zl = zlabel('Cross-track, $z$, m');
                    leg1 = legend('$x(t)$','$x_0$','$x_f$','Thrust Plume','location','best');
                    set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
                    view([-90,90])
                    
                    figure(2)
                    subplot(211)
                    hold on
                    grid on
                    plot(hcw.Time(1:end-1),hcw.U(1,:),plotColor,'linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    title1 = title('Under-actuated Minimum Fuel Transfer');
                    yl1 = ylabel('In-track, $u_y$, N');
                    subplot(212)
                    hold on
                    grid on
                    plot(hcw.Time(1:end-1),hcw.U(2,:),plotColor,'linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl2 = ylabel('Cross-track, $u_z$, N');
                    xl = xlabel('Time, $t$, s');
                    set([title1, yl1, yl2, xl],'interpreter','latex','fontsize',11)
                    
                case 3
                    
                    figure(1)
                    hold on
                    grid on
                    plot3(hcw.X(1,:),hcw.X(2,:),hcw.X(3,:),plotColor,'linewidth',2)
                    axis tight
                    plot3(hcw.X0(1),hcw.X0(2),hcw.X0(3),'b.','MarkerSize',25)
                    plot3(hcw.Xf(1),hcw.Xf(2),hcw.Xf(3),'r.','MarkerSize',25)
                    quiver3(hcw.Xq(:,1),hcw.Xq(:,2),hcw.Xq(:,3),hcw.Uq(:,1),hcw.Uq(:,2),hcw.Uq(:,3),'r','LineWidth',2);
                    title1 = title('Trajectory in Space');
                    xl = xlabel('Radial, $x$, m');
                    yl = ylabel('In-track, $y$, m');
                    zl = zlabel('Cross-track, $z$, m');
                    leg1 = legend('$x(t)$','$x_0$','$x_f$','Thrust Plume','location','best');
                    set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
                    view([-90,90])
                    
                    figure(2)
                    subplot(311)
                    hold on
                    grid on
                    plot(hcw.Time(1:end-1),hcw.U(1,:),plotColor,'linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    title1 = title('Fully-actuated Minimum Fuel Transfer');
                    yl1 = ylabel('Radial, $u_x$, N');
                    subplot(312)
                    hold on
                    grid on
                    plot(hcw.Time(1:end-1),hcw.U(2,:),plotColor,'linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl2 = ylabel('In-track, $u_y$, N');
                    subplot(313)
                    hold on
                    grid on
                    plot(hcw.Time(1:end-1),hcw.U(3,:),plotColor,'linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl3 = ylabel('Cross-track, $u_z$, N');
                    xl = xlabel('Time, $t$, s');
                    set([title1,yl1,yl2,yl3,xl],'interpreter','latex','fontsize',11)
            end
            
        end
        
        function hcw = plotTransferFeasibility(hcw)
            % This function plots the relative trajectory and the control
            % forces from the feasibility problem
            switch hcw.nu
                case 2
                    figure
                    hold on
                    grid on
                    plot3(hcw.Xs(1,:),hcw.Xs(2,:),hcw.Xs(3,:),'k','linewidth',2)
                    axis tight
                    plot3(hcw.X0(1),hcw.X0(2),hcw.X0(3),'b.','MarkerSize',25)
                    plot3(hcw.Xf(1),hcw.Xf(2),hcw.Xf(3),'r.','MarkerSize',25)
                    %                     quiver3(hcw.Xq(:,1),hcw.Xq(:,2),hcw.Xq(:,3),zeros(size(hcw.Uq(:,1))),hcw.Uq(:,1),hcw.Uq(:,2),'r','LineWidth',2);
                    title1 = title('Trajectory in Space');
                    xl = xlabel('Radial, $x$, m');
                    yl = ylabel('In-track, $y$, m');
                    zl = zlabel('Cross-track, $z$, m');
                    leg1 = legend('$x(t)$','$x_0$','$x_f$','location','best');
                    set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
                    view([-90,90])
                    
                    figure
                    subplot(211)
                    hold on
                    grid on
                    plot(hcw.tu,hcw.Us(1,:),'k','linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    title1 = title('Under-actuated Minimum Fuel Transfer');
                    yl1 = ylabel('In-track, $u_y$, N');
                    subplot(212)
                    hold on
                    grid on
                    plot(hcw.tu,hcw.Us(2,:),'k','linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl2 = ylabel('Cross-track, $u_z$, N');
                    xl = xlabel('Time, $t$, s');
                    set([title1, yl1, yl2, xl],'interpreter','latex','fontsize',11)
                    
                case 3
                    
                    figure
                    hold on
                    grid on
                    plot3(hcw.Xs(1,:),hcw.Xs(2,:),hcw.Xs(3,:),'k','linewidth',2)
                    axis tight
                    plot3(hcw.X0(1),hcw.X0(2),hcw.X0(3),'b.','MarkerSize',25)
                    plot3(hcw.Xf(1),hcw.Xf(2),hcw.Xf(3),'r.','MarkerSize',25)
                    %                     quiver3(hcw.Xq(:,1),hcw.Xq(:,2),hcw.Xq(:,3),hcw.Uq(:,1),hcw.Uq(:,2),hcw.Uq(:,3),'r','LineWidth',2);
                    title1 = title('Trajectory in Space');
                    xl = xlabel('Radial, $x$, m');
                    yl = ylabel('In-track, $y$, m');
                    zl = zlabel('Cross-track, $z$, m');
                    leg1 = legend('$x(t)$','$x_0$','$x_f$','location','best');
                    set([title1 xl yl zl leg1],'interpreter','latex','fontsize',11)
                    view([-90,90])
                    
                    figure
                    subplot(311)
                    hold on
                    grid on
                    plot(hcw.tu,hcw.Us(1,:),'k','linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    title1 = title('Fully-actuated Minimum Fuel Transfer');
                    yl1 = ylabel('Radial, $u_x$, N');
                    subplot(312)
                    hold on
                    grid on
                    plot(hcw.tu,hcw.Us(2,:),'k','linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl2 = ylabel('In-track, $u_y$, N');
                    subplot(313)
                    hold on
                    grid on
                    plot(hcw.tu,hcw.Us(3,:),'k','linewidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Ub, hcw.Ub],'k--','LineWidth',2)
                    plot([hcw.Time(1), hcw.Time(end)],[hcw.Lb, hcw.Lb],'k--','LineWidth',2)
                    axis([hcw.Time(1),hcw.Time(end),2*hcw.Lb,2*hcw.Ub])
                    yl3 = ylabel('Cross-track, $u_z$, N');
                    xl = xlabel('Time, $t$, s');
                    set([title1,yl1,yl2,yl3,xl],'interpreter','latex','fontsize',11)
            end
            
        end
    end
    
end

function [Xout,Uout] = quivThrust(t,x,u,nu,step)

X = x(:,1:3); U = u(:,1:nu);
for ii = 1:length(t)
    if mod(ii,step) == 0
        Xout(ii,:) = X(ii,:);
        Uout(ii,:) = -U(ii,:);
    else
    end
end
end

function [A,B] = HCW_Matrices(n,nu,mass)

a = 3*n^2;
b = 2*n;
c = -2*n;
d = n^2;

A11 = zeros(3,3); A12 = eye(3);
A21 = [a 0 0; 0 0 0; 0 0 d];
A22 = [0 b 0; c 0 0; 0 0 0];
A = [A11,A12; A21,A22];

switch nu
    case 2
        B = 1/mass.*[zeros(4,2); eye(2)];
    case 3
        B = 1/mass.*[zeros(3,3); eye(3)];
end

end

function [T,X,U] = minTimeFeasibility(A,B,varhi,varlo,err,X0,Xf)
Bt = [B,-B];
nu = size(B,2);
N = varhi;
while (varhi ~= varlo)
    Xt = zeros(6,N+1);
    Ut = zeros(nu,N);
    cvx_begin
        variable Xt(6,N+1);
        variable Ut(nu,N);
    subject to
        Xt(:,2:N+1) == A*Xt(:,1:N) + Bt*Ut;
        Xt(:,1) == X0;
        Xt(:,N+1) == Xf;
        Ut >= 0;
        Ut <= Ub;
    cvx_end
    if ~isempty(strfind(cvx_statU,'Solved'))
        varhi = N;
        feas = 1;
    else
        varlo = N + 1;
        feas = 0;
    end
    if (varhi - varlo > 1)
        N = varlo + round((varhi - varlo)/2);
    else
        N = varlo;
    end
end

if feas == 1
    N = varhi;
    cvx_begin
        variable X(n,N+1);
        variable U(nu,N);
    subject to
        X(:,2:N+1) == A*X(:,1:N) + Bt*U;
        X(:,1) == X0;
        X(:,N+1) == Xf;
        U >= 0;
        U <= Ub;
    cvx_end
else
    X = [];
    U = [];
end
end