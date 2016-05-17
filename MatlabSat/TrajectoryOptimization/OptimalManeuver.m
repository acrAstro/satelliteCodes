classdef OptimalManeuver < handle
    
    properties
        X
        U
        T
        Xq
        Uq
        optimalObjective
        
        samples
        t0
        dt
        tf
        
        umax
        umin
        B
        numInput
        numState
        
        X0
        Xf
        Nsim
        
        motionModel
        descriptor
    end
    
    methods
        function obj = OptimalManeuver(initStruct)
            obj.descriptor = initStruct.descriptor;
            obj.motionModel = eval(strcat(obj.descriptor,'(initStruct)'));
            obj.t0 = initStruct.timeParams.t0;
            obj.dt = initStruct.timeParams.dt;
            obj.tf = initStruct.timeParams.tf;
            
            if isempty(initStruct.maneuverParams{1})
                obj.samples = [];
            else
                obj.samples = initStruct.maneuverParams{1};
            end
            obj.B = initStruct.maneuverParams{2};
            obj.numInput = size(obj.B,2);
            obj.numState = size(obj.B,1);
            obj.umax = initStruct.maneuverParams{3};
            obj.umin = initStruct.maneuverParams{4};
            obj.motionModel.makeTimeVector();
            obj.Nsim = length(obj.motionModel.time);
            obj.motionModel.makeDiscreteMatrices();
            obj.T = obj.motionModel.time;
            obj.X0 = obj.motionModel.initialCondition;
            obj.Xf = obj.motionModel.terminalCondition;
        end
        
        function obj = minimizeFuel(obj)
            Bkk = [obj.motionModel.Bk,-obj.motionModel.Bk];
            u = sdpvar(2*obj.numInput,obj.Nsim);
            uslack = sdpvar(2*obj.numInput,obj.Nsim);
            x = sdpvar(obj.numState,obj.Nsim +1);
            x1 = sdpvar(obj.numState,1);
            x2 = sdpvar(obj.numState,1);
            constraints = [];
            objective = 0;
            constraints = [constraints, x(:,1) == obj.X0];
            for kk = 1:obj.Nsim
                for jj = 1:2*obj.numInput
                    objective = objective + uslack(jj,kk);
                end
                if strcmp(obj.descriptor,'HCW') == 1
                    constraints = [constraints, x(:,kk+1) == ...
                        obj.motionModel.Ak*x(:,kk) + Bkk*u(:,kk)];
                else
                    constraints = [constraints, ...
                        x(:,kk+1) == obj.motionModel.Ak(:,:,kk)*x(:,kk) + ...
                        Bkk(:,:,kk)*u(:,kk)];
                end
                for ii = 1:2*obj.numInput
                    constraints = [constraints, u(ii,kk) <= uslack(ii,kk)];
                    constraints = [constraints, -u(ii,kk) <= uslack(ii,kk)];
                    constraints = [constraints, 0 <= u(ii,kk) <= obj.umax];
                end
            end
            constraints = [constraints, x(:,obj.Nsim+1) == obj.Xf];
            options = sdpsettings('solver','gurobi','verbose',2);
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
            [solutions,~] = controller{{obj.X0,obj.Xf}};
            Uint = solutions{1};
            obj.X = solutions{2};
            obj.optimalObjective = solutions{3};
            if obj.numInput == 2
                obj.U = [Uint(1,:)-Uint(3,:); Uint(2,:)-Uint(4,:)];
            elseif obj.numInput == 3
                obj.U = [Uint(1,:)-Uint(4,:); Uint(2,:)-Uint(5,:); Uint(3,:)-Uint(6,:)];
            else
            end
            [obj.Xq,obj.Uq] = quivThrust(obj.T(1:end-1),transpose(obj.X),transpose(obj.U),obj.numInput,10);
        end
        
        function obj = minimizeEnergy_Constrained(obj)
            u = sdpvar(obj.numInput,obj.Nsim);
            x = sdpvar(obj.numState,obj.Nsim +1);
            x1 = sdpvar(obj.numState,1);
            x2 = sdpvar(obj.numState,1);
            constraints = [];
            objective = 0;
            constraints = [constraints, x(:,1) == obj.X0];
            for kk = 1:obj.Nsim
                for jj = 1:obj.numInput
                    objective = objective + u(jj,kk)^2;
                end
                if strcmp(obj.descriptor,'HCW') == 1
                    constraints = [constraints, x(:,kk+1) == ...
                        obj.motionModel.Ak*x(:,kk) + obj.motionModel.Bk*u(:,kk)];
                else
                    constraints = [constraints, ...
                        x(:,kk+1) == obj.motionModel.Ak(:,:,kk)*x(:,kk) + ...
                        obj.motionModel.Bk(:,:,kk)*u(:,kk)];
                end
                for ii = 1:obj.numInput
                    constraints = [constraints, obj.umin <= u(ii,kk) <= obj.umax];
                end
            end
            constraints = [constraints, x(:,obj.Nsim+1) == obj.Xf];
            options = sdpsettings('solver','gurobi','verbose',2);
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
            [solutions,~] = controller{{obj.X0,obj.Xf}};
            obj.U = solutions{1};
            obj.X = solutions{2};
            obj.optimalObjective = solutions{3};
            
            [obj.Xq,obj.Uq] = quivThrust(obj.T(1:end-1),transpose(obj.X),transpose(obj.U),obj.numInput,10);
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
