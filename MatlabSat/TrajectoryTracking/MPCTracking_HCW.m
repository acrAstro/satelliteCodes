classdef MPCTracking_HCW < handle
    
    properties
    X
    U
    T
    
    Ref
    URef
    TRef
    
    optimalObjective
    
    N  % Horizon
    dt
    ny
    Q
    R
    M
    umax
    umin
    
    n
    mass
    Ac
    Bc
    Cc
    A
    B
    C
    
    controller
    
    baseBody
    
    chiefOrbit
    end
    
    methods
        function obj = MPCTracking_HCW(initStruct,chiefStruct)
            obj.Ref         = initStruct.refTrajectory;
            obj.URef        = initStruct.refControl;
            obj.TRef        = initStruct.refTrajectoryTime;
            
            obj.n           = initStruct.meanMotion;
            obj.mass        = initStruct.mass;
            obj.ny          = initStruct.controlParams.numOutput;
            [obj.Ac,obj.Bc,obj.Cc] = makeHCWMatrices(obj.n,obj.mass,obj.ny);
            sysc            = ss(obj.Ac,obj.Bc,obj.Cc,[]);
            sysd            = c2d(sysc,obj.dt,'zoh');
            obj.A           = sysd.a;
            obj.B           = sysd.b;
            obj.C           = sysd.c;
%             obj.baseBody    = initStruct.baseBody;
            obj.N           = initStruct.controlParams.horizon;
            obj.dt          = initStruct.controlParams.timeStep;
            obj.Q           = initStruct.controlParams.statePenalty;
            obj.R           = initStruct.controlParams.controlPenalty;
            [~,obj.M,~]     = dlqr(obj.A,obj.B,obj.Q,obj.R);
            obj.umax        = initStruct.controlParams.upperBound;
            obj.umin        = initStruct.controlParams.lowerBound;
            obj.chiefOrbit  = TwoBodyOrbit(chiefStruct);
        end
        
        function obj = buildConvexMPC_OneNorm(obj)
            u       = sdpvar(3,obj.N);
            ur      = sdpvar(3,obj.N);
            x       = sdpvar(6,obj.N+1);
            r       = sdpvar(obj.ny,obj.N+1);
            e       = sdpvar(6,obj.N+1);
            
            constraints = [];
            objective = 0;
            for kk = 1:obj.N
                e(:,kk) = x(:,kk) - r(:,kk);
                ubar = u(:,kk) - ur(:,kk);
                objective = objective + norm(obj.Q*e(:,kk),1) + norm(obj.R*ubar,1);
                constraints = [constraints, e(:,kk+1) == obj.A*e(:,kk) + obj.B*ubar];
                constraints = [constraints, 0 <= u(:,kk) <= obj.umax];
            end
            objective = objective + norm(obj.M*e(:,obj.N+1),1);
            params_in = {x(:,1),[r],[ur]};
            solutions_out = {[u],objective};
            obj.controller = optimizer(constraints,objective,options,params_in,solutions_out);
        end
        
        function obj = buildConvexMPC_Quadratic(obj)
            u       = sdpvar(3,obj.N);
            ur      = sdpvar(3,obj.N);
            x       = sdpvar(6,obj.N+1);
            r       = sdpvar(obj.ny,obj.N+1);
            e       = sdpvar(6,obj.N+1);
            
            constraints = [];
            objective = 0;
            for kk = 1:obj.N
                e(:,kk) = x(:,kk) - r(:,kk);
                ubar = u(:,kk) - ur(:,kk);
                objective = objective + transpose(e(:,kk))*obj.Q*e(:,kk) + transpose(ubar)*obj.R*ubar;
                constraints = [constraints, e(:,kk+1) == obj.A*e(:,kk) + obj.B*ubar];
                constraints = [constraints, 0 <= u(:,kk) <= obj.umax];
            end
            objective = objective + transpose(e(:,obj.N+1))*obj.M*e(:,obj.N+1);
            params_in = {x(:,1),[r],[ur]};
            solutions_out = {[u],objective};
            obj.controller = optimizer(constraints,objective,options,params_in,solutions_out);
        end
        
        function obj = buildMIMPC_OneNorm(obj)
            u       = semivar(3,obj.N);
            ur      = sdpvar(3,obj.N);
            x       = sdpvar(6,obj.N+1);
            r       = sdpvar(obj.ny,obj.N+1);
            e       = sdpvar(6,obj.N+1);
            
            constraints = [];
            objective = 0;
            for kk = 1:obj.N
                e(:,kk) = x(:,kk) - r(:,kk);
                ubar = u(:,kk) - ur(:,kk);
                objective = objective + norm(obj.Q*e(:,kk),1) + norm(obj.R*ubar,1);
                constraints = [constraints, e(:,kk+1) == obj.A*e(:,kk) + obj.B*ubar];
                constraints = [constraints, obj.umin <= u(:,kk) <= obj.umax];
            end
            objective = objective + norm(obj.M*e(:,obj.N+1),1);
            params_in = {x(:,1),[r],[ur]};
            solutions_out = {[u],objective};
            obj.controller = optimizer(constraints,objective,options,params_in,solutions_out);
        end
        
        function obj = buildMIMPC_Quadratic(obj)
            u       = semivar(3,obj.N);
            ur      = sdpvar(3,obj.N);
            x       = sdpvar(6,obj.N+1);
            r       = sdpvar(obj.ny,obj.N+1);
            e       = sdpvar(6,obj.N+1);
            
            constraints = [];
            objective = 0;
            for kk = 1:obj.N
                e(:,kk) = x(:,kk) - r(:,kk);
                ubar = u(:,kk) - ur(:,kk);
                objective = objective + transpose(e(:,kk))*obj.Q*e(:,kk) + transpose(ubar)*obj.R*ubar;
                constraints = [constraints, e(:,kk+1) == obj.A*e(:,kk) + obj.B*ubar];
                constraints = [constraints, obj.umin <= u(:,kk) <= obj.umax];
            end
            objective = objective + transpose(e(:,obj.N+1))*obj.M*e(:,obj.N+1);
            params_in = {x(:,1),[r],[ur]};
            solutions_out = {[u],objective};
            obj.controller = optimizer(constraints,objective,options,params_in,solutions_out);
        end
        
        function obj = simulateSystem(obj)
            
        end
        
    end
end

function [A,B,C] = makeHCWMatrices(n,mass,ny)
a = 3*n^2;
b = 2*n;
c = -2*n;
d = -n^2;

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     a 0 0 0 b 0;
     0 0 0 c 0 0;
     0 0 d 0 0 0];
 
B = 1/mass.*[zeros(3,3); eye(3)];

switch ny
    case 3
        C = [eye(3), zeros(3,3)];
    case 6
        C = eye(6);
end
end