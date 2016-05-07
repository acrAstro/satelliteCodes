classdef HCW < handle
    properties
        mu
        a
        initialConditions
        n
        period
        
        samples
        
        t0
        dt
        tf
        time
        
        Phi
        A
        B
        X
        Ak
        Bk
    end
    
    methods
        function obj = HCW(initStruct)
            obj.mu = initStruct.params{1};
            obj.a = initStruct.params{2};
            obj.initialConditions = initStruct.params{3};
            obj.n = sqrt(obj.mu/obj.a^3);
            obj.period = 2*pi/obj.n;
            obj.samples = initStruct.maneuverParams{1};
            obj.B = initStruct.maneuverParams{2};
            obj.A = HCW_Matrix(obj.n);
            
            obj.t0 = initStruct.timeParams.t0;
            obj.dt = initStruct.timeParams.dt;
            obj.tf = initStruct.timeParams.tf;
            
            obj.makeTimeVector();
        end
        
        function obj = makeTimeVector(obj)
            obj.time = obj.t0:obj.dt:obj.tf;
        end
        
        function obj = propagateModel(obj)
            obj.A = HCW_Matrix(obj.n);
            Ad = expm(obj.A*obj.dt);
            for ii = 1:length(obj.time)-1
                obj.Phi(:,:,ii) = Ad^ii;
            end
        end
        
        function obj = makeDiscreteMatrices(obj)
            sysc = ss(obj.A,obj.B,[],[]);
            sysd = c2d(sysc,obj.dt,'zoh');
            obj.Ak = sysd.A;
            obj.Bk = sysd.B;
        end
        
        function obj = propagateState(obj)
            obj.propagateModel();
            obj.X = zeros(6, length(obj.time));
            obj.X(:,1) = obj.initialConditions;
            for ii = 1:length(obj.time)-1
                obj.X(:,ii+1) = obj.Phi(:,:,ii)*obj.initialConditions;
            end
        end
    end
end

function A = HCW_Matrix(n)
A = [0     0 0     1     0 0;
     0     0 0     0     1 0;
     0     0 0     0     0 1;
     3*n^2 0 0     0   2*n 0;
     0     0 0    -2*n   0 0;
     0     0 -n^2  0     0 0];
end
