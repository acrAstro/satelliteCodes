classdef HCW < handle
    properties
        mu
        a
        initialConditions
        n
        period
        numInputs
        
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
            obj.numInputs = initStruct.params{4};
            obj.B = inputMatrix(obj.numInputs);
            
            obj.t0 = initStruct.timeParams{1};
            obj.dt = initStruct.timeParams{2};
            obj.tf = initStruct.timeParams{3};
            
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

function B = inputMatrix(numInput)

method = numInput;
switch method
    case 0
        B = zeros(6,1);
    case 2
        B = [zeros(4,2); eye(2)];
    case 3
        B = [zeros(3,3); eye(3)];
    case 4
        B = [zeros(4,4); eye(2), -eye(2)];
    case 6
        B = [zeros(3,6); eye(3), -eye(3)];
end

end
