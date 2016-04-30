classdef HCW < handle
    properties
        mu
        a
        initialConditions
        n
        period
        
        t0
        dt
        tf
        time
        
        Phi
        X
    end
    
    methods
        function obj = HCW(initStruct)
            obj.mu = initStruct.params{1};
            obj.a = initStruct.params{2};
            obj.initialConditions = initStruct.params{3};
            obj.n = sqrt(obj.mu/obj.a^3);
            obj.period = 2*pi/obj.n;
            
            obj.t0 = initStruct.timeParams{1};
            obj.dt = initStruct.timeParams{2};
            obj.tf = initStruct.timeParams{3};
            
            obj.makeTimeVector();
            
            A = HCW_Matrix(obj.n);
            obj.Phi = expm(A*obj.dt);
        end
        
        function obj = makeTimeVector(obj)
            obj.time = obj.t0:obj.dt:obj.tf;
        end
        
        function obj = propagateModel(obj)
            obj.X = zeros(6, length(obj.time));
            obj.X(:,1) = obj.initialConditions;
            for ii = 1:length(obj.time)-1
                obj.X(:,ii+1) = obj.Phi^ii*obj.initialConditions;
            end
        end
        
        function obj = plotOrbit(obj)
            x = obj.X(1,:);
            y = obj.X(2,:);
            z = obj.X(3,:);
            
            figure
            hold on
            grid on
            plot3(x,y,z,'k-','LineWidth',2)
            axis tight
            title1 = title('HCW Relative Motion');
            xl = xlabel('Radial, $x$, km');
            yl = ylabel('In-track, $y$, km');
            zl = zlabel('Cross-track, $z$, km');
            set([title1,xl,yl,zl],'interpreter','latex','fontsize',12);
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
