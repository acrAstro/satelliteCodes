classdef TwoBodyOrbit < handle
    
    properties
        initialConditions
        kepElems
        J2
        dragCoefficient
        time
        t0
        period
        numPeriod
        mu
        states
        OsculatingElements
        mass
        Req
    end
    
    methods
        function TBP = TwoBodyOrbit(kepElems,J2,mu,Req,t0,numPeriod)
            TBP.kepElems = kepElems;
            TBP.J2 = J2;
            TBP.mu = mu;
            TBP.Req = Req;
            TBP.t0 = t0;
            TBP.numPeriod = numPeriod;
        end
        
        function TBP = makeTimeVector(TBP)
            n = sqrt(TBP.mu/TBP.kepElems(1)^3);
            if isempty(TBP.t0)
                TBP.t0 = 0;
            else
            end
            TBP.period = 2*pi/n;
            tf = TBP.numPeriod*TBP.period;
            TBP.time = linspace(TBP.t0,tf,100*TBP.numPeriod);
        end
        
        function TBP = propagateOrbit(TBP)
            options = odeset('RelTol',1e-9,'AbsTol',1e-12);
            TBP.makeTimeVector();
            [~,X] = ode45(@orbitEquation,TBP.time,TBP.initialConditions,options,TBP.mu,TBP.J2,TBP.Req);
            TBP.states = X;
            TBP.getKepElems();
        end
        
        function TBP = setInitialConditions(TBP)
            [R,V] = oe2rv(TBP.kepElems,TBP.mu,'r');
            TBP.initialConditions = [R(:); V(:)];
        end
        
        function TBP = getKepElems(TBP)
            oscElems = zeros(6,length(TBP.time));
            for ii = 1:length(TBP.time)
                oscElems(:,ii) = rv2oe(TBP.states(ii,:),TBP.mu);
            end
            TBP.OsculatingElements = oscElems;
        end
        
        function TBP = plotOrbit(TBP)
            [xs,ys,zs] = sphere(30);
            xs = TBP.Req.*xs; ys = TBP.Req.*ys; zs = TBP.Req.*zs;
            figure
            hold on
            grid on
            plot3(TBP.states(:,1),TBP.states(:,2),TBP.states(:,3),'k','linewidth',2);
            surf(xs,ys,zs);
            title1 = title('$J_2$-Perturbed Orbit');
            xl = xlabel('$X$, km');
            yl = ylabel('$Y$, km');
            zl = zlabel('$Z$, km');
            leg1 = legend('Orbit','Location','Best');
            axis equal
            set([title1,xl,yl,zl,leg1],'interpreter','latex','fontsize',12);
        end

    end
end