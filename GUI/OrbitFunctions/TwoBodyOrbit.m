classdef TwoBodyOrbit
    
    properties
        initialConditions
        kepElems
        J2
        dragCoefficient
        time
        mu
        states
        mass
%         xArea
    end
    methods
        function TBP = TwoBodyOrbit(kepElems,J2,Time,mu)
            TBP.kepElems = kepElems;
            TBP.initialConditions = setInitialConditions(kepElems);
            TBP.J2 = J2;
            TBP.time = Time;
            TBP.mu = mu;
        end
        
        function TBP = propagateOrbit(TBP)
            options = odeset('RelTol',1e-9,'AbsTol',1e-12);
            tspan = TBP.time;
            X0 = TBP.initialConditions;
            [~,X] = ode45(@orbitEquation,tspan,X0,options,TBP.mu,TBP.J2);
            TBP.states = X;
        end
        
        function TBP = setInitialConditions(TBP)
            [R,V] = oe2rv(TBP.kepElems,TBP.mu,'r');
            TBP.initialConditions = [R(:); V(:)];
        end
    end
end