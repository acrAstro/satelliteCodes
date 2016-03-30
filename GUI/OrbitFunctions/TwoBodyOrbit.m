classdef TwoBodyOrbit < handle
% This is a basic two-body orbit perturbed by J2. The user is required only
% to input the initial orbit based on Kepler elements, a few physical
% parameters, and the total time of flight.
%
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%
%% Fundamental properties of class
    properties
        initialConditions   % Derived from Kepler elements, IC's for ode45
        kepElems            % Initial Kepler elements of orbit
        J2                  % J2 zonal harmonic coefficient
        dragCoefficient     % Drag coefficient of satellite, not used yet
        time                % Derived from Kepler elements and TOF
        t0                  % Initial time
        period              % Derived from Kepler elements
        numPeriod           % Time of flight
        mu                  % Gravitational parameter of Earth
        states              % R and V histories from ode45
        OsculatingElements  % Kepler element histories from R and V
        mass                % Spacecraft mass, not used yet
        Req                 % Equatorial radius of Earth
    end

%% Functionality of class
    methods
        function TBP = TwoBodyOrbit(kepElems,J2,mu,Req,t0,numPeriod)
            % Constructor object
            TBP.kepElems = kepElems;
            TBP.J2 = J2;
            TBP.mu = mu;
            TBP.Req = Req;
            TBP.t0 = t0;
            TBP.numPeriod = numPeriod;
        end
        
        function TBP = makeTimeVector(TBP)
            % Take constructor and derive time vector
            n = sqrt(TBP.mu/TBP.kepElems(1)^3); % Mean motion
            if isempty(TBP.t0)
                TBP.t0 = 0;
            else
            end
            TBP.period = 2*pi/n;                % Orbital period
            tf = TBP.numPeriod*TBP.period;      % Time of flight
            TBP.time = linspace(TBP.t0,tf,100*TBP.numPeriod);
        end
        
        function TBP = propagateOrbit(TBP)
            % Take constructor and propagate orbit
            options = odeset('RelTol',1e-9,'AbsTol',1e-12); % ode45 options
            TBP.makeTimeVector(); % Instantiates the time vector
            [~,X] = ode45(@orbitEquation,TBP.time,TBP.initialConditions,...
                options,TBP.mu,TBP.J2,TBP.Req); % Integrates equations
            TBP.states = X;
            TBP.getKepElems(); % Derives Kepler elements from R and V
        end
        
        function TBP = setInitialConditions(TBP)
            % Takes Kepler elements and returns R and V for ode45
            [R,V] = oe2rv(TBP.kepElems,TBP.mu,'r');
            TBP.initialConditions = [R(:); V(:)];
        end
        
        function TBP = getKepElems(TBP)
            % Take R and V and return Kepler elements over time
            oscElems = zeros(6,length(TBP.time));
            for ii = 1:length(TBP.time)
                oscElems(:,ii) = rv2oe(TBP.states(ii,:),TBP.mu);
            end
            TBP.OsculatingElements = oscElems;
        end
        
        function TBP = plotOrbit(TBP)
            % Plots the orbit in R^3, makes a sphere as well to represent
            % Earth
            [xs,ys,zs] = sphere(30);
            xs = TBP.Req.*xs; ys = TBP.Req.*ys; zs = TBP.Req.*zs;
            figure
            hold on
            grid on
            plot3(TBP.states(:,1),TBP.states(:,2),TBP.states(:,3),'k',...
                'linewidth',2);
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