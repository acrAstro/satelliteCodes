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
        initialKeplerElems  % Initial Kepler elements of orbit
        J2                  % J2 zonal harmonic coefficient
        dragCoefficient     % Drag coefficient of satellite, not used yet
        safetyAltitude      % Altitude below which orbit is 'degraded'
        time                % Derived from Kepler elements and TOF
        t0                  % Initial time
        period              % Derived from Kepler elements
        numPeriod           % Time of flight
        mu                  % Gravitational parameter of Earth
        RVStates            % R and V histories from ode45
        OsculatingElements  % Kepler element histories from R and V
        mass                % Spacecraft mass, not used yet
        Req                 % Equatorial radius of Earth
        parameterization    % Parameterization, R and V or Kepler, string
    end
    
    %% Functionality of class
    methods
        function TBP = TwoBodyOrbit(initStruct)
            % Constructor object
            TBP.initialKeplerElems = initStruct.kepElems;
            TBP.J2 = initStruct.params{1};
            TBP.mu = initStruct.params{2};
            TBP.Req = initStruct.params{3};
            TBP.t0 = initStruct.params{4};
            TBP.numPeriod = initStruct.params{5};
            TBP.safetyAltitude = initStruct.params{6};
            TBP.parameterization = initStruct.Parameterization;
        end
        
        function TBP = checkParameterization(TBP)
            method = TBP.parameterization;
            switch method
                case 'RV'
                    TBP.setInitialConditions()
                case 'OE'
                    TBP.setInitialConditions()
            end
        end
        
        function TBP = makeTimeVector(TBP)
            % Take constructor and derive time vector
            n = sqrt(TBP.mu/TBP.initialKeplerElems(1)^3); % Mean motion
            if isempty(TBP.t0)
                TBP.t0 = 0;
            else
            end
            TBP.period = 2*pi/n;                % Orbital period
            tf = TBP.numPeriod*TBP.period;      % Time of flight
            TBP.time = linspace(TBP.t0,tf,100*TBP.numPeriod);
        end
        
        function TBP = propagateOrbit(TBP)
            TBP.makeTimeVector(); % Instantiates the time vector
            % Take constructor and propagate orbit
            method = TBP.parameterization;
            switch method
                case 'RV'
                    options = odeset('RelTol',1e-12,'AbsTol',1e-15); % ode45 options
                    [~,X] = ode45(@orbitEquation,TBP.time,TBP.initialConditions,...
                        options,TBP.mu,TBP.J2,TBP.Req); % Integrates equations
                    TBP.RVStates = X;
                    TBP.getKepElems(); % Derives Kepler elements from R and V
                case 'OE'
                    options = odeset('RelTol',1e-12,'AbsTol',1e-15); % ode45 options
                    [~,X] = ode45(@GaussVariationalEquations,TBP.time,TBP.initialConditions,options,TBP.mu,TBP.Req,TBP.J2);
                    TBP.OsculatingElements = X;
                    TBP.getRV();
            end
        end
        
        function TBP = setInitialConditions(TBP)
            method = TBP.parameterization;
            switch method
                case 'RV'
                    % Takes Kepler elements and returns R and V for ode45
                    [R,V] = oe2rv(TBP.initialKeplerElems,TBP.mu,'r');
                    TBP.initialConditions = [R(:); V(:)];
                    TBP.isInsideEarth();
                case 'OE'
                    TBP.initialConditions = TBP.initialKeplerElems;
                    TBP.isInsideEarth();
            end
        end
        
        function TBP = getKepElems(TBP)
            % Take R and V and return Kepler elements over time
            oscElems = zeros(6,length(TBP.time));
            for ii = 1:length(TBP.time)
                oscElems(:,ii) = rv2oe(TBP.RVStates(ii,:),TBP.mu);
            end
            TBP.OsculatingElements = oscElems;
        end
        
        function TBP = getRV(TBP)
            states = zeros(6,length(TBP.time));
            for ii = 1:length(TBP.time)
                states(:,ii) = oe2rv(TBP.OsculatingElements(ii,:),TBP.mu,'r');
            end
            TBP.RVStates = states';
        end
        
        function TBP = plotOrbit(TBP)
            % Plots the orbit in R^3, makes a sphere as well to represent
            % Earth
            [xs,ys,zs] = sphere(30);
            xs = TBP.Req.*xs; ys = TBP.Req.*ys; zs = TBP.Req.*zs;
            figure
            hold on
            grid on
            plot3(TBP.RVStates(:,1),TBP.RVStates(:,2),TBP.RVStates(:,3),'k',...
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
        
        function TBP = plotOrbitalElements(TBP)
            SMA      = TBP.OsculatingElements(:,1);
            Ecc      = TBP.OsculatingElements(:,2);
            Inc      = TBP.OsculatingElements(:,3);
            Raan     = TBP.OsculatingElements(:,4);
            argPer   = TBP.OsculatingElements(:,5);
            F        = TBP.OsculatingElements(:,6);
            Time     = 1/TBP.period.*TBP.time;
            
            figure
            hold on
            grid on
            plot(Time, SMA,'k','LineWidth',2)
            axis tight
            title1 = title('Semi-major Axis vs Time');
            xl = xlabel('Time, $n$-orbits');
            yl = ylabel('Semi-major Axis, $a$, km');
            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
            figure
            hold on
            grid on
            plot(Time, Ecc,'k','LineWidth',2)
            axis tight
            title1 = title('Eccentricity vs Time');
            xl = xlabel('Time, $n$-orbits');
            yl = ylabel('Eccentricity, $e$');
            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
            figure
            hold on
            grid on
            plot(Time, 180/pi.*Inc,'k','LineWidth',2)
            axis tight
            title1 = title('Inclination vs Time');
            xl = xlabel('Time, $n$-orbits');
            yl = ylabel('Inclination, $i$, $\deg$');
            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
            figure
            hold on
            grid on
            plot(Time, 180/pi.*Raan,'k','LineWidth',2)
            axis tight
            title1 = title('Right-Ascension vs Time');
            xl = xlabel('Time, $n$-orbits');
            yl = ylabel('Right-Ascension, $\Omega$, $\deg$');
            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
            figure
            hold on
            grid on
            plot(Time, 180/pi.*argPer,'k','LineWidth',2)
            axis tight
            title1 = title('Argument of Perigee vs Time');
            xl = xlabel('Time, $n$-orbits');
            yl = ylabel('Argument of Perigee, $\omega$, $\deg$');
            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
            %            figure
            %            hold on
            %            grid on
            %            plot(Time, 180/pi.*F,'k','LineWidth',2)
            %            axis tight
            %            title1 = title('True Anomaly vs Time');
            %            xl = xlabel('Time, $n$-orbits');
            %            yl = ylabel('True Anomaly, $f$, $\deg$');
            %            set([title1,xl,yl],'interpreter','latex','fontsize',10);
            
        end
        
        function isInsideEarth(TBP)
            if TBP.Req/1e3 < 10
                % Radius at periapsis
                rp = TBP.initialKeplerElems(1)*(1 - TBP.initialKeplerElems(2));
                periapsisAltitude = rp - TBP.Req;
                %                 fprintf(['Periapsis altitude = ' num2str(periapsisAltitude) '\n']);
                
                safeAltRatio = (TBP.safetyAltitude + TBP.Req)/TBP.Req;
                rpRatio = (rp)/TBP.Req;
                if rpRatio <= safeAltRatio
                    error(['Orbit periapsis altitude is ' num2str(periapsisAltitude) ' km, check eccentricity and semi-major axis']);
                else
                end
            elseif TBP.Req/1e3 > 10
                rp = TBP.initialKeplerElems(1)*(1 - TBP.initialKeplerElems(2));
                periapsisAltitude = rp - TBP.Req;
                fprintf(['Periapsis altitude = ' num2str(periapsisAltitude) '\n']);
                
                safeAltRatio = (TBP.safetyAltitude + TBP.Req)/TBP.Req;
                rpRatio = (rp)/TBP.Req;
                if rpRatio <= safeAltRatio
                    error(['Orbit periapsis altitude is ' num2str(periapsisAltitude) ' m, check eccentricity and semi-major axis']);
                else
                end
            else
            end
        end
    end
    
end


%% Daughter functions
function DX = orbitEquation(t,x,mu,J2,Req)
% Main function to integrate the orbital equation of motion with the
% influence of J2
% Author: Andrew Rogers, Ph.D.
% Date:   29 March 2016
%
% State vector
r = x(1:3,1);
v = x(4:6,1);
R = norm(r); % Magnitude of R

% Return J2 acceleration
aj2 = J2Vector(mu,r,J2,Req);

% Differential equation in state space form
rdot = v;
vdot = -mu/R^3.*r + aj2;
DX = [rdot(:); vdot(:)];
end

function DX = GaussVariationalEquations(t,X,mu,Req,J2)

% This function contains the necessary procedures for computing Gauss'
% Variation of Parameters. Equations derived from Battin (1987), Schaub and
% Junkins (2009) and Vallado (2007). This function also contains the
% differential equation used to describe the relative motion in cartesian
% coordinates of a chief and deputy. GVE give the chief's state, and the
% formation flying equations give the deputy's. J2 is added to provide the
% disturbance vector field but can be turned off using the J2flag option.
% The flag option calls angle arguments in either radians ('r') or degrees
% ('d'). The vector 'u' will take the place of any control vector field.

% Global options to be used throughout analysis
% global mu;

% Extract the Orbital Elements from the input vector to make typing/
% reading the equations easier.
a = X(1); ecc = X(2); inc = X(3);
W = X(4);  w = X(5); f = X(6);
% More parameters calculated for convenience
p = a*(1-ecc^2);
r = p/(1+ecc*cos(f));
h = sqrt(mu*p);
th = w+f;

% The acceleration vector is written in terms of the position of the
% satellite, so convert the orbital elements into the R and V vectors using
% oe2rv subroutine.
els = [a ecc inc W w f];
State = oe2rv(els,mu,'r');
R = State(1:3); V = State(4:6);
Rot = ijk_to_LVLH([R V]);
% Return J2 acceleration
aj2 = Rot*J2Vector(mu,R,J2,Req);
ar = aj2(1); at = aj2(2); ah = aj2(3);
% Gauss' Variational Equations
de(1) = (2*a^2/h)*(ecc*sin(f)*ar+(p/r)*at);
de(2) = (1/h)*(p*sin(f)*ar+((p+r)*cos(f)+r*ecc)*at);
de(3) = (r*cos(th)/h)*ah;
de(4) = r*sin(th)/(h*sin(inc))*ah;
de(5) = (1/(h*ecc))*(-p*cos(f)*ar+(p+r)*sin(f)*at)-(r*sin(th)...
    *cos(inc))/(h*sin(inc))*ah;
de(6) = h/r^2+1/(ecc*h)*(p*cos(f)*ar-(p+r)*sin(f)*at);

DX = de(:);
end

function aj2 = J2Vector(mu,r,J2,req)
% J2 vector
x = r(1); y = r(2); z = r(3);
R = norm(r);
gamma = -3/2*J2*(mu/R^2)*(req/R)^2;
j2x = (1 - 5*(z/R)^2)*x/R;
j2y = (1 - 5*(z/R)^2)*y/R;
j2z = (3 - 5*(z/R)^2)*z/R;

aj2 = gamma.*[j2x; j2y; j2z];
end

function [KepElems] = rv2oe(ic,mu)
% Function imports a 1x6 state vector of ephemeris data and computes the
% corresponding initial Keplerian state vector for use with Gauss'
% Equations, follows algorithm 9 (pp 120-121) in Fundamentals of ...
% Astrodynamics with Applications, David Vallado, 2007.

% Inputs:   ic         = Cartesian state vector [x,y,z,xd,yd,zd]
% Outputs:  KepElems   = 1x6 vector of Keplerian Elements
%
% Inertial unit vectors, IJK
KK = [0 0 1];
% Radial and velocity vectors
rvec = [ic(1) ic(2) ic(3)]; vvec = [ic(4) ic(5) ic(6)];
r = norm(rvec); v = norm(vvec);
% Angular momentum
hvec = cross(rvec,vvec); h = norm(hvec);
nvec = cross(KK,hvec); n = norm(nvec);
% Laplace-Runge-Lenz vector (eccentricity vector)
eccvec = (1/mu)*((v^2-mu/r).*rvec - (dot(rvec,vvec)).*vvec);
ecc = norm(eccvec);
% Specific mechanical energy
Energy = v^2/2 - mu/r;
% Compute semi-major axis, semi-minor axis and semi-latus rectum
if ecc ~= 1.0
    a = -mu/(2*Energy);
else
    a = inf;
end
b = a*sqrt(1-ecc^2);
% Inclination
inc = acosd((hvec(3))/h);
% Right-ascension of the ascending node
if nvec(2) > 0
    W = acosd((nvec(1))/n);
else
    W = 360 - acosd((nvec(1))/n);
end
% Argument of periapsis
if eccvec(3) > 0
    w = acosd((dot(nvec,eccvec))/(n*ecc));
else
    w = 360 - acosd((dot(nvec,eccvec))/(n*ecc));
end
% True anomaly
if dot(rvec,vvec) > 0
    f = acosd(dot(eccvec,rvec)/(ecc*r));
else
    f = 360 - acosd(dot(eccvec,rvec)/(ecc*r));
end

% Special cases
% Elliptic Equatorial
if inc == 0 && ecc ~= 0
    if eccvec(2) > 0
        w_true = acosd(eccvec(1)/ecc);
    else
        w_true = 360 - acosd(eccvec(1)/ecc);
    end
    KepElems = [a,ecc,inc,W,w_true,f]';
    % Circular Inclined
elseif inc ~= 0 && ecc == 0
    if rvec(3) > 0
        u_l = acosd(dot(nvec,rvec)/(n*r));
    else
        u_l = 360 - acosd(dot(nvec,rvec)/(n*r));
    end
    KepElems = [a,ecc,inc,W,w,u_l]';
    % Circular Equatorial
elseif inc == 0 && ecc == 0
    if rvec(2) > 0
        lam = acosd(rvec(1)/r);
    else
        lam = 360 - acosd(rvec(1)/r);
    end
    KepElems = [a,ecc,inc,W,w,lam]';
    % All other orbits
else
    KepElems = [a,ecc,inc,W,w,f]';
end
end

function states = oe2rv(elements,mu,flag)
% This function converts the Keplerian elements to the Cartesian state
% vector. From Fundamentals of Astrodynamics with Applications,...
% David Vallado, 2007., pp.126-127, Algorithm 10.
%
% AC Rogers 14 Feb 2013
%
% Inputs:
% p   = parameter
% ecc = eccentricity
% inc = inclination (in radians)
% O   = right-ascension of the ascending node
% w   = argument of periapsis
% nu  = true anomaly
%
% Outputs:
% r   = 3x1 vector of current Cartesian positions (ECI frame)
% v   = 3x1 vector of current Cartesian velocities (ECI frame)

% Set default value of mu to that of Earth if no value is given
if nargin < 2
    mu = 3.986e5;
end
% Determine if vector of orbit elements is the correct size and return an
% error if it is not
if length(elements) ~= 6
    error('Orbit elements vector is not the correct size.')
end
a = elements(1); ecc = elements(2); inc = elements(3);
O = elements(4); w = elements(5);  f = elements(6);
p = a*(1-ecc^2);

if strcmp(flag,'r') == 1
    rpqw = [p*cos(f)/(1+ecc*cos(f));
        p*sin(f)/(1+ecc*cos(f));
        0];
    vpqw = [-sqrt(mu/p)*sin(f);
        sqrt(mu/p)*(ecc+cos(f));
        0];
    map = rotation3(-O,flag)*rotation1(-inc,flag)*rotation3(-w,flag);
elseif strcmp(flag,'d') == 1
    rpqw = [p*cosd(f)/(1+ecc*cosd(f));
        p*sind(f)/(1+ecc*cosd(f));
        0];
    vpqw = [-sqrt(mu/p)*sind(f);
        sqrt(mu/p)*(ecc+cosd(f));
        0];
    map = rotation3(-O,flag)*rotation1(-inc,flag)*rotation3(-w,flag);
else
end
r = map*rpqw; v = map*vpqw;
if length(r)~= 3, error('Wrong size Position vector')
end
if length(v)~= 3, error('Wrong size Velocity vector')
end

states = [r(:); v(:)];

end

function R1 = rotation1(theta,flag)
% This function performs a 1-rotation, takes arguments in radians
if strcmp(flag,'d') == 1
    theta = theta*180/pi;
    R1 = [1 0            0         ;
        0 cosd(theta)  sind(theta);
        0 -sind(theta)   cosd(theta)];
elseif strcmp(flag,'r') == 1
    R1 = [1 0            0         ;
        0 cos(theta)  sin(theta);
        0 -sin(theta)   cos(theta)];
else
end
end

function R3 = rotation3(theta,flag)
% This function performs a 3-rotation, takes arguments in radians
if strcmp(flag,'d') == 1
    theta = theta*180/pi;
    R3 = [cosd(theta) sind(theta) 0;
        -sind(theta) cosd(theta)  0;
        0                 0    1];
elseif strcmp(flag,'r') == 1
    R3 = [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta)  0;
        0                 0    1];
else
end
end