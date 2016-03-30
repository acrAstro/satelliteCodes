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
            
        end
        
        function TBP = setInitialConditions(TBP)
            % Takes Kepler elements and returns R and V for ode45
            [R,V] = oe2rv(TBP.kepElems,TBP.mu,'r');
            TBP.initialConditions = [R(:); V(:)];
            
            function [r,v] = oe2rv(elements,mu,flag)
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
                
            end
            
            function [R,V] = HillFrameToInertial(R1,V1,X)
                X = X(:);
                Rot = ijk_to_LVLH([R1(:);V1(:)]);
                
                R = Rot'*X(1:3);
                V = Rot'*X(4:6);
            end
            
            function RSW = ijk_to_LVLH(ic)
                
                % Converts the inertial frame to the orbit frame (Local-Vertical,
                % Local-Horizontal). Page 163, 172 of Fundamentals of Astrodynamics with
                % Applications, David Vallado, 2007.
                
                rvec = [ic(1) ic(2) ic(3)]'; vvec = [ic(4) ic(5) ic(6)]';
                r_hat = (1/norm(rvec)).*rvec;
                w_hat = (1/(norm(cross(rvec,vvec)))).*cross(rvec,vvec);
                s_hat = cross(w_hat,r_hat);
                
                RSW = [r_hat(:), s_hat(:), w_hat(:)]';
                
                % Make sure it's a good rotation matrix. Needs full rank and unity
                % determinant
                if round(det(RSW)) ~= 1 || rank(RSW) ~= 3
                    error('Inappropriately-sized matrix')
                else
                end
            end
            
            function R2 = rotation2(theta,flag)
                % This function performs a 2-rotation, takes arguments in radians
                if strcmp(flag,'d') == 1
                    theta = theta*180/pi;
                    R2 = [cosd(theta) 0 -sind(theta)
                        0 1 0
                        sind(theta) 0 cosd(theta)];
                elseif strcmp(flag,'r') == 1
                    R2 = [cos(theta) 0 -sin(theta)
                        0 1 0
                        sin(theta) 0 cos(theta)];
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
        end
        
        function TBP = getKepElems(TBP)
            % Take R and V and return Kepler elements over time
            oscElems = zeros(6,length(TBP.time));
            for ii = 1:length(TBP.time)
                oscElems(:,ii) = rv2oe(TBP.states(ii,:),TBP.mu);
            end
            TBP.OsculatingElements = oscElems;
            
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