function isInsideEarth(a,ecc,Req,safetyAltitude)

% Function to test if an orbit will pass inside safetyAltitude

if Req/1e3 < 10
    % Radius at periapsis
    rp = a*(1-ecc);
    periapsisAltitude = rp - Req;
    fprintf(['Periapsis altitude = ' num2str(periapsisAltitude) '\n']);
    
    safeAltRatio = (safetyAltitude + Req)/Req;
    rpRatio = (rp)/Req;
    if rpRatio <= safeAltRatio
        error(['Orbit periapsis altitude is ' num2str(periapsisAltitude) ' km, check eccentricity and semi-major axis']);
    else
    end
elseif Req/1e3 > 10
    rp = a*(1-ecc);
    periapsisAltitude = rp - Req;
    fprintf(['Periapsis altitude = ' num2str(periapsisAltitude) '\n']);
    
    safeAltRatio = (safetyAltitude + Req)/Req;
    rpRatio = (rp)/Req;
    if rpRatio <= safeAltRatio
        error(['Orbit periapsis altitude is ' num2str(periapsisAltitude) ' m, check eccentricity and semi-major axis']);
    else
    end
else
end
end