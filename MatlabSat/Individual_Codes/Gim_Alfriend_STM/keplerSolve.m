function [E, f] = keplerSolve(e, M, Tol1);

% Kepler's Equation
E = M;
FF = 1;
while abs(FF) > Tol1
    FF = M - (E - e*sin(E));
    dFFdE = -(1 - e*cos(E));
    del_E = -FF / dFFdE;
    E = E + del_E;
end
while (E < 0)
    E = E + (2*pi);
end
while (E >= (2*pi))
    E = E - (2*pi);
end

%
kk_plus = 0;
while (M < 0)
    kk_plus = kk_plus + 1;
    M = M + (2*pi);
end
kk_minus = 0;
while (M >= (2*pi))
    kk_minus = kk_minus + 1;
    M = M - (2*pi);
end

% True Anomaly
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
if (E >= 0) && (E <= (pi))
    f = abs(f);
else
    f = (2*pi) - abs(f);
end;

f = f - kk_plus*(2*pi) + kk_minus*(2*pi);

return;
