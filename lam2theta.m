function [theta, F] = lam2theta(lambda, q1, q2, Tol)

% Calculation of true longitude theta = f + w 
% from mean longitude lambda = M + w
%
% input :
%    lambda = mean longitude = M + w
%    q1 = e * cos(w)
%    q2 = e * sin(w)
%
% output :
%    theta = true longitude = f + w
%    F = eccentric longitude = E + w

eta = sqrt( 1- q1^2 - q2^2 );
beta = 1 / ( eta*(1+eta) );

% Modified Kepler's Equation
F = lambda;
FF = 1;
while abs(FF) > Tol
  FF = lambda - (F-q1*sin(F)+q2*cos(F));
  dFFdF = -(1-q1*cos(F)-q2*sin(F));
  del_F = -FF / dFFdF;
  F = F + del_F;
end;
%while (F < 0)
%  F = F + (2*pi);
%end;
%while (F >= (2*pi))
%  F = F - (2*pi);
%end;

% True Longitude
num = (1+eta)*(eta*sin(F)-q2) + q2*(q1*cos(F)+q2*sin(F));
den = (1+eta)*(eta*cos(F)-q1) + q1*(q1*cos(F)+q2*sin(F));
theta = atan2( num, den );
while (theta < 0)
  theta = theta + (2*pi);
end;
while (theta >= (2*pi))
  theta = theta - (2*pi);
end;

%
if (lambda < 0)
   kk_plus = 0;
   quad_plus = 0;
   while (lambda < 0)
     kk_plus = kk_plus + 1;
     lambda = lambda + (2*pi);
   end;
   if (lambda < (pi/2)) & (theta > (pi))
      quad_plus = 1;
   elseif (theta < (pi/2)) & (lambda > (pi))
      quad_plus = -1;
   end;
   theta = theta - (kk_plus+quad_plus)*(2*pi);
else
   kk_minus = 0;
   quad_minus = 0;
   while (lambda >= (2*pi))
     kk_minus = kk_minus + 1;
     lambda = lambda - (2*pi);
   end;
   if (lambda < (pi/2)) & (theta > (pi))
      quad_minus = -1;
   elseif (theta < (pi/2)) & (lambda > (pi))
      quad_minus = 1;
   end;
   theta = theta + (kk_minus+quad_minus)*(2*pi);
end;
%theta = theta - (kk_plus+quad_plus)*(2*pi) + (kk_minus+quad_minus)*(2*pi);

end
