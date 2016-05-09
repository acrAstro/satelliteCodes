function [lambda] = theta2lam(a, theta, q1, q2)

% Calculation of mean longitude lambda = M + w 
% from true longitude theta = f + w
%
% input :
%    a = semi major axis
%    theta = true longitude
%    q1 = e * cos(w)
%    q2 = e * sin(w)
%
% output :
%    lambda = mean longitude

eta = sqrt(1- q1^2 - q2^2);
beta = 1/(eta*(1+eta));
R = (a*eta^2) / (1+q1*cos(theta)+q2*sin(theta));

numerator = R*(1+beta*q1^2)*sin(theta) - beta*R*q1*q2*cos(theta) + a*q2;
denominator = R*(1+beta*q2^2)*cos(theta) - beta*R*q1*q2*sin(theta) + a*q1;

F = atan2(numerator, denominator);
lambda = F - q1*sin(F) + q2*cos(F);
while (lambda < 0)
  lambda = lambda + (2*pi);
end;
while (lambda >= (2*pi))
  lambda = lambda - (2*pi);
end;

if (theta < 0)
   kk_plus = 0;
   quad_plus = 0;
   while (theta < 0)
     kk_plus = kk_plus + 1;
     theta = theta + (2*pi);
   end;
   if (theta < (pi/2)) && (lambda > (pi))
      quad_plus = 1;
   elseif (lambda < (pi/2)) && (theta > (pi))
      quad_plus = -1;
   end;
   lambda = lambda - (kk_plus+quad_plus)*(2*pi);

else
   kk_minus = 0;
   quad_minus = 0;
   while (theta >= (2*pi))
     kk_minus = kk_minus + 1;
     theta = theta - (2*pi);
   end;
   if (theta < (pi/2)) && (lambda > (pi))
      quad_minus = -1;
   elseif (lambda < (pi/2)) && (theta > (pi))
      quad_minus = 1;
   end;
   lambda = lambda + (kk_minus+quad_minus)*(2*pi);
end;


end