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