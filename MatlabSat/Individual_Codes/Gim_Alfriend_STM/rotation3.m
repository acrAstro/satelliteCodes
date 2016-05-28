function R3 = rotation3(theta,flag)
% This function performs a 3-rotation, takes arguments in radians
if strcmp(flag,'d')
    theta = theta*180/pi;
    R3 = [cosd(theta) sind(theta) 0;
        -sind(theta) cosd(theta)  0;
        0                 0    1];
elseif strcmp(flag,'r')
    R3 = [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta)  0;
        0                 0    1];
else
end
end