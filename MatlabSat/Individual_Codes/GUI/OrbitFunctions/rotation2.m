function R2 = rotation2(theta,flag)
% This function performs a 2-rotation, takes arguments in radians
if strcmp(flag,'d')
    theta = theta*180/pi;
    R2 = [cosd(theta) 0 -sind(theta)
        0 1 0
        sind(theta) 0 cosd(theta)];
elseif strcmp(flag,'r')
    R2 = [cos(theta) 0 -sin(theta)
        0 1 0
        sin(theta) 0 cos(theta)];
else
end
end