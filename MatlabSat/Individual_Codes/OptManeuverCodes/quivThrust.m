function [Xout,Uout] = quivThrust(t,x,u,nu,step)

X = x(:,1:3); U = u(:,1:nu);
for ii = 1:length(t)
    if mod(ii,step) == 0
        Xout(ii,:) = X(ii,:);
        Uout(ii,:) = -U(ii,:);
    else
    end
end
end