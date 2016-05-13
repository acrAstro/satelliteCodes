function makeSpaceShuttle
plotShuttle(0,0,0,0,-pi/2,pi/2,0.0487,1e-3,[1,1,0.5])
end

function plotShuttle(x,y,z,pitch,roll,yaw,scale_factor,step,cv)
load shuttle;
V = [-V(:,2) V(:,1) V(:,3)];
V(:,1) = V(:,1)-round(sum(V(:,1))/size(V,1));
V(:,2) = V(:,2)-round(sum(V(:,2))/size(V,1));
V(:,3) = V(:,3)-round(sum(V(:,3))/size(V,1));

correction = max(abs(V(:,1)));
V = V./(scale_factor*correction);
ii = length(x);
resto = mod(ii,step);

y = y;
z = z;
pitch = pitch;
roll = roll;
yaw = -yaw;

for jj = 1:step:(ii-resto)
    theta = pitch(jj);
    phi = -roll(jj);
    psi = yaw(jj);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tbe = [cos(psi)*cos(theta), -sin(psi)*cos(theta), sin(theta);
        cos(psi)*sin(theta)*sin(phi)+sin(psi)*cos(phi) ...
        -sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) ...
        -cos(theta)*sin(phi);
        -cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi) ...
        sin(psi)*sin(theta)*cos(phi)+cos(psi)*sin(phi) ...
        cos(theta)*cos(phi)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vnew = V*Tbe;
    rif = [x(jj) y(jj) z(jj)];
    X0 = repmat(rif,size(Vnew,1),1);
    Vnew = Vnew + X0;
    p = patch('faces', F, 'vertices' ,Vnew);
    set(p, 'facec', cv);
    set(p, 'EdgeColor','none');
    H1 = light('Position',[-100 0 0],'Style','local');
    hold on
%     lighting phong
    daspect([1 1 1]) ;
    
end
end