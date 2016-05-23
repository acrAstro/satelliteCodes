classdef OrbitPlotter < handle
    
    properties
        States
        StatesQ
        Controls
        ControlsQ
        Points
        Time
        
        plotLinesStates
        plotLineStatesMods
        plotLineStatesSizes
        
        plotLinesStatesQ
        plotLineStatesQMods
        plotLineStatesQSizes
        
        plotLegends
        plotIdentifier
        plotTitle
        plotLabelX
        plotLabelY
        plotLabelZ
        axisBounds
        
        numInput
        spaceshuttle
        fig3d
    end
    
    methods
        function obj = OrbitPlotter(inputStruct)
            for ii = 1:length(inputStruct.times)
                obj.Time{ii}            = inputStruct.times{ii};
            end
            for ii = 1:length(inputStruct.states.states)
                obj.States{ii}          = inputStruct.states.states{ii};
                obj.plotIdentifier{ii}  = inputStruct.id{ii};
            end
            for ii = 1:length(inputStruct.states.statesq)
                obj.StatesQ{ii}         = inputStruct.states.statesq{ii};
            end
            for ii = 1:length(inputStruct.controls.controls)
                obj.Controls{ii}        = inputStruct.controls.controls{ii};
            end
            for ii = 1:length(inputStruct.controls.controlsq)
                obj.ControlsQ{ii}       = inputStruct.controls.controlsq{ii};
            end
            for ii = 1:length(inputStruct.lines.linestates)
                obj.plotLinesStates{ii}         = inputStruct.lines.linestates{ii};
                obj.plotLineStatesMods{ii}      = inputStruct.lines.linemods{ii};
                obj.plotLineStatesSizes(ii)     = inputStruct.lines.linesizes(ii);
            end
            for ii = 1:length(inputStruct.lines.linestatesq)
                obj.plotLinesStatesQ{ii}        = inputStruct.lines.linestatesq{ii};
                obj.plotLineStatesQMods{ii}     = inputStruct.lines.linemodsq{ii};
                obj.plotLineStatesQSizes(ii)    = inputStruct.lines.linesizesq(ii);
            end
            for ii = 1:length(inputStruct.legends)
                obj.plotLegends{ii}             = inputStruct.legends{ii};
            end
            if ~isempty(inputStruct.controls.controls)
                for ii = 1:length(inputStruct.controls.controls)
                    obj.numInput{ii}            = size(inputStruct.controls.controls{ii},1);
                end
            else obj.numInput{ii} = 0;
            end
            
            obj.plotTitle  = inputStruct.title;
            obj.axisBounds = inputStruct.bounds;
            obj.plotLabelX = inputStruct.labels{1};
            obj.plotLabelY = inputStruct.labels{2};
            obj.plotLabelZ = inputStruct.labels{3};
            if strcmp(inputStruct.shuttleFlag,'yes')
                obj.spaceshuttle = inputStruct.shuttleFlag;
            else
                obj.spaceshuttle = [];
            end
        end
        
        function obj = plot3DOrbit(obj)
            obj.fig3d = figure;
            hold on
            grid on
            for ii = 1:length(obj.States)
                X = obj.States{ii};
                if strcmp(obj.plotIdentifier{ii},'GASTM') == 1
                    plot3(X(1,:),X(3,:),X(5,:),obj.plotLinesStates{ii},obj.plotLineStatesMods{ii},obj.plotLineStatesSizes(ii));
                else
                    plot3(X(1,:),X(2,:),X(3,:),obj.plotLinesStates{ii},obj.plotLineStatesMods{ii},obj.plotLineStatesSizes(ii));
                end
            end
            for ii = 1:length(obj.StatesQ)
                X = obj.StatesQ{ii};
                U = obj.ControlsQ{ii};
                if strcmp(obj.plotIdentifier{ii},'GASTM') == 1
                    quiver3(X(1,:),X(3,:),X(5,:),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                    switch obj.numInput{ii}
                        case isempty(obj.numInput{ii})
                            nothing();
                        case 2
                            quiver3(X(:,1),X(:,3),X(:,5),zeros(size(X(:,1))),U(:,1),U(:,2),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                        case 3
                            quiver3(X(:,1),X(:,3),X(:,5),U(:,1),U(:,2),U(:,3),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                    end
                else
                    switch obj.numInput{ii}
                        case isempty(obj.numInput{ii})
                            nothing();
                        case 2
                            quiver3(X(:,1),X(:,2),X(:,3),zeros(size(X(:,1))),U(:,1),U(:,2),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                        case 3
                            quiver3(X(:,1),X(:,2),X(:,3),U(:,1),U(:,2),U(:,3),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                    end
                    plot3(X(1,:),X(2,:),X(3,:),obj.plotLinesStatesQ{ii},obj.plotLineStatesQMods{ii},obj.plotLineStatesQSizes(ii));
                end
            end
            
            legString = cell(length(obj.plotLegends),1);
            for ii = 1:length(obj.plotLegends)
                legString{ii} = [obj.plotLegends{ii}];
            end
            leg = legend(legString);
            title1 = title(obj.plotTitle);
            xl = xlabel(obj.plotLabelX);
            yl = ylabel(obj.plotLabelY);
            zl = zlabel(obj.plotLabelZ);
            if ~isempty(obj.spaceshuttle)
                makeSpaceShuttle();
            else
            end
            if strcmp(obj.axisBounds,'tight') == 1
                axis tight
            else
                axis(obj.axisBounds)
            end
            set([title1,xl,yl,zl,leg],'interpreter','latex','fontsize',12)
        end
    end
end

function nothing(~,~)
end

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