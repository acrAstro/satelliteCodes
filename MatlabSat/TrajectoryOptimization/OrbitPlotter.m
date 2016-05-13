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
        
        plotLinePoints
        plotLinePointsMods
        plotLinePointsSizes
        
        plotLegends
        plotIdentifier
        plotTitle
        plotLabelX
        plotLabelY
        plotLabelZ
        axisBounds
        
        numInput
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
            for ii = 1:length(inputStruct.states.points)
                obj.Points{ii}          = inputStruct.states.points{ii};
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
            for ii = 1:length(inputStruct.lines.points)
                obj.plotLinePoints{ii}          = inputStruct.lines.points{ii};
                obj.plotLinePointsMods{ii}      = inputStruct.lines.pointmods{ii};
                obj.plotLinePointsSizes(ii)     = inputStruct.lines.pointsizes(ii);
            end
            for ii = 1:length(inputStruct.legends)
                obj.plotLegends{ii}     = inputStruct.legends{ii};
            end
            if ~isempty(inputStruct.controls.controls)
                for ii = 1:length(inputStruct.controls.controls)
                    obj.numInput{ii}  = size(inputStruct.controls.controls{ii},1);
                end
            else obj.numInput{ii} = 0;
            end
            
            obj.plotTitle  = inputStruct.title;
            obj.axisBounds = inputStruct.bounds;
            obj.plotLabelX = inputStruct.labels{1};
            obj.plotLabelY = inputStruct.labels{2};
            obj.plotLabelZ = inputStruct.labels{3};
        end
        
        function obj = plot3DOrbit(obj)
            figure
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
            for ii = 1:length(obj.Points)
                X = obj.Points{ii};
                plot3(X(1),X(2),X(3),obj.plotLinePoints{ii},obj.plotLinePointsMods{ii},obj.plotLinePointsSizes(ii));
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
            if strcmp(obj.axisBounds,'tight') == 1
                axis tight
            else
                axis(obj.axisBounds)
            end
            set([title1,xl,yl,zl,leg],'interpreter','latex','fontsize',12)
        end
        
        function obj = plotControls(obj)
            % I've standardized how you plot control signals; if you're
            % fully actuated, the control plot looks for three signals, if
            % you're under-actuated, the control plotter looks for only y
            % and z directions
            method = obj.numInput;
            switch method
                case 2
                    figure
                    subplot(212)
                    hold on
                    grid on
                    for ii = 1:length(obj.Controls)
                        U = obj.Controls{ii};
                        T = obj.Time{ii};
                        plot(T,U(1,:),obj.plotLinesStates{ii},obj.plotLineStatesMods{ii},obj.plotLineStatesSizes{ii});
                    end
                case 3
            end
        end
    end
end

function nothing(~,~)
end