classdef OrbitPlotter < handle
    
    properties
        States
        Controls
        Time
        
        plotLegends
        plotLines
        plotLineMods
        plotLineSizes
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
            if length(inputStruct.states) ~= length(inputStruct.legends)
                error('Each state history should have a matching legend!')
            else
            end
            obj.States          = {};
            obj.Controls        = {};
            obj.Time            = {};
            obj.plotLines       = {};
            obj.plotIdentifier  = {};
            obj.plotLegends     = {};
            obj.plotLineMods    = cell(size(inputStruct.lineMods));
            obj.numInput        = {};
            for ii = 1:length(inputStruct.states)
                obj.Time{ii}            = inputStruct.times{ii};
                obj.States{ii}          = inputStruct.states{ii};
                obj.Controls{ii}        = inputStruct.controls{ii};
                obj.plotLines{ii}       = inputStruct.lines{ii};
                obj.plotLegends{ii}     = inputStruct.legends{ii};
                obj.plotIdentifier{ii}  = inputStruct.id{ii};
                obj.plotLineMods{ii}    = inputStruct.lineMods{ii};
                obj.plotLineSizes(ii)   = inputStruct.lineSizes(ii);
                if ~isempty(inputStruct.controls{ii})
                    obj.numInput{ii}        = size(inputStruct.controls{ii},1);
                else obj.numInput{ii} = 0;
                end
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
                U = obj.Controls{ii};
                if strcmp(obj.plotIdentifier{ii},'GASTM') == 1
                    plot3(X(1,:),X(3,:),X(5,:),obj.plotLines{ii},obj.plotLineMods{ii},obj.plotLineSizes(ii));
                elseif strcmp(obj.plotIdentifier{ii},'0') == 1
                    plot3(X(1),X(2),X(3),obj.plotLines{ii},obj.plotLineMods{ii},obj.plotLineSizes(ii));
                elseif strcmp(obj.plotIdentifier{ii},'q') == 1
                    switch obj.numInput{ii}
                        case isempty(obj.numInput{ii})
                            nothing();
                        case 2
                            quiver3(X(:,1),X(:,2),X(:,3),zeros(size(X(:,1))),U(:,1),U(:,2),obj.plotLines{ii},obj.plotLineMods{ii},obj.plotLineSizes(ii))
                        case 3
                            quiver3(X(:,1),X(:,2),X(:,3),U(:,1),U(:,2),U(:,3),obj.plotLines{ii},obj.plotLineMods{ii},obj.plotLineSizes(ii))
                    end
                else
                    plot3(X(1,:),X(2,:),X(3,:),obj.plotLines{ii},obj.plotLineMods{ii},obj.plotLineSizes(ii));
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
            if strcmp(obj.axisBounds,'tight') == 1
                axis tight
            else
                axis(obj.axisBounds)
            end
            set([title1,xl,yl,zl,leg],'interpreter','latex','fontsize',12)
        end
        
        function obj = plotControls(obj)
            
        end
    end
end

function nothing(~,~)
end